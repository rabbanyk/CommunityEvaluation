/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Copyright © 2009  Peter Ronhovde                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 *  This file is part of the RNMRA program.                                    *
 *                                                                             *
 *  RNMRA is free software: you can redistribute it and/or modify              *
 *  it under the terms of the GNU General Public License as published by       *
 *  the Free Software Foundation, either version 3 of the License, or          *
 *  (at your option) any later version.                                        *
 *                                                                             *
 *  This program is distributed in the hope that it will be useful,            *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 *  GNU General Public License for more details.                               *
 *                                                                             *
 *  You should have received a copy of the GNU General Public License          *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 *  Created by Peter Ronhovde on 10/15/09 (email: ronhovde@hbar.wustl.edu)     *
 *	Modified on 11/12/09                                                       *
 *  Location: Washington University in St. Louis, St. Louis, MO 63101          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
/*
Peter Ronhovde
August 2005
Optimized memory handling done November 2005
Additional functions added November 2005
MSL - My (or Math) Scientific Library - Array1D
See MSL_Array.h for general comments and some documentation.

Array1D Specific Notes:
These notes handle only behavior specific to the Array1D implementation.

1. If StrictErrorChecking != 1 then if the objects are of different sizes,
the assignment operator destroys the original object and recreates one with the
new 'correct' size to make the assignment.

2. The 'Assign' function is a (useful) function that equates the elements of the
Array1D from Initial up to but not including the element arrayName[L]
(i.e. it fills a total of (L-Initial) elements of the array). The changed
array is returned by reference as well as having been passed by reference, so it
can be used either way.


To-Do List:
1. Add a conversion from Array1D back to T*.
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TARRAY1D_H
#define MSL_TARRAY1D_H

//#define DEBUG_MODE

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
#ifndef MSL_TARRAY_H
#include "MSL_Array.h"
#endif

using namespace std;

namespace MSL {

// forward declarations
template <typename T> class Array1D;
template <typename T> Array1D<T> operator*(T x, Array1D<T>& b);
template <typename T> istream& operator>>(istream& fin, Array1D<T>& a);
// prototype some classes
template <typename T> class TList;
//template <typename T> class  TListArray;

#ifndef MSL_TVECTOR_T
#define MSL_TVECTOR_T
// define some common vector types
typedef Array1D<double>         TVectorFloat; // double  array type definition
typedef Array1D<int>            TVectorInt;   // integer array type definition
typedef Array1D<unsigned char>  TVectorChar;  // char array type definition
#endif

//------------------------------------------------------------------------------
template <typename T>
class Array1D : public TVector {
// Memory allocation is dynamic, and this Array1D uses the standard C++ array
// indexing with brackets as in arrayName[index] with zero indexing.
 friend Array1D<T> operator*<T>(T x, Array1D<T>& b);
 friend istream&   operator>><T>(istream& fin, Array1D<T>& a);

 public:
 // overloaded operators
         Array1D& operator=(const Array1D& b);
         Array1D& operator=(TList<T>& b); 
         Array1D& operator+=(const Array1D& b);
         Array1D& operator-=(const Array1D& b);
         Array1D& operator*=(T x);        // *= by scalar
         Array1D& operator/=(T x);        // /= by scalar
         Array1D  operator-();            // '-' negation
         Array1D  operator+(const Array1D& a);
         Array1D  operator-(const Array1D& a);
         Array1D  operator*(T x);         // * by scalar
         Array1D  operator/(T x);         // / by scalar
  //inline T&     operator[](unsigned i); // standard C/C++ array reference
  inline T&       operator[](unsigned i) const;
  // Public member functions
         void     init(T Value = 0);      // Reinitialize the Array to a value
         // initialize with an increment
         void     initStep(T begin = 0, T step = 1, bool bRandomizeOrder = 0);
  // vector-like functions         
         T        dot(Array1D& b) const;
         double   getMag() const;
         T        getSum() const;
         T        getMax(bool bAbsoluteValue = 0) const;
         T        getMin(bool bAbsoluteValue = 0) const;
         double   getSigma(bool SigmaN = 0) const;
         double   getAvg()  const;
  inline int      getSize() const { return L; };
  inline int      getL()    const { return L; };
  // utility functions
  inline  bool    isValid(bool bEmptyOK = 1) const;
  inline  bool    isValid(bool bEmptyOK = 1);
          void    randomizeOrder();       // randomize the order of elements
          // bFakeResize adjusts the current vector length without returning
          // memory.  This is useful for insert, drop and add operations.
          void    resize(unsigned Length, bool Fill = 0, T FillValue = 0, 
                         bool bFakeResize = 0);
          bool    add(T val); // add a value to the end of the list (if L<LMax)
  inline  bool    drop() { L--; };  // drop the last value on end of list
          bool    drop(unsigned element);  // drop the last value on end of list
          bool    insert(T val, bool bCheckDuplicates=0);
          //bool  isMember(T val, bool bBinarySearch=0);
          bool    isMember(T val);  // performs linear search
          void    assign(Array1D& b, unsigned Length, unsigned Initial = 0);
          void    assign(T* a, unsigned Length, unsigned Initial=0, 
                         bool Construct=0);
         ostream& display(ostream& fout, bool bColor=1);
  // Standard public member functions
  Array1D(unsigned Length = 0);                         // bare constructor
  Array1D(unsigned Length, T Value);                    // init constructor
  Array1D(unsigned Length, T* a, unsigned Initial = 0); // init with *array
  Array1D(const Array1D &b);                            // copy constructor
  ~Array1D();                                           // destructor
 protected:
  // Data elements
  unsigned  L, LMax; // number of elements in the array indexed 0..(L-1)
  T*        data;    // actual data array: a standard 1D dynamic array
}; // end class Array1D;

// Take the unusual step of including TList<T> after the class definition since
// TList<T> needs the Array1D<T> functions defined, but the Array1D<T> functions
// need a full declaration of the TList<T> type for the operator= overload.
} // must end MSL namespace temporarily otherwise a sub-namespace is declared
#ifndef MSL_TLIST_H
#include "MSL_List_Template.h"
#endif
namespace MSL {  // restart namespace

//---------------------Array1D<T> Member Function Declarations------------------

template <typename T> 
Array1D<T>::Array1D(unsigned Length) { // constructor without an init value
  // Raw constructor that skips the initialization
  // Without the initialization, the user takes the chance of accessing junk,
  // but it is provided for "more efficient" applications.
  // If used as the default constructor, it generates a Array1D object that is 
  // a zero length but invalid Array1D object.
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Constructing TArray1D<T> object with ");
  #endif
  L = Length;  LMax = Length;  // LMax is used for 'fast' resizes
  if(L>0) {
    #ifdef DEBUG_MODE_HIGH
    debugMsg("size "+itos(L)+"... ");
    #endif
    data = new T[L];
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } // end un-initialized array creation
  else if(L==0) {
    #ifdef DEBUG_MODE_HIGH
    debugMsg("default size = 0... ");
    #endif
    #ifdef DEBUG_MODE
    ValidArray = 0;  
    #endif
    data = NULL;
  } // end zero length if
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting\n");
  #endif
}; // end Array1D default and no-initialization constructor

template <typename T>
Array1D<T>::Array1D(unsigned Length, T Value) {// constructor with an init value
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Constructing TArray1D<T> object with init... ");
  #endif
  L = Length;  LMax = Length;
  if(L>0) {
    data = new T[L];
    T* pElement = data;
    for(unsigned i=0; i<L; i++) { *pElement = Value;  pElement++; } // initialize
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } // end un-initialized array creation
  else if(L==0) {
    data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    warningMsg("Attempting to create a size 0 Array1D object with an init value");
    #endif
  } // end zero length if
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting");
  #endif
}; // end Array1D constructor

template <typename T>
Array1D<T>::Array1D(const Array1D &b) { // copy constructor
  #ifdef DEBUG_MODE
  debugMsg("Copying TArray1D<T> object... ");
  if(!b.isValid())  
    errorMsg("Invalid Array1D<T> source object in copy constructor"); 
  #endif
  if(b.L>0) { // check for valid array
    L = b.L;
    data = new T[L];  // allocate copy's own dynamic memory
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=0; i<b.L; i++) { 
      *pElement = (*pbElement);  pElement++;  pbElement++; } // initialize
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } else { 
    data = NULL; 
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else for L=0
  #ifdef DEBUG_MODE
  debugMsg("exiting\n");
  #endif
}; // end Array1D copy constructor

template <typename T>
Array1D<T>::Array1D(unsigned Len, T* a, unsigned Initial) { // initializing with T* array
// This constructor requires a little trust of the regular C pointer array
// object *a and Len(gth) being passed to it.  Sadly, this is somewhat
// unavoidable (and part of the reason for classes being created in the first
// place), but it is a constructor that should be included.
// Len is the *Length* of the Array1D (as opposed to the final index) where

// Initial is the initial *index* value from which to start the copy which will
// normally be zero so that the entire array is copied.
// i.e. Construct an Array1D starting at element Initial with a Length of Len.
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Constructing TArray1D<T> object with C array init... ");
  #endif
  L = Len-Initial;
  data = new T[L];  // allocate new object's own dynamic memory
  T* pElement  = data;
  T* paElement = a;
  for(unsigned i=Initial; i<L; i++) { 
    *pElement = (*paElement);  pElement++;  paElement++; } // initialize
  #ifdef DEBUG_MODE
  ValidArray = Valid1DArrayFlag;
  #endif
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting\n");
  #endif
};

template <typename T>
Array1D<T>::~Array1D() {
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Destructing TArray1D<T> object with size = "+itos(L)
           //+" and data pointer = "+itos((int)data)+"... ",brown);
           +" and data pointer = "+itos(-1)+"... ",brown);
  #endif
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("Error Destructing Array1D");
  #endif
  delete[] data;  data = NULL;  // just to be safe
  L = -1;  LMax = -1;
  #ifdef DEBUG_MODE
  ValidArray = 0;
  #endif
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting\n",brown);
  #endif
}; // end Array1D destructor


template <typename T>
inline bool Array1D<T>::isValid(bool bEmptyOK) { 
  // bEmptyOK is for use in the destructor and resize() functions where 
  // destructing an (consistent) empty array is still a valid operation
  #ifdef DEBUG_MODE
  bool bResult = 0;
  if(bEmptyOK && ValidArray==0 && L==0 && data==NULL)  return 1;
  else if(ValidArray!=Valid1DArrayFlag)  
        warningMsg("Invalid Array1D<T> object with flag = "+itos(ValidArray)); 
  else if(((data==NULL && L!=0) || (data!=NULL && L==0)) && !(L<LMax && L==0))  
        // fake size invalidates this check
        warningMsg("Inconsistent Array1D<T> object data pointer = "
                   //+itos((int)data)+" with a size of "+itos(L)); 
                   +itos(-1)+" with a size of "+itos(L)); 
  else if( (!(ValidArray==Valid1DArrayFlag && data!=NULL && L!=0) ) && !(L<LMax && L==0))  
        warningMsg("Invalid Array1D<T> object flag.  Flag = "+itos(ValidArray)+
                //", data pointer = "+itos((int)data)+" with a size of "+itos(L)); 
                ", data pointer = "+itos(-1)+" with a size of "+itos(L)); 
  else  bResult = 1;  // if we make it to here, it is ok
  return bResult;
  #else
  return 1;
  #endif
};     // valid array check
// const version of isValid
template <typename T>
inline bool Array1D<T>::isValid(bool bEmptyOK) const {
  // bEmptyOK is for use in the destructor and resize() functions where 
  // destructing an (consistent) empty array is still a valid operation
  #ifdef DEBUG_MODE
  bool bResult = 0;
  if(bEmptyOK && ValidArray==0 && L==0 && data==NULL)  return 1;
  else if(ValidArray!=Valid1DArrayFlag)  
        warningMsg("Invalid Array1D<T> object with flag = "+itos(ValidArray)); 
  else if((data==NULL && L!=0) || (data!=NULL && L==0))  
        warningMsg("Inconsistent Array1D<T> object data pointer = "
                   //+itos((int)data)+" with a size of "+itos(L)); 
                   +itos(-1)+" with a size of "+itos(L)); 
  else if( !(ValidArray==Valid1DArrayFlag && data!=NULL && L!=0) )  
        warningMsg("Invalid Array1D<T> object flag.  Flag = "+itos(ValidArray)+
                //", data pointer = "+itos((int)data)+" with a size of "+itos(L)); 
                ", data pointer = "+itos(-1)+" with a size of "+itos(L)); 
  else  bResult = 1;  // if we make it to here, it is ok
  return bResult;
  #else
  return 1;
  #endif
}; // end Array1D<T>::isValid


// Other useful array functions
template <typename T>
void Array1D<T>::init(T Value) {
  // Reinitialize the Array.  Though it doesn't really matter too much, it is
  // assumed that this is a reinitialization since one of the constructors can be
  // used for the first initialization.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in init()");
  #endif
  T* pElement = data;
  for(unsigned i=0; i<L; i++) { *pElement = Value;  pElement++; }
  return;
};

template <typename T>
void Array1D<T>::initStep(T begin, T step, bool bRandomOrder) {
  //initialize given an starting number and a step size
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in initStep()"); 
  #endif
  T* pElement = data;
  for(unsigned i=0; i<L; i++) { *pElement = begin + step*i;  pElement++; }
  if(bRandomOrder)  randomizeOrder();
  return;
}; // end init range

template <typename T>
void Array1D<T>::randomizeOrder() {
  //initialize given an integer range
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in randomizeOrder()"); 
  #endif
  unsigned  iS, jS;
  T    temp;
  T   *pk, *piS, *pjS;  // micro-optimization pointers
  // first, ensure that every spin is moved at least once
  for(unsigned k=0; k<L; k++) {      // swap two nodes at a time
    jS   = randomInt(0,L-1);   // not worried about duplicates
    pjS  = &(data[jS]);  pk  = &(data[k]);
    temp = *pjS;      *pjS = *pk;      *pk = temp;    // swap nodes
  } // end for k
  // then mix things up some more choosing two random nodes at a time
  // I'm not sure if this has any real significance.
/*
  for(unsigned k=0; k<L; k++) {  // swap two nodes at a time
    iS = randomInt(0,L-1);  jS = randomInt(0,L-1);  // duplicates don't matter
    piS = &(data[iS]);  pjS  = &(data[jS]);
    temp = *pjS;      *pjS = *piS;       *piS = temp; // swap nodes
  } // end for k
*/
  return;
}; // end init range
  
template <typename T>
bool  Array1D<T>::add(T val) {
  // end add a value to the end of the list
  //cout << "L = " << L << " and LMax = " << LMax << endl; // debugging
  if(L<LMax) {
    L += 1;
    data[L-1] = val;
  } else
    errorMsg("Array1D object size is maximized.  Cannot add a new value!");
}; // end add a value

template <typename T>
T  Array1D<T>::dot(Array1D& b) const { // dot product of 2 Array1D's
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in dot()"); 
  #endif
  T absum = 0;
  if(L==b.L) {
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=0; i<L; i++) { 
      absum += (*pElement)*(*pbElement);  pElement++;  pbElement++; }
  } // end if isValid
  else  errorMsg("Invalid dot product for Array1D object.\n");
  return absum;
}; // end dot product

template <typename T>
TFloat  Array1D<T>::getMag() const {        // magnitude of current Array1D
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getMag()"); 
  #endif
  TFloat a2sum = 0.0;
  T* pElement  = data;
  for(unsigned i=0; i<L; i++) { a2sum += (*pElement)*(*pElement);  pElement++; }
  return sqrt(a2sum);
}; // end getMag

template <typename T>
T  Array1D<T>::getSum() const {        // magnitude of current Array1D
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getSum()"); 
  #endif
  T  asum = 0;
  T* pElement  = data;
  for(unsigned i=0; i<L; i++) { asum += (*pElement);  pElement++; }
  return asum;
}; // end getMag

template <typename T>
T  Array1D<T>::getMax(bool bAbsoluteValue) const { // get maximum value
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getSum()"); 
  #endif
  T  *pElement = data;
  T  max;
  if(bAbsoluteValue) {
    max = abs(data[0]);
    for(unsigned i=1; i<L; i++) { 
      if(abs(*pElement)>max)  max = abs(*pElement);  
      pElement++;
    } // end for i
  } else{
    max = data[0];
    for(unsigned i=1; i<L; i++) { 
      if((*pElement)>max)  max = (*pElement);  
      pElement++;
    } // end for i
  } // end else
  return max;
}; // end getMax

template <typename T>
T  Array1D<T>::getMin(bool bAbsoluteValue) const { // get minimum value
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getSum()"); 
  #endif
  T  *pElement = data;
  T  min;
  if(bAbsoluteValue) {
    min = abs(data[0]);
    for(unsigned i=1; i<L; i++) { 
      if(abs(*pElement)<min)  min = abs(*pElement);  
      pElement++;
    } // end for i
  } else {
    min = data[0];
    for(unsigned i=1; i<L; i++) { 
      if((*pElement)<min)  min = (*pElement);  
      pElement++;
    } // end for i
  } // end else
  return min;
}; // end getMin

template <typename T>
TFloat  Array1D<T>::getAvg() const { 
  // average of current Array1D<T>
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getAvg()"); 
  #endif
  TFloat a2sum = 0.0;
  T* pElement = data;
  for(unsigned i=0; i<L; i++) { a2sum  += (*pElement);  pElement++; }
  return ( a2sum/(TFloat)L );
}; // end getAvg

template <typename T>
TFloat  Array1D<T>::getSigma(bool useSigmaN) const { 
  // Standard deviation of current Array1D
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in getSigma()"); 
  #endif
  TFloat asum = 0.0, sigsum = 0.0, dev, avg;
  T* pElement = data;
  for(unsigned i=0; i<L; i++) { asum += (*pElement);  pElement++; }
  avg = asum/(TFloat)L;
  pElement = data;
  for(unsigned i=0; i<L; i++) {dev = (*pElement)-avg; sigsum += dev*dev; pElement++;}
  if(useSigmaN)  return sqrt(sigsum/(TFloat)L);       // Sigma_N
  else           return sqrt(sigsum/(TFloat)(L-1));   // Sigma_N-1
}; // end getSigma

template <typename T>
void Array1D<T>::resize(unsigned Len, bool bFill, T FillValue, bool bFake) {
  // bFake performs only a numerical change in the length of the vector
  // leaving the system to handle the actual size upon destruct. 
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("Invalid Array1D<T> object in resize()");
  #endif
  if(bFake && Len<=LMax && Len>=0) { 
    //cout << "Fake resizing vector to " << Len << endl;  // debugging
    L = Len;  return; 
  } // end if bFake

  //if(Len>LMax)  errorMsg("Length exceeds LMax in Array1D?");  // debugging
  //if(Len>L+1)   errorMsg("Length change exceeds +1?");        // debugging

  // now, the logic for a real resize operation
  if(Len>0) {
    T *cdata, *pd, *pc; // declare all pointers used
    cdata = new T[Len];
    pd = data;  pc = cdata;
    unsigned  minL = min(L,Len);
    for(unsigned i=0; i<minL; i++) { *pc = *pd;  pd++; pc++; } // copy current data
    delete[] data;  // delete old data - still safe if data==NULL
    data = cdata;   // now set new data pointer to new array
    L = Len;  LMax = Len;
/*
    T *cdata, *pd, *pc; // declare all pointers used
    cout << "Resizing vector to " << Len << flush;  // debugging
    cout << "a " << flush;  // debugging
    delete[] data;
    cout << "b " << flush;  // debugging
    data = new T[Len];
    cout << "c " << flush;  // debugging
*/
    // now fill the remainder if appropriate
    if(bFill && L<Len) {
      pd = &data[L];
      for(unsigned i=L; i<Len; i++) { *pd = FillValue;  pd++; }
    } // end if Fill
    //cout << "done." << endl; // debugging
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } else { 
    // just get rid of current data and set pointer to NULL
    delete[] data;
    data = NULL;
    L = 0;  LMax = 0;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else for L = 0
  return;
};


template <typename T>
bool Array1D<T>::drop(unsigned element) {
  // swap element with the last value on end of list and drop it from the array
  T  temp;
  temp = data[L-1];  data[L-1] = data[element];  data[element] = temp;
  drop();  // inline drop (in case we change logic later
};  


template <typename T>
bool Array1D<T>::isMember(T val) { //, bool bBinarySearch=0) {
  // bBinarySearch implements a ... binary search, otherwise we search linearly
  unsigned i=0;
  bool     bFound = 0;

  T        *pi = &(data[0]);
  while(i<L && !bFound) {
    if((*pi)==val)  bFound = 1;
    pi++;  i++;
  } // end for i
  return bFound;
}; // end if isMember


template <typename T>
bool Array1D<T>::insert(T val, bool bCheckDuplicates) {
  // Uses an optimized insertion sort logic to insert a single node into *this
  // *sorted* Array1D<T> object.  The intention is that the function is used to 
  // `build' the integer list as it is created rather than to be an iterated 
  // function call in an insertion sort (although that is possible).
  int   j;
  T    *pJ, *pPrevJ;       // pointers to j'th, previous j'th #s
  bool bInserted;
  
  if(L==0) { // just add it to the list
    resize(1,0,0,1);      // fake resize the vector
    data[0] = val;
    return 1;
  } // end length check

  if(bCheckDuplicates) {
    // do this as a separate case to avoid redundant loops
    //cout << "A " << flush;  // debugging
    pPrevJ = &(data[L-1]);
    //cout << "B " << flush;  // debugging
    j = L;                            // j starts at candidate location
    //cout << kij[iN] << endl;  // debugging
    while(j>0 && (*pPrevJ)>val) { pPrevJ--;  j--; } // end while j
    //cout << "C with pPrevJ = " << (*pPrevJ) << flush;  // debugging
    // j ends with the location of the addition

    if((*pPrevJ)==val) {
      //cout << "D " << flush;  // debugging
      bInserted = 0;        // ignore addition of val
    } else {                                 // copy the data now
      resize(L+1,0,0,1);      // fake resize the vector
      pJ = &(data[L]);
      pPrevJ = pJ - 1;
      //cout << "E " << flush;  // debugging
      // move data element up to next spot element by element
      for(int k=L; k>j; k--) {
       *pJ = *pPrevJ;                      // move element up to next spot
        pJ  = pPrevJ;  pPrevJ--;           // decrement j'th, (j-1)'th pointers
      } // end for k
      //cout << "F " << flush;  // debugging
      *pJ = val;
      bInserted = 1;
      //cout << "G " << flush;  // debugging
    } // end else
    //cout << "H " << flush;  // debugging
  } else { // else not bCheckDuplicates
    resize(L+1,0,0,1);  // fake resize of this vector
    pJ = &(data[L]);    // dummy starting location is incremented in loop
    pPrevJ = pJ - 1;    // we need j and j-1 pointer index locations stored
    j = L;              // j and i counters start at same number index
    while(j>0 && (*pPrevJ)>val) {
     *pJ = *pPrevJ;          // move data element up to next spot
      pJ = pPrevJ;  pPrevJ--;// decrement j'th and (j-1)'th pointers
      j--;
    } // end while j
    *pJ = val;  // ok, finally add the integer to the list
    bInserted = 1;
  } // end else bCheckDuplicates
  
  //cout << "I " << flush;  // debugging
  return bInserted;  // was the node actually inserted
}; // end insert (Sorted)


template <typename T>
void Array1D<T>::assign(Array1D& b, unsigned Final, unsigned Initial) {
// Assignment function that will equate the elements of the Array1D from Initial
// up to but not including the element arrayName[L] (i.e. it fills a total
// of (IFinal-Initial) elements of the array).

// Final defaults to the entire length of *this.  INT_MAX is shown because of
// the required constant default parameter.
  #ifdef DEBUG_MODE
  if(!isValid())      errorMsg("Invalid Array1D<T> object in assign()");
  else  if(!b.isValid())  
          errorMsg("Invalid Array1D<T> source object in assign()");
  #endif
  if(Final>=0 && Initial>=0) {
    Final = min(Final,L);
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=Initial; i<Final; i++) {
      (*pElement) = (*pbElement);  pElement++;  pbElement++; }
  } // end if isValid
  return;
}; // end variable range assignment

template <typename T>
void Array1D<T>::assign(T* a, unsigned Final, unsigned Initial, bool Construct_a) {
// NEEDS TO BE REIMPLEMENTED AS A NON-MEMBER FUNCTION in another form
// potentially since the current form could be confused with 
// Assign(Array1D b,...) in that the previous copies *from* b and the current 
// copies *to* a.
                       
// Assignment function that exports *this to a C pointer array *a.
// Final defaults to the entire length of *this.  INT_MAX is shown because of
// the required constant default parameter.
// If *a is usually an existing pointer array, it will be used unless
// Construct_a is set to true.  If Construct_a is true, (*a will be deleted and)
// the memory for *a will be (re)allocated.
  #ifdef DEBUG_MODE
  if(!isValid()) errorMsg("Invalid Array1D<T> object in assign() init version");
  #endif
  if(Final>=0 && Initial>=0) {
    unsigned  aSize = Final-Initial+1;
    Final = min(Final,L);
    if(Construct_a) { delete[] a;  a = new T[aSize]; } // end if
    T* pElement  = &(data[Initial]);
    T* paElement = a;
    for(unsigned i=Initial; i<Final; i++) {
      (*paElement) = (*pElement);  pElement++;  paElement++; }
  } // end if isValid
  return;
}; // end variable range assignment


template <typename T>
inline T&  Array1D<T>::operator[](unsigned i) const { 
  #ifdef DEBUG_MODE
  if(i>=L)  errorMsg("Array1D index is out of bounds! i = "+itos(i)+
                     " in a size "+itos(L)+" vector");// debugging
  #endif
  return data[i]; 
}; // end reference operator


template <typename T>
Array1D<T>& Array1D<T>::operator=(const Array1D& b) { // assignment for 2 arrays
  // check for valid arrays first
  #ifdef DEBUG_MODE
  //if(L!=b.L)  
  //  warningMsg("Array1D '=' sizes do not match. "+itos(L)+" vs "+itos(b.L));
  if(!isValid(1))  errorMsg("Invalid Array1D<T> object in operator=()");
  else if(!b.isValid())  
         errorMsg("Invalid Array1D<T> source object in operator=()");
  #endif
  if(L==b.L && L>0) {
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=0; i<L; i++) { 
      *pElement = (*pbElement);  pElement++;  pbElement++; }
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } // end if
  else if(b.L>0) { // delete current data and copy from b
    // this version of the assignment operator allows size of *this to change
    L = b.L;  LMax = b.LMax;
    if(data!=NULL) delete[] data;
    data = new T[L];  // allocate new size of data
    for(unsigned i=0; i<L; i++) { data[i] = b.data[i]; }
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } // end mutable assignment
  else { 
    if(data!=NULL)  delete[] data;
    L = 0;  LMax = 0;
    data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else for L==0
  return *this;
}; // end =


template <typename T>
Array1D<T>& Array1D<T>::operator+=(const Array1D& b) {   // additive assignment
  #ifdef DEBUG_MODE
  if(L!=b.L)  
    errorMsg("Array1D '=' sizes do not match. "+itos(L)+" vs "+itos(b.L));
  if(!isValid())   errorMsg("Invalid Array1D<T> object in operator+=()");
  else  if(!b.isValid())  
                   errorMsg("Invalid Array1D<T> source object in operator+=()");
  #endif
  // check for valid arrays first
  if(isValid() && b.isValid() && L==b.L) {
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=0; i<L; i++) { 
      *pElement += (*pbElement);  pElement++;  pbElement++; }
  } // end if isValid
  else cerr << "Invalid '+=' assignment for Array1D.\n";
  return *this;
}; // end +=


template <typename T>
Array1D<T>& Array1D<T>::operator-=(const Array1D& b) {   // subtraction assignment
  #ifdef DEBUG_MODE
  if(L!=b.L)  
    errorMsg("Array1D '=' sizes do not match. "+itos(L)+" vs "+itos(b.L));
  if(!isValid())   errorMsg("Invalid Array1D<T> object in operator-=()");
  else  if(!b.isValid())  
                   errorMsg("Invalid Array1D<T> source object in operator-=()");
  #endif
  // check for valid arrays first
  if(L==b.L) {
    T* pElement  = data;
    T* pbElement = b.data;
    for(unsigned i=0; i<L; i++) { 
      *pElement -= (*pbElement);  pElement++;  pbElement++; }
  } // end if isValid
  return *this;
}; // end -=


template <typename T>
Array1D<T>& Array1D<T>::operator*=(T x) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator*=(T)");
  #endif
  T* pElement = data;
  for(unsigned i=0; i<L; i++) { *pElement *= x;  pElement++; }
  return *this;
}; // end *=

template <typename T>
Array1D<T>& Array1D<T>::operator/=(T x) {   // T division assignment
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator/=(T)");
  #endif
  if(x!=0.0) { 
    T* pElement = data;
    for(unsigned i=0; i<L; i++) { *pElement /= x;  pElement++; }
  } // end if x!=0
  else  errorMsg("Division by zero in Array1D scale operation!\n");
  return *this;
}; // end /=

/*
// iffdef STRICT_ERROR_CHECKING_P
template <typename T>
inline T& Array1D<T>::operator[](unsigned i) const {
  // check for valid array first
  if(StrictErrorChecking==1 && isValid() && (i<L) && (i>=0)) return data[i];
  else if(i>=L || i<0) cerr << "Array1D index out of bounds.\n";
  return data[L];  // error return value
}; // end array reference
#else
// standard C/C++ array reference
template <typename T>
inline T& Array1D<T>::operator[](unsigned i) const { return data[i]; };
#endif
*/

template <typename T>
Array1D<T>  Array1D<T>::operator-() {             // unary negative
  Array1D c(L);
  if(isValid()) { 
    T* pElement  = data;
    T* pcElement = c.data;
    for(unsigned i=0; i<L; i++) { 
      *pcElement = -(*pElement);  pElement++;  pcElement++; }
  } // end if isValid
  else cerr << "Invalid array in unuary '-' for Array1D object.\n";
  return c;
}; // end unary -


template <typename T>
Array1D<T> Array1D<T>::operator+(const Array1D& b) { // array addition
  #ifdef DEBUG_MODE
  if(L!=b.L)  
    errorMsg("Array1D '=' sizes do not match. "+itos(L)+" vs "+itos(b.L));
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator+()");
  #endif
  Array1D c(L);
  if(L==b.L) {
    T* pElement  = data;
    T* pbElement = b.data;
    T* pcElement = c.data;
    for(unsigned i=0; i<b.L; i++) { 
      *pcElement = (*pElement) + (*pbElement);  
      pElement++;  pbElement++;  pcElement++;
    } // end for i
  } // end if
  return c;
}; // end +

template <typename T>
Array1D<T> Array1D<T>::operator-(const Array1D& b) { // array subtraction
  #ifdef DEBUG_MODE
  if(L!=b.L)  
    errorMsg("Array1D '=' sizes do not match. "+itos(L)+" vs "+itos(b.L));
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator+()");
  #endif
  Array1D c(L);
  if(L==b.L) {
    T* pElement  = data;
    T* pbElement = b.data;
    T* pcElement = c.data;
    for(unsigned i=0; i<b.L; i++) { 
      *pcElement = (*pElement) - (*pbElement);  
      pElement++;  pbElement++;  pcElement++; 
    } // end for i
  } // end if 
  return c;
}; // end binary -


template <typename T>
Array1D<T> Array1D<T>::operator*(T x) {   // vector scaling vec*x
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator*(T)");
  #endif
  Array1D c(L);
  T* pElement  = data;
  T* pcElement = c.data;
  for(unsigned i=0; i<L; i++) {
    *pcElement = (*pElement)*x;  pElement++;  pcElement++; }
  return c;
}; // end scalar multiplication


template <typename T>
Array1D<T> Array1D<T>::operator/(T x) {   // scalar division
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array1D<T> object in operator/(T)");
  #endif
  Array1D c(L);
  if(x!=0.0) { 
    T* pElement  = data;
    T* pcElement = c.data;
    for(unsigned i=0; i<L; i++) { 
      *pcElement = (*pElement)/x;  pElement++;  pcElement++; }
  } else  errorMsg("Scaling division by zero error on Array1d object.\n");
  return c;
}; // end scalar division
// x/vector is not a valid operation, so it is not defined

template <class T>
ostream&  Array1D<T>::display(ostream& fout, bool bColor) {
  // Outputs the Array1D in MatrixMarket format.  Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  if(isValid()) {
    //fout << MMHeader;
    // must somehow allow comment lines here at some point
    if(bColor)  fout << brown << "Size " << getSize() << ":  " << green;
    else        fout << "Size " << getSize() << ":  ";
    T* paElement = &(data[0]);
    for(int i=0; i<getSize(); i++) { fout << (*paElement) << " "; paElement++; }
    if(bColor)  fout << normText;
  } // end if isValid
  return fout;
};

//------------------------------------------------------------------------------
// End Array1D method definitions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Overloaded iostream operators
template <class T>
ostream& operator<<(ostream& fout, Array1D<T>& a) {
  // Outputs the Array1D in MatrixMarket format.  Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  if(a.isValid()) {
    //fout << MMHeader;
    // must somehow allow comment lines here at some point
    fout << "Size " << a.getSize() << ":  ";
    T* paElement = &(a[0]);
    for(int i=0; i<a.getSize(); i++) {
      fout << (*paElement) << " ";  paElement++; }
  } // end if isValid
  return fout;
};

template <typename T>
void ArrayOut(ostream& fo, Array1D<T> &v, unsigned length,
              string label, int W1, int W2, char commentChar) {
  // output the vector horizontally rather than in MatrixMarket format
  // There is a bug with the streaming return operator.  The return ofstream&
  // is having problems.
  fo << commentChar << setw(W1) << label.c_str();
  for (int i=0; i<min(v.getSize(),length); i++) { fo << setw(W2) << v[i]; }
  fo << endl;
  return;
}; // end ArrayOut

template <typename T>
istream& operator>>(istream& fin, Array1D<T>& a) {
  unsigned M, N;  // cols and rows - here M = 1 since it is an Array1D
  unsigned nComments=0;
  string sInput;
  bool inputErrorFlag = 0;
  // Inputs the Array1D from MatrixMarket format. Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  if(a.isValid()) {
    // Skip all comment lines (which should all be at the top of the output file
    // All comment lines should begin with a '%'.  Comments cannot be on an
    // existing data line even if preceded by a '%' after the data.
    // 1024 is technically the maximum length in MM format, but we ignore a line
    // much larger than that just in case.

    // First, check first format line.  It must conform to MatrixMarket
    // specification by.  We only import the type.
    // "%%MatrixMarket matrix array real general"
    fin >> sInput; if(sInput!="%%MatrixMarket") { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="matrix")         { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="array")          { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="real")           { inputErrorFlag=1; }
    std::getline(fin,sInput,'\n'); if(sInput!="general") { inputErrorFlag=1; }
    if(inputErrorFlag!=0) {
      cerr << "Array1D MatrixMarket header incorrect.\n"; return fin; }
    // Now proceed to ignore all comment lines after the header.
    char checkCommentLine = fin.peek();
    while (checkCommentLine == '%') {
      if (checkCommentLine == '%') {
        fin.ignore(MaxLineLength,'\n');
        nComments++;
      } // end taking care of a comment line
      checkCommentLine = fin.peek();  // get next line's first char
    } // end comment line skipping
    // Import all actual data starting with the number of rows and columns.
    fin >> N >> M; // discard M here since it must be 1
    if(M!=1) {
      cerr << "Error importing Array1D. Number of columns > 1.\n"; return fin; }
    a.~Array1D<T>(); // explicit destructor call since it is being redefined
    a.L = N;
    T* paElement = &(a[0]);
    for(unsigned i=0; i<a.L; i++) { fin >> (*paElement);  paElement++; }
    if(a.VerboseErrorReports!=0) {
      cout << "Array Length = " << N
           << "Number of comment lines " << nComments << ".\n";
    } // end if VerboseErrorReports
  } // end if
  return fin;
};

// Now, code the commuted x*vector operation as a NON-MEMBER function.
template <typename T>
Array1D<T> operator*(T x, const Array1D<T>& b) {
  Array1D<T> c(b.getSize());
  if(b.isValid()) {
    T* pbElement = &(b[0]);
    T* pcElement = &(c[0]);
    for(unsigned i=0; i<b.getSize(); i++) {
      *pcElement = (*pbElement)*x;  pbElement++;  pcElement++; }
  } // end if isValid
  else cerr << "Invalid non-method scalar '*' operation for Array1D object.\n";
  return c;
};


template <typename T>
ostream& TVectorOut(ostream& fo, Array1D<T> &v, 
             bool bHiLiteMax = 0, bool bHiLiteMin = 0, bool bStopOnNegative = 0,
             const char colorCode[] = normText, unsigned length=INT_MAX, 
             string label = "", int W1=0, int W2=18, char commentChar = '\0') {
  // output the vector horizontally rather than in MatrixMarket format
  // There is a bug with the streaming return operator.
  // The return ofstream& is having problems.

  // find min and max values
  T  maxVal = v[0], minVal = maxVal;  // will have problems for T double
  T  vi;                              // temporary value
  int  maxIndex = -1, minIndex = -1, minLength = min(v.getSize(),length);
  int i;
  bool bFinished = 0;
  if(bHiLiteMax || bHiLiteMin) {
    maxIndex = 0;     minIndex = 0;
    i = 1;
    while(i<minLength && !bFinished) {
      // if bStopOnNegative, we only track for first positive values in array
      if( !(bStopOnNegative && v[i]<0) ) {
        vi = v[i];
        if(vi > maxVal)      { maxVal = vi;  maxIndex = i; }
        else if(vi < minVal) { minVal = vi;  minIndex = i; }
      } else bFinished = 1;
      
      i++;
    } // end while i
  } // end if bHiLite...

  // now output results
  //fo << commentChar << setw(W1) << colorCode << label.c_str();
  fo << colorCode;
  i = 0;
  bFinished = 0;
  while(i<minLength && !bFinished) {
    if( !(bStopOnNegative && v[i]<0) ) {
      //if(bHiLiteMax && i==maxIndex)       cout << green;
      //else if(bHiLiteMin && i==minIndex)  cout << brown;
      if(i==maxIndex)       cout << green;
      else if(i==minIndex)  cout << brown;
      //fo << setw(W2) << v[i];
      fo << v[i] << " ";
      if((bHiLiteMax && i==maxIndex) || (bHiLiteMin && i==minIndex)) 
        cout << colorCode;
    } else bFinished = 1; 

    i++;
  } // end while i
  //fo << endl;
  
  return fo;  // streaming usage not working, not sure how to implement it?
} // end TVectorOut


template <typename T>
Array1D<T>& Array1D<T>::operator=(TList<T>& b) { 
  // assignment for copying a list to an Array1D<T>
  // check for valid arrays first
  #ifdef DEBUG_MODE
  if(!isValid(1)) errorMsg("Invalid Array1D<T> object in operator=(TList<T>&)");
  else if(!b.isValid())
        errorMsg("Invalid TList<T> object in Array1D<T>::operator=(TList<T>&)");
  #endif
  if(L==b.getSize() && b.getSize()>0) {
    T* pElement = data;
    L = b.getSize();
    // use TList overloaded op[] to copy list data
    b.begin();  // begin manual iteration
    for(unsigned i=0; i<L; i++) { 
      *pElement = b.getCurrent();  pElement++;  b.next(); }
    //for(unsigned i=0; i<L; i++) { *pElement = b.get(i);  pElement++;  }
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  }  // end if
  else if(b.getSize()>0) {
    // this version of the assignment operator allows size of *this to change
    #ifdef DEBUG_MODE
    if(L>0)  debugMsg("Array1D<T> object changed size in op=(TList<T>&) from "
                      +itos(L));
    #endif
    if(data!=NULL)  delete[] data;
    L = b.getSize();  LMax = b.getSize();
    data = new T[L];  // allocate new size of data
    T* pElement = data;
    b.begin();  // begin manual iteration
    for(unsigned i=0; i<L; i++) { 
      *pElement = b.getCurrent();  pElement++;  b.next(); }
    //for(unsigned i=0; i<L; i++) { *pElement = b.get(i);  pElement++;  }
    #ifdef DEBUG_MODE
    ValidArray = Valid1DArrayFlag;
    #endif
  } // end mutable assignment
  else {
    if(data!=NULL)  delete[] data;
    L = 0;  LMax = 0;  data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    warningMsg("Array1D<T> operator=(TList&) is 'copying' a size 0 TList<T>");
    #endif
  } // end else for L==0
  return *this;
}; // end =


} // end namespace MSL
#endif
