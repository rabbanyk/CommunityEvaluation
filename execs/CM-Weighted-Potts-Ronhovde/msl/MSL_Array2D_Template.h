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
MSL - My (or Math) Scientific Library - Array2D<T>
See MSL_Array.h for general comments and some documentation.

Array2D<T> Specific Notes:
These notes handle only behavior specific to the Array2D<T> implementation.

1. If StrictErrorChecking != 1:  removed functionality 5/22/08 - 
Define DEBUG_MODE as a pre-processor parameter instead.  This method only 
inserts the additional error checking code if actually debugging rather than 
trying to rely on a 'smart' compiler to eliminate the dead code (except for 
dead if statements).

2. operator[] is defined only for the T* return value since once we have a T*, 
it is a intrinsic type.  Operator[] must be defined as a member function which 
we cannot do for double* type.  Therefore, we are forced to let the compiler 
handle out of bounds and error checking after this point.
However, full error checking is done at the T* operator[] stage here; therefore, 
at least up to there, it is assured that we have a valid object.
An alternate implementation would be to create an array of Array1D<T> objects 
which might prove useful if I implemented a full TMatrix class.

Known problems:
1. For some reason operator<< is trying to call the pure virtual isValid() 
function; therefore, the error check is commented out for now.
*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TARRAY2D_H
#define MSL_TARRAY2D_H

//#define DEBUG_MODE

#include <string>

#ifndef MSL_TARRAY_H
#include "MSL_Array.h"
#endif

namespace MSL {

template <typename T> class  Array2D;
template <typename T> Array2D<T>  operator*(T x, const Array2D<T>& a);
template <typename T> istream&    operator>>(istream& fin, Array2D<T>& a);

#ifndef MSL_TMATRIX_T
// define some common vector types
typedef Array2D<TFloat> TMatrixFloat;  // TFloat  array type definition
typedef Array2D<int>    TMatrixInt;    // integer array type definition
//typedef Array2D<bool> TMatrixBool;   // integer array type definition
#define MSL_TMATRIX_T
#endif

//------------------------------------------------------------------------------
template <typename T> 
class Array2D : public TMatrix {
 // Memory allocation is dynamic, and this Array2D<T> uses the standard C++ array
 // indexing with brackets as in arrayName[col][row] with zero indexing.
 //friend class Sparse2D<T>;
 friend Array2D<T>  operator*<T>(T x, const Array2D<T>& a);
 friend istream&    operator>><T>(istream& fin, Array2D<T>& a);
 //friend class Array3D;
 public:
  // overloaded operators
         Array2D<T>& operator= (const Array2D<T>& b); // = will resize if needed
         Array2D<T>& operator+=(const Array2D<T>& b); // arrays must be same
         Array2D<T>& operator-=(const Array2D<T>& b); //   size for *= and /=
         Array2D<T>& operator*=(T x);                 // *= by scalar
         Array2D<T>& operator/=(T x);                 // /= by scalar
         // standard C array reference error checking only on column reference
  inline T*          operator[](unsigned col) const;  // standard C array ref
  inline T*          operator[](unsigned col);
         // alternate reference operator that allows full error checking
  inline T&          operator()(unsigned col, unsigned row) const ;
  inline T&          operator()(unsigned col, unsigned row);
         // the following operations return by value to they are inefficient
         Array2D<T>  operator-();                    // unary negative
         Array2D<T>  operator-(const Array2D<T>& a); // binary +
         Array2D<T>  operator+(const Array2D<T>& a); // binary -
         Array2D<T>  operator*(T x);                 // binary * by scalar
         Array2D<T>  operator/(T x);                 // binary / by scalar
  // Public member functions
  inline int         getCols() const { return Cols; };
  inline int         getRows() const { return Rows; };
  inline T           getMax(bool bSymm, bool bAbsoluteValue=0) const;
  inline T           getMin(bool bSymm, bool bAbsoluteValue=0) const;
         // if one also wants the matrix indices for the max or min value...
         T           getMax(unsigned &col, unsigned &row, bool bSymm=0, bool bAbsVal=0) const;
         T           getMin(unsigned &col, unsigned &row, bool bSymm=0, bool bAbsVal=0) const;
         // utility functions
         void        resize(unsigned Cols, unsigned Rows, bool bCopyData=0);
         void        init(T Value);                  // reinitialize
         void        initDiagonal(T Value);          // reinitialize
         void        erase();                        // return dynamic memory
         // eliminate an array column and row by summing the column and row data
         // Requires a square matrix at the moment. Returns an N-1 x N-1 matrix.
         // Stores the collapsed data in the minimum index of ColOne or ColTwo.
         void        collapse(unsigned ColOne, unsigned ColTwo);
  // for debugging... is the current object a valid one? Methods perform regular
  // checks if DEBUG_MODE is defined as a pre-processor parameter, otherwise
  // it simply returns a boolean true.
  inline virtual bool isValid(bool bEmptyOK=0) const;       // valid array check
  // standard member functions
          Array2D(unsigned Cols = 0, unsigned Rows = 0);  // default constructor
          Array2D(unsigned Cols, unsigned Rows, T Value); // constructor w/init
          Array2D(const Array2D<T> &b);                   // copy constructor
  virtual ~Array2D();                                     // destructor
 protected:
  // protected data elements
  unsigned  Cols, Rows;       // number of columns and rows in the array
  T**       data;             // actual data array
}; // end class Array2D<T>;

//----------------------Array2D<T> Member Function Declarations--------------------

//Assign default static parameters
//Array2D<T>::StrictErrorChecking = 1;

template <typename T> 
Array2D<T>::Array2D(unsigned Cs, unsigned Rs) {
  // Raw constructor that skips the initialization.  Of course, it takes the 
  // chance that the user will access junk, but it is provided for "more 
  // efficient" applications where the user wants to initialize himself.
  // It is also the default constructor, it generates a Array2D<T> object that 
  // is a zero length real Array2D<T> object but officially an invalid one since
  // it holds no data.
  Cols = Cs;  Rows = Rs;
  if(Cols>0 && Rows>0) {
    data = new T*[Cols];                                // set up rows
    for(unsigned i=0; i<Cols; i++)  data[i] = new T[Rows];  // set up columns
    #ifdef DEBUG_MODE
    ValidArray = Valid2DArrayFlag;
    #endif
  } // end un-initialized array creation
  else if(Cols==0 || Rows==0) {
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
    data = NULL;
    if(VerboseErrorReports!=0)
      errorMsg("Error creating a zero length Array2D<T> object!  Rows x Cols = "
               +itos(Rows)+" x "+itos(Cols));
  } // end zero length if
}; // end Array2D<T> constructor


template <typename T> 
Array2D<T>::Array2D(unsigned Cs, unsigned Rs, T Value) { // constructor with init value
  Cols = Cs;  Rows = Rs;
  if(Cols>0 && Rows>0) {
    data = new T*[Cols];
    for(unsigned j=0; j<Cols; j++) { data[j] = new T[Rows]; } // Set up columns
    // Now, initialize array
    T  *pa;
    for(unsigned j=0; j<Cols; j++) { 
      pa = &(data[j][0]);  (*pa) = Value;
      for(unsigned i=1; i<Rows; i++) { pa++;  (*pa) = Value; }
    } // end for j
    #ifdef DEBUG_MODE
    ValidArray = Valid2DArrayFlag;
    #endif
  } // end un-initialized array creation
  else if(Cols==0 && Rows==0) {
    data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    if(VerboseErrorReports!=0)
      errorMsg("Error creating zero length Array2D<T> object!  Rows x Cols = "
               +itos(Rows)+" x "+itos(Cols));
    #endif
  } // end zero length if
}; // end Array2D<T> constructor


template <typename T> 
Array2D<T>::Array2D(const Array2D<T> &b) { // copy constructor
  #ifdef DEBUG_MODE
  if( !(isValid() && b.isValid()) )  
    errorMsg("Invalid Array2D<T> object in copy constructor");
  #endif
  Cols = b.Cols;  Rows = b.Rows;
  if(b.Cols>0 && b.Rows>0) { 
    data = new T*[Cols];                                // set up rows
    for (int j=0; j<Cols; j++)  data[j] = new T[Rows];  // set up columns
    // Now, initialize array
    T  *pa, *pb;
    for(unsigned j=0; j<Cols; j++) { 
      pa = &(data[j][0]);  pb = &(b.data[j][0]);
      (*pa) = (*pb);
      for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  (*pa) = (*pb); }
    } // end for j
    #ifdef DEBUG_MODE
    ValidArray = b.ValidArray;
    #endif
  } // end if
  else {
    data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    warningMsg("Copying a size zero Array2D<T> object?");
    #endif
  } // end else
}; // end Array2D<T> copy constructor


template <typename T> 
inline Array2D<T>::~Array2D() {
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Entering Array2D<T> destructor on a size "+itos(Rows)+" x "
           +itos(Cols)+" object... ");
  #endif
  #ifdef DEBUG_MODE
  if(!isValid(1))  warningMsg("Destructed Array2D<T> object is not valid?");
  #endif
  erase();
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting Array2D<T> destructor\n");
  #endif
}; // end Array2D<T> destructor


template <typename T> 
void Array2D<T>::erase() {
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("Invalid Array2D<T> object in erase()");
  #endif
  for(unsigned j=0; j<Cols; j++)  delete[] data[j];  // Delete columns
  delete[] data;  data = NULL;
  Rows = -1;  Cols = -1;  // set to invalid value just in case
  #ifdef DEBUG_MODE
  ValidArray = 0;
  #endif
  return;
}; // end erase


template <typename T> 
inline bool Array2D<T>::isValid(bool bEmptyOK) const { 
  #ifdef DEBUG_MODE
  //if(StrictErrorChecking==1) { // deprecated by me
  bool bResult = 0;
  if(bEmptyOK && ValidArray==0 && Cols==0 && Rows==0 && data==NULL)  return 1;
  else if(ValidArray==0)  warningMsg("Array2D<T> object is invalid-empty?");
  else if(ValidArray!=Valid2DArrayFlag) 
    warningMsg("Invalid Array2D<T> object signature field");
  // perform a consistency check on the 'unsigned' indices in case we point to
  // a random memory address
  // perform a consistency check data pointer data member
  else if(Rows==0 && Cols==0 && data!=NULL) 
    warningMsg("Array2D<T> is empty by data pointer is not NULL?");
  else if(!(Rows==0 && Cols==0) && data==NULL) 
    warningMsg("Array2D<T> data is NULL but sizes are not zero?");
  // appears to be a valid object
  else   bResult = 1;  // if we make it to here, it is ok as best we can tell
  return bResult;
  //} // if StrictErrorChecking
  #else
  return 1;
  #endif
}; // end Array2D<T>::isValid


template <typename T> 
Array2D<T>& Array2D<T>::operator=(const Array2D<T>& b) {
  #ifdef DEBUG_MODE
  if( !(isValid(1) && b.isValid()) )  
    errorMsg("Invalid Array2D<T> object in assignment operator=()");
  #endif
  if( !(Cols==b.Cols || Rows==b.Rows) ) {
    #ifdef DEBUG_MODE
    cout << red << "Resizing Array2D<T> object for '=' operation... " << flush;
    #endif
    erase();
    resize(b.Cols,b.Rows,0);
    #ifdef DEBUG_MODE
    cout << "done with resize.\n" << normText << flush;
    #endif
  } // end if resize check
  T  *pa, *pb;
  for(unsigned j=0; j<Cols; j++) { 
    pa = &(data[j][0]);  pb = &(b.data[j][0]);
    (*pa) = (*pb);
    for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  (*pa) = (*pb); }
  } // end for j
  #ifdef DEBUG_MODE
  ValidArray = b.ValidArray;
  #endif
  return *this;
}; // end =


template <typename T> 
Array2D<T>& Array2D<T>::operator+=(const Array2D<T>& b) {   // additive assignment
  #ifdef DEBUG_MODE
  if( !(isValid() && b.isValid()) ) 
    errorMsg("Invalid Array2D<T> object in operator+=()");
  #endif
  if(Cols==b.Cols && Rows==b.Rows) {
    T  *pa, *pb;
    for(unsigned j=0; j<Cols; j++) { 
      pa = &(data[j][0]);  pb = &(b.data[j][0]);
      (*pa) += (*pb);
      for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  (*pa) += (*pb); }
    } // end for j
  } // end if
  else  errorMsg("Array2D<T> objects are difference sizes in operator+=()");
  return *this;
}; // end +=


template <typename T> 
Array2D<T>& Array2D<T>::operator-=(const Array2D<T>& b) {   // subtraction assignment
  #ifdef DEBUG_MODE
  if( !(isValid() && b.isValid()) ) 
    errorMsg("Invalid Array2D<T> object in operator+()");
  #endif
  if(Cols==b.Cols && Rows==b.Rows) {
    T  *pa, *pb;
    for(unsigned j=0; j<Cols; j++) {
      pa = &(data[j][0]);  pb = &(b.data[j][0]);
      (*pa) -= (*pb);
      for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  (*pa) -= (*pb); }
    } // end for j
  } // end if
  else  errorMsg("Array2D<T> objects are difference sizes in operator-=()");
  return *this;
}; // end -=


template <typename T> 
Array2D<T>& Array2D<T>::operator*=(T x) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in scalar operator*=()");
  #endif
  T  *pa;
  for(unsigned j=0; j<Cols; j++) { 
    pa = &(data[j][0]);  (*pa) *= x;
    for(unsigned i=1; i<Rows; i++) { pa++;  (*pa) *= x; }
  } // end for j
  return *this;
}; // end *=


template <typename T> 
Array2D<T>& Array2D<T>::operator/=(T x) {   // T division assignment
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in scalar operator/=()");
  #endif
  if(x!=0) {
    T  *pa;
    for(unsigned j=0; j<Cols; j++) {
      pa = &(data[j][0]);  (*pa) /= x;
      for(unsigned i=1; i<Rows; i++) { pa++;  (*pa) /= x; }
    } // end for j
  } else  errorMsg("Array2D<T> division by zero error in operator/=()");
  return *this;
}; // end /=


template <typename T> 
inline T* Array2D<T>::operator[](unsigned i) const {
#ifdef DEBUG_MODE
  // check for valid array first
  if(!isValid())    errorMsg("Array2D<T> object is invalid in op[i] const reference");
  else if(i>=Cols)  errorMsg("Array2D<T> index out of bounds");
  else              return data[i];  // column reference appears valid
#else
  return data[i];
#endif
}; // end array reference const
template <typename T> 
inline T* Array2D<T>::operator[](unsigned i) {
  #ifdef DEBUG_MODE
  // check for valid array first
  if(!isValid())    errorMsg("Array2D<T> object is invalid in op[i] reference");
  else if(i>=Cols)  errorMsg("Array2D<T> index out of bounds");
  else              return data[i];  // column reference appears valid
  #else
  return data[i];
  #endif
}; // end array reference
// include an alternate reference operator that allows for full error checking
// on both i and j indices
template <typename T> 
inline T& Array2D<T>::operator()(unsigned i, unsigned j) const {
  #ifdef DEBUG_MODE
  // check for valid array first
  if(!isValid())  errorMsg("Array2D<T> object is invalid in (i,j) const reference");
  else if(i>=Cols)  errorMsg("Array2D<T> col index out of bounds");
  else if(j>=Rows)  errorMsg("Array2D<T> row index out of bounds");
  else return data[i][j];  // finally return a valid element within array bounds
  #else
  return data[i][j];
  #endif
}; // end array reference const
template <typename T> 
inline T& Array2D<T>::operator()(unsigned i, unsigned j) {
  #ifdef DEBUG_MODE
  // check for valid array first
  if(!isValid())    errorMsg("Array2D<T> object is invalid in (i,j) reference");
  else if(i>=Cols)  errorMsg("Array2D<T> col index out of bounds");
  else if(j>=Rows)  errorMsg("Array2D<T> row index out of bounds");
  else return data[i][j];  // finally return a valid element within array bounds
  #else
  return data[i][j];
  #endif
}; // end array reference


template <typename T> 
Array2D<T>  Array2D<T>::operator-() {             // unary negative
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in unary operator-()");
  #endif
  Array2D<T> c(Cols,Rows);
  T  *pa;
  for(unsigned j=0; j<Cols; j++) { 
    pa = &(c.data[j][0]);  (*pa) *= -1;
    for(unsigned i=1; i<Rows; i++) { pa++;  (*pa) *= -1; }
  } // end for j
  return c;
}; // end unary negative


template <typename T> 
Array2D<T> Array2D<T>::operator+(const Array2D<T>& b) { // array addition
  #ifdef DEBUG_MODE
  if( !(isValid() && b.isValid()) )  
    errorMsg("Invalid Array2D<T> object in binary operator+()");
  #endif
  Array2D<T> c(Cols,Rows);
  if(Rows==b.Rows && Cols==b.Cols) {
    T  *pa, *pb, *pc;
    for(unsigned j=0; j<Cols; j++) { 
      pa = &(data[j][0]);  pb = &(b.data[j][0]);  pc = &(c.data[j][0]);
      (*pc) = (*pa)+(*pb);
      for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  pc++;  (*pc) = (*pa)+(*pb); }
    } // end for j
  } // end if error check
  else  errorMsg("Array2D<T> objects are difference sizes in operator+()");
  return c;
}; // end +


template <typename T> 
Array2D<T> Array2D<T>::operator-(const Array2D<T>& b) { // array subtraction
  #ifdef DEBUG_MODE
  if( !(isValid() && b.isValid()) )  
    errorMsg("Invalid Array2D<T> object in binary operator-()");
  #endif
  Array2D<T> c(Cols,Rows);
  if(Rows==b.Rows && Cols==b.Cols) {
    T  *pa, *pb, *pc;
    for(unsigned j=0; j<Cols; j++) {
      pa = &(data[j][0]);  pb = &(b.data[j][0]);  pc = &(c.data[j][0]);
      (*pc) = (*pa)-(*pb);
      for(unsigned i=1; i<Rows; i++) { pa++;  pb++;  pc++;  (*pc) = (*pa)-(*pb); }
    } // end for j
  } // end if
  else  errorMsg("Array2D<T> objects are difference sizes in operator+()");
  return c;
}; // end binary -


template <typename T> 
Array2D<T> Array2D<T>::operator*(T x) {   // vector scaling vec*x
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in scalar operator*()");
  #endif
  Array2D<T> c(Cols,Rows);
  T  *pa, *pc;
  for(unsigned j=0; j<Cols; j++) { 
    pa = &(data[j][0]);  pc = &(c.data[j][0]);;
    (*pc) = (*pa)*x;
    for(unsigned i=1; i<Rows; i++) { pa++;  pc++;  (*pc) = (*pa)*x; }
  } // end for j
  return c;
}; // end scalar multiplication


template <typename T> 
Array2D<T> Array2D<T>::operator/(T x) {   // scalar division
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in scalar operator/()");
  #endif
  Array2D<T> c(Cols,Rows);
  if(x!=0) {
    T  *pa, *pc;
    for(unsigned j=0; j<Cols; j++) { 
      pa = &(data[j][0]);  pc = &(c.data[j][0]);
      (*pc) = (*pa)/x;
      for(unsigned i=1; i<Rows; i++) { pa++;  pc++;  (*pc) = (*pa)/x; }
    } // end for j
  } // end if
  else  errorMsg("Division by zero error on Array2D object.\n");
  return c;
}; // end scalar division
// x/Array2D<T> is not a valid operation, so it is not defined


// Other useful array functions
template <typename T> 
void  Array2D<T>::collapse(unsigned c1, unsigned c2) {
  // this is an expensive operation that collapses two given rows/columns
  // and returns the reduced size matrix by reference
  // at the moment it is implemented only on square matrices
  // not currently optimized to allow testing changes more easily
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in collapse()");
  #endif
  Array2D<T>  b;
  if(Cols!=Rows) {
    errorMsg("Array2D collapse(...) is only allowed on square matrices");
    //return *this;
    return;
  } // end error test
  else if(c1==c2 || c1>Rows || c2>Rows || c1>Cols || c2>Cols) {
    errorMsg("Array2D collapse parameter error: c1="+itos(c1)+"c2="+itos(c2));
    //return *this;
    return;
  } // end error test
  else {
    unsigned  minC, maxC, i, j;  
    if(c1>c2) { minC = c2;  maxC = c1; }
    else      { minC = c1;  maxC = c2; }
    // now create the collapsed matrix
    b.resize(Cols-1,Rows-1);
    // copy before maxC
    for(i=0; i<maxC; i++) {
      for(j=0; j<maxC; j++)       b[i][j]   = data[i][j];
      for(j=maxC+1; j<Rows; j++)  b[i][j-1] = data[i][j];
      b[minC][i] += data[maxC][i];           // merge the row data
      b[i][minC] += data[i][maxC];           // merge the column data
    } // end for i
    // copy after maxC
    for(i=maxC+1; i<Cols; i++) {
      for(j=0; j<maxC; j++)       b[i-1][j]   = data[i][j];
      for(j=maxC+1; j<Rows; j++)  b[i-1][j-1] = data[i][j];
      b[minC][i-1] += data[maxC][i];         // merge the row data
      b[i-1][minC] += data[i][maxC];         // merge the column data
    } // end for i
    b[minC][minC]  += data[maxC][maxC];      // merge the omitted element
  } // end if Cols, Rows
  //return b;
  // now copy the data back to current one
  (*this) = b;
  return;
}; // end collapse


template <typename T> 
void Array2D<T>::resize(unsigned Cs, unsigned Rs, bool bCopyData) {
  // Crude resize function to "get things working" at the moment
  // At the moment it DELETES ALL DATA!!!
  // first delete existing array data
  #ifdef DEBUG_MODE
  if(!isValid(1)) { // check whether this object is otherwise valid
    if( !(ValidArray==0 && data==NULL) ) 
      errorMsg("Array2D resize() error.  Invalid non-empty object in resize.");
    // else // OK to resize an empty object, do nothing
  } // end isValid check
  #endif
  if(bCopyData) {
    // declare and fill a copy to minimum values and set data pointer to copy
    // the rest of the array could be uninitialized
    T **dataCopy;  dataCopy = new T*[Cs];
    T *pa, *pb;  // optimization pointers
    unsigned minR = min(Rs,Rows), minC = min(Cs,Cols);
    for(unsigned i=0; i<Cs; i++)  dataCopy[i] = new T[Rs];
    // now copy data over to dataCopy
    for(unsigned i=0; i<minC; i++) {
      pa = &(data[i][0]);  pb = &(dataCopy[i][0]);
      (*pa) = (*pb);
      for(unsigned j=0; j<minR; j++) { pa++;  pb++;  (*pa) = (*pb); }
    } // end for i
    erase();          // erase all data in *this object
    data = dataCopy;  // set data pointer to copy
    Cols = Cs;  Rows = Rs;
    #ifdef DEBUG_MODE
    ValidArray = Valid2DArrayFlag;
    #endif
  } else {
    // return memory and reallocate uninitialized new Array2<T>
    erase();
    // now create the new image
    data = new T*[Cs];
    for(unsigned i=0; i<Cs; i++)  data[i] = new T[Rs];
    Cols = Cs;  Rows = Rs;
    #ifdef DEBUG_MODE
    ValidArray = Valid2DArrayFlag;
    #endif
  } // end if bCopyData
  return;
}; // end resize


template <typename T> 
void Array2D<T>::init(T Value) {
  // reinitialize the Array
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in init()");
  #endif
  T  *pa;
  for(unsigned j=0; j<Cols; j++) {
    pa = &(data[j][0]);  (*pa) = Value;
    for(unsigned i=1; i<Rows; i++) { pa++;  (*pa) = Value; }
  } // end for j
  #ifdef DEBUG_MODE
  ValidArray = Valid2DArrayFlag;
  #endif
  return;
}; // end init


template <typename T> 
void Array2D<T>::initDiagonal(T Value) {
  // reinitialize the array diagonal
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in init()");
  #endif
  T  *pa;
  int iMax = min(Cols,Rows);
  for(unsigned i=0; i<iMax; i++)  data[i][i] = Value;
  #ifdef DEBUG_MODE
  ValidArray = Valid2DArrayFlag;
  #endif
  return;
}; // end init


template <typename T> 
T Array2D<T>::getMin(unsigned &col, unsigned &row, 
                     bool bSymmetric, bool bAbsoluteValue) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  // col and row are the first index location where the minimum is located
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in getMin()");
  #endif
  T  *pa = &(data[0][0]), minVal;
  unsigned jEnd;  // for symmetric optimization
  if(bAbsoluteValue) {
    minVal = abs(*pa);
    for(unsigned i=0; i<Cols; i++) { 
      pa = &(data[i][0]);
      if(bSymmetric) jEnd = i;  else jEnd = Rows;
      for(unsigned j=0; j<jEnd; j++) {
        if(abs(*pa)<minVal) { minVal = abs(*pa);  col = i;  row = j; }
        pa++; // must increment pa regardless of minimum result
      } // end for j
    } // end for i
  } else {
    minVal = (*pa);
    for(unsigned i=0; i<Cols; i++) { 
      pa = &(data[i][0]);
      if(bSymmetric) jEnd = i;  else jEnd = Rows;
      for(unsigned j=0; j<jEnd; j++) {
        if((*pa)<minVal) { minVal = (*pa);  col = i;  row = j; }
        pa++; // must increment pa regardless of minimum result
      } // end for j
    } // end for i
  } // end 
  return minVal;
}; // end getMin with row,col
template <typename T> 
inline T  Array2D<T>::getMin(bool bSymmetric, bool bAbsoluteValue) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  unsigned iDummy, jDummy;
  return getMin(iDummy,jDummy,bSymmetric,bAbsoluteValue);
}; // end getMin


template <typename T> 
T Array2D<T>::getMax(unsigned &col, unsigned &row, 
                     bool bSymmetric, bool bAbsoluteValue) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  // col and row are the first index location where the maximum is located
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("Invalid Array2D<T> object in getMax()");
  #endif
  T   *pa, maxVal;
  int  jStart;
  if(bAbsoluteValue) {
    maxVal = abs(*pa);
    for(unsigned i=0; i<Cols; i++) { 
      pa = &(data[i][0]);
      if(bSymmetric)  jStart = i+1;
      else            jStart = 0;
      for(unsigned j=jStart; j<Rows; j++)
        if(abs(*pa)>maxVal) { maxVal = abs(*pa);  col = i;  row = j;  pa++; }
    } // end for i
  } else {
    maxVal = (*pa);
    for(unsigned i=0; i<Cols; i++) { 
      pa = &(data[i][0]);
      if(bSymmetric)  jStart = i+1;
      else            jStart = 0;
      for(unsigned j=jStart; j<Rows; j++)
        if((*pa)>maxVal) { maxVal = (*pa);  col = i;  row = j;  pa++; }
    } // end for i
  } // end 
  return maxVal;
}; // end getMax with row,col
template <typename T> 
inline T  Array2D<T>::getMax(bool bSymmetric, bool bAbsoluteValue) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  unsigned iDummy, jDummy;
  return getMax(iDummy,jDummy,bSymmetric,bAbsoluteValue);
}; // end getMax
//------------------------------------------------------------------------------
// End Array2D<T> method definitions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// the commuted x*matrix operation must be a non-member function
template <typename T> 
inline Array2D<T> operator*(T x, const Array2D<T>& a) { return (a*x); };

// Overloaded iostream operators
template <typename T> 
ostream& operator<<(ostream& fout, Array2D<T>& a) {
  // Outputs the Array1D in MatrixMarket format.  Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  
  // There is a bug where the insertion operator << is calling the pure virtual 
  // a.isValid(). The extraction op >> does not have this problem though.
  //if(a.isValid()) { 
    //fout << MMHeader;  // disabled for general << - write a different function
    // must somehow allow comment lines here at some point
    fout << a.getRows() << " " << a.getCols() << endl;
    for(int i=0; i<a.getRows(); i++) {
      for(int j=0; j<a.getCols(); j++)  fout << a[j][i] << "\t"; 
      cout << "\n";  // line break
    } // end for j
  //} // end if
  return fout;
}; // end stream output


template <typename T> 
istream& operator>>(istream& fin, Array2D<T>& a) {
  int M, N;  // cols and rows - here M = 1 since it is an Array1D
  int nComments=0;
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
      cerr << "Array2D<T> MatrixMarket header incorrect.\n"; return fin; }
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
    //cout << "Before destructor in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    a.erase();
    //cout << "After  destructor in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    a.Rows = N;  a.Cols = M;
    a.data = new T*[a.Cols];
    for (int i=0; i<a.Cols; i++) { a.data[i] = new T[a.Rows]; }
    a.ValidArray = Valid2DArrayFlag;
    for(unsigned j=0; j<a.getCols(); j++) { for(unsigned i=0; i<a.getRows(); i++) {
      fin >> a.data[j][i]; }} // end ij
    //cout << "After redefinition in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    if(a.VerboseErrorReports!=0) {
      cout << "Array Rows x Cols = " << N << " " << M << "  "
           << "Number of comment lines " << nComments << ".\n";
    } // end if VerboseErrorReports
  } // end if
  return fin;
}; // end stream input

} // end namespace MSL
#endif

