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
MSL - My (or Math) Scientific Library - TCMSparse
See MSL_Array.h for general comments and some documentation.

TCMSparse Specific Notes:
These notes handle only behavior specific to the TCMSparse implementation.

1. If StrictErrorChecking != 1 then if the objects are of different sizes,
the assignment operator destroys the original object and recreates one with the
new 'correct' size to make the assignment.

4. Operator[] is defined only for the T* return value since once we have a
T* (usually a double*), it is a intrinsic type.  Operator[] must be defined
as a member function which we cannot do for double* type.  Therefore, we are
forced to let the compiler handle out of bounds and error checking after this
point.  However, full error checking is done at the pointer return stage;
therefore, from there at least, it is assured that we have valid objects after
this point.
If needed, an alternate implementation at some point would be to create an
array of Array1D objects.  At this point though, I do not need any of the 
added functionality of such an implementation.  It might prove useful if I 
implemented a full TMatrix class.
*/
/* Known problems:
1. For some reason operator<< is trying to call the pure virtual isValid() function;
therefore, the error check is commented out for now.
*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TCMSPARSE_UNWEIGHTED_H
#define MSL_TCMSPARSE_UNWEIGHTED_H

//#define DEBUG_MODE

#include <iostream>
#include <iomanip>
#include <string>

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
//#ifndef MSL_TLIST_H
#include "MSL_List_Template.h"
//#endif
#ifndef MSL_TLISTARRAY_H
#include "MSL_ListArray_Template.h"
#endif
#ifndef MSL_TCMATRIX_H
#include "MSL_CMatrix.h"
#endif

namespace MSL {

class TCluster;
class TClusterList;
template <typename T> class Array1D;

// prototype some classes
//template <typename T> class TCMSparse;
//template <typename T> class  TList;
//template <typename T> class  TListArray;
// prototype a few functions from other classes???
//template <typename T> inline int TListArray<T>::getSize() const;
//template <typename T> TCMSparse  operator*(T x, const TCMSparse& a);
//template <typename T> istream&     operator>>(istream& fin, TCMSparse& a);

#ifndef MSL_TSPARSE_T
#define MSL_TSPARSE_T
// define some common vector types
//typedef TCMSparse<TFloat> TSparseFloat;// TFloat  array type definition
//typedef TCMSparse<int>    TCMSparseInt;     // integer array type definition
#endif


//------------------------------------------------------------------------------
//template <typename T> 
class TCMSparse : public TCMatrix {
 // Memory allocation is dynamic
 //friend istream&    operator>><T>(istream& fin, TCMSparse& a);
 public:
          // overloaded operators
          TCMSparse& operator=(const TCMSparse& b);
          // a = b for TListArray defines a square matrix of size b.getSize()
          // using default values of storedValue=1 and unstoredValues=0
          // use set(...) for a non-square matrix
          TCMSparse& operator=(const TListArray<int>& b);
  inline  TCMSparse& operator*=(TCMData x);          // *= by scalar
  inline  TCMSparse& operator/=(TCMData x);          // /= by scalar
         // operator() is a constant reference - cannot change elements
         // does not use op[][] due to structure constraints for binary search
  //          TCMData  operator()(unsigned col, unsigned row) const;
  //          TCMData  operator()(unsigned col, unsigned row);  // error here
  virtual TCMData  operator()(unsigned col, unsigned row) const;
  virtual TCMData  operator()(unsigned col, unsigned row);  // error here
         // set() sets the specified element to the storedValue adding the 
         // element location if it does not exist (expensive operation)
         //bool           set(unsigned col, unsigned row);
         // copy a set of indices to the matrix, that is set a whole kij list
  //inline void           set(unsigned col, Array1D<int> &a);
  //inline void           set(unsigned col, TList<int> &a);
         // returns true if unset or false if element was not being stored
         //bool           unset(unsigned col, unsigned row);
  // other functions
  virtual inline void   scale(TCMData a, TCMData b, TCMData aOld=1, TCMData bOld=1);
         // initialize the sparse matrix in various ways
                 void   initByClusters(TClusterList &c, double pin, double pout);
                 void   initHeterogenous(TClusterList &c, Array1D<int> &v,
                                         double p2, double p1, double p0);
  void                  initHierarchy(unsigned N, unsigned n2, unsigned n1, unsigned n0,
                                      double p3,  double p2,   double p1,   double p0,
                                      bool bRandomizeNodes=1);
         // initialize a 'staggered' hierarchy
  void                  initSHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
                                      double p0,  double p1a, double p2a, double p3a, 
                                                  double p1b, double p2b, double p3b, 
                                      bool bRandomOrder=1);
  void                  initRandom(unsigned N, double p);
  int                   input(string    fname, int soffset=0);// not yet working
  int                   inputGML(string fname, int soffset=0);
  //int                 inputDIMACS(string fname, int soffset=0);
          // utility functions
  virtual inline  void  init(TCMData stored=1, TCMData unstored=-1);
                  void  display(int bOutputLarge=0, bool bShowZeroes=1) const;
          // data member functions
  virtual TCMData       getMax(bool bAbs=0) const;
  virtual TCMData       getMin(bool bAbs=0) const;
          // mostly internal functions
  virtual inline bool   isValid(bool bEmptyOK=0) const;         // valid array check
  virtual        int    output_dimacs(string fname, bool useBinary = 0) const;

  // Standard public member functions
            TCMSparse(unsigned Cols=0, TCMData stored=1, TCMData unstored=1, 
                      string descr="");                     // constructor
            TCMSparse(TCMSparse &b);                        // copy constructor
  virtual  ~TCMSparse();                                    // destructor
  // public data member - only modify for non-standard sparse matrix
  // generalized constant for our cluster implementation.  It is the value of 
  // the non-stored member elements, i.e. zero for real sparse matrices.
  // It is bad form, but significantly simplifies the implementation for the 
  // cluster programs because I do not have to implement a derived class to 
  // satisfy etiquet.
            TCMData   storedValue;    // stored element value (no weights)
            TCMData   unstoredValue;  // return this value for unstored data
 protected:
  static const unsigned iMinBinarySearchSize = 10;  // not yet implemented
                      // fill the edges array given a TListArray
            void      copyEdges(TListArray<int> &edges, bool bCheckDupl=0);
}; // end class TCMSparse;

//------------------ TCMSparse Member Function Declarations -----------------

//Assign default static parameters
//TCMSparse::StrictErrorChecking = 1;

/*
//template <typename T> 
inline TCMData TCMSparse::operator()(unsigned i, unsigned j) const {
  // binary search on [][] operator for Array1D data structure
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in operator()()");
  #endif
  Array1D<int>  *pi = &(kij[i]);
  if(i==j)  return 0;
  else if(pi->getSize()<iMinBinarySearchSize || j<iMinBinarySearchSize) {  
    // use a linear search
    int k = 0, *pik = &((*pi)[0]);
    while((*pik)<j && k< pi->getSize()-1) { k++;  pik++; } // end while
    if((*pik)==j)  return storedValue;
    else           return unstoredValue; // not found if we get to here
  } else {
    if(i==j)  return 0;
    int  high = pi->getSize()-1, low = 0, mid;
    while(low<high) {
      mid = (low+high)/2;
      if((*pi)[mid] < j)
        low = mid + 1;
        else  high = mid; // can't be high = mid-1: here A[mid] >= value,
                          // so high can't be < mid if A[mid] == value
     } // end while
     // high == low, using high or low depends on taste 
     //if(low<N && (A[low]==j))  return low // found
     //else                      return -1 // not found        
     if(low<pi->getSize() && ((*pi)[low]==j)) return storedValue;   // found
     else                                     return unstoredValue; // not found        
   } // end else
}; // end array reference

//template <typename T> 
inline TCMData  TCMSparse::operator()(unsigned i, unsigned j) {
  // taken from wikipedia entry on binary searches since my original *textbook*
  // source had problems
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in operator()()");
  #endif
  Array1D<int>  *pi = &(kij[i]);
  if(i==j)  return 0;
  else if(pi->getSize()<iMinBinarySearchSize || j<iMinBinarySearchSize) {  
    // use a linear search
    int k = 0, *pik = &((*pi)[0]);
    while((*pik)<j && k< pi->getSize()-1) { k++;  pik++; } // end while
    if((*pik)==j)  return storedValue;
    else           return unstoredValue; // not found if we get to here
  } else {
    if(i==j)  return 0;
    int  high = pi->getSize()-1, low = 0, mid;
    while(low<high) {
      mid = (low+high)/2;
      if((*pi)[mid] < j)
        low = mid + 1;
        else  high = mid; // can't be high = mid-1: here A[mid] >= value,
                          // so high can't be < mid if A[mid] == value
     } // end while
     // high == low, using high or low depends on taste 
     //if(low<N && (A[low]==j))  return low // found
     //else                      return -1 // not found        
     if(low<pi->getSize() && ((*pi)[low]==j)) return storedValue;   // found
     else                                     return unstoredValue; // not found        
   } // end else
}; // end array reference
*/

//template <typename T> 
inline bool TCMSparse::isValid(bool bEmptyOK) const {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid(bEmptyOK);
  // Now TCMSparse specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid
/*
//template <typename T> 
inline bool TCMSparse::isValid() {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid();
  // Now TCMSparse specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid const
*/
//template <typename T> 
inline void TCMSparse::scale(TCMData a, TCMData b, TCMData aOld, TCMData bOld) {
  // scale the stored and unstored elements differently
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in scale()");
  #endif
  storedValue   = (storedValue/aOld)  *a;
  unstoredValue = (unstoredValue/bOld)*b;
  gammaa = a;   gammab = -b;
  return;
}; // end scale


//template <typename T> 
inline TCMSparse& TCMSparse::operator*=(TCMData x) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in operator*=()");
  #endif
  storedValue *= x;  unstoredValue *= x;
  return *this;
}; // end *=

//template <typename T> 
inline TCMSparse& TCMSparse::operator/=(TCMData x) {   // T division assignment
  #ifdef DEBUG_MODE
  if(!isValid())
    errorMsg("TCMSparse object is not valid in TListArray operator/=()");
  #endif
  if(x!=0) { storedValue /= x;  unstoredValue /= x; }
  else errorMsg("Division by zero in TCMSparse scale operation!\n");
  return *this;
}; // end /=



// Other useful array functions
//template <typename T> 
inline void TCMSparse::init(TCMData stored, TCMData unstored) {
  // reinitialize the Array but not the kijes
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in init()");
  #endif
  storedValue = stored;  unstoredValue = unstored;
  gammaa      = stored;  gammab        = -unstored;
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
  return;
};


// non-member functions
inline ostream& coutNodeConnections(ostream &fo, TCMatrix &CM, unsigned wNode) {
  cout << "Node " << wNode << " is connected to:  ";
  // now output CMatrix
  if(wNode<CM.getSize())  cout << CM.kij[wNode] << "\n";
  return fo;
} // end coutCMatrix



} // end namespace MSL
#endif

