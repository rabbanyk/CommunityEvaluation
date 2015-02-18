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
MSL - My (or Math) Scientific Library - TCMSparseW
See MSL_Array.h for general comments and some documentation.

TCMSparseW Specific Notes:
These notes handle only behavior specific to the TCMSparseW implementation.

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

#ifndef MSL_TCMSPARSEW_H
#define MSL_TCMSPARSEW_H

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
class TCMDataNode;
template <typename T> class Array1D;  // just need a forward declaration of type

// prototype some classes
//template <typename T> class TCMSparseW;
//template <typename T> class  TList;
//template <typename T> class  TListArray;
// prototype a few functions from other classes???
//template <typename T> inline int TListArray<T>::getSize() const;
//template <typename T> TCMSparseW  operator*(T x, const TCMSparseW& a);
//template <typename T> istream&     operator>>(istream& fin, TCMSparseW& a);

#ifndef MSL_TSPARSE_T
#define MSL_TSPARSE_T
// define some common vector types
//typedef TCMSparseW<TFloat> TSparseFloat;// TFloat  array type definition
//typedef TCMSparseW<int>    TCMSparseWInt;     // integer array type definition
#endif


//------------------------------------------------------------------------------
//template <typename T> 
class TCMSparseW : public TCMatrix {
 // Memory allocation is dynamic
 //friend istream&    operator>><T>(istream& fin, TCMSparseW& a);
 public:
          // overloaded operators
          TCMSparseW& operator=(const TCMSparseW& b);
          // a = b for TListArray defines a square matrix of size b.getSize()
          // using default values of storedValue=1 and unstoredValues=0
          // use set(...) for a non-square matrix
          TCMSparseW& operator=(const TListArray<int>& b);
          TCMSparseW& operator*=(TCMData x);          // *= by scalar
          TCMSparseW& operator/=(TCMData x);          // /= by scalar
          // operator() is a constant reference - cannot change elements
          // does not use op[][] due to structure constraints for binary search
          //TCMData  operator()(unsigned col, unsigned row) const;
          //TCMData  operator()(unsigned col, unsigned row);
  virtual TCMData  operator()(unsigned col, unsigned row) const;
  virtual TCMData  operator()(unsigned col, unsigned row);
         // set() sets the specified element to the storedValue adding the 
         // element location if it does not exist (expensive operation)
         //bool           set(unsigned col, unsigned row);
         // copy a set of indices to the matrix, that is set a whole kij list
  //inline void           set(unsigned col, Array1D<int> &a);
  //inline void           set(unsigned col, TList<int> &a);
         // returns true if unset or false if element was not being stored
         //bool           unset(unsigned col, unsigned row);
  // other functions
  virtual void   scale(TCMData a, TCMData b, TCMData aOld=1, TCMData bOld=1);
         // initialize the sparse matrix in various ways
          void   initByClusters(TClusterList &c, int N, double pin, double pout, 
                              TCMData lowW = 1, TCMData highW = 1);
          void   initByClusters2(TClusterList &c, int N, double pin, double pout, 
            int &nEdgesIn, int &nEdgesOut, TCMData lowW=1, TCMData highW=1);
          void   initByClustersSparse(TClusterList &c, double pin, double pout, 
                              TCMData lowWeight = 1, TCMData highWeight = 1);
          void   initByClustersRS(TClusterList &c, unsigned N,  
                              double pin, double pout, unsigned long long ISeed, 
                              TCMData lowW = 1, TCMData highW = 1);
          void   initByClustersPowerRS(TClusterList &c, unsigned N, 
                            double pin, double kMin, double kMax, double alpha, 
                            int &nEdgesIn, int &nEdgesOut, unsigned long long ISeed, 
                            TCMData  bottomWeight = 1, TCMData topWeight = 1);
          void   initHeterogeneous(TClusterList &c, TClusterList &merged, int N, 
                              Array1D<int> &v, double p2, double p1, double p0, 
                              TCMData lowWeight = 1, TCMData highWeight = 1);
          void   initHierarchy(unsigned N, unsigned n2, unsigned n1, unsigned n0,
                              double p3,  double p2,   double p1,   double p0,
                              bool bRandomizeNodes=1, TCMData lowW = 1, TCMData highW = 1);
         // initialize a 'staggered' hierarchy
          void   initSHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
                              double p0,  double p1a, double p2a, double p3a, 
                              double p1b, double p2b, double p3b, bool bRandOrder=1, 
                              TCMData lowWeight = 1, TCMData highWeight = 1);
  inline  void   initRandom(unsigned N, double p, 
                 TCMData lowW = 1, TCMData highW = 1, bool bN2 = 0);
          void   initRandomLL(unsigned N, double p, 
                 TCMData lowW = 1, TCMData highW = 1, bool bN2 = 0);
          void   initRandomRS(unsigned N, double p, unsigned long long ISeed, 
                 TCMData lowW = 1, TCMData highW = 1);
   inline bool   insertNode(int iN, int jN, bool checkDuplicates = 1);
          void   initRandomWeights(TCMData lowW = 1, TCMData highW = 1);
          // for the input routines, specialized matrix types are not finished
          int    input(string fname, int soffset=0, bool bDirected=0);
          int    inputGML(string fname, int soffset=0, 
                          bool bDirected=0, bool bSymmetrize=0); 
          int    inputLFRDAT(string fname, int N);
          int    exportGML(string fname, bool bDirected=0);
  //      int    inputDIMACS(string fname, int soffset=0);
          // utility functions
          // initializing stored data overwrites the current weights
          void   destroy();
  virtual void   init(TCMData stored=1, TCMData unstored=-1);
  virtual void   display(int bOutputLarge=0, bool bShowZeroes=1) const;
          void   disp(int bOutputLarge, bool bShowZeroes, bool bForceOutput) const;
          // data member functions
  virtual TCMData getMax(bool bAbs=0) const;
  virtual TCMData getMin(bool bAbs=0) const;
          // mostly internal functions
  virtual inline bool   isValid(bool bEmptyOK=0) const;    // valid array check
  virtual        int    output_dimacs(string fname, bool useBinary = 0) const;

  // Standard public member functions
           TCMSparseW(unsigned Cols=0, TCMData unstored=-1, string describe="");
           TCMSparseW(TCMSparseW &b);                     // copy constructor
  virtual ~TCMSparseW();                                  // destructor
  // public data member - only modify for non-standard sparse matrix
  // generalized constant for our cluster implementation.  It is the value of 
  // the non-stored member elements, i.e. zero for real sparse matrices.
  // It is bad form, but significantly simplifies the implementation for the 
  // cluster programs because I do not have to implement a derived class to 
  // satisfy etiquet.
         TCMData   unstoredValue; // return this value for unstored data
         int       copyEdges(TListArray<int> &edges, int checkDupl=0);
         int       copyEdges(TListArray<TCMDataNode> &edges, int checkDupl, 
                             bool bIgnoreThisParameter);
 //protected:
            Array1D<TCMData>  *data; // weights (indices kij are in base class)
  static const unsigned iMinBinarySearchSize = 10;  // not yet implemented
                      // fill the edges array given a TListArray
}; // end class TCMSparseW;

//------------------ TCMSparseW Member Function Declarations -----------------

//Assign default static parameters
//TCMSparseW::StrictErrorChecking = 1;

inline void TCMSparseW::initRandom(unsigned N, double p, 
                 TCMData lowW, TCMData highW, bool bN2) {
  // it just makes things easier to have a 'generic' initRandom
  initRandomLL(N,p,lowW,highW,bN2);  
};

inline bool TCMSparseW::insertNode(int iN, int jN, bool bCheckDuplicates) {
  // Uses an optimized insertion sort logic to insert a single node into a 
  // sorted integer list.  The intention is that the function is used to 
  // `build' the integer list as it is created rather than to be an iterated 
  // function call (although it is inline).
  int  j, length = kij[iN].getSize();      // current length of the iN vector
  int *pJ, *pPrevJ;                        // pointers to j'th, previous j'th #s
  bool bInserted;
  
  if(length==0) { // just add it to the list
    kij[iN].resize(1,0,0,1);      // fake resize the vector
    kij[iN][0] = jN;
    return 1;
  } // end length check

  if(bCheckDuplicates) {
    // do this as a separate case to avoid redundant loops
    //cout << "A " << flush;  // debugging
    pPrevJ = &(kij[iN][length-1]);
    //cout << "B " << flush;  // debugging
    j = length;                            // j starts at candidate location
    //cout << kij[iN] << endl;  // debugging
    while(j>0 && (*pPrevJ)>jN) { pPrevJ--;  j--; } // end while j
    //cout << "C with pPrevJ = " << (*pPrevJ) << flush;  // debugging
    // j ends with the location of the addition

    if((*pPrevJ)==jN) {
      //cout << "D " << flush;  // debugging
      bInserted = 0;        // ignore addition of jN
    } else {                                 // copy the data now
      kij[iN].resize(length+1,0,0,1);      // fake resize the vector
      pJ = &(kij[iN][length]);
      pPrevJ = pJ - 1;
      //cout << "E " << flush;  // debugging
      // move data element up to next spot element by element
      for(int k=length; k>j; k--) {
       *pJ = *pPrevJ;                      // move element up to next spot
        pJ  = pPrevJ;  pPrevJ--;           // decrement j'th, (j-1)'th pointers
      } // end for k
      //cout << "F " << flush;  // debugging
      *pJ = jN;
      bInserted = 1;
      //cout << "G " << flush;  // debugging
    } // end else
    //cout << "H " << flush;  // debugging
  } else { // else not bCheckDuplicates
    kij[iN].resize(length+1,0,0,1);  // fake resize of this vector
    pJ = &(kij[iN][length]); // dummy starting location is incremented in loop
    pPrevJ = pJ - 1;         // we need j and j-1 pointer index locations stored
    j = length;              // j and i counters start at same number index
    while(j>0 && (*pPrevJ)>jN) {
     *pJ = *pPrevJ;          // move data element up to next spot
      pJ = pPrevJ;  pPrevJ--;// decrement j'th and (j-1)'th pointers
      j--;
    } // end while j
    *pJ = jN;  // ok, finally add the integer to the list
    bInserted = 1;
  } // end else bCheckDuplicates
  
  //cout << "I " << flush;  // debugging
  return bInserted;  // was the node actually inserted
} // end insertNode


//template <typename T> 
inline bool TCMSparseW::isValid(bool bEmptyOK) const {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid(bEmptyOK);
  // Now TCMSparseW specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid
/*
//template <typename T> 
inline bool TCMSparseW::isValid() {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid();
  // Now TCMSparseW specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid const
*/


// --------------- End TCMSparseW Member Function Declarations ----------------
// ----------------------------------------------------------------------------

// ------------------------- TCMDataNode class --------------------------------
// TCMDataNode for the initialization routines.  The TCMSparseW c(i,j) reference
// operator requires a sort by index due to the binary search; but for input
// data files, the weights need to be sorted attached to the indices also
// (assuming we do not want to read the weighted data file twice).
// The easiest way to accomplish this is to use a simple encapsulated data
// structure for the TList and define the corresponding boolean operators.
class TCMDataNode {
 public:
         // overloaded operators
  inline TCMDataNode& operator=(const TCMDataNode& b);
  // the boolean operators are defined on the index rather than the weight
  // because of the binary search in the reference TCMSparse::operator(i,j)
  inline bool operator> (const TCMDataNode& b) { return index>b.index;  };
  inline bool operator< (const TCMDataNode& b) { return index<b.index;  };
  inline bool operator>=(const TCMDataNode& b) { return index>b.index || 
                                                        index==b.index; };
  inline bool operator<=(const TCMDataNode& b) { return index<b.index || 
                                                        index==b.index; };
  inline bool operator==(const TCMDataNode& b) { return index==b.index; };
  inline TCMDataNode operator+ (const TCMDataNode& b);
  inline TCMDataNode operator+ (const int    i);
  inline TCMDataNode operator+ (const double d);
  // Standard public member functions
  inline      TCMDataNode(int theIndex=0, TCMData val=1); // constructor
  inline      TCMDataNode(TCMDataNode &b);                // copy constructor
  //         ~TCMDataNode() {};                           // destructor
  // keep the data public for the implementation classes
         int      index;
         TCMData  value;
}; // end class TCMDataNode

inline TCMDataNode::TCMDataNode(int theIndex, TCMData val)  { 
  index = theIndex; value = val; };     // constructor
inline TCMDataNode::TCMDataNode(TCMDataNode &b) { 
  index = b.index;  value = b.value; }; // copy const
inline TCMDataNode& TCMDataNode::operator=(const TCMDataNode& b) { 
  index = b.index;  value = b.value;  return *this; };
// the operators on TCMDataNode only change the weight (value) not the indices
inline TCMDataNode TCMDataNode::operator+(const TCMDataNode& b) { 
  TCMDataNode c = (*this);
  c.value += b.value;
  return c;  
};
// the operator functions are needed mostly to maintain basic compatibility
// with TList operations such as display(...)
inline TCMDataNode TCMDataNode::operator+(const int i) { 
  TCMDataNode c = (*this);
  c.value += (TCMData)i;
  return c;  
};
inline TCMDataNode TCMDataNode::operator+(const double i) { 
  TCMDataNode c = (*this);
  c.value += (TCMData)i;
  return c;  
};

// additional TCMDataNode functions
inline ostream& operator<<(ostream& fout, TCMDataNode& a) {
  // Outputs the weight of the TCMData node only
  fout << a.index << "(" << a.value << ")";
  return fout;
};


// ----------------------- End TCMDataNode class ------------------------------

} // end namespace MSL
#endif

