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
MSL - My (or Math) Scientific Library - TCMDense
See MSL_Array.h for general comments and some documentation.

TCMDense Specific Notes:
These notes handle only behavior specific to the TCMDense implementation.

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

#ifndef MSL_TCMDENSE_H
#define MSL_TCMDENSE_H

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
#ifndef MSL_TARRAY2D_H
#include "MSL_Array2D_Template.h"
#endif

namespace MSL {

class TCluster;
class TClusterList;
//template <typename T> class Array2D;
// prototype of init potential for initVij
double Vij(double x1, double y1, double x2, double y2);  

// prototype some classes
//template <typename T> class TCMDense;
//template <typename T> class  TList;
//template <typename T> class  TListArray;
// prototype a few functions from other classes???
//template <typename T> inline int TListArray<T>::getSize() const;
//template <typename T> TCMDense  operator*(T x, const TCMDense& a);
//template <typename T> istream&     operator>>(istream& fin, TCMDense& a);

#ifndef MSL_TSPARSE_T
#define MSL_TSPARSE_T
// define some common vector types
//typedef TCMDense<TFloat> TSparseFloat;// TFloat  array type definition
//typedef TCMDense<int>    TCMDenseInt;     // integer array type definition
#endif


//------------------------------------------------------------------------------
//template <typename T> 
class TCMDense : public TCMatrix {
 // Memory allocation is dynamic
 //friend istream&    operator>><T>(istream& fin, TCMDense& a);
 public:
          // overloaded operators
          TCMDense& operator=(const TCMDense& b);
          // a = b for TListArray defines a square matrix of size b.getSize()
          // using default values of storedValue=1 and unstoredValues=0
          // use set(...) for a non-square matrix
          TCMDense& operator=(const TListArray<int>& b);
  inline  TCMDense& operator*=(TCMData x) { data *= x; }; // *= by scalar
  inline  TCMDense& operator/=(TCMData x) { data /= x; }; // /= by scalar
          // operator() is a constant reference - cannot change elements
          // does not use op[][] due to structure constraints for binary search
  inline  TCMData  operator()(unsigned col, unsigned row) const;
  inline  TCMData  operator()(unsigned col, unsigned row);  // error here
  // other functions
  virtual void   scale(TCMData a, TCMData b, TCMData aOld=1, TCMData bOld=1);
          // initialize matrix
          int    input(string fname, int soffset=0, bool bDirected=0);
          // directed matrices are automatically symmetrized due to legacy issues
          int    inputTextMatrix(string fname, int soffset=0, bool bDirected=0);
                                 //bool bDirected=0, bool bSymmetrize=0);
          int    inputGML(string fname, int soffset=0, bool bDirected=0);
                                 //bool bDirected=0, bool bSymmetrize=0);
          int    exportGML(string fname, bool bDirected=0);
          void   initRandom(unsigned N, double p, 
                            TCMData lowW = 1, TCMData highW = 1, bool bN2 = 0);
          void   initHeterogeneous(TClusterList &c, TClusterList &merged, int N, 
                                  Array1D<int> &v, double p2, double p1, double p0, 
                                  TCMData lowWeight = 1, TCMData highWeight = 1);
          void   initVij(TMatrixFloat &rn, double t=0.0, TCMData a=1, TCMData b=-1);
          void   initByClusters2(TClusterList &c, int N, double pin, double pout,
                   int &nEdgesIn, int &nEdgesOut, TCMData bottomW=1, TCMData topW=1);
          void   initByClusters(TClusterList &c, int N, double pin, double pout,
                   TCMData bottomW=1, TCMData topW=1);
          void   initByClustersPower(TClusterList &c, unsigned N, 
                          double pin, double kMin, double kMax, double alpha, 
                          int &nEdgesIn, int &nEdgesOut);
          void   initByCMArray(Array2D<int> &a, int N, string d = "");
          void   initNoiseExp(unsigned N, double kMean, bool bWeighted=0, double wMean=1.0);

          int    addNoiseExp(TClusterList& s, double pout, double kMean);
          int    addNoisep(TClusterList& s, double pout);
          int    addNoisePower(TClusterList& s, double kMin, double kMax, double alpha);
          // utility functions
          // initializing stored data overwrites the current weights
  virtual void   init(TCMData stored=1, TCMData unstored=-1);
  virtual void   display(int bOutputLarge=0, bool bShowZeroes=1) const;
          // data member functions
  virtual inline TCMData getMax(bool bAbs=0) const;
  virtual inline TCMData getMin(bool bAbs=0) const;
          // mostly internal functions
  virtual inline bool    isValid(bool bEmptyOK=0) const;         // valid array check
          int            copyEdges(TListArray<int> &edges);
                 void    copyEdges(Array1D<int> &degreeList);

  // Standard public member functions
            TCMDense(unsigned Cols=0, string descr="");   // constructors
            TCMDense(Array2D<TCMData> &a, unsigned Cs, string descr="");
            TCMDense(TCMDense &b);                        // copy constructor
  virtual  ~TCMDense();                                   // destructor

  protected:
            Array2D<TCMData>  data;
}; // end class TCMDense;

//------------------ TCMDense Member Function Declarations -----------------

//Assign default static parameters
//TCMDense::StrictErrorChecking = 1;

inline TCMData TCMDense::operator()(unsigned col, unsigned row) const {
  return data[col][row]; 
};
inline TCMData TCMDense::operator()(unsigned col, unsigned row) {
  return data[col][row]; 
};

inline TCMData TCMDense::getMax(bool bAbs) const { return data.getMax(bAbs); };
inline TCMData TCMDense::getMin(bool bAbs) const { return data.getMin(bAbs); };

//template <typename T> 
inline bool TCMDense::isValid(bool bEmptyOK) const {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid(bEmptyOK);
  // Now TCMDense specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid
/*
//template <typename T> 
inline bool TCMDense::isValid() {
  #ifdef DEBUG_MODE
  bool result = TCMatrix::isValid();
  // Now TCMDense specific validation check
  return result;
  #else
  return 1;
  #endif
}; // end isValid const
*/
// --------------- End TCMDense Member Function Declarations ----------------
// ----------------------------------------------------------------------------

} // end namespace MSL
#endif

