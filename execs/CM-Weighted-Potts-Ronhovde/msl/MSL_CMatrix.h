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
June 2008
MSL - My (or Math) Scientific Library - TCMatrix base class

Comments:
-
  
**************************************************************************************/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TCMATRIX_H
#define MSL_TCMATRIX_H

//#define DEBUG_MODE

#include <string>

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "MSL_Array1D_Template.h"
#endif

namespace MSL {

// define the type of data used in the CMatrix.  This is done as a typedef
// rather than using templates because of compiler problems (not understanding
// the errors) and the fact that it is unlikely that I will need to instance 
// both an int and a double version of the matrix in the same program.
typedef int TCMData;  

// define constants
const int ValidCMatrixFlag = 54321;

// need to prototype template class Array1D<T> for TCMatrix data member
template <typename T> class Array1D;

class TCMatrix {
 public:
  // data access methods
  inline  double      getDensity() const { return p;      };
  inline  unsigned    getNEdges()  const { return nEdges; };
  inline  unsigned    getNNodes()  const { return nNodes; };
  inline  unsigned    getSize()    const { return nNodes; };
  // convenience functions
  inline  unsigned    getN()       const { return nNodes; };
  inline  unsigned    getL()       const { return nEdges; };
  inline  double      getND()      const { return (double)nNodes; };
  inline  double      getLD()      const { return (double)nEdges; };
  //inline  unsigned  getCols()    const { return nNodes; };  // compatibility
  //inline  unsigned  getRows()    const { return nNodes; };  // compatibility

  // the reference operators are implemented as overloaded functions rather
  // than virtual function because there was a *huge* performance hit
  virtual TCMData     operator()(unsigned col, unsigned row) const = 0;
  virtual TCMData     operator()(unsigned col, unsigned row)       = 0;
  //      TCMData     operator()(unsigned col, unsigned row) const {return 0.0;};
  //      TCMData     operator()(unsigned col, unsigned row)       {return 0.0;};
  // ki and kij return the degree of the node and the connected nodes
  inline  int         ki(unsigned i) const;  // return the degree of node i
                      // return the jth node connected to i
  inline  int         ki(unsigned i, unsigned j) const;

  // virtual functions 
  virtual void        scale(TCMData a, TCMData b, TCMData c=1, TCMData d=1) = 0;
  virtual void        init(TCMData connected=1, TCMData unconnected=-1) = 0;
  virtual TCMData     getMax(bool bAbs=0) const = 0;
  virtual TCMData     getMin(bool bAbs=0) const = 0;
  //virtual int       output_dimacs(string fname, bool useBinary=0) const = 0;
  virtual void        display(int bOutputLarge=0, bool bShowZeroes=1) const = 0;

  // utility functions
  //      void        display(int bOutputLarge=0, bool bShowZeroes=1) const;
          void        resize(unsigned Cols);               // DELETES ALL DATA!
  inline  bool        isValid(bool bEmptyOK=0) const;

  // standard member functions
                      TCMatrix(unsigned Cs, TCMData a=1, TCMData b=1, string desc="");
                      TCMatrix(TCMatrix& a);
  virtual            ~TCMatrix();
  // other data elements
  static const bool   VerboseErrorReports = 0;  // not used much at the moment
  string              description;
  Array1D<int>       *kij; // the node connection list is public but functions
                           // are provided and should be used most of the time

  // these are the Potts model weights almost equivalent to the stored or 
  // unstored values (except that unstored value = -gammab) for the sparse 
  // matrices.  They are the current scale factors being used in the matrix.
  TCMData             gammaa, gammab;   // not used at the moment
 protected:
  void           destroykij();
  void           createkij(unsigned size);
  unsigned       nNodes, nEdges; // number of nodes and edges in the matrix
  double         p;              // edge density
  //int            nodeOffset;   // the integer offset for displayed information
  // some functions are off limits for later possible usage
  // Data elements
  unsigned       Cols;    // number of columns - a CMatrix is always square
  int            ValidArray;
}; // end virtual class TCMatrix


//template <typename T> 
inline bool TCMatrix::isValid(bool bEmptyOK) const {
  #ifdef DEBUG_MODE
  if(p<0.0) {
    warningMsg("TCMatrix object has an edge density < 0.0?");  return 0; }
  else if(Cols!=nNodes) {
    warningMsg("TCMatrix object has an inconsistent number of nodes w/nNodes = "
             +itos(nNodes)+" and Cols = "+itos(Cols));    return 0; }
  else if(nEdges>(nNodes*(nNodes-1))) {
    warningMsg("TCMatrix object has too many edges.\nnEdges = "+itos(nEdges)
             +" but the maximum number is "+itos(nNodes*(nNodes-1))
             +" for "+itos(nNodes)+" nodes");  return 0; } 
  // after above genuine inconsistencies, check for simply an empty array
  else if(ValidArray==0 && Cols==0 && nNodes==0 && kij==NULL)  return 1;
  // else continue with error checks
  else if(ValidArray!=ValidCMatrixFlag) {
    warningMsg("TCMatrix object has an incorrect signature flag of "
             +itos(ValidArray));  return 0; }
  //else if(Cols!=Rows)  errorMsg("TCMSparse object is not a square matrix?");
  return 1;  // if we make it to here, it is ok
  #else
  return 1;
  #endif
}; // end isValid const


/*
//template <typename T> 
inline void scale(TCMData a, TCMData b, TCMData aOld=1, TCMData bOld=1) {
  #ifdef DEBUG_MODE
  warningMsg("TCMatrix scale() is a dummy function.  Use child version.");
  #endif
}; // end scale
*/


//template <typename T> 
inline int TCMatrix::ki(unsigned i) const { 
  // returns the degree of node i
  #ifdef DEBUG_MODE
  if(!isValid())    errorMsg("Invalid TCMatrix object in ki(i). ");
  else if(i>=Cols)  errorMsg("Index out of bounds in TCMatrix ki() with i = "
                        +itos(i)+" with a maximum of "+itos(Cols));
  #endif
  return kij[i].getSize();
}; // end ki()


//template <typename T> 
inline int TCMatrix::ki(unsigned i, unsigned j) const { 
  // returns the jth node connected to i
  #ifdef DEBUG_MODE
  if(!isValid())    errorMsg("Invalid TCMatrix object in ki(i,j)");
  else if(i>=Cols)  errorMsg("Index out of bounds in TCMatrix ki(i) with i = "
                        +itos(i)+" with a maximum of "+itos(Cols));
  else if(j>=kij[i].getSize()) 
         errorMsg("Index out of bounds in TCMatrix ki(i,j) with j = "
                  +itos(j)+" with a maximum of "+itos(kij[i].getSize()));
  #endif
  return kij[i][j];
}; // end ki()


} // end namespace MSL
#endif
