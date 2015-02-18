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

//#ifndef MSL_TCMATRIX_H

//#define DEBUG_MODE

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "MSL_Array1D_Template.h"
#endif
#ifndef MSL_TCMATRIX_H
#include "MSL_CMatrix.h"
#endif

namespace MSL {

// define the type of data used in the CMatrix.  This is done as a typedef
// rather than using templates because of compiler problems (not understanding
// the errors) and the fact that it is unlikely that I will need to instance 
// both an int and a double version of the matrix in the same program.
typedef int TCMData;

// need to prototype template class Array1D<T> for TCMatrix data member
template <typename T> class Array1D;


//template <typename T> 
TCMatrix::TCMatrix(unsigned Cs, TCMData a, TCMData b, string descr) {
  #ifdef DEBUG_MODE
  debugMsg("Constructing TCMatrix object with ");
  #endif
  Cols   = Cs;  // a square matrix
  nNodes = Cs;  nEdges = 0;
  p = 0.0;
  description = descr;
  gammaa = a;  gammab = b;
  if(Cols>0) {
    #ifdef DEBUG_MODE
    debugMsg("size "+itos(Cs)+"... ");
    #endif
    //kij = new Array1D<int>[Cols];   // Set up kij columns
    // allow initializing routines to set up kij, this is almost always required
    // anyhow, and it prevents an allocation and de-allocation of an empty array
    kij = NULL;  
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end un-initialized array creation
  else {
    kij = NULL;
    #ifdef DEBUG_MODE
    debugMsg("default size = 0... ");
    ValidArray = 0;
    warningMsg("Creating a size zero TCMatrix object");
    #endif
  } // end zero length if
  #ifdef DEBUG_MODE
  debugMsg("exiting with array flag = "+itos(ValidArray)+"... ");
  #endif
}; // end constructor


//template <typename T> 
TCMatrix::TCMatrix(TCMatrix &b) { // copy constructor
  #ifdef DEBUG_MODE
  if(!isValid(1)) 
    warningMsg("TCMatrix object is not valid in copy constructor");
  debugMsg("Entering TCMSparse<T> copy constructor... ");
  #endif
  if(b.getSize()>0) {
    kij = new Array1D<int>[b.getSize()];   // Set up kij columns
    // Now, initialize array and the rest of the data
    for(int i; i<b.getSize(); i++)  kij[i] = b.kij[i];  // use Array1D<T> '='
    Cols   = b.Cols;
    nNodes = b.nNodes;  nEdges = b.nEdges;  p = b.p;
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
    description = b.description;
    gammaa = b.gammaa;  gammab = b.gammab;
  } // end un-initialized array creation
  else { 
    kij  = NULL;
    Cols   = b.Cols;
    nNodes = b.nNodes;  nEdges = b.nEdges;  p = b.p;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    errorMsg("Attempting to copy a zero length TCMatrix object?");
    #endif
  } // end else
  #ifdef DEBUG_MODE
  debugMsg("exiting TCMSparse<T> copy constructor... ");
  #endif
}; // end TCMSparse<T> copy constructor


//template <typename T> 
TCMatrix::~TCMatrix() {
  #ifdef DEBUG_MODE
  debugMsg("Entering TCMatrix destructor...\n",brown);
  if(!isValid(1))  warningMsg("TCMatrix object is not valid in destructor");
  #endif
  destroykij();
  #ifdef DEBUG_MODE
  debugMsg("exiting TCMatrix destructor.\n",brown);
  #endif
}; // end TCMSparse<T> destructor


//template <typename T> 
void  TCMatrix::destroykij() {
  #ifdef DEBUG_MODE
  //if(!isValid(1))  warningMsg("TCMatrix object is not valid in destroykij()");
  #endif
  delete[] kij;  kij = NULL;
  //Rows = -1;  Cols = -1;  p = -1.0;
  #ifdef DEBUG_MODE
  ValidArray = 0;
  #endif
}; // end TCMSparse<T> destructor


void  TCMatrix::createkij(unsigned size) {
  kij = new Array1D<int>[size];
  //Rows = -1;  Cols = -1;  p = -1.0;
  #ifdef DEBUG_MODE
  ValidArray = 0;  
  #endif
}; // end 

//template <typename T> 
void TCMatrix::resize(unsigned Cs) {
  // Crude resize function to "get things working" at the moment
  // At the moment it DELETES ALL DATA!!!
  // first delete existing array data
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMatrix object is not valid in resize()");
  #endif
  Cols   = Cs; // square matrix
  nNodes = Cs;  nEdges = 0;
  p      = 0.0;
  if(Cs>0) {
    delete[] kij;
    kij = new Array1D<int>[Cs];
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end if
  else {
    kij = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return;
}; // end resize()

/*
//template <typename T> 
void TCMatrix::display(int bOutputLarge, bool bShowZeroes) const {
  int i, j, cij, N = nNodes;

  //cout << blueText << " CM " << CMatrix[1][0] << normText;
  if(N>96 || !bShowZeroes)  return;  // debugging
  //if(N>64 && !(bool)bOutputLarge)  return;  // debugging
  //if(bOutputLarge < 2) return;  // debugging

  if(bOutputLarge>0) { // scan for largest integer
    int maxInt = 0;
    for(i=0; i<N; i++)  for(j=0; j<N; j++) 
      if(abs((*this)(i,j))>maxInt)  maxInt = abs((*this)(i,j));
    if(maxInt==0) { 
      warningMsg("Matrix is a zero matrix!  Exiting matrix output.");  
      return; 
    } // end if maxInt
    // now define spacers for text
    int spacingSize = (int)log10((double)maxInt) + 1;
    char *spacer;  spacer = new char[spacingSize+1];
    for(i=0; i<spacingSize; ++i)  spacer[i] = ' ';
    spacer[spacingSize] = '\0';
    
    // finally output matrix
    int iSize = (int)log10((double)N) + 1;
    int nSize;
    // print the line at top
    cout << normText << " " << itos(N,0,',',1,iSize,' ') << " |";
    if(spacingSize==1)  // need to account for too large values sometime!!!
      for(j=0; j<N; ++j)  cout << spacer << (j%10);
    else { // allow two or more digit column headers
      for(j=0; j<N; ++j) {
        if(j<10)  nSize = 0;
        else      nSize = (int)log10((double)j);
        cout << &spacer[nSize] << j;
      } // end for j
    } // end else 
    cout << "\n";
    for(j=0; j<(spacingSize+1)*N+iSize+5; ++j)  cout << "-";  
    cout << "\n";

    // now output CMatrix
    for(i=0; i<N; ++i) {
      cout << " " << itos(i,0,',',1,iSize,' ') << " |";
      for(j=0; j<N; ++j) {
        cij = (*this)(j,i);
        //if(CM(i,j)==1) cout << (i%10) << " ";
        if(cij>0) {
          if(cij<10)  nSize = 0;
          else        nSize = (int)log10((double)cij);  
          cout << &spacer[nSize] << cij; 
        } // end if
        else if(cij==0) {
               if(bShowZeroes) cout << spacer << greyText << '0' << normText;
               else            cout << spacer << ' ';
             } // end if cij==0
        else if(cij<0) {
               if(abs(cij)<10)  nSize = 0;
               else             nSize = (int)log10((double)abs(cij));  
               cout << redText << &spacer[nSize] << -cij << normText;
             } // end if cij<0
      } // end for j
      cout << "\n";
    } // end for i
    delete[] spacer;

  } // end if bOutputLarge
  else {  // just use original version - kept just in case
    int maxInt = 0;
    for(i=0; i<N; i++)  for(j=0; j<N; j++) 
      if(abs((*this)(i,j))>maxInt)  maxInt = (*this)(i,j);
    // now define spacers for text
    // ...
    int iSize = (int)log10((double)N) + 1;
    // print the line at top
    cout << normText << " " << itos(N,0,',',1,iSize,' ') << " | ";
    for(j=0; j<N; ++j)            cout << (j%10) << " ";  cout << "\n";
    for(j=0; j<2*N+iSize+5; ++j)  cout << "-";            cout << "\n";
    // now output CMatrix
    for(i=0; i<N; ++i) {
      cout << " " << itos(i,0,',',1,iSize,' ') << " | ";
      for(j=0; j<N; ++j) {
        cij = (*this)(j,i);
        //if(CM(i,j)==1) cout << (i%10) << " ";
        if(cij>0)        cout << cij << " ";
        else if(cij==0) {
               if(bShowZeroes) cout << greyText << '0' << normText << ' ';
               else            cout << "  ";
             } // end if cij==0
        else if(cij<0)   cout << redText << -cij << " " << normText;
      } // end for j
      cout << "\n";
    } // end for i
  } // end else
  return;
}; // end display
*/



} // end namespace MSL
//#endif
