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
1. For some reason operator<< is trying to call the pure virtual isValid() 
function; therefore, the error check is commented out for now.
*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//#ifndef MSL_TCMSPARSEW_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

//#define DEBUG_MODE

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
#ifndef MSL_TCMSparseW_H
#include "MSL_CMatrix_SparseW.h"
#endif
#ifndef MSL_TCLUSTERLIST_H
#include "../clusterclasses.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "./msl/MSL_Array1D_Template.h"
#endif

namespace MSL {

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
//------------------ TCMSparseW Member Function Declarations -----------------

//Assign default static parameters
//TCMSparseW::StrictErrorChecking = 1;

//template <typename T> 
TCMSparseW::TCMSparseW(unsigned Cs, TCMData unstoredWeight, string d) 
           : TCMatrix(Cs,1,-unstoredWeight,d) {
  // Raw constructor that skips the initialization
  // As the default constructor, it generates a TCMSparseW object that is a 
  // zero length real TCMSparseW object but officially an invalid one since it
  // holds no data.
  #ifdef DEBUG_MODE
  debugMsg("Constructing TCMSparseW object... ");
  #endif
  // now declare the sparse data
  // allow initializing routines to set up data, this is almost always required
  // anyhow, and it prevents an allocation and de-allocation of an empty array
  //data = new Array1D<TCMData>[Cs];  // create a null version of the data array
  data = NULL;  

  unstoredValue = unstoredWeight;   // set the value of the unstored elements
  #ifdef DEBUG_MODE
  debugMsg("done.\n");
  #endif
}; // end TCMSparseW constructor


//template <typename T> 
TCMSparseW::TCMSparseW(TCMSparseW &b) : TCMatrix(b){ // copy constructor
  #ifdef DEBUG_MODE
  debugMsg("Copying TCMSparseW object\n");
  if(!isValid(1)) 
    errorMsg("TCMSparseW object is not valid in copy constructor");
  cout << "Entering TCMSparseW destructor... " << endl;
  #endif
  data = new Array1D<TCMData>[Cols];  // create a null version of the data array
  for(unsigned i=0; i<Cols; i++)  data[i] = b.data[i];  // use optimized Array1D '='
  unstoredValue = b.unstoredValue;
}; // end TCMSparseW copy constructor


//template <typename T> 
TCMSparseW::~TCMSparseW() {
  #ifdef DEBUG_MODE
  //debugMsg("Entering TCMSparseW destructor... ",brown);
  if(!isValid(1))  errorMsg("TCMSparseW object is not valid in destructor");
  #endif
  debugMsg("Entering TCMSparseW destructor... ",brown);
  delete[] data;  data = NULL;
  //Rows = 0;  Cols = 0;  p = 0.0; // i.e. do nothing here
  //#ifdef DEBUG_MODE
  //ValidArray = 0; let TCMatrix destructor set this due to reversed call order
  debugMsg("exiting TCMSparseW destructor.\n",brown);
  //#endif
}; // end TCMSparseW destructor


void TCMSparseW::destroy() {
  #ifdef DEBUG_MODE
  //debugMsg("Entering TCMSparseW destructor... ",brown);
  if(!isValid(1))  errorMsg("TCMSparseW object is not valid in destroy()");
  debugMsg("Entering TCMSparseW destructor... ",brown);
  #endif
  destroykij();
  delete[] data;  data = NULL;
  nNodes = 0, nEdges = 0, Cols = 0, p = 0.0;
  #ifdef DEBUG_MODE
  ValidArray = 0;
  debugMsg("exiting TCMSparseW destructor.\n",brown);
  #endif
}; // end TCMSparseW destroy()

//template <typename T> 
void TCMSparseW::scale(TCMData a, TCMData b, TCMData aOld, TCMData bOld) {
  // scale the stored and unstored elements differently
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in scale()");
  #endif
  // this current implementation could be inefficient; but we allow integer 
  // weights, and we would have problems if we use integer division
  // try to use use optimized Array1D routines
  if(aOld!=1)  for(unsigned i=0; i<Cols; i++)  data[i] /= aOld; 
  for(unsigned i=0; i<Cols; i++)               data[i] *= a; 
  unstoredValue = (unstoredValue/bOld)*b;
  gammaa = a;    gammab = -b;
  return;
}; // end scale


//template <typename T> 
inline TCMSparseW& TCMSparseW::operator*=(TCMData x) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in operator*=()");
  #endif
  for(unsigned i=0; i<Cols; i++)  data[i] *= x; // use Array1D optimized routine
  unstoredValue *= x;
  return *this;
}; // end *=

//template <typename T> 
inline TCMSparseW& TCMSparseW::operator/=(TCMData x) {   // T division assignment
  #ifdef DEBUG_MODE
  if(!isValid())
    errorMsg("TCMSparseW object is not valid in TListArray operator/=()");
  #endif
  if(x!=0) { 
    for(unsigned i=0; i<Cols; i++)  data[i] /= x; // use Array1D optimized routine
    unstoredValue /= x; 
  } // end if
  else errorMsg("Division by zero in TCMSparseW scale operation!\n");
  return *this;
}; // end /=



// Other useful array functions
//template <typename T> 
inline void TCMSparseW::init(TCMData stored, TCMData unstored) {
  // reinitialize the Array but not the kijes
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in init()");
  #endif
  // use Array1D optimized routines
  for(unsigned i=0; i<Cols; i++)  data[i].init(stored);
  unstoredValue = unstored;
  gammaa = stored;  gammab = -unstored;
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
  return;
};


//template <typename T> 
TCMSparseW& TCMSparseW::operator=(const TCMSparseW& b) {
  // assignment for 2 arrays
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("TCMSparseW object is not valid in operator=()");
  #endif
  if(b.getSize()>Cols)                      // assumes a square matrix here
    errorMsg("TCMSparseW row too large in '=' exceeds specified Rows size.");
  delete[] kij;  delete[] data;
  unstoredValue = b.unstoredValue;
  Cols          = b.Cols;
  nNodes        = b.nNodes;       nEdges = b.nEdges;
  p             = b.p;
  //gammaa        = b.gammaa;       gammab = b.gammab;
  description   = b.description;
  if(Cols>0) {
    kij  = new Array1D<int>[Cols];                // Set up kij columns
    data = new Array1D<TCMData>[Cols];            // Set up data columns
    for(unsigned i=0; i<Cols; i++) {    // use Array1D<T> optimized routines
      kij[i] = b.kij[i];
      data[i] = b.data[i];
    } // end for i
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end if
  else {
    kij = NULL;  data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return *this;
}; // end =


//template <typename T> 
TCMData TCMSparseW::operator()(unsigned i, unsigned j) const {
  // binary search on [][] operator for Array1D data structure
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in operator()");
  #endif
  Array1D<int>  *pi = &(kij[i]);
  if(i==j)  return 0;
  else if(pi->getSize()<iMinBinarySearchSize || j<iMinBinarySearchSize) {  
    // use a linear search
    int k = 0, *pik = &((*pi)[0]);
    while((*pik)<j && k< pi->getSize()-1) { k++;  pik++; }
    if((*pik)==j)  return data[i][k];
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
     if(low<pi->getSize() && (*pi)[low]==j)  return data[i][low];  // found
     else                                    return unstoredValue; // not found        
   } // end else
}; // end array reference

//template <typename T> 
TCMData  TCMSparseW::operator()(unsigned i, unsigned j) {
  // taken from wikipedia entry on binary searches since my original *textbook*
  // source had problems
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in operator()");
  #endif
  Array1D<int>  *pi = &(kij[i]);
  if(i==j)  return 0;
  else if(pi->getSize()<iMinBinarySearchSize || j<iMinBinarySearchSize) {  
    // use a linear search
    int k = 0, *pik = &((*pi)[0]);
    while((*pik)<j && k< pi->getSize()-1) { k++;  pik++; }
    if((*pik)==j)  return data[i][k];
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
     if(low<pi->getSize() && (*pi)[low]==j)  return data[i][low];  // found
     else                                    return unstoredValue; // not found        
   } // end else
}; // end array reference


//template <typename T> 
TCMData TCMSparseW::getMax(bool bAbs) const {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in getMax()");
  #endif
  TCMData max, tMax;
  // start with unstored value for simplicity
  if(bAbs)  max = abs(unstoredValue);
  else      max = unstoredValue;
  // now iterate over the Array1D data using the optimized routines
  for(unsigned i=0; i<Cols; i++) {
    tMax = data[i].getMax(bAbs); // use Array1D optimized routine
    if(tMax>max)  max = tMax;
  } // end for i
  return max;
}; // end getMax


//template <typename T> 
TCMData TCMSparseW::getMin(bool bAbs) const {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in getMin()");
  #endif
  TCMData min, tMin;
  // start with unstored value for simplicity
  if(bAbs)  min = abs(unstoredValue);
  else      min = unstoredValue;
  // now iterate over the Array1D data using the optimized routines
  for(unsigned i=0; i<Cols; i++) {
    tMin = data[i].getMin(bAbs); // use Array1D optimized routine
    if(tMin<min)  min = tMin;
  } // end for i
  return min;
}; // end getMin

/*
void TCMSparseW::initEven(unsigned N, int q, bool bRandomNodes) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initEven()");
  #endif
  //  N - is the total size of the system
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  int iS, jS;

  msg("Assigning even communities... ");
  TVectorInt s(N);
  s.initStep(0,1,bRandomNodes);
    
  cout << red  << "Here c with kij pointer value as " << (int)kij << "\n" 
       << normText << endl; // debugging
  delete[] kij;
  cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[N];
  cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);

  cout << red << "Here e "  << endl;

  // loop over all pairs of nodes
  for(unsigned i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    for(unsigned j=i+1; j<N; j++) {
      jS = s[j];
      if((i/n0)==(j/n0)) {
        if(randomDouble()<p3) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n1)==(j/n1)) {
        if(randomDouble()<p2) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n2)==(j/n2)) {
        if(randomDouble()<p1) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else {
        if(randomDouble()<p0) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;} 
      } // end else
    } // end for j
  } // end for i
  
  //cout << redText << "Here f " << endl;
  copyEdges(edges);
  //cout << redText << "Here g " << endl;

  msg("done\n");
  return;
}; // end initEven
*/

int TCMSparseW::copyEdges(TListArray<TCMDataNode> &edges, int checkDuplicates,
                           bool bIgnoreThisParameter) {
  // copy a TListArray to the kij edge list and update class variables
  // ignored parameter is due to type conflicts for certain TCMDataNode types
  // returns the number of erased edges due to duplicates
  //cout << red << "Here 0"  << endl;
  int N = edges.getSize();
  //cout << red << "Here 1 "  << endl;

  destroykij();
  //cout << red << "data = " << (int)data << endl;  // debugging
  delete[] data;
  unstoredValue = -1;  // unweighted missing edge
  //cout << red << "Here 2 "  << endl;
  createkij(N);
  data = new Array1D<TCMData>[N];
  //cout << red << "Here 3 "  << endl;

  //msg("copying to sparse matrix data structure... ");
  nEdges = 0;
  // we copy in reverse order so that we can return memory as data is copied
  // this prevents us from having to allocate 2x the memory for init stage
  // (the erase routine swaps pointers for optimized speed, but even if constant
  // order is enforced, it needlessly adds O(q^2) pointer copies which competes
  // for the dominant order of the hierarchical algorithm)
  //cout << endl;  // debugging
  int tMod = max(1,N/50);  // user update interval
  int eraseCount = 0;  // track how many duplicates were deleted (or added)
 
  //cout << red << "Here 4 "  << endl;

  localInitTime();  // debugging
  TList<TCMDataNode> *pk;  // easy access to current TList
  int tMax = -1;
  for(int k=N-1; k>=0; k--) {
    //if((k%tMod)==0)  timeMsg("before sort "," step "+itos(k)+"...\n");
    //if((k%tMod)==0)  localInitTime();  // debugging
    pk = &(edges[k]);
    if(pk->getSize()>tMax)  tMax = pk->getSize();
/*
    if((k%tMod)==0)  cout << magenta << "\nCurrent max size = " <<  tMax << endl;
    //if((k%tMod)==0) cout << "." << flush;  // update user on progress
    if((k%tMod)==0) {
      cout << ".(" << k << ", size = " <<  edges[k].getSize()
           << ")\n";  // update user on progress
      cout << brown << "Before edges for node " << k << ":  " << edges[k] 
           << endl; // debugging
    }
*/
    //cout << "\nSorting list k = " << k << endl << edges[k] << flush; // debugging
    pk->sort();   // sort each list - required for binary search in kij
                  // sort is O(N) if already in order or is reverse ordered
/*
    if((k%tMod)==0)  timeMsg("sort "," step "+itos(k)+"...\n");
    if((k%tMod)==0) {
      cout << brown << "After edges for node " << k << ":  " << edges[k] 
           << endl; // debugging
    }
*/
    //cout << " done with sort... " << endl; // debugging

    if(checkDuplicates==1) { // check for duplicate entries end delete them
      // this uses the 'current' pointer logic for the TList class
      // there is a *very* difficult to find bug if I attempt to use the getNext
      // proceedure logic that is commented out below.  If the deleted node 
      // appears in the first position (very repeatable), there is a fault.  
      // I avoid this by deleting the second copy instead.  
      // This is only a workaround because the memory corruption could be 
      // occuring from somewhere else.
      //cout << "  copy with remove edges" << endl; // debugging

      // track previous rather than next node due to bug mentioned above
      TCMDataNode        *pn;
      pk->begin();
      pn = &(pk->getCurrent());  // set initial 'previous' since we start at 2
      // we start with the second node, this is fine even for a duplicate node 0
      pk->next();
                  
      int i = 0, nTests = pk->getSize()-1;  // getSize() can change
      //cout << "Starting Loop i (size " << (nTests+1) << "):  " << flush; // debugging
      while(i<nTests) {
        //cout << brown << i << " (" << flush << pk->getCurrent().index << ") " << flush;
        
        if(pk->getCurrent().index == pn->index) {
          pk->eraseCurrent();
          eraseCount++;
        } // end if
        // store current as new previous before iterating unless we deleted 
        // current, then iPrevious is unchanged
        else  pn = &(pk->getCurrent());  // update previous TCMDataNode

        pk->next();
        i++;
      } // end while i
      //cout << "done with step " << (i-1) << endl;  // debugging
    } // end if check duplicates
    else if(checkDuplicates>=2) { // check for duplicate entries and *add* them
      // this uses the 'current' pointer logic for the TList class
      // there is a *very* difficult to find bug if I attempt to use the getNext
      // proceedure logic that is commented out below.  If the deleted node 
      // appears in the first position (very repeatable), there is a fault.  
      // I avoid this by deleting the second copy instead.  
      // This is only a workaround because the memory corruption could be 
      // occuring from somewhere else.
      //cout << "  copy with remove edges" << endl; // debugging
      //if((k%tMod)==0) cout << "." << flush;  // update user on progress

      // track previous rather than next node due to bug mentioned above
      TCMDataNode        *pn;
      pk->begin();
      pn = &(pk->getCurrent());  // set initial 'previous' since we start at 2
      // we start with the second node, this is fine even for a duplicate node 0
      pk->next();
                  
      int i = 0, nTests = pk->getSize()-1;  // getSize() can change
      //cout << "Starting Loop i (size " << (nTests+1) << "):  " << flush; // debugging
      while(i<nTests) {
        //cout << brown << i << " (" << flush << pk->getCurrent().index << " and "
        //     << pn->index << ") " << flush;
        
        if(pk->getCurrent().index == pn->index) {
          //cout << red << " adding weights on indices " << pk->getCurrent().index
          //     << " and " << pn->index << brown << endl;  // debugging
          pn->value += pk->getCurrent().value;  // add the weights on 2
          pk->eraseCurrent();
          eraseCount++;
        } // end if
        // store current as new previous before iterating unless we deleted 
        // current, then iPrevious is unchanged
        else  pn = &(pk->getCurrent());  // update previous TCMDataNode

        pk->next();
        i++;
      } // end while i
      //cout << "done with step " << (i-1) << endl;  // debugging
    } // end if check duplicates
/*
    if(checkDuplicates==1) { // check for duplicate entries end delete them
      // this uses the 'current' pointer logic for the TList class
      TList<TCMDataNode> *pk = &(edges[k]);
      //cout << "Start duplicate test at k = " << k << "... " << flush; // debugging
      pk->begin();
      //cout << "begin with starting value = " << pk->getStart() 
      //     << "... " << endl; // debugging
      int i = 0;
      while(i<pk->getSize()-1) {
        cout << "Loop i=" << i << ":  " << flush << pk->getCurrent() << " and " 
             << flush << pk->getNext() << " ---> " << flush; // debugging
        if((pk->getCurrent()) == (pk->getNext())) {
          //cout << magenta << "Erasing edge i = " << flush << i << " with a value of " 
          //     << flush << pk->getCurrent().value << flush << "... " << endl; // debugging
          pk->eraseCurrent();
          eraseCount++;
          //cout << magenta << "after Erasing edge i = " << flush << i << " with a value of " 
          //     << flush << pk->getCurrent().value << flush << "... " << flush; // debugging
          //cout << "done." << endl;  // debugging
          //cout << "done." << endl << edges[k] << endl;  // debugging
          i--;  // decrement counter
        } // end if
        //cout << "ne" << flush; // debugging
        pk->next();
        //cout << "xt... " << flush << " with the next size as " << pk->getSize()
        //     << " " << endl; // debugging
        i++;
      } // end while i
    } // end if check duplicates
    else if(checkDuplicates>=2) { // check for duplicate entries and *add* them
      // this uses the 'current' pointer logic for the TList class
      TList<TCMDataNode> *pk = &(edges[k]);
      //cout << "Start duplicate test at k = " << k << "... " << flush; // debugging
      pk->begin();
      //cout << "begin with starting value = " << pk->getStart() 
      //     << "... " << endl; // debugging
      int i = 0;
      while(i<pk->getSize()-1) {
        //cout << "Loop i=" << i << ":  " << flush << pk->getCurrent() << " and " 
        //     << flush << pk->getNext() << " ---> " << flush; // debugging
        if((pk->getCurrent()) == (pk->getNext())) {
          //cout << "Erasing edge i = " << flush << i << " with a value of " 
          //     << flush << pk->getCurrent() << flush << "... " << flush; // debugging
          pk->getNext().value += pk->getCurrent().value;
          pk->eraseCurrent();
          //cout << "done." << endl << edges[k] << endl;  // debugging
          i--;  // decrement counter
        } // end if
        //cout << "ne" << flush; // debugging
        pk->next();
        //cout << "xt... " << flush << " with the next size as " << pk->getSize()
        //     << " " << endl; // debugging
        i++;
      } // end while i
    } // end if check duplicates
*/
    
    //cout << "Creating data array of size (" << edges[k].getSize() 
    //     << ")... " << flush;  // debugging
    kij[k].resize(pk->getSize());
    //cout << " kij " << flush;
    data[k].resize(pk->getSize());
    //if((k%tMod)==0)  timeMsg("resize "," step "+itos(k)+"... ");
    //cout << " data " << flush;
    pk->begin();
    for(int j=0; j<edges[k].getSize(); j++) {
      kij[k][j]  = pk->getCurrent().index;
      data[k][j] = pk->getCurrent().value;
      //kij[k][j]  = edges[k][j].index;
      //data[k][j] = edges[k][j].value;
      //cout << "Copying edge connecting nodes " << k << " to " << edges[k][j].index 
      //     << " with weight of " << edges[k][j].value << " or " << k << " to " 
      //     << kij[k][j] << " with weight of " << data[k][j] << "\n"; // debugging
      pk->next();
    } // end for j
    //cout << blue << "kij[" << k << "]:  " << normText << data[k] << endl;  // debugging
    //if((k%tMod)==0)  timeMsg("copy "," step "+itos(k)+"... ");
    nEdges += kij[k].getSize();
    edges.erase(k);    // now return the memory of this TList
    //if((k%tMod)==0)  timeMsg("erase "," step "+itos(k)+"... ");
    //cout << "done with data array for k = " << k << endl;  // debugging
    //cout << brown << " exiting " << endl;
  } // end for k
  
  nEdges /= 2;  // correct for double counted edges
  nNodes = N;  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1)); // calculate actual density
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  cout << "done." << endl;
  #endif

  return eraseCount;
}; // end copyEdges with TCMDataNode


int TCMSparseW::copyEdges(TListArray<int> &edges, int checkDuplicates) {
  // copy a TListArray to the kij edge list and update class variables
  // returns the number of erased edges due to duplicates
  //cout << red << "Here 1 "  << endl;
  if(checkDuplicates>=2)  
    errorMsg("Cannot add edge weights with copyEdges with TListArray<int>");

  int N = edges.getSize();
  //cout << red << "Here 1 "  << endl;
  destroykij();
  //cout << red << "Here 1 "  << endl;
  createkij(N);
  //cout << red << "Here 1 "  << endl;
  unstoredValue = -1;  // unweighted missing edge

  //msg("copying to sparse matrix data structure... ");
  nEdges = 0;
  // we copy in reverse order so that we can return memory as data is copied
  // this prevents us from having to allocate 2x the memory for init stage
  // (the erase routine swaps pointers for optimized speed, but even if constant
  // order is enforced, it needlessly adds O(q^2) pointer copies which competes
  // for the dominant order of the hierarchical algorithm)
  //cout << endl;  // debugging
  int tMod = max(1,N/50);  // user update interval
  int eraseCount = 0;  // track how many duplicates were deleted (or added)

  for(int k=N-1; k>=0; k--) {
    //cout << brown << "Nodes " << k << ":  " << green << edges[k] << endl; // debugging
    //cout << "\nSorting list k = " << k << endl << edges[k] << flush; // debugging
    //cout << k  << flush;
    //cout << red << "\nStart test at k = " << k << " with " << edges[k] 
    //     << "... " << endl; // debugging
    edges[k].sort();   // sort each list - required for binary sort in kij
                       // sort is O(N) if already in order or is reverse ordered
    //cout << " done with sort..." << flush << edges[k] << endl; // debugging

    if(checkDuplicates==1) { // check for duplicate entries end delete them
      // this uses the 'current' pointer logic for the TList class
      // there is a *very* difficult to find bug if I attempt to use the getNext
      // proceedure logic that is commented out below.  If the deleted node 
      // appears in the first position (very repeatable), there is a fault.  
      // I avoid this by deleting the second copy instead.  
      // This is only a workaround because the memory corruption could be 
      // occuring from somewhere else.
      if((k%tMod)==0) cout << "." << flush;  // update user on progress

      TList<int> *pk = &(edges[k]);
      //pk->begin();
      edges[k].begin();
      int iPrevious = pk->getCurrent();
      // we start with the second node, this is fine even for a duplicate node 0
      pk->next();  
                  
      //cout << magenta << "begin with starting value = " << pk->getStart()
      //     << "... " << endl; // debugging
      int i = 0, nTests = pk->getSize()-1;  // getSize() can change
      //while(i<pk->getSize()-1) {
      while(i<nTests) {
        //cout << brown << "Loop i=" << i << ":  " << flush << pk->getCurrent()
        //     << " and " << flush;
        //if(i<nTests-1)
        //  cout << pk->getNext() << " ---> " << flush; // debugging
        
        if(edges[k].getCurrent() == iPrevious) {
        //if(edges[k].getCurrent() == edges[k].getNext()) {
        //if(pk->getCurrent() == pk->getNext()) {  // very difficult to find bug
          //cout << normText << "Erasing edge i = " << flush << i << " with a value of "
          //     << flush << pk->getCurrent() << " for k = " << k << endl;
          //cout << grey << edges[k-1] << flush; // debugging
          //cout << edges[k]   << flush; // debugging
          edges[k].eraseCurrent();
          eraseCount++;  // track how many duplicates were deleted (or added)
          //cout << "done with erase." << endl;
          //pk->eraseCurrent();
          //cout << flush << magenta << "after erasing edge i = " << flush << i
          //     << flush << " at " << pk->getCurrent() << flush;
              // << " and " << pk->getNext() << flush 
              // << "... " << flush; // debugging
          //cout << "done." << endl << edges[k] << endl;  // debugging
          //i--;  // decrement counter for the deleted node - needed for Next version
        } // end if
        // store current before iterating
        // unless we deleted current, then iPrevious is unchanged
        else  iPrevious = pk->getCurrent();  

        //cout << "ne" << flush; // debugging
        edges[k].next();
        //pk->next();
        //cout << "xt... " << flush << " with the next size as " 
        //     << pk->getSize() << " " << endl; // debugging
        i++;
      } // end while i
    } // end if check duplicates

    kij[k]  = edges[k]; // use overloaded '='
    nEdges += kij[k].getSize();
    // return memory of TList k - fast doesn't matter since we delete at the end
    edges.erase(k,1);
    //cout << brown << " exiting " << endl;
  } // end for k
  
  nEdges /= 2;  // correct for double counting edges
  nNodes = N;  Cols = N;
  p = (double)nEdges*2.0/((double)N*(double)(N-1)); // calculate actual density
  //cout << red << "Ok, setting valid array flag for SparseW with nEdges = " 
  //     << nEdges << " and a density of " << p << " and N = " << N << endl;
  //errorMsg("Done");
  #ifdef DEBUG_MODE
  msg("done");
  //cout << red << "Ok, setting valid array flag for SparseW." << green << endl;
  ValidArray = ValidCMatrixFlag;
  #endif
  //cout << "done." << green << endl;
  return eraseCount;
}; // end copyEdges with int


void TCMSparseW::initHeterogeneous(TClusterList &c, TClusterList &merged,
       int N, Array1D<int> &v, double p2, double p1, double p0, 
       TCMData bottomWeight, TCMData topWeight) {
  // initialize a heterogeneous hierarchy but with equal density propabilities
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  // v - specifies which clusters to merge.  i.e. if v = {1,3,2}, then 
  //     leave cluster 1 alone (the 1); merge clusters 2, 3, and 4 (the 3);
  //     and merge clusters 5 and 6 (the 2)
  // bRandomNodes - specifies whether the nodes will the randomized or 
  //                sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is invalid in initHeterogenous()");
  #endif
  if(v.getSum()!=c.getSize())  
    errorMsg("Merge vector sum != c.nClusters in initHeterogenous()");
    
  int iS, jS;
  unstoredValue = -1;  // unweighted missing edge

  //cout << "Just before assign... " << flush;
  msg("Assigning heterogeneous hierachical communities... ");
    
  //cout << red  << "Here c with kij pointer value as " << (int)kij 
  //     << "\n" << endl; // debugging
  //delete[] kij;
  //cout << red << "Here c2 " << endl;
  //kij = new Array1D<int>[N];
  //cout << red << "Here d "  << endl;
  //cout << "after delete with N = " << N << "... " << flush;  // debugging

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);

  //cout << red << "Here e "  << endl;                             // debugging

  // work with super-clusters first for level 1 communities
  //TClusterList merged;
  merged = c;
  // merge original clusters by cluster merge list v.  We merge backwards in
  // order to not change the individual cluster order during the move(I,To_J) 
  // function calls although the ending merged order is actually mixed due to 
  // the move call speed optimizations.  Since we use merged internally here, 
  // this should not be a problem.
  //cout << "before merges with pre-merge as:\n"; // << merged << flush; // debugging
  int tSize, bDisplay = 0;
  cout << endl;
  for(int j=0; j<v.getSize(); j++) { // loop over all merges
/*
    //cout << magenta << j << ":  " << normText << flush;          // debugging
    tSize = 0;  // crude detection for special communities
    for(int k=0; k<v[v.getSize()-1-j]; k++) { // only merge if v[j]>1
      //cout << red << k << " " << normText << flush;              // debugging
      tSize += merged[merged.getSize()-j-1-k].getSize();
    } // end for k
    if(tSize==27 || tSize==180)  bDisplay = 1;
    if(bDisplay) {
      //cout << magenta << "Size " << tSize << " merged community is:\n";
      for(int k=0; k<v[v.getSize()-1-j]; k++) { // only merge if v[j]>1
        //cout << red << k << " " << normText << flush;            // debugging
        merged[merged.getSize()-j-1-k].display();
        //cout << "\n";
      } // end for k
      //cout << "Stepwise generation of the community:\n";
    } // end if bDisplay
*/
    for(int k=0; k<v[v.getSize()-1-j]-1; k++) { // only merge if v[j]>1
/*
      //cout << red << k << " " << normText << flush;              // debugging
      if(bDisplay) {
        merged[merged.getSize()-j-1].display();
        cout << endl;
      } // end if bDisplay
*/
      merged.move(merged.getSize()-j-1,merged.getSize()-j-2,1);
    } // end for k
/*
    if(bDisplay) {
      cout << normText << endl;                                    // debugging
      bDisplay = 0;
    } // end if bDisplay
*/
  } // end for j
/*
  tSize = 0;
  for(int i=0; i<c.getSize(); i++)  tSize += c[i].getSize();
  cout << "Average original list size is:" << (double)tSize/(double)c.getSize() 
       << "\n";
  tSize = 0;
  for(int i=0; i<merged.getSize(); i++)  tSize += merged[i].getSize();
  cout << "Average merged list size is:" << (double)tSize/(double)merged.getSize() 
       << "\n";
  //errorMsg("Done");  // debugging
*/
  //merged.display(0,0,"merged inside ");                          // debugging

  //int nNodes0 = 0, nNodes1 = 0, nNodes2 = 0;  // count number of nodes at each level
  int nEdges0 = 0, nEdges1 = 0, nEdges2 = 0, nodeEdgeCount;  // count number of nodes at each level
  double nodeDensity;
  int NEdgesMax0 = 0, NEdgesMax1 = 0, NEdgesMax2 = 0;  // count max number of nodes at each level
  // fill matrix - this should eliminate duplicate edges and make the algorithm
  // fill with the density as advertized (instead of Level 2 being p2, 
  // it was (1 - (1-p1)*(1-p2)) or about 0.93 for p1 = 0.3 and p2 = 0.9).
  int icStart = 0;  // start of which clusters are we looping over
  for(int n=0; n<v.getSize(); n++) {
    for(int m=icStart; m<icStart+v[n]; m++) {
      //cout << "Looking at cluster size " << c[m].getSize() << " with: ";
      for(int k=m+1; k<icStart+v[n]; k++) {
        NEdgesMax1 += c[m].getSize()*c[k].getSize();
        //cout << c[k].getSize() << " ";

        for(int i=0; i<c[m].getSize(); i++) {
          iS = c[m].getNode(i);
        
          // loop over all pairs of nodes in cluster k
          for(int j=0; j<c[k].getSize(); j++) {
            jS = c[k].getNode(j);
            if(randomDouble()<p1) { 
              edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdges1++; }
          } // end for j
        } // end for i
      
      } // end for n
      //cout << endl; // debugging
    } // end for k
    
    icStart += v[n];  // step up to the next set of sub-clusters based on v[]
  } // end for m
  
  //errorMsg("Stop here!");
  cout << "Initializing level 2... " << flush;  // debugging
  // now work with original clusters for level 2 most dense communities
  for(int k=0; k<c.getSize(); k++) {
    NEdgesMax2 += c[k].getSize()*(c[k].getSize()-1)/2;
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iS = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jS = c[k].getNode(j);
        // this could add some duplicate links which are removed during the copy
        // stage below.  Otherwise, we have to do an O(k) search for each link
        // added to check whether it is a duplicate.
        if(randomDouble()<p2) { 
          edges[iS].add(jS);  edges[jS].add(iS);
          nEdges++;  nEdges2++;
        } // end if randomDouble
      } // end for j
    } // end for i
  } // end for k
  
  // now fill the rest of the system super-cluster to super-cluster at p0
  for(int k=0; k<merged.getSize(); k++) {
    // loop over all other clusters m (past k) for outside connections
    for(int m=k+1; m<merged.getSize(); m++) {
      NEdgesMax0 += merged[k].getSize()*merged[m].getSize();

      // loop over all nodes in cluster k
      for(int i=0; i<merged[k].getSize(); i++) {
        iS = merged[k].getNode(i);

        // loop over nodes j in cluster m
        for(int j=0; j<merged[m].getSize(); j++) {
          jS = merged[m].getNode(j);
          if(randomDouble()<p0) {
            edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdges0++; }
        } // end for j
      } // end for i
    } // end for m
  } // end for k
  
  for(int i=0; i<edges.getSize(); i++)  edges[i].sort();
  //cout << "Just before copy... " << flush;
  copyEdges(edges,1);  // copy and erase duplicate entries from step p2 above
  //cout << "done with copy" << endl;
  initRandomWeights(bottomWeight,topWeight);
/*
  // consistency check on the number of edges for all level 2 clusters
  for(int k=0; k<c.getSize(); k++) {
    for(int i=0; i<c[k].getSize(); i++) {
      iS = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      nodeEdgeCount = 0;  // reset the edge count for this node
      for(int j=0; j<c[k].getSize(); j++) {
        jS = c[k].getNode(j);
        if((*this)(iS,jS)>0)  nodeEdgeCount++;
      } // end for j

      nodeDensity = (double)nodeEdgeCount/(double)c[k].getSize();
      if(0.75*p2>nodeDensity)
        warningMsg("Node "+itos(iS)+" only has a level 3 density of "
                   +ftos(nodeDensity));

    } // end for i
  } // end for k
*/
  // write some initialization information to the matrix description string
  description += 
      "\"Initialized to a size " + itos(N) + " heterogeneous hierarchy:\","
      + "\"p1 = " + dtos(p0) + "\"," + "\"p2 = " + dtos(p1) + "\","
      + "\"p3 = " + dtos(p2)
      + "\"\n\"Actual densities:\","
      + "\"p1a = " + dtos((double)nEdges0/(double)NEdgesMax0) + "\"," 
      + "\"p2a = " + dtos((double)nEdges1/(double)NEdgesMax1) + "\","
      + "\"p3a = " + dtos((double)nEdges2/(double)NEdgesMax2)
      + "\"\n\"New interior edges at each level (no sub-cluster edges):\","
      + "\"L1 = " + itos(nEdges0,1) + "\"," + "\"L2 = " + itos(nEdges1,1) +"\","
      + "\"L3 = " + itos(nEdges2,1) + ",\"Total = " + itos(nEdges,1)
      + "\"\n\"Maximum new interior edges at each level:\","
      + "\"M1 = " + itos(NEdgesMax0) + "\"," 
      + "\"M2 = " + itos(NEdgesMax1) + "\","
      + "\"M3 = " + itos(NEdgesMax2)
      + "\"\n\"Number of clusters:\","
      + "\"q1 = " + itos(1) + "\"," + "\"q2 = " + itos(merged.getSize()) + "\","
      + "\"q3 = " + itos(c.getSize())
      + "\"\n\"Number of q3 clusters merged to make 2:\"";
  if(v.getSize()<901)  // approximate column limitation in OpenOffice
    for(int i=0; i<v.getSize(); i++)  description += "," + itos(v[i]);
  description += "\n\"Merged cluster size list for level 2 (Total " 
               + itos(merged.getSize()) + "):  \"";
  if(merged.getSize()<901)  // approximate column limitation in OpenOffice
    for(int i=0; i<merged.getSize(); i++)  
      description += "," + itos(merged[i].getSize());
  description += "\n\"Original cluster size list for level 3 (Total " 
               + itos(c.getSize()) + "):\"";
  if(c.getSize()<901)  // approximate column limitation in OpenOffice
    for(int i=0; i<c.getSize(); i++)  
      description += "," + itos(c[i].getSize());
  description += "\n";  // finish description information

    // output the maximum and minimum sizes of the constructed lists
    int  minClusterSize = BigInteger, maxClusterSize = -BigInteger;  
    // work on the merged list
    for(int i=0; i<c.getSize(); i++)  {
      tSize = c[i].getSize();
      if(tSize<minClusterSize)  minClusterSize = tSize;
      if(tSize>maxClusterSize)  maxClusterSize = tSize;
    } // end for i
    description += "\"max/min cluster sizes for level 3 are:\"," 
                   +itos(maxClusterSize)+","+itos(minClusterSize)+"\n";
    // work on the testing list
    minClusterSize = BigInteger;  maxClusterSize = -BigInteger;  
    for(int i=0; i<merged.getSize(); i++)  {
      tSize = merged[i].getSize();
      if(tSize<minClusterSize)  minClusterSize = tSize;
      if(tSize>maxClusterSize)  maxClusterSize = tSize;
      // crude identification of some communities for the paper
      if(tSize==27)   merged[i].display();
      if(tSize==180)  merged[i].display();
    } // end for i
    description += "\"max/min cluster sizes for level 2 are:\"," 
                   +itos(maxClusterSize)+","+itos(minClusterSize)+"\n";

  msg("done\n");
  return;
}; // end initHeterogeneous


void TCMSparseW::initRandomWeights(TCMData bottomWeight, TCMData topWeight) {
  // with an already defined kij, initialize the weights randomly
  // now define random weights for the random edges
  TCMData range = topWeight - bottomWeight;
  //cout << "Deleting data array... " << flush;
  //cout << red << "data = " << (int)data << endl;  // debugging
  delete[] data;  
  unstoredValue = -1;  // unweighted missing edge
  //cout << "Allocating new data array... " << flush;
  data = new Array1D<TCMData>[Cols];
  //cout << "Filling data array... " << flush;
  for(unsigned i=0; i<Cols; i++) {
    data[i].resize(kij[i].getSize());
    for(int j=0; j<kij[i].getSize(); j++)
      data[i][j] = (TCMData)randomDouble((double)range) + bottomWeight;
    //cout << magenta << "i = " << i << ":  " << blue << kij[i] << endl;  // debugging
  } // end for i
  //cout << "Done." << endl;
  return;
}; // end copyEdges with int


//template <typename T> 
void TCMSparseW::initByClusters2(TClusterList &c, int N, double pin, double pout,
     int &nEdgesIn, int &nEdgesOut, TCMData  bottomW, TCMData topW) {
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid()) errorMsg("TCMSparseW object is not valid in initByClusters()");
  #endif
  if(c.getNNodes()==0)  
    errorMsg("In initByClusters(), cluster list says it has zero nodes?");

  //msg("Assigning communities by a specified set of clusters... ");
    
  //destroykij();
  //createkij(N);

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);
  //cout << red << "Here e with N = " << N << endl;

  //c.display();  // debugging

  int iS, jS;
  nEdgesIn = 0, nEdgesOut = 0;
  for(int k=0; k<c.getSize(); k++) {
    // loop over all nodes in cluster k
    c[k].begin();
    //cout << red << " i = " << k << flush;  // debugging

    for(int i=0; i<c[k].getSize(); i++) {
      //cout << red << " j = " << i << flush;  // debugging
      iS = c[k].getCurrentNode();
      //iS = c[k].getNode(i);  // redunant debugging
      c[k].getNode(i);  // redunant debugging - bug with op[] below - workaround

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        //cout << red << " j = " << i << flush;  // debugging
        jS = c[k].getNode(j);
        if(randomDouble()<pin) { 
          edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdgesIn++; }
      } // end for j

      // loop over all other clusters m (past k) for outside connections
      for(int m=k+1; m<c.getSize(); m++) {
        // loop over nodes j in cluster m
        for(int j=0; j<c[m].getSize(); j++) {
          jS = c[m].getNode(j);
          if(randomDouble()<pout) { 
            edges[iS].add(jS); edges[jS].add(iS); nEdges++;  nEdgesOut++; }
        } // end for j
      } // end for m

    c[k].next();
    } // end for i
  } // end for k
  
  //cout << "k = " << k << ":  " << edges[k] << endl;  // debugging

  copyEdges(edges);
  //cout << green << "\nInitializing weights with flag = " << ValidArray << endl;
  initRandomWeights(bottomW,topW);

  //msg("done\n");
  return;
}; // end initByClusters2
void TCMSparseW::initByClusters(TClusterList &c, int N, double pin, double pout,
     TCMData  bottomWeight, TCMData topWeight) {
  // a version that does not require the reference edge return values
  int ignore1, ignore2;
  initByClusters2(c,N,pin,pout,ignore1,ignore2,bottomWeight,topWeight); 
}; // end non-reference version


//template <typename T> 
void TCMSparseW::initByClustersSparse(TClusterList &c, double pin, double pout,
     TCMData  bottomWeight, TCMData topWeight) {
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initByClusters()");
  #endif
  if(c.getNNodes()==0)  
    errorMsg("In initByClusters(), cluster list says it has zero nodes?");
  int iS, jS;
  int N = c.getNNodes();

  //msg("Assigning communities by a specified set of clusters... ");
    
  //cout << red  << "Here c " << endl; // debugging
  delete[] kij;
  //cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[N];
  //cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);
  //cout << red << "Here e with N = " << N << endl;

  for(int k=0; k<c.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iS = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jS = c[k].getNode(j);
        if(randomDouble()<pin) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } // end for j

    } // end for i
  } // end for k
  
  //cout << blue << "No. of edges before random additions = " << nEdges;

  // now randomly select nodes outside between nodes each mutually outside of
  // their own clusters.  Add the edge, but we rely on the copy edges routine
  // to catch duplicates.  We do this by looping over all nodes and adding an
  // evenly distributed number of edges per node.
  double nExteriorPossibleEdges = 0.0; // use double because of likely overflow
  // loop over all cluster pairs, calculate the total number of possible edges
  for(int i=0; i<c.getSize(); i++)  for(int j=i+1; j<c.getSize(); j++)
    nExteriorPossibleEdges += (double)c[i].getSize()*(double)c[j].getSize();

  // calculate the target number of exterior edges, here, we are relying on the 
  // sparse structure to limit interior and exteriror duplicates since it is 
  // time consuming and difficult to catch and eliminate them.
  // nEdges is invalid after this; however, it will be updated in copyEdges()
  int nEdgesTarget = (int)( pout*nExteriorPossibleEdges );
  //edges.sort();  // need a sorted list to count double edges in logic below
  //cout << "Before random additions:\n" << edges << endl;
  //copyEdges(edges); // copy edges and eliminate duplicates // debugging
  //initRandomWeights(bottomWeight,topWeight);  // debugging
  //display(); // debugging
  //cout << "Back in init " << endl;

  for(int i=0; i<nEdgesTarget; i++) {
    iS = randomInt(0,N-1);
    jS = randomInt(0,N-1);
    while(iS==jS)  jS = randomInt(0,N-1);        // no self-connections allowed
    //cout << "iS = " << iS << "\tjS = " << jS << "   ";           // debugging
    edges[iS].add(jS);  edges[jS].add(iS);
  } // end for i
  nEdges += nEdgesTarget;  // assuming no duplicates at the moment
  //cout << "After random additions:\n" << edges << endl;

  // count number of duplicates that were introduced by the random assignments
  int dCount = 0;
  TList<int> *pk;
  for(int k=0; k<edges.getSize(); k++) {
    edges[k].sort();  // need a sorted list to count double edges in logic below
    pk = &(edges[k]);
    pk->begin();

    for(int i=0; i<pk->getSize()-1; i++) {
      if((pk->getCurrent()) == (pk->getNext()))  dCount++;
      pk->next();
    } // end for i
  } // end for k
  //cout << "After duplicate count:\n" << edges << endl;

  //cout << blue << "\nTarget No. of exterior edges = " << nEdgesTarget << " and "
  //     << "No. of possible edges = " << nExteriorPossibleEdges
  //     << magenta << " (with " << dCount << " duplicates)" << normText << "\n";
  //cout << blue << "No. of edges after random additions = " << nEdges;

  //cout << "Here f " << endl;
  copyEdges(edges,1); // copy edges and eliminate duplicates
  //cout << green << "Initializing weights " << endl;
  initRandomWeights(bottomWeight,topWeight);
  //cout << blue << "Final L = " << nEdges;

  msg("done\n");
  return;
}; // end initByClustersSparse


void TCMSparseW::initByClustersRS(TClusterList &c, unsigned N, 
           double pin, double pout, unsigned long long ISeed, 
           TCMData  bottomWeight, TCMData topWeight) {
  //  N - is the total size of the system
  //cout << red << "here a in initRandom " << normText << flush;   // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initRandom()");
  #endif
  nNodes = N;
  //msg("Assigning random communities...");
  cout << "Assigning random communities... " << flush;
  cout << " with a seed of " << rngSeedString() << endl; 
  
  destroykij();  // just in case - still valid if NULL

  //cout << red << "here c in initRandom " << normText << flush;   // debugging
  // use O(L) version
  // tEdges is the number of expected edges
  int tEdges = (int)( pout*(double)N*(double)(N-1)/2.0 + 0.5 );
  //cout << red << "Number of random edges expected is " << tEdges << endl;
  int iN, jN, tMod = tEdges/50;
  nEdges = 0;
  TVectorInt  degrees(N,0);  // store a count of the number of edges

  rngReseed(ISeed);  // start with initial seed here and just count edges
  //r.Reseed(ISeed);   // debugging
  //rngReseed(rngSeedString());  // start with initial seed here and just count edges
  //cout << red << "first 4 random numbers are "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " " << endl;
    
  double maxcCount = 0.0, cCount, totalCount;
  for(int k=0; k<c.getSize(); k++) {
    maxcCount += (double)c[k].getSize()*(double)(c[k].getSize()-1)/2.0;
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iN = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jN = c[k].getNode(j);
        if(randomDouble()<pin) { degrees[iN] += 1;  degrees[jN] += 1; }
      } // end for j
    } // end for i
  } // end for k
  
  cCount = maxcCount*pin;  totalCount = cCount + tEdges;
  //cout << "Number of cluster edges expected is " << cCount 
  //     << " for a total possible of " << totalCount << endl;

  for(unsigned i=0; i<tEdges; i++) {
    if((i%tMod)==0) cout << "." << flush;  // update user on progress
    iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
    // assumes symmetric
    while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
    // checking for duplicates at this point is expensive, but we should not
    // have many duplicates (for large systems), so increment nEdges and allow 
    // copyEdges to catch duplicates
    degrees[iN] += 1;  degrees[jN] += 1;
  } // end for i

  //cout << "Allocating degree data... " << flush;
  //cout << endl << degrees << endl << " with i = " << flush;
  // now create kij at the (almost, save for duplicates) correct sizes
  createkij(N);  // allocate the rows (columns)
  data = new TVectorInt[N];
  for(unsigned i=0; i<N; i++) {
    //cout << i << "," << flush;
    kij[i].resize((unsigned)degrees[i]);   // actual memory allocation
    //cout << " " << flush;
    kij[i].resize(0,0,0,1);    // fake resize to zero to prepare insertion sort
  } // end for i

  // I think this will benefit cache to allocate in a separate 'chunk'
  for(unsigned i=0; i<N; i++)  data[i].resize((unsigned)degrees[i]);  

  // now regenerate the entire list of random numbers to fill the list
  rngReseed(ISeed);  // start with initial seed here and just count edges
  //cout << "Filling degree data... " << flush;
  bool bInserted;
  for(int k=0; k<c.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iN = c[k].getNode(i);
      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jN = c[k].getNode(j);
        if(randomDouble()<pin) { // these are guaranteed to not be duplicated
          bInserted = insertNode(iN,jN,1);
          bInserted = insertNode(jN,iN,1);
          nEdges++;
        } // end if random
      } // end for j
    } // end for i
  } // end for k


  for(unsigned i=0; i<tEdges; i++) {
    if((i%tMod)==0) cout << "." << flush;  // update user on progress
    iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
    // assumes symmetric
    while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
    // now try to add the node to an existing list using insertion sort logic
    // duplicates are checked in the insertion (sort) where the parameter
    // is set to ignore duplicates (duplicate memory is not returned since it 
    // is assumed to be small relative to the size of kij)
    bInserted = insertNode(iN,jN,1);
    if(bInserted) {
      insertNode(jN,iN,1);  // assumes symmetric edges
      nEdges++;
    } // end if bInserted
  } // end for i
  cout << "\n";  // finish user updates

  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1));  // calculate actual density

  //initRandomWeights(bottomWeight,topWeight);
  for(unsigned i=0; i<N; i++)  data[i].init(1);

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
  msg("done.");
  return;
}; // end initByClustersRS


void TCMSparseW::initByClustersPowerRS(TClusterList &c, unsigned N, 
       double pin, double kMin, double kMax, double alpha, 
       int &nEdgesIn, int &nEdgesOut,
       unsigned long long ISeed, TCMData  bottomWeight, TCMData topWeight) {
  //  N - is the total size of the system
  //cout << red << "here a in initRandom " << normText << flush;   // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initRandom()");
  #endif
  nEdgesIn = 0;  nEdgesOut = 0;
  nNodes = N;
  #ifdef DEBUG_MODE
  //msg("Assigning random communities...");
  cout << "Assigning random communities with a seed of " << rngSeedString()
       << "\nBeginning kMin " << kMin << ") and kMax = " << kMax << endl;
  #endif
  
  destroykij();  // just in case - still valid if NULL

  //cout << red << "here c in initRandom " << normText << flush;   // debugging
  // use O(L) version
  // tEdges is the number of expected edges
  //int tEdges = (int)( pout*(double)N*(double)(N-1)/2.0 + 0.5 );
  int iN, jN, tMod = N/50;
  nEdges = 0;
  TVectorInt  degrees(N,0);  // store a count of the number of edges
  TVectorInt  kStored(N,0);  // store a count of the number of edges

  rngReseed(ISeed);  // start with initial seed here and just count edges
  //r.Reseed(ISeed);   // debugging
  //cout << red << "first 4 random numbers are "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " " << endl;

  #ifdef DEBUG_MODE
  cout << magenta << "Counting edges for degree distribution... " 
       << endl;  // debugging
  #endif
  int tEdges;
  double maxcCount = 0.0, cCount, totalCount;
  for(int k=0; k<c.getq(); k++) {
    //cout << cyan << "cluster " << k << blue << " (size = " << c[k].getn() << ")"
    //     << cyan << ":  " << green << flush;  // debugging
    maxcCount += c[k].getnD()*(c[k].getnD()-1.0)/2.0;
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getn(); i++) {
      iN = c[k].getNode(i);
      //cout << iN << " " << flush;  // debugging

      // assign internal edges
      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getn(); j++) {
        jN = c[k].getNode(j);
        if(randomDouble()<pin) { 
          degrees[iN] += 1;  degrees[jN] += 1;  nEdgesIn++; }
      } // end for j
      kStored[i] = degrees[iN];

      // assign external edges
      //cout << "." << flush;  // debugging
      //tEdges = (int)( randomPower(kMin,kMax,alpha) + 0.5 ) - kStored[i];
      tEdges = (int)( randomPower(kMin,kMax,alpha) + 0.5 ) - degrees[i];
      if(tEdges<0)  tEdges = 0;
      //cout << " (tEdges = " << tEdges << ") " << endl;  // debugging
      for(unsigned j=0; j<tEdges; j++) {
        //if((i%tMod)==0) cout << "." << flush;  // update user on progress
        iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
        // assumes symmetric
        while(iN==jN)  jN = randomInt(0,N-1); // try again on duplicate
        // checking for duplicates at this point is expensive, but we should 
        // not have many duplicates (for large systems), so increment nEdges 
        // and allow copyEdges to catch duplicates
        degrees[iN] += 1;  degrees[jN] += 1;
        nEdgesOut++;
      } // end for j
      //cout << ". " << flush;  // debugging
    } // end for i
    //cout << endl;  // debugging
  } // end for k

  #ifdef DEBUG_MODE
  cout << "Allocating degree data... " << flush;
  //cout << endl << degrees << endl << " with i = " << flush;
  #endif
  // now create kij at the (almost, save for duplicates) correct sizes
  createkij(N);  // allocate the rows (columns)
  data = new TVectorInt[N];
  for(unsigned i=0; i<N; i++) {
    kij[i].resize((unsigned)degrees[i]);   // actual memory allocation
    kij[i].resize(0,0,0,1);    // fake resize to zero to prepare insertion sort
  } // end for i

  // I think this will benefit cache to allocate in a separate 'chunk'
  for(unsigned i=0; i<N; i++)  data[i].resize((unsigned)degrees[i]);  

  // now regenerate the entire list of random numbers to fill the list
  rngReseed(ISeed);  // start with initial seed here and just count edges

  #ifdef DEBUG_MODE
  cout << "Filling degree data... " << flush;
  #endif
  bool bInserted;
  for(int k=0; k<c.getq(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getn(); i++) {
      iN = c[k].getNode(i);
      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getn(); j++) {
        jN = c[k].getNode(j);
        if(randomDouble()<pin) { // these are guaranteed to not be duplicated
          bInserted = insertNode(iN,jN,1);
          bInserted = insertNode(jN,iN,1);
          nEdges++;
        } // end if random
      } // end for j

      //tEdges = (int)( randomPower(kMin,kMax,alpha) + 0.5 ) - kStored[i];
      tEdges = (int)( randomPower(kMin,kMax,alpha) + 0.5 ) - degrees[i];
      if(tEdges<0)  tEdges = 0;
      for(unsigned j=0; j<tEdges; j++) {
        //if((i%tMod)==0) cout << "." << flush;  // update user on progress
        iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
        // assumes symmetric
        while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
        // now try to add the node to an existing list using insertion sort 
        // logic duplicates are checked in the insertion (sort) where the 
        // parameter is set to ignore duplicates (duplicate memory is not 
        // returned since it is assumed to be small relative to the size of kij)
        bInserted = insertNode(iN,jN,1);
        if(bInserted) {
          insertNode(jN,iN,1);  // assumes symmetric edges
          nEdges++;
        } // end if bInserted
      } // end for j
    } // end for i
  } // end for k
  #ifdef DEBUG_MODE
  cout << "\n";  // finish user updates
  #endif

  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1));  // calculate actual density

  //initRandomWeights(bottomWeight,topWeight);
  for(unsigned i=0; i<N; i++)  data[i].init(1);  // temporary

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  msg("done.");
  #endif
  return;
}; // end initByClustersPowerRS


//template <typename T> 
void TCMSparseW::initSHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
     double p0, double p1a, double p2a, double p3a, 
                double p1b, double p2b, double p3b, bool bRandomNodes,
     TCMData  bottomWeight, TCMData topWeight) {
  // initialize a 'staggered' hierarchy'.  That is, half of the hiearchy is 
  // initialized by p1a, p2a, and p3a for its levels and the other half is 
  // initialized by p1b, p2b, and p3b for its levels.
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.

  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initHierarchy()");
  #endif
  //  N - is the total size of the system
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  int iS, jS;

  msg("Assigning random staggered hierarchical communities... ");
  cout << red << "Here 1 " << endl;
  TVectorInt s(N);
  cout << red << "Here 2 " << endl;
  s.initStep(0,1,bRandomNodes);
  cout << red << "Here 3 " << endl;
    
  delete[] kij;
  cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[N];
  cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);

  //cout << red << "Here e "  << endl;

  // loop over all pairs of nodes
  for(unsigned i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    double pI, pII, pIII;
    for(unsigned j=i+1; j<N; j++) {
      jS = s[j];
      if((i/n3)==(j/n3)) {
        if(i/(N/2)==0) pIII = p3a;  else pIII = p3b;  // pick side of hiearchy
        if(randomDouble()<pIII){edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n2)==(j/n2)){
        if(i/(N/2)==0) pII = p2a;   else pII = p2b;   // pick side of hiearchy
        if(randomDouble()<pII) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n1)==(j/n1)){
        if(i/(N/2)==0) pI = p1a;    else pI = p1b;    // pick side of hiearchy
        if(randomDouble()<pI)  {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else {
        if(randomDouble()<p0)  {edges[iS].add(jS); edges[jS].add(iS); nEdges++;} 
      } // end else
    } // end for j
  } // end for i
  
  //cout << redText << "Here f " << endl;
  copyEdges(edges);
  //cout << redText << "Here g " << endl;
  initRandomWeights(bottomWeight,topWeight);

  //msg("done\n");
  return;
}; // end initSHierarchy


void TCMSparseW::initHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
     double p0, double p1, double p2, double p3, bool bRandomNodes,
     TCMData  bottomWeight, TCMData topWeight) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initHierarchy()");
  #endif
  //  N - is the total size of the system
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  int iS, jS;

  msg("Assigning random hierarchical communities... ");
  //cout << red << "Here 1 " << endl;
  TVectorInt s(N);
  //cout << red << "Here 2 " << endl;
  s.initStep(0,1,bRandomNodes);
  //cout << red << "Here 3 " << endl;
    
  delete[] kij;
  unstoredValue = -1;  // unweighted missing edge
 //cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[N];
  //cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(N);

  //cout << red << "Here e "  << endl;

  // loop over all pairs of nodes
  for(unsigned i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    for(unsigned j=i+1; j<N; j++) {
      jS = s[j];
      if((i/n3)==(j/n3)) {
        if(randomDouble()<p3) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n2)==(j/n2)) {
        if(randomDouble()<p2) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else if((i/n1)==(j/n1)) {
        if(randomDouble()<p1) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } else {
        if(randomDouble()<p0) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;} 
      } // end else
    } // end for j
  } // end for i
  
  //cout << red << "Here f "  << endl;

  //cout << redText << "Here f " << endl;
  copyEdges(edges);
  //cout << red << "Here g "  << endl;

  //cout << redText << "Here g " << endl;
  initRandomWeights(bottomWeight,topWeight);

  //cout << red << "Here h "  << endl;

  //msg("done\n");
  return;
}; // end initHierarchy


void TCMSparseW::initRandomRS(unsigned N, double randomDensity, 
           unsigned long long ISeed, TCMData  bottomWeight, TCMData topWeight) {
  //  N - is the total size of the system
  //cout << red << "here a in initRandom " << normText << flush;   // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initRandom()");
  #endif
  nNodes = N;
  //msg("Assigning random communities...");
  cout << "Assigning random communities... " << flush;
  cout << " with a seed of " << rngSeedString() << endl; 
  
  destroykij();  // just in case - still valid if NULL

  //cout << red << "here c in initRandom " << normText << flush;   // debugging
  // use O(L) version
  // tEdges is the number of expected edges
  int tEdges = (int)( randomDensity*(double)N*(double)(N-1)/2.0 + 0.5 );
  cout << red << "Number of edges expected is " << tEdges << endl;
  int iN, jN, tMod = tEdges/50;
  TVectorInt  degrees(N,0);  // store a count of the number of edges

  rngReseed(ISeed);  // start with initial seed here and just count edges
  //r.Reseed(ISeed);   // debugging
  //rngReseed(rngSeedString());  // start with initial seed here and just count edges
  //cout << red << "first 4 random numbers are "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " " << endl;
    
  for(unsigned i=0; i<tEdges; i++) {
    if((i%tMod)==0) cout << "." << flush;  // update user on progress
    iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
    // assumes symmetric
    while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
    // checking for duplicates at this point is expensive, but we should not
    // have many duplicates (for large systems), so increment nEdges and 
    // allow copyEdges to catch duplicates
    //edges[iN].add(jN);  edges[jN].add(iN);  nEdges++;
    degrees[iN] += 1;  degrees[jN] += 1;
  } // end for i

  cout << "Allocating degree data... " << flush;
  //cout << endl << degrees << endl << " with i = " << flush;
  // now create kij at the (almost, save for duplicates) correct sizes
  createkij(N);  // allocate the rows (columns)
  data = new TVectorInt[N];
  for(unsigned i=0; i<N; i++) {
    //cout << i << "," << flush;
    kij[i].resize((unsigned)degrees[i]);   // actual memory allocation
    //cout << " " << flush;
    kij[i].resize(0,0,0,1);    // fake resize to zero to prepare insertion sort
  } // end for i

  for(unsigned i=0; i<N; i++) {
    // I think this will benefit cache to allocate in a separate 'chunk'
    data[i].resize((unsigned)degrees[i]);  
  } // end for i

  // now regenerate the entire list of random numbers to fill the list
  rngReseed(ISeed);  // start with initial seed here and just count edges
  //r.Reseed(ISeed);   // debugging
  //rngReseed(rngSeedString());  // start with initial seed here and just count edges
  //cout << red << "first 4 random numbers are "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " "
  //     << randomInt(0,N-1) << " " << randomInt(0,N-1) << " " << endl;

  cout << endl << "Filling degree data... " << flush;
  bool bInserted;
  for(unsigned i=0; i<tEdges; i++) {
    if((i%tMod)==0) cout << "." << flush;  // update user on progress
    //cout << "a " << flush;  // debugging
    iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
    //cout << "b with iN = " << iN << " and jN = " << jN << flush;  // debugging
    //cout << kij[iN] << endl;  // debugging
    //cout << kij[jN] << endl;  // debugging
    // assumes symmetric
    while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
    //cout << "c " << flush;  // debugging
    // now try to add the node to an existing list using insertion sort logic
    // duplicates are checked in the insertion (sort) where the parameter
    // is set to ignore duplicates (duplicate memory is not returned since it 
    // is assumed to be small relative to the size of kij[][])
    //cout << "d " << flush;  // debugging
    bInserted = insertNode(iN,jN,1);
    //cout << "e " << flush;  // debugging
    if(bInserted) {
      insertNode(jN,iN,1);  // assumes symmetric edges
      nEdges++;
    } // end if bInserted
    //cout << "r " << flush;  // debugging
  } // end for i
  cout << "\n";  // finish user updates

  //nEdges /= 2;  // correct for double counting edges
  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1));  // calculate actual density

  //for(unsigned i=0; i<N; i++)  
  //  cout << blue << "kij[" << i << "] = " << green << kij[i] << endl;

  //initRandomWeights(bottomWeight,topWeight);
  for(unsigned i=0; i<N; i++) {
    data[i].init(1);
  } // end for i
  unstoredValue = -1;  // unweighted missing edge


  //msg("done.");
  return;
}; // end initRandomRS


//template <typename T> 
void TCMSparseW::initRandomLL(unsigned N, double randomDensity, 
                            TCMData  bottomWeight, TCMData topWeight, bool bN2){
  //  N - is the total size of the system
  //  p - is a short vector describing both the levels and probabilities of the 
  //      hierachies,  For example, [0.9, 0.5, 0.1] would define a 3 level 
  //      hiearchy with the most dense (smallest) structures at 0.9 probability
  //      for the moment it is crudely implemented as 
  //  s - is the corresponding vector to p that describes the sizes of the 
  //      communities and sub-communities.  Note that they must be evenly 
  //      divisible upwards.
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  //cout << red << "here a in initRandom " << normText << flush;      // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in initRandom()");
  #endif
  nNodes = N;
  msg("Assigning random communities...");

  //cout << red << "here b in initRandom " << normText << flush;      // debugging
  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file.  We then copy it
  // over to the more constant sparse matrix structure within this class.
  TListArray<int>  edges(N);
  
  //cout << red << "here c in initRandom " << normText << flush;   // debugging
  if(bN2) { // use N^2 version (for whatever reason)
    // loop over all pairs of nodes
    for(unsigned i=0; i<N; i++) 
      for(unsigned j=i+1; j<N; j++) 
        if(randomDouble()<randomDensity) {
          edges[i].add(j);  edges[j].add(i);  nEdges++; }
  } else { // use O(L) version
    // use a more efficient way O(L) rather than O(N^2)
    // number of expected edges
    int tEdges = (int)( randomDensity*(double)N*(double)(N-1)/2.0 + 0.5 );  
    cout << red << "Number of edges expected is " << tEdges << endl;
    int iN, jN, tMod = tEdges/50;
    for(unsigned i=0; i<tEdges; i++) {
      if((i%tMod)==0) cout << "." << flush;  // update user on progress
      iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
      // assumes symmetric
      while(iN==jN)  jN = randomInt(0,N-1);  // try again on duplicate
      // checking for duplicates at this point is expensive, but we should not
      // have many duplicates (for large systems), so increment nEdges and 
      // allow copyEdges to catch duplicates
      edges[iN].add(jN);  edges[jN].add(iN);  nEdges++;
    } // end for i
    //cout << "\n";  // end user update
  } // end else bN2

  //cout << edges << endl; // debugging

  if(bN2)  copyEdges(edges);    // no duplicates possible in N2 version
  else     copyEdges(edges,1);  // catch duplicates with random version
  initRandomWeights(bottomWeight,topWeight);

  //msg("done.");
  return;
}; // end initRandomLL


//template <typename T> 
int TCMSparseW::input(string fname, int sOffset, bool bDirected) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in input()");
  #endif

  int errorCode;
  // text matrix format
  //if(fname[fname.size()-4]=='.' && fname[fname.size()-3]=='t' && 
  //   fname[fname.size()-2]=='x' && fname[fname.size()-1]=='t')
  //  errorCode = Cp.inputTextMatrix(Cp,nodeOffset,bDirected);
  //else 
  // GML format
  if(fname[fname.size()-4]=='.' && fname[fname.size()-3]=='g' && 
     fname[fname.size()-2]=='m' && fname[fname.size()-1]=='l')
    errorCode = inputGML(fname,sOffset,bDirected);
  else if(fname[fname.size()-4]=='.' && fname[fname.size()-3]=='d' && 
          fname[fname.size()-2]=='a' && fname[fname.size()-1]=='t')
    // LFR data file format - pure edge list (with no header information)
    //errorCode = inputLFRDAT(fname,sOffset);  
    errorMsg("LFRDAT input currently disabled");
  else  
    errorMsg("TCMSparseW::input() did not recognize the data file type for \""
             +fname+"\".  Known file type is currently only gml.");

  if(errorCode<=0)  errorMsg("There was a problem reading the input file.");

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
}; // end input


//template <typename T> 
int TCMSparseW::inputGML(string fname, int soffset, 
                         bool bDirected, bool bSymmetrize) {
  // Input the CMatrix undirected graph in gml format
  // NOTE:  It is not thoroughly implemented!!!
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in inputGML()");
  #endif
  string sWord;  // temporary line variable
  char   cNext;  // temporary character variable - good for checking line type
  const int MaxLineLength = 1024;
  char   buffer[MaxLineLength];  // char* variable for getline() function
  int    p1, p2;              // particles 1 and 2 of the defined edge
  int    weight=1;            // weight for weighted graphs

  if(fname=="")  errorMsg("No filename specified in TCMSpare input()!");

  nNodes = 0;  nEdges = 0;

  cout << green << "Importing file: \"" << fname 
       << "\" assuming weighted edges " << normText << "... ";

  ifstream din(fname.c_str());
  if(din.bad()) {
    errorMsg("Input filename "+fname+" does not exist!");
    return -1;  // Error filename does not exist
  } // if din.bad
  cout << "  Reading GML data...\n";
  
  // Ascii input stuff here - skip comment lines that begin with a 'c' character
  bool bStart = 0;
  while(!din.eof() && !bStart) {
    din >> sWord;                    // read first character of this line
    cNext = (char)( sWord[0] );      // look at the first character a case test
    // need to properly read the initial info and skip it...
    // for the moment we ignore everything except node and edge designations
    if(cNext=='n')  bStart = 1;      // cut out to read data next
    else  din.getline(buffer,MaxLineLength); // otherwise skip rest of line
    //cout << redText << "Line read as: (" << sWord << ") " << buffer << endl; // debugging
  } // end while !eof
  din.getline(buffer,MaxLineLength); // get rest of first node line

  cout << magenta << "Reading nodes... "; // debugging
  // Ascii input stuff here
  //cout << redText << "In node read " << normText << endl; // debugging
  bool bStartEdges = 0;
  //din >> sWord;    // read first character of this line
  //din.getline(buffer,MaxLineLength); // get rest of line
  while(!din.eof() && !bStartEdges) {
      din >> sWord;  // read first word of this line
      cNext = (char)( sWord[0] );  // look at the first character a case test
      //cout << redText << sWord << " Line read as: " << buffer
      //     << " and cNext = " << cNext << " with nNodes = "
      //     << nNodes << normText << endl; // debugging
      switch(cNext) {

        case 'i': {  // skip whole comment line
          nNodes++;  // count the current edge
          //cout << nNodes << " ";  // debugging
          din.getline(buffer,MaxLineLength); // get rest of line after first char
          break;
        } // end case c

        case 'v': { // value line
          //cout << redText << (char)din.get() << (char)din.get() 
          //                  << (char)din.get() << " ";  // debugging
/*
          din.get();  din.get();
          //cout << "\'" << (char)din.peek() << "\'"; 
          tChar = (char)din.get();  
          //din >> tChar;
          cout << "\'" << tChar << "\'"; 
          if(tChar=='c')       a[0].add(nNodes-1);
          else if(tChar=='n')  a[1].add(nNodes-1);
          else if(tChar=='l')  a[2].add(nNodes-1);
*/
          din.getline(buffer,MaxLineLength); // get rest of line after chars
          break;
        } // end case c
 
        case 'e': { // begin edge readings
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          bStartEdges = 1;
          break;
        } // end case c

        case 'l':  case 'n':  case '[':  case ']': { // end of edge section
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          break;
        } // end case c

        default: {
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          // error - close file and exit
          din.close();
          return -2;  // Corrupted DIMACS data file
          break; // redundant break for default
        } // end default
      } // end switch for node input
  } // end while !eof

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file.  We then copy it
  // over to the more constant sparse matrix structure within this class.
  TListArray<TCMDataNode>  edges(nNodes);
  TCMDataNode d;  // temporary node while reading edges
  
  bool  bEdgeDone = 0;
  // Ascii input stuff here
  cout << "Reading edges... " << flush; // debugging
  //din >> sWord;    // read first character of this line
  //if(sWord[0]=='e') {    
  //din.getline(buffer,MaxLineLength); // get rest of line
  p1 = -1;  p2 = -1;  weight = 1;
   while(!din.eof()) {
      din >> sWord;  // read first word of this line
      cNext = (char)( sWord[0] );  // look at the first character a case test
      //cout << redText << sWord << " Edge read as: " << buffer 
      //     << " and cNext = " << cNext << " with nNodes = " 
      //     << nNodes << normText << endl; // debugging
      switch(cNext) {

        case 's': { // source line
          din >> p1;  // count the current edge
          din.getline(buffer,MaxLineLength); // get rest of line after first char
          break;
        } // end case c

        case 't': { // target line
          din >> p2;  // count the current edge
          din.getline(buffer,MaxLineLength); // get rest of line after first char
          break;
        } // end case c
 
        case 'v': { // value line
          //cout << "Weight for edge " << nEdges << " connecting nodes " 
          //     << p1-soffset << " to " <<p2-soffset; // debugging
          din >> weight;  // count the current edge
          //cout << " is " << weight << grey << " (assuming nodes defined first) " 
          //     << normText << "\n"; // debugging
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          break;
        } // end case c
 
        case 'e':  case '[': {   // end of edge section
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          break;
        } // end case c

        case ']': {   // end of edge section
          bEdgeDone = 1;
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          break;
        } // end case c

        default: {
          din.getline(buffer,MaxLineLength);// get rest of line after first char
          // error - close file and exit
          din.close();
          return -2;  // Corrupted DIMACS data file
          break; // redundant break for default
        } // end default
      } // end switch for node input

      // now assign the edge - hack fix
      // NOTE:  This assumes that weight has to be >=0 integer
      if(p1>=0 && p2>=0 && bEdgeDone) {
        nEdges++;
        d.index = p2-soffset;  d.value = weight;
        edges[p1-soffset].add(d);
        //if(!bDirected) {  // we are currently automatically symmetrizing
        // by crudely just summing (rather than averaging) the weights
          d.index = p1-soffset;
          edges[p2-soffset].add(d);  // add symmetric edge
        //}
        //cout << green << "Weight for edge " << nEdges << " connecting nodes " 
        //     << p1-soffset << " to " <<p2-soffset << " is " << weight << "\n"; // debugging
        //CM[p1-soffset][p2-soffset] = weight;
        // for the moment we are only using undirected
        //CM[p2-soffset][p1-soffset] = weight;  
        weight = 1;     // return weight to default value for next read
        bEdgeDone = 0;
        p1 = -1;  p2 = -1;
        //cout << nEdges << " ";  // debugging
      } // end if p1, p2, weight
    } // end while !eof
  //} // end if sWord

  
  int nDuplicates = 0;
  // third index is a dummy index, possible ambiguity class if TCMData is int
  // add duplicates - crude symmetrize without the /2 afterward - did not want
  // to lose data on integer division.  This should be changed after TCMData
  // can be a float or double in the main code
  if(bDirected)  nDuplicates = copyEdges(edges,2,0);  
  // copy edges but check for duplicates as a consistency check, should be zero
  else           nDuplicates = copyEdges(edges,0,0);
  cout << "Done.\n  The number of GML nodes was " << nNodes << " and edges was "
       << nEdges << ".  ";
  if(bDirected)  
    cout << "There were " << nDuplicates << " added edges (directed graph).";
  else  if(nDuplicates>0) 
          warningMsg("There were "+itos(nDuplicates)+" duplicate edges");
  cout << "done." << endl;
/*       
  // perform a final consistency check on the number of edges and what is
  // specified in the problem definition
  if(edges.getSize()!=nNodes) {
        cout << red << "  Warning:  Number of nodes from edge definitions "
             << "does not match the problem definition.\n  File specified " 
             << nNodes << " nodes, but the number of nodes were counted as " 
             << edges.getSize() << normText << endl;
        din.close();
        return 2;  // number of nodes does not match problem definition
  } // end edges check
*/
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif

  // successfully finished reading GML data file - close and exit
  din.close();
  return 1;
}; // end inputGML


//template <typename T> 
int TCMSparseW::inputLFRDAT(string fname, int N) {
  // Input the CMatrix undirected graph in Lancichinetti et al. network.dat 
  // format.  Note that they use a node offset of 1 which we automatically
  // correct internally to be 0 based.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparseW object is not valid in inputLFRDAT()");
  #endif
  string sWord;  // temporary line variable
  char   cNext;  // temporary character variable - good for checking line type
  const int MaxLineLength = 1024;
  char   buffer[MaxLineLength]; // char* variable for getline() function
  int    p1, p2;                // particles 1 and 2 of the defined edge
  int    weight = 1;            // LFR DAT is only defined for unweighted graphs
  int    sOffset = 1;   // Note that they use a node offset of 1 which we 
                        // automatically correct internally to be 0 based.
  bool   bDirected = 0; // assume symmetric graph
  if(fname=="")  errorMsg("No filename specified in TCMSparseW inputLFRDAT()!");

  nNodes = N;  nEdges = 0;
  int maxNode = -1;  // for consistency track max node number and check at end
  unstoredValue = -1;

  //cout << green << "Importing file: \"" << fname 
  //     << "\" (LFR DAT assumes unweighted edges) ... ";

  ifstream din(fname.c_str());
  if(din.bad()) {
    errorMsg("Input filename "+fname+" does not exist!");
    return -1;  // Error filename does not exist
  } // if din.bad
  cout << "  Reading LFR DAT file... ";
  
  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file.  We then copy it
  // over to the more constant sparse matrix structure within this class.
  TListArray<TCMDataNode>  edges(nNodes);
  TCMDataNode d;  // temporary node while reading edges
  
  //int eCount = 0;
  d.value = weight;  // unweighted graph, so initialize outside loop
  while(!din.eof()) {
    din >> p1 >> p2;                    // read first character of this line
    din.getline(buffer,MaxLineLength);  // skip rest of edge line

    if(p1>=sOffset && p2>=sOffset && p1<(N+sOffset) && p2<(N+sOffset)) {
      nEdges++;  //eCount++;
      d.index = p2-sOffset;  edges[p1-sOffset].add(d);
      // omit symmetric edge fo LFR data because they double declare edges
      //d.index = p1-sOffset;  edges[p2-sOffset].add(d);  // add symmetric edge
      //weight = 1;     // LFR DAT is always unweighted
      // now track the max node number for consistency check
      if(p1>maxNode)  maxNode = p1;
      if(p2>maxNode)  maxNode = p1;
      p1 = -1;  p2 = -1;  // reset to invalid values
    } // end if p1, p2, weight
    else if(!din.eof())
           errorMsg("Invalid edge read in inputLFRDAT? ("
                    +itos(p1)+" and "+itos(p2)+" on edge "+itos(nEdges)+")");
  } // end while !eof
  //cout << red << "L = " << nEdges << normText << endl; // debugging
  if((nEdges%2)==1)  errorMsg("Invalid edge count in LFR DAT file (L%2 = 1)?");
  else               nEdges /= 2;  // correct for double counted edges

  // now copy edges to actual internal data structure
  int nDuplicates = 0;
  // third index is a dummy index, possible ambiguity class if TCMData is int
  // add duplicates - crude symmetrize without the /2 afterward - did not want
  // to lose data on integer division.  This should be changed after TCMData
  // can be a float or double in the main code
  if(bDirected)  nDuplicates = copyEdges(edges,2,0);  
  // copy edges but check for duplicates as a consistency check, should be zero
  else           nDuplicates = copyEdges(edges,0,0);
  cout << "LFR nodes = " << nNodes << " and edges = "
       << nEdges << ".  ";

  if(bDirected)  
    cout << "There were " << nDuplicates << " added edges (directed graph).";
  else if(nDuplicates>0) 
         warningMsg("There were "+itos(nDuplicates)+" duplicate edges");
  //cout << "done.\n";
       
  // perform a final consistency check on the number of edges and what is
  // specified in the problem definition
  if(maxNode!=nNodes) {
        cout << red << "  Warning:  Number of nodes from edge definitions "
             << "does not match the problem definition.\n  User specified " 
             << nNodes  << " nodes, but the number of nodes wa identified as " 
             << maxNode << normText << endl;
        din.close();
        return 2;  // number of nodes does not match problem definition
  } // end edges check

  #ifdef DEBUG_MODE
  for(int i=0; i<N; i++)  for(int j=0; j<ki(i); j++) 
    if( (i==j && (*this)(i,j)!=0) || (i!=j && abs((*this)(i,j))!=1) )
      errorMsg("Node weight is invalid in LFR import.  Cp("
               +itos(i)+","+itos(j)+") = "+itos((*this)(i,j))+"?");
  ValidArray = ValidCMatrixFlag;
  #endif

  cout << grey << "done reading file.\n" << normText;
  // successfully finished reading LFR data file - close and exit
  din.close();
  return 1;
}; // end inputLFRDAT


int TCMSparseW::exportGML(string fname, bool bDirected) {
  int N = nNodes;
  //ofstream  fout(fname.c_str(),ios::out | ios::app);
  ofstream  fout(fname.c_str(),ios::out);
  fout << "Created by Peter Ronhovde\n" << description;  // output CMatrix notes
  fout << "graph\n[\n";  // begin graph section

  // output node identities
  for(int i=0; i<N; i++)  
    fout << "  node\n  [\n    id " << i //<< "\n    degree " << kij[i].getSize()
         << "\n  ]\n";
  
  // output edges identities
  int jNode, j;
  for(int i=0; i<N; i++) {
    j = 0;
    jNode = kij[i][0];  // jNode starts at first neighbor node
    while(j<kij[i].getSize() && jNode<i) {
      jNode = kij[i][j];  // get next jNode
      if(jNode<i) 
        fout << "  edge\n  [\n    source " << i << "\n    target " << jNode
             << "\n    value " << data[i][j] << "\n  ]\n";
      j++;
    } // end while j
  } // end for i

  fout << "]\n";  // end graph section
  fout.close();
  return 1;
}; // end exportGML



// other non-member functions....?????????????????????
//template <typename T> 
TCMSparseW& TCMSparseW::operator=(const TListArray<int>& b) { 
  // assignment for a TListArray to an unweighted sparse matrix
  #ifdef DEBUG_MODE
  if(!isValid(1))  
    errorMsg("TCMSparseW object is invalid in operator=(TListArray&)");
  else if(!b.isValid())  
    errorMsg("TListArray object is invalid in TCMSparseW op=(TListArray&)");
  #endif
  Cols = b.getSize();  nNodes = b.getSize();   // assumes a square matrix
  if(Cols>0) {
    if(kij!=NULL)  delete[] kij;
    kij = new Array1D<int>[Cols];              // Set up kij columns
    int maxRows = -1;
    nEdges = 0;
    delete[] data;
    for(unsigned i=0; i<Cols; i++) {
      nEdges += b[i].getSize();         // sum number of edges
      kij[i]  = b[i];      // use overloaded Array1D<T> operator=(TList&)
      data[i].resize(b[i].getSize());
      for(unsigned j=0; j<Cols; j++)  data[i][j] = 1; // use a default value of 1
      if(b[i].getSize()>maxRows)  maxRows = b[i].getSize();
    } // end for i
    nEdges /= 2;
    p = (double)nEdges/(double)( nNodes*(nNodes-1) );
    if(maxRows>Cols)                        // assumes a square matrix
      errorMsg("TCMSparseW TList '=' row is too large. Exceeds max size.");
    // no data values can be specified in the passed parameter, so we set the 
    // default stored and unstored data values
    unstoredValue = 0;
    description = "";
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end if Cols>0
  else {
    kij = NULL;  data = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return *this;
}; // end =


int TCMSparseW::output_dimacs(string fname, bool useBinary) const {
  // Output the CM undirected graph in binary format.
  // The routine actually exports in ascii mode not binary mode since the edge
  // connection definitions are easily calculated as 'char's rather than
  // converting to some binary format.
  
  int i, j;
  int preLen;
  int NPar = Cols;  // legacy definition
  // use string variables to stored preamble so that I can readily calculate
  // the preamble length required for the binary format
  string sLine1("c Graph definition in DIMACS format\n");
  string sLine2("c Output from dclique.cpp code by Brian Albright and PR\n");
  //string sLine3("c RandomLib seed is " + ISeed + "\n");

  if(fname=="") {
    if(useBinary)  fname = "cm_dimacs.b";
    else           fname = "cm_dimacs.dat";
  } // end if fname
  cout << greenText << "Creating DIMACS data file \"" << fname << "\"" << endl;
  ofstream dout(fname.c_str(), ios::out|ios::trunc);
  // need to count the edges for the dimacs problem format
  int eCount = nEdges; // legacy definition
  //for(i=0; i<NPar; ++i) { for(j=i+1; j<NPar; ++j) {
  //  if((*this)(i,j)==1)  eCount++; }} // end for ij

  // now calculate the preamble length for dimacs binary format requirement
  preLen = sLine1.size() + sLine2.size() +
           (int)log10(NPar)+1 + (int)log10(eCount)+1
         //+ "p edge " + # + " " # + \n // now add problem definition space
           + 7             +  1    + 1;

  if(useBinary)  dout << preLen << "\n"; // preamble length required in binary
  // write comment lines for both formats
  dout << sLine1 << sLine2;// << sLine3;
  // write "problem" definition ("edge" is just part of the definition)
  dout << "p edge " << NPar << " " << eCount << "\n";

  if(useBinary) {
    int bitCounter = 0, currentBitCount = 0;
    unsigned char bitChar = 0;  // need an 8 bit data type for logic below
    // now write edge descriptors for CM in DIMACS binary format
    // the output converts the lower triangular section of CM to bits
    // in the order encountered by shifting bits in an unsigned char variable.
    // Note that we correct for a number of bits that is not divisible by 8
    // during the bitChar output section.
    for(i=0; i<NPar; i++) {  
      bitCounter = 0;
      for(j=0; j<=i; j++) { 
        // The order that I am scanning the array is logic oriented by the 
        // binary format definition rather than the most efficient way.  Since
        // the format here deals with undirected graphs, I could invert i and j
        //    for(j=0; j<NPar; j++) {  
        //      for(i=0; i<=j; i++) { 
        //        ... 
        //        if((bitCounter%8==0) || (i==NPar-1)) {  
        //          ...
        // loops to be much more efficient with properly optimized code.
        // The order the bits are encountered is important.
        bitCounter++;  currentBitCount = bitCounter%8;
      
        bitChar = bitChar << 1;  // shift bits to left to clear least significant bit
        // now, set left-most (most significant) bit to current CM value
        //bitChar = (unsigned char)( (int)bitChar + CM[j][i] );
        // Bitwise or effectively adds CM[j][i] to bitChar here since
        // CM is 0 or 1 and we have already cleared the 1 bit above.
        bitChar = bitChar | (unsigned char)(*this)(j,i);  
        // output char to buffer if we have read 8 bits since last write
        // The j==i check accounts for reading the last bits of the row where
        // the set of bits may not be a full char's worth of data.  This is
        // required due to what seems a clunky binary format specification.
        if((currentBitCount==0) || (j==i)) {  
          // if last bit, correct for the fact that we add bits on the right
          if((j==i) || (currentBitCount!=0)) 
             bitChar = bitChar << (8-currentBitCount);
          
          dout << bitChar;  
          bitChar = 0;  // reset bit counter for next 8 bits
        } // end if bitCounter
      } // end for j  
    } // end for i
  } // end if useBinary edge output
  else { 
    // ok, output ascii format since user specified not to use binary
    // write edge descriptors for CM in DIMACS ascii format
    // Note that the DIMACS formate evidently specifies numbering from 1 not 0
    for(i=0; i<NPar; i++)  for(j=0; j<=i; j++)
      if((*this)(j,i)==1) dout << "e " << (i+1) << " " << (j+1) << "\n";
  } // end else useBinary
  

  // finished - write file and exit
  dout.close();
  return 1;
} // end outputCMatrix


//template <typename T> 
void TCMSparseW::display(int bOutputLarge, bool bShowZeroes) const {
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
      maxInt = 1;
      //warningMsg("Matrix is a zero matrix!  Exiting matrix output.");  
      //return; 
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
  cout << "done with display." << endl;
  return;
}; // end display


//template <typename T> 
void TCMSparseW::disp(int bOutputLarge, bool bShowZeroes, 
                         bool bForceOutput) const {
  int i, j, cij, N = nNodes;

  //cout << blueText << " CM " << CMatrix[1][0] << normText;
  //if(N>96 && !bForceOutput) || !bShowZeroes)  return;  // debugging
  //if(N>64 && !(bool)bOutputLarge)  return;  // debugging
  //if(bOutputLarge < 2) return;  // debugging

  if(bOutputLarge>0) { // scan for largest integer
    int maxInt = 0;
    for(i=0; i<N; i++)  for(j=0; j<N; j++) 
      if(abs((*this)(i,j))>maxInt)  maxInt = abs((*this)(i,j));
    if(maxInt==0) { 
      maxInt = 1;
      //warningMsg("Matrix is a zero matrix!  Exiting matrix output.");  
      //return; 
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
  cout << "done with display." << endl;
  return;
}; // end display

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Overloaded iostream operators
//template <typename T>
ostream& operator<<(ostream& fout, TCMSparseW& a) {
  a.display(1,1);
  return fout;
};

} // end namespace MSL
//#endif

