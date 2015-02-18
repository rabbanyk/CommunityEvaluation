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

//#ifndef MSL_TCMSPARSE_UNWEIGHTED_H

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
#ifndef MSL_TCMSPARSE_UNWEIGHTED_H
#include "MSL_CMatrix_Sparse_Unweighted.h"
#endif
#ifndef MSL_TCLUSTERLIST_H
#include "../clusterclasses.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "./msl/MSL_Array1D_Template.h"
#endif

namespace MSL {

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
//------------------ TCMSparse Member Function Declarations -----------------

//Assign default static parameters
//TCMSparse::StrictErrorChecking = 1;

//template <typename T> 
TCMSparse::TCMSparse(unsigned Cs, TCMData storedWeight, TCMData unstoredWeight,
                     string d) : TCMatrix(Cs,storedWeight,unstoredWeight,d) {
  // Raw constructor that skips the initialization
  // As the default constructor, it generates a TCMSparse object that is a 
  // zero length real TCMSparse object but officially an invalid one since it
  // holds no data.
  #ifdef DEBUG_MODE
  debugMsg("Constructing TCMSparse object... ");
  #endif
  storedValue = storedWeight;  unstoredValue = unstoredWeight;
  #ifdef DEBUG_MODE
  debugMsg("done.\n");
  #endif
}; // end TCMSparse constructor


//template <typename T> 
TCMSparse::TCMSparse(TCMSparse &b) : TCMatrix(b){ // copy constructor
  #ifdef DEBUG_MODE
  debugMsg("Copying TCMSparse object\n");
  if(!isValid(1)) errorMsg("TCMSparse object is not valid in copy constructor");
  cout << "Entering TCMSparse destructor... " << endl;
  #endif
  unstoredValue = b.unstoredValue;  storedValue = b.storedValue;
}; // end TCMSparse copy constructor


//template <typename T> 
TCMSparse::~TCMSparse() {
  #ifdef DEBUG_MODE
  //debugMsg("Entering TCMSparse destructor... ",brown);
  if(!isValid(1))  errorMsg("TCMSparse object is not valid in operator=()");
  #endif
  debugMsg("Entering TCMSparse destructor... ",brown);
  //Rows = 0;  Cols = 0;  p = 0.0; // i.e. do nothing here
  //#ifdef DEBUG_MODE
  //ValidArray = 0; let TCMatrix destructor set this due to reversed call order
  debugMsg("exiting TCMSparse destructor.\n",brown);
  //#endif
}; // end TCMSparse destructor


//template <typename T> 
TCMSparse& TCMSparse::operator=(const TCMSparse& b) {
  // assignment for 2 arrays
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("TCMSparse object is not valid in operator=()");
  #endif
  if(b.getSize()>Cols)                      // assumes a square matrix here
    errorMsg("TCMSparse row too large in '=' exceeds specified Rows size.");
  delete[] kij;
  Cols          = b.Cols;
  unstoredValue = b.unstoredValue;  storedValue = b.storedValue;
  nNodes        = b.nNodes;         nEdges      = b.nEdges;
  p             = b.p;
  gammaa        = b.gammaa;         gammab      = b.gammab;
  description   = b.description;
  if(Cols>0) {
    kij = new Array1D<int>[Cols];                // Set up kij columns
    for(int i=0; i<Cols; i++)  kij[i] = b.kij[i]; // use Array1D<T> optimized '=' 
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end f
  else {
    kij = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return *this;
}; // end =


//template <typename T> 
TCMData TCMSparse::operator()(unsigned i, unsigned j) const {
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
TCMData  TCMSparse::operator()(unsigned i, unsigned j) {
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


//template <typename T> 
TCMData TCMSparse::getMax(bool bAbs) const {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in getMax()");
  #endif
  if(bAbs) { if(abs(storedValue)>abs(unstoredValue))  
             return abs(storedValue);  else return abs(unstoredValue); }
  else     { if(storedValue>unstoredValue)  return storedValue;  
             else return unstoredValue; }
}; // end getMax


//template <typename T> 
TCMData TCMSparse::getMin(bool bAbs) const {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in getMin()");
  #endif
  if(bAbs) { if(abs(storedValue)<abs(unstoredValue))  
             return abs(storedValue);  else return abs(unstoredValue); }
  else     { if(storedValue<unstoredValue)  return storedValue;  
             else return unstoredValue; }
}; // end getMin

/*
void TCMSparse::initEven(unsigned N, int q, bool bRandomNodes) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in initEven()");
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
  for(int i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    for(int j=i+1; j<N; j++) {
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

void TCMSparse::copyEdges(TListArray<int> &edges, bool bCheckDuplicates) {
  // copy a TListArray to the kij edge list and update class variables
  int N = edges.getSize();
  delete[] kij;
  kij = new Array1D<int>[N];

  msg("copying to sparse matrix data structure... ");
  nEdges = 0;
  // we copy in reverse order so that we can return memory as data is copied
  // this prevents us from having to allocate 2x the memory for init stage
  // (the erase routine swaps pointers for optimized speed, but even if constant
  // order is enforced, it needlessly adds O(q^2) pointer copies which competes
  // for the dominant order of the hierarchical algorithm)
  //cout << endl;  // debugging
  for(int k=N-1; k>=0; k--) {
    //cout << brown << "Nodes " << k << ":  " << green << edges[k] << endl; // debugging
    //cout << "\nSorting list k = " << k << endl << edges[k] << flush; // debugging
    edges[k].sort();   // sort each list - required for binary sort in kij
                       // sort is O(N) if already in order or is reverse ordered
    //cout << " done with sort..." << flush << edges[k] << endl; // debugging

    if(bCheckDuplicates) { // check for duplicate entries end delete them
      // this uses the 'current' pointer logic for the TList class
      TList<int> *pk = &(edges[k]);
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

    //cout << edges[k];  // debugging
    kij[k] = edges[k]; // use overloaded equals operator of Array1D<T>
    nEdges += kij[k].getSize();
    edges.erase(k);    // now return the memory of this TList
    //cout << brown << " exiting " << endl;
  } // end for k
  
  nNodes = N;  Cols = N;
  p = (double)nEdges/(double)( N*(N-1) );  // calculate actual density
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
}; // end copyEdges


//template <typename T> 
void TCMSparse::initByClusters(TClusterList &c, double pin, double pout) {
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in initByClusters()");
  #endif
  int iS, jS;

  msg("Assigning communities by a specified set of clusters... ");
    
  cout << red  << "Here c with kij pointer value as " << (int)kij 
       << "\n" << endl; // debugging
  delete[] kij;
  cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[c.getNNodes()];
  cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(c.getNNodes());

  cout << red << "Here e "  << endl;

  for(int k=0; k<c.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iS = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jS = c[k].getNode(j);
        if(randomDouble()<pin) { 
          edges[iS].add(jS); edges[jS].add(iS); nEdges++; }
      } // end for j

      // loop over all other clusters m (past k) for outside connections
      for(int m=k+1; m<c.getSize(); m++) {
        // loop over nodes j in cluster m
        for(int j=0; j<c[m].getSize(); j++) {
          jS = c[m].getNode(j);
          if(randomDouble()<pout) { 
            edges[iS].add(jS); edges[jS].add(iS); nEdges++; }
        } // end for j
      } // end for m
  
    } // end for i
  } // end for k
  
  cout << redText << "Here f " << endl;
  copyEdges(edges);
  cout << redText << "Here g " << endl;

  msg("done\n");
  return;
}; // end initByClusters


void TCMSparse::initHeterogenous(TClusterList &c, Array1D<int> &v,
                                 double p2, double p1, double p0) {
  // initialize a heterogeneous hierarchy but with equal density propabilities
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  // v - specifies which clusters to merge.  i.e. if v = {1,3,2}, then 
  //     leave cluster 1 alone (the 1); merge clusters 2, 3, and 4 (the 3);
  //     and merge clusters 5 and 6 (the 2)
  // bRandomNodes - specifies whether the nodes will the randomized or 
  //                sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is invalid in initHeterogenous()");
  #endif
  if(v.getSum()!=c.getSize())  
    errorMsg("Merge vector sum != c.nClusters in initHeterogenous()");

  int iS, jS;

  msg("Assigning heterogeneous hierachical communities... ");
    
  //cout << red  << "Here c with kij pointer value as " << (int)kij 
  //     << "\n" << endl; // debugging
  delete[] kij;
  //cout << red << "Here c2 " << endl;
  kij = new Array1D<int>[c.getNNodes()];
  //cout << red << "Here d "  << endl;

  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file, and additionally
  // eliminates concerns over allocating enough memory for the edges.  We then 
  // copy it over to the more constant sparse matrix structure within this 
  // class.  The TList sort is optimized to handle ordered or reverse-ordered
  // elements in O(N).
  TListArray<int>  edges(c.getNNodes());

  //cout << red << "Here e "  << endl;

  // work with super-clusters first for level 1 communities
  TClusterList merged;
  merged = c;
  // merge original clusters by cluster merge list v.  We merge backwards in
  // order to not change the individual cluster order during the move(I,To_J) 
  // function calls although the ending merged order is actually mixed due to 
  // the move call speed optimizations.  Since we use merged internally here, 
  // this should not be a problem.
  for(int j=0; j<v.getSize(); j++)  // loop over all merges
    for(int k=0; k<v[v.getSize()-1-j]-1; k++)     // only merge if v[j]>1
      merged.move(merged.getSize()-j-1,merged.getSize()-j-2);

  // fill matrix
  for(int k=0; k<merged.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<merged[k].getSize(); i++) {
      iS = merged[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<merged[k].getSize(); j++) {
        jS = merged[k].getNode(j);
        if(randomDouble()<p1) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } // end for j

    } // end for i
  } // end for k
  
  // now work with original clusters for level 2 most dense communities
  for(int k=0; k<c.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<c[k].getSize(); i++) {
      iS = c[k].getNode(i);

      // loop over all pairs of nodes in cluster k (past i) for inside connect.
      for(int j=i+1; j<c[k].getSize(); j++) {
        jS = c[k].getNode(j);
        // this could add some duplicate links which are removed during the copy
        // stage below.  Otherwise, we have to do an O(k) search for each link
        // added to check whether it is a duplicate.
        if(randomDouble()<p2) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
      } // end for j

    } // end for i
  } // end for k
  
  // now fill the rest of the system super-cluster to super-cluster at p0
  for(int k=0; k<merged.getSize(); k++) {
    // loop over all nodes in cluster k
    for(int i=0; i<merged[k].getSize(); i++) {
      iS = merged[k].getNode(i);

      // loop over all other clusters m (past k) for outside connections
      for(int m=k+1; m<merged.getSize(); m++) {
        // loop over nodes j in cluster m
        for(int j=0; j<merged[m].getSize(); j++) {
          jS = merged[m].getNode(j);
          if(randomDouble()<p0) {edges[iS].add(jS); edges[jS].add(iS); nEdges++;}
        } // end for j
      } // end for m
  
    } // end for i
  } // end for k
  
  //cout << red << "Here f " << endl;
  for(int i=0; i<edges.getSize(); i++) {
    edges[i].sort();
    //cout << edges[i] << endl; // debugging
  } // end for i
  copyEdges(edges,1);  // copy and erase duplicate entries from step p2 above
  //cout << red << "Here g " << endl;

  msg("done\n");
  return;
}; // end initHeterogenous


//template <typename T> 
void TCMSparse::initSHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
     double p0, double p1a, double p2a, double p3a, 
                double p1b, double p2b, double p3b, bool bRandomNodes) {
  // initialize a 'staggered' hierarchy'.  That is, half of the hiearchy is 
  // initialized by p1a, p2a, and p3a for its levels and the other half is 
  // initialized by p1b, p2b, and p3b for its levels.
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.

  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in initHierarchy()");
  #endif
  //  N - is the total size of the system
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  int iS, jS;

  msg("Assigning random hierarchical communities... ");
  TVectorInt s(N);
  s.initStep(0,1,bRandomNodes);
    
  //cout << red  << "Here c with kij pointer value as " << (int)kij 
  //     << "\n" << endl; // debugging
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

  //cout << red << "Here e "  << endl;

  // loop over all pairs of nodes
  for(int i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    double pI, pII, pIII;
    for(int j=i+1; j<N; j++) {
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

  msg("done\n");
  return;
}; // end initSHierarchy


void TCMSparse::initHierarchy(unsigned N, unsigned n1, unsigned n2, unsigned n3,
     double p0, double p1, double p2, double p3, bool bRandomNodes) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in initHierarchy()");
  #endif
  //  N - is the total size of the system
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  int iS, jS;

  msg("Assigning random hierarchical communities... ");
  TVectorInt s(N);
  s.initStep(0,1,bRandomNodes);
    
  //cout << red  << "Here c with kij pointer value as " << (int)kij << "\n" 
  //     << endl; // debugging
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

  //cout << red << "Here e "  << endl;

  // loop over all pairs of nodes
  for(int i=0; i<N; i++) {
    if((i%(N/10))==0) cout << "." << flush;
    iS = s[i];
    for(int j=i+1; j<N; j++) {
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
  
  //cout << redText << "Here f " << endl;
  copyEdges(edges);
  //cout << redText << "Here g " << endl;

  msg("done\n");
  return;
}; // end initHierarchy


//template <typename T> 
void TCMSparse::initRandom(unsigned N, double randomDensity) {
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
  cout << red << "here a in initRandom " << normText << flush;      // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in initRandom()");
  #endif
  nNodes = N;
  msg("Assigning random communities...");

  cout << red << "here b in initRandom " << normText << flush;      // debugging
  // use a temporary holding structure of a TListArray.  This eliminates any
  // concerns over out of order edge listings in the gml file.  We then copy it
  // over to the more constant sparse matrix structure within this class.
  TListArray<int>  edges(N);
  
  //cout << red << "here c in initRandom " << normText << flush;      // debugging
  // loop over all pairs of nodes
  for(int i=0; i<N; i++) 
    for(int j=i+1; j<N; j++) 
      if(randomDouble()<randomDensity) {
        edges[i].add(j);  edges[j].add(i);  nEdges++; }
  // use a more efficient way O(L) rather than O(N^2)
/* // errror - causes duplicate connections
  int tEdges = (int)(p*(double)( N*(N-1) ) + 0.5 );  // number of expected edges
  int iN, jN;
  for(int i=0; i<tEdges; i++) {
    if((i%(tEdges/100))==0) cout << "." << flush;
    iN = randomInt(0,N-1);
    jN = randomInt(0,N-1);
    if(iN==jN)  jN = randomInt(0,N-1);  // try again on a match - crude
    // it is expensive to check for duplicates, so we ignore them
    if(iN!=jN) { edges[iN].add(jN);  edges[jN].add(iN); }
  } // end for i
*/
  copyEdges(edges);
  msg("done.");
  return;
}; // end initRandom

/*
//template <typename T> 
bool  TCMSparse::set(unsigned i, unsigned j) {
  // uses a binary search on Array1D data structure to 
  // use operator() if data is not being written within the array
  // returns a boolean value of whether the element already existed
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in set()");
  #endif
  Array1D<int>  *pii = &(kij[i]);
  int  high = pii->getSize()-1, low = 0, guess, jkij, lastStep;
  while(low<high) {
    guess = ((high+low) >> 1);
    jkij = kij[guess];
    if(jkij==j)       return 1;  // element already exists
    else if(jkij<j) { low  = guess + 1;  lastStep =  1; }
    else              { high = guess - 1;  lastStep = -1; }
  } // end while
  // if we make it here the element does not exists, so add it in sorted order
  // low==high with lastStep gives the location that the new element is to be 
  // added at that kij location
  Array1D<int>  *pb;  pb = new Array1D<int>; // create default size zero 
  pb->resize(pii->getSize()+1);              // enlarge to new size
  jkij = low+(lastStep+1)/2;  // where are we inserting the new value?
  for(int k=0; k<jkij; k++)              (*pb)[k] = (*pii)[k];
  (*pb)[jkij] = j; 
  for(int k=jkij+1; k<a.getSize(); k++)  (*pb)[k] = (*pii)[k-1];
  // now finally copy data with new member back to class data members
  //delete kij;
  kij[i] = (*pb);  // an extra O(Z) operation here due to encapsulation
                   // moving on for time constraints
  return 0;        // element had to be added
}; // end set


//template <typename T> 
bool  TCMSparse::unset(unsigned i, unsigned j) {
  // uses a binary search on Array1D data structure to find the element by kij
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in unset()");
  #endif
  bool bFound = 0;
  int  high = pii->getSize()-1, low = 0, guess, jkij;
  while(low<high) {
    guess = ((high+low) >> 1);
    jkij = kij[guess];
    if(jkij==j)      bFound = 1;
    else if(jkij<j)  low    = guess + 1;
    else             high   = guess - 1;
  } // end while
  // delete the stored element by copy and replace
  if(bFound) {
    Array1D<int>  *pb;  pb = new Array1D<int>; // create default size zero 
    pb->resize(pii->getSize()-1);              // enlarge to new size
    Array1D<int>  *pii = &(kij[i]);
    for(int k=0; k<jkij; k++)                 (*pb)[k] = (*pii)[k];
    for(int k=jkij; k<pii->getSize()-1; k++)  (*pb)[k] = (*pii)[k+1];
    // now finally copy data with new member back to class data members
    //delete kij;
    kij[i] = (*pb);  // an extra O(Z) operation here due to encapsulation
                       // moving on for time constraints
    return 1;          // element was deleted
  } // end if
  else return 0;       // element did not exist, so no deletion was made
}; // end set


//template <typename T> 
inline void TCMSparse::set(unsigned i, Array1D<int> &a) {
  // copies the passed kij array at kij location i
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in Array1D set()");
  #endif
  if(i<Cols)  kij[i] = a; // use the Array1D overloaded = operator for a TList
  else  errorMsg("Attempting to insert an kij Array1D that is too large!");
  return;  
}; // end set


//template <typename T> 
inline void TCMSparse::set(unsigned i, TList<int> &a) {
  // takes a TList object and copies the values to the kij array at kij i
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in TList set()");
  #endif
  if(i<Cols)  kij[i] = a;  // use defined Array1D<T> operator
  else        errorMsg("Attempting to insert an kij TList that is too large");
  return;  
}; // end set
*/
/*
//template <typename T> 
void initFixedValency(TClusterList &a, double zin, double zout) {
  // setup system if a fixed outside valency is specified by user
  // At the moment we assume a total number of connect edges per node of 16 
  // since that is the standard used in the comparison graphs in the literature
  double pin, pout;
  if(useFixedValency) {
    warningMsg("  User specified to use a fixed valency.");
    warningMsg("    -pout and -pin flags are ignored for fixed valency.",greyText);
    double n = (double)N/(double)QMax;
    pin  = (16.0 - zout)/(n - 1.0);
    pout = (16.0 - (n-1.0)*pin)/((double)N-n);
    cout << normText << "    Fixed valency values pin = " << pin 
                     << " and pout = " << pout << "\n";
  } // end if useFixedValency

  return;
}; // end initRandom


//template <typename T> 
void initByClusters(TClusterList &a, double pin, double pout) {
  // setup system with a random density within and without the clusters
  double pin, pout;
  if(useFixedValency) {
    warningMsg("  User specified to use a fixed valency.");
    warningMsg("    -pout and -pin flags are ignored for fixed valency.",greyText);
    double n = (double)N/(double)QMax;
    pin  = (16.0 - zout)/(n - 1.0);
    pout = (16.0 - (n-1.0)*pin)/((double)N-n);
    cout << normText << "    Fixed valency values pin = " << pin 
                     << " and pout = " << pout << "\n";
  } // end if useFixedValency

  return;
}; // end initRandom
*/

//template <typename T> 
int TCMSparse::input(string fname, int soffset) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in input()");
  #endif
/*
  bool   isBinary;
  string binOutputFile;
  if(useInputFile && NRunsMax==1) {  // only import external file once
    if(inputFileName[inputFileName.size()-2]=='.'  && 
       inputFileName[inputFileName.size()-1]=='b'  )  isBinary = 1;
    else                                              isBinary = 0;
      errorCode = Cp.inputGML(Cp,nodeOffset);
    if(errorCode<=0)  errorMsg("There was a problem reading the input file.");
  }
  inputGML(fname,soffset);
*/
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
}; // end input

//template <typename T> 
int TCMSparse::inputGML(string fname, int soffset) {
  // Input the CMatrix undirected graph in gml format
  // NOTE:  It is not thoroughly implemented!!!
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMSparse object is not valid in inputGML()");
  #endif
  string sWord;  // temporary line variable
  char   cNext;  // temporary character variable - good for checking line type
  const int MaxLineLength = 1024;
  char   buffer[MaxLineLength];  // char* variable for getline() function
  int    p1, p2;              // particles 1 and 2 of the defined edge
  int    weight=1;            // weight for weighted graphs

  if(fname=="")  errorMsg("No filename specified in TCMSpare input()!");

  nNodes = 0;  nEdges = 0;

  cout << greenText << "Importing file: \"" << fname << "\" " 
       << normText << "\n";

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
    cout << redText << "Line read as: (" << sWord << ") " << buffer << endl; // debugging
  } // end while !eof
  din.getline(buffer,MaxLineLength); // get rest of first node line

  cout << magentaText << "Reading nodes... "; // debugging
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
  TListArray<int>  edges(nNodes);
  
  bool  bEdgeDone = 0;
  // Ascii input stuff here
  cout << "Reading edges " << normText << endl; // debugging
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
          cout << "Weight for edge " << nEdges; // debugging
          din >> weight;  // count the current edge
          cout << " is " << weight << endl;     // debugging
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
        if(weight!=1)  
          warningMsg("Ignoring weight for edge "+itos(p1-soffset)+" to "
                     +itos(p1-soffset));
        edges[p1-soffset].add(p2-soffset);
        edges[p2-soffset].add(p1-soffset);  // add symmetric edge
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

  
  cout << "  The number of GML nodes was " << nNodes << " and edges was "
       << nEdges << ".  ";
  copyEdges(edges);
  msg("done.");
       
  // perform a final consistency check on the number of edges and what is
  // specified in the problem definition
  if(edges.getSize()!=nNodes) {
        cout << redText << "  Warning:  Number of nodes from edge definitions "
             << "does not match the problem definition.\n  File specified " 
             << nNodes << " nodes, but the number of nodes were counted as " 
             << edges.getSize() << normText << endl;
        din.close();
        return 2;  // number of nodes does not match problem definition
  } // end edges check
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif

  // successfully finished reading GML data file - close and exit
  din.close();
  return 1;
}; // end inputGML


//template <typename T> 
void TCMSparse::display(int bOutputLarge, bool bShowZeroes) const {
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



// other non-member functions....?????????????????????
//template <typename T> 
TCMSparse& TCMSparse::operator=(const TListArray<int>& b) { 
  // assignment for a TListArray to an unweighted sparse matrix
  #ifdef DEBUG_MODE
  if(!isValid(1))  
    errorMsg("TCMSparse object is invalid in operator=(TListArray&)");
  else if(!b.isValid())  
    errorMsg("TListArray object is invalid in TCMSparse op=(TListArray&)");
  #endif
  Cols = b.getSize();  nNodes = b.getSize();   // assumes a square matrix
  if(Cols>0) {
    if(kij!=NULL)  delete[] kij;
    kij = new Array1D<int>[Cols];              // Set up kij columns
    int maxRows = -1;
    nEdges = 0;
    for(int i=0; i<Cols; i++) {
      nEdges += b[i].getSize();         // sum number of edges
      kij[i]  = b[i];      // use overloaded Array1D<T> operator=(TList&)
      if(b[i].getSize()>maxRows)  maxRows = b[i].getSize();
    } // end for i
    nEdges /= 2;
    p = (double)nEdges/(double)( nNodes*(nNodes-1) );
    if(maxRows>Cols)                        // assumes a square matrix
      errorMsg("TCMSparse TList '=' row is too large. Exceeds max size.");
    // no data values can be specified in the passed parameter, so we set the 
    // default stored and unstored data values
    storedValue = 1;  unstoredValue = 0;
    description = "";
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end if Cols>0
  else {
    kij = NULL;
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return *this;
}; // end =


int TCMSparse::output_dimacs(string fname, bool useBinary) const {
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





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Overloaded iostream operators
//template <typename T>
ostream& operator<<(ostream& fout, TCMSparse& a) {
  a.display(1,1);
  return fout;
};
/*
//template <typename T> 
istream& operator>>(istream& fin, TCMSparse& a) {
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
      cerr << "TCMSparse MatrixMarket header incorrect.\n"; return fin; }
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
    a.~TCMSparse();  // explicit destructor call since it is being redefined
    //cout << "After  destructor in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    a.Rows = N;  a.Cols = M;
    a.data = new T*[a.Cols];
    for (int i=0; i<a.Cols; i++) { a.data[i] = new T[a.Rows]; }
    a.ValidArray = ValidCMatrixFlag;
    for(int j=0; j<a.getCols(); j++) { for(int i=0; i<a.getRows(); i++) {
      fin >> a.data[j][i]; }} // end ij
    //cout << "After redefinition in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    if(a.VerboseErrorReports!=0) {
      cout << "Array Rows x Cols = " << N << " " << M << "  "
           << "Number of comment lines " << nComments << ".\n";
    } // end if VerboseErrorReports
  } // end if
  return fin;
};
*/


} // end namespace MSL
//#endif

