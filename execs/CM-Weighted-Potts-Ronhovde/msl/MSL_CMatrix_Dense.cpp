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

//#ifndef MSL_TCMDense_WEIGHTED_H

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
#ifndef MSL_TCMDENSE_H
#include "MSL_CMatrix_Dense.h"
#endif
#ifndef MSL_TCLUSTERLIST_H
#include "../clusterclasses.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "MSL_Array2D_Template.h"
#endif
#include "MSL_Stats1D.h"  // debugging

namespace MSL {

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
//typedef TCMDense<TFloat> TSparseFloat; // TFloat  array type definition
//typedef TCMDense<int>    TCMDenseInt;  // integer array type definition
#endif


//------------------------------------------------------------------------------
//------------------ TCMDense Member Function Declarations -----------------

//Assign default static parameters
//TCMDense::StrictErrorChecking = 1;

//template <typename T> 
TCMDense::TCMDense(unsigned Cs, string d) : TCMatrix(Cs,1,-1,d) {
  // Raw constructor that skips the initialization
  // As the default constructor, it generates a TCMDense object that is a 
  // zero length real TCMDense object but officially an invalid one since it
  // holds no data.
  #ifdef DEBUG_MODE
  debugMsg("Constructing TCMDense object... ");
  #endif
  // now declare the sparse data
  data.resize(Cs,Cs);
  #ifdef DEBUG_MODE
  debugMsg("done.\n");
  #endif
}; // end TCMDense constructor


TCMDense::TCMDense(Array2D<TCMData> &a, unsigned Cs, string d) : 
          TCMatrix(Cs,1,-1,d) {
  #ifdef DEBUG_MODE
  debugMsg("Constructing TCMDense object from Array2D<int> object... ");
  #endif
  data.resize(a.getCols(),a.getCols());
  data = a;  // declare the data
  //debugMsg("with array flag = "+itos(ValidArray)+" a... ");
  // the mergeMatrix is sometimes passed in upper-triangular form
  // (with zero fill values).  We need to catch this.
  bool bZeroij, bZeroji;
  for(int i=0; i<nNodes; i++) {
    for(int j=i+1; j<nNodes; j++) {
      data[i][j] *= -1;  data[j][i] *= -1;  // is opposite sign
      bZeroij = (data[i][j] == 0);  bZeroji = (data[j][i] == 0);
      // catch a single zero value
      if((bZeroij || bZeroji) && !(bZeroij && bZeroji)) { // exclusive or
        if(bZeroij)  data[i][j] = data[j][i];
        else         data[j][i] = data[i][j];
      } // end if single zero catch
    } // end for j
  } // end for i
  //debugMsg("with array flag = "+itos(ValidArray)+" b... ");
  
  // set the number of edges and the density
  nEdges = 0;
  int degreeCount;  
  TVectorInt nDegrees(Cols,0);
  TCMData *pj;
  //debugMsg("with array flag = "+itos(ValidArray)+" c... ");
  createkij(nNodes);
  //debugMsg("with array flag = "+itos(ValidArray)+" d... ");
  for(int i=0; i<nNodes; i++) {

    pj = &(data[i][0]);
    for(int j=0; j<nNodes; j++) { 
      if((*pj)>0) { nEdges++;  nDegrees[i] += 1; }
      pj++;
    } // end for j

    // now create the degree list
    kij[i].resize(nDegrees[i]);
    pj = &(data[i][0]);
    degreeCount = 0;
    for(int j=0; j<nNodes; j++) { 
      if((*pj)>0)  { kij[i][degreeCount] = j;  degreeCount++; }
      pj++;
    } // end for j
  } // end for i
  //debugMsg("with array flag = "+itos(ValidArray)+" e... ");
  // correct for double counted edges
  nEdges /= 2;

  p = nEdges*2.0/( (double)nNodes*(double)(nNodes-1) );

  // now fill the rest of the information about the CMatrix
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag; // have to redefine since we destroyed kij
  debugMsg("done with array flag = "+itos(ValidArray)+"\n");
  #endif
}; // end TCMDense constructor


//template <typename T> 
TCMDense::TCMDense(TCMDense &b) : TCMatrix(b){ // copy constructor
  #ifdef DEBUG_MODE
  debugMsg("Copying TCMDense object\n");
  if(!isValid(1)) 
    errorMsg("TCMDense object is not valid in copy constructor");
  cout << "Entering TCMDense destructor... " << endl;
  #endif
  data = b.data;  // use optimized Array2D '='
}; // end TCMDense copy constructor


//template <typename T> 
TCMDense::~TCMDense() {
  #ifdef DEBUG_MODE
  debugMsg("Entering TCMDense destructor... ",brown);
  if(!isValid(1))  errorMsg("TCMDense object is not valid in operator=()");
  #endif
  //Rows = -1;  Cols = -1;  p = -1.0; // i.e. do nothing here
  // data automatically destructs with TArray2D destructor
  #ifdef DEBUG_MODE
  ValidArray = 0;
  debugMsg("exiting TCMDense destructor.\n",brown);
  #endif
}; // end TCMDense destructor


//template <typename T> 
void TCMDense::scale(TCMData a, TCMData b, TCMData aOld, TCMData bOld) {
  // scale the stored and unstored elements differently
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMDense object is not valid in scale()");
  #endif
  // this current implementation could be inefficient; but we allow integer 
  // weights, and we would have problems if we use integer division
  // try to use use optimized Array1D routines
  TCMData *pj;
  for(int i=0; i<Cols; i++) {
    pj = &(data[i][0]);
    for(int i=0; i<Cols; i++) { 
      if((*pj)>0)  (*pj) = ( (*pj)/aOld )*a; 
      else         (*pj) = ( (*pj)/bOld )*b; 
      pj++;
    } // end for j
  } // end for i
  return;
}; // end scale


// Other useful array functions
//template <typename T> 
TCMDense& TCMDense::operator=(const TCMDense& b) {
  // assignment for 2 arrays
  #ifdef DEBUG_MODE
  if(!isValid(1))  errorMsg("TCMDense object is not valid in operator=()");
  #endif
  if(b.getSize()>Cols)                      // assumes a square matrix here
    errorMsg("TCMDense row too large in '=' exceeds specified Rows size.");
  delete[] kij;
  Cols          = b.Cols;
  nNodes        = b.nNodes;         nEdges      = b.nEdges;
  p             = b.p;
  description   = b.description;
  if(Cols>0) {
    kij = new Array1D<int>[Cols];                 // Set up kij columns
    for(int i=0; i<Cols; i++)  kij[i] = b.kij[i]; // use Array1D<T> optimized routines
    data = b.data;
    #ifdef DEBUG_MODE
    ValidArray = ValidCMatrixFlag;
    #endif
  } // end if
  else {
    kij = NULL;  data.resize(0,0);
    #ifdef DEBUG_MODE
    ValidArray = 0;
    #endif
  } // end else
  return *this;
}; // end =


void TCMDense::initHeterogeneous(TClusterList &c, TClusterList &merged,
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

  //cout << "Just before assign... " << flush;
  msg("Assigning heterogeneous hierachical communities... ");

  data.resize(nNodes,nNodes);
  data.init(-1);
  for(int i=0; i<N; i++)  data[i][i] = 0;  // clear diagonal
  TVectorInt degreeCount(N,0);
    
  //cout << red  << "Here c with kij pointer value as " << (int)kij 
  //     << "\n" << endl; // debugging
  //delete[] kij;
  //cout << red << "Here c2 " << endl;
  //kij = new Array1D<int>[N];
  //cout << red << "Here d "  << endl;
  //cout << "after delete with N = " << N << "... " << flush;  // debugging

  // work with super-clusters first for level 1 communities
  //TClusterList merged;
  merged = c;
  // merge original clusters by cluster merge list v.  We merge backwards in
  // order to not change the individual cluster order during the move(I,To_J) 
  // function calls although the ending merged order is actually mixed due to 
  // the move call speed optimizations.  Since we use merged internally here, 
  // this should not be a problem.
  //cout << "before merges with pre-merge as:\n"; // << merged << flush; // debugging
  for(int j=0; j<v.getSize(); j++) { // loop over all merges
    //cout << magenta << j << ":  " << normText << flush;            // debugging
    for(int k=0; k<v[v.getSize()-1-j]-1; k++) { // only merge if v[j]>1
      //cout << red << k << " " << normText << flush;                // debugging
      merged.move(merged.getSize()-j-1,merged.getSize()-j-2,1);
    } // end for k
    //cout << endl;                                                  // debugging
  } // end for j

  //merged.display(0,0,"merged inside ");                            // debugging

  //int nNodes0 = 0, nNodes1 = 0, nNodes2 = 0;  // count number of nodes at each level
  int nEdges0 = 0, nEdges1 = 0, nEdges2 = 0;  // count number of nodes at each level
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
              //edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdges1++; }
              data[iS][jS] = 1;  data[jS][iS] = 1;  nEdges++;  nEdges1++; 
              degreeCount[iS] += 1;  degreeCount[jS] += 1;
            } // end if
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
          //edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdges2++;}
          data[iS][jS] = 1;  data[jS][iS] = 1;  nEdges++;  nEdges1++; 
          degreeCount[iS] += 1;  degreeCount[jS] += 1;
        } // end if
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
            //edges[iS].add(jS);  edges[jS].add(iS);  nEdges++;  nEdges0++; }
            data[iS][jS] = 1;  data[jS][iS] = 1;  nEdges++;  nEdges1++; 
            degreeCount[iS] += 1;  degreeCount[jS] += 1;
          } // end if
        } // end for j
      } // end for i
    } // end for m
  } // end for k
  
  // now assign the edge connect/degree data - I could make this more efficient
  // by incorpating the read above a vector at a time to get the size
  destroykij();
  createkij(nNodes);
  int edgeCount = 0;
  for(int j=0; j<nNodes; j++) {
    kij[j].resize(degreeCount[j]);
    edgeCount = 0;
    for(int i=0; i<nNodes; i++) 
      if(data[i][j]>0) { kij[j][edgeCount] = i;  edgeCount++; }
    // consistency check
    if(edgeCount!=degreeCount[j]) 
      errorMsg("Degree assignments are inconsistent?");
  } // end for j

  p = nEdges*2.0/( (double)nNodes*(double)(nNodes-1) );

  //initRandomWeights(bottomWeight,topWeight);  // temporarily removed

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
  for(int i=0; i<v.getSize(); i++)  description += "," + itos(v[i]);
  description += "\n\"Merged cluster size list for level 2 (Total " 
               + itos(merged.getSize()) + "):  \"";
  for(int i=0; i<merged.getSize(); i++)  
    description += "," + itos(merged[i].getSize());
  description += "\n\"Original cluster size list for level 3 (Total " 
               + itos(c.getSize()) + "):\"";
  for(int i=0; i<c.getSize(); i++)  
    description += "," + itos(c[i].getSize());
  description += "\n";  // finish description information

  msg("done\n");
  return;
}; // end initHeterogeneous


//template <typename T> 
int TCMDense::input(string fname, int sOffset, bool bDirected) {
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMDense object is not valid in input()");
  #endif

  int errorCode;
  // text matrix format
  if(fname[fname.size()-4]=='.' && fname[fname.size()-3]=='t' && 
     fname[fname.size()-2]=='x' && fname[fname.size()-1]=='t')
    errorCode = inputTextMatrix(fname,sOffset,bDirected);
  else 
  // GML format
  if(fname[fname.size()-4]=='.' && fname[fname.size()-3]=='g' && 
     fname[fname.size()-2]=='m' && fname[fname.size()-1]=='l')
    errorCode = inputGML(fname,sOffset,bDirected);
  else  
    errorMsg("TCMDense::input() did not recognize the data file type for \""+fname
             +"\".  Known file types are currently gml and text matrix (txt).");

  if(errorCode<=0)  errorMsg("There was a problem reading the input file.");

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
}; // end input


//template <typename T> 
int TCMDense::inputTextMatrix(string fname, int soffset, 
                              //bool bDirected, bool bSymmetrize) {
                              bool bDirected) {
  // Input the CMatrix a basic ascii text matrix delimited by spaces
  // the top of the file can have comment lines starting by '#'
  // the data is not assumed to be symmetric and this is not checked
  // after the comment data, two lines defining the number of nodes and edges
  // is expected, such as:
  // N = 64
  // L = 1192
  // 0 1 5 0 0 ...
  // NOTE:  It is not thoroughly implemented!!!
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMDense object is not valid in inputTextMatrix()");
  #endif
  char   cNext;  // temporary character variable - good for checking line type
  const  int MaxLineLength = 1024;
  char   buffer[MaxLineLength];  // char* variable for getline() function

  if(fname=="")  
    errorMsg("No filename specified in TCMDense inputTextMatrix()!");

  nNodes = 0;  nEdges = 0;

  msg("Importing file: \""+fname,green);

  ifstream din(fname.c_str());
  if(din.bad()) {
    errorMsg("Input filename "+fname+" does not exist!");
    return -1;  // Error filename does not exist
  } // if din.bad
  cout << "  Reading text matrix data...\n";
  
  // Ascii input stuff here - skip comment lines that begin with a 'c' character
  bool bStart = 0;
  while(!din.eof() && !bStart) {
    //din >> sWord;                    // read first character of this line
    //cNext = (char)( sWord[0] );      // look at the first character a case test
    cNext = (char)( din.peek() );      // look at the first character a case test
    // need to properly read the initial info and skip it...
    // for the moment we ignore everything except node and edge designations
    if(cNext=='#')  din.getline(buffer,MaxLineLength); // skip rest of line
    else            bStart = 1;                        // begin data read
    cout << redText << "Comment line read as:  " << buffer << endl; // debugging
  } // end while !eof
  //din.getline(buffer,MaxLineLength); // get rest of first node line

  cout << magentaText << "Reading nodes... "; // debugging
  // Ascii input for file parameters
  // the number of nodes and edges as defined in the header of the file
  // used mostly as a consistency check
  int NInput = -1, LInput = -1;
  //din >> sWord;    // read first character of this line
  //din.getline(buffer,MaxLineLength); // get rest of line
  bStart = 0;
  while(!din.eof() && !bStart) {
    //din >> sWord;  // read first word of this line
    //cNext = (char)( sWord[0] );  // look at the first character a case test
    cNext = (char)( din.peek() );      // look at the first character a case test
    //cout << redText << sWord << " Edge read as: " << buffer 
    //     << " and cNext = " << cNext << " with nNodes = " 
    //     << nNodes << normText << endl; // debugging
    switch(cNext) {

      case 'N': { // number of nodes definition
        din.getline(buffer,MaxLineLength,'='); // get rest of line up to '='
        din >> NInput;  // count the current edge
        din.getline(buffer,MaxLineLength); // get rest of line after first char
        break;
      } // end case c

      case 'L': { // number of edges definition
        din.getline(buffer,MaxLineLength,'='); // get rest of line up to '='
        din >> LInput;  // count the current edge
        break;
      } // end case c

      default: {
        bStart = 1;  // now start matrix read
        break; // redundant break for default
      } // end default
    } // end switch for node input

  } // end while !eof

  cout << "Input text parameters are:  N = " << NInput 
       << " and L = " << LInput << endl;                           // debugging
  
  if(LInput==-1 || NInput ==-1)  
    errorMsg("Input text matrix is missing N or L value line(s)!  (Read as N = "
             +itos(NInput)+" and L = "+itos(LInput)+")\n");
    
  data.resize(NInput,NInput);
  //if(bDirected)  data.init(0);
  data.init(0);
  
  int         i, j=0;
  TCMData     element;
  TVectorInt  degreeCount(NInput,0);
  while(!din.eof() && j<NInput) {

    i = 0;
    while(!din.eof() && i<NInput) {
      din >> element;
      //cout << element << " ";  // debugging
      if(bDirected) { 
        data[i][j] += element; // sum the directed edge to symmetrize it
        data[j][i] += element;
      } else {
        data[i][j]  = element; // matrix is read in full, so only one assignment
      } // end else

      if(element>0) { nEdges++;  degreeCount[i] += 1; } // keep track of edges
      i++;
    } // end while !eof
    
    din.getline(buffer,MaxLineLength); // get rest of line after this data line
    //cout << endl;  // debugging
    j++;
  } // end while !eof

  // The above automatically symmetrizes the matrix data if bDirected.  This
  // presents a small problem since we use integer data and the matrix element
  // can be odd after the sum; therefore, we do not divide by two on the end.
  //if(bDirected)  data /= 2;  // correct for symmetrization on directed graph
  
  // some consistency checks
  nEdges /= 2;                   // correct for overcounted edges
  if(nEdges!=LInput)             
    warningMsg("Number of edges is inconsistent!  Counted as "+itos(nEdges)+
               " and input as L = "+itos(LInput));
  if(!(i==NInput && j==NInput))  errorMsg("Input text matrix is missing data!");
  nNodes = NInput;  Cols = NInput;

  // now assign the edge connect/degree data - I could make this more efficient
  // by incorpating the read above a vector at a time to get the size
  destroykij();
  createkij(nNodes);
  int edgeCount = 0;
  for(j=0; j<nNodes; j++) {
    kij[j].resize(degreeCount[j]);
    edgeCount = 0;
    for(i=0; i<nNodes; i++) 
      if(data[i][j]>0) { kij[j][edgeCount] = i;  edgeCount++; }
    // consistency check
    if(edgeCount!=degreeCount[j]) 
      errorMsg("Degree assignments are inconsistent?");
  } // end for j

  p = nEdges*2.0/( (double)nNodes*(double)(nNodes-1) );

  cout << "  N was read as " << nNodes << " and L as " << nEdges << ".  ";
  msg("done.\n");
  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif

  // successfully finished readinrg text matrix data file - close and exit
  din.close();
  return 1;
}; // end inputTextMatrix


void TCMDense::init(TCMData stored, TCMData unstored) {
  errorMsg("TCMDense::init(...) is not yet implemented");
}; // end init


void TCMDense::initVij(TMatrixFloat &rn, double t, TCMData a, TCMData b) {
  // We take a potential V(r_i, r_j) and use it to define a connection matrix
  // nxy is a 2D floating point vector that defines the position of all N nodes.
  //   nxy[node_i][0] = n_ix and nxy[node_i][1] = n_iy
  // Because we have a legacy integer data type for the matrix values, the
  // units are in microUnits.  That is, we scale all values by 10^6.
  int N = rn.getCols(), i, j;
  // we use the absolute values of the scales for the system definition
  b = abs(b);  a = abs(a);  
  nNodes = N;  Cols = N;
  nEdges = 0;

  // reset the edge list and array data structures
  destroykij();
  createkij(N);
  data.resize(N,N);
  data.init(0);
  double   vijmtVal, rnix, rniy;
  TVectorInt  degreeCount(N,0);  // track how many edges are added to each node
  for(i=0; i<N; i++) {
    rnix = rn[i][0];  rniy = rn[i][1];
    
    for(j=i+1; j<N; j++) {

      vijmtVal = Vij(rnix,rniy,rn[j][0],rn[j][1]) - t;
      if(-vijmtVal>0.0) {
        // we have a valid edge.  Include data and track information
        data[i][j] = (TCMData)( vijmtVal*a );
        data[j][i] = (TCMData)( vijmtVal*a );
        nEdges++;
        degreeCount[i] += 1;  // assumes symmetric right now
        degreeCount[j] += 1;
      } else {
        // we have a "missing" edge.  Do not track edge information.
        data[i][j] = (TCMData)( vijmtVal*b );
        data[j][i] = (TCMData)( vijmtVal*b );
      } // end else vijVal
      
    } // end for j
  } // end for i
  
  // now assign the edge connect/degree data
  //degreeCount /= 2;  // we do not double count above, so do not divide by 2
  int edgeCount = 0;
  for(j=0; j<nNodes; j++) {
    kij[j].resize(degreeCount[j]);
    edgeCount = 0;
    for(i=0; i<nNodes; i++) 
      if(data[i][j]>0) { kij[j][edgeCount] = i;  edgeCount++; }

    // consistency check
    if(edgeCount!=degreeCount[j]) 
      errorMsg("Degree assignments are inconsistent in initVij?  (count is "
               +itos(edgeCount)+" and stored is "+itos(degreeCount[j])
               +" for "+itos(j)+" and "+itos(i)+")");
  } // end for j

  p = (double)nEdges*2.0/( (double)nNodes*(double)(nNodes-1) );

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
  return;
}; // end initVij

//template <typename T> 
void TCMDense::display(int bOutputLarge, bool bShowZeroes) const {
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


int TCMDense::inputGML(string fname, int soffset, 
                              //bool bDirected, bool bSymmetrize) {
                              bool bDirected) {
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
  if(bDirected)  
    warningMsg("Symmetrizing directed matrix by b_ij = (a_ij + a_ji).");

  nNodes = 0;  nEdges = 0;

  cout << green << "Importing file: \"" << fname 
       << "\" assuming weighted edges " << normText << "\n" << flush;

  ifstream din(fname.c_str());
  if(din.bad()) {
    errorMsg("Input filename "+fname+" does not exist!");
    return -1;  // Error filename does not exist
  } // if din.bad
  cout << "  Reading GML data... ";
  
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

  cout << magentaText << "nodes... "; // debugging
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
          return -2;  // Corrupted GML data file
          break; // redundant break for default
        } // end default
      } // end switch for node input
  } // end while !eof

  bool  bEdgeDone = 0;
  // Ascii input stuff here
  TVectorInt  degrees(nNodes,0);
  data.resize(nNodes,nNodes);
  if(bDirected)  data.init(0);
  else { // normal symmetric matrix
    data.init(-1);
    for(int k=0; k<nNodes; k++)  data[k][k] = 0;  // clear diagonal
  } // end else

  cout << "edges... " << normText << endl; // debugging
  TVectorInt  degreeCount(nNodes,0);
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
        //cout << green << "Weight for edge " << nEdges << " connecting nodes " 
        //     << p1-soffset << " to " <<p2-soffset << " is " << weight << endl; // debugging
        if(bDirected) {
          // we automatically symmetrize a directed graph (for now)
          // note that this assumes a zero initialized data matrix from above
          if(data[p1-soffset][p2-soffset]==0) {
            nEdges++;
            degreeCount[p1-soffset] += 1;  // assumes symmetric right now
            degreeCount[p2-soffset] += 1;
          } // end if existing edge check for symmetrized matrices
          // now sum the weights of the edges.  Note that we do not average 
          // them because of integer rounding with the legacy integer weights.
          data[p1-soffset][p2-soffset] += weight;
          data[p2-soffset][p1-soffset] += weight;  
        } else { // else assume a regular symmetric matrix
          nEdges++;
          data[p1-soffset][p2-soffset]  = weight;
          data[p2-soffset][p1-soffset]  = weight;  
          degreeCount[p1-soffset] += 1;  // assumes symmetric right now
          degreeCount[p2-soffset] += 1;
        } // end else 
        weight = 1;     // return weight to default value for next read
        bEdgeDone = 0;
        p1 = -1;  p2 = -1;
        //cout << nEdges << " " << flush;  // debugging
      } // end if p1, p2, weight
    } // end while !eof
  //} // end if sWord
  
  // The above automatically symmetrizes the matrix data if bDirected.  This
  // presents a small problem since we use integer data and the matrix element
  // can be odd after the sum; therefore, we do not divide by two on the end.
  //if(bDirected)  data /= 2;  // correct for symmetrization on directed graph
  // because of the symmetrized assignements above, we now have to crudely set 
  // all non-edges to a value of -1
  if(bDirected) {
    for(int i=0; i<nNodes; i++)  
      for(int j=i+1; j<nNodes; j++) 
        if(data[i][j]==0) { data[i][j] = -1;  data[j][i] = -1; }
  } // end if bDirected

  // now assign the edge connect/degree data - I could make this more efficient
  // by incorpating the read above a vector at a time to get the size
  //cout << "Initializing kij in TCMDense inputGML with L = " << nEdges << endl; // debugging
  //cout << "Degree list:  " << degreeCount << endl;
  //cout << "nNodes = " << nNodes << endl;
  //degreeCount /= 2;  // temporary workaround for double counting with Mousumi
  if(bDirected)  warningMsg("Skipping degree check in TCMDense::inputGML()");
  destroykij();
  createkij(nNodes);
  int edgeCount = 0, i;
  for(int j=0; j<nNodes; j++) {
    kij[j].resize(degreeCount[j]);
    edgeCount = 0;
    for(i=0; i<nNodes; i++) 
      if(data[i][j]>0) { kij[j][edgeCount] = i;  edgeCount++; }
    // consistency check which we skip if importing a directed graph
    if(edgeCount!=degreeCount[j] && !bDirected)
      errorMsg("Degree assignments are inconsistent? (count is "+itos(edgeCount)
               +" and stored is "+itos(degreeCount[j])+" for "+itos(j)+" and "
               +itos(i)+")");
  } // end for j

  p = nEdges*2.0/( (double)nNodes*(double)(nNodes-1) );

  cout << "  The number of (TCMDense) GML nodes in was " << nNodes 
       << " and edges was " << nEdges << ".  " << flush;
  msg("done.\n");

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif

  // successfully finished reading GML data file - close and exit
  din.close();
  return 1;
}; // end inputGML


void TCMDense::initRandom(unsigned N, double randomDensity, 
                          TCMData  bottomWeight, TCMData topWeight, bool bN2) {
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
  data.resize(N,N);
  data.init(-1);
  //for(unsigned i=0; i<N; i++)  // ??? unsure of what this was for???

  //cout << red << "here c in initRandom " << normText << flush;   // debugging
  // loop over all pairs of nodes
  TCMData  wWeight;
  if(bN2) { // use N^2 version (for whatever reason)
    for(unsigned i=0; i<N; i++) 
      for(unsigned j=i+1; j<N; j++) 
        if(randomDouble()<randomDensity) { 
          wWeight = randomInt(bottomWeight,topWeight);
          data[i][j] = wWeight;  data[j][i] = wWeight;
          nEdges++; 
        } // end if
  } else { // use O(L) version
    // use a more efficient way O(L) rather than O(N^2)
    // number of expected edges
    int tEdges = (int)( randomDensity*(double)N*(double)(N-1)/2.0 + 0.5 );  
    int iN, jN, tMod = tEdges/1000;
    for(unsigned i=0; i<tEdges; i++) {
      if((i%tMod)==0) cout << "." << flush;  // update user on progress
      iN = randomInt(0,N-1);   jN = randomInt(0,N-1);
      // assumes symmetric
      while(iN==jN || data[iN][jN]<=0)  jN = randomInt(0,N-1);  // try again
      // should not have any duplicates at this point, so just increment nEdges
      nEdges++;   
      wWeight = randomInt(bottomWeight,topWeight);
      data[iN][jN] = wWeight;  data[jN][iN] = wWeight;  
    } // end for i
  } // end else bN2

  msg("done.");
  return;
}; // end initRandom


int  TCMDense::exportGML(string fname, bool bDirected) {
  warningMsg("exportGML() not yet implemented for TCMDense!");
};


int TCMDense::addNoiseExp(TClusterList& s, double pout, double kMean) {
  // loop over all edges and add edges with a probability of pout 
  // existing edges are ignorned.  Currently assumes a symmetric matrix.
  // s is the "signal" list - noise is not added *within* existing clusters
  // Cp is changed to include the specified level of additional noise between
  // clusters in s.  Note that this can be an expensive operation due to the
  // current structure of the CMatrix class.  The neighbor list matrix must be
  // re-written even for a dense 2D matrix structure.
  // nzout is the (average) number of additional edges to add to each node.
  // The noise edges are added randomly based on a resulting probability density
  // calculated based on nzout and the maximum number of possible new edges.
  int  nAdded = 0, N = nNodes, L = nEdges, qi, qj, iDegree, iDegreeMax;
  TListArray<int>  edges(N);
  double MaxK = (double)N/4.0;

  double pouti, douti;  // the connection probability for *this* node
  TVectorInt  degrees(N,0);
  for(int i=0; i<N; i++)  degrees[i] = ki(i);  // store this degree list
  // loop over all edges for a roughly exponential distribution
  for(int i=0; i<N; i++) {
    qi = s.nodes[i];
    // get external connection probability for node i by a random exponential
    // but we set a rejection level
    douti = 9.9E99;
    douti = -kMean*log(randomDouble()); // base e log, exponential variate
    if(douti>MaxK)  douti = MaxK;  // crude rejection test
    iDegreeMax = (int)( douti + 0.5 );

    // loop over all other nodes without repetition - We add external edges  
    // if the exponential degree is greater than the existing number of edges.
    // Thus it is not explicitly an exponential distribution.
    iDegree = degrees[i];
    if(iDegreeMax>iDegree) {
      pouti = (double)(iDegreeMax - iDegree)/(double)(N - s[qi].getn() - i);
      for(int j=i+1; j<N; j++) {
        qj = s.nodes[j];

        // check for adding a noise edge - existing edges are ignored
        if(data[i][j]<0 && qi!=qj && randomDouble()<pouti) {
          //cout << "Adding edge " << i << " to " << j << "\n";  // debugging
          edges[i].add(j);  edges[j].add(i);
          degrees[i] += 1;  degrees[j] += 1;
          nAdded++;
        } // end if
      } // end for j

      //degrees[i] = iDegree;  // write final value of degrees for node i
    } // end if
  } // end for i

  #ifdef DEBUG_MODE
  // check stats (not sure if we will modify output
  //if((double)nAdded/(double)N > zout)  
  //  cout << "Actual pout is " << ( (double)nAdded/(double)N ) 
  //       << " compared to the specified " << pout << endl;  // debugging
  #endif

  // now add the existing edges to the list for easy sorting and reconstruction
  int jNode;
  for(int i=0; i<N; i++) {
    iDegree = ki(i);
    for(int j=0; j<iDegree; j++)  edges[i].add(ki(i,j));
  } // end for i
  
  copyEdges(edges);

  return nAdded;
}; // end addNoiseExp


void TCMDense::initNoiseExp(unsigned N, double kMean, bool bWeighted, double wMean) {
  //  N - is the total size of the system
  //cout << red << "here a in initRandom " << normText << flush;   // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMDense object is not valid in initNoiseExp()");
  //msg("Assigning random communities...");
  cout << "kMean = " << kMean << ")" << endl;
  #endif

  destroykij();  // just in case - still valid if NULL
  if(N!=nNodes)  data.resize(N,N);
  data.init(-1);                           // init to unconnected
  for(int i=0; i<N; i++)  data[i][i] = 0;  // zero the diagonal
  nNodes = N;

  int iN, jN;
  nEdges = 0;
  #ifdef DEBUG_MODE
  cout << magenta << "Counting edges for degree distribution... " 
       << endl;  // debugging
  #endif
  TVectorInt  degrees(N,0);       // store a count of the number of edges
  TVectorInt  expDegrees(N,0);    // store target degree for each node
  TVectorInt  expdCount(N,0);     // count current target degree for each node
  // determine target degree for each node
  for(int i=0; i<N; i++)  expDegrees[i] = (int)( randomExp(kMean) + 0.5 );
  expdCount = expDegrees;
  TVectorInt  nodesLeft(N);       // which nodes have "completed" # of edges
  nodesLeft.initStep(0,1,1);
  int iLeft, nLeft = N;           // how many nodes are left to assign edges

  // assign external edges - assumes symmetric edges
  // this is not an efficient implementation
  // we fill the edges in order of decreasing degree to avoid any residual
  // unintended clustering
  //cout << "Adding exponential edges (nLeft = " << nodesLeft.getSize() 
  //     << ")... " << flush;  // debugging
  int  jLeft, iLeftMax, iVal;
  while(nodesLeft.getSize()>1) {
    nodesLeft.randomizeOrder();
    // get (current) max node degree
    iLeftMax = -1;  // start with invalid value
    for(int i=0; i<nodesLeft.getSize(); i++) {
      iVal = expdCount[nodesLeft[i]];
      if(iVal>iLeftMax) { iLeftMax = iVal;  iLeft = i; }
    } // end for i
    iN = nodesLeft[iLeft];  // get this node number

    // error checks
    if(iLeftMax<1 || iLeft<0 || iLeft>nodesLeft.getSize()) {
      cout << "nodesLeft vector is:  " << nodesLeft << "\n";
      cout << "expdCount values are:  ";
      for(int i=0; i<nodesLeft.getSize(); i++)
        cout << expdCount[nodesLeft[i]] << " ";
      cout << "\n";
      errorMsg("Current max node degree in initNoiseExp() is invalid ("
               +itos(iLeftMax)+")?");
    } // end if error check

    // Traverse nodes list in sequential order, but list is randomly ordered
    // Assign edges - start with zeroth node unless iLeft = 0
    int w = 1;  // default is unweighted edge
    if(iLeft==0) jLeft = 1;  else jLeft = 0;
    while(expdCount[iN]>0 && jLeft<nodesLeft.getSize()) {
      jN = nodesLeft[jLeft];

      if(data[iN][jN]<0) {
        if(bWeighted)  w = max((int)( randomExp(wMean) + 0.5 ), 1);
        degrees[iN]   += 1;  degrees[jN]   += 1;
        expdCount[iN] -= 1;  expdCount[jN] -= 1;
        data[iN][jN]   = w;  data[jN][iN]   = w;
        // check if we have "filled" the specified degree. If yes, drop the node
        if(expdCount[jN]<1) {
          if(iLeft==nodesLeft.getSize()-1) { // is iLeft the last is nodesLeft[]
            nodesLeft[jLeft] = nodesLeft[iLeft];
            iLeft = jLeft;
            nodesLeft.drop();  // decrement size without copying any elements
          } // end if
          else  nodesLeft.drop(jLeft);  // just drop the node as normal
        } // end if filled degree check
        //cout << "(" << iN << "," << jN << ") ";  // debugging
      } // end if

      jLeft++;  // data[iN][jN]<0 above catches special case of iN == jN
    } // end while jLeft

    // drop node at position 0 if it is empty or we cannot satisfy the power
    // law degree distribution (the only two possible cases here)
    nodesLeft.drop(iLeft);  // get rid of this node now even if not filled
  } // end while iLeft
  //cout << expDegrees << "\n" << expdCount << endl;  // debugging

  // debugging
  //cout << "Generate degrees excess... " << endl;  // debugging
  double zLeftSum = 0.0;
  TStats1D  zLeftVals("zLeft");
  TVectorInt  degreesLeft(N);
  for(int i=0; i<N; i++) {
    degreesLeft[i] = expDegrees[i] - degrees[i];
    zLeftSum      += degreesLeft[i];
    zLeftVals.store(degreesLeft[i]);
  } // end for i
  //cout << cyan << degreesLeft << " with a sum of " << zLeftSum 
  //     << " and avgerage of " << (zLeftSum/(double)N) << "\n" 
  //     << brown << expDegrees << "\n"
  //     << magenta << degrees << normText << endl;  // debugging
  // debugging
  cout << cyan << "degreesLeft has a sum of " << zLeftSum << ", an average of " 
       << brown << zLeftVals.avg() << cyan << ", and a sigma of " 
       << magenta << zLeftVals.sigma() << normText << "\n";

  // now assign the neighbors of all nodes
  copyEdges(degrees);
  
  return;
}; // end initNoiseExp


int TCMDense::addNoisePower(TClusterList& s, double kMin, double kMax, double alpha) {
  // loop over all edges and add edges with a probability of pout 
  // existing edges are ignorned.  Currently assumes a symmetric matrix.
  // s is the "signal" list - noise is not added *within* existing clusters
  // Cp is changed to include the specified level of additional noise between
  // clusters in s.  Note that this can be an expensive operation due to the
  // current structure of the CMatrix class.  The neighbor list matrix must be
  // re-written even for a dense 2D matrix structure.
  // nzout is the (average) number of additional edges to add to each node.
  // The noise edges are added randomly based on a resulting probability density
  // calculated based on nzout and the maximum number of possible new edges.
  int  nAdded = 0, N = nNodes, L = nEdges, qi, qj, iDegree, iDegreeMax;
  TListArray<int>  edges(N);
  //cout << red << "kMin = " << kMin << "  kMax = " << kMax 
  //     << "  alpha = " << alpha << "\n" << normText;  // debugging
  //TStats1D rs(N,0), ps(N,0), ds(N,0), ms(N,0);

  double pouti;  // the connection probability for *this* node
  TVectorInt  degrees(N,0);
  for(int i=0; i<N; i++)  degrees[i] = ki(i);  // store this degree list
  // loop over all edges for a roughly power-law distribution
  for(int i=0; i<N; i++) {
    qi = s.nodes[i];
    // get degree size for node i by a random power-law, but set rejection level
    iDegreeMax = (int)( randomPower(kMin,kMax,alpha) + 0.5 );

    // loop over all other nodes without repetition - We add external edges  
    // if the power-law degree is greater than the existing number of edges.
    // Thus it is not explicitly a power-law distribution.
    iDegree = degrees[i];
    //ds.store(iDegree);  // debugging
    if(iDegreeMax>iDegree) {
      pouti = (double)(iDegreeMax - iDegree)/(double)(N - s[qi].getn());
      //ps.store(pouti);  // debugging
      for(int j=0; j<N; j++) {
        qj = s.nodes[j];

        // check for adding a noise edge - existing edges are ignored
        if(data[i][j]<0 && qi!=qj && randomDouble()<pouti) {
          //cout << "Adding edge " << i << " to " << j << "\n";  // debugging
          edges[i].add(j);  edges[j].add(i);
          degrees[i] += 1;  degrees[j] += 1;  // assumes symmetric
          // data assignments are not sufficient, so we use copyEdges below
          // these just stop overwrites of same elements
          data[i][j]  = 1;  data[j][i]  = 1;  // assumes symmetric
          nAdded++;
        } // end if
      } // end for j
    } // end if
  } // end for i

  #ifdef DEBUG_MODE
  // check stats (not sure if we will modify output
  //if((double)nAdded/(double)N > zout)  
  //  cout << "Actual pout is " << ( (double)nAdded/(double)N ) 
  //       << " compared to the specified " << pout << endl;  // debugging
  #endif

  // now add the existing edges to the list for easy sorting and reconstruction
  int jNode;
  for(int i=0; i<N; i++) {
    iDegree = ki(i);
    for(int j=0; j<iDegree; j++)  edges[i].add(ki(i,j));
  } // end for i
  
  copyEdges(edges);

  return nAdded;
}; // end addNoisePower


void TCMDense::initByClustersPower(TClusterList &c, unsigned N, 
       double pin, double kMin, double kMax, double alpha, 
       int &nEdgesIn, int &nEdgesOut) {
  //  N - is the total size of the system
  //cout << red << "here a in initRandom " << normText << flush;   // debugging
  #ifdef DEBUG_MODE
  if(!isValid())  errorMsg("TCMDense object is not valid in initRandom()");
  //msg("Assigning random communities...");
  cout << "Assigning random communities with a seed of " << rngSeedString()
       << "\nBeginning kMin " << kMin << ") and kMax = " << kMax << endl;
  #endif
  
  destroykij();  // just in case - still valid if NULL
  if(N!=nNodes)  data.resize(N,N);
  data.init(-1);                           // init to unconnected
  for(int i=0; i<N; i++)  data[i][i] = 0;  // zero the diagonal
  nNodes = N;
  nEdgesIn = 0;  nEdgesOut = 0;

  int iN, jN;
  nEdges = 0;
  #ifdef DEBUG_MODE
  cout << magenta << "Counting edges for degree distribution... " 
       << endl;  // debugging
  #endif
  TVectorInt  degrees(N,0);       // store a count of the number of edges
  TVectorInt  powerDegrees(N,0);  // store target degree for each node
  TVectorInt  powerdCount(N,0);   // count current target degree for each node
  for(int i=0; i<N; i++)          // determine target degree for each node
    powerDegrees[i] = (int)( randomPower(kMin,kMax,alpha) + 0.5 );
  powerdCount = powerDegrees;
  TVectorInt  nodesLeft(N);       // which nodes have "completed" # of edges
  nodesLeft.initStep(0,1,1);
  int iLeft, nLeft = N;           // how many nodes are left to assign edges

  // assign external edges - assumes symmetric edges
  // this is not an efficient implementation
  //cout << "Adding power law edges (nLeft = " << nodesLeft.getSize() 
  //     << ")... " << flush;  // debugging
  int  jLeft, qi, qj, iLeftMax, iVal;
  while(nodesLeft.getSize()>1) {
    nodesLeft.randomizeOrder();
    // get (current) max node degree 
    iLeftMax = -1;  // start with invalid value
    for(int i=0; i<nodesLeft.getSize(); i++) {
      iVal = powerdCount[nodesLeft[i]];
      if(iVal>iLeftMax) { iLeftMax = iVal;  iLeft = i; }
    } // end for i
    iN = nodesLeft[iLeft];  // get this node number
    qi = c.nodes[iN];

    // error checks
    if(iLeftMax<1 || iLeft<0 || iLeft>nodesLeft.getSize()) {
      cout << red << "Warning:" << lightred << "  nodesLeft vector is:  " 
           << nodesLeft << "\npowerdCount values are:  " << normText;
      for(int i=0; i<nodesLeft.getSize(); i++)  
        cout << powerdCount[nodesLeft[i]] << " ";
      cout << "\n";
      errorMsg("Current max node degree in initByClustersPowerRS() is invalid ("
               +itos(iLeftMax)+")?");
    } // end if error check

    // At this point, all nodes iN have are expected to add at least some edges
    // (i.e. powerdCount[iN]>0).
    //cout << "Starting iN = " << iN << " (iLeft = " << iLeft 
    //     << "):  expected edges = " << powerdCount[iN]
    //     << " with nLeft = " << nodesLeft.getSize() << endl;  // debugging

    // Traverse nodes list in sequential order, but list is randomly ordered
    // Assign edges - start with zeroth node unless iLeft = 0
    if(iLeft==0) jLeft = 1;  else jLeft = 0;
    while(powerdCount[iN]>0 && jLeft<nodesLeft.getSize()) {
      jN = nodesLeft[jLeft];
      qj = c.nodes[jN];
      //cout << "Trying jN = " << jN << " (jLeft = " << jLeft 
      //     << "):  expected edges = " << powerdCount[jN] 
      //     << " with nLeft = " << nodesLeft.getSize() << endl;  // debugging

      if(qi!=qj && data[iN][jN]<0) {
        degrees[iN] += 1;      degrees[jN] += 1;
        data[iN][jN] = 1;      data[jN][iN] = 1;
        powerdCount[iN] -= 1;  powerdCount[jN] -= 1;
        nEdgesOut++;
        // check if we have "filled" the specified degree. If yes, drop the node
        if(powerdCount[jN]<1) {
          if(iLeft==nodesLeft.getSize()-1) { // is iLeft the last is nodesLeft[]
            nodesLeft[jLeft] = nodesLeft[iLeft];
            iLeft = jLeft;
            nodesLeft.drop();  // decrement size without copying any elements
          } // end if
          else  nodesLeft.drop(jLeft);  // just drop the node as normal
        } // end if filled degree check
        //cout << "(" << iN << "," << jN << ") ";  // debugging
      } // end if

      jLeft++;  // data[iN][jN]<0 above catches special case of iN == jN
    } // end while jLeft

    // drop node at position 0 if it is empty or we cannot satisfy the power
    // law degree distribution (the only two possible cases here)
    nodesLeft.drop(iLeft);  // get rid of this node now even if not filled
  } // end while iLeft
  //cout << powerDegrees << "\n" << powerdCount << endl;  // debugging

  // debugging
  //cout << "Generate degrees excess... " << endl;  // debugging
  double zLeftSum = 0.0;
  TStats1D  zLeftVals("zLeft");
  TVectorInt  degreesLeft(N);
  for(int i=0; i<N; i++) {
    degreesLeft[i] = powerDegrees[i] - degrees[i];
    zLeftSum      += degreesLeft[i];
    zLeftVals.store(degreesLeft[i]);
  } // end for i
  //cout << cyan << degreesLeft << " with a sum of " << zLeftSum 
  //     << " and avgerage of " << (zLeftSum/(double)N) << "\n" 
  //     << brown << powerDegrees << "\n"
  //     << magenta << degrees << normText << endl;  // debugging
  // debugging
  cout << cyan << "degreesLeft has a sum of " << zLeftSum << ", an average of " 
       << brown << zLeftVals.avg() << cyan << ", and a sigma of " 
       << magenta << zLeftVals.sigma() << normText << "\n";

  //cout << "internal edges... " << flush;  // debugging
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
          // create a new edge if it does not already exist
          if(data[iN][jN]<0) {
            // these edges may "overfill" the degrees of nodes, but our 
            // benchmark requires a specified average density
            degrees[iN] += 1;      degrees[jN] += 1;
            data[iN][jN] = 1;      data[jN][iN] = 1;
            powerdCount[iN] -= 1;  powerdCount[jN] -= 1;
            nEdgesIn++;
          } // end if
        } // end if
        else if(data[iN][jN]>0) {
          // remove existing edge so that we maintain the proper ratio 
          // of internal connected edges based on the density pin
          degrees[iN] -= 1;      degrees[jN] -= 1;
          data[iN][jN] = -1;     data[jN][iN] = -1;
          powerdCount[iN] += 1;  powerdCount[jN] += 1;
          nEdgesIn--;
        } // end if
      } // end for j
    } // end for i
  } // end for k


  #ifdef DEBUG_MODE
  cout << "Allocating degree data... " << flush;
  //cout << endl << degrees << endl << " with i = " << flush;
  #endif
  // create kij with the known correct sizes - an N^2 operation as implemented
  createkij(N);  // allocate the rows (columns)
  int dCount;
  nEdges = 0;
  for(unsigned i=0; i<N; i++) {
    // memory allocation
    kij[i].resize((unsigned)degrees[i]);

    // add the new edges to the degree list
    dCount = 0;
    for(unsigned j=0; j<N; j++) {
      //cout << data[i][j] << " ";  // debugging
      if(data[i][j]>0) { kij[i][dCount] = j;  dCount++; }
    } // end for j
    nEdges += dCount;

    if(dCount!=degrees[i])
      errorMsg("Degree error in TCMDense::initByClustersPower()?  (Counted "
              +itos(dCount)+" edges with "+itos(degrees[i])+" edges expected)");
    //cout << kij[i] << endl;  // debugging
  } // end for i
  nEdges /= 2; // correct for double counting
  
  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1));  // calculate actual density

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  msg("done.");
  #endif
  //cout << "done with init." << endl;
  
  return;
}; // end initByClustersPower


void TCMDense::copyEdges(Array1D<int> &degrees) {
  // create the degree list from a degree vector and the existing data matrix
  #ifdef DEBUG_MODE
  cout << "Allocating degree data... " << flush;
  //cout << endl << degrees << endl << " with i = " << flush;
  #endif
  int N = degrees.getSize();

  // create kij with the known correct sizes - an N^2 operation as implemented
  createkij(N);  // allocate the rows (columns)
  int dCount;
  nEdges = 0;    // use class variable as a summing variable
  for(unsigned i=0; i<N; i++) {
    // memory allocation
    kij[i].resize((unsigned)degrees[i]);

    // add the new edges to the degree list
    dCount = 0;
    for(unsigned j=0; j<N; j++) {
      //cout << data[i][j] << " ";  // debugging
      if(data[i][j]>0) { kij[i][dCount] = j;  dCount++; }
    } // end for j
    nEdges += dCount;

    if(dCount!=degrees[i])
      errorMsg("Degree error in TCMDense::copyEdges()?  (Counted "
              +itos(dCount)+" edges with "+itos(degrees[i])+" edges expected)");
    //cout << kij[i] << endl;  // debugging
  } // end for i
  nEdges /= 2; // correct for double counting
  
  Cols = N;
  p = 2.0*(double)nEdges/((double)N*(double)(N-1));  // calculate actual density

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  msg("done.");
  #endif
  //cout << "done with init." << endl;
  return;
}; // end copyEdges for TCMDense matrix


void TCMDense::initByCMArray(Array2D<int> &a, int N, string d) {
  // initialize the connection matrix with weighted array
  //    'a' should be a symmetric integer matrix
  //    a[i][j] is positive - defines and edge and its weight
  //    a[i][j] is negative - defines a missing edge
  //    a[i][j] is zero     - error except for a[i][i] diagonal, no self edges
 
  #ifdef DEBUG_MODE
  debugMsg("Initializing TCMDense object from Array2D<int> object... ");
  #endif
  if(a.getRows()!=a.getCols())
    errorMsg("Connection matrix is not square in initByCMArray().");

  //cout << "Size of parameter matrix a is " << N << endl;// debugging
  description = d; // assign description string
  //int N = a.getCols();
  if(a.getRows()!=N || a.getCols()!=N || data.getCols()!=N) {
    //cout << "Resizing data matrix... " << flush;// debugging
    data.resize(N,N);
    //cout << "done." << flush;// debugging
  } // end if
  //cout << "Copying data matrix... " << flush;// debugging
  nNodes = N;
  data   = a;  // declare the data matrix
  //cout << "done." << flush;// debugging
  TVectorInt degrees(N,0);  // track assigned degrees according to 'a'

  // scan the connection matrix to determine the degrees of each node
  //cout << "Counting node neighbors in initByCMArray()... " << flush;// debugging
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      if(data[i][j]>0) { degrees[i] += 1;  degrees[j] += 1; }

      // catch some errors
      if(data[i][j]!=data[j][i])
        errorMsg("Array is not symmetric in initByCMArray().");
      if(data[i][j]==0 || data[j][i]==0)
        errorMsg("Array has zero off-diagonal entry in initByCMArray().");
    } // end for j

    if(data[i][i]!=0)
      errorMsg("Array has non-zero diagonal entry in initByCMArray().");
      #ifdef DEBUG_MODE
      if(degrees[i]==0)  
        errorMsg("Node "+itos(i)+" has no neighbors in initByCMArray()");
      #else
      //if(degrees[i]==0)  
        //warningMsg("Node "+itos(i)+" has no neighbors in initByCMArray()");
      #endif
  } // end for i
  //cout << "done with counting node neighbors.\n" << flush;// debugging

  // now define the neighbor list
  //cout << "Copying edges in initByCMArray()... " << flush;// debugging
  copyEdges(degrees);
  //cout << "done with copying edges.\n" << flush;// debugging
  return;
}; // end init with connection matrix


int TCMDense::addNoisep(TClusterList& s, double pout) {
  // loop over all edges and add edges with a probability of pout 
  // existing edges are ignorned.  Currently assumes a symmetric matrix.
  // s is the "signal" list - noise is not added *within* existing clusters
  // Cp is changed to include the specified level of additional noise between
  // clusters in s.  Note that this can be an expensive operation due to the
  // current structure of the CMatrix class.  The neighbor list matrix must be
  // re-written even for a dense 2D matrix structure.
  // nzout is the (average) number of additional edges to add to each node.
  // The noise edges are added randomly based on a resulting probability density
  // calculated based on nzout and the maximum number of possible new edges.
  int  nAdded = 0, N = nNodes, L = nEdges, qi, qj;
  TListArray<int>  edges(N);

  // loop over all edges
  for(int i=0; i<N; i++) {
    qi = s.nodes[i];

    // loop over all other nodes without repetition
    for(int j=i+1; j<N; j++) {
      qj = s.nodes[j];

      // check for adding a noise edge - existing edges are ignored
      if(data[i][j]<0 && qi!=qj && randomDouble()<pout) {
        //cout << "Adding edge " << i << " to " << j << "\n";  // debugging
        edges[i].add(j);  edges[j].add(i);
        nAdded++;
      } // end if
    } // end for j
  } // end for i
  
  #ifdef DEBUG_MODE
  // check stats (not sure if we will modify output
  //if((double)nAdded/(double)N > zout)  
  //  cout << "Actual pout is " << ( (double)nAdded/(double)N ) 
  //       << " compared to the specified " << pout << endl;  // debugging
  #endif

  // now add the existing edges to the list for easy sorting and reconstruction
  int iDegrees, jNode;
  for(int i=0; i<N; i++) {
    iDegrees = ki(i);
    for(int j=0; j<iDegrees; j++)  edges[i].add(ki(i,j));
  } // end for i
  
  copyEdges(edges);

  return nAdded;
}; // end addNoisep


//template <typename T> 
void TCMDense::initByClusters2(TClusterList &c, int N, double pin, double pout,
     int &nEdgesIn, int &nEdgesOut, TCMData  bottomW, TCMData topW) {
  // This structure poses a more difficult identification than the symmetric
  // hiearchical structure.
  //  bRandomNodes - specifies whether the nodes will the randomized or 
  //      sequentially assigned.
  #ifdef DEBUG_MODE
  if(!isValid()) 
    errorMsg("TCMSparseW object is not valid in initByClusters2()");
  #endif
  if(c.getNNodes()==0)  
    errorMsg("In initByClusters2(), cluster list says it has zero nodes?");

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
          edges[iS].add(jS);  edges[jS].add(iS);
          nEdges++;           nEdgesIn++;
        } // end if
      } // end for j

      // loop over all other clusters m (past k) for outside connections
      for(int m=k+1; m<c.getSize(); m++) {
        // loop over nodes j in cluster m
        for(int j=0; j<c[m].getSize(); j++) {
          jS = c[m].getNode(j);
          if(randomDouble()<pout) { 
            edges[iS].add(jS);  edges[jS].add(iS);
            nEdges++;           nEdgesOut++; 
          } // end if
        } // end for j
      } // end for m

    c[k].next();
    } // end for i
  } // end for k
  
  // edges.sort(); is broken!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //for(int i=0; i<N; i++)  edges[i].sort(); // debugging
  //cout << "edges list:  " << edges << flush;  // debugging
  for(int i=0; i<N; i++)  
    if(edges[i].getSize()==0) { 
      c.display(0,1);
      errorMsg("Node "+itos(i)+" in cluster "+itos(c.nodes[i])
               +" has zero edges in initByClusters2() with pin = "
               +ftos(pin)+" and pout = "+ftos(pout)+"?");
    } // end if
  copyEdges(edges);
  //cout << brown << "done." << flush;  // debugging
  
  //msg("done\n");
  return;
}; // end initByClusters2
void TCMDense::initByClusters(TClusterList &c, int N, double pin, double pout,
     TCMData  bottomWeight, TCMData topWeight) {
  // a version that does not require the reference edge return values
  int ignore1, ignore2;
  initByClusters2(c,N,pin,pout,ignore1,ignore2,bottomWeight,topWeight); 
}; // end non-reference version



int TCMDense::copyEdges(TListArray<int> &edges) {
  // copy a TListArray to the kij edge list and update class variables
  // does not check for duplicates at the moment
  // return value is nonsense at the moment
  int N = edges.getSize(), jNode;

  destroykij();
  createkij(N);
  data.resize(N,N);
  data.init(-1);  // the default is unconnected and unweighted
  for(int i=0; i<N; i++)  data[i][i] = 0;  // no self edges
  
  //cout << "a: i = " << flush; // debugging
  int iDegree;
  nEdges = 0;
  for(int i=0; i<N; i++) {
    //cout << i << " " << flush; // debugging
    iDegree = edges[i].getSize();
    edges[i].sort();
    kij[i].resize(iDegree);
    nEdges += iDegree;
    
    edges[i].begin();  // begin manual iteration
    for(int j=0; j<iDegree; j++) {
      jNode = edges[i].getCurrent();
      kij[i][j] = jNode;
      data[i][jNode] = 1;  // assume unweighted edge at the moment
      edges[i].next();
    } // end for j
  } // end for i
  //cout << "a " << flush; // debugging
  
  // new density
  nEdges /= 2;  // correct for double counted edges
  nNodes = N;  Cols = N;  // just in case
  p = (double)nEdges*2.0/( (double)N*(double)(N-1));

  #ifdef DEBUG_MODE
  ValidArray = ValidCMatrixFlag;
  #endif
  return 0;  
}; // end copyEdges for dense structure


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Overloaded iostream operators
//template <typename T>
ostream& operator<<(ostream& fout, TCMDense& a) {
  a.display(1,1);
  return fout;
};

// distance cutoff functions
// need to declare this before the dense functions
double Vij(double x1, double y1, double x2, double y2) {
  //specify form of V_ij, the generalized particle interaction function
  double vijVal;
  //vijVal = ...???
  return vijVal;  // function here 
}


} // end namespace MSL
//#endif


