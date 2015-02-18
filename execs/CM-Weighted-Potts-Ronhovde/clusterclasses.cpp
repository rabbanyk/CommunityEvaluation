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
January 2008

Specific Notes:
As a science application, we do *not* perform extraneous redundant checks for
a valid input parameter.

To-Do List:
1. 
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//#define STRICT_ERROR_CHECKING_P // uncomment when there are problems
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#ifndef ML_UTILS_H
#include "./msl/ML_Utils.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "./msl/MSL_Array1D_Template.h"
#endif
#ifndef ML_TCLUSTERLIST_H
#include "clusterclasses.h"
#endif

using namespace std;
using namespace MSL;

namespace MSL {

//typedef TBool unsigned short;  // define my own boolean type

//------------------------------------------------------------------------------
//--------------------- TNode Member Function Declarations----------------------
TNode::TNode(int s) { 
  node = s;      clean = 0;       flag = 0;      //pCl = NULL;  //fixed = 0; 
};
TNode::~TNode() {};  // end TNode // no action needed
// following integer assignments clears the clean parameter???
//TNode&  TNode::operator=(int   s) { 
//  node = s;    clean = 1;       flag = 0;      pCl = NULL;  //fixed = 0;    
//  return *this;
//};
TNode&  TNode::operator=(const TNode &s) {
  node = s.node; clean = s.clean; flag = s.flag; //pCl = s.pCl; //fixed = s.fixed;
  return *this;
};
void    TNode::init(int s, TBool b, TBool f) {  //, TBool fi = 0) { 
  node = s;      clean = b;       flag = f;      //pCl = NULL;  //fixed = fi;      
  return;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------- TClusterNode Member Function Declarations-----------------
TClusterNode::TClusterNode(int s) { 
  n.node = s;  n.clean = 0;  n.flag = 0;  //n.fixed = 0; 
  //n.unsetSuperNode();
  pNext = NULL;
};
TClusterNode::~TClusterNode() {};  // no action is needed

/*
TClusterNode&  TClusterNode::operator=(const TNode &s) {
  n.node = s.node;  n.clean = s.clean;  n.flag = s.flag;  //n.fixed = s.fixed;
  n.setSuperNode(s.superNodeCluster());
  pNext = NULL;
  return *this;
};
*/

//void  TClusterNode::init(int s, TBool c, TBool f, TCluster *pc, TClusterNode *pn) {
void  TClusterNode::init(int s, TBool c, TBool f, TClusterNode *pn) {
  n.node = s;  n.clean = c;  n.flag = f;  //fixed = 0;    
  //n.setSuperNode(*pc);
  pNext = pn;
  return;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------- TCluster Member Function Declarations-----------------
TCluster::TCluster() {
  pBegin    = NULL;  pEnd = NULL;
  nNodes    = 0;     nEdges = 0;  //nMaxEdges = 0;
  energy    = 0;     partialQ = 0.0;
  pCurrent  = NULL;  pPrevious = NULL;  // user pointer iteration parameters
  piCurrent = NULL;  iCurrent = 0;      // operator[] parameters
  bFlag     = 0;
  //expectedEnergy = 0.0;
};  // default constructor
/*
TCluster::TCluster(int nodeList[], int size) {
  pBegin = NULL;     pEnd = NULL;
  nNodes = size;     nEdges = 0;  //nMaxEdges = 0;  
  energy = 0;        partialQ = 0.0;
  pCurrent  = NULL;  pPrevious = NULL;  // user pointer iteration parameters
  piCurrent = NULL;  iCurrent = 0;      // operator[] parameters
  if(size>0) {
    pBegin = new TClusterNode;
    pBegin->init(nodeList[0]);
    TClusterNode *pn = pBegin;
    for(int i=1; i<size; i++) {
      pn->pNext = new TClusterNode;
      pn->init(nodeList[i]);
      pn = pn->pNext;
    } // end for i
    pEnd = pn; // finally set end pointer to last nodes in list
  } else { 
    pBegin = NULL;  pEnd = NULL; 
    warningMsg("Attempting to construct cluster with size zero array");
  } // end else
};  // integer list constructor
*/

TCluster::~TCluster() {
  TClusterNode *pn, *pc = pBegin;
  // traverse and delete linked list
  for(int i=0; i<nNodes; i++) { 
    pn = pc->pNext;  
    #ifdef DEBUG_MODE
    if(pc==NULL)  errorMsg("~TCluster has premature terminating NULL pointer?");
    #endif
    delete pc;
    pc = pn; 
  } // end for i
  //pBegin = NULL;    pEnd = NULL;       // just in case
  //pCurrent = NULL;  pPrevious = NULL;  // just in case
  #ifdef DEBUG_MODE
  if(pc!=NULL)  errorMsg("~TList<T> exits on non-NULL pointer?");
  #endif
};  // destructor


void TCluster::initClean() {
  if(nNodes>0) {
    TClusterNode *pn = pBegin;
    for(int i=0; i<nNodes; i++) { pn->setClean();  pn = pn->pNext; }
  } // end if nNodes
  return;
};  // initClean

void TCluster::clearFlags() {
  if(nNodes>0) {
    TClusterNode *pn = pBegin;
    for(int i=0; i<nNodes; i++) { pn->n.flag = 0;  pn = pn->pNext; }
  } // end if nNodes
  bFlag = 0;
  return;
};  // clearFlags


int TCluster::calckSum(TCMatrix &Cp) const {
  int iNode, dSum = 0;  // the degree summing variable
  if(nNodes>0) {
    TClusterNode *pn = pBegin;
    for(int i=0; i<nNodes; i++) { 
      iNode = pn->n.node;
      #ifdef DEBUG_MODE
      if(iNode>=Cp.getSize() || iNode<0)  
        errorMsg("Invalid value of iNode = "+itos(iNode)+" in calckSum()?");
      #endif
      dSum += Cp.ki(iNode);  
      pn = pn->pNext; 
    } // end for i
    #ifdef DEBUG_MODE
    if(pn!=NULL)  errorMsg("Invalid non-NULL terminator in calckSum()?");
    #endif
  } // end if nNodes
  return dSum;
}; // end calckSum


int TCluster::calcEdges(TCMatrix &Cp, int jNode) const {
  // assumes symmetric Cp
  int jEdges = 0, iNode;
  if(nNodes>0) {
    TClusterNode *pn = pBegin;
    for(int i=0; i<nNodes; i++) {
      //cout << "Trying iNode = " << flush;  // debugging
      iNode = pn->n.node;
      //cout << iNode << endl;  // debugging
      #ifdef DEBUG_MODE
      if(iNode>=Cp.getSize() || iNode<0)  
        errorMsg("Invalid value of iNode = "+itos(iNode)+" in calcEdges()?");
      #endif
      //cout << "Trying Cp(" << iNode << "," << jNode << ") = " << flush;  // debugging
      if(Cp(iNode,jNode)>0)  jEdges++;  
      //cout << Cp(iNode,jNode) << endl;  // debugging
      //cout << "Stepping to next node... " << flush;  // debugging
      pn = pn->pNext; 
      //cout << "done." << endl;  // debugging
    } // end for i
    #ifdef DEBUG_MODE
    if(pn!=NULL)  errorMsg("Invalid non-NULL terminator in calcEdges()?");
    #endif
  } // end if nNodes
  return jEdges;
}; // end calcEdges


void TCluster::copyToVector(TVectorInt &v) {
  v.resize(nNodes);
  if(nNodes>0) {
    TClusterNode *pn = pBegin;
    for(int i=0; i<nNodes; i++) { v[i] = pn->n.node;  pn = pn->pNext; }
  } // end if nNodes
  return;
};  // copy to vector int


void TCluster::init(int nodeList[], int size) {
  nNodes   = size;  nEdges    = 0;  //nMaxEdges = 0;  
  energy   = 0;     partialQ  = 0.0;
  pCurrent = NULL;  pPrevious = NULL;
  if(size>0) {
    pBegin = new TClusterNode;
    pBegin->init(nodeList[0]);
    pEnd = pBegin;
    TClusterNode *pn = pBegin;
    for(int i=1; i<size; i++) {
      pn->pNext = new TClusterNode;
      pn->init(nodeList[i]);
      pn = pn->pNext;
    } // end for i
    pEnd = pn; // finally set end pointer to last nodes in list
  } else { 
    pBegin = NULL;  pEnd = NULL; 
    warningMsg("Attempting to initialize cluster with size zero array");
  } // end else
  bFlag = 0;

  return;
};  // init

void  TCluster::init(TVectorInt &v, TBool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  // the integers can be general
  int N = v.getSize();  // use given size of v
  TVectorInt *pv = &v;
  
  // now randomize the node order
  if(bRandomOrder) {
    TVectorInt u(N);  
    u = v;
    u.randomizeOrder();
    pv = &u;   // re-assign pointer to randomized object
  } // end if bRandomOrder

  // now finally add the vector to the cluster
  for(int k=0; k<N; k++)  add((*pv)[k]);
  //cout << s << endl; // debugging
  bFlag = 0;

  return;
}; // end init


TCluster& TCluster::operator=(const TCluster& c) {
  if(nNodes>0)  erase();  // delete current data
  nNodes   = c.nNodes;  nEdges    = c.nEdges;  //nMaxEdges = 0;
  energy   = c.energy;  partialQ  = c.partialQ;
  pPrevious = NULL;
  bFlag = c.bFlag;

  if(c.nNodes>0) {
    TClusterNode  *pc = c.pBegin;
    pBegin      = new TClusterNode;
    pBegin->n   = c.pBegin->n;          // copy node data
    TClusterNode  *pn = pBegin;
    for(int i=1; i<c.nNodes; i++) {
      pn->pNext = new TClusterNode;
      pn = pn->pNext;  pc = pc->pNext;  // step to next node on both
      pn->n = pc->n;                    // copy node data
    } // end for i
    pn->pNext = NULL;
    pEnd      = pn;
    piCurrent = pBegin;
  } // end if nNodes
  else { pBegin = NULL;  pEnd = NULL; }
  
  return *this;
};  // equals operator


TNode&  TCluster::operator[](int i) {
  // an array[] reference operator for the linked-list implementation
  // obviously, it is not efficient if i>0.
  // Use the iteration pointers for an array-like iterative implementation 
  // as long as the iteration is relatively sequential through the list.
  // Random data access will be O(n).
  // Things can get messed up if there are deletions at the point of the 
  // iCurrent pointer which often occurs on a node move... be careful.
  // I did not want to bury more special cases in the rest of the code.
  // It is best used in a closed short loop starting at node zero.
  //if(nNodes==0)  errorMsg("cluster op["+itos(i)+"] is a size zero cluster?");
  //cout << "iCurrent = " << iCurrent << endl;  // debugging
  if(i>=nNodes)  errorMsg("cluster op["+itos(i)+"] reference is out of bounds");
  else if(i==0)            { piCurrent = pBegin;            iCurrent = 0; }
  // up until the final case, we are assuming now that piCurrent is valid
  else if(iCurrent==i-1)   { piCurrent = piCurrent->pNext;  iCurrent++;   }
  else if(iCurrent==i)     {} // do nothing here
  else if(iCurrent<i)  
         while(iCurrent<i) { piCurrent = piCurrent->pNext;  iCurrent++;   }
  // now we just start over which is much safer
  else { 
    // since we start over, this case catches all other cases including iCurrent
    // out of bounds
    iCurrent = 0;  piCurrent = pBegin;
    while(iCurrent<i)      { piCurrent = piCurrent->pNext;  iCurrent++;   }
  } // end else
  // do some error detection prior to exiting - this still will not catch most 
  // cases of moved or deleted nodes at the point of iCurrent however
  //#ifdef DEBUG_MODE
  if(pBegin==NULL)
    errorMsg("Cluster operator["+itos(i)+"] detected pBegin==NULL!? " +
             "(with iCurrent = "+itos(iCurrent)+")");
  else if(piCurrent==NULL)  
    errorMsg("Cluster operator["+itos(i)+"] is returning a NULL!?  " + 
             "(with iCurrent = "+itos(iCurrent)+")");
  //#endif
  return piCurrent->n;
}; // end operator[] for linked-list implementation

/*
TBool TCluster::operator==(const TCluster& c) {
  TBool isEqual = 1;  // default state is true

  return isEqual;
};  // constructor
*/


void TCluster::sortNodes() {
  // Uses an optimized insertion sort to sort a generic integer data array.
  // implementation is a linked list, so it simplifies the insertion logic
  // we pick the next node in the list and scan previous nodes.  Then 
  // swap pointers at the correct location.
  int     i, j;
  int     nI, nJ;                 // temporary integer value
  TClusterNode  *pJ, *pI;         // pointers to i'th, j'th
  TClusterNode  *pPrevI, *pPrevJ; // pointer to previous i'th and j'th numbers
  TBool   bFound = 0;

  pI = pBegin;                    // starting location
  for(i=1; i<nNodes; i++) {
    // pick up next to test insert in earlier part of list
    pPrevI = pI;                  // pointer to previous i'th node for insertion
    pI     = pI->pNext;  nI = pI->n.node;
    if(nI<pPrevI->n.node) {       // only check for move if out of place
      pJ     = pBegin;   nJ = pJ->n.node;
      j      = 0;        bFound = 0;
      while(j<i && !bFound) {
        if(nJ>=nI)  bFound = 1;
        else { 
          pPrevJ = pJ;            // store previous position for insertion
          pJ = pJ->pNext;    nJ = pJ->n.node;
          j++;             // if not incrementing j, loop is stopped by bFound
        } // end else
      } // end while j
    
      // finally swap cluster nodes if needed
      if(bFound) {      // only swap if we found a larger node
        if(j==0) {      // j==0 is a special case because of pBegin
          pPrevI->pNext = pI->pNext;  // shortcut pI
          pBegin        = pI;
          pI->pNext     = pJ;
          //pPrevJ      = NULL;  // no need to update pPrevJ
          // this decrements pI trivially which is made up in the follow up loop       
          pI            = pPrevI;     // because pointer data changes location
        } // end if j == 0
        else if(j<i) {  // swap any nodes if after j==0
          pPrevI->pNext = pI->pNext;
          pPrevJ->pNext = pI;
          pI->pNext     = pJ;
          // this decrements pI trivially which is made up in the follow up         
          pI            = pPrevI;     // because pointer data changes location
        } // end if j<0
        if(i==nNodes-1) {
          pEnd = pPrevI;  // moved last node, so reset pEnd
          //pEnd->pNext = NULL;  // redundant with above logic
        } // if end of list
      } // end if bFound
    } // end if check i first
  } // end for i
  
  return;
} // end sortNodes


bool TCluster::isValid(string s) const {
  // perform some standard validation checks on a cluster - these are not
  // comprehensive and do not cover function specific errors
  // the string is to identify the calling function
  #ifdef DEBUG_MODE
  string isv = "via isValid()";
  if(nNodes>0) {
    if(pEnd==NULL)
      errorMsg("pEnd = NULL on a non-empty list in "+s+isv+"?");
    else if(pBegin==NULL)
      errorMsg("pBegin = NULL on a non-empty list in "+s+isv+"?");
    else if(pEnd->pNext!=NULL)   
     errorMsg(s+" detected non-NULL pointer at end of list "+isv+"?");
  } // end error checks on non-empty list 
  else if(nNodes<0) errorMsg(s+" detected a negative number of nodes "+isv+"?");
  #endif
  return 1;  // if we make it to here, we are ok as far as isValid is concerned
}; // end isValid


//int TCluster::add(int s, TBool bClean, TBool bFlag, TCluster *pc) {
int TCluster::add(int s, TBool bClean, TBool bFlag) {
  // perform a pure addition by parameters including allocation of a new 
  // cluster node, but cluster properties are *not* updated
  #ifdef DEBUG_MODE
  isValid("add(int n)");  // debugging - standard error detection code
  #endif

  //cout << "TCluster::adding node " << s << endl;  // debugging

  if(nNodes>0) {
    if(pEnd==NULL)
      errorMsg("pEnd = NULL on a non-empty list in add(int n)?");
    else if(pBegin==NULL)
      errorMsg("pBegin = NULL on a non-empty list in add(int n)?");
    else if(pEnd->pNext!=NULL)   
     errorMsg("add(int n) detected non-NULL pointer at end of list?");
  } // end error checks on non-empty list 
  //#endif

  TClusterNode *pn;  pn = new TClusterNode;
  pn->n.node  = s;   pn->n.clean = bClean;  pn->n.flag = bFlag;
  //if(pc!=NULL)  errorMsg("add(int,...) does not allow a superNode to be set");
  //pn->n.setSuperNode(*pc);  //pn->n.flag = s.flag;
  pn->pNext = NULL;  // update new end pointer to NULL
  if(nNodes==0) {
    pBegin    = pn;  // if cluster is currently empty, point pBegin
    piCurrent = pn;  iCurrent = 0;
  } // end if nNodes
  else  pEnd->pNext = pn;   // set previous end to point to added node
  pEnd = pn;         // update pEnd to point to added node
  nNodes++;
  return nNodes;
};  // add with no energy updates

int TCluster::add(TNode &s) {
  // perform a pure addition including allocation of a new cluster node
  // other cluster properties are *not* updated
  TClusterNode *pn;  pn = new TClusterNode;
  pn->n     = s;          // copy node data using overloaded '='
  pn->pNext = NULL;       // update new end pointer to NULL
  if(nNodes==0) {
    pBegin    = pn;       // if cluster is currently empty, point pBegin
    piCurrent = pn;
  } // end if nNodes
  else  pEnd->pNext = pn; // set previous end to point to added node
  pEnd = pn;              // update pEnd to point to added node
  nNodes++;
  return nNodes;
};  // add with no energy updates


int TCluster::add(TCluster &c) {
  // add a whole cluster at a time by allocating memory and copying all data
  // energy is invalid after this due to larger number of edge changes possible
  TClusterNode *pn = c.pBegin;
  for(int i=0; i<c.nNodes; i++) {
    add(pn->n);      // add this node
    pn = pn->pNext;  // move to next node in list, still ok if last node
  } // end for i
  if(pn!=NULL)  errorMsg("TCluster add exit pointer is not NULL?");

  return nNodes;
};  // add with no energy updates


int TCluster::moveCurrent(TCluster &c) {
  // moves current node in iteration to cluster c without deleting or 
  // allocating anything.  Other cluster parameters such as energy are *not* 
  // updated.  Current node must be a valid cluster node.
  // The iteration loop is not invalidated and a subsequent next() will 
  // correctly increment to the next node in the list.
  // There is the error case where a previous update has set pCurrent equal 
  // to pPrevious without a subsequent next() call for a new iteration.
  // In this case, this function can produce an error and/or jump lists.
  // We attempt to catch this error case with some debug code up front.
  //#ifdef DEBUG_MODE
  // debugging - error detection
  if(nNodes==0)
   errorMsg("Attempting to move a node from an empty cluster?");
  else if(pBegin==NULL)
   errorMsg("pBegin = NULL on moveCurrent attempt?  Empty list?");
  else if(pCurrent==NULL)
   errorMsg("Attempting moveCurrent with pCurrent = NULL?  End of list?");
  else if(pEnd->pNext!=NULL)   
   errorMsg("moveCurrent detected non-NULL pointer at end of list?");
  else if(pCurrent==pPrevious)   
   errorMsg("Attempting moveCurrent with pCurrent = pPrevious? Forgot next()?");
  //#endif

  TClusterNode *pn = pCurrent;  // store current cluster node for addition below
  // remove it from old list
  if(nNodes==1) {
    pBegin    = NULL;  pEnd     = NULL; 
    nNodes    = 0;     //energy   = 0;  partialQ = 0.0;
    piCurrent = NULL;  iCurrent = -1;  // set op[] variables to invalid values
    // finally set pCurrent to pPrevious as an indicator which also maintains
    // consistency with the iteration loop
    pCurrent = NULL;  pPrevious = NULL;// pPrevious is supposed to be redundant
  } else if(pCurrent==pBegin) { // pCurrent is at the beginning of the list
    pBegin    = pBegin->pNext;
    piCurrent = pBegin;  iCurrent = 0;
    // finally set pCurrent to pPrevious as an indicator which also maintains
    // consistency with the iteration loop
    pCurrent = NULL;  pPrevious = NULL;// pPrevious is supposed to be redundant
    nNodes--;
  } else if(pCurrent==pEnd) {   // pCurrent is at the end of the list
    //cout << red << "ok, moving end pointer now" << normText;  // debugging
    pEnd = pPrevious;
    pPrevious->pNext = NULL;
    pCurrent = pPrevious;
    piCurrent = pBegin;  iCurrent = 0;  // reset op[] iterators
    nNodes--;  // decrement number of nodes in either case
  } else {
    // pCurrent is somewhere in the middle of the list, 
    pPrevious->pNext = pCurrent->pNext; 
    pCurrent         = pPrevious;
    piCurrent = pBegin;  iCurrent = 0;  // reset op[] iterators
    nNodes--;          // decrement number of nodes in either case
  } // end else

  // add it to new list
  // we are not updating the iteration pointers on the target cluster since we
  // are adding the node to the end of the list.
  pn->pNext = NULL;  // pn is will be the end node on the target cluster c
  if(c.nNodes==0) { c.pBegin = pn;       c.pEnd = pn; } // if c is empty
  else            { c.pEnd->pNext = pn;  c.pEnd = pn; } // end else
  // no need to updata op[] stuff for target list since we add to end of c list
  c.nNodes++;

  return nNodes;
};  // end moveCurrent with no energy updates


int TCluster::moveStart(TCluster &c) {
  // moves start node in iteration to cluster c without deleting or 
  // allocating anything.Other cluster parameters such as energy are *not* 
  // updated.  Current node must be a valid cluster node.
  // The iteration loop is not invalidated and a subsequent next() will 
  // correctly increment to the next node in the list.
  //#ifdef DEBUG_MODE
  // debugging - error detection
  if(nNodes==0)
    errorMsg("Attempting to move a node from an empty cluster in moveStart?");
  else if(pBegin==NULL)
    errorMsg("pBegin = NULL on moveStart attempt?  Empty list?");
  else if(pEnd->pNext!=NULL)   
    errorMsg("moveStart detected non-NULL pointer at end of list?");
  //#endif

  TClusterNode *pn = pBegin;  // store current cluster node for addition below
  // remove it from old list
  if(nNodes==1) {
    // this case is the intended application, moving the sole member of a
    // cluster elsewhere, though the next else makes it completely general
    pBegin = NULL;     pEnd = NULL;
    piCurrent = NULL;  iCurrent = -1;  // set op[] variables to invalid values
    // finally set pCurrent to pPrevious as an indicator which also maintains
    // consistency with the iteration loop
    pPrevious = NULL;  pCurrent = NULL;
    nNodes = 0;  //energy = 0;  partialQ = 0.0;  // does not update for speed
  } else if(nNodes>1) {
    // pEnd points to a valid node other than pBegin
    pBegin    = pBegin->pNext;
    piCurrent = pBegin;  iCurrent = 0;  
    // keep iteration pointers consistent - pn is original pBegin
    if(pCurrent==pn)   pCurrent  = pBegin;
    // pPrevious is a separate case, not an 'else' since pCurrent can equal
    // pPrevious in when other code happens to set them equal for consistency.
    if(pPrevious==pn)  pPrevious = NULL;  
    nNodes--;   //energy  = 0;  partialQ = 0.0;  // does not update for speed
  } // end of valid cases
  else errorMsg("Attempting to use moveStart(c) node on an empty cluster?");

  // add it to new list
  // we are not updating the iteration pointers on the target cluster since we
  // are adding the node to the end of the list.
  pn->pNext = NULL;  // pn is will be the end node on the target cluster c
  if(c.nNodes==0) { c.pBegin = pn;       c.pEnd = pn; } // if c is empty
  else            { c.pEnd->pNext = pn;  c.pEnd = pn; } // end else
  // no need to updata op[] stuff for target list since we add to end of c list
  c.nNodes++;

  return nNodes;
};  // end moveStart with no energy updates


void  TCluster::move(TCluster &c) {
  // moves this cluster to cluster 'c' without deleting or allocating anything.
  // Cluster parameters such as energy are *not* updated.
  //cout << "Moving cluster:     " << flush;  display();    cout << endl;  // debugging
  //cout << "Moving to cluster:  " << flush;  c.display();  cout << endl;  // debugging
  //debugPause();
  if(nNodes==0)  errorMsg("Attempting to move a size zero cluster?");
  if(c.nNodes>0) {
    #ifdef DEBUG_MODE
    if(c.pEnd->pNext!=NULL)
      errorMsg("Current ending cluster does not end in pNext = NULL in move()?");
    #endif
    c.pEnd->pNext = pBegin;  c.pEnd = pEnd;  c.nNodes += nNodes; // move to c
  } else {
    // c is size zero
    //warningMsg("Moving to a size zero cluster?");
    c.pBegin = pBegin;  c.pEnd = pEnd;  c.nNodes = nNodes;
  } // end else
  pBegin = NULL;     pEnd = NULL;       nNodes = 0;        // remove from this
  pCurrent = NULL;   pPrevious = NULL;  piCurrent = NULL;  iCurrent = 0;
  // ending error check
  #ifdef DEBUG_MODE
  if(c.pEnd->pNext!=NULL)
    errorMsg("New ending cluster does not end in pNext = NULL in move()?");
  #endif
  // nEdges += would be incorrect since there could be connecting edges
  // that are not accounted for in the direct move call
  // no need to updata op[] stuff for target list since we add to end of c list
  return;
};  // add with no energy updates


int  TCluster::eraseCurrent() {
  // delete current node in linked list and return memory, but does not destruct
  // the object.  pCurrent must be a valid node.
  TClusterNode *pn = pCurrent;
  if(nNodes>1) {
    // In each case we set pCurrent backwards (or to NULL) as an indicator of
    // the change.  It also maintains consistency with the iteration loop.
    if(pCurrent==pBegin) { // pCurrent is the beginning of the list
      pBegin = pCurrent->pNext;
      pCurrent = NULL;  pPrevious = NULL;          // pCurrent is now invalid
    } // end if at start
    else if(pCurrent==pEnd) { // pCurrent is the end of the list
      pEnd = pPrevious;  pPrevious->pNext = NULL; // pPrevious is now the end
      pCurrent = pPrevious;
    } // end if at end
    else { // pCurrent is somewhere in the middle
      pPrevious->pNext = pCurrent->pNext;     // shortcut around current node
      pCurrent = pPrevious;
    } // end else
    piCurrent = pBegin;  iCurrent = 0;        // update array[] 'iterator' data
    nNodes--;
    delete pn;                                // finally return the memory
  } // end if more than one node
  else if(nNodes==1) {
    pBegin = NULL;  pEnd = NULL;  pCurrent = NULL;  pPrevious = NULL;
    piCurrent = NULL;  iCurrent = 0;          // update array[] 'iterator' data
    nNodes = 0;
    delete pn;                                // finally return the memory
  } // end else if nNodes
  else  warningMsg("Attempting to erase a node from a size zero cluster?");
  // now finish with upkeep and memory
  //piCurrent = pBegin;  iCurrent = 0;        // update array[] 'iterator' data
  //delete pn;  // finally return the actual memory

  return nNodes;
};  // end eraseCurrent


void TCluster::erase() {
  TClusterNode *pn, *pc = pBegin;
  // delete linked list
  //display();
  for(int i=0; i<nNodes; i++) { pn = pc->pNext;  delete pc;  pc = pn; }
  // must invalidate the list to indicate empty status
  pBegin   = NULL;  pEnd      = NULL;
  pCurrent = NULL;  pPrevious = NULL;  piCurrent = NULL;  iCurrent  = 0;
  nNodes   = 0;     energy    = 0;
  bFlag    = 0;
  return;
}; // end erase

/*
void TCluster::destroy() {
  TClusterNode *pn, *pc = pBegin;
  // delete linked list
  for(int i=0; i<nNodes; i++) { pn = pc->pNext;  delete pc;  pc = pn; }
  // must invalidate the list to indicate empty status
  pBegin   = NULL;  pEnd      = NULL;
  pCurrent = NULL;  pPrevious = NULL;  piCurrent = NULL;  iCurrent  = 0;
  nNodes   = 0;
  return;
}; // end destroy
*/

void TCluster::display(int nodeOffset) const {
  // crude community list omitting edge connections (no graphical capability)
  int MaxEdges = nNodes*(nNodes-1)/2;

  // output row label
  cout << brownText << "Size ";
  if(nNodes<10) cout << " ";  // crude alignment for less than 100 nodes/group
  cout << nNodes << ":  " << green;

  if(nNodes==0)  return;

  // output cluster nodes - special nodes are color-coded
  TClusterNode *pn = pBegin;
  for(int i=0; i<nNodes; i++) {
    if(pn==NULL)  errorMsg("NULL pointer for next node in TCluster display?");
    //cout << "[ with pn = " << (int)pn << flush;  // debugging
    if(pn->n.flag)                     cout << red;   // marked flagged values
    cout << (pn->n.node + nodeOffset) << " " << flush;
    if(pn->n.flag)                     cout << green; // reset to standard color
    pn = pn->pNext;  // increment to next node in list
    //cout << "] with pn = " << (int)pn << flush;  // debugging
  } // end for i
  if(pn!=NULL)  errorMsg("TCluster display exit pointer is not NULL?");

  // output some cluster parameters - not generalized for weighted edges
  if(energy == -MaxEdges && MaxEdges>0)  cout << blue; // full clique
  else                                   cout << normText;
  cout << ": with Es = " << energy;
  if(energy == -MaxEdges && energy<0)    cout << normText;
  //cout << "  dQ = " << partialQ << "  edges = " << nEdges << "\n";
  cout << ",  Ls = " << nEdges << ", and  " << magenta << "ps = ";
  if(nNodes>1)  cout << ((double)nEdges/(double)MaxEdges) << normText;
  else          cout << "N/A" << normText;
  
  return;
};  // display cluster


/*
TFloat TCluster::getModChange(int i, TCMatrix &CM) const {
  // calculate modularity change *if* particle i were to be added to cluster i
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  //TFloat mSum = 0;
  //for(int j=0; j<nNodes; j++)  mSum += CM[i][nodes[j].node]; // old version
  // try optimized version - doesn't appear to make much difference
  //int    *pI = &(CM[i][0]);
  double  mSum = 0.0;
  int     jNode;
  for(int j=0; j<nNodes; j++) {
    jNode = nodes[j].node;
    //mSum += (double)CM[i][jNode] - p[i][jNode];
  } // end for j
  return mSum;
}; // end energy change with test node
TFloat TCluster::getModChangeRemoved(int i, TCMatrix &CM) const {
  // calculate modularity change *if* particle i removed from cluster i
  // i is the number identity not the index.  It also assumes that node
  // i is *in* the cluster to be valid and this is *not* checked.
  // The i==i diagonal element is supposed to be zero for i==i.
  TFloat mSum = 0.0;
  int    jNode;
  for(int j=0; j<nNodes; j++) {
    jNode = nodes[j].node;
    //mSum += p[i][jNode] - (double)CM[i][jNode];
  } // end for j
  return mSum;
}; // end energy change with test node
*/

int TCluster::calcEnergy(TCMatrix &Cp, TBool bSymmetric) {
  // calculate the energy (and everything else) of the cluster from scratch
  // Cp is the connection matrix including connected *and* unconnected nodes
  // ki is the degree vector for each of the nodes
  // bSymmetric states whether we are dealing with a symmetric Cp (default is 1)
  TClusterNode *pa = pBegin, *pb = pBegin;
  int Cij, iNode, jNode;
  //int nMaxEdges = nNodes*(nNodes-1)/2;  // maximum # possible for cluster
  // calculate L - need to do this here because it is needed in the null model 
  // comparison.  This is inefficient since it is calculated again for each 
  // cluster in a cluster list, but it keeps the data consistent with objects.
  double Ld = (double)Cp.getNEdges();  // total number of edges in the system
  //int  L = Cp.getNEdges();           // total number of edges in the system
  cEdgeSum = 0;  uEdgeSum = 0;  // zero connected and unconnected edge sums

  if(bSymmetric) { // symmetric Cp case
    // following class variables are temporarily used as a summing variables
    nEdges   = 0;      // number of connecting edges *within* the cluster
    //partialQ = 0.0;
    energy   = 0;
    // sum over all pairs of 
    for(int i=0; i<nNodes; i++) {
      iNode = pa->n.node;
      partialQ += (double)Cp.ki(iNode);// sum partialQ only for nodes in cluster

      pb = pa->pNext;  // skips lower half of Cp just like for arrays j=i+1
      for(int j=i+1; j<nNodes; j++) {
        jNode = pb->n.node;

        Cij = Cp(iNode,jNode);
        energy -= Cij;                // weighted energy sum
        // for a weighted Cp, it is connected if Cij>0 and unconnected if Cij<=0
        //if(Cij>0)  nEdges++;        // sum number of connected edges
        if(Cij>0) { 
          nEdges++;                   // sum number of connected edges
          cEdgeSum += Cij;            // sum connected edge weight
        //} else  if(Cij<0)  uEdgeSum += -Cij;// sum abs unconnected edge weight
        // the zero terms are irrelevant for the next sum
        } else  uEdgeSum -= Cij;      // sum abs unconnected edge weight

        pb = pb->pNext;               // increment to next node in list
      } // end for j

      pa = pa->pNext;                 // increment to next node in list
    } // end for i
    //energy = -cEdgeSum + uEdgeSum;  // weighted energy sum

    double W  = (double)cEdgeSum;  // total number of edges in the system
    // calculate modularity
    //partialQ  = ( (double)nEdges - partialQ*partialQ/(4.0*Ld) )/Ld;
    //partialQW = ( (double)cEdgeSum - partialQW*partialQW/(4.0*W) )/W;// is wrong

  } else { // non-symmetric case
    // following class variables are temporarily used as a summing variables
    nEdges   = 0;      // number of connecting edges *within* the cluster
    partialQ = 0.0;
    energy   = 0;
    // The following loops are complete and will handle non-symmetric matrices.
    // This is partially required by the linked list implementation since it 
    // would require and O(n) operation to find the j=i+1 node anyhow.
    // This does assume that self-energy terms are zero or at least accounted.
    for(int i=0; i<nNodes; i++) {
      iNode = pb->n.node;
      // Q calculation may not be valid for non-symmetric matrices???
      partialQ += (double)Cp.ki(iNode);  // sum degrees of nodes in cluster
      for(int j=0; j<nNodes; j++) {
        jNode   = pb->n.node;
        Cij     = Cp(iNode,jNode);
        energy -= Cij;                // weighted energy sum
        // for a weighted Cp, it is connected if Cij>0 and unconnected if Cij<=0
        //if(Cij>0)  nEdges++;          // sum number of connected edges
        if(Cij>0) { 
          nEdges++;          // sum number of connected edges
          cEdgeSum += Cij;   // sum connected edge weight
        } else  if(Cij<0)  uEdgeSum += -Cij; // sum abs unconnected edge weight
        pb = pb->pNext;               // increment to next node in list
      } // end for j
      pa = pa->pNext;                 // increment to next node in list
    } // end for i
    if((energy%2)==1)  
      warningMsg("Odd energy detected in calcEnergy. This may be OK for non-symmetric matrices");
    energy /= 2;  // correct for double counting edges (for symmetric matrices)
    //nEdges /= 2;  // should I divide by 2 here???
    // calculate modularity as given by Sales-Pardo
    //partialQ = (double)nEdges/Ld - pow( partialQ/(2.0*Ld) ,2);
  } // end non-symmetric matrix

  return energy;
}; // end energy change with test node


void TCluster::clearEnergy() {
  // clear the energy and other variables for the cluster - this is useful and
  // faster for some initializations where the connection matrix is not needed
  //int nMaxEdges = nNodes*(nNodes-1)/2;  // maximum # possible for cluster
  cEdgeSum = 0;  uEdgeSum = 0;  // zero connected and unconnected edge sums
  nEdges   = 0;  partialQ = 0.0;  energy = 0;
  return;
}; // end clear energy
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//----------------- TClusterList Member Function Declarations------------------
TClusterList::TClusterList(int MaxSize, int NTotalNodes, int NTotalEdges) {
  // a proper constructor definition requires a filled cluster list which we
  // don't have yet, so...
  NEdges       = NTotalEdges;
  NNodes       = NTotalNodes;
  NMaxClusters = MaxSize;
  nodes.resize(NMaxClusters);
  clusters = new TCluster*[NMaxClusters];
  //nodes.resize(NNodes);
  //clusters = new TCluster*[NNodes];
  //TCluster **pc = clusters;
  //for(int i=0; i<NMaxClusters; i++) { (*pc) = NULL;  pc++; }
  //for(int i=0; i<NMaxClusters; i++)  clusters[i] = NULL;
  for(int i=0; i<NMaxClusters; i++)  clusters[i] = new TCluster;
  nClusters   = 0;    // nothing there yet
  energyTotal = 0;    // not yet able to calculate the energy
  modularity  = 0.0;  // trivial value to start
};  // constructor

TClusterList::~TClusterList() { 
  //TCluster **pc = clusters;
  //for(int i=0; i<NMaxClusters; i++) { delete (*pc);  pc++; }
  #ifdef DEBUG_MODE
  int errorCode = 0;
  #endif
  for(int i=0; i<NMaxClusters; i++) {
    #ifdef DEBUG_MODE
    // this is invalid with the new empty cluster list
    //if(i>=nClusters && clusters[i]!=NULL)  errorCode -= 1;
    if(i>=nClusters && clusters[i]->getSize()>0)  errorCode -= 1;
    #endif
    delete clusters[i];
  } // end for i
  delete[] clusters;
  #ifdef DEBUG_MODE
  if(errorCode<0)
    warningMsg("~TClusterList is non-empty for i>=nClusters? (errorCode = "
               +itos(errorCode)+")");
  #endif
}; // destructor


void  TClusterList::initClean() {
  //TCluster **pc = clusters;
  //for(int i=0; i<NMaxClusters; i++) { (*pc)->initClean();  pc++; }
  for(int i=0; i<nClusters; i++)  clusters[i]->initClean();
};  // initClean


void  TClusterList::clearFlags() {
  //TCluster **pc = clusters;
  //for(int i=0; i<NMaxClusters; i++) { (*pc)->clearFlags();  pc++; }
  for(int i=0; i<nClusters; i++)  clusters[i]->clearFlags();
};  // clearFlags


void  TClusterList::resize(unsigned size, bool bSaveData) {
  // resize the cluster list and create a new set of empty clusters
  // currently it DELETES ALL DATA
  TClusterList *pcTemp;
  if(bSaveData) { // a crude data save routine
    pcTemp = new TClusterList;
    (*pcTemp) = (*this);
  } // end bSaveData
  destroy();
  NMaxClusters = size;
  nodes.resize(NMaxClusters);
  clusters = new TCluster*[NMaxClusters];
  // create the *array of empty clusters
  for(int i=0; i<NMaxClusters; i++)  clusters[i] = new TCluster;
  if(bSaveData && NMaxClusters>=pcTemp->NMaxClusters) { 
    for(int i=0; i<pcTemp->NMaxClusters; i++)  (*clusters[i]) = (*pcTemp)[i];
    nClusters   = pcTemp->nClusters;    // nothing there yet
    energyTotal = pcTemp->energyTotal;  // not yet able to calculate the energy
    modularity  = pcTemp->modularity;   // trivial value to start
    NNodes      = pcTemp->NNodes;
    NEdges      = pcTemp->NEdges;
  } else {
    NEdges = 0;  NNodes = 0;
    nClusters   = 0;    // nothing there yet
    energyTotal = 0;    // not yet able to calculate the energy
    modularity  = 0.0;  // trivial value to start
  } // end else
}; // end resize


void TClusterList::display(int nodeOffset, TBool bSort, string s, int verbosity){
  // s   - is a string that is displayed just prior to the energy output
  //       which would typically hold the name of the cluster set for display
  // verbosity controls how much information is displayed about the system
  //   0 - nothing but the clusters themselves
  //   1 - displays energy and modularity information afterwards
  //   2 - displays... not implemented yet
  // just display each cluster in sequence
  // Note that the calculations are *not* called prior to output!
  if(bSort)  sortClusters();
  // now output information
  for(int i=0; i<nClusters; i++) {
    if(clusters[i]->getSize()>0) {
     cout << blue << "Cluster " << i << ":  ";
     clusters[i]->display(nodeOffset);
     cout << "\n";
    } // end if
  } // end for i
  if(nClusters>0)
    cout << magenta << s << "energy is " << energyTotal << " with a "
         << "modularity of " << modularity << " and " << NEdges << " edges"
         << " with " << nClusters << " clusters" << normText << "\n" << flush;
  else 
    warningMsg(s+"cluster list has size 0:  energy = 0 and modularity = 0.",
               magenta);
  return;
};  // display cluster list


void  TClusterList::initNodesList() {
  // now fill node to cluster list for the sparse matrix iteration
  TCluster  *pc;
  for(int i=0; i<nClusters; i++) {
    pc = clusters[i];
    pc->begin();
    for(int j=0; j<pc->getSize(); j++) {  // should be same as number of nodes
      nodes[pc->getCurrentNode()] = i;
      pc->next();
    } // end for j
  } // end for i
  return;
};  // end initNodesList


int TClusterList::calcNNodes() {
  // now fill node to cluster list for the sparse matrix iteration
  NNodes = 0;
  for(int i=0; i<nClusters; i++)  NNodes += clusters[i]->getn();
  return NNodes;
};  // end calcNNodes


void  TClusterList::invalidateNodesList(int i) {
  // invalidate node to cluster list for cluster i
  TCluster  *pc = clusters[i];
  pc->begin();
  for(int j=0; j<pc->getSize(); j++) {  // should be same as number of nodes
    nodes[pc->getCurrentNode()] = -1;
    pc->next();
  } // end for i
  return;
};  // end initNodesList


void  TClusterList::updateNodesList(int i) {
  // update node to cluster list for cluster i
  TCluster  *pc = clusters[i];
  pc->begin();
  for(int j=0; j<pc->getSize(); j++) {  // should be same as number of nodes
    nodes[pc->getCurrentNode()] = i;
    pc->next();
  } // end for i
  return;
};  // end updateNodesList


void  TClusterList::updateNodesList(int i, int target) {
  // update node to cluster list for cluster i
  TCluster  *pc = clusters[i];
  pc->begin();
  for(int j=0; j<pc->getSize(); j++) {  // should be same as number of nodes
    nodes[pc->getCurrentNode()] = target;
    pc->next();
  } // end for i
  return;
};  // end updateNodesList with target


void  TClusterList::initEven(int N, int &q, TBool bRandomOrder, int iBegin) {
  // initializes the vector either as sequential integers from 0 to N-1
  // bRandomOrder will randomize the order of the integers
  //cout << red << "here in initEven with N = " << N << " and q = " << q 
  //     << normText << endl; // debugging
  if(q>N) {
    warningMsg("In TClusterList initEven(N,q...) q = "+itos(q)
               +" is larger than N = "+itos(N)+".  Resetting q=N.");
    q = N;
  } // end if q
  
  cout << green << "b " << endl; // debugging
  TVectorInt v(N,-1);                   // initialize a working vector
  cout << green << "b " << endl; // debugging
  v.initStep(iBegin,1);                 // initialize nodes in consecutive order
  cout << green << "b " << endl; // debugging
  if(bRandomOrder)  v.randomizeOrder(); // now randomize the node order
  cout << green << "b " << endl; // debugging
  // random is 0 below because we already randomized the vector if desired
  initEven(v,q,0);                      // call basic function
  cout << green << "b " << endl; // debugging
  return;
}; // end initEven given integer range


void  TClusterList::initEven(TVectorInt &v, int &q, bool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  // the integers can be general
  int  N = v.getSize();  // use given size of v

  TVectorInt *pv = &v;  // define pointer in case we randomize
  if(q>N) {
    warningMsg("In TClusterList initEven(TVectorInt...) q = "+itos(q)
               +" is larger than N = "+itos(N)+".  Resetting q = N.");
    q = N;
  } // end if q

  // now randomize the node order if indicated
  if(bRandomOrder) {
    TVectorInt u(N);  // we don't want to disturb the passed vector
    u = v;
    u.randomizeOrder();
    pv = &u;          // reset v pointer for below loop to vector u
  } // end if bRandomOrder
  
  if(N==NMaxClusters)  erase();        // erase current cluster list
  else                 resize(N);
  NNodes = N; 
  nClusters = q;

  bool qdivN = !( (bool)(N%q) );      // does q divide N? logic is easier
  if(qdivN) {
    int n = N/q;                  // number of nodes per cluster
    int nCount = 0;               // nCount - track total number all nodes added
    for(int i=0; i<q; i++)        // add nodes to current cluster
      for(int k=0; k<n; k++) { clusters[i]->add((*pv)[nCount]);  nCount++; }
    
    // consistency check before exit
    if(nCount!=N)  errorMsg("Did not correctly add even nodes in initEven!");
  } else {
    // q does not divide N so initialize as evenly as possible
    int NModq = N%q, nMin = N/q;  // minimum # of nodes per cluster
    // need some variables to try to evenly distribute the 'extra' nodes evenly
    // among the clusters
    // qProb  - probability of adding an 'extra' spin to a given cluster
    // nExtra - track total number of 'extra' nodes added during process, 
    //          this should be (N mod q) by the end of the initialization
    // nCount - track total number all nodes added during process
    TFloat qProb = (TFloat)(N%q)/(TFloat)(N/q);  
    int    nCount = 0, nExtra = 0, i = 0;
    while(i<q) {
      // add minimum number of nodes to current cluster
      for(int k=0; k<nMin; k++) { clusters[i]->add((*pv)[nCount]);  nCount++; }
      // for uneven clusters, do we add another single node to this cluster?
      if(randomDouble()<qProb && nExtra<NModq) {
        clusters[i]->add((*pv)[nCount]);  
        nCount++;  nExtra++; 
      } // end if extra
      
      i++;
    } // end while i
    // because of probabilistic inclusion, we now have to ensure that we have 
    // included all 'extra' nodes.  Here we place 'left-over' nodes randomly
    // among the q clusters whereas above, we were trying to maintain 'even'
    // clusters by adding at most a single extra node to each cluster.
    int wCluster;
    while(nExtra<NModq) {
      wCluster = randomInt(0,q-1);             // choose random cluster
      clusters[wCluster]->add((*pv)[nCount]);  // add nodes to random cluster
      nCount++;  nExtra++;
    } // end while i

    // consistency checks before exit
    if(nCount!=N)  errorMsg("Did not correctly add nodes in initEven!");
    else if(nExtra!=NModq)  
                   errorMsg("Did not correctly add extra nodes in initEven!");
  } // end else (q does not divide N)

  initNodesList();
  return;
}; // end initEven given integer vector


void  TClusterList::initRandomS(int q) {
  // bRandomOrder will randomize the order of the integers
  // the integers can be general
  int  N = NNodes, wCluster;
  if(q>N)  warningMsg("In TClusterList initRandomS(q) q = "+itos(q)
                      +" is larger than N = "+itos(N)+".  Resetting q = N.");
  if(NNodes==0)  warningMsg("In TClusterList initRandomS(q) N = "+itos(N)+"?");

  // start by initializing the current cluster list to one-node-per-cluster
  // randomize the node order is redundant, so we set to zero always
  initSymmetric(0);

  // more efficiently randomize the first q clusters (rather than the list)
  //cout << "b " << endl; // debugging
  for(int i=0; i<q; i++) {
    wCluster = randomInt(q,N-1);
    //cout << "Swapped clusters are " << i << " and " << wCluster << endl; // debugging
    swapClusters(i,wCluster);
  } // end for i
  //display();  // debugging
  //cout << "a inside initRandomS with N = " << N << endl; // debugging


  // each cluster starts with at least one node.  At this point there is only 
  // one node in each cluster so we count backwards and move each node to fill
  for(int i=N-1; i>=q; i--) {
    wCluster = randomInt(0,q-1);
    //cout << "Moving node " << clusters[i]->getStartNode() << " in cluster " 
    //     << i << " to cluster " << wCluster << " " << flush; // debugging
    //clusters[i]->begin();  // begin iteration
    //clusters[i]->moveCurrent(*(clusters[wCluster]));
    clusters[i]->moveStart(*(clusters[wCluster]));
    //cout << "done with move." << endl; // debugging
  } // end for i
  nClusters = q;
  //cout << "b " << endl; // debugging

  initNodesList();
  //cout << "e " << endl; // debugging
  return;
}; // end initRandom cluster list


void  TClusterList::initEvenS(int q, bool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  // the integers can be general
  int  N = NNodes, wCluster;
  if(q>N)  warningMsg("In TClusterList initEvenS(q) q = "+itos(q)
                      +" is larger than N = "+itos(N)+".  Resetting q = N.");
  if(NNodes==0)  warningMsg("In TClusterList initEvenS(q) N = "+itos(N)+"?");

  // start by initializing the current cluster list to one-node-per-cluster
  // and randomize the node order if indicated
  initSymmetric(bRandomOrder);
  
  bool qdivN = ( (N%q)==0 );      // does q divide N? logic is easier
  if(qdivN) {
    int n = N/q;                  // number of nodes per cluster
    int cCluster = N-1;           // track last cluster in the list
    for(int i=0; i<q; i++)
      // add nodes to current cluster - each cluster starts with one node
      for(int j=1; j<n; j++) {
        //clusters[cCluster]->begin();  // begin iteration
        //clusters[cCluster]->moveCurrent(*(clusters[i]));
        clusters[cCluster]->moveStart(*(clusters[i]));
        cCluster--;
      } // end for j
    nClusters = q;

    // consistency check before exit
    if(cCluster!=q-1)  
      errorMsg("Did not correctly add even nodes in initEvenS(q)!");
  } else {
    errorMsg("initEvenS(q) is currently only implemented for N mod q = 0.");
/*
    // q does not divide N so initialize as evenly as possible
    int NModq = N%q, nMin = N/q;  // minimum # of nodes per cluster
    // need some variables to try to evenly distribute the 'extra' nodes evenly
    // among the clusters
    // qProb  - probability of adding an 'extra' spin to a given cluster
    // nExtra - track total number of 'extra' nodes added during process, 
    //          this should be (N mod q) by the end of the initialization
    // nCount - track total number all nodes added during process
    TFloat qProb = (TFloat)(N%q)/(TFloat)(N/q);  
    int    nCount = 0, nExtra = 0, i = 0;
    while(i<q) {
      // add minimum number of nodes to current cluster
      for(int k=0; k<nMin; k++) { clusters[i]->add((*pv)[nCount]);  nCount++; }
      // for uneven clusters, do we add another single node to this cluster?
      if(randomDouble()<qProb && nExtra<NModq) {
        clusters[i]->add((*pv)[nCount]);  
        nCount++;  nExtra++; 
      } // end if extra
      
      i++;
    } // end while i
    // because of probabilistic inclusion, we now have to ensure that we have 
    // included all 'extra' nodes.  Here we place 'left-over' nodes randomly
    // among the q clusters whereas above, we were trying to maintain 'even'
    // clusters by adding at most a single extra node to each cluster.
    int wCluster;
    while(nExtra<NModq) {
      wCluster = randomInt(0,q-1);             // choose random cluster
      clusters[wCluster]->add((*pv)[nCount]);  // add nodes to random cluster
      nCount++;  nExtra++;
    } // end while i

    // consistency checks before exit
    if(nCount!=N)  errorMsg("Did not correctly add nodes in initEven!");
    else if(nExtra!=NModq)  
                   errorMsg("Did not correctly add extra nodes in initEven!");
*/
  } // end else (q does not divide N)

  initNodesList();
  return;
}; // end initEven given integer vector

/*
void  TClusterList::initEven(int q, TBool bRandomOrder) {
  // initializes the cluster list given a single cluster
  NNodes = 0;
  for(int i=0; i<nClusters; i++)  NNodes += clusters[i]->getSize();

  if(q>NNodes) {
    warningMsg("In initEven(q...) q = "+itos(q)+" is larger than N = "
              +itos(NNodes)+" in initEven(...) cluster list.  Resetting q=N!");
    q = NNodes;
  } // end if q
  
  TCluster b;                   // declare a working cluster
  for(int i=0; i<nClusters; i++)  b.add(*(clusters[i]));
  // finally initialize with the temporary cluster
  initEven(b,q,bRandomOrder);
  return;
}; // end initEven given integer range


void  TClusterList::initEven(TCluster &b, int &q, TBool bRandomOrder) {
  // initialize to an even set of clusters given a single cluster
  // bRandomOrder will randomize the order of the integers
  int  N = b.getSize();  // use given size of b
  NNodes = b.getSize();  // now set the number of nodes
  //TVectorInt *pv = &b;  // define pointer in case we randomize
  if(q>N) {
    warningMsg("In initEven(cluster...) q = "+itos(q)+" is larger than N = "
               +itos(N)+" in initEven(...) cluster list.  Resetting q = N!");
    q = N;
  } // end if q

  // declare a working node *array* - not efficient but it makes the function 
  // a little safer, keeps us from modifying b if the order is randomized, 
  // and it allows us to avoid declaring TClusterList as a friend of TCluster
  TNode **a;  a = new TNode*[b.getSize()];
  for(int i=0; i<b.getSize(); i++) {
    a[i]    = new TNode;
    *(a[i]) = b[i];  // uses TCluster operator[] for linked-list cluster b
  } // end for i

  // now randomize the node order if indicated which is the default action
  if(bRandomOrder) {
    int    jS;
    TNode *pTemp;
    // ensure that every node is moved at leahttp://www.cnn.com/US/st once to a random location
    for(int k=0; k<b.getSize(); k++) {      // swap two nodes at a time
      jS    = randomInt(0,b.getSize()-1);  // not worried about duplicates
      pTemp = a[jS];  a[jS] = a[k];  a[k] = pTemp;  // swap node pointers
    } // end for k
  } // end if bRandomOrder

  erase();                        // erase current cluster list
  TBool qdivN = !((TBool)(N%q));    // does q divide N? logic is easier
  if(qdivN) {
    int n = N/q;                  // number of nodes per cluster
    int nCount = 0;               // nCount - track total number all nodes added
    TCluster c;                   // declare the temporary cluster
    for(int i=0; i<q; i++) {
      c.erase();                  // erase current cluster
      // add nodes to current cluster
      for(int k=0; k<n; k++) { c.add(*(a[nCount]));  nCount++; }
      add(c);                                // add full cluster to cluster list
    } // end for i
    // end if qdivN
    
    // consistency check before exit
    if(nCount!=N)  errorMsg("Did not correctly add even nodes in initEven!");
  } else {
    // q does not divide N so initialize as evenly as possible
    int NModq = N%q, nMin = N/q;             // minimum # of nodes per cluster
    TCluster c;                
    // need some variables to try to evenly distribute the 'extra' nodes evenly
    // among the clusters
    // qProb  - probability of adding an 'extra' spin to a given cluster
    // nExtra - track total number of 'extra' nodes added during process, 
    //          this should be (N mod q) by the end of the initialization
    // nCount - track total number all nodes added during process
    TFloat qProb = (TFloat)(N%q)/(TFloat)(N/q);  
    int    nCount = 0, nExtra = 0, i = 0;
    while(i<q) {
      c.erase();                             // erase current temporary cluster
      // add minimum number of nodes to current cluster
      for(int k=0; k<nMin; k++) { c.add(*(a[nCount]));  nCount++; }
      // for uneven clusters, do we add another single node to this cluster?
      if(randomDouble()<qProb && nExtra<NModq) {
        c.add(*(a[nCount]));  
        nCount++;  nExtra++; 
      } // end if extra
      // now add current cluster to the overall cluster list
      add(c);
      
      i++;
    } // end while i

    // because of probabilistic inclusion, we now have to ensure that we have 
    // included all 'extra' nodes.  Here we place 'left-over' nodes randomly
    // among the q clusters whereas above, we were trying to maintain 'even'
    // clusters by adding at most a single extra node to each cluster.
    int wCluster;
    while(nExtra<NModq) {
      wCluster = randomInt(0,q-1);            // choose random cluster
      clusters[wCluster]->add(*(a[nCount]));   // add nodes to random cluster
      nCount++;  nExtra++;
    } // end while i
    
    // consistency checks before exit
    if(nCount!=N)  errorMsg("Did not correctly add nodes in initEven!");
    else if(nExtra!=NModq)  
           errorMsg("Did not correctly add extra nodes in initEven!");
  } // end else (q does not divide N)

  initNodesList();
  
  for(int i=0; i<b.getSize(); i++)  delete a[i];
  delete[] a; // return local cluster array
  return;
}; // end initEven given integer range
*/

void  TClusterList::createSymmetric(TVectorInt &v, TBool bRandomOrder) {
  int N = v.getSize();
  // bRandomOrder will randomize the order of the integers
  #ifdef DEBUG_MODE_HIGH
  cout << green << "Entering createSymmetric with N = " << NNodes 
       << ", nClusters = " << nClusters << ", NMax = " << NMaxClusters 
       << ", and a target size of " << N << normText << endl;  // debugging
  #endif

  // now randomize the node order if indicated
  if(bRandomOrder)  v.randomizeOrder();
  
  // should I traverse the existing list and just move nodes rather than erasing 
  // and recreating the list?
  if(NMaxClusters==N && N>0)  erase();
  else                        resize(N);
  NNodes    = N;  // manually set the *total* number of nodes
  nClusters = N;  // and clusters

  int *pv = &(v[0]); // optimize, this ties it to the array structure of Array1D
  nodes.resize(NNodes);
  for(int i=0; i<NNodes; i++) {
    clusters[i] = new TCluster;
    clusters[i]->add(*pv);
    nodes[*pv] = i;
    pv++;  // step to next item in the node order array
  } // end for i

  clearEnergy();  // initializes the new system to zero energies
  // consistency check
  if(NNodes!=nClusters)  
    warningMsg("createSymmetric(v) did not end with correct number of clusters. nClusters = "
               +itos(nClusters)+" and NNodes = "+itos(NNodes));
  #ifdef DEBUG_MODE_HIGH
  cout << green << "Exiting createSymmetric(v,q) with N = " << NNodes 
       << ", nClusters = " << nClusters << ", and NMax = " << NMaxClusters 
       << normText << endl;  // debugging
  #endif
  return;
}; // end createSymmetric given integer vector


void  TClusterList::initSymmetric(bool bRandomOrder) {
  // Traverse the existing list and just move nodes rather than erasing and
  // recreating the entire cluster list.  This version is faster.
  // bRandomOrder will randomize the order of the clusters after all moves 
  // have been completed
  // NNodes is supposed to be unchanged through the course of the function
  // need to add some maximum cluster checks...

  //cout << green << "Entering initSymmetric with N = " << NNodes 
  //     << " and nClusters = " << nClusters << normText << endl;    // debugging

  // stop at the current max cluster though nClusters itself changes below
  int stopCluster = nClusters, nNodeChecks;

  //cout << "Here inside init A " << flush;                          // debugging
  //display(0,0,"Inside init ");
  //cout << "pq = "   << (int)clusters[nClusters]   << " and "
  //     << "pq+1 = " << (int)clusters[nClusters+1] << " and "  
  //     << "pq+2 = " << (int)clusters[nClusters+2] << endl;       // debugging

  if(nClusters<NNodes) {  // stops a trivial re-initialization
    for(int i=0; i<stopCluster; i++) {
      clusters[i]->begin();                   // begin manual pointer iteration
      nNodeChecks = clusters[i]->getSize()-1; // need to store since it changes

      // now iterate up to the last node in the cluster, but leave it there
      for(int j=0; j<nNodeChecks; j++) {
        //cout << "pq = "   << (int)clusters[nClusters]   << " and "
        //     << "pq+1 = " << (int)clusters[nClusters+1] << " and "  
        //     << "pq+2 = " << (int)clusters[nClusters+2] << endl;   // debugging
        //cout << " q" << flush;       // debugging
        addEmpty(); // increment nClusters one at a time to track last position
        //cout << " = " << flush;       // debugging
        clusters[i]->moveCurrent(*(clusters[nClusters-1]));
        //clusters[i]->display();  cout << " ";
        //clusters[nClusters-1]->display();
        //cout << nClusters << endl;       // debugging

        clusters[i]->next(); // manually iterate to the next node in cluster[i]
      } // end for j
    } // end for i
  } // end if
  //cout << "B " << flush; // debugging

  // now additionally randomize the cluster order if indicated by moving each
  // cluster to a new randomized location
  if(bRandomOrder) {
    //cout << "Ok randomizing the order in initSymmetric now" << endl;
    //display(0,0,"Before ");
    for(int i=0; i<nClusters; i++)  swapClusters(i,randomInt(0,nClusters-1));
    //display(0,0,"After ");
    //errorMsg("Stop now!");
  } // end if 
  
  //cout << "C " << flush; // debugging
  clearEnergy();  // (almost) trivially update new energy state

  //cout << "D " << endl; // debugging
  // now initialize the nodes to cluster list and set clean flag
  for(int i=0; i<nClusters; i++)  nodes[clusters[i]->getStartNode()] = i;

  //cout << "Here inside init E  " << endl; // debugging
  // consistency check
  if(NNodes!=nClusters)
    warningMsg("initSymmetric() ended with incorrect q. q = "
              + itos(nClusters) + " and NNodes = " + itos(NNodes));
  //cout << green << "Exiting initSymmetric() with N = " << NNodes
  //     << " and nClusters = " << nClusters << normText << endl;  // debugging
  return;
}; // end initSymmetric given integer vector


void  TClusterList::createPower(int N, int nMin, int nMax, double alpha, 
                                TBool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  cout << green << "Entering createPower() with N = " << NNodes 
       << " and nClusters = " << nClusters << " and NMax = " << NMaxClusters 
       << normText << endl;  // debugging

  //createSymmetric(N);  // debugging
  cout << "Allocating storage... " << flush;  // debugging
  // set up data storage
  if(NMaxClusters==N && N>0)  erase();
  else                        resize(N);

  cout << "Allocating node vector... " << flush;  // debugging
  // randomize the node order if indicated
  TVectorInt v(N);  
  v.initStep(0,1,bRandomOrder);
  //cout << "Node vector:  " << v << endl;  // debugging

  cout << "Filling easy clusters... " << flush;  // debugging
  int    nLeft = N, ix, iq, nAdded = 0, nCount;
  nClusters = 0;
  // consider case where nLeft > nMax.  
  // We are guaranteed to be able to create another complete cluster.
  while(nLeft > nMax) {
    ix = (int)( randomPower(nMin,nMax,alpha) + 0.5 );  // get next cluster size
    nCount = 0;
    while(nCount < ix) {
      //cout << "Creating cluster " << nClusters << " with n = " 
      //     << clusters[nClusters]->getn() << "... " << flush;  // debugging
      clusters[nClusters]->add(v[nAdded]);
      nCount++;
      nAdded++;
    } // end while nAdded
    nClusters++;
    nLeft -= nCount;
  } // end while nLeft

  cout << "Filling questionable clusters... " << flush;  // debugging
  // consider case where nMin <= nLeft <= nMax 
  // We are not sure if we can create a complete new cluster.
  while(nLeft >= nMin) {
    cout << "Generating random number " << flush;  // debugging
    ix = (int)( randomPower(nMin,nMax,alpha) + 0.5 );  // get next cluster size
    cout << "... " << flush;  // debugging
    // since iLeft >= nMin, but ix > nLeft we'll just create a new cluster 
    // with the remaining nodes if ix > nLeft
    if(ix > nLeft)  ix = nLeft;
    // now add the new cluster
    nCount = 0;
    while(nCount < ix) {
      //cout << "Creating cluster " << nClusters << " with n = " 
      //     << clusters[nClusters]->getn() << "... " << flush;  // debugging
      clusters[nClusters]->add(v[nAdded]);
      nCount++;
      nAdded++;
    } // end while nAdded
    nClusters++;
    nLeft  -= nCount;
  } // end while nLeft
  
  cout << "Dispersing remaining " << red << nLeft << normText << " nodes... " 
       << flush;  // debugging
  // since the number of remaining nodes is less than nMin, spread them
  // randomly among the existing communities subject to maximum size limit
  while(nLeft > 0) {
    iq = randomInt(0,nClusters-1);
    if(clusters[iq]->getn() < nMax) {
      clusters[iq]->add(v[nAdded]);
      nAdded++;
      nLeft--;
    } // end if
  } // end while nLeft
  NNodes = nAdded;
  initNodesList();

  // consistency check
  if(NNodes!=N)  
    warningMsg("createPower() did not end with correct number of nodes. N = "
               +itos(N)+" and NNodes = "+itos(NNodes));

  cout << green << "Exiting createPower() with N = " << NNodes 
       << " and nClusters = " << nClusters << normText << endl;    // debugging
  return;
}; // end createPower given integer vector


void  TClusterList::initPowerS(int N, int nMin, int nMax, double alpha, 
                               TBool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  //cout << green << "Entering initPowerS() with N = " << NNodes 
  //     << " and nClusters = " << nClusters << " and NMax = " << NMaxClusters 
  //     << " and alpha = " << alpha << normText << endl;  // debugging

  initSymmetric(1);

  //cout << "Filling easy clusters starting with nClusters = " << nClusters 
  //     << "... " << flush;  // debugging
  int    nLeft = N, ix, iq, nAdded = 0, nCount, clusterCount = 0;
  // consider case where nLeft > nMax.  
  // We are guaranteed to be able to create another complete cluster.
  while(nLeft > nMax) {
    ix = (int)( randomPower(nMin,nMax,alpha) + 0.5 );  // get next cluster size
    nCount = 0;
    // with a symmetric init, each cluster already has one node so use (ix-1)
    while(nCount < ix-1) {
      //cout << "Adding node " << clusters[nClusters-1]->getStartNode() 
      //     << " to cluster " << clusterCount << " with n = " 
      //     << clusters[clusterCount]->getn() << "... " << flush;  // debugging
      move(nClusters-1,clusterCount);
      nCount++;
      nAdded++;
    } // end while nAdded
    clusterCount++;
    nLeft -= nCount+1;  // must add extra node already in cluster from init
    nAdded++;           // must add extra node already in cluster from init
  } // end while nLeft

  //cout << "Filling questionable clusters starting with nAdded = " << nAdded 
  //     << " and nLeft = " << nLeft << "... " << flush;  // debugging
  // consider case where nMin <= nLeft <= nMax 
  // We are not sure if we can create a complete new cluster.
  while(nLeft >= nMin) {
    //cout << "Generating random number " << flush;  // debugging
    ix = (int)( randomPower(nMin,nMax,alpha) + 0.5 );  // get next cluster size
    //cout << "... " << flush;  // debugging
    // since iLeft >= nMin, but ix > nLeft we'll just create a new cluster 
    // with the remaining nodes if ix > nLeft
    if(ix > nLeft)  ix = nLeft;
    // now add the new cluster
    nCount = 0;
    // with a symmetric init, each cluster already has one node so use (ix-1)
    while(nCount < ix-1) {
      //cout << "Creating cluster " << nClusters << " with n = " 
      //     << clusters[nClusters]->getn() << "... " << flush;  // debugging
      move(nClusters-1,clusterCount);
      nCount++;
      nAdded++;
    } // end while nAdded
    clusterCount++;
    nLeft -= nCount+1;  // must add extra node already in cluster from init
    nAdded++;           // must add extra node already in cluster from init
  } // end while nLeft
  clearZeroClusters();  // debugging
  
  //cout << "Dispersing starting with nAdded = " << nAdded  << " with " 
  //     << red << nLeft << normText << " remaining nodes and " 
  //     << nClusters << " clusters" << endl;  // debugging
  // since the number of remaining nodes is less than nMin, spread them
  // randomly among the existing communities subject to maximum size limit
  while(nLeft > 0) {
    iq = randomInt(0,nClusters-1-nLeft);
    if(clusters[iq]->getn() < nMax) {
      move(nClusters-1,iq);
      nAdded++;
      nLeft--;
    } // end if
  } // end while nLeft
  NNodes = nAdded;
  initNodesList();
  //display(0,1);  // debugging

  // consistency check
  if(NNodes!=N)  
    errorMsg("initPowerS() did not end with correct number of nodes. N = "
             +itos(N)+" and NNodes = "+itos(NNodes));

  //display();  // debugging
  //cout << green << "Exiting initPowerS() with N = " << NNodes 
  //     << " and nClusters = " << nClusters << normText << endl;    // debugging
  return;
}; // end initPowerS


void  TClusterList::randomizeOrder() {
  // Traverse the existing list and simply randomizes the order of the clusters.
  // This is essentially an ad-hoc function for the 'NZMod' algorithm.
  for(int i=0; i<nClusters; i++)  swapClusters(i,randomInt(0,nClusters-1));
  //cout << green << "Exiting randomizeOrder() with N = " << NNodes 
  //     << " and nClusters = " << nClusters << normText  << endl; // debugging
  return;
}; // end randomize order


void  TClusterList::initRandom(int N, int &q, int iBegin) {
  // bRandomOrder will randomize the order of the integers
  cout << "A " << endl; // debugging
  TVectorInt v(N,-1);              // initialize a working vector
  cout << "B " << endl; // debugging
  v.initStep(iBegin,1);            // initialize nodes in consecutive order
  cout << "C " << endl; // debugging
  // random is 0 below because we already randomized the vector if desired
  initRandom(v,q);                 // call basic function
  cout << "D " << endl; // debugging
  return;
}; // end initRandom cluster list


void  TClusterList::initRandom(TVectorInt &v, int &q) {
  errorMsg("There is a bug in initRandom... temporarily use initRandomS");

  // initialize by randomly choosing clusters to add nodes
  int  N = v.getSize();              // use given size of v
  int  wCluster;                     // add node to which cluster?
  if(q>N) {
    warningMsg("In initRandom q = "+itos(q)+" is larger than N = "
               +itos(N)+" in initEven(...) cluster list.  Resetting q = N!");
    q = N;
  } // end if q

  // randomize the node order - partially redundant with q..N part loop below
  TVectorInt u(N);       // we don't disturb the passed vector
  u = v; 
  u.randomizeOrder();

  // completely get rid of current cluster list either by erasing or resizing
  cout << "a " << endl; // debugging
  if(N==NMaxClusters) {
    cout << "erase " << flush;
    erase();   
  } else {
    cout << "resize " << flush;
    resize(N);
  }
  NNodes = N; 
  nClusters = q;
  cout << "b " << endl; // debugging

  // since q is specified by constraint, give each cluster at least one node
  for(int i=0; i<q; i++)  clusters[i]->add(u[i]);
  cout << "c " << endl; // debugging

  //display(); // debugging
  // now add rest of nodes randomly
  for(int i=q; i<N; i++) {
    wCluster = randomInt(0,q-1);    // choose random cluster
    clusters[wCluster]->add(u[i]);  // add nodes to random cluster
  } // end for i - random fill

  cout << "d " << endl; // debugging
  initNodesList();
  cout << "e " << endl; // debugging

  return;
}; // end initRandom cluster list


int  TClusterList::createEmpty() {
  // add an empty cluster to the end of the list
  #ifdef DEBUG_MODE
  if(nClusters==NMaxClusters)
    errorMsg("Cannot add another cluster to the list!  nClusters = "+
              itos(nClusters));
  #endif
  clusters[nClusters] = new TCluster;  // new version does not re-allocate
  nClusters++;
  return nClusters;
}; // add a single cluster to the end of the list


int  TClusterList::add(TCluster &c) {
  // add a cluster to this list
  //cout << red << "here in add(cluster) with cluster size = " << c.getSize() 
  //     << " and nClusters of " << nClusters << endl; // debugging
  if(nClusters==NMaxClusters)
    errorMsg("Cannot add another cluster to the list!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero cluster to the list.");
  else {
    //clusters[nClusters] = new TCluster;// new method always has empty clusters
    *(clusters[nClusters]) = c;
    updateNodesList(nClusters);  // inefficient - updates entire cluster
    nClusters++;
  } // end else

  return nClusters;
}; // add a single cluster to the end of the list


int  TClusterList::add(TClusterList &c) {
  // add a whole cluster at a time by copying data
  if(nClusters>NMaxClusters-c.getSize())
    errorMsg("Appended cluster list is too large to add to this list!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero cluster list to the list.");
  else { // add the list clusters
    for(int i=0; i<c.getSize(); i++)  add(*(c.clusters[i]));
  } // end else
  return nClusters;
}; // add a single cluster to the end of the list

/*
int  TClusterList::move(TCluster &c) {
  // add a cluster to this list by pointer assignment
  // note that this function has no mechanism to unset other pointers to 'c'
  // In particular, local objects could be destructed and this pointer would
  // point to unallocated memory.
  if(nClusters==MaxClusters)  
    errorMsg("Cannot add another cluster to the list!");
  else if(c.getSize()<1)  
    warningMsg("Attempting to add a size zero cluster to the list.");
  else { 
    clusters[nClusters] = &c;
    nClusters++;
  } // end else
  return nClusters;
}; // add cluster
*/

int TClusterList::move(TClusterList &c) {
  // move a whole cluster at a time by pointer assignment
  if(nClusters>NMaxClusters-c.getSize())
    errorMsg("Appended cluster list is too large to add to this list!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero cluster list to the list.");
  else { 
    for(int i=0; i<c.getSize(); i++) {
      c.destroy(i,1);  // need to return memory even if it is empty
      c.clusters[nClusters+i] = clusters[i];  
      updateNodesList(nClusters+i);
    } // end for i
    // finally update parameters of moved list
    c.nClusters += nClusters;
    nClusters = 0;    NNodes = 0;       NEdges = 0;
    energyTotal = 0;  modularity = 0.0;
  } // end else

  return c.nClusters;
}; // move cluster list


void TClusterList::offsetNodes(int offset) {
  // add 'i' to each node index ID
  // Note that this directly changes the indices of the cluster and will 
  // implicitly invalidate the connection matrix without additional work
  int nCount = 0;  // debugging
  //cout << red << "offset = " << offset << "\n" << normText << flush;  // debugging
  for(int i=0; i<nClusters; i++) {
    //cout << "Start " << nCount << ".\nNode i = " << endl;
    clusters[i]->begin();   // begin manual iteration
    for(int j=0; j<clusters[i]->getn(); j++) {
      // add i to each node
      //cout << nCount << ": " << clusters[i]->getCurrent().node << "+" << offset 
      //     << "=" << flush;  // debugging
      clusters[i]->getCurrent().node += offset;
      //cout << clusters[i]->getCurrent().node << "?  " << flush; // debugging
      clusters[i]->next();  // step to next node in cluster
      nCount++; // debugging
    } // end for j
    //cout << " end" << nCount << "." << endl;  // debugging
  } // end for i
  //cout << "exiting offset." << endl;  // debugging
  return;
};

int TClusterList::clearZeroClusters(TBool bUpdateNodesList) {
  // Eliminate size zero clusters, but do not return the memory.
  // A small optimization of the base NZ community detection algorithm is to 
  // leave size zero clusters in the list after a few iterations which saves 
  // many node-cluster list updates.
  // Start at end of list to avoid having to check for non-empty clusters.
  // Order is not preserved. 
  int nErased = 0, startCluster = nClusters-1;
  
  for(int k=startCluster; k>=0; k--)
    if(clusters[k]->getSize() == 0) { 
      // swap i with last cluster, last cluster is not empty now
      if(k<nClusters-1)  swapClusters(k,nClusters-1); 
      nClusters--;  // simply decrement cluster counter, do not return memory
      nErased++;
    } // end if
  // now update the nodes list from scratch, there could be lots of changes
  if(bUpdateNodesList)  initNodesList();
  return nErased;
}; // end clearZeroClusters


void  TClusterList::erase() {
  // delete all clusters but don't return cluster pointer array (just clusters)
  // leave empty clusters intact
  for(int i=0; i<nClusters; i++) { clusters[i]->erase(); }
  nClusters   = 0;  NNodes      = 0;   NEdges = 0;
  energyTotal = 0;  modularity = 0.0;
  nodes.init(-1);
  return;
}; // delete all clusters


int  TClusterList::erase(int i, TBool bFast) {
  // for bFast speed, the order of clusters is *not* preserved
  if(bFast) {  // more efficient
    if(i<nClusters-1) {
      swapClusters(i,nClusters-1);  // swap i with last
      updateNodesList(i);           // only update if swapping clusters
    } // end if

    if(clusters[nClusters-1]->getSize()>0) {
      // invalidate nodes in nodes-cluster list before erasing
      invalidateNodesList(nClusters-1);
      clusters[nClusters-1]->erase();  // erase cluster but do not return memory
    } // end if size check

  } else {
    if(clusters[i]->getSize()>0)  invalidateNodesList(i);
    TCluster *pc = clusters[i];
    pc->erase();  // erase cluster but do not return memory
    // now copy the pointers down to preserve the order - can easily optimize
    for(int j=i; j<nClusters-1; j++) {
      clusters[j] = clusters[j+1];
      updateNodesList(j);               // inefficient
    } // end for j
    clusters[nClusters-1] = pc;         // now 'move' i to last cluster
  } // end else
  // do not invalidate the empty cluster in new method
  nClusters--;

  //cout << "done." << endl;
  return nClusters;
}; // erase cluster i


void  TClusterList::destroy() {
  // delete all clusters but don't return cluster pointer array (just clusters)
  for(int i=0; i<nClusters; i++) { delete clusters[i];  clusters[i] = NULL; }
  delete[] clusters;
  nClusters   = 0;  NNodes      = 0;  NEdges = 0;
  energyTotal = 0;  modularity = 0.0;
  nodes.resize(0);
  return;
}; // destroy all clusters


int  TClusterList::destroy(int i, TBool bFast) {
  // for bFast speed, the order of clusters is *not* preserved
  if(bFast) {  // more efficient
    if(i<nClusters-1) {
      swapClusters(i,nClusters-1);  // swap i with last
      updateNodesList(i);           // only update if swapping clusters
    } // end if
    // invalidate nodes in nodes list before erasing
    if(clusters[nClusters-1]->getSize()>0)  invalidateNodesList(nClusters-1);
    clusters[nClusters-1]->erase();
  } else {
    if(clusters[i]->getSize()>0)  invalidateNodesList(i);
    clusters[i]->erase();
    // now copy the pointers down to preserve the order - can easily optimize
    for(int j=i; j<nClusters-1; j++) {
      clusters[j] = clusters[j+1];
      updateNodesList(j);                        // inefficient
    } // end for j
  } // end else
  // in both cases we return the memory for the empty cluster at the end
  delete clusters[nClusters-1];  clusters[nClusters-1] = NULL;
  nClusters--;

  return nClusters;
}; // destroy cluster i


TClusterList& TClusterList::operator=(const TClusterList& b) {
  // copy cluster b to the current cluster, resizing if necessary
  //display(0,0,"best inside operator=");
  //b.display();
  //cout << "erasing... " << flush; // debugging
  //cout << "a with " << nClusters << " and " << NMaxClusters  
  //     << " and " << b.NMaxClusters << "... "    << endl; // debugging
  if(NMaxClusters==b.NMaxClusters)  erase();
  else                              resize(b.NMaxClusters);

  //cout << "copying clusters... " << flush;  // debugging
  //cout << "b... " << endl; // debugging
  // copy cluster by cluster to the empty list
  for(int i=0; i<b.nClusters; i++) {
    //cout << i << " " << *(b.clusters[i]) << endl; // debugging
    *(clusters[i]) = *(b.clusters[i]);  
  } // end for i
  //cout << "c... " << endl; // debugging

  //cout << "copying extras... " << flush;    // debugging
  // now finish with specific TClusterList info
  NMaxClusters = b.NMaxClusters;  nClusters   = b.nClusters;
  NNodes       = b.NNodes;        NEdges      = b.NEdges;
  energyTotal  = b.energyTotal;   modularity  = b.modularity;
  nodes        = b.nodes;  // copy nodes list with op= for TVectorInt
  //cout << "d... " << endl; // debugging

  //cout << "done with operator=. " << flush; // debugging
  return *this;
}; // end operator=


void TClusterList::sortClusters() {
  // sorts clusters by the first node in the list.  This will, of course,   
  // be a problem when the first node is the missing or extra node.
  // Uses an insertion sort
  int        i, j;
  TCluster  *pTemp;  // pointer to temp cluster

  //cout << "Here a in sort clusters... " << endl;  // debugging
  for(i=0; i<nClusters; i++) { 
    //clusters[i]->display();  // debugging
    //cout << "Here a1 in sort clusters with i = " << i << " and cluster size of "
    //     << clusters[i]->getSize();  // debugging
    if(clusters[i]->getSize()>1)  clusters[i]->sortNodes();
    //cout << " after... " << endl;  // debugging
    //clusters[i]->display();  // debugging
  } // end for i

  //cout << "Here b in sort clusters... " << endl;  // debugging

  for(i=1; i<nClusters; i++) {
    j = i;                 // j and i counters start at same number index
    pTemp = clusters[i];   // store next item in list to drop in correct place
    while(j>0 && 
          //clusters[j-1]->getStart().getNode()>pTemp->getStart().getNode() ) {
          (clusters[j-1]->getSize()>pTemp->getSize() // sort by size
          // and sub-sort by particle number if size is the same
           || (clusters[j-1]->getSize()==pTemp->getSize() && 
               clusters[j-1]->getStart().getNode()>pTemp->getStart().getNode()))
         ) {
    // while( j>0 && clusters[j-1]->getSize()>pTemp->getSize() ) {
      clusters[j] = clusters[j-1]; // move cluster pointer up to next spot 
      j--;
    } // end while j
    clusters[j] = pTemp;   // move cluster pointer up to next spot
  } // end for i
 
  initNodesList();  // update node-cluster list - innefficient
  //cout << "Here at end of sort clusters... " << endl;  // debugging

  return;
} // end sortClusters


int  TClusterList::calcEnergy(TCMatrix &Cp, TBool bFastCalc) {
  // requires energies of clusters to already be updated
  TCluster *pc;
  // class variables temporarily used a summing variables
  NNodes      = 0;  NEdges     = 0;
  energyTotal = 0;  
  //modularity = 0.0;  modularityW = 0.0;
  int cWeight = 0,  uWeight    = 0;
  cEdgeSum    = 0;  uEdgeSum   = 0;

  if(bFastCalc) {
    // Now do the calculation.  This fast NZ version requires a constant 'b' 
    // (unconnected edge weight).
    errorMsg("Fast O(NZ) energy calculation is not yet implemented.");
  } else {
    // Do the calculation for the 'slow' version.  Currently it is an O(Nn) 
    // algorithm, but it can be used with an arbitrary Cp.
    //cout << magenta << "nClusters = " << nClusters << endl;     // debugging
    for(int i=0; i<nClusters; i++) {
      pc = clusters[i];      // pointer to current cluster
      //n  = pc->nNodes;     // number of nodes in current cluster
      //cout << magenta << i << flush;          // debugging
      pc->calcEnergy(Cp); // debugging
      //cout << "." << flush;                   // debugging
      energyTotal += pc->getEnergy();
      //cout << "." << flush;                   // debugging
      //modularity  += pc->getPartialQ();
      //cout << "." << flush;                   // debugging
      //cEdgeSum  += pc->getCEdgeSum();
      //uEdgeSum  += pc->getUEdgeSum();
      // following are redundant after first call but a complete constructor 
      // defintion requires a filled cluster list so...
      NNodes += pc->getSize();
      NEdges += pc->getNEdges();
    } // end for interior Edges
/*
    // now calculate the normalized energy
    // find weights for normalized energy calculation.  Anecdotal data suggests 
    // that normalized energy shows resolution limit effects.
    // currently assumes a symmetric case!
    // Currently is an N^2 algorithm, but it can be used with an arbitrary Cp.
    int cij;
    for(unsigned i=0; i<Cp.getCols(); i++) {
      for(unsigned j=i+1; j<Cp.getRows(); j++) {
        cij = Cp(i,j);
        if(cij>0)       cWeight += cij;
        else if(cij<0)  uWeight += -cij;
    }} // end ij
    //eNorm = ((double)cEdgeSum/(double)cWeight - (double)uEdgeSum/(double)uWeight);
    eNorm = (double)energyTotal/(double)cWeight;
*/
  } // end else

  return energyTotal;
}; // end calcEnergy


double TClusterList::calcQ(TCMatrix &Cp, double gamma, double r) { 
  //, TBool bFast) {
  // initialize data for sparse cases
  //TVectorInt  edgeCount(N,0);  // count number of nodes connected to a cluster
  if(NNodes==0)  errorMsg("calcQ() has zero nodes?");
  if(!floatEq(r,0.0))  errorMsg("AFG r<>0 is not yet implemented.");
  
  #ifdef MSL_TCMWEIGHTED
  //errorMsg("TClusterList::calcQ() is not yet implemented for weighted graphs.");
  #endif

  if(nClusters==1)  return 0.0;  // trivial case
  int N = Cp.getSize();
  //cout << brown << "N = " << NNodes << " in calcQ()"; // debugging

  //TVectorFloat  wss(N,0.0);      // connected weight sum for each cluster
  TVectorFloat  ws(N,0.0);       // total weight sum for each cluster
  //TVectorFloat  wsk(N,0.0);      // connected weight sum for node k to cluster s
  TVectorFloat  Wk(N,0.0);       // weight sum for each node - constant
  double        W2 = 0.0;        // 2*constant total weight of the system edges
  double        g2WNr;           // scale factor optimization variable
  double        dQ, maxdQ;       // modularity changes
  int    *pNeighborStart, *pNeighbor; // optimzing neighbor pointers
  
  // -------------------------------------------------------------------------
  // calculate starting modularity
  // we manually calculate the modularity since it is secondary to the cluster
  // class implementations and is not tracked consistently in the class
  // Most of the time, we will call this function with a 'symmetric' 
  // initialized system, but this initialization code is completely general.
  modularity = 0.0;  // clear modularity to use as a summing variable
  // start by initializing the weight (or degree) sum matrix for each cluster
  double   dCjk, wSum, wssiCluster;
  int      jNode, jDegrees, kNode, kCluster;
  for(int iCluster=0; iCluster<nClusters; iCluster++) {
    clusters[iCluster]->begin();  // begin manual iteration
    wssiCluster = 0.0;  // clear current cluster connected weight sum

    for(int j=0; j<clusters[iCluster]->getn(); j++) {
      jNode    = clusters[iCluster]->getCurrentNode();
      jDegrees = Cp.ki(jNode);
      wSum = 0.0;       // clear current node weight sum

      // now we sum the actual internal weights
      // scan edge connection list for iNode and sum edge weight contributions
      pNeighborStart = &(Cp.kij[jNode][0]);
      pNeighbor = pNeighborStart;
      for(int k=0; k<jDegrees; k++) {
        kNode = (*pNeighbor);
        kCluster = nodes[kNode];
        #ifdef MSL_TCMWEIGHTED
        dCjk = (double)Cp(jNode,kNode);
        #else
        dCjk = 1.0;   // unweighted version
        #endif
        // increment internal cluster weight sums - This sum tracks the actual 
        // connected weights of jNode for iCluster.  For this initialization,
        // we do not need to track the other clusters.
        wSum += dCjk; // increment weight sum for this node
        // increment internally connected edge weights for this node
        if(kCluster==iCluster)  wssiCluster += dCjk; // double counts

        pNeighbor++;       // increment neighbor pointer
      } // end for k - neighbor list search

      // increment weight sums
      ws[iCluster] += wSum;   // sum total weight for cluster i
      //Wk[jNode]  = wSum;    // sum total weight for jNode
      // sum total weight of system (doubles counts edges on purpose)
      W2           += wSum;   // sum total weight of the system - double counts

      clusters[iCluster]->next();
    } // end for j

    // sum starting modularity and scale by system weight
    // this update covers both the RB (gamma) and the AFG (r) methods using
    // *pQ as a temporary summing variable for the internal edge weights
    //wss[i] = wssiCluster/2.0;  // not actually needed below - why define it?
    //*pQ += 2.0*wss[i] + clusters[i]->getn()*r;
    //*pQ += wssiCluster + clusters[iCluster]->getn()*r;
    // without r
    modularity += wssiCluster;  // edge contribution of iCluster - double counts
  } // end for iCluster

  // correct for double counted edge weights
  //W /= 2.0;  // we omit the correction since W is always used as 2*W
  //set overall scale factor (used in whole most of the time)
  //g2WNr = gamma/( W2 + (double)N*r );
  // without r
  g2WNr = gamma/W2;
  // Finish starting modularity absent the overall system scale factor.
  // We do not perform this overall scale until the end of the function
  // That is, we use a non-normalized modularity for the actual iteration and
  // scale the overall result at the end.  This saves a lot of operations.
  // This update covers both the RB (gamma) and the AFG (r) methods.
  // We used *pQ as a temporary sum variable for all internal edge weights.
  for(int jCluster=0; jCluster<nClusters; jCluster++)
    modularity -= g2WNr*ws[jCluster]*ws[jCluster];
    //modularity -= g2WNr*pow( ws[iCluster]+clusters[iCluster]->getnD()*r ,2);
  // end get starting modularity

  // before exiting, finally scale the working modularity to be consistent with 
  // the normalized definition.  This scaling incorporates both the AFG scaling.
  //modularity /= (W2 + (double)N*r);
  modularity /= W2;
/*
  // debugging
  double degreeSqSum = 0.0, symmQ;
  for(int i=0; i<N; i++)  degreeSqSum += pow( (double)Cp.ki(i),2 );
  symmQ = -degreeSqSum/ pow( 2.0*(double)Cp.getL(),2 );
  cout << "The symmetric Q is " << symmQ << endl;
  debugPause("Waiting...");
  // debugging
*/
  return modularity;
}; // end calcQ


void TClusterList::clearEnergy() {
  NEdges      = 0;  cEdgeSum    = 0;  uEdgeSum   = 0;
  energyTotal = 0;  modularity = 0.0;
  // class variables temporarily used a summing variables
  //NNodes      = 0;  NNodes does not change

  // Do the calculation for the 'slow' version.  Currently it is an O(Nn) 
  // algorithm, but it can be used with an arbitrary Cp.
  for(int i=0; i<nClusters; i++) {
    clusters[i]->clearEnergy(); 
    //NNodes   += clusters[i]->getSize();
    //cEdgeSum += pc->getCEdgeSum();
    //uEdgeSum += pc->getUEdgeSum();
  } // end for interior Edges

  return;
}; // end clearEnergy


void TClusterList::calcPZParams(TCMatrix &Cp, 
                                double &zin, double &zout, double &zAvg, 
                                double &pin, double &pout, double &p) {
  #ifdef DEBUG_MODE
  if(NNodes<1)  warningMsg("Number of nodes is not set in calcZParams(.)?");
  #endif
  // for the moment, we assume that calcEnergy has been called to initialize
  // various variables such as nEdges, etc.
  // use as a summing variable, use double locally to avoid max int problems
  double  NEdgesOutMax = 0.0, NEdgesInMax = 0.0;
  double  nEdgesOut = 0.0,    nEdgesIn = 0.0, degreeSum = 0.0,
          L = (double)Cp.getL();
  
  // trivial case
  for(int i=0; i<nClusters; i++) {
    NEdgesInMax += 0.5*clusters[i]->getnD()*(clusters[i]->getnD()-1.0);
    // sum total number of internal edges
    nEdgesIn += (double)clusters[i]->getl();
    for(int j=i+1; j<nClusters; j++)
      NEdgesOutMax += clusters[i]->getnD()*clusters[j]->getnD();

    // now get the degree sum for cluster i to easily calculate nEdgesOut
    degreeSum = 0.0;
    // begin manual iteration via easy-to-use function.  I probably need to 
    // change this for the internal implementation at some point.
    clusters[i]->begin();  
    for(int k=0; k<clusters[i]->getn(); k++) {
      //degreeSum += clusters[i]->getCurrent().ki(k);
      degreeSum += (double)Cp.ki(clusters[i]->getCurrentNode());
      clusters[i]->next();  // go to next node
    } // end for k
    nEdgesOut += degreeSum/2.0 - (double)clusters[i]->getl();
  } // end for i

  // now calculate the de-referenced variables
  pin  = nEdgesIn/NEdgesInMax;
  p    = 2.0*L/( (double)NNodes*(double)(NNodes-1) );
  zAvg = 2.0*L/(double)NNodes;
  zin  = 2.0*nEdgesIn/(double)NNodes;
  if(nClusters>1 && NEdgesOutMax>0.9999) {
    pout = nEdgesOut/NEdgesOutMax;
    zout = 2.0*nEdgesOut/(double)NNodes;
  } else {  // problem case of one cluster or no external edges anywhere
    pout = 0.0;  zout = 0.0;
  } // end else

  #ifdef DEBUG_MODE  
  cout << brown 
       << "NEdgesInMax = " << NEdgesInMax << ", NEdgesOutMax = " << NEdgesOutMax
       << ", nEdgesIn = " << nEdgesIn << ", nEdgesOut = " << nEdgesOut 
       << red << endl
       << "zin = " << zin << ", zout = " << zout << ", Z = " << zAvg
       << ", pin = " << pin << ", pout = " << pout << ", p = " << p 
       << normText << endl;  // debugging
  #endif
  
  return;
}; // end calcPZParams
void TClusterList::calcZParams(TCMatrix &Cp, 
       double &zIn, double &zOut, double &zAvg) {
  double pinDummy, poutDummy, pDummy;
  calcPZParams(Cp,zIn,zOut,zAvg,pinDummy,poutDummy,pDummy);
  return;
}; // end calcPZParams without density reference variables


//template <typename T> 
int TClusterList::inputLFRANS(string fname, int N) {
  // Input the cluster list identities in Lancichinetti et al. communities.dat 
  // format.  Note that they use a cluster and node offset of 1 which we 
  // automatically correct internally to be 0 based in both cases.
  #ifdef DEBUG_MODE
  //if(!isValid())  errorMsg("TClusterList object is not valid in inputLFRANS()");
  #endif
  string sWord;  // temporary line variable
  char   cNext;  // temporary character variable - good for checking line type
  const int MaxLineLength = 1024;
  char   buffer[MaxLineLength]; // char* variable for getline() function
  int    node, cluster;         // particles 1 and 2 of the defined edge
  int    nodeCount = 0;         // number of nodes counted
  int    sOffset = 1;  // Note that they use node and cluster offsets of 1 which
                       // we automatically correct internally to be 0 based.
  // LFR does not put information up front on the number of clusters q, so we
  // must determine q from the partition memberships via qMax.  We then set q
  // and nClusters at the end of the function.
  int    qMax = -1; // actual and current number of clusters

  if(fname=="")  errorMsg("No filename specified in TCMSparseW inputLFRDAT()!");
  //cout << green << "Importing file: \"" << fname << "\" w/N = " << N << "... ";

  erase();
  if(NMaxClusters!=N)  resize(N);
  nClusters = N; // temporarily assume N clusters as max possible (just in case)
  NNodes = N;

  ifstream din(fname.c_str());
  if(din.bad()) {
    errorMsg("Input filename "+fname+" does not exist!");
    return -1;  // Error filename does not exist
  } // if din.bad
  //cout << "  Reading LFR cluster list file... ";

  while(!din.eof() && nodeCount<N) {
    din >> node >> cluster;             // input data
    din.getline(buffer,MaxLineLength);  // get rest of node line
    //cout << "  node = " << node << "   cluster = " << cluster << endl; // debugging

    if(cluster>qMax)  qMax = cluster;   // track the best guess on the system q

    // NOTE:  This assumes that weight has to be >0 integer
    if(node>=sOffset && cluster>=sOffset && cluster<(N+sOffset)) {
      nodeCount++;
      //cout << "before " << (node-sOffset) << endl;  // debugging
      clusters[cluster-sOffset]->add(node-sOffset);  // add node to the cluster
      //clusters[cluster-sOffset]->display();  // debugging
      //cout << "after " << (cluster-sOffset) << endl;  // debugging
      //node = -1;  cluster = -1;                    // reset to invalid values
    } // end if node, cluster, weight
    else if(!din.eof())
           errorMsg("Invalid node membership in inputLFRANS? (node "
                    +itos(node)+" and cluster "+itos(cluster)+")");
                    
    //cout << "Current nNodes = " << nodeCount << endl;  // debugging
    //clusters[cluster-sOffset]->display();  cout << endl; // debugging
    //din.getline(buffer,MaxLineLength); // otherwise skip rest of line
  } // end while !eof
  // now set the actual determined number of clusters q
  nClusters = qMax;

  // perform consistency check
  if(nodeCount!=N) {
        cout << red << "  Warning:  Number of nodes does not match the problem "
             << "definition in LFR cluster list read.\n  File specified " 
             << nodeCount << " nodes, but the number of nodes was given as " 
             << N << normText << endl;
        din.close();
        return 2;  // number of nodes does not match problem definition
  } // end edges check
  #ifdef DEBUG_MODE
  // test actual number of nodes within established qMax
  nodeCount = 0;
  for(int i=0; i<qMax; i++)  nodeCount += clusters[i]->getn();
    cout << "N_calc = " << nodeCount << "... ";
  #endif

  //cout << grey << "done reading file with q = " << nClusters << " clusters." 
  //     << normText << endl;
  //display();
  //sortClusters();
  initNodesList();  // update node-cluster list - innefficient
  // successfully finished reading LFR cluster list data file - close and exit
  din.close();
  return 1;
}; // end inputLFRANS


int TClusterList::inputRaw(string fname, int N, int nodeOffset) {
  int  iSize, iNode, MaxLength = 10001;
  // track current number of nodes and the max node label 
  // (remember to correct for a file node offset of 1)
  int nTotal = 0, nodeCount = 0, maxNode = -1, q = 0;  
  char buffer[MaxLineLength];  // char* variable for getline() function

  if(nodeOffset<0)  
    errorMsg("Node offset not implemented for negative values in inputRaw().");
  ifstream fin(fname.c_str());
  if(fin.bad()) {
    errorMsg("Input filename "+fname+" does not exist in inputRaw()!");
    return -1;  // Error filename does not exist
  } // if din.bad

  //cout << red << "Reading raw cluster data in inputRaw()...\n";  // debugging
  erase();
  if(NMaxClusters!=N)  resize(N);
  nClusters = N;  // current q is set to max and corrected at end
  
  while(!fin.eof() && q<N) {
    nodeCount = 0;
    //cout << cyan << "q_i = " << q << ":  " << green;  // debugging

    // the '\r' check is for windows text files with a "\r\n" line termination
    while((fin.peek()!='\n' && fin.peek()!='\r') && nodeCount<MaxLength) {
      fin >> iNode;  // read integer node
      //cout << iNode << " ";  // debugging
      if(iNode>maxNode)  maxNode = iNode; // track max node label

      clusters[q]->add(iNode-nodeOffset); // correct node starting at 0

      //iNode = -1;  //set iNode to invalid control value
      // skip trailing spaces, this makes the newline detection consistent
      while(fin.peek()==' ')  fin.get();           
      nTotal++;                           // track number of nodes (error check)
      nodeCount++;                        // track n for this cluster
    } // end while eol check

    //fin.getline(buffer,MaxLineLength);    // skip to next line (just reads
    // skip trailing spaces, this makes the eof detection consistent
    fin.get();  // skip to next line (just reads the ending \n)
    fin.peek(); // crude peek for eof, I think this triggers eof state at end
    //cout << "\n"; // debugging
    q++;
  } // end while eof check
  
  // max final size assignments
  NNodes = nTotal;
  nClusters = q;   // communities after this are zero
  initNodesList(); // now set internal node-to-cluster list

  // calculate associated system energy 
  //display(1,0,"Raw input from inputRaw");  
  //cout << normText << "done.\n";  // debugging
  fin.close();
  
  // do some error checking before exiting
  int nZero;
  nZero = clearZeroClusters();
  if(nZero!=0)  warningMsg("Empty clusters were found after inputRaw()?");
  if(nTotal!=N)
    errorMsg("Counted Nc = "+itos(nTotal)+" does not match expected N = "
             +itos(N)+" in inputRaw()?");
  if((maxNode-nodeOffset+1)!=N)  
    errorMsg("Max node does not match expected N in inputRaw()?");
  
  return 1;
}; // end inputRaw cluster list


int TClusterList::outputRaw(string fname, int nodeOffset, bool bSort) {
  int nTotal = 0, nodeCount = 0, maxNode = -1, q = 0;  
  char buffer[MaxLineLength];  // char* variable for getline() function

  if(nodeOffset<0)  
    errorMsg("Node offset not implemented for negative values.");

  if(bSort)  sortClusters();
  ofstream fout(fname.c_str());  // open and overwrite existing file
  //if(fout.bad()) {
  //  errorMsg("Input filename "+fname+" does not exist in inputRaw()!");
  //  return -1;  // Error filename does not exist
  //} // if din.bad

  cout << red << "Writing raw cluster data in outputRaw... ";  // debugging

  for(int i=0; i<nClusters; i++) {
    clusters[i]->begin();  // begin manual iteration
    
    for(int j=0; j<clusters[i]->getn(); j++) {
      // output this node
      fout << (clusters[i]->getCurrentNode()+nodeOffset) << " ";  
      clusters[i]->next();
    } // end for j
    fout << "\n";  // debugging
  } // end for i

  cout << red << "done." << endl;  // debugging
  fout.close();
  return 1;
}; // end outputRaw cluster list


// --- TClusterListArray functions --------------------------------------------
TClusterListArray::TClusterListArray(int MaxLists, int listSize) {
  NMaxLists = MaxLists;
  lists = new TClusterList*[NMaxLists];
  //for(int i=0; i<NMaxLists; i++)  lists[i] = NULL;
  // fill the array with empty lists using default constructor
  for(int i=0; i<NMaxLists; i++) {
    lists[i] = new TClusterList;
    if(listSize>0)  lists[i]->resize(listSize);
  } // end for i
  // if the user specifies a size for the empty lists, then we assume that he
  // wants random access to them as valid empty lists
  if(listSize==-1)  nLists = MaxLists;  
  else              nLists = listSize;
}; // end default constructor

/*
TClusterListArray::~TClusterListArray() { 
  for(int i=0; i<NMaxLists; i++)  delete lists[i];
  delete[] lists;
}; // destructor
*/
TClusterListArray::~TClusterListArray() { 
  destroy(); 
}; // end destructor


int TClusterListArray::getMinEIndex() const {
  int minE = BigInteger;
  int minEIndex = -1;
  for(int i=0; i<nLists; i++) 
    if(lists[i]->getE()<minE) { minE = lists[i]->getE();  minEIndex = i; }
  #ifdef DEBUG_MODE
  if(minEIndex<0)
    warningMsg("Attempting to get the minimum energy list on an empty array?");
  #endif
  return minEIndex; 
};  // end 'fast' is NZ calculation


TClusterList* TClusterListArray::getMinEList() const {
  // returns a pointer to the list with the lowest energy
  int minE = BigInteger;
  int minEIndex = -1;
  for(int i=0; i<nLists; i++)
    if(lists[i]->getE()<minE) { minE = lists[i]->getE();  minEIndex = i; }
  #ifdef DEBUG_MODE
  if(minEIndex<0)
    warningMsg("Attempting to get the minimum energy list on an empty array?");
  #endif
  if(minEIndex>=0)  return lists[minEIndex]; 
  return NULL;  // otherwise it is an invalid result
};  // end 'fast' is NZ calculation


int TClusterListArray::getMaxQIndex() const {
  double maxQ = -1.0;
  int maxQIndex = -1;
  for(int i=0; i<nLists; i++)
    if(lists[i]->getQ()>maxQ) { maxQ = lists[i]->getQ();  maxQIndex = i; }
  #ifdef DEBUG_MODE
  if(maxQIndex<0)
    warningMsg("Attempting to get the minimum energy list on an empty array?");
  #endif
  return maxQIndex; 
};  // end 'fast' is NZ calculation


TClusterList* TClusterListArray::getMaxQList() const {
  // returns a pointer to the list with the lowest energy
  double maxQ = -1.0;
  int maxQIndex = -1;
  for(int i=0; i<nLists; i++)
    if(lists[i]->getQ()>maxQ) { maxQ = lists[i]->getQ();  maxQIndex = i; }
  #ifdef DEBUG_MODE
  if(maxQIndex<0)
    warningMsg("Attempting to get the minimum energy list on an empty array?");
  #endif
  if(maxQIndex>=0)  return lists[maxQIndex]; 
  return NULL;  // otherwise it is an invalid result
};  // end 'fast' is NZ calculation

// ---------------------------------------------------------------------------
} // end namespace MSL
// ---------------------------------------------------------------------------
