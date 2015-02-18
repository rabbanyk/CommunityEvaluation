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

#ifndef MSL_TLIST_H
#define MSL_TLIST_H

#include <iostream>
#include <iomanip>
#include <string>

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "MSL_Array1D_Template.h"
#endif

// need to prototype class Array1D since it is used in assignments
//template <typename T> class Array1D;

using namespace std;

namespace MSL {

template <typename T> class TListNode;
template <typename T> class TList;
template <typename T> istream& operator>>(istream& fin, TList<T>& a);

#ifndef MSL_TLIST_T
// define some common vector types
typedef TList<double> TListFloat;// double  array type definition
typedef TList<int>    TListInt;  // integer array type definition
//typedef Array1D<bool>   TVectorBool; // integer array type definition
#define MSL_TLIST_T
#endif

template <typename T> 
class TListNode {
 friend class TList<T>;
 public:
  // overloaded operators
  // Public member functions
  // Standard public member functions
  TListNode<T>(T node=0) { n = node; };            // bare constructor
  ~TListNode<T>() {};                                // destructor
  // everything is public here
  // Data elements
 protected:
  T             n;      // current number Lists stored
  TListNode<T> *pNext;  // maximum number stored Lists allowed
}; // end class TList


//------------------------------------------------------------------------------
template <typename T> 
class TList {
// Memory allocation is dynamic and can use standard C++ array indexing with
// brackets as in arrayName[index] with zero indexing.
// friend istream&   operator>>(istream& fin, TList& a);

 public:
 // overloaded operators
           TList<T>&  operator=(const TList& b);
           //bool     operator==(const TList& b);
  // C array reference for a linked-list, not efficient!
           T&         operator[](int i);
  // Public member functions
  // initialize a List given a List or integer list
           void       init(TVectorInt &v, bool bRandom=1);
           void       init(int nodeList[], int size);
  // add a node to the end, only the dE/dQ version updates properties
  // add node by detailed parameters, node changes return number of nodes left
           int        add(T node);      // add node to this List by value
           //int      add(T &node);     // add node to this List
           int        add(TList &List); // add a List to this List
  // 'move' functions only manipulate pointers; they do not allocate memory
           int        moveCurrent(TList &c); // move current node *to* c
           int        move(TList &c);        // move this List *to* c List
  // erase function delete the object and return the memory
           int        eraseCurrent();  // erase current node and return memory
           void       erase();         // erase all data and return memory
  // return node integer value for pCurrent
  inline   T&         getStart() const   { return pBegin->n;   };
  inline   T&         get(int i)         { return (*this)[i];  };
  // return node by reference for pCurrent
  inline   T&         getCurrent() const { return pCurrent->n; };
  // peak at next element, but do not increment pCurrent pointer
  inline   T&         getNext()    const { return pCurrent->pNext->n; };
  // other functions on the List
           void       sort();   // performs an insertion sort on the nodes
           bool       isMember(T n)    const; // returns whether n is in this TList
           // erases member if it exists and returns true, otherwise returns false
           bool       eraseMember(T n); 
           T*         getMember(T n) const; // returns pointer to n if it is in TList
           void       display(int offset=0) const; // display List to console
           void       copyToVector(TVectorInt &v);
  // iteration loop functions simulate array-like iteration through the list
  // begin() is required which starts a pCurrent pointer which will step through 
  // the list as next() is called.  It is only modified by user calls except for
  // some consistency updates.  There is only one iteration loop at a time, but
  // interior functions do not directly utilize the variables.
  inline   T&         begin();       // start an iteration
  inline   T&         next();        // step to next node in an iteration
  inline   bool       isValid() const { return 1; } // debugging function
  inline   bool       isValid() { return 1; } // debugging function
  // data acquisition functions
  inline   int        getSize() const    { return nNodes;               };

  // Standard public member functions
  TList();                          // bare constructor
  TList(int nodeList[], int size);  // constructor with integer array
  ~TList();                         // destructor
  
 protected:                 // data elements are protected
  int            nNodes;    // current number Lists stored
  //int          MaxNodes;  // maximum number stored Lists allowed
  TListNode<T>  *pBegin;    // pointer to beginning of data - NULL if no data
  TListNode<T>  *pEnd;      // pointer to end of data list - NULL if no data
  TListNode<T>  *pCurrent;  // pointer to the current node in an iteration
  TListNode<T>  *pPrevious; // pointer to the previous node in an iteration
  TListNode<T>  *piCurrent; // pointer to the current node in op[] iteration
  int            iCurrent;  // internal index for operator[] function which provides
                            //  simulated array referencing as an alternative to
                            //  the single pointer iterator
  // restricted due to internal changes, but it needs to be understood better
  //T&             operator[](int i) const;  
  bool           eraseElement(TListNode<T> *pPrevious, TListNode<T> *pElement); 
}; // end class TList;


// --- inline List functions ------------------------------------------------
template <typename T> 
inline T& TList<T>::begin() {
  // starts an iteration loop, or alternatively just returns the starting node
  // The iteration loop is an 'array' equivalent mechanism to provide current 
  // node in a loop.
  // The iteration loops are for user usage and are not incremented within the 
  // classes.  They are only updated on user changes in appropriate function.
  pPrevious = NULL;
  pCurrent  = pBegin;  // pBegin is NULL if TList is empty
  if(nNodes<1 || pBegin==NULL)
    //warningMsg("Attempting to begin iteration on a size zero List.");
    errorMsg("Attempting to begin iteration on a size zero TList.");
  return pCurrent->n;
};  // begin pointer iteration


template <typename T> 
inline T& TList<T>::next() {
  // starts an iteration loop, or alternatively just returns the starting node
  if(pCurrent!=NULL) {
    // also handles case of moved node where pPrevious is set equal to pCurrent
    pPrevious = pCurrent;    // simply redundant if a node has been moved
    pCurrent  = pCurrent->pNext;
  } // end if
  // if pCurrent is null and List size > 0, reset to beginning
  else if(nNodes>0) { 
   errorMsg("pCurrent!=NULL on a non-zero (n = "+itos(nNodes)+")TList?");
   pPrevious = NULL;  pCurrent = pBegin; 
  } else errorMsg("Attempting to use iterate with next() on size zero TList");
  //else  warningMsg("Attempting to iterate with pCurrent past end of list.");
  return pCurrent->n;
};  // iterate to next node on pointer iteration


//------------------------------------------------------------------------------
//--------------------- TList Member Function Declarations-----------------
template <typename T> 
TList<T>::TList() {
  pBegin    = NULL;  pEnd = NULL;
  nNodes    = 0;
  pCurrent  = NULL;  pPrevious = NULL;  // user pointer iteration parameters
  piCurrent = NULL;  iCurrent = 0;      // operator[] parameters
};  // default constructor

template <typename T> 
TList<T>::TList(int nodeList[], int size) {
  pBegin = NULL;     pEnd = NULL;
  nNodes = size;
  pCurrent  = NULL;  pPrevious = NULL;  // user pointer iteration parameters
  piCurrent = NULL;  iCurrent = 0;      // operator[] parameters
  if(size>0) {
    pBegin = new TList;
    pBegin->init(nodeList[0]);
    TListNode<T> *pn = pBegin;
    for(int i=1; i<size; i++) {
      pn->pNext = new TList;
      pn->init(nodeList[i]);
      pn = pn->pNext;
    } // end for i
    pEnd = pn; // finally set end pointer to last nodes in list
  } else { 
    pBegin = NULL;  pEnd = NULL; 
    warningMsg("Attempting to construct List with size zero array");
  } // end else
};  // integer list constructor


template <typename T> 
TList<T>::~TList() {
  TListNode<T>  *pn, *pc = pBegin;
  // traverse and delete linked list
  for(int i=0; i<nNodes; i++) { 
    pn = pc->pNext;  
    #ifdef DEBUG_MODE
    if(pc==NULL)  errorMsg("~TList<T> has premature terminating NULL pointer?");
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


template <typename T> 
void TList<T>::copyToVector(TVectorInt &v) {
  v.resize(nNodes);
  if(nNodes>0) {
    TListNode<T> *pn = pBegin;
    for(int i=0; i<nNodes; i++) { v[i] = pn->n;  pn = pn->pNext; }
  } // end if nNodes
  return;
};  // copy to vector int


template <typename T> 
void TList<T>::init(int nodeList[], int size) {
  nNodes   = size;
  pCurrent = NULL;  pPrevious = NULL;
  if(size>0) {
    pBegin = new TList;
    pBegin->init(nodeList[0]);
    pEnd = pBegin;
    TListNode<T>  *pn = pBegin;
    for(int i=1; i<size; i++) {
      pn->pNext = new TListNode<T>;
      pn->init(nodeList[i]);
      pn = pn->pNext;
    } // end for i
    pEnd = pn; // finally set end pointer to last nodes in list
  } else { 
    pBegin = NULL;  pEnd = NULL; 
    warningMsg("Attempting to initialize TList with size zero array");
  } // end else

  return;
};  // init

template <typename T> 
void  TList<T>::init(TVectorInt &v, bool bRandomOrder) {
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

  // now finally add the vector to the List
  for(int k=0; k<N; k++)  add((*pv)[k]);
  //cout << s << endl; // debugging

  return;
}; // end init


template <typename T> 
bool  TList<T>::isMember(T v) const {
  // returns a boolean value indicating whether n is a member of this list
  bool bFound = 0;  
  if(nNodes>0) {
    //cout << red << "\nStarting a member test for cluster = " << v 
    //            << "... i = 0, "; // debugging
    TListNode<T>  *pn = pBegin;
    bFound = (pn->n == v);  
    int i=1;
    while(i<nNodes && !bFound) {
      //cout << i << ", " << flush; // debugging
      pn = pn->pNext;
      bFound |= (pn->n == v);
      i++;  // increment i
    } // end while i
    //cout << "done." << normText << endl; // debugging
  } // end if nNodes

  return bFound;
}; // end isMember


template <typename T> 
bool  TList<T>::eraseMember(T v) {
  // returns a boolean value indicating whether n is a member of this list
  bool bFound = 0;  
  if(nNodes>0) {
    //cout << red << "\nStarting a member test for cluster = " << v 
    //            << "... i = 0, "; // debugging
    TListNode<T>  *pn = pBegin, *pp = pn;
    //bFound = (pn->n == v);  
    int i = 0;
    do {
      bFound |= (pn->n == v);
      // now erase it if the member was found and exit
      if(bFound)  eraseElement(pp,pn);

      pp = pn;  // track previous node for erase operation
      //cout << i << ", " << flush; // debugging
      pn = pn->pNext;
      i++;  // increment i
    } while(i<nNodes && !bFound); // end while i
    //cout << "done." << normText << endl; // debugging
  } // end if nNodes

  return bFound;
}; // end eraseMember


template <typename T> 
bool  TList<T>::eraseElement(TListNode<T> *pPrev, TListNode<T> *pCurr) {
  // based upon the passes previous and current element pointers, erase the
  // pCurr element from the list
  // NOTE that the pCurrent, etc. iteration pointers and counters are *not*
  // updated during the function call

  // some error and consistency checks
  if(pPrev==NULL || pCurr==NULL)
    errorMsg("Passed null pointer in TList::eraseElement()?");
  // if pCurr points at the beginning of the TList, require pPrev==pCurr or
  // pPrev==NULL to be a valid function call
  if(pCurr==pBegin && !(pPrev==pCurr || pPrev==NULL))
    errorMsg("Inconsistent beginning passed pointers in TList::eraseElement()");
  if(pPrev->pNext!=pCurr && !(pPrev==pCurr || pPrev==NULL))
    errorMsg("pPrev next node doesn't point at pCurr in TList::eraseElement()");

  //T *pn;  pn = pCurr;  // save a pointer to current element
  // now erase the element identified by the passed pointers
  if(nNodes>1) {
    // In each case we set pCurrent backwards (or to NULL) as an indicator of
    // the change.  It also maintains consistency with the iteration loop.
    if(pCurr==pBegin) {
      // pCurrent is the beginning of the list
      pBegin = pCurr->pNext;
    } // end if at start
    else if(pCurr==pEnd) {
      // pCurrent is the end of the list
      pEnd = pPrev;  pPrev->pNext = NULL; // pPrevious is now the end
    } // end if at end
    else {
      // pCurrent is somewhere in the middle
      pPrev->pNext = pCurr->pNext;        // shortcut around current node
    } // end else
    nNodes--;
    delete pCurr;                         // finally return the memory
    //pCurr = NULL;  pPrev = NULL;        // just in case - pCurr is invalid now
  } // end if more than one node
  else if(nNodes==1) {
    pBegin = NULL;  pEnd = NULL;
    nNodes = 0;
    delete pCurr;                         // finally return the memory
    //pCurr = NULL;  pPrev = NULL;        // just in case - pCurr is invalid now
  } // end else if nNodes
  else  warningMsg("Attempting to erase a node from a size zero TList?");
}; // erase element

  
template <typename T> 
T*  TList<T>::getMember(T v) const {
  // returns a boolean value indicating whether n is a member of this list
  bool bFound = 0;  
  if(nNodes>0) {
    //cout << red << "\nStarting a member test for cluster = " << v 
    //            << "... i = 0, "; // debugging
    TListNode<T>  *pn = pBegin;
    bFound = (pn->n == v);  
    int i=1;
    while(i<nNodes && !bFound) {
      //cout << i << ", " << flush; // debugging
      pn = pn->pNext;
      bFound |= (pn->n == v);
      i++;  // increment i
    } // end while i
    //cout << "done." << normText << endl; // debugging
    if(bFound)  return &(pn->n);
    else        return NULL;
  } // end if nNodes
  else return NULL;
}; // end getMember


template <typename T> 
TList<T>& TList<T>::operator=(const TList& c) {
  //errorMsg("Oh, we use this");  // debugging
  if(nNodes>0)  erase();  // delete current data
  nNodes   = c.nNodes;
  //pPrevious = NULL;

  if(c.nNodes>0) {
    TListNode<T>  *pc = c.pBegin;
    pBegin    = new TListNode<T>;
    pBegin->n = c.pBegin->n;            // copy node data
    TListNode<T>  *pn = pBegin;
    for(int i=1; i<c.nNodes; i++) {
      pn->pNext = new TListNode<T>;
      pn = pn->pNext;  pc = pc->pNext;  // step to next node on both
      pn->n = pc->n;                    // copy node data
    } // end for i
    pn->pNext = NULL;                   // set end pointer to NULL
    pEnd      = pn;
  } // end if nNodes
  else { pBegin = NULL;  pEnd = NULL; }

  piCurrent = pBegin;  iCurrent = 0;
  return *this;
};  // equals operator

/*
template <typename T> 
T&  TList<T>::operator[](int i) const {
  // an array[] reference operator for the linked-list implementation
  // obviously, it is not efficient if i>0.
  // Use the iteration pointers for an array-like iterative implementation 
  // as long as the iteration is relatively sequential through the list.
  // Random data access will be O(n).
  // Things can get messed up if there are deletions at the point of the 
  // iCurrent pointer which often occurs on a node move... be careful.
  // I did not want to bury more special cases in the rest of the code.
  // It is best used in a closed short loop starting at node zero.
  if(i>=nNodes)  errorMsg("List op["+itos(i)+"] reference is out of bounds");
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
  if(pBegin==NULL)
    errorMsg("List operator["+itos(i)+"] detected pBegin==NULL!? " +
             "(with iCurrent = "+itos(iCurrent)+")");
  else if(piCurrent==NULL)  
    errorMsg("List operator["+itos(i)+"] is returning a NULL!?  " + 
             "(with iCurrent = "+itos(iCurrent)+")");
  return piCurrent->n;
}; // end operator[] for linked-list implementation
*/
template <typename T> 
T&  TList<T>::operator[](int i) {
  // an array[] reference operator for the linked-list implementation
  // obviously, it is not efficient if i>0.
  // Use the iteration pointers for an array-like iterative implementation 
  // as long as the iteration is relatively sequential through the list.
  // Random data access will be O(n).
  // Things can get messed up if there are deletions at the point of the 
  // iCurrent pointer which often occurs on a node move... be careful.
  // I did not want to bury more special cases in the rest of the code.
  // It is best used in a closed short loop starting at node zero.
  if(i>=nNodes)  errorMsg("List op["+itos(i)+"] reference is out of bounds");
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
  #ifdef DEBUG_MODE
  if(pBegin==NULL)
    errorMsg("List operator["+itos(i)+"] detected pBegin==NULL!? " +
             "(with iCurrent = "+itos(iCurrent)+")");
  else if(piCurrent==NULL)  
    errorMsg("List operator["+itos(i)+"] is returning a NULL!?  " + 
             "(with iCurrent = "+itos(iCurrent)+")");
  #endif
  return piCurrent->n;
}; // end operator[] for linked-list implementation

/*
bool TList<T>::operator==(const TList& c) {
  bool isEqual = 1;  // default state is true

  return isEqual;
};  // constructor
*/


template <typename T> 
void TList<T>::sort() {
  // Uses an optimized insertion sort to sort a generic integer data array.
  // implementation is a linked list, so it simplifies the insertion logic
  // we pick the next node in the list and scan previous nodes.  Then 
  // swap pointers at the correct location.
  int     i, j;
  T       nI, nJ;                 // temporary integer value
  TListNode<T>  *pJ, *pI;         // pointers to i'th, j'th
  TListNode<T>  *pPrevI, *pPrevJ; // pointer to previous i'th and j'th numbers
  bool    bFound = 0;

  pI = pBegin;                    // starting location
  for(i=1; i<nNodes; i++) {
    // pick up next to test insert in earlier part of list
    pPrevI = pI;                  // pointer to previous i'th node for insertion
    pI     = pI->pNext;  nI = pI->n;
    if(nI<pPrevI->n) {            // only check for move if out of place
      pJ   = pBegin;     nJ = pJ->n;
      j    = 0;          bFound = 0;
      while(j<i && !bFound) {
        if(nJ>=nI)  bFound = 1;
        else { 
          pPrevJ = pJ;            // store previous position for insertion
          pJ = pJ->pNext;    nJ = pJ->n;
          j++;           // if not incrementing j, loop is stopped by bFound
        } // end else
      } // end while j
    
      // finally swap List nodes if needed
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
} // end sort


template <typename T> 
int TList<T>::add(T v) {
  // perform a pure addition including allocation of a new TList node
  // other List properties are *not* updated
  // no need to updata op[] stuff for target list since we add to end
  //cout << "Here in add with v = " << v << endl;  // debugging
  TListNode<T> *pn;  
   pn = new TListNode<T>(v);
  //pn->n       = v;          // copy node data
  pn->pNext   = NULL;       // update new end pointer to NULL
  if(nNodes==0) {
    pBegin    = pn;         // if List is currently empty, point pBegin
    piCurrent = pn;  iCurrent = 0;
  } // end if nNodes
  else  pEnd->pNext = pn;   // set previous end to point to added node
  pEnd = pn;                // update pEnd to point to added node
  nNodes++;
  return nNodes;
};  // add node

/*
template <typename T> 
int TList<T>::add(T &s) {
  // perform a pure addition including allocation of a new List node
  // other List properties are *not* updated
  TList *pn;  pn = new TList;
  pn->n       = s;    // copy node data
  pn->pNext   = NULL; // update new end pointer to NULL
  if(nNodes ==0) {
    pBegin    = pn;   // if List is currently empty, point pBegin
    piCurrent = pn;
  } // end if nNodes
  else            pEnd->pNext = pn;   // set previous end to point to added node
  pEnd = pn;      // update pEnd to point to added node
  nNodes++;
  return nNodes;
};  // add node
*/

template <typename T> 
int TList<T>::add(TList &c) {
  // add a whole List at a time by allocating memory and copying all data
  //if(c.getSize()>0) { // is redundant even if c.pBegin is NULL (or anything)
  TListNode<T> *pn = c.pBegin;
  for(int i=0; i<c.getSize(); i++) {
    add(pn->n);      // add this node
    pn = pn->pNext;  // move to next node in list, still ok if last node
  } // end for i
  //} // end if
  return nNodes;
};  // add a List


template <typename T> 
int TList<T>::moveCurrent(TList &c) {
  // moves current node in iteration to List c without deleting or 
  // allocating anything.  Current node must be a valid List node.
  // The iteration loop is not invalidated and a subsequent next() will 
  // correctly increment to the next node in the list.

  TListNode<T> *pn = pCurrent;  // store current List node for addition below
  // remove it from old list
  if(nNodes==1)   {
    pBegin    = NULL;  pEnd     = NULL; 
    nNodes    = 0;
    piCurrent = NULL;  iCurrent = 0;
    pPrevious = NULL;  pCurrent = NULL;     // just in case
  } else if(pCurrent==pBegin) { // pCurrent is at the beginning of the list
    pBegin    = pCurrent->pNext;
    pPrevious = NULL;    pCurrent = NULL;
    piCurrent = pBegin;  iCurrent = 0;
    nNodes--;                   // decrement number of nodes in either case
  } else if(pCurrent==pEnd) {   // pCurrent is at the end of the list
    //cout << redText << "ok, moving end pointer now" << normText;  // debugging
    pEnd             = pPrevious;
    pPrevious->pNext = NULL;
    pCurrent         = pPrevious;  // indicator
    piCurrent = pBegin;  iCurrent = 0;
    nNodes--;          // decrement number of nodes in either case
    if(pCurrent->pNext!=NULL)   // debugging - error detection
      errorMsg("List move does not show pNext=NULL pointer at end of list.");
  } else if(nNodes==0) {  // debugging
      errorMsg("Attempting to move a node from a size zero TList???");
  } else {                      // somewhere in the middle of the list
    //pCurrent         = pCurrent->pNext;  
    //pPrevious->pNext = pCurrent->pNext; 
    //cout << redText << "ok, moving the pointer now" << normText;  // debugging
    pPrevious->pNext = pCurrent->pNext; 
    pCurrent         = pPrevious;  
    piCurrent = pBegin;  iCurrent = 0;
    nNodes--;          // decrement number of nodes in either case
  } // end else

  // add it to new list
  pn->pNext = NULL;
  if(c.nNodes==0) { c.pBegin = pn;       c.pEnd = pn; } // if c is empty
  else            { c.pEnd->pNext = pn;  c.pEnd = pn; } // end else
  // no need to updata op[] stuff for target list since we add to end of c list
  c.nNodes++;

  // finally set pCurrent to pPrevious as an indicator which also maintains
  // consistency with the iteration loop
  pCurrent = pPrevious;  
  return nNodes;
};  // end moveCurrent with no energy updates


template <typename T> 
int TList<T>::move(TList &c) {
  // moves this List to List 'c' without deleting or allocating anything.
  // no need to updata op[] stuff for target list since we add to end of c list
  if(c.getSize()>0) {
    c.pEnd->pNext = pBegin;  c.pEnd = pEnd;  c.nNodes += nNodes; // move to c
  } else {
    // c is size zero
    //warningMsg("Moving to a size zero list");
    c.pBegin = pBegin;  c.pEnd = pEnd;  c.nNodes = nNodes;
  } // end else
  pBegin = NULL;     pEnd = NULL;   nNodes = 0;     // remove from this
  piCurrent = NULL;  iCurrent = 0;
  return c.nNodes;
}; // move to another list


template <typename T> 
void TList<T>::erase() {
  TListNode<T> *pn, *pc = pBegin;
  // delete linked list
  for(int i=0; i<nNodes; i++) { pn = pc->pNext;  delete pc;  pc = pn; }
  pBegin    = NULL;   pEnd      = NULL;
  pCurrent  = NULL;   pPrevious = NULL;
  piCurrent = NULL;   iCurrent  = 0;
  nNodes    = 0;
  return;
}; // end erase


template <typename T> 
int  TList<T>::eraseCurrent() {
  // delete current node in linked list and return memory, but does not destruct
  // the object.  pCurrent must be a valid node.
  TListNode<T> *pn = pCurrent;
  if(nNodes>1) {
    // In each case we set pCurrent backwards (or to NULL) as an indicator of
    // the change.  It also maintains consistency with the iteration loop.
    if(pCurrent==pBegin) { // pCurrent is the beginning of the list
      //cout << "Ok, here we are... " << (int)pBegin << " " << (int)pCurrent 
      //     << " " << (int)(pCurrent->pNext) << " " 
      //     << (int)(pCurrent->pNext->pNext) 
      //     << " " << pCurrent->n << " " << pCurrent->pNext->n << " "
      //     << pCurrent->pNext->pNext->n << " " << flush;
      pBegin = pCurrent->pNext;
      //cout << "1... " << flush;
      pCurrent = NULL;  //pPrevious = NULL;     // pCurrent is now invalid
      //cout << "2... " << flush;
    } // end if at start
    else if(pCurrent==pEnd) { // pCurrent is the end of the list
      pEnd = pPrevious;  pPrevious->pNext = NULL;  // pPrevious is now the end
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
    pBegin = NULL;     pEnd = NULL;
    pCurrent = NULL;   pPrevious = NULL;
    piCurrent = NULL;  iCurrent = 0;          // update array[] 'iterator' data
    nNodes = 0;
    delete pn;                                // finally return the memory
  } // end else if nNodes
  else  warningMsg("Attempting to erase a node from a size zero TList?");
  // now finish with upkeep and memory
  //cout << "3... " << flush;
  //cout << "3b " << flush << (int)pBegin << " " << (int)pCurrent << " " 
  //         << (int)(pCurrent->pNext) << " " 
  //         << (int)(pCurrent->pNext->pNext) << " " << flush;
  //piCurrent = pBegin;  iCurrent = 0;        // update array[] 'iterator' data
  //cout << "4... " << flush;
  //delete pn;  // finally return the actual memory
  //cout << "done with nNodes = " << flush << nNodes << endl;  // debugging

  return nNodes;
}; // end eraseCurrent


template <typename T> 
void TList<T>::display(int nodeOffset) const {
  // crude community list omitting edge connections (no graphical capability)
  //int MaxEdges = nNodes*(nNodes-1)/2;

  // output row label
  cout << brown << "Size ";
  if(nNodes<10) cout << " ";  // crude alignment for less than 100 nodes/group
  cout << nNodes << ":  " << green;

  if(nNodes==0)  return;
  // output List nodes
  TListNode<T> *pn = pBegin;
  cout << green;
  for(int i=0; i<nNodes; i++) {
    cout << pn->n << " ";
    pn = pn->pNext;  // increment to next node in list
    //cout << "here e" << endl;
  } // end for i
  cout << normText; // << "\n";

  return;
};  // display List


// now define other non-member related functions
// --- related functions -------------------------------------------------------
template <typename T> 
inline ostream& operator<<(ostream& fout, TList<T>& a) {
  // Outputs the List utilizing the display function
  a.display();
  return fout;
};

// ---------------------------------------------------------------------------
} // end namespace MSL
// ---------------------------------------------------------------------------

#endif
