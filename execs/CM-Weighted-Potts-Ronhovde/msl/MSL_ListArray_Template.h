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

Specific Notes:
As a science application, we do *not* perform extraneous redundant checks for
a valid input parameter, but there is a pre-processor parameter DEBUG_MODE
that will enable substantial error checking.

To-Do List:
1. 
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TLISTARRAY_H
#define MSL_TLISTARRAY_H

#include <iostream>
#include <iomanip>
#include <string>

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
#ifndef MSL_TLIST_H
#include "MSL_List_Template.h"
#endif

using namespace std;

namespace MSL {

template <typename T> class TListNode;
template <typename T> class TList;
template <typename T> class TListArray;

#ifndef MSL_TLISTARRAY_T
#define MSL_TLISTARRAY_T
// define some common types
typedef TListArray<double> TListArrayFloat;// double  array type definition
typedef TListArray<int>    TListArrayInt;  // integer array type definition
#endif

//------------------------------------------------------------------------------
template <typename T>
class TListArray {
// Memory allocation is dynamic and uses standard C++ array indexing with 
// brackets as in arrayName[index] with zero indexing.

 public:
 // overloaded operators
           TListArray<T>& operator=(const TListArray<T>& b);
  inline   TList<T>&      operator[](int i) const { return *(lists[i]); };
  inline   TList<T>&      operator[](int i)       { return *(lists[i]); };
  // Public member functions
           // initialize a list list as evenly as possible given N and q
           void   initEven  (int N, int &q, bool bRandomOrder=1, int iBegin=0);
           // initialize a list list given a numerical vector
           void   initEven  (Array1D<T> &v, int &q, bool bRandomOrder=1);
           // take a list and evenly initialize a new list with it into q groups
           void   initEven  (TList<T> &c, int &q, bool bRandomOrder=1);
           // take the list and evenly re-initialize it into q groups
           void   initEven  (int q, bool bRandomOrder=1);
           // initialize a list completely random within q groups
           void   initRandom(Array1D<T> &v, int &q);
           void   initRandom(int N, int &q, int iBegin=0);
           // add and move functions
  inline   int    addEmpty();                 // add an empty list to the end
           int    add(TList<T> &c);           // add a list to the end by copy
           int    add(TListArray<T> &c);      // add a list to end by copy
           int    move(TListArray<T> &c);     // move a list to end by pointer
  inline   void   move(int i, int j);   // move list i to j and remove i
           // erase all lists, do not return memory, keep empty lists accessible?
           void   erase(bool bKeepEmpty=0); 
           void   resize(int size);     // move list i to j and remove i
           void   destroy();            // destroy all lists, return memory
           int    erase(int i, bool bFast=1); // erase list i, return memory
           //void   resize();                 // resize the number of lists
                                              //   (at the moment it deletes everything!!!)
           void   swap(int i, int j);         // swap lists i and j
  // other list functions
           // display list list to console
           void   display(int offset=0, bool bSort=0, string s="", int verbosity=1);
           //void resize(int Length, bool Fill, T FillValue);
           void   sort();         // pointer based insertion sort on the lists
  inline   bool   isValid(bool bEmptyOK=0) const;   // valid check on object
  // data acquisition functions
  inline   int    getSize()               const { return nLists; };
  inline   int    getnNodes(bool bFast=0) const; // if bFast=1, no recalculation

  // Standard public member functions
  TListArray<T>(int NLists = 0);
  ~TListArray<T>();
  
 //protected:
  // Data elements
  int        NMaxLists;
  int        nLists;       // current number lists stored
  int        nNodes;       // constant total number of nodes in graph
  //int      NodeOffset;   // constant setting at what number the nodes begin
  // the List data is stored as pointer list in order to prevent a lot of data
  //   movement during moving operations
  TList<T>  **lists;      // the array of lists
}; // end class TListArray<T>;

// --- inline functions for listList ----------------------------------------
template <typename T>
inline  void  TListArray<T>::move(int i, int j) {
  // move list i to j and remove i without deleting anything except the 
  // empty list i
  lists[i]->move(*(lists[j]));
  //erase(i);  // redundant or even wrong?
  return;
}; // end inline move

template <typename T>
inline bool TListArray<T>::isValid(bool bEmptyOK) const {
  #ifdef DEBUG_MODE
  bool bResult = 0;
  // need to flesh this out some more and implement within the class!!!
  if(nLists>nNodes)
    warningMsg("Number of lists exceeds number of nodes in TListArray object");
  else   bResult = 1; // if we make it to here, it passed all basic checks
  return bResult;
  #else
  return 1;
  #endif
}; // end isValid


//------------------------------------------------------------------------------
//--------------------- TListArray Member Function Declarations-----------------
template <typename T>
TListArray<T>::TListArray(int size) {
  // a proper constructor definition requires a filled list list which we
  // don't have yet, so...
  // this was a passed parameter, now it is deprecated for the moment
  bool bAddAllEmpty = 1;  
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Entering TListArray constructor... ");
  #endif
  NMaxLists = size;
  nNodes    = 0;
  //cout << red << "here a " << normText << flush;      // debugging
  lists     = new TList<T>*[NMaxLists];
  //cout << red << "here b " << normText << flush;      // debugging
  if(bAddAllEmpty) {
    for(int i=0; i<NMaxLists; i++)  lists[i] = new TList<T>;
    // conceptual problem here - Do we define nLists = 0 or NMaxLists
    // Some applications go either way.  Current applications prefer the latter.
    nLists = NMaxLists;  
  } else {
    for(int i=0; i<NMaxLists; i++)  lists[i] = NULL;
    nLists = 0;    // nothing there yet
  } // end else
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting TListArray constructor.");
  #endif
};  // constructor


template <typename T>
TListArray<T>::~TListArray() { 
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Entering TListArray destructor... ");
  #endif
  for(int i=0; i<nLists; i++) {
    #ifdef DEBUG_MODE
    if(i>=nLists && lists[i]!=NULL)
      errorMsg("~TListArray<T> has non-NULL pointer for i>=nLists?");
    #endif
    delete lists[i];
  } // end for i
  delete[] lists;
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting TListArray destructor.\n");
  #endif
}; // destructor


template <typename T>
void TListArray<T>::destroy() {
  for(int i=0; i<nLists; i++) {
    #ifdef DEBUG_MODE
    if(i>=nLists && lists[i]!=NULL)
      errorMsg("TListArray<T> destroy() has non-NULL pointer for i>=nLists?");
    #endif
    delete lists[i];
  } // end for i
  delete[] lists;
  return;
}; // destroy


template <typename T>
void TListArray<T>::resize(int size) {
  destroy();
  // this was a passed parameter, now it is deprecated for the moment
  bool bAddAllEmpty = 1;  
  #ifdef DEBUG_MODE_HIGH
  debugMsg("Entering TListArray resize()... ");
  #endif
  NMaxLists = size;
  nNodes    = 0;
  //cout << red << "here a " << normText << flush;      // debugging
  lists     = new TList<T>*[NMaxLists];
  //cout << red << "here b " << normText << flush;      // debugging
  if(bAddAllEmpty) {
    for(int i=0; i<NMaxLists; i++)  lists[i] = new TList<T>;
    // conceptual problem here - Do we define nLists = 0 or NMaxLists
    // Some applications go either way.  Current applications prefer the latter.
    nLists = NMaxLists;  
  } else {
    for(int i=0; i<NMaxLists; i++)  lists[i] = NULL;
    nLists = 0;    // nothing there yet
  } // end else
  #ifdef DEBUG_MODE_HIGH
  debugMsg("exiting TListArray resize().");
  #endif
  return;
};  // resize


template <typename T>
inline int TListArray<T>::addEmpty() {
  // add an empty list to the end of the list
  if(nLists==NMaxLists) 
    errorMsg("Cannot add another list in TListArray::addEmpty()!");
  //else { lists[nLists] = new TList<T>;  nLists++; }
  // we eliminated the possibility of NULL pointers rather than empty TLists
  else  nLists++;
  
  return nLists;
}; // add a single list to the end of the list


template <typename T>
void TListArray<T>::display(int nodeOffset, bool bSort, string s, int verbosity){
  // s   - is a string that is displayed just prior to the output
  //       which would typically hold the name of the list set for display
  // verbosity controls how much information is displayed about the system
  //   0 - nothing but the lists themselves
  //   1 - displays... not implemented yet
  //   2 - displays... not implemented yet
  // just display each list in sequence
  if(bSort)  sort();
  // now output information
  if(verbosity>=0) {
    for(int i=0; i<nLists; i++) {
      if(lists[i]->getSize()>0) {
        cout << blueText << "List " << i << ":  ";
        lists[i]->display(nodeOffset);
        cout << "\n";
       } // end if
    } // end for i
  } // end if verbosity
  if(nLists==0) warningMsg(s+"ListArray has size"+itos(nLists),magentaText);
  return;
};  // display list array


template <typename T>
void  TListArray<T>::initEven(int N, int &q, bool bRandomOrder, int iBegin) {
  // initializes the vector either as sequential integers from 0 to N-1
  // bRandomOrder will randomize the order of the integers
  //cout << redText << "here in initEven a" << endl; // debugging
  Array1D<T> v(N);                 // initialize a working vector
  //cout << redText << "here in initEven b" << endl; // debugging
  v.initStep(iBegin,1);            // initialize nodes in consecutive order
  //cout << redText << "here in initEven c" << endl; // debugging
  if(bRandomOrder)  v.randomizeOrder(); // now randomize the node order
  //cout << redText << "here in initEven d" << endl; // debugging
  // random is 0 below because we already randomized the vector if desired
  //cout << redText << "here in initEven e" << endl; // debugging
  initEven(v,q,0);                 // call basic function
  //cout << redText << "here in initEven f" << endl; // debugging
  return;
}; // end initEven given integer range


template <typename T>
void  TListArray<T>::initEven(int q, bool bRandomOrder) {
  // initializes the list list given a single list
  getnNodes();
  if(q>nNodes) {
    q = nNodes;
    warningMsg("q is larger than the number of nodes in the list list!");
  } // end if q
  
  TList<T> b;                   // declare a working list
  for(int i=0; i<nLists; i++)  b.add(*(lists[i]));
  // finally initialize with the temporary list
  initEven(b,q,bRandomOrder);
  return;
}; // end initEven given integer range


template <typename T>
void  TListArray<T>::initEven(TList<T> &b, int &q, bool bRandomOrder) {
  // initialize to an even set of lists given a single list
  // bRandomOrder will randomize the order of the integers
  int N = b.getSize();  // use given size of v
  //Array1D<T> *pv = &b;  // define pointer in case we randomize
  if(q>N) {
    q = N;
    warningMsg("q is larger than N in initRandom(...) list list.  Resetting q = N!");
  } // end if q

  // declare a working node *array* - not efficient but it makes the function 
  // a little safer, keeps us from modifying b if the order is randomized, 
  // and it allows us to avoid declaring TListArray as a friend of TList
  TListNode<T> **a;  a = new TListNode<T>*[b.getSize()];
  for(int i=0; i<b.getSize(); i++) {
    a[i]    = new TListNode<T>;
    *(a[i]) = b[i];  // uses TList operator[] for linked-list list b
  } // end for i

  // now randomize the node order if indicated which is the default action
  if(bRandomOrder) {
    int        jS;
    TListNode<T>  *pTemp;
    // ensure that every node is moved at least once to a random location
    for(int k=0; k<b.getSize(); k++) {      // swap two nodes at a time
      jS    = randomInt(0,b.getSize()-1);  // not worried about duplicates
      pTemp = a[jS];  a[jS] = a[k];  a[k] = pTemp;  // swap node pointers
    } // end for k
  } // end if bRandomOrder

  erase();                        // erase current list list
  bool qdivN = !((bool)(N%q));    // does q divide N? logic is easier
  if(qdivN) {
    int n = N/q;                  // number of nodes per list
    int nCount = 0;               // nCount - track total number all nodes added
    TList<T> c;                   // declare the temporary list
    for(int i=0; i<q; i++) {
      c.erase();                  // erase current list
      // add nodes to current list
      for(int k=0; k<n; k++) { c.add(*(a[nCount]));  nCount++; }
      add(c);                                // add full list to list list
    } // end for i
    // end if qdivN
    
    // consistency check before exit
    if(nCount!=N)  errorMsg("Did not correctly add even nodes in initEven!");
  } else {
    // q does not divide N so initialize as evenly as possible
    int NModq = N%q, nMin = N/q;             // minimum # of nodes per list
    TList<T> c;                
    // need some variables to try to evenly distribute the 'extra' nodes evenly
    // among the lists
    // qProb  - probability of adding an 'extra' spin to a given list
    // nExtra - track total number of 'extra' nodes added during process, 
    //          this should be (N mod q) by the end of the initialization
    // nCount - track total number all nodes added during process
    double qProb = (double)(N%q)/(double)(N/q);  
    int    nCount = 0, nExtra = 0, i = 0;
    while(i<q) {
      c.erase();                             // erase current temporary list
      // add minimum number of nodes to current list
      for(int k=0; k<nMin; k++) { c.add(*(a[nCount]));  nCount++; }
      // for uneven lists, do we add another single node to this list?
      if(randomDouble()<qProb && nExtra<NModq) {
        c.add(*(a[nCount]));  
        nCount++;  nExtra++; 
      } // end if extra
      // now add current list to the overall list list
      add(c);
      
      i++;
    } // end while i
    // because of probabilistic inclusion, we now have to ensure that we have 
    // included all 'extra' nodes.  Here we place 'left-over' nodes randomly
    // among the q lists whereas above, we were trying to maintain 'even'
    // lists by adding at most a single extra node to each list.
    int wList;
    while(nExtra<NModq) {
      wList = randomInt(0,q-1);            // choose random list
      lists[wList]->add(*(a[nCount]));   // add nodes to random list
      nCount++;  nExtra++;
    } // end while i
    
    // consistency checks before exit
    if(nCount!=N)  errorMsg("Did not correctly add nodes in initEven!");
    else if(nExtra!=NModq)  
           errorMsg("Did not correctly add extra nodes in initEven!");
  } // end else (q does not divide N)

  for(int i=0; i<b.getSize(); i++)  delete a[i];
  delete[] a; // return local list array
  return;
}; // end initEven given integer range


template <typename T>
void  TListArray<T>::initEven(Array1D<T> &v, int &q, bool bRandomOrder) {
  // bRandomOrder will randomize the order of the integers
  // the integers can be general
  int N = v.getSize();  // use given size of v
  //cout << redText << "here in initEven a" << endl; // debugging
  Array1D<T> *pv = &v;  // define pointer in case we randomize
  //cout << redText << "here in initEven b" << endl; // debugging
  if(q>N) {
    q = N;
    warningMsg("q is larger than N in initRandom(...) list list. Resetting q = N!");
  } // end if q
  //cout << redText << "here in initEven c" << endl; // debugging

  // now randomize the node order if indicated
  if(bRandomOrder) {
    Array1D<T> u(N);  // we don't want to disturb the passed vector
    u = v;
    u.randomizeOrder();
    pv = &u;          // reset v pointer for below loop to vector u
  } // end if bRandomOrder
  
  //cout << redText << "here in initEven d" << endl; // debugging
  erase();                        // erase current list list
  //cout << redText << "here in initEven e0" << endl; // debugging
  bool qdivN = !((bool)(N%q));      // does q divide N? logic is easier
  //cout << redText << "here in initEven e1 with qdivN = " << qdivN << endl; // debugging
  if(qdivN) {
    //cout << redText << "here in initEven e2" << endl; // debugging
    int n = N/q;                  // number of nodes per list
    //cout << redText << "here in initEven e3" << endl; // debugging
    int nCount = 0;               // nCount - track total number all nodes added
    //cout << redText << "here in initEven e4" << endl; // debugging
    TList<T> c;                   // declare the temporary list
    //cout << redText << "here in initEven e5" << endl; // debugging
    for(int i=0; i<q; i++) {
      //cout << redText << "here in initEven f0" << endl; // debugging
      c.erase();                  // erase current list
      // add nodes to current list
      //cout << redText << "here in initEven f1" << endl; // debugging
      for(int k=0; k<n; k++) { c.add((*pv)[nCount]);  nCount++; }
      //cout << redText << "here in initEven f2" << endl; // debugging
      add(c);                                // add full list to list list
      //cout << redText << "here in initEven f3" << endl; // debugging
    } // end for i
    // end if qdivN
    //cout << redText << "here in initEven g" << endl; // debugging
    
    // consistency check before exit
    if(nCount!=N)  errorMsg("Did not correctly add even nodes in initEven!");
  } else {
    //cout << redText << "here in initEven A" << endl; // debugging
    // q does not divide N so initialize as evenly as possible
    int NModq = N%q, nMin = N/q;  // minimum # of nodes per list
    //cout << redText << "here in initEven B" << endl; // debugging
    TList<T> c;                
    //cout << redText << "here in initEven C" << endl; // debugging
    // need some variables to try to evenly distribute the 'extra' nodes evenly
    // among the lists
    // qProb  - probability of adding an 'extra' spin to a given list
    // nExtra - track total number of 'extra' nodes added during process, 
    //          this should be (N mod q) by the end of the initialization
    // nCount - track total number all nodes added during process
    double qProb = (double)(N%q)/(double)(N/q);  
    //cout << redText  << "here in initEven D0 with qProb = " << qProb 
    //     << normText << endl; // debugging
    int    nCount = 0, nExtra = 0, i = 0;
    while(i<q) {
      //cout << redText << "here in initEven D1 with i = " << i 
      //     << normText << endl; // debugging
      c.erase();                             // erase current temporary list
      //cout << redText << "here in initEven D2 with i = " << i 
      //     << normText << endl; // debugging
      // add minimum number of nodes to current list
      for(int k=0; k<nMin; k++) { c.add((*pv)[nCount]);  nCount++; }
      //cout << redText  << "here in initEven D3 with nCount = " << nCount 
      //     << normText << endl; // debugging
      // for uneven lists, do we add another single node to this list?
      if(randomDouble()<qProb && nExtra<NModq) {
        c.add((*pv)[nCount]);  
        nCount++;  nExtra++; 
      } // end if extra
      //cout << redText << "here in initEven D4 with nExtra = " << nExtra 
      //     << normText << endl; // debugging
      // now add current list to the overall list list
      //c.display();  // debugging
      add(c);
      //cout << redText << "here in initEven D5 with i = " << i 
      //     << normText << endl; // debugging
      
      i++;
    } // end while i
    //cout << redText << "here in initEven E with qProb = " << qProb << endl; // debugging
    // because of probabilistic inclusion, we now have to ensure that we have 
    // included all 'extra' nodes.  Here we place 'left-over' nodes randomly
    // among the q lists whereas above, we were trying to maintain 'even'
    // lists by adding at most a single extra node to each list.
    int wList;
    while(nExtra<NModq) {
      wList = randomInt(0,q-1);             // choose random list
      lists[wList]->add((*pv)[nCount]);  // add nodes to random list
      nCount++;  nExtra++;
    } // end while i
    cout << redText << "here in initEven F with qProb = " << qProb << endl; // debugging

    // consistency checks before exit
    if(nCount!=N)  errorMsg("Did not correctly add nodes in initEven!");
    else if(nExtra!=NModq)  
      errorMsg("Did not correctly add extra nodes in initEven!");
  } // end else (q does not divide N)

  return;
}; // end initEven given integer vector


template <typename T>
void  TListArray<T>::initRandom(int N, int &q, int iBegin) {
  // bRandomOrder will randomize the order of the integers
  //cout << "here a in initRandom" << endl;  // debugging
  Array1D<T> v(N);                 // initialize a working vector
  //cout << "here b in initRandom" << endl;  // debugging
  v.initStep(iBegin,1);            // initialize nodes in consecutive order
  //cout << "here c in initRandom" << endl;  // debugging
  // random is 0 below because we already randomized the vector if desired
  initRandom(v,q);                 // call basic function
  //cout << "here d in initRandom" << endl;  // debugging
  return;
}; // end initRandom list list


template <typename T>
void  TListArray<T>::initRandom(Array1D<T> &v, int &q) {
  // initialize by randomly choosing lists to add nodes
  int N = v.getSize();              // use given size of v
  int wList;                     // add node to which list?
  if(q>N) {
    q = N;
    warningMsg("q is larger than N in initRandom(...) list list.  Resetting q = N!");
  } // end if q
  //cout << "here a in initRandom with N = " << N << " and q = " << q << endl;  // debugging

  // randomize the node order - partially redundant with q..N part loop below
  Array1D<T> u(N);     // we don't want to disturb the passed vector
  //cout << "here b in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  u = v; 
  //cout << "here c in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  u.randomizeOrder();

  cout << u << endl;  // debugging

  //cout << "here d in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  erase();             // get rid of current list list
  //cout << "here e in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  // since q is specified by constraint, give each list at least one node
  for(int i=0; i<q; i++) {
    //cout << "here e1 in initRandom with N = " << N << " and q = " << q << endl;  // debugging
    lists[i] = new TList<T>;     // declare a new valid list
    //cout << "here e2 in initRandom with N = " << N << " and q = " << q << endl;  // debugging
    lists[i]->add(u[i]);         // add one random node to the new list
    //cout << "here e3 in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  } // end for i - declaration and one node addition
  //cout << "here f in initRandom with N = " << N << " and q = " << q << endl;  // debugging
  // now add rest of nodes randomly
  for(int i=q; i<N; i++) {
    wList = randomInt(0,q-1);   // choose random list
    lists[wList]->add(u[i]);  // add nodes to random list
  } // end for i - random fill
  cout << "here g in initRandom with N = " << N << " and q = " << q << endl;  // debugging

  return;
}; // end initRandom list list


template <typename T>
int  TListArray<T>::add(TList<T> &c) {
  // add a list to this list
  //cout << redText << "here in add(list) with list size = " << c.getSize() 
  //     << " and nLists of " << nLists << endl; // debugging
  if(nLists==NMaxLists)
    errorMsg("Cannot add another list to the TListArray!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero list to the TListArray.  Use addEmpty()");
  else {
    //cout << redText << "here in initEven b0 " << endl; // debugging
    //lists[nLists] = new TList<T>;  // eliminated non-empty TLists
    //cout << redText << "here in initEven b1 " << endl; // debugging
    *(lists[nLists]) = c;
    //cout << redText << "here in initEven b2 " << endl; // debugging
    nLists++;
    //cout << redText << "here in initEven b3 " << endl; // debugging
  } // end else

  return nLists;
}; // add a single list to the end of the list


template <typename T>
int  TListArray<T>::add(TListArray<T> &c) {
  // add a whole list at a time by copying data
  if(nLists>NMaxLists-c.getSize())
    errorMsg("Appended list list is too large to add to this list!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero list list to the list.");
  else  for(int i=0; i<c.getSize(); i++)  add(*(c.lists[i]));

  return nLists;
}; // add a single list to the end of the list


template <typename T>
int TListArray<T>::move(TListArray<T> &c) {
  // move a whole list at a time by pointer assignment
  if(nLists>NMaxLists-c.getSize())
    errorMsg("Appended list is too large to add to this list array!");
  else if(c.getSize()<1)
    warningMsg("Attempting to add a size zero list to the list array.");
  else { 
    for(int i=0; i<c.getSize(); i++) {
      c.lists[nLists+i] = lists[i];  
      //lists[i] = NULL;  // no longer removing empty TLists
      lists[i] = new TList<T>;  // create a new empty TList for this
    } // end for i
    // finally update parameters of moved list
    c.nLists += nLists;  c.nNodes += nNodes;
    nLists = 0;  nNodes = 0;
  } // end else

  return c.nLists;
}; // move list list


template <typename T>
void  TListArray<T>::erase(bool bKeepEmpty) {
  // delete all lists but don't return list pointer array (just lists)
  //for(int i=0; i<nLists; i++) { lists[i]->erase();  lists[i] = NULL; }
  // redo to assume that we maintain an set of empty TLists
  //cout << "nLists = " << nLists << " in TListArray::erase()" << endl; // debugging
  for(int i=0; i<nLists; i++)  lists[i]->erase();
  if(!bKeepEmpty)  nLists = 0;
  nNodes = 0;
  return;
}; // delete all lists


template <typename T>
int  TListArray<T>::erase(int i, bool bFast) {
  // for bFast speed, the order of lists is *not* preserved
  //cout << redText << "Deleting list " << i << "!" << normText << endl;
  // redo assuming that we maintain an set if empty TLists (not NULL pointers)
  //cout << magenta << "k-1 " << (int)lists[i-1] << " " << flush
  //                << "k   "   << (int)lists[i]   << " " << flush
  //                << "k+1 " << (int)lists[i+1] << " " << endl; // debugging

  if(bFast) {
    if(i<nLists-1)  swap(i,nLists-1);  // swap i with last list
    if(lists[nLists-1]->getSize()>0)  lists[nLists-1]->erase();
  } else {
    //if(lists[nLists-1]==NULL)  errorMsg("Yea");  // debugging
    //if(lists[i]==NULL)  errorMsg("Nay"); // debugging
    if(i<nLists-1) { // copy TLists down preserving the order
      TList<T> *pn = lists[i];  // store to 'move' to end of list as empty
      if(pn->getSize()>0)  pn->erase();
      // now copy the pointers down to preserve the order
      for(int j=i; j<nLists-1; j++)  lists[j] = lists[j+1];
      lists[nLists-1] = pn;
    } // end if i
    else  // else just erase last list
      if(lists[nLists-1]->getSize()>0)  lists[nLists-1]->erase();
  } // end else bFast
  nLists--;

  return nLists;
}; // erase list i


template <typename T>
void  TListArray<T>::swap(int i, int j) {
  // perform a fast pointer swap to change the positions of lists i and j
  TList<T>* pc = lists[i];  lists[i] = lists[j];  lists[j] = pc;
  return;
}; // swap lists i and j


template <typename T>
TListArray<T>& TListArray<T>::operator=(const TListArray<T>& b) {
  // assumes same size at the moment!pij
  if(NMaxLists<b.NMaxLists)  
    errorMsg("TListArrays are different sizes for '=' operator!");

  TList<T> *pc;
  nLists = b.nLists;
  for(int i=0; i<b.nLists; i++) {
    pc = lists[i];
    //cout << " a = with i = " << i << endl;
    //if(lists[i]!=NULL && )  lists[i]->erase();  // no longer allow non-empty
    if(pc->getSize()>0)  pc->erase();
    //cout << " c = with i = " << i << endl;
    //lists[i] = new TList<T>; // no longer allow non-empty TLists
    //cout << " d = with i = " << i << endl;
    (*pc) = *(b.lists[i]);  // copy list by list
    //cout << " e = with i = " << i << endl;
  } // end for i
  // now finish with specific TListArray info
  nNodes = b.nNodes;

  return *this;
}; // end operator=


template <typename T>
void TListArray<T>::sort() {
  // sorts lists by the first node in the list.  This will, of course,   
  // be a problem when the first node is the missing or extra node.
  // Uses an insertion sort
  int     i, j;
  TList<T> *pTemp;  // pointer to temp list

  errorMsg("Bug in TListArray<T>::sort().  Temp. sort by individual TLists[i]");

  //cout << "Here a in sort lists... " << endl;  // debugging
  for(i=0; i<nLists; i++) { 
    //lists[i]->display();  // debugging
    //cout << "Here a1 in sort lists with i = " << i << " and list size of "
    //     << lists[i]->getSize();  // debugging
    if(lists[i]->getSize()>1)  lists[i]->sort();
    //cout << " after... " << endl;  // debugging
    //lists[i]->display();  // debugging
  } // end for i

  //cout << "Here b in sort lists... " << endl;  // debugging

  for(i=1; i<nLists; i++) {
    j = i;                 // j and i counters start at same number index
    pTemp = lists[i];   // store next item in list to drop in correct place
    while(j>0 && 
          //lists[j-1]->getStart().getListNode()>pTemp->getStart().getListNode() ) {
          (lists[j-1]->getSize()>pTemp->getSize() // sort by size
          // and sub-sort by particle number if size is the same
           || (lists[j-1]->getSize()==pTemp->getSize() && 
               //lists[j-1]->getStart().getListNode()>pTemp->getStart().getListNode()))
               lists[j-1]->getStart()>pTemp->getStart()))
         ) {
    ///while( j>0 && lists[j-1]->getSize()>pTemp->getSize() ) {
      lists[j] = lists[j-1]; // move list pointer up to next spot 
      j--;
    } // end while j
    lists[j] = pTemp;   // move list pointer up to next spot
  } // end for i
 
  return;
} // end sortLists


template <typename T>
inline int TListArray<T>::getnNodes(bool bFast) const { 
  if(bFast)  return nNodes;
  else {
    nNodes = 0;
    for(int i=1; i<nLists; i++)  nNodes += lists[i].getSize();
  } // end else
  return nNodes;              
}; // end getnNodes


// --- related functions -------------------------------------------------------
template <typename T>
inline ostream& operator<<(ostream& fout, TListArray<T>& a) {
  // Outputs the List utilizing the display function
  a.display();
  return fout;
};



// ---------------------------------------------------------------------------
} // end namespace MSL
// ---------------------------------------------------------------------------

#endif
