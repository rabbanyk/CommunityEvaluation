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

#ifndef MSL_TCLUSTERLIST_H
#define MSL_TCLUSTERLIST_H

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
#ifndef MSL_TCMSPARSEW_H
// the edges and data data structure
#include "./msl/MSL_CMatrix_SparseW.h"
#endif

using namespace std;
using namespace MSL;

namespace MSL {

class TCluster;
class TClusterList;
class TNodeSphere;
typedef unsigned short TBool;  // define my own boolean type

class TNode {
// Memory allocation is dynamic and uses standard C++ array indexing with 
// brackets as in arrayName[index] with zero indexing.
 //friend class TCluster;      // not required now
 //friend class TClusterNode;  // don't need since basic data are public anyhow
 
 public:
  // Public member functions
  // overloaded operators
           //TNode&   operator=(const int    s);
           TNode&     operator=(const TNode& s);
           void       init(int s, TBool c, TBool f);  // initialize w/raw data
           //void     display();                      // display data to console
  // data functions
  // unsetSuperNode does *not* expand supernode since that requires a cluster
  //inline void       setSuperNode(TCluster &c) { pCl = &c;    return;  }; 
  //inline void       unsetSuperNode()          { pCl = NULL;  return;  }; 
  //inline TBool      isSuperNode()       const { return (pCl != NULL); };
  //inline TCluster&  superNodeCluster()  const { return *pCl;          };
  inline   int        getNode()           const { return node;          };
  //inline TBool      isClean()           const { return clean;         };
  //inline TBool      isFlagged()         const { return flag;          };
  //inline void       setClean(TBool c=1)       { clean = c;  return;   };
  // insert placeholders for physical system generation functions
  inline   bool       isOverlap(double x, double y, double r);
  inline   bool       isOverlap(TNodeSphere &s);

  // standard public member functions
  TNode(int node = 0);                          // bare constructor
  ~TNode();                                     // destructor
  
  // data elements    // basic data elements are public
  int        node;    // current number clusters stored
  TBool      clean;   // has the node been modified or moved?
  TBool      flag;    // is the node flagged by user?
  //TBool    fixed;   // flag to fix this node as unmovable (if needed)
  
 protected:           // protected data elements
  //TCluster*  pCl;   // a pointer to an *external* cluster whose memory
                      //  is *not* returned on destruct - just hierarchical
}; // end class TNode;

inline bool TNode::isOverlap(double x, double y, double r) {
  errorMsg("isOverlap() is not implemented for TCluster (use TNodeSphere...)");
  return 0;
};  // end isOverlap for TNode
inline bool TNode::isOverlap(TNodeSphere &s) {
  errorMsg("isOverlap() is not implemented for TCluster (use TNodeSphere...)");
  return 0;
};  // end isOverlap for TNode


class TClusterNode {
 friend class TCluster;
 
 public:
  // overloaded operators
           //TNode&  operator=(const int    s);
           //TNode&  operator=(const TNode& s);
  // Public member functions
  //           void    init(int s, TBool c=0, TBool f=0, // initialize the node
  //                        TCluster *pc = NULL, TClusterNode *pn = NULL);  
             void    init(int s, TBool c=0, TBool f=0, // initialize the node
                          TClusterNode *pn = NULL);  
//           void       display();          // display this cluster to console
  //inline   TBool    isSuperNode() const { return n.isSuperNode();      };
  inline   void    setClean()          { n.clean = 1;  return;        };
  // unsetSuperNode does *not* expand supernode since that requires a cluster
  //inline   void    unsetSuperNode()    { n.unsetSuperNode();  return; }; 
  // Standard public member functions
  TClusterNode(int node = 0);                          // bare constructor
  ~TClusterNode();                                     // destructor
  // everything is public here
  // Data elements
 protected:
  TNode         n;      // current number clusters stored
  TClusterNode *pNext;  // maximum number stored clusters allowed
}; // end class TClusterNode;

//------------------------------------------------------------------------------
class TCluster {
// Memory allocation is dynamic and can use standard C++ array indexing with
// brackets as in arrayName[index] with zero indexing.
// friend istream&   operator>>(istream& fin, TCluster& a);

 public:
 // overloaded operators
           TCluster&  operator=(const TCluster& b);
           //TBool       operator==(const TCluster& b);
  // C array reference for a linked-list, not efficient!
           TNode&     operator[](int i);
  // Public member functions
           // initialize a cluster given a cluster or integer list
           void       init(int nodeList[], int length);  // reinitialize
           void       init(TVectorInt &v, TBool bRandom=1);
           void       initClean();                  // reset all nodes to clean
           void       clearFlags();                 // reset all flags to false
  inline   void       clearFlag() { bFlag = 0;  return; } // reset cluster flag to false
  // add a node to the end, only the dE/dQ version updates properties
  // add node by detailed parameters, node changes return number of nodes left
  //       int        add(int node, TBool bClean=0, TBool bFlag=0, TCluster *pCl = NULL);
  //inline int        add(int node, int dE, TBool bClean, TBool bFlag, TCluster *pCl = NULL);
           int        add(int node, TBool bClean=0, TBool bFlag=0);
  inline   int        add(int node, int dE, TBool bClean, TBool bFlag);
           int        add(TNode    &node);      // add node to this cluster
           int        add(TCluster &cluster);   // add cluster to this cluster
  // add and update with a known energy and/or modularity change
  inline   int        add(TNode &n, int dE, TFloat dQ=0.0);
  // 'move' functions only manipulate pointers; they do not allocate memory
           int        moveCurrent(TCluster &c); // move current node *to* c
           int        moveStart(TCluster &c);   // move start node *to* c
  inline   int        moveCurrent(TCluster &c, int deltaEAdded, int deltaERemoved);
           void       move(TCluster &c);        // move this cluster *to* c
  inline   TNode&     setCurrent(int wNode);    // set manual iteration to wNode
                      // move this cluster *to* c with a known energy change
  inline   void       move(TCluster &c, int deltaE);
  // erase function delete the object and return the memory
           int        eraseCurrent();  // erase current node and return memory
           void       erase();         // erase all data and return memory
  // energy calculation functions - get... functions do not make energy or 
  // modularity changes.  They simply return a 'what-if' value.
           // return node integer value for pCurrent - for when only the integer
           // node value is needed
  inline   TNode&     getStart()       const { return pBegin->n;             };
  inline   int        getStartNode()   const { return pBegin->n.node;        };
  inline   int        getNode(int i)         { return (*this)[i].node;       };
  inline   bool       getFlag(int i)         { return (bool)(*this)[i].flag; };
  inline   bool       getFlag()              { return (bool)bFlag;           };
  inline   void       setFlag(int i, TBool b){ (*this)[i].flag = b; return;  };
  inline   void       setFlag(TBool b = 1)   { bFlag = b;  return;           };
           // return node by reference for pCurren t
  inline   TNode&     getCurrent()     const { return pCurrent->n;           };
  inline   int        getCurrentNode() const { return pCurrent->n.node;      };
  inline   int        getNextNode()    const; // with no step iteration
  inline   int        getEnergyChangeAdded  (TNode &n,  TCMatrix &Cp) const;
  inline   int        getEnergyChangeAdded  (TNode &n,  TCMatrix &Cp, int &nConnected) const;
  inline   int        getEnergyChangeRemoved(TNode &n,  TCMatrix &Cp) const;
  inline   int        getMergeEnergy(TCluster &b, TCMatrix &Cp) const;
  inline   int        getMergeEnergy(TCluster &b, TCMatrix &Cp, int &nCon)const;
  inline   int        getMergeEnergies(TCluster &b, TCMatrix &Cp, 
                                       int &eCon, int &eUncon, int &nCon) const;
  //inline int        getkSum() const { return kSum;      };
           int        calckSum(TCMatrix &Cp) const;
           int        calcEdges(TCMatrix &Cp, int iNode) const;
  // calculate energy from scratch - required at least once for some data
           int        calcEnergy(TCMatrix &Cp, TBool bSymm=1);
           void       clearEnergy();
  // other functions on the cluster
           void       sortNodes();   // performs an insertion sort on the nodes
           void       display(int offset=0) const; // display cluster to console
           void       copyToVector(TVectorInt &v);
  // iteration loop functions simulate array-like iteration through the list
  // begin() is required which starts a pCurrent pointer which will step through 
  // the list as next() is called.  It is only modified by user calls except for
  // some consistency updates.  There is only one iteration loop at a time, but
  // interior functions do not directly utilize the variables.
  inline   TNode&     begin();       // start an iteration
  inline   TNode&     next();        // step to next node in an iteration
  // data acquisition functions
  inline   int        getSize()        const { return nNodes;              };
  inline   int        getn()           const { return nNodes;              };
  inline   double     getnD()          const { return (double)nNodes;      };
  // these data aquisition functions may require a calcEnergy(...) call
  inline   int        getEnergy()      const { return energy;              };
  inline   int        getE()           const { return energy;              };
  inline   int        getCEdgeSum()    const { return cEdgeSum;            }; 
  inline   int        getUEdgeSum()    const { return uEdgeSum;            }; 
  inline   TBool      isClique()       const; 
  inline   int        getNEdges()      const { return nEdges;              }; 
  inline   unsigned   getl()           const { return nEdges; };
  inline   void       setNEdges(int e)       { nEdges = e;  return;        }; 
  inline   int        getMaxEdges()    const { return nNodes*(nNodes-1)/2; };
  inline   double     getPartialQ()    const { return partialQ;            };
  inline   double     getPartialQW()   const { return partialQW;           };
  inline   double     getDensity()     const;

  // Standard public member functions
  TCluster();                          // bare constructor
  TCluster(int nodeList[], int size);  // constructor with integer array
  ~TCluster();                         // destructor
  
  // utility public member functions
           bool       isValid(string callingFunction) const;
 // temporarily disabled protected status of variables
 //protected:               // data elements are protected
  int            nNodes;    // current number clusters stored
  //int          MaxNodes;  // maximum number stored clusters allowed
  int            nEdges;    // number of edges in the cluster from calcEnergy()
  double         partialQ;  // partial contribution of modularity for cluster
  double         partialQW; // partial contribution of weighted modularity
  int            energy;    // energy of each cluster
  int            cEdgeSum;  // sum of connected edge weights in cluster
  int            uEdgeSum;  // abs sum of unconnected edge weights in cluster
  //int          kSum;      // stores the degree sum in the cluster
  TClusterNode  *pBegin;    // pointer to beginning of data - NULL if no data
  TClusterNode  *pEnd;      // pointer to end of data list - NULL if no data
  TClusterNode  *pCurrent;  // pointer to the current node in an iteration
  TClusterNode  *pPrevious; // pointer to the previous node in an iteration
  TClusterNode  *piCurrent; // pointer to the current node in op[] iteration
  int            iCurrent;  // an internal index for operator[] functions which
                            //   provides simulated array referencing as an 
                            //   alternative to the pointer iterators above
  TBool          bFlag;     // a generic flag for use by the application
}; // end class TCluster;


inline int TCluster::getNextNode() const { 
  if(pCurrent->pNext!=NULL)  return pCurrent->pNext->n.node;
  else { // invalid attempt to retrieve next node number
    if(nNodes==0)
      errorMsg("Attempting to retrieve next node on a size zero cluster.");
    else 
      errorMsg("Attempting to retrieve next node at the end of the cluster.");
  } // end else
  return -1;  // internal invalid value
};  // end getNextNode

inline double TCluster::getDensity() const { 
  if(nNodes>1)
        return (double)nEdges*2.0/( (double)nNodes*(double)(nNodes-1) ); 
  else  return 0.0;
};  // end getDensity

inline TNode& TCluster::begin() {
  // starts an iteration loop, or alternatively just returns the starting node
  // The iteration loop is an 'array' equivalent mechanism to provide current 
  // node in a loop.
  // The iteration loops are for user usage and are not incremented within the classes.
  // They are only updated on user changes in appropriate function.
  pPrevious = NULL;
  pCurrent  = pBegin;
  //if(nNodes<1 || pBegin==NULL)  
  //  warningMsg("Attempting to begin iteration on a size zero cluster.");
  return pCurrent->n;
};  // add with no energy updates

inline TBool TCluster::isClique() const { 
  return nEdges==( nNodes*(nNodes-1)/2 );              
}; // end isClique

// --- related functions -------------------------------------------------------
inline ostream& operator<<(ostream& fout, TCluster& a) {
  // Outputs the Cluster utilizing the display function
  a.display();
  return fout;
};

// --- inline cluster functions ------------------------------------------------
inline int TCluster::getEnergyChangeAdded(TNode &s, TCMatrix &Cp) const {
  // calculate energy change *if* node s were to be added to another cluster
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  TClusterNode *pn = pBegin;
  int eSum = 0;
  for(int j=0; j<nNodes; j++) {
    eSum -= Cp(s.node,pn->n.node);
    pn    = pn->pNext;  // increment to next node in list - end of list is NULL
  } // end for j
  #ifdef DEBUG_MODE
  if(pn!=NULL)
    errorMsg("Pointer iteration in getdEAdded exits with a non-NULL value!?");
  #endif
  return eSum;
}; // end energy change with test node

inline int TCluster::getEnergyChangeAdded(TNode &s, TCMatrix &Cp, 
                                          int &nConnected) const {
  // calculate energy change *if* node s were to be added to another cluster
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  // This version also calculates the number of edges that s is connected to
  // inside this cluster.
  TClusterNode *pn = pBegin;
  int eSum = 0, eVal;
  nConnected = 0;
  for(int j=0; j<nNodes; j++) {
    eVal  = Cp(s.node,pn->n.node);
    eSum -= eVal;
    if(eVal>0)  nConnected++;
    pn    = pn->pNext;  // increment to next node in list - end of list is NULL
  } // end for j
  #ifdef DEBUG_MODE
  if(pn!=NULL)
    errorMsg("Pointer iteration in dE added exits with a non-NULL value!?");
  #endif
  return eSum;
}; // end energy change with test node

inline int TCluster::getEnergyChangeRemoved(TNode &s, TCMatrix &Cp) const {
  // calculate energy change *if* particle i removed from the cluster.
  // i is the number identity not the index.  It also assumes that node
  // i is *in* the cluster to be valid and this is *not* checked.
  // The i==i diagonal element is supposed to be zero, so this element is *not*
  // skipped since it simply adds zero.  There is a debug mode check however.
  TClusterNode *pn = pBegin;
  int eSum = 0;
  for(int j=0; j<nNodes; j++) {
    eSum += Cp(s.node,pn->n.node);
    pn    = pn->pNext;  // increment to next node in list - end of list is NULL
  } // end for j

  #ifdef DEBUG_MODE
  if(Cp(s.node,s.node)!=0)  
    errorMsg("Cp "+itos(Cp(s.node,s.node))+" diagonal element is not zero!?");
  if(pn!=NULL)
    errorMsg("Pointer iteration in dE removed exits with a non-NULL value!?");
  #endif
  return eSum;
}; // end energy change with test node

inline int TCluster::getMergeEnergy(TCluster &b, TCMatrix &Cp) const {
  // calculate energy change *if* node s were to be added to another cluster
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  // This version also calculates the number of edges that s is connected to
  // inside this cluster.
  TClusterNode *pa = pBegin, *pb;
  //int iNode;
  int eSum = 0, eVal;
  for(int i=0; i<nNodes; i++) {
    pb = b.pBegin;  // restart b iteration
    //iNode = pa->n.node;
    
    for(int j=0; j<b.nNodes; j++) {
      eVal  = Cp(pa->n.node,pb->n.node);
      eSum -= eVal;
      pb = pb->pNext;  // increment to next node in list - end of list is NULL
    } // end for j
    
    pa = pa->pNext;  // increment to next node in list - end of list is NULL
  } // end for i
  #ifdef DEBUG_MODE
  if(pa!=NULL || pb!=NULL)
    errorMsg("Pointer in getMergeEnergy exits with a non-NULL value!?");
  #endif
  return eSum;
}; // end energy change with test node

inline int TCluster::getMergeEnergy(TCluster &b, TCMatrix &Cp, 
                                    int &nConnected) const {
  // calculate energy change *if* node s were to be added to another cluster
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  // This version also calculates the number of edges that s is connected to
  // inside this cluster.
  TClusterNode *pa = pBegin, *pb;
  //int iNode;
  int eSum = 0, eVal;
  nConnected = 0;
  for(int i=0; i<nNodes; i++) {
    pb = b.pBegin;  // restart b iteration
    //iNode = pa->n.node;
    
    for(int j=0; j<b.nNodes; j++) {
      eVal  = Cp(pa->n.node,pb->n.node);
      eSum -= eVal;
      if(eVal>0)  nConnected++;
      pb = pb->pNext;  // increment to next node in list - end of list is NULL
    } // end for j
    
    pa = pa->pNext;  // increment to next node in list - end of list is NULL
  } // end for i
  #ifdef DEBUG_MODE
  if(pa!=NULL || pb!=NULL)
    errorMsg("Pointer in getMergeEnergy exits with a non-NULL value!?");
  #endif
  return eSum;
}; // end energy change with test node

inline int TCluster::getMergeEnergies(TCluster &b, TCMatrix &Cp, 
             int &eConnected, int &eUnconnected, int &nConnected) const {
  // calculate energy change *if* node s were to be added to another cluster
  // Note, i is the number identity not the index.  If particle i happens to 
  // be in the test cluster, the i==i diagonal element is supposed to be zero.
  // This version also calculates the number of edges that s is connected to
  // inside this cluster.
  TClusterNode *pa = pBegin, *pb;
  int eVal;
  eConnected = 0;  eUnconnected = 0;  nConnected = 0;
  for(int i=0; i<nNodes; i++) {
    pb = b.pBegin;  // restart b iteration
    //iNode = pa->n.node;
    
    for(int j=0; j<b.nNodes; j++) {
      eVal = Cp(pa->n.node,pb->n.node);
      if(eVal>0) { eConnected   -= eVal;  nConnected++; }
      else       { eUnconnected -= eVal;                }
      pb = pb->pNext;  // increment to next node in list - end of list is NULL
    } // end for j
    
    pa = pa->pNext;  // increment to next node in list - end of list is NULL
  } // end for i
  #ifdef DEBUG_MODE
  if(pa!=NULL || pb!=NULL)
    errorMsg("Pointer in getMergeEnergy exits with a non-NULL value!?");
  #endif
  return (eConnected+eUnconnected);
}; // end energy change with test node

inline TNode& TCluster::next() {
  // starts an iteration loop, or alternatively just returns the starting node
  if(pCurrent!=NULL) {
    // this also handles the case of a moved node where pPrevious is set equal to pCurrent
    pPrevious = pCurrent;    // simply redundant if a node has been moved
    pCurrent  = pCurrent->pNext;
  } // end if
  // if pCurrent is null and cluster size > 0, reset to beginning
  else if(nNodes>0) { pPrevious = NULL;  pCurrent = pBegin; } 
  //else if(pCurrent==NULL) 
  //warningMsg("Iteration error in TCluster next().  pCurrent = NULL, but n=0");
  //     errorMsg("Attempting to iterate with pCurrent on size zero cluster");
  //else  warningMsg("Attempting to iterate with pCurrent past end of list.");
  return pCurrent->n;
};  // add with no energy updates


inline int TCluster::add(TNode &s, int dE, TFloat dQ) {
  // add a node with a known change in energy and modularity
  // this updates the cluster properties for the passed parameters
  add(s);
  // now update other parameters according to this call
  energy   += dE;
  partialQ += dQ;
  return nNodes;
};  // add with known energy change


inline int TCluster::add(int node, int dE, TBool bClean, TBool bFlag) {
  //add(node,bClean,bFlag,pCl);  // temporarily removed pCl parameter
  add(node,bClean,bFlag);
  energy += dE;  // increment the energy
};

  
inline  int  TCluster::moveCurrent(TCluster &c, int deltaEAdd, int deltaERem) {
  // move current node *to* c given an energy change
  moveCurrent(c);
  c.energy += deltaEAdd;
  energy   += deltaERem;
  return nNodes;
};

  
inline  TNode&  TCluster::setCurrent(int wNode) {
  // move current node *to* c given an energy change
  begin();
  while(pCurrent->n.node != wNode)  next();
  return (pCurrent->n);
};

  
inline void TCluster::move(TCluster &c, int deltaE) {
  // move this cluster *to* c and update the energy with a known energy change
  move(c);
  c.energy += deltaE;
}; 

//------------------------------------------------------------------------------
class TClusterList {
// Memory allocation is dynamic and uses standard C++ array indexing with 
// brackets as in arrayName[index] with zero indexing.

 public:
 // overloaded operators
           TClusterList& operator=(const TClusterList& b);
//         TClusterList& operator+=(const TClusterList& b);
//         TClusterList  operator+(const TClusterList& a);
  inline   TCluster&     operator[](int i) const { return *clusters[i]; };
  // Public member functions
           void   initClean();            // reset all cluster nodes to clean
           void   clearFlags();           // reset all cluster node flags to 0
           // initialize a cluster list as evenly as possible given N and q
           void   initEven(int N, int &q, TBool bRandomOrder=1, int iBegin=0);
           // initialize a cluster list given an integer list
           void   initEven(TVectorInt &v, int &q, bool bRandomOrder=1);
           // take a cluster and evenly initialize a cluster list with it
           void   initEven(TCluster &c,   int &q, bool bRandomOrder=1);
           // take the cluster list and evenly re-initialize it
           void   initEven(int q, TBool bRandomOrder=1);
  inline   void   createSymmetric(int N, TBool bRandomOrder=1);
           void   createSymmetric(TVectorInt &v, TBool bRandomOrder=1);
           void   createPower(int N, int nMin, int nMax, double alpha = -1.0, 
                              TBool bRandomOrder = 0);
           void   initPowerS(int N, int nMin, int nMax, double alpha = -1.0, 
                             TBool bRandomOrder = 0);
           void   randomizeOrder();
/*
           // Here is an effort make an all purpose initialization routine
  inline   void   init(TInitClusterList &iType);  // abstracted type - invalid
  inline   void   init(TInitRandom &iType);
  inline   void   init(TInitSymmetric &iType);
  inline   void   init(TInitEven &iType);
  inline   void   init(TInitPower &iType);
*/
           void   initRandom(TVectorInt &v, int &q);
           void   initRandom(int N, int &q, int iBegin=0);
           // initializations that just reorder the existing cluster list 
           void   initEvenS(int q, bool bRandomOrder=1);   // evenly sized
           void   initRandomS(int q); // randomly placed
           void   initSymmetric(bool bRandomOrder=1);      // one node/cluster
           // add and move functions
           int    createEmpty();          // create an empty cluster at the end
  inline   int    addEmpty(); // allow access to an existing empty cluster at the end
  inline   int    eraseEmpty(); // allow access to an existing empty cluster at the end
           int    add(TCluster &c);       // add a cluster to the end by copy
           int    add(TClusterList &c);   // add a clusterList to end by copy
           //int  move(TCluster &c);      // move cluster to end by pointer
           int    move(TClusterList &c);  // move clusterList to end by pointer
           // note that cluster order is not preserved for the move(i,j) call
           // nor the erase(i) call unless bFast is set to 0
           // move cluster i to j and remove i
  inline   void   move(int i, int j, TBool bFast=1);
  inline   void   move(int i, int j, int  deltaE, TBool bFast);
           void   erase();  // delete all clusters, return nodes but leave clusters empty
           int    erase(int i, TBool bFast=1);// delete cluster i, return memory
           void   destroy();          // delete all clusters, return memory
           int    destroy(int i, TBool bFast=1);// delete cluster i, return memory
  inline   void   swapClusters(int i, int j);     // swap clusters i and j
           int    calcEnergy(TCMatrix &Cp, TBool bFast=0);  // 'fast' is NZ calculation
           // in calcQ, r (AFG parameter) is not currently implemented
           double calcQ(TCMatrix &Cp, double gamma=1.0, double r=0.0);
                      //TBool bFastCalc=0);
           void   clearEnergy();
           int    calcNNodes();  // calculate the number of nodes in the list
           // resize with empty lists, currently it DELETES ALL DATA
           void   resize(unsigned newsize, bool bSaveData = 0);   
  // other cluster functions
           // display cluster list to console
           void   display(int offset=0, TBool bSort=0, string s="", int verbosity=1); 
           //void resize(int Length, TBool Fill, T FillValue);
           void   initNodesList();            // initialize the node to cluster list
           void   invalidateNodesList(int i); // invaldate the node to cluster list
                                              //   for cluster i
           void   updateNodesList(int i); // update the node-cluster list on one cluster
           // update the node-cluster list on one cluster with a known change of target
           void   updateNodesList(int i, int target); 
           void   offsetNodes(int offset = 1);  // add 'i' to each node index ID
           void   sortClusters();         // insertion sort on the clusters
           int    clearZeroClusters(TBool bUpdateNodesList=1);
           int    inputLFRANS(string fname, int N);
           int    inputRaw(string fname, int N, int nodeOffset = 0);
           int    outputRaw(string fname, int nodeOffset = 0, bool bSort = 0);
           // calculate average of some parameters
           void   calcZParams (TCMatrix &Cp, double &zIn, double &zOut, double &zAvg);
           void   calcPZParams(TCMatrix &Cp, double &zin, double &zout, double &zavg,
                               double &pin, double &pout, double &p);
  // data acquisition functions
  inline   int    getSize()           const { return nClusters;           };
  inline   int    getq()              const { return nClusters;           };
  inline   int    getqi(int wNode)    const { return nodes[wNode];        };
  // these data functions may require a calcEnergy(...) call depending on what
  // data was passed to the constructor and what changes have been made
  inline   int    getNNodes()         const { return NNodes;              };
  inline   int    getNEdges()         const { return NEdges;              };
  inline   int    getNMaxClusters()   const { return NMaxClusters;        };
  inline   int    getNMaxEdges()      const { return NNodes*(NNodes-1)/2; };
  // these data functions require a calcEnergy(...) call
  inline   int    getEnergy()         const { return energyTotal;         };
  inline   int    getE()              const { return energyTotal;         };
  //inline double getENorm()          const { return eNorm;               };
  // note that getq() above is for the number of clusters (i.e., case sensitive)
  inline   double getQ()              const { return modularity;          };
  inline   double getModularity()     const { return modularity;          };
  // RB modularity is calculated by calcQ(...) and is stored in modularity
  //inline double getQW()             const { return modularityW;         };
  // return average community edge density (unweighted at the moment)
  // size 1 communities are excluded from the average
  inline   double getAvgDensity()     const;
 
  TVectorInt  nodes;         // cluster data arranged by nodes in cluster

  // Standard public member functions
  TClusterList(int MaxClusters=0, int NNodes=0, int NEdges=0);  
  ~TClusterList();
  
 //protected:
  // Data elements
  int       NMaxClusters;
  int       nClusters;       // current number clusters stored
  int       NNodes;          // constant total number of nodes in graph
  int       NEdges;          // constant total number of edges in graph
  //int     NodeOffset;      // constant setting at what number the nodes begin
  // the next values can only be referenced after a call to calcEnergy()
  int       energyTotal;     // sum of energies from all clusters
  int       cEdgeSum;        // sum of connected edge weights in cluster
  int       uEdgeSum;        // abs sum of unconnected edge weights in cluster
  TFloat    modularity;      // modularity of the clustering
  double    eNorm;           // normalized energy on [-1,+1] range - not used
  TFloat    modularityW;     // modularity of the clustering       - not used
  // the clustering data is stored as "2D" dynamic array to prevent a lot of 
  //   data movement during cluster move operations
  TCluster  **clusters;      // actual cluster data
}; // end class TClusterList;

// --- related functions -----------------------------------------------------
inline ostream& operator<<(ostream& fout, TClusterList& a) {
  // Outputs the Cluster utilizing the display function
  a.display();
  return fout;
};

inline double TClusterList::getAvgDensity() const { 
  double  pSum   = 0.0;
  int     cCount = 0;
  for(int i=0; i<nClusters; i++) 
    if(clusters[i]->getSize()>1) { 
      pSum += clusters[i]->getDensity();  cCount++; }
  // now finally return the average density excluding size one clusters
  return pSum/(double)cCount;
};  // end getAvgDensity



// --- inline functions for clusterList --------------------------------------
inline void TClusterList::createSymmetric(int N, TBool bRandomOrder) {
  TVectorInt v(N);                          // order vector
  v.initStep(0,1);                          // initialize order vector
  //if(bRandomOrder)  v.randomizeOrder();
  createSymmetric(v,bRandomOrder);          // init symmetric, no randomize
  //createSymmetric(v,0);          // init symmetric, no randomize
  return;
};


inline void  TClusterList::swapClusters(int i, int j) {
  // perform a fast pointer swap to change the positions of clusters i and j
  TCluster* pc = clusters[i];  clusters[i] = clusters[j];  clusters[j] = pc;
  return;
}; // swap clusters i and j


inline  void  TClusterList::move(int i, int j, TBool bFast) {
  // move cluster i to j and remove i without deleting anything except the 
  // empty cluster i
  // note that the cluster order is not preserved for the erase call
  //cout << "before move... " << flush;  // debugging
  updateNodesList(i,j);  // update with known change (more efficient)
  clusters[i]->move(*(clusters[j]));
  //cout << "done... " << flush;  // debugging
  erase(i,bFast);  // erases the now empty cluster
  //cout << "really done " << endl;  // debugging
  return;
}; // end inline move

inline  void  TClusterList::move(int i, int j, int deltaE, TBool bFast) {
  // move cluster i to j and remove i without deleting anything except the 
  // empty cluster i
  // note that the cluster order is not preserved for the erase call
  updateNodesList(i,j);  // update with known change (more efficient)
  clusters[i]->move(*(clusters[j]),deltaE);
  erase(i,bFast);  // erases the now empty cluster
  return;
}; // end inline move

inline int TClusterList::addEmpty() {
  // add an empty cluster to the end of the list
  #ifdef DEBUG_MODE
  if(nClusters==NMaxClusters)
    errorMsg("Cannot add another cluster to the list!  nClusters = "+
              itos(nClusters));
  #endif
  //cout << "Adding an empty cluster... " << flush;  // debugging
  //clusters[nClusters] = new TCluster;  // new version does not re-allocate
  nClusters++;
  //cout << "done." << endl;  // debugging
  return nClusters;
}; // add a single cluster to the end of the list


inline int TClusterList::eraseEmpty() {
  // add an empty cluster to the end of the list
  #ifdef DEBUG_MODE
  if(nClusters==0)
    errorMsg("Cannot delete a cluster from the list!  nClusters = "+
              itos(nClusters));
  #endif
  //clusters[nClusters] = new TCluster;  // new version does not re-allocate
  nClusters--;
  return nClusters;
}; // add a single cluster to the end of the list

/*
inline void TClusterList::init(TInitClusterList &iType) {
  //int N = iType.N;
  //if(NNodes!=N)  calcNNodes();
  //if(NMaxClusters!=N || NNodes!=N) { destroy();  createSymmetric(N); }
  if(iType.type()==Symmetric)    init(iType.bRandomize);
  else if(iType.type()==Power)   
    initPowerS(N,iType.nMin,iType.nMax,iType.beta,iType.bRandomize);
  else if(iType.type()==Random)  initRandomS(iType.q);
  else if(iType.type()==Even)    initEvenS(iType.q,iType.bRandomize);
  else errorMsg("Invalid init type in TClusterList::init().");
  //errorMsg("TInitClusterList is an abstract type in TClusterList::init().");
  return;
}; // end init 
inline void TClusterList::init(TInitRandom &iType) {
  // cluster list should actually exist before calling this function
  int N = iType.N;
  if(NNodes!=N)  calcNNodes();
  if(NMaxClusters!=N || NNodes!=N) { destroy();  createSymmetric(N); }
  initRandomS(iType.q);
  return;
}; // end init 
inline void TClusterList::init(TInitEven &iType) {
  // cluster list should actually exist before calling this function
  int N = iType.N;
  if(NNodes!=N)  calcNNodes();
  if(NMaxClusters!=N || NNodes!=N) { destroy();  createSymmetric(N); }
  initEvenS(iType.q,iType.bRandomize);
  return;
}; // end init 
inline void TClusterList::init(TInitSymmetric &iType) {
  // cluster list should actually exist before calling this function
  int N = iType.N;
  if(NNodes!=N)  calcNNodes();
  if(NMaxClusters!=N || NNodes!=N) {
    destroy();  createSymmetric(iType.N,iType.bRandomize); 
  } else        initSymmetric(iType.bRandomize);
  return;
}; // end init 
inline void TClusterList::init(TInitPower &iType) {
  // cluster list should actually exist before calling this function
  int N = iType.N;
  if(NNodes!=N)  calcNNodes();
  if(NMaxClusters!=N || NNodes!=N) { destroy();  createSymmetric(N); }
  initPowerS(N,iType.nMin,iType.nMax,iType.beta,iType.bRandomize);
  return;
}; // end init 
*/

//------------------------------------------------------------------------------
class TClusterListArray {
// Memory allocation is dynamic and uses standard C++ array indexing with 
// brackets as in arrayName[index] with zero indexing.

 public:
 // overloaded operators
           //TClusterList& operator=(const TClusterList& b);
  inline   TClusterList& operator[](int i) const { return *lists[i]; };
  // Public member functions
  inline   void   initClean();            // reset all cluster nodes to clean
  inline   void   clearFlags();           // reset all cluster node flags to 0
           // add and move functions
           //int    createEmpty();        // create an empty cluster at the end
  //inline   int    addEmpty(); // allow access to an existing empty cluster at the end
  //inline   int    eraseEmpty(); // allow access to an existing empty cluster at the end
  inline   int    add(TCluster &c);       // add a cluster to the end as a cluster list
  inline   int    add(TClusterList &c);   // add a clusterList to end by copy
  inline   int    addEmpty();             // add a cluster to the end as a cluster list
  inline   int    move(TClusterList &c);  // move clusterList to end by pointer
  // return index of or a pointer to the cluster list with the lowest energy
           int    getMinEIndex() const;         
           int    getMaxQIndex() const;         
   TClusterList*  getMinEList()  const;      
   TClusterList*  getMaxQList()  const;      
           // note that cluster list order is not preserved for the move(i,j) 
           // call, nor the erase(i) call unless bFast is set to 0
           // move cluster i to j and remove i
  inline   void   move(int i, int j, TBool bFast=1);
  inline   void   erase();  // delete all cluster lists, return nodes but leave clusters empty
  inline   int    erase(int i, TBool bFast=1);// delete cluster list i, return memory
  inline   void   destroy();          // delete all cluster lists, return memory
  inline   int    destroy(int i, TBool bFast=1);// delete cluster list i, return memory
  inline   void   swap(int i, int j);     // swap cluster lists i and j
  inline   void   calcEnergy(TCMatrix &Cp, TBool bFastCalc=0);  // 'fast' is NZ calculation
  // other cluster functions
           // display cluster list to console
           //void resize(int Length, TBool Fill, T FillValue);
           //void sort();         // insertion sort on the clusters
  // data acquisition functions
  inline   int    getSize()           const { return nLists;    };
  inline   int    getMaxLists()       const { return NMaxLists; };
 
  // Standard public member functions
  // the listSize = -1 is a flag that tell the constructor to assume that all
  // lists are valid and random access is expected
  TClusterListArray(int MaxLists = 20, int listSize = -1);
  ~TClusterListArray();
  
 //protected:
  // Data elements
  int       nLists;          // the total number of cluster lists stored
  int       NMaxLists;       // the maximum total number of cluster lists

  // the cluster list data is stored as "2D" dynamic array to prevent a lot of 
  //   data movement during move operations
  TClusterList  **lists;         // actual cluster list array
}; // end class TClusterList;

// --- inline functions ------------------------------------------------------
inline void TClusterListArray::initClean() {
  // reset all cluster lists to clean
  for(int i=0; i<nLists; i++)  lists[i]->initClean();
  return; 
}; // end initClean


inline   void   TClusterListArray::clearFlags() {
  // reset all cluster node flags to 0
  for(int i=0; i<nLists; i++)  lists[i]->clearFlags();
  return; 
}; // end clearFlags


inline int TClusterListArray::add(TCluster &c) {
  // add a cluster to the end as an independent cluster list
  if(nLists==NMaxLists)  errorMsg("Cannot add another cluster list");
  TClusterList *pa; 
  pa = new TClusterList;  // add a size one cluster list
  (*pa).add(c);           // use cluster copy routine
  lists[nLists] = pa;
  nLists++;
  return 1; 
}; // end add cluster as list

/*
inline   int    TClusterListArray::addEmpty() {
  // add a cluster to the end as a cluster list
  if(nLists==NMaxLists)  errorMsg("Cannot add another cluster list");
  lists[nLists] = new TClustersList;
  nLists++;
  return 1; 
}; // end add cluster as list
*/

inline int TClusterListArray::add(TClusterList &c) {
  // add clusterList
  if(nLists==NMaxLists)  errorMsg("Cannot add another cluster list");
  (*lists[nLists]) = c;
  nLists++;
  return nLists; 
}; // end add clusterList


inline int TClusterListArray::move(TClusterList &c) {
  errorMsg("TClusterListArray move(c) needs to be fixed.");
  if(nLists==NMaxLists)  errorMsg("Cannont add another cluster list");
  lists[nLists] = &c;
  nLists++;
};  // end move clusterList to end by pointer


inline void TClusterListArray::move(int i, int j, TBool bFast) {
  errorMsg("TClusterListArray move(i,j) needs to be fixed.");
  lists[j] = lists[i];
  lists[i] = NULL;
  return;
}; // end move


inline void TClusterListArray::erase() {
  for(int i=0; i<nLists; i++)  lists[i]->erase();
  nLists = 0;
  return; 
}; // end delete all cluster lists, return nodes but leave clusters empty


inline int TClusterListArray::erase(int i, TBool bFast) {
  if(bFast) {  // more efficient
    if(i<nLists-1)  this->swap(i,nLists-1);  // swap i with last
    // if size...??? CHECK THIS ???  erase cluster but do not return memory
    if(lists[nLists-1]->getq()>0)  lists[nLists-1]->erase();  
  } else {
    TClusterList *pc = lists[i];
    pc->erase();  // erase cluster but do not return memory
    // now copy the pointers down to preserve the order - can easily optimize
    for(int j=i; j<nLists-1; j++)  lists[j] = lists[j+1];
    lists[nLists-1] = pc;         // now 'move' i to last cluster
  } // end else
  // do not invalidate the empty cluster list in new method
  nLists--;

  return nLists;
}; // end delete cluster list i, return memory


inline void TClusterListArray::destroy() {
  // end delete all cluster lists, including empty ones
  for(int i=0; i<NMaxLists; i++)  lists[i]->destroy();
  delete[] lists;
  nLists = -1;  NMaxLists = -1; // invalidate just in case
  return; 
}; // end delete all cluster lists, return memory


inline int TClusterListArray::destroy(int i, TBool bFast) {
  if(bFast) {  // more efficient
    if(i<nLists-1)  swap(i,nLists-1);  // swap i with last
    if(lists[nLists-1]->getSize()>0)  lists[nLists-1]->destroy();  
  } else {
    TClusterList *pc = lists[i];
    pc->destroy();  // erase cluster and return memory
    // now copy the pointers down to preserve the order - can easily optimize
    for(int j=i; j<nLists-1; j++)  lists[j] = lists[j+1];
    lists[nLists-1] = NULL;
  } // end else
  // do not invalidate the empty cluster list in new method
  nLists--;

  return nLists;
}; // end delete cluster list i, return memory


inline void TClusterListArray::swap(int i, int j) {
  if(i<nLists && j<nLists && i>-1 && j>-1) {
    TClusterList *ps;
    ps = lists[i];  lists[i] = lists[j];  lists[j] = ps;
  } else errorMsg("Cluster list array swap is out of range");
  return;
}; // end swap cluster lists i and j


inline void TClusterListArray::calcEnergy(TCMatrix &Cp, TBool bFastCalc) {
  for(int i=0; i<nLists; i++)  lists[i]->calcEnergy(Cp,bFastCalc);
  return; 
};  // end energy calc for array


// --- related functions -----------------------------------------------------
inline ostream& operator<<(ostream& fout, TClusterListArray& a) {
  // Outputs the Cluster utilizing the display function
  for(int i=0; i<a.getSize(); i++) { a[i].display();  fout << "\n"; }
  return fout;
};



// ---------------------------------------------------------------------------
// init classes are an attempt to allow easily used multi-purpose functions
// in main code base
// ---------------------------------------------------------------------------
enum TInitType { Symmetric, Even, Random, Power };
class TInitClusterList {
 protected:
  TInitType  initType;
 public:
  bool       bRandomize;
  int        N;
  TInitType  type() { return initType; };
  virtual    void    init(TClusterList &c) = 0;
  virtual    string  initString(int verbosity=0) = 0;
  TInitClusterList() {};
};
class TInitEven : public TInitClusterList {
 public:
  int        q;
  TInitEven(int NNodes, int nq, bool  bRandom = 1) 
    { initType = Even;  N = NNodes;  q = nq;  bRandomize = bRandom; };
  inline  virtual  void    init(TClusterList &c);
  inline  virtual  string  initString(int verbosity=0);
};

class TInitRandom : public TInitClusterList {
 public:
  int        q;
  TInitRandom(int NNodes, int nq, bool  bRandom = 1) 
    { initType = Random;  N = NNodes;  q = nq;  bRandomize = bRandom; };
  inline  virtual  void    init(TClusterList &c);
  inline  virtual  string  initString(int verbosity=0);
};

class TInitSymmetric : public TInitClusterList {
 public:
  int        q;
  TInitSymmetric(int NNodes, bool  bRandom = 1)
    { initType = Symmetric;  N = NNodes;  bRandomize = bRandom; };
  inline  virtual  void    init(TClusterList &c);
  inline  virtual  string  initString(int verbosity=0);
};

class TInitPower : public TInitClusterList {
 public:
  int        nMin, nMax;
  double     beta;
  TInitPower(int NNodes, int NMin, int NMax, double Beta=-1.0, bool bRandom=1)
    { initType = Power;  N = NNodes;  nMin = NMin;  nMax = NMax;  beta = Beta;
      bRandomize = bRandom; };
  inline  virtual  void    init(TClusterList &c);
  inline  virtual  string  initString(int verbosity=0);
};
inline void TInitRandom::init(TClusterList &c) {
  // cluster list should actually exist before calling this function
  if(c.getNNodes()!=N)  c.calcNNodes();
  if(c.getNMaxClusters()!=N || c.getNNodes()!=N) { 
    if(c.getNMaxClusters()>0)  c.destroy(); 
    c.createSymmetric(N); 
  } // end if
  c.initRandomS(q);
  return;
}; // end init 
inline void TInitEven::init(TClusterList &c) {
  // cluster list should actually exist before calling this function
  if(c.getNNodes()!=N)  c.calcNNodes();
  if(c.getNMaxClusters()!=N || c.getNNodes()!=N) {
    if(c.getNMaxClusters()>0)  c.destroy(); 
    c.createSymmetric(N); 
  } // end if
  c.initEvenS(q,bRandomize);
  return;
}; // end init 
inline void TInitPower::init(TClusterList &c) {
  // cluster list should actually exist before calling this function
  if(c.getNNodes()!=N)  c.calcNNodes();
  if(c.getNMaxClusters()!=N || c.getNNodes()!=N) { 
    if(c.getNMaxClusters()>0)  c.destroy(); 
    c.createSymmetric(N); 
  } // end if
  c.initPowerS(N,nMin,nMax,beta,bRandomize);
  return;
}; // end init 
inline void TInitSymmetric::init(TClusterList &c) {
  // cluster list should actually exist before calling this function
  //cout << "check number of nodes currently at N = " << c.getNNodes() << "... " 
  //     << flush;  // debugging
  if(c.getNNodes()!=N)  c.calcNNodes();
  //cout << "calculated N = " << c.getNNodes() << " with init N = " << N 
  //     << endl << "check max number of nodes currently at " 
  //     << c.getNMaxClusters() << "... " << flush;  // debugging
  if(c.getNMaxClusters()!=N || c.getNNodes()!=N) {
    //cout << "destroy old cluster list? " << flush;  // debugging
    if(c.getNMaxClusters()>0) {
      //cout << "destroying... " << flush;  // debugging
      c.destroy();
    } // end if
    else cout << "did not destroy the current list\n" << flush;  // debugging
    //cout << "create new cluster list with N = " << N << "... " << endl;  // debugging
    c.createSymmetric(N,bRandomize);   // debugging
    //c.display();  // debugging
    //cout << "done with create with NMax = " << c.getNMaxClusters() << endl; // debugging
  } else {
    //cout << "init symmetric... " << endl;  // debugging
    c.initSymmetric(bRandomize);
  } // end else
  return;
}; // end init 
inline string TInitRandom::initString(int verbosity) {
  string s = "random init";
  if(verbosity>0) s += " with q = " + itos(q);
  if(verbosity>1) s += " (n_avg = " + dtos((double)N/(double)q) + ")";
  return s;
}; // end initString 
inline string TInitEven::initString(int verbosity) {
  string s = "even init";
  if(verbosity>0) s += " with q = " + itos(q);
  if(verbosity>1) {
    if((N%q)==0)  s += " (n_avg = " + itos(N/q) + ")";
    else          s += " (n_avg = " + dtos((double)N/(double)q) + ")";
  } // end if
  return s;
}; // end initString 
inline string TInitPower::initString(int verbosity) {
  string s = "power law init";
  if(verbosity>0) s += " with nMin = " + itos(nMin) + " nMax = " + itos(nMax)
                     + " and beta = "  + dtos(beta);
  if(verbosity>1) s += " (n_avg expected = " 
                     + dtos(randomPowerMean(nMin,nMax,beta)) + ")";
  return s;
}; // end initString 
inline string TInitSymmetric::initString(int verbosity) {
  string s = "symmetric init";
  return s;
}; // end initString 
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
} // end namespace MSL
// ---------------------------------------------------------------------------

#endif
