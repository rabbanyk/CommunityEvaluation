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
#include <iostream>
#include <iomanip>
#include <cmath>
#ifndef MSL_ARRAY_H
#include "./msl/MSL_Array.h"
#include "./msl/MSL_Array1D_Template.h"
#endif
//#ifndef _ML_UTILS
//#include "ML_Utils.h"
//#endif
#include "./msl/MSL_Stats1D.h"
*/
using namespace MSL;


/* ***************************************************************************
  community detection functions
*************************************************************************** */
// ---------------------------------------------------------------------------
// main community detection functions
// ---------------------------------------------------------------------------
int communityDetectNZ(TClusterList &clusters, TCMatrix &Cp, 
       TCMData a, TCMData b,   int NStepMax = 1000, bool useZeroE = 0, 
       bool bqConstrained = 0, int verbosity = 0,   int  debugMode = 0) {
  // implements a neighbor-only search of the system.  Note that it uses 
  // optimizations that depend on the specific structure of the connection data 
  // structures.  Care should be used when changing the classes.
  // all int's should be positive (and non-zero) values
  // TryHard is unused but left due to legacy calls

  // -------------------------------------------------------------------------
  // BEGIN MAIN LOOP
  // -------------------------------------------------------------------------
  int     BigInteger = (2<<29);
  int     N = Cp.getSize();  // the number of nodes
  int     iNode, iCluster, jNode, jCluster, minCluster, jDegrees;
  TNode   wNode, *pwCurrent;
  bool    useZeroEnergy = 0;    // do we allow zero energy moves?
  //bool  bqConstrained = 0;    // do we contrain the number of clusters?
  int     lastNodeMoved = -2,   thisNodeMoved = -1; // invalid initial values
  int     lastCluster   = -2,   thisCluster   = -1; // invalid initial values
  //int   lastEnergy    = -1,   thisEnergy    = 0; 
  //int     constantEnergyCount = 0, 
  int     zeroEnergyMoveCount = 0; 
  int     energyChange  = 0, deltaEAdded, deltaERemoved;
  int     equalEnergyCount = 0, nStepMoves;
  unsigned long long int moveCount = 0;
  bool    converged = 0,        moved = 0;
  int    *pNeighborStart,      *pNeighbor; // optimzing neighbor pointers
  int    *pjClustereSum;        // optimzing pointer for energy sum
  int     systemEnergy = 0;     // sum of current system energy

  // initialize data for sparse cases
  TVectorInt  edgeCount(N,0); // count number of nodes connected to a cluster
  TVectorInt  eSum(N,0);      // energy sum for each cluster for current node
  TVectorInt  bcChecked(N,0); // was cluster 'x' checked?
  //clusters.initNodesList(); // initialize the node-in-which-cluster array

  //Cp.display();
  //if(useZeroE)
  //  errorMsg("Zero energy moves are currently disabled in communityDetectNZ!");

  // ************* while nStep ***********************************************
  if(verbosity>0) {
    cout << green << "Starting community detect with ";
    #ifdef MSL_TCMWEIGHTED
    cout << "a " << brown << "weighted";
    #else
    cout << "an " << brown << "unweighted";
    #endif
    cout << green << " neighbor iteration:\n";
  } // end if verbosity

  int nStep = 0;
  while(nStep<NStepMax && !converged) {
     clusters.initClean();  // can optimize this by storing moved nodes
     moved = 0; // indicates whether node moves have occured (not converged)
     nStepMoves = 0; // indicates whether node moves have occured (not converged)

    #ifdef DEBUG_MODE
    if(debugMode>0)                                              // debugging
      cout << magenta << "Begin step loop " << nStep << normText << endl 
           << "nClusters = " << clusters.getSize() << endl;
    #endif

    // loop over particles and move each one in turn
    int i = 0;
    while(i<clusters.getq()) {
      clusters[i].begin();  // begin user iteration
      iCluster = i;         // just to keep clear code

      // ----------------------------------------------------------------
      for(int k=0; k<clusters[i].getn(); k++) {
        //wNode = clusters[i].getCurrent();
        pwCurrent = &( clusters[i].getCurrent() );
        iNode = clusters[i].getCurrentNode();

        // only test node if it hasn't already been moved this step
        //if(wNode.clean) {
        if(pwCurrent->clean) {
          minCluster = i;               // default cluster is its own

          // I think I can optimize the triple neighbor search by storing the
          // node pointer addresses on the first pass and then iterating 
          // directly over the node addresses on the two following lists.
          // Since this is the inner loop, it may have an impact.

          // scan edge connection list for iNode and sum energy contributions
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7.5% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            //jCluster = clusters.nodes[jNode];
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            // count number of nodes jNode is connected to in jCluster
            edgeCount[jCluster] += 1;
            // sum (subtract) the energy for node iNode merge into jCluster for
            // the unweighted version, we can skip the energy sum and just use
            // the count with the stored and unstored values
            // I use pre-processor if's to avoid an extra inner loop if() call
            #ifdef MSL_TCMWEIGHTED
            eSum[jCluster] -= Cp(iNode,jNode); 
            #else
            eSum[jCluster] -= a; 
            #endif

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and update for unconnected edges
          // first, set the energy change for removing iNode from iCluster
          // which is slightly different from the added case
          eSum[iCluster] += (clusters[iCluster].getn()-1-edgeCount[iCluster])
                            *(-b);
          deltaERemoved   = eSum[iCluster];
          bcChecked[iCluster] = 1;
          pNeighbor = pNeighborStart;    // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
              // sum (subtract) the energy for node iNode merge into jCluster
              // note that the unstored value is a negative number
              eSum[jCluster] += (clusters[jCluster].getSize()-edgeCount[jCluster])
                                *(-b);
              bcChecked[jCluster] = 1;
            } // end if bcChecked

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and clear changed values
          deltaEAdded = BigInteger;
          pNeighbor   = pNeighborStart;  // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
            // this if statement catches iCluster!=jCluster by default but has 
            // a lot of extraneous assignments.  Not sure which is better.
            //if((*pjClustereSum)<deltaEAdded) {
            // find the lowest energy change but skip the current cluster
            if(iCluster!=jCluster && (*pjClustereSum)<deltaEAdded) {
              deltaEAdded = (*pjClustereSum);
              minCluster  = jCluster;
            } // end if eSum
            // clear the changed parameters - assignements can be redundant
            (*pjClustereSum)    = 0;
            edgeCount[jCluster] = 0;
            bcChecked[jCluster] = 0;

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search
          energyChange = deltaEAdded - deltaERemoved;

          // -------------------------------------------------------------
          // now we have to check the energy and possible moves
          // -------------------------------------------------------------
          // move if energy is lowered, everything is already set if the energy 
          // is lowered, here we check the special cases
          if(energyChange<0) {
            // we have a candidate move, but we must still check to see if 
            // deltaE>0 and if q is unconstrained since we can lower the energy
            // even more in that case
            if(deltaERemoved>0 && !bqConstrained) {  // open new cluster instead...
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead to lower the energy more
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              clusters.addEmpty();                  // declare an empty cluster
              minCluster = clusters.getq()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Opening an empty cluster at position " 
                   << minCluster << normText << endl;  // debugging
              #endif
            } // end empty cluster case
          } else if(energyChange==0 && useZeroEnergy && clusters[i].getn()>1) {
              // allow zero energy moves as long as we have a  
              // this case is ok - allow the move
              zeroEnergyMoveCount++;
              #ifdef DEBUG_MODE
              cout << green << "Attempting to execute a zero energy move of "
                   << "node " << iNode << " to cluster " << minCluster 
                   << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else if(deltaERemoved>0 && !bqConstrained) {
              // if q is unconstrained, we can still potentially move iNode
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              clusters.addEmpty();               // declare an empty cluster
              minCluster = clusters.getq()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Adding a preferred empty cluster at position "
                   << minCluster << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else {
            #ifdef DEBUG_MODE
            cout << red << "Invalidating move for node " << iNode 
                 << " in cluster " << iCluster << normText << endl; // debugging
            #endif
            minCluster = i;  // otherwise invalidate this move
            //return 0; // debugging
          } // end else

          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          if(minCluster!=i) {
            moveCount++;
            nStepMoves++; // indicates whether node moves have occured (not converged)
            systemEnergy += energyChange;  // track tot al energy incrementally?

            #ifdef DEBUG_MODE
            if(debugMode>0) 
              cout << red << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            //if(debugMode>1)
              cout << "Known energy change for move " << moveCount << " is:  "
                   << energyChange << " = " 
                   << (-deltaERemoved) << " + " << deltaEAdded
                   << " for node " << iNode << " moving from cluster "
                   << i << " (size " << clusters[i].getSize() << ") to cluster " 
                   << minCluster << " (size " << clusters[minCluster].getSize() 
                   << ") on step " << nStep << ".  The total energy is supposedly " 
                   << systemEnergy << endl; // debugging
            #endif

            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            clusters.nodes[iNode] = minCluster;
            // move from the old to new cluster, after following move, the 
            // 'current' iteration pointer is invalid (points to 'previous')
            clusters[i].moveCurrent(clusters[minCluster],deltaEAdded,deltaERemoved);

            //cout << blueText << "The node to cluster list is:  " << redText 
            //     << clusters.nodes << normText << endl;  // debugging

            moved = 1;  // changes were made, so converged is false, unless...
            pwCurrent->clean = 0;  // invalidate further moves of wNode on this step

            // unless we are moving the same node between the same clusters
            if(lastNodeMoved==thisNodeMoved) {
            //if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << blue  << "Moving node " << iNode << " again on iteration " 
                   << nStep << " from cluster " << i << " to " << minCluster
                   << " with an energy of " << systemEnergy
                   << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

            // unless we have reached the zero energyChange move limit
            if(energyChange<0)        equalEnergyCount = 0;
            else if(energyChange==0)  equalEnergyCount++;
            if(equalEnergyCount>1000 ) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << red << "Ended trial due to excessive consecutive zero "
                   << " energy moves " << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

          } // end if minCluster
          #ifdef DEBUG_MODE
          else  cout << grey  << "Rejected a move of node " << iNode 
                     << " to cluster " << minCluster << normText << endl;
          #endif
          // -------------------------------------------------------------
          // -------------------------------------------------------------
        } // end if wNode.clean

        clusters[i].next();  // increment to next node for user iteration
      } // end for k - particle loop
      // -----------------------------------------------------------------

      //if(clusters[i].getq()==0 && !bqConstrained) {  // debugging
      // nStep<3 is an optimization that limits the number of updates of the
      // nodes-cluster list.  Need to clear zero clusters at end of function.
      if(clusters[i].getSize()==0 && nStep<3 && !bqConstrained) {
        #ifdef DEBUG_MODE
        if(debugMode>0)  // debugging
          cout << red << "deleting empty cluster " << i << ":  " << clusters[i]
               << normText << endl;
        #endif
        clusters.erase(i);      // does not return memory
      } // end delete empty cluster
      else  i++;  // if no erase, increment as normal
    } // end while i - cluster loop

    converged = !moved;

    if(converged && useZeroE && !useZeroEnergy) {
      useZeroEnergy = 1; // turn on zero energy moves
      converged = 0;     // unset converged status
      #ifdef DEBUG_MODE
      cout << red << "Switching on zero energy moves at energy " << systemEnergy 
           << normText << "\n";
      #endif
    } // end if

    nStep++;
    if(verbosity>1) 
      cout << blue << nStepMoves << " moves on " << nStep << ", " << normText;
  } // end while nStep - iteration loop

  if(verbosity>1) {
    cout << blue << "Done after " << nStep << " steps "  << magenta << "(Made " 
         << moveCount << " moves)" << normText << endl;
  } // end if

  // eliminate size zero clusters before exiting, order is not preserved
  int nErased;  nErased = clusters.clearZeroClusters(1);
  if(verbosity>1) {
    cout << blue    << "Number of erased clusters at end was " << magenta 
         << nErased << normText << endl;
  } // end if

  return nStep;
// ---------------------------------------------------------------------------
} // end communityDetectNZ
// ---------------------------------------------------------------------------


int communityDetectNZRedo(TClusterList &clusters, TCMatrix &Cp, 
       TCMData a, TCMData b, int NStepMax = 1000, 
       bool useZeroE = 0, int verbosity = 0, int  debugMode = 0) {
  // implements a neighbor-only search of the system.  It is significanlty 
  // faster.  Note that the implementation uses optimizations that depend on the
  // specific structure of the connection data structures.  Care should be used
  // when changing the classes.
  // all int's should be positive (and non-zero) values

  // -------------------------------------------------------------------------
  // BEGIN MAIN LOOP
  // -------------------------------------------------------------------------
  int     BigInteger = (2<<29);
  int     N = Cp.getSize();  // the number of nodes
  int     iNode, iCluster, jNode, jCluster, minCluster, jDegrees;
  TNode   wNode, *pwCurrent;
  bool    useZeroEnergy = 0;    // do we allow zero energy moves?
  bool    bqConstrained = 0;    // do we contrain the number of clusters?
  int     lastNodeMoved = -2,   thisNodeMoved = -1; // invalid initial values
  int     lastCluster   = -2,   thisCluster   = -1; // invalid initial values
  //int   lastEnergy    = -1,   thisEnergy    = 0; 
  //int     constantEnergyCount = 0, 
  int     zeroEnergyMoveCount = 0; 
  int     energyChange  = 0, deltaEAdded, deltaERemoved;
  int     equalEnergyCount = 0, nStepMoves;
  unsigned long long int moveCount = 0;
  bool    converged = 0,        moved = 0;
  int    *pNeighborStart,      *pNeighbor; // optimzing neighbor pointers
  int    *pjClustereSum;        // optimzing pointer for energy sum
  int     systemEnergy = 0;     // sum of current system energy

  // initialize data for sparse cases
  TVectorInt  edgeCount(N,0); // count number of nodes connected to a cluster
  TVectorInt  eSum(N,0);      // energy sum for each cluster for current node
  TVectorInt  bcChecked(N,0); // was cluster 'x' checked?
  //clusters.initNodesList(); // initialize the node-in-which-cluster array
  TVectorInt  storedEdges(N,0); // track the number of edges of each node within
                                // its own cluster on the previous iteration

  // ************* while nStep ***********************************************
  if(verbosity>0) {
    cout << green << "Starting community detect with ";
    #ifdef MSL_TCMWEIGHTED
    cout << "a " << brown << "weighted";
    #else
    cout << "an " << brown << "unweighted";
    #endif
    cout << green << " neighbor iteration:\n";
  } // end if verbosity
  
  int nStep = 0;
  while(nStep<NStepMax && !converged) {
     clusters.initClean();  // can optimize this by storing moved nodes
     moved = 0; // indicates whether node moves have occured (not converged)
     nStepMoves = 0; // indicates whether node moves have occured (not converged)

    #ifdef DEBUG_MODE
    if(debugMode>0)                                              // debugging
      cout << magenta << "Begin step loop " << nStep << normText << endl 
           << "nClusters = " << clusters.getSize() << endl;
    #endif

    // loop over particles and move each one in turn
    int i = 0;
    while(i<clusters.getSize()) {
      clusters[i].begin();  // begin user iteration
      iCluster = i;         // just to keep clear code

      // ----------------------------------------------------------------
      for(int k=0; k<clusters[i].getSize(); k++) {
        //wNode = clusters[i].getCurrent();
        pwCurrent = &( clusters[i].getCurrent() );
        iNode = clusters[i].getCurrentNode();

        // only test node if it hasn't already been moved this step
        if( !(pwCurrent->clean) && //nStepMoves<MoveThreshold && 
            storedEdges[iNode]>(int)( 0.7*(double)Cp.ki(iNode)+0.5 ) )
          pwCurrent->clean = 0;
        //if(wNode.clean) {
        if(pwCurrent->clean) {
          minCluster = i;               // default cluster is its own

          // I think I can optimize the triple neighbor search by storing the
          // node pointer addresses on the first pass and then iterating 
          // directly over the node addresses on the two following lists.
          // Since this is the inner loop, it may have an impact.

          // scan edge connection list for iNode and sum energy contributions
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7.5% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            //jCluster = clusters.nodes[jNode];
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            // count number of nodes jNode is connected to in jCluster
            edgeCount[jCluster] += 1;
            // sum (subtract) the energy for node iNode merge into jCluster for
            // the unweighted version, we can skip the energy sum and just use
            // the count with the stored and unstored values
            // I use pre-processor if's to avoid an extra inner loop if() call
            #ifdef MSL_TCMWEIGHTED
            eSum[jCluster] -= Cp(iNode,jNode); 
            #else
            eSum[jCluster] -= a; 
            #endif

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and update for unconnected edges
          // first, set the energy change for removing iNode from iCluster
          // which is slightly different from the added case
          eSum[iCluster] += (clusters[iCluster].getSize()-1-edgeCount[iCluster])
                            *(-b);
          deltaERemoved   = eSum[iCluster];
          bcChecked[iCluster] = 1;
          pNeighbor = pNeighborStart;    // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
              // sum (subtract) the energy for node iNode merge into jCluster
              // note that the unstored value is a negative number
              eSum[jCluster] += (clusters[jCluster].getSize()-edgeCount[jCluster])
                                *(-b);
              bcChecked[jCluster] = 1;
            } // end if bcChecked

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // keep track of interior edges
          storedEdges[iNode] = edgeCount[i]; 

          // scan edge connection list and clear changed values
          deltaEAdded = BigInteger;
          pNeighbor   = pNeighborStart;  // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
            // this if statement catches iCluster!=jCluster by default but has 
            // a lot of extraneous assignments.  Not sure which is better.
            //if((*pjClustereSum)<deltaEAdded) {
            // find the lowest energy change but skip the current cluster
            if(iCluster!=jCluster && (*pjClustereSum)<deltaEAdded) {
              deltaEAdded = (*pjClustereSum);
              minCluster  = jCluster;
            } // end if eSum
            // clear the changed parameters - assignements can be redundant
            (*pjClustereSum)    = 0;
            edgeCount[jCluster] = 0;
            bcChecked[jCluster] = 0;

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search
          energyChange = deltaEAdded - deltaERemoved;

          // -------------------------------------------------------------
          // now we have to check the energy and possible moves
          // -------------------------------------------------------------
          // move if energy is lowered, everything is already set if the energy 
          // is lowered, here we check the special cases
          if(energyChange<0) {
            // we have a candidate move, but we must still check to see if 
            // deltaE>0 and if q is unconstrained since we can lower the energy
            // even more in that case
            if(deltaERemoved>0 && !bqConstrained) {  // open new cluster instead...
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead to lower the energy more
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              clusters.addEmpty();                  // declare an empty cluster
              minCluster = clusters.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Opening an empty cluster at position " 
                   << minCluster << normText << endl;  // debugging
              #endif
            } // end empty cluster case
          } else if(energyChange==0 && useZeroEnergy && clusters[i].getSize()>1) {
              // allow zero energy moves as long as we have a  
              // this case is ok - allow the move
              zeroEnergyMoveCount++;
              #ifdef DEBUG_MODE
              cout << green << "Attempting to execute a zero energy move of "
                   << "node " << iNode << " to cluster " << minCluster 
                   << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else if(deltaERemoved>0 && !bqConstrained) {
              // if q is unconstrained, we can still potentially move iNode
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              clusters.addEmpty();                  // declare an empty cluster
              minCluster = clusters.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Adding a preferred empty cluster at position "
                   << minCluster << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else {
            #ifdef DEBUG_MODE
            cout << red << "Invalidating move for node " << iNode 
                 << " in cluster " << iCluster << normText << endl; // debugging
            #endif
            minCluster = i;  // otherwise invalidate this move
            //return 0; // debugging
          } // end else

          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          if(minCluster!=i) {
            moveCount++;
            nStepMoves++; // indicates whether node moves have occured (not converged)
            systemEnergy += energyChange;  // track tot al energy incrementally?

            #ifdef DEBUG_MODE
            if(debugMode>0) 
              cout << red << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            //if(debugMode>1)
              cout << "Known energy change for move " << moveCount << " is:  "
                   << energyChange << " = " 
                   << (-deltaERemoved) << " + " << deltaEAdded
                   << " for node " << iNode << " moving from cluster "
                   << i << " (size " << clusters[i].getSize() << ") to cluster " 
                   << minCluster << " (size " << clusters[minCluster].getSize() 
                   << ") on step " << nStep << ".  The total energy is supposedly " 
                   << systemEnergy << endl; // debugging
            #endif

            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            clusters.nodes[iNode] = minCluster;
            // move from the old to new cluster, after following move, the 
            // 'current' iteration pointer is invalid (points to 'previous')
            clusters[i].moveCurrent(clusters[minCluster],deltaEAdded,deltaERemoved);

            //cout << blueText << "The node to cluster list is:  " << redText 
            //     << clusters.nodes << normText << endl;  // debugging

            moved = 1;  // changes were made, so converged is false, unless...
            pwCurrent->clean = 0;  // invalidate further moves of wNode on this step

            // unless we are moving the same node between the same clusters
            if(lastNodeMoved==thisNodeMoved) {
            //if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << blue  << "Moving node " << iNode << " again on iteration " 
                   << nStep << " from cluster " << i << " to " << minCluster
                   << " with an energy of " << systemEnergy
                   << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

            // unless we have reached the zero energyChange move limit
            if(energyChange<0)        equalEnergyCount = 0;
            else if(energyChange==0)  equalEnergyCount++;
            if(equalEnergyCount>1000 ) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << red << "Ended trial due to excessive consecutive zero "
                   << " energy moves " << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

          } // end if minCluster
          #ifdef DEBUG_MODE
          else  cout << grey  << "Rejected a move of node " << iNode 
                     << " to cluster " << minCluster << normText << endl;
          #endif
          // -------------------------------------------------------------
          // -------------------------------------------------------------
        } // end if wNode.clean

        clusters[i].next();  // increment to next node for user iteration
      } // end for k - particle loop
      // -----------------------------------------------------------------

      //if(clusters[i].getSize()==0 && !bqConstrained) {
      // nStep<3 is an optimization - need to clear zero clusters at the end
      if(clusters[i].getSize()==0 && nStep<3 && !bqConstrained) {
        #ifdef DEBUG_MODE
        if(debugMode>0)  // debugging
          cout << red << "deleting empty cluster " << i << ":  " << clusters[i]
               << normText << endl;
        #endif
        clusters.erase(i);      // does not return memory
        //clusters.destroy(i);  // returns memory
      } // end delete empty cluster
      else  i++;  // if no erase, increment as normal
    } // end while i - cluster loop

    converged = !moved;
    //clusters.display();  // debugging
/*
    if(converged && useZeroE && !useZeroEnergy) {
      useZeroEnergy = 1; // turn on zero energy moves
      converged = 0;     // unset converged status
      #ifdef DEBUG_MODE
      cout << red << "Switching on zero energy moves at energy " << systemEnergy 
           << normText << "\n";
      #endif
    } // end if
*/
    nStep++;
    if(verbosity>1) 
      cout << blue << nStepMoves << " moves on " << nStep << ", " << normText;
  } // end while nStep - iteration loop

  if(verbosity>1) {
    cout << blue << "Done after " << nStep << " steps "  << magenta << "(Made " 
         << moveCount << " moves)" << normText << endl;
  } // end if

  // eliminate size zero clusters before exiting, order is not preserved
  int nErased;  nErased = clusters.clearZeroClusters(1);
  clusters.clearFlags();
  if(verbosity>1) {
    cout << blue    << "Number of erased clusters at end was " << magenta 
         << nErased << normText << endl;
  } // end if

  return nStep;
// ---------------------------------------------------------------------------
} // end communityDetectNZRedo
// ---------------------------------------------------------------------------


int communityDetectNZOpt(TClusterList &c, TCMatrix &Cp, TCMData a, TCMData b, 
       int NStepMax  = 100, bool useZeroEnergy = 0, bool bqConstrained = 0, 
       int verbosity = 0,   int debugMode = 0) {
  // implements a neighbor-only search of the system
  // This optimized version improves with some additional shortcuts including
  // cutting off individual node checks when it is sufficiently `stable.'
  // Note that the implementation uses optimizations that depend on the specific
  // structure of the connection data structures.  
  // Care should be exercised when changing the classes.

  // -------------------------------------------------------------------------
  // BEGIN MAIN LOOP
  // -------------------------------------------------------------------------
  int     N = Cp.getSize();  // the number of nodes
  int     iNode, iCluster, jNode, jCluster, minCluster, jDegrees;
  TNode   wNode, *pwCurrent;
  int     lastNodeMoved = -2,   thisNodeMoved = -1; // invalid initial values
  int     lastCluster   = -2,   thisCluster   = -1; // invalid initial values
  //int   lastEnergy    = -1,   thisEnergy    = 0; 
  int     zeroEnergyMoveCount = 0; 
  int     energyChange  = 0, deltaEAdded, deltaERemoved;
  int     equalEnergyCount = 0, nStepMoves;
  unsigned long long int moveCount = 0;
  bool    converged = 0,      moved = 0;
  int    *pNeighborStart,    *pNeighbor; // optimzing neighbor pointers
  int    *pjClustereSum;      // optimzing pointer for energy sum
  int     systemEnergy = 0;   // sum of current system energy

          b = abs(b);         // use absolute value of b
  int     MoveThreshold = (int)(0.02*(double)N);
  int     addCount = 0;  // debugging
  // initialize data for sparse cases
  TVectorInt  edgeCount(N,0);   // count number of nodes connected to a cluster
  TVectorInt  eSum(N,0);        // energy sum for each cluster for current node
  TVectorInt  bcChecked(N,0);   // was cluster 'x' checked?
  //c.initNodesList();   // initialize the node-in-which-cluster array
  // optimization testing
  TVectorInt  storedEdges(N,0); // track the number of edges of each node within
                                // its own cluster on the previous iteration
  c.clearFlags();

  // ************* while nStep ***********************************************
  if(verbosity>0) {
    cout << green << "Starting community detect with ";
    #ifdef MSL_TCMWEIGHTED
    cout << "a " << brown << "weighted";
    #else
    cout << "an " << brown << "unweighted";
    #endif
    cout << green << " neighbor iteration using ";
    if(bqConstrained)  cout << red   << "constrained";
    else               cout << brown << "unconstrained";
    cout << green << " communities" << normText << "\n";
  } // end if verbosity

  //c.display(0,0,"inside communityDetectNZOpt ");  // debugging

  int nStep = 0;
  while(nStep<NStepMax && !converged) {
     c.initClean();  // can optimize this by storing moved nodes
     moved = 0; // indicates whether node moves have occured (not converged)
     nStepMoves = 0; // indicates how many node moves have occured on this step

    #ifdef DEBUG_MODE
    if(debugMode>1)                                                // debugging
      cout << magenta << "Begin step loop " << nStep << normText
           << "\nnClusters = " << c.getSize() << endl;
    #endif

    // loop over particles and move each one in turn
    int i = 0;
    while(i<c.getSize()) {
      c[i].begin();  // begin user iteration
      //iCluster = i;       // just to keep clear code
      //cout << blue << i << ": " << normText << flush;  // debugging

      // ----------------------------------------------------------------
      for(int k=0; k<c[i].getSize(); k++) {
        //wNode = c[i].getCurrent();
        pwCurrent = &( c[i].getCurrent() );
        iNode = c[i].getCurrentNode();

        //cout << blue << " " << iNode << normText << flush;  // debugging

        //if( !(pwCurrent->flag) && nStepMoves<MoveThreshold && 
        //    storedEdges[iNode]>(int)(0.9*(double)Cp.ki(iNode)+0.5) )
        //  pwCurrent->flag = 1;// mark as having a well-established membership
        // only test node if it hasn't already been moved this step or if it has
        //   significant exterior edge connections
        //if(wNode.clean) {
        //if(pwCurrent->clean && !(pwCurrent->flag)) {
        if(pwCurrent->clean) {  // debugging
          minCluster = i;               // default cluster is its own

          // I think I can optimize the triple neighbor search by storing the
          // node pointer addresses on the first pass and then iterating 
          // directly over the node addresses on the two following lists.
          // Since this is the inner loop, it may have an impact.

          // scan edge connection list for iNode and sum energy contributions
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7.5% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            //jCluster = c.nodes[jNode];
            jNode = (*pNeighbor);
            jCluster = c.nodes[jNode];
            // count number of nodes jNode is connected to in jCluster
            edgeCount[jCluster] += 1;
            // sum (subtract) the energy for node iNode merge into jCluster for
            // the unweighted version, we can skip the energy sum and just use
            // the count with the stored and unstored values
            // I use pre-processor if's to avoid an extra inner loop if() call
            #ifdef MSL_TCMWEIGHTED
            eSum[jCluster] -= Cp(iNode,jNode);
            #else
            eSum[jCluster] -= a;
            #endif

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and update for unconnected edges
          // first, set the energy change for removing iNode from iCluster
          // which is slightly different from the added case
          eSum[i] += (c[i].getSize()-1-edgeCount[i])*b;
          deltaERemoved = -eSum[i];
          deltaEAdded = BigInteger;
          bcChecked[i]  = 1;
          pNeighbor = pNeighborStart;    // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = c.nodes[jNode];
            pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
            if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
              // sum (subtract) the energy for node iNode merge into jCluster
              // note that the unstored value is a negative number
              (*pjClustereSum) += (c[jCluster].getSize()-edgeCount[jCluster])*b;
              bcChecked[jCluster] = 1;
              // this if statement catches iCluster!=jCluster by default but has 
              // a lot of extraneous assignments.  Not sure which is better.
              //if((*pjClustereSum)<deltaEAdded) {
              // find the lowest energy change but skip the current cluster
              //if(i!=jCluster && (*pjClustereSum)<deltaEAdded) {
              if((*pjClustereSum)<deltaEAdded) {
                deltaEAdded = (*pjClustereSum);
                minCluster  = jCluster;
              } // end if eSum
            } // end if bcChecked

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // keep track of interior edges
          storedEdges[iNode] = edgeCount[i]; 

          // scan edge connection list and clear changed values
          pNeighbor = pNeighborStart;  // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = c.nodes[jNode];
            // clear the changed parameters - assignements can be redundant
            eSum[jCluster]      = 0;
            edgeCount[jCluster] = 0;
            bcChecked[jCluster] = 0;

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search
          energyChange = deltaEAdded + deltaERemoved;

          // -------------------------------------------------------------
          // now we have to check the energy and possible moves
          // -------------------------------------------------------------
          // move if energy is lowered, everything is already set if the energy 
          // is lowered, here we check the special cases
          if(energyChange<=0) {
            // we have a candidate move, but we must still check to see if 
            // deltaE>0 and if q is unconstrained since we can lower the energy
            // even more in that case
            if(deltaERemoved<0 && !bqConstrained) {  // open new cluster instead...
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead to lower the energy more
              energyChange = deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              c.addEmpty();                  // declare an empty cluster
#if 0
              addCount++;
              cout << red << "Adding an empty cluster (" << addCount << ")... " 
                   << flush;      // debugging
              cout << "Current size " << c[c.getSize()-1].getSize() 
                   << " cluster at position " << (c.getSize()-1) 
                   << ": " << endl << c[c.getSize()-1] << endl; // debugging
              c.addEmpty();                  // declare an empty cluster
              cout << "Old size " << c[c.getSize()-2].getSize() 
                   << " cluster at position " << (c.getSize()-2) 
                   << ": " << endl << c[c.getSize()-2] << endl; // debugging
              cout << red << "New size " << flush;
              cout << c[c.getSize()-1].getSize() 
                   << " cluster at position " << (c.getSize()-1) 
                   << ": " << endl; 
              cout << c[c.getSize()-1] << endl; // debugging
#endif
              minCluster = c.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              if(debugMode>3) 
                cout << green << "Opening an empty cluster at position " 
                     << minCluster << normText << endl;  // debugging
              #endif
            } // end empty cluster case
          }
          else if(energyChange==0 && useZeroEnergy && c[i].getSize()>1) {
              // allow zero energy moves as long as we have a  
              // this case is ok - allow the move
              zeroEnergyMoveCount++;
              #ifdef DEBUG_MODE
              if(debugMode>3) 
                cout << green << "Attempting to execute a zero energy move of "
                     << "node " << iNode << " to cluster " << minCluster 
                     << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else if(deltaERemoved<0 && !bqConstrained) {
              // if q is unconstrained, we can still potentially move iNode
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead
              energyChange = deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              // in the optimized version, all clusters up to maximum are empty
              c.addEmpty();                  // declare an empty cluster
#if 0
              addCount++;
              cout << red << "Adding an empty cluster (" << addCount << ")... " 
                   << flush;      // debugging
              cout << "Current size " << c[c.getSize()-1].getSize() 
                   << " cluster at position " << (c.getSize()-1) 
                   << ": " << endl << c[c.getSize()-1] << endl; // debugging
              c.addEmpty();                  // declare an empty cluster
              cout << "Old size " << c[c.getSize()-2].getSize() 
                   << " cluster at position " << (c.getSize()-2) 
                   << ": " << endl << c[c.getSize()-2] << endl; // debugging
              cout << red << "New size " << flush;
              cout << c[c.getSize()-1].getSize() 
                   << " cluster at position " << (c.getSize()-1) 
                   << ": " << endl; 
              cout << c[c.getSize()-1] << endl; // debugging
#endif
              minCluster = c.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              if(debugMode>3) 
                cout << green << "Adding a preferred empty cluster at position "
                     << minCluster << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else {
            #ifdef DEBUG_MODE
              if(debugMode>3) 
                cout << red << "Invalidating move for node " << iNode 
                     << " in cluster " << i << normText << endl; // debugging
            #endif
            minCluster = i;  // otherwise invalidate this move
            //return 0; // debugging
          } // end else

          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          if(minCluster!=i) {
            nStepMoves++; // indicates whether node moves have occured (not converged)
            systemEnergy += energyChange;  // track total energy incrementally?

            #ifdef DEBUG_MODE
            if(debugMode>3) 
              cout << red << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            if(debugMode>3) 
              cout << "Known energy change for move " << moveCount << " is:  "
                   << energyChange << " = " 
                   << deltaERemoved << " + " << deltaEAdded
                   << " for node " << iNode << " moving from cluster "
                   << i << " (size " << c[i].getSize() << ") to cluster " 
                   << minCluster << " (size " << c[minCluster].getSize() 
                   << ") on step " << nStep << ".  The total energy is supposedly " 
                   << systemEnergy << endl; // debugging
            #endif

            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            c.nodes[iNode] = minCluster;
            // move from the old to new cluster, after following move, the 
            // 'current' iteration pointer is invalid (points to 'previous')
            c[i].moveCurrent(c[minCluster],deltaEAdded,-deltaERemoved);

            //cout << blueText << "The node to cluster list is:  " << redText 
            //     << c.nodes << normText << endl;  // debugging

            moved = 1;  // changes were made, so converged is false, unless...
            pwCurrent->clean = 0;  // invalidate further moves of wNode on this step

            // unless we are moving the same node between the same clusters
            //if(lastNodeMoved==thisNodeMoved) {
#if 0
            if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              if(debugMode>3) 
                cout << blue  << "Moving node " << iNode << " again on iteration " 
                     << nStep << " from cluster " << i << " to " << minCluster
                     << " with an energy of " << systemEnergy
                     << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved
#endif
            // unless we have reached the zero energyChange move limit
            if(energyChange<0)        equalEnergyCount = 0;
            else if(energyChange==0)  equalEnergyCount++;
            if(equalEnergyCount>1000 ) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              if(debugMode>2) 
                cout << red << "Ended trial due to excessive consecutive zero "
                     << " energy moves " << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

          } // end if minCluster
/*
// debugging
          // -------------------------------------------------------------
          // now we have to check the energy and possible moves
          // -------------------------------------------------------------
          // move if energy is lowered, everything is already set if the energy 
          // is lowered, here we check the special cases
          if(energyChange<0) {
            // we have a candidate move, but we must still check to see if 
            // deltaE>0 and if q is unconstrained since we can lower the energy
            // even more in that case
            if(deltaERemoved>0 && !bqConstrained) {  // open new cluster instead...
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead to lower the energy more
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              c.addEmpty();                  // declare an empty cluster
              minCluster = c.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Opening an empty cluster at position " 
                   << minCluster << normText << endl;  // debugging
              #endif
            } // end empty cluster case
          } else if(energyChange==0 && useZeroEnergy && c[i].getSize()>1) {
              // allow zero energy moves as long as we have a  
              // this case is ok - allow the move
              zeroEnergyMoveCount++;
              #ifdef DEBUG_MODE
              cout << green << "Attempting to execute a zero energy move of "
                   << "node " << iNode << " to cluster " << minCluster 
                   << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else if(deltaERemoved>0 && !bqConstrained) {
              // if q is unconstrained, we can still potentially move iNode
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              c.addEmpty();                  // declare an empty cluster
              minCluster = c.getSize()-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Adding a preferred empty cluster at position "
                   << minCluster << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else {
            #ifdef DEBUG_MODE
            cout << red << "Invalidating move for node " << iNode 
                 << " in cluster " << iCluster << normText << endl; // debugging
            #endif
            minCluster = i;  // otherwise invalidate this move
            //return 0; // debugging
          } // end else
// debugging
// debugging
          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          if(minCluster!=i) {
            moveCount++;
            nStepMoves++; // indicates whether node moves have occured (not converged)
            systemEnergy += energyChange;  // track tot al energy incrementally?

            #ifdef DEBUG_MODE
            if(debugMode>0) 
              cout << red << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            //if(debugMode>1)
              cout << "Known energy change for move " << moveCount << " is:  "
                   << energyChange << " = " 
                   << (-deltaERemoved) << " + " << deltaEAdded
                   << " for node " << iNode << " moving from cluster "
                   << i << " (size " << c[i].getSize() << ") to cluster " 
                   << minCluster << " (size " << c[minCluster].getSize() 
                   << ") on step " << nStep << ".  The total energy is supposedly " 
                   << systemEnergy << endl; // debugging
            #endif

            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            c.nodes[iNode] = minCluster;
            // move from the old to new cluster, after following move, the 
            // 'current' iteration pointer is invalid (points to 'previous')
            c[i].moveCurrent(c[minCluster],deltaEAdded,deltaERemoved);

            //cout << blueText << "The node to cluster list is:  " << redText 
            //     << clusters.nodes << normText << endl;  // debugging

            moved = 1;  // changes were made, so converged is false, unless...
            pwCurrent->clean = 0;  // invalidate further moves of wNode on this step

            // unless we are moving the same node between the same clusters
            if(lastNodeMoved==thisNodeMoved) {
            //if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << blue  << "Moving node " << iNode << " again on iteration " 
                   << nStep << " from cluster " << i << " to " << minCluster
                   << " with an energy of " << systemEnergy
                   << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

            // unless we have reached the zero energyChange move limit
            if(energyChange<0)        equalEnergyCount = 0;
            else if(energyChange==0)  equalEnergyCount++;
            if(equalEnergyCount>1000 ) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << red << "Ended trial due to excessive consecutive zero "
                   << " energy moves " << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

          } // end if minCluster
// debugging
*/
          #ifdef DEBUG_MODE
          else if(debugMode>3) 
                 cout << grey  << "Rejected a move of node " << iNode 
                      << " to cluster " << minCluster << normText << endl;
          #endif
          // -------------------------------------------------------------
          // -------------------------------------------------------------
        } // end if wNode.clean

        c[i].next();  // increment to next node for user iteration
      } // end for k - particle loop
      // -----------------------------------------------------------------
      moveCount += nStepMoves;  // keep track of total number of moves

      //if(c[i].getSize()==0 && nStep<3 && !bqConstrained) {  // debugging omit
      if(c[i].getSize()==0 && !bqConstrained) {
        // erase (effectively swapping with last cluster in the list) for the 
        // first iteration(s) which saves some updates of the node to cluster 
        // list.  This is counter productive if it is too small.  Empirically,
        // it appears that it is at a minimum around a few (say <=2) iterations
        // although the effect depends upon at least n and the system confusion.
        // This effect is internal to the function, and the nodes list and 
        // cluster list are cleaned and re-initialized before exiting.
        //if(c[i].getSize()==0 && !bqConstrained) {
        #ifdef DEBUG_MODE
        if(debugMode>3)  // debugging
          cout << red << "deleting empty cluster " << i << ":  " << c[i]
               << normText << endl;
        #endif
        c.erase(i);      // does not return memory
        //c.destroy(i);  // returns memory
      } // end delete empty cluster
      else  i++;  // if no erase, increment as normal
    } // end while i - cluster loop

    converged = !moved;
    //c.display();  // debugging

    if(!moved && !useZeroEnergy) {
      useZeroEnergy = 1; // turn on zero energy moves
      converged = 0;     // unset converged status
      //nStep--;
      if(verbosity>1)
        cout << red << "Switching on zero energy moves at energy " 
             << systemEnergy << normText << "\n";
    } // end if

    nStep++;
    if(verbosity>2)
      cout << blue << nStepMoves << " moves on " << nStep << ", " << normText;
  } // end while nStep - iteration loop

  if(verbosity>1) {
    cout << blue << "Done after " << nStep << " steps "  << magenta << "(Made " 
         //<< itos(moveCount,1) << " moves)" << normText << endl;
         << moveCount << " moves)" << normText << endl;
  } // end if

  // eliminate size zero clusters before exiting, order is not preserved
  int nErased;  nErased = c.clearZeroClusters(1);
  c.clearFlags();
  if(verbosity>1) {
    cout << blue    << "Number of erased clusters at end was " << magenta 
         << nErased << normText << endl;
  } // end if

  return nStep;
// ---------------------------------------------------------------------------
} // end communityDetectNZOpt
// ---------------------------------------------------------------------------


int communityDetectNZOpt2(TClusterList &clusters, TCMatrix &Cp, 
       TCMData a, TCMData b, int NStepMax = 1000, //int CheckInterval = 100, 
       bool useEnergy = 1,   int TryHard = 3,     int verbosity = 0,   
       int  debugMode = 0) {
  // implements a neighbor-only search of the system.  It is significanlty 
  // faster for 'smaller' systems, but oddly the scaling is worse for large
  // systems even when the clusters are not shuffled as they are moved.
  // Note that the implementation uses optimizations that depend on the 
  // specific structure of the connection data structures.  Care should be
  // used when changing the classes.
  // Experience showed that it actually scales worse on the high size end.
  // This actually does not make sense because it appears to happen even when
  // I allow arbitrary memory usage for the empty clusters (eliminating the
  // nodes-cluster list update.
  // all int's should be positive (and non-zero) values

  // -------------------------------------------------------------------------
  // BEGIN MAIN LOOP
  // -------------------------------------------------------------------------
  int     BigInteger = (2<<29);
  int     N = Cp.getSize();  // the number of nodes
  int     iNode, iCluster, jNode, jCluster, minCluster, jDegrees;
  TNode   wNode, *pwCurrent;
  bool    useZeroEnergy = 0;    // do we allow zero energy moves?
  bool    bqConstrained = 0;    // do we contrain the number of clusters?
  int     lastNodeMoved = -2,   thisNodeMoved = -1; // invalid initial values
  int     lastCluster   = -2,   thisCluster   = -1; // invalid initial values
  //int   lastEnergy    = -1,   thisEnergy    = 0; 
  //int     constantEnergyCount = 0, 
  int     zeroEnergyMoveCount = 0; 
  int     energyChange  = 0, deltaEAdded, deltaERemoved;
  unsigned  equalEnergyCount = 0, nStepMoves;
  unsigned long moveCount = 0;
  bool    converged = 0,        moved = 0;
  int    *pNeighborStart,      *pNeighbor; // optimzing neighbor pointers
  int    *pjClustereSum;        // optimzing pointer for energy sum
  int     systemEnergy = 0;     // sum of current system energy

          b = abs(b);           // use absolute value of b
  int     MoveThreshold = (int)(0.01*(double)N);
  int     addCount = 0;  // debugging
  // initialize data for sparse cases
  TVectorInt  edgeCount(10*N,0);   // count number of nodes connected to a cluster
  TVectorInt  eSum(10*N,0);        // energy sum for each cluster for current node
  TVectorInt  bcChecked(10*N,0);   // was cluster 'x' checked?
  //clusters.initNodesList();   // initialize the node-in-which-cluster array
  // optimization testing
  TVectorInt  storedEdges(N,0); // track the number of edges of each node within
                                // its own cluster on the previous iteration
  TVectorInt  nodeFlag(N,0);
  //TVectorInt  nodeClean(N);     // initialized at beginning of the loop
  TVectorInt  nodeCount(10*N,1);   // tracks the number of nodes in each cluster
  int         nClusters = clusters.getSize(); // start with symmetric list

  // for the moment, warn user if optimizing by modularity since the modularity
  // functions were not implemented in the new linked-list classes
  if(!useEnergy)
    errorMsg("communityDetectNZ is not currently implemented for direct modularity optimization");
  // note that this is not implemented for sub-cluster runs at the moment!!!!!
  //warningMsg("Note:  communityDetectNZ is not currently implemented for sub-cluster trials.",
  //            brown);
  //cout << green << "The Cp gammaa = " << a << " and gammab = " << b
  //     << normText << "\n"; // debugging

  // ************* while nStep ***********************************************
  //if(verbosity>0) {
    cout << green << "Starting community detect with ";
    #ifdef MSL_TCMWEIGHTED
    cout << "a " << brown << "weighted";
    #else
    cout << "an " << brown << "unweighted";
    #endif
    cout << green << " neighbor iteration:\n";
  //} // end if verbosity

  int nStep = 0, iSearch;
  while(nStep<NStepMax && !converged) {
    moved = 0;      // indicates whether node moves have occured (not converged)
    nStepMoves = 0; // indicates how many node moves have occured on this step

    //clusters.initClean();  // can optimize this by storing moved nodes
    //nodeClean.init(1);  // can optimize this by storing moved nodes
    //iSearchStart = randomInt(0,N-1);

    clusters.randomizeOrder();  // randomizes the node order for each step
    //clusters.display(); // debugging

    #ifdef DEBUG_MODE
    cout << magenta << "Begin step loop " << nStep << normText << endl 
         << "nClusters = " << nClusters << endl;          // debugging
    cout << magenta << "Iteration nStep = " << nStep << " with " << nClusters 
         << "total clusters:  " << flush;  // debugging
    #endif
    // loop over particles and move each one in turn

    // iterate from end to prevent a lot updating on the first step
    int i = N-1;
    //while(i<N) {
    while(i>=0) {
      //cout << green << "Starting (i = " << i << flush; // debugging
      //clusters[i].begin();  // begin user iteration
      //cout << blue << i << ": " << normText << flush;  // debugging

      // ----------------------------------------------------------------
      //for(int k=0; k<clusters[i].getSize(); k++) {
        //iSearch = (iSearch+1)%N;
        // for the moment, we assume a symmetric clusters parameter, i.e. one 
        // node per cluster
        iNode = clusters[i].getStartNode(); // ad-hoc randomize search order
        //cout << ") with node " << flush;  // debugging
        iCluster = clusters.nodes[iNode];  // just to keep clear code
        //cout << iNode << normText << endl;  // debugging

        //if(wNode.clean) {
        if(nodeFlag[iNode]==0 && nStepMoves<MoveThreshold && 
           storedEdges[iNode] > (int)(0.8*(double)Cp.ki(iNode)+0.5) )
          nodeFlag[iNode] = 1; // mark as having a well-established membership
        // only test node if it hasn't already been moved this step or if it has
        //   significant exterior edge connections
        //if(nodeClean[iNode]==1 && nodeFlag[iNode]==0) {
        if(nodeFlag[iNode]==0) {
          minCluster = iCluster;               // default cluster is its own

          // scan edge connection list for iNode and sum energy contributions
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7.5% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            //jCluster = clusters.nodes[jNode];
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            // count number of nodes jNode is connected to in jCluster
            edgeCount[jCluster] += 1;
            // sum (subtract) the energy for node iNode merge into jCluster for
            // the unweighted version, we can skip the energy sum and just use
            // the count with the stored and unstored values
            // I use pre-processor if's to avoid an extra inner loop if() call
            #ifdef MSL_TCMWEIGHTED
            eSum[jCluster] -= Cp(iNode,jNode); 
            #else
            eSum[jCluster] -= a; 
            #endif

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and update for unconnected edges
          // first, set the energy change for removing iNode from iCluster
          // which is slightly different from the added case
          //eSum[i] += (clusters[i].getSize()-1-edgeCount[i])*b;
          eSum[iCluster] += (nodeCount[iCluster]-1-edgeCount[iCluster])*b;
          deltaERemoved = eSum[iCluster];
          bcChecked[iCluster]  = 1;
          pNeighbor = pNeighborStart;    // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
              // sum (subtract) the energy for node iNode merge into jCluster
              // note that the unstored value is a negative number
              //eSum[jCluster] += (clusters[jCluster].getSize()-edgeCount[jCluster])*b;
              eSum[jCluster] += (nodeCount[jCluster]-edgeCount[jCluster])*b;
              bcChecked[jCluster] = 1;
            } // end if bcChecked

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // keep track of interior edges
          storedEdges[iNode] = edgeCount[iCluster]; 

          // scan edge connection list and clear changed values
          deltaEAdded = BigInteger;
          pNeighbor   = pNeighborStart;  // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = clusters.nodes[jNode];
            pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
            // this if statement catches iCluster!=jCluster by default but has 
            // a lot of extraneous assignments.  Not sure which is better.
            //if((*pjClustereSum)<deltaEAdded) {
            // find the lowest energy change but skip the current cluster
            if(i!=jCluster && (*pjClustereSum)<deltaEAdded) {
              deltaEAdded = (*pjClustereSum);
              minCluster  = jCluster;
            } // end if eSum
            // clear the changed parameters - assignements can be redundant
            (*pjClustereSum)    = 0;
            edgeCount[jCluster] = 0;
            bcChecked[jCluster] = 0;

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search
          energyChange = deltaEAdded - deltaERemoved;

          // -------------------------------------------------------------
          // now we have to check the energy and possible moves
          // -------------------------------------------------------------
          // move if energy is lowered, everything is already set if the energy 
          // is lowered, here we check the special cases
          if(energyChange<0) {
            // we have a candidate move, but we must still check to see if 
            // deltaE>0 and if q is unconstrained since we can lower the energy
            // even more in that case
            if(deltaERemoved>0 && !bqConstrained) {  // open new cluster instead...
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead to lower the energy more
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              //clusters.addEmpty();                  // declare an empty cluster
              //minCluster = clusters.getSize()-1;    // move to the new cluster
              nodeCount[nClusters] = 0;
              nClusters++;
              minCluster = nClusters-1;    // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Opening an empty cluster at position " 
                   << minCluster << normText << endl;  // debugging
              #endif
            } // end empty cluster case
          } else if(energyChange==0 && useZeroEnergy && nodeCount[iCluster]>1) {
              // allow zero energy moves as long as we have a  
              // this case is ok - allow the move
              zeroEnergyMoveCount++;
              #ifdef DEBUG_MODE
              cout << green << "Attempting to execute a zero energy move of node " 
                   << iNode << " to cluster " << minCluster << normText << endl; // debugging
              #endif
          } // end empty cluster case
          else if(deltaERemoved>0 && !bqConstrained) {
              // if q is unconstrained, we can still potentially move iNode
              // if we have unconstrained communities with a positive energy
              // move, move to an empty cluster instead
              energyChange = -deltaERemoved;
              deltaEAdded  = 0;
              // declare an empty cluster for the next move statement
              // in the optimized version, all clusters up to maximum are empty
              //clusters.addEmpty();                 // declare an empty cluster
              //minCluster = clusters.getSize()-1;   // move to the new cluster
              nodeCount[nClusters] = 0;
              nClusters++;
              minCluster = nClusters-1;              // move to the new cluster
              #ifdef DEBUG_MODE
              cout << green << "Adding a preferred empty cluster at position "
                   << minCluster << normText << endl;  // debugging
              #endif
          } // end empty cluster case
          else {
            #ifdef DEBUG_MODE
            cout << red << "Invalidating move for node " << iNode 
                 << " in cluster " << i << normText << endl; // debugging
            #endif
            minCluster = iCluster;  // otherwise invalidate this move
            //return 0; // debugging
          } // end else

          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          //if(minCluster!=i) {
          if(minCluster!=iCluster) {
            nStepMoves++; // indicates whether node moves have occured (not converged)
            systemEnergy += energyChange;  // track tot al energy incrementally?

            #ifdef DEBUG_MODE
            if(debugMode>0) 
              cout << red << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            //if(debugMode>1)
              cout << "Known energy change for move " << moveCount << " is:  "
                   << energyChange << " = " 
                   << (-deltaERemoved) << " + " << deltaEAdded
                   << " for node " << iNode << " moving from cluster "
                   << iCluster << " (size " << nodeCount[iCluster] << ") to cluster " 
                   << minCluster << " (size " << nodeCount[minCluster]
                   << ") on step " << nStep << ".  The total energy is supposedly " 
                   << systemEnergy << endl; // debugging
            #endif

            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            clusters.nodes[iNode] = minCluster;
            // move from the old to new cluster, after following move, the 
            // 'current' iteration pointer is invalid (points to 'previous')
            //clusters[i].moveCurrent(clusters[minCluster],deltaEAdded,deltaERemoved);
            nodeCount[minCluster] += 1;
            nodeCount[iCluster]   -= 1;      
/*
            // check if we should delete iCluster
            // this adjustment is very expensive, need to fix it
            if(nodeCount[iCluster]==0 && nStep==0 && nClusters== && !bqConstrained) {
              // erase (effectively swap with last cluster in the list) for the 
              // first iteration(s) which saves some updates of the node to 
              // cluster list. This version in 'NZMod' is more expensive since 
              // we do not explicitly know the members of cluster[nClusters-1] 
              // for the swap.
              #ifdef DEBUG_MODE
              //if(debugMode>0)  // debugging
                cout << red << "deleting empty cluster " << iCluster << " size "
                     << nodeCount[iCluster] << " with particle number "
                     << iNode << "... " << normText << flush;
              #endif
              // have to traverse the nodes list to reset cluster assignments 
              // for the last cluster since it is implicitly swapped to the i'th 
              // position.  Obviously the order is not preserved
              // traverse list only if deleting a cluster not already at the end
              if(iCluster<nClusters-1) {
                int *pc = &(clusters.nodes[0]);  // optimizing pointer
                for(int k=0; k<N; k++) {
                  //if(clusters.nodes[k]==nClusters-1)  
                  //  clusters.nodes[k] = iCluster;
                  if((*pc)==nClusters)  (*pc) = iCluster;
                  pc++;
                } // end for k
                nodeCount[iCluster] = nodeCount[nClusters-1];
              } // end if iCluster
              nClusters--;
              //clusters.erase(i);      // does not return memory
              //clusters.destroy(i);  // returns memory
            } // end delete empty cluster
*/
            //cout << blueText << "The node to cluster list is:  " << redText 
            //     << clusters.nodes << normText << endl;  // debugging

            moved = 1;  // changes were made, so converged is false, unless...
            //nodeClean[iNode] = 0;  // invalidate further moves of wNode on this step
/*
            // unless we are moving the same node between the same clusters
            //if(lastNodeMoved==thisNodeMoved) {
            if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
              moved = 0;  // repeated moves of same node do not count
              #ifdef DEBUG_MODE
              cout << blue  << "Moving node " << iNode << " again on iteration " 
                   << nStep << " from cluster " << i << " to " << minCluster
                   << " with an energy of " << systemEnergy
                   << normText << endl;  // debugging
              #endif
            } // end if lastNodeMoved

            // unless we have reached the zero energyChange move limit
            if(energyChange<0)        equalEnergyCount = 0;
            else if(energyChange==0)  equalEnergyCount++;
            if(equalEnergyCount>1000 ) {
              moved = 0;  // repeated moves of same node do not count
              //#ifdef DEBUG_MODE
              cout << red << "Ended trial due to excessive consecutive zero "
                   << " energy moves " << normText << endl;  // debugging
              //#endif
            } // end if lastNodeMoved
*/
          } // end if minCluster
          #ifdef DEBUG_MODE
          else  cout << grey  << "Rejected a move of node " << iNode 
                     << " to cluster " << minCluster << normText << endl;
          #endif
          // -------------------------------------------------------------
          // -------------------------------------------------------------
        } // end if wNode.clean

        //clusters[i].next();  // increment to next node for user iteration
      //} // end for k - particle loop
      // -----------------------------------------------------------------

      moveCount += nStepMoves;  // keep track of total number of moves
      i--;  // decrement to next node in the list
      //cout << green << "done with node " << iNode << normText << endl; // debugging
    } // end while i - main node loop

    converged = !moved;
    //clusters.display();  // debugging
/*
    if(!moved && !useZeroEnergy) {
      useZeroEnergy = 1; // turn on zero energy moves
      converged = 0;     // unset converged status
      //nStep--;
      #ifdef DEBUG_MODE
      cout << red << "Switching on zero energy moves at energy " << systemEnergy 
           << normText << "\n";
      #endif
    } // end if
*/
    nStep++;
    //if(verbosity>1) {
      cout << blue << nStepMoves << " moves on " << nStep << ", " << normText;
    //} // end if
  } // end while nStep - iteration loop

  //if(verbosity>1) {
    cout << blue << "Done after " << nStep << " steps "  << magenta << "(Made " 
         << itos(moveCount,1) << " moves)" << normText << endl;
  //} // end if
/*
  // now use the nodes list to build the answer in the clusters list
  // reset the node list to sequential order
  //for(int k=0; k<N; k++)  clusters[k].getStart().node = k; // cheating here
  // now reorder the list according to the nodes solution, do this without
  // moving clusters around to make the logic easier
  TVectorInt nodesList;  nodesList = clusters.nodes; // wasteful
  clusters.erase();  // erase current list - memory intensive unfortunately
  clusters.nClusters = N;  // bad form
  int *pc = &(nodesList[0]);  // optimizing pointer
  for(int i=0; i<N; i++) {
    clusters[(*pc)].add(i);
    pc++;
  } // end for k
  // now get rid of empty clusters
  clusters.clearZeroClusters();
  // energy changes were not stored, so calculate energy from scratch
  clusters.calcEnergy(Cp);
*/
  return nStep;
// ---------------------------------------------------------------------------
} // end communityDetectNZOpt2
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// main community detection functions
// ---------------------------------------------------------------------------
//int communityDetect(TClusterList &clusters, TCMDense &Cp, 
// I has having trouble with compilation with the dense highland iteration 
// for the next line. It was due to purposefully omitting the virtual function
// calls and the compiler could not resolve the function call properly.
int communityDetect(TClusterList &clusters, TCMatrix &Cp, 
       int  NStepMax = 1000, int CheckInterval = 100, 
       bool bDense = 1,      bool useEnergy = 1,      int TryHard = 3, 
       int  verbosity = 0,   int debugMode = 0) {
       // all int's should be positive (and non-zero) values

  // -------------------------------------------------------------------------
  // BEGIN MAIN LOOP
  // -------------------------------------------------------------------------
  int NMax = Cp.getSize(), N = NMax; // N can be less than NMax on sub-cluster trials
  const bool HasMoved = 0;
  const int  BigInteger = 2 << 29;
  int     nStep = 0, iNode, minCluster;
  TNode   wNode;
  int     lastNodeMoved = -2,   thisNodeMoved = -1; // invalid initial values
  int     lastCluster   = -2,   thisCluster   = -1; // invalid initial values
  int     lastEnergy    = -1,   thisEnergy    = 0; 
  //double  lastMod       = -1,   thisMod       = 0; 
  int     constantEnergyCount = 0; 
  int     energyChange  = 0,    energyAdded, energyRemoved;
  double  modChange     = 0.0,  modAdded,    modRemoved;
  int     equalEnergyCount = 0, lowEnergy = BigInteger, moveCount = 0;
  int     lowEnergyAdded,       lowEnergyRemoved;
  double  highMod = -10000.0,   highModAdded, highModRemoved;
  bool    converged = 0,        moved = 0;
  bool    useZeroEnergy = 0;    // bool    useZeroEnergy = 0;
  bool    bqConstrained = 0;    // do we fix q or allow it to vary

  //Cp.display();  // debugging
  //errorMsg("Done.");

  // initialize data for sparse cases
  TVectorInt  bClusterChecked(NMax,0);// how big is enough? NMax seems too big
  clusters.initNodesList();

  // for the moment, warn user if optimizing by modularity since the modularity
  // functions were not implemented in the new linked-list classes
  if(!useEnergy)
    errorMsg("communityDetect is not currently implemented for direct modularity optimization");

  // ************* while nStep ***********************************************
  if(verbosity>1) {
    cout << greenText << "Starting community detect with a " << brownText;
    if(bDense)  cout << "dense ";  else  cout << "sparse "; 
    cout << greenText << "cluster iteration:  ";
  } // end if verbosity

  //cout << "clusters size = " << clusters.getSize() << flush; // debugging
  if(verbosity>1)  cout << "Periodic energies are:  " << flush;
  while(nStep<NStepMax && !converged) {
    if(debugMode>0)  cout << brownText << "Beginning loop " << normText << endl;

    nStep++;
    if(nStep%CheckInterval==0)  if(verbosity>0)  cout << nStep << ": " << flush;
    clusters.initClean();
    moved = 0;

    #ifdef DEBUG_MODE
    if(debugMode>0 && nStep&CheckInterval==1) {                    // debugging
      clusters.display(0);                                         // debugging
      cout << magentaText << "For step " << nStep << normText << "\n";
    } // end if debugMode                                          // debugging
    //cout << magentaText << "Original clusters " << normText << "\n";
    //clusters.display();  // debugging

    //cout << nStep << " " << flush;       // debugging
    #endif

    // loop over particles and move each one in turn
    int i = 0;
    while(i<clusters.getSize()) {

      #ifdef DEBUG_MODE
      if(debugMode>0)                                              // debugging
        cout << magentaText << "Begin clusters loop " << normText << endl;
      if(bClusterChecked.getSize()!=N) {                           // debugging
        cout << redText << "bClusterChecked size = "               // debugging
          << bClusterChecked.getSize() << " and i = " << i << endl;// debugging
        return 0;                                                  // debugging
      } // end if                                                  // debugging
      #endif

      clusters[i].begin();  // begin user iteration

      //cout << magentaText << "cluster[i] " << normText << "\n";
      //clusters[i].display();  // debugging

      for(int j=0; j<clusters[i].getSize(); j++) {
        wNode = clusters[i].getCurrent();
        iNode = clusters[i].getCurrentNode();

        // only test node if it hasn't already been moved this step
        if(wNode.clean) {
          lowEnergy      = BigInteger;  // set as extremely high value
          lowEnergyAdded = BigInteger;  lowEnergyRemoved = BigInteger;  
          minCluster = i;               // default cluster is its own
          energyRemoved  = clusters[i].getEnergyChangeRemoved(wNode,Cp);
          //highMod      = -10000.0;    // set as invalid values
          //highModAdded = -10000.0;    highModRemoved = -10000.0;  
          //modRemoved   = clusters[i].getModChangeRemoved(wNode,Cp,ki);

          // check energy change for each cluster
          int k, nClustersToCheck, wCluster;  
          if(bDense) { // if graph is dense, search all clusters sequentially
            if(i==0) k = 1;  else k = 0;            // if i=0 start at k=1
            nClustersToCheck = clusters.getSize();  // end search w/last cluster
            wCluster = k;       // use wCluster as proxy variable for bDense
          } else {              // search only clusters of connected nodes
            // sparse case where we search only connected nodes 
            #ifdef DEBUG_MODE
            if(debugMode>0)  cout << brown << "A " << flush; // debugging
            clusters.initNodesList();  // debugging
            if(bClusterChecked.getSize()!=N) {                   // debugging
              cout << red << "bClusterChecked size = "       // debugging
                   << bClusterChecked.getSize() << " and k = "   // debugging
                   << k << endl;                                 // debugging
              return 0;                                          // debugging
            } // end if                                          // debugging
            #endif

            // -----------------------------------------------------------------
            // implemented using a constant connection list
            nClustersToCheck = Cp.ki(iNode);
            // k==0==nClustersToCheck skips the case of no connected edges
            k = 0;                  // otherwise start counting
            if(nClustersToCheck>0)  wCluster = clusters.nodes[Cp.ki(iNode,k)];

            while(k<nClustersToCheck && (wCluster<0 || wCluster==i)) {
              k++;
              if(k<nClustersToCheck)  wCluster = clusters.nodes[Cp.ki(iNode,k)];
            } // end while

            #ifdef DEBUG_MODE
            if(debugMode>0)                                        // debugging
              cout << brownText << "D with iNode = " << iNode      // debugging
                   << ", wCluster = " << wCluster  << " k = " << k // debugging 
                   << ", and nClustersToCheck = " << nClustersToCheck << endl;
                   // << blueText << Cp.ki(iNode) << endl;         // debugging
            #endif
            // -----------------------------------------------------------------
          } // end if else bDense

          #ifdef DEBUG_MODE
          if(debugMode>0 && wCluster==-1) {                        // debugging
            cout << redText << "wCluster is -1 again with k = " << k << endl;
            return 0;                                              // debugging
          } // end if                                              // debugging
          if(debugMode>0 && nStep&CheckInterval==0)  cout << "\n"; // debugging
          if(debugMode>1)  cout << brownText << "E " << flush;     // debugging
          if(bClusterChecked.getSize()!=N) {                       // debugging
            cout << redText << "bClusterChecked size = "           // debugging
                 << bClusterChecked.getSize() << " and k = "       // debugging
                 << k << endl;                                     // debugging
            return 0;                                              // debugging
          } // end if                                              // debugging
          #endif

          while(k<nClustersToCheck) {
            energyAdded = clusters[wCluster].getEnergyChangeAdded(wNode,Cp);
            energyChange = energyAdded + energyRemoved;
            // modularity functions are not currently implemented in new classes
            //modAdded     = clusters[k].getModChange(iNode,Cp,ki);
            //modChange    = modAdded + modRemoved;

            #ifdef DEBUG_MODE
            if(debugMode>0 && nStep&CheckInterval==0 && energyChange<=0)
              cout << "Energy change is " << energyChange          // debugging
                   << " for node " << iNode  << " moving from cluster "
                   << i << " to " << wCluster << " and wCluster has " 
                   << clusters[wCluster].getSize() << " nodes." << endl; // debugging
            #endif

            // now test versus current minimum energy change so far
            // The equal energy move appears to help convergence at times.
            if(useEnergy) {
              if(energyChange<lowEnergy && (
                 energyChange<0 || 
                 // allow user to turn off zero energy moves until the first 
                 // local minimum has been reached which speeds convergence in 
                 // some cases
                 ( energyChange==0 && useZeroEnergy 
                 // prevent single node from bouncing between between empty 
                 // communities
                 //&& !(clusters[k].getSize()==0 && clusters[i].getSize()==1)
                ))) {

                #ifdef DEBUG_MODE
                if(debugMode>1)
                  cout << "Found lower energy of " << energyChange // debugging
                       << " compared to " << lowEnergy << " for "  // debugging
                       << " cluster " << k << endl;                // debugging
                #endif

                lowEnergy      = energyChange; minCluster = wCluster;
                lowEnergyAdded = energyAdded;  lowEnergyRemoved = energyRemoved;
                //highMod      = modChange;
                //highModAdded = modAdded;     highModRemoved = modRemoved;
              } // end if energyChange
            } // end if useEnergy
/*
            // try a modularity optimization solution for comparison instead
            else if(modChange>highMod && modChange>0.0) {
              lowEnergy      = energyChange;  minCluster = wCluster;
              lowEnergyAdded = energyAdded;   lowEnergyRemoved = energyRemoved;
              //highMod        = modChange;     
              //highModAdded   = modAdded;      highModRemoved = modRemoved;
            } // end else useEnergy (use modularity measure instead of energy)
*/
            #ifdef DEBUG_MODE
            if(debugMode>1)  cout << brownText << "H " << flush;   // debugging
            #endif

            if(bDense) {          // search all clusters
              k++;  if(k==i) k++; // increment k, but skip k==i
              wCluster = k;       // use wCluster as proxy variable
            } else {              // search only clusters of connected nodes
              // sparse case where we only search connected nodes

              #ifdef DEBUG_MODE
              if(wCluster>=bClusterChecked.getSize()) {            // debugging
                cout << red << "wCluster at " << wCluster      // debugging
                     << " is too big with k = " << k << endl;      // debugging
                return 0;                                          // debugging
              } // end if                                          // debugging
              if(debugMode>0)  cout << brownText << "I " << flush; // debugging
              if(wCluster<0)                                       // debugging
                cout << red << "wCluster<0 again " << wCluster // debugging
                     << " at k = " << k << " and i = " << endl;    // debugging
              if(wCluster>N)                                       // debugging
                cout << red << "wCluster<0 again " << wCluster // debugging
                     << " at k = " << k << " and i = " << endl;    // debugging
              if(bClusterChecked.getSize()!=N) {                   // debugging
                cout << red << "bClusterChecked size = "       // debugging
                     << bClusterChecked.getSize() << " and k = " << k << endl;
                return 0;                                          // debugging
              } // end if                                          // debugging
              #endif

              if(wCluster>=0)  bClusterChecked[wCluster] = 0;

              #ifdef DEBUG_MODE
              if(bClusterChecked.getSize()!=N) {                   // debugging
                cout << red << "bClusterChecked size = "       // debugging
                     << bClusterChecked.getSize() << " and k = " << k << endl;
                return 0;                                          // debugging
              } // end if                                          // debugging
              if(debugMode>0)            
                          // debugging
                cout << blue << "wCluster = " << wCluster << ", k = " << k 
                     << " and minCluster = " << minCluster         // debugging
                     << " with ki[iNode] = " << Cp.ki(iNode) << endl; // debugging
              #endif

              do {
                k++;

                #ifdef DEBUG_MODE
                if(debugMode>0)  cout << red << "."            // debugging
                   << "wCluster = " << wCluster << flush;          // debugging
                #endif

                if(k<nClustersToCheck)  wCluster=clusters.nodes[Cp.ki(iNode,k)];

                #ifdef DEBUG_MODE
                if(debugMode>0)  
                  cout << brownText << "wCluster = " << wCluster   // debugging
                       << " and i = " << i << ", bClusterChecked"  // debugging
                       << " size = "  << bClusterChecked.getSize() << endl;
                if(debugMode>0)  cout << clusters[i] << endl;      // debugging
                #endif

              } while(k<nClustersToCheck && (
                      wCluster<0  || // hack fix to handle sub-cluster instances
                      wCluster==i || // don't run on same cluster
                      bClusterChecked[wCluster]==1)); // don't repeat a cluster

              #ifdef DEBUG_MODE
              if(bClusterChecked.getSize()!=N) {                   // debugging
                cout << red << "bClusterChecked size = "       // debugging
                     << bClusterChecked.getSize() << " and k = " << k << endl; 
                return 0;                                          // debugging
              } // end if                                          // debugging
              if(debugMode>1 && wCluster==-1) {                    // debugging
                cout << red << "wCluster is -1 again with k = " << k << endl
                     << "and " << bClusterChecked << endl;         // debugging
                return 0;                                          // debugging
              } // end if                                          // debugging
              #endif
            } // end else for if bDense

            #ifdef DEBUG_MODE
            if(debugMode>0)  cout << brownText << " J " << flush;  // debugging
            if(debugMode>1)                                        // debugging
              cout << blueText << "wCluster = " << wCluster << ", k = " << k 
                   << " and minCluster = " << minCluster << " with ki[iNode] = "
                   << Cp.ki(iNode) << endl;                           // debugging
            if(debugMode>1)                                        // debugging
              cout << "Known energies for move " << moveCount << " are:  " 
                   << energyChange  << " =? " << energyAdded << " + " 
                   << energyRemoved << " for node " << iNode << " moving from "
                   << "cluster " << i << " to " << minCluster << endl;// debugging
            #endif
          } // end while k - cluster loop

          // clear checked nodes from loop for sparse version
          if(!bDense) { 
            for(int m=0; m<nClustersToCheck; m++) {
              wCluster = clusters.nodes[Cp.ki(iNode,m)];  // reuse wCluster var
              if(wCluster>=0)  bClusterChecked[wCluster] = 0;
            } // end for m
          } // end if !bDense

          #ifdef DEBUG_MODE
          if(debugMode>0)                                          // debugging
            cout << blueText << "wCluster = " << wCluster << ", k = " << k 
                 << " and minCluster = " << minCluster << " with ki[iNode] = "
                 << Cp.ki(iNode) << endl;                             // debugging
          if(debugMode>0)                                          // debugging
            cout << redText << "Attempting to move node " << iNode << " from "
                 << "cluster " << i << " to " << minCluster << normText << endl;
          #endif

          // now add the node to the cluster where it had the lowest energy 
          // change (minCluster) and mark it moved
          if(minCluster!=i) {
            moveCount++;

            #ifdef DEBUG_MODE
            if(debugMode>0 && nStep&CheckInterval==0) 
              cout << redText << "Moving node " << iNode    << " from cluster "
                   << i << " to " << minCluster << normText << endl;// debugging
            #endif

            // keep track of last node moved, if these are the same after two 
            // complete loops, the node is bouncing between different clusters
            // and we should terminate the run???
            // if changes are made that allow a single node to move twice in
            // in one loop, this logic would have to be changed.
            lastNodeMoved = thisNodeMoved;  thisNodeMoved = iNode;
            lastCluster   = thisCluster;    thisCluster   = minCluster;

            #ifdef DEBUG_MODE
            //if(debugMode>=0 && energyChange<=0) 
            if(verbosity>3) 
              cout << "Known energy change for move " << moveCount << " is:  "
                   << lowEnergy << " for node " << iNode << " moving from "
                   << "cluster " << i << " to " << minCluster 
                   << " which has " << clusters[minCluster].getSize() 
                   << " nodes." << endl; // debugging
            if(verbosity>2) 
              cout << "Low energyies are:  added = " << lowEnergyAdded 
                   << "for node " << iNode << " moving from "
                   << "cluster " << i << " to " << minCluster 
                   << " which has " << clusters[minCluster].getSize() 
                   << " nodes." << endl; // debugging
            if(verbosity>3) 
              cout << magentaText << "i before with " << clusters.getSize() 
                   <<  " total clusters" << normText << endl;  // debugging
            #endif

            // move from the old to new cluster
            // after following move, the 'current' iteration pointer is invalid
            // (actually afterwards points to 'previous')
            //cout << brownText << " a " << flush; // debugging
            if(lowEnergy>0 && !bqConstrained) {
              // if we have unconstrained communities with a positive energy 
              // move, move to an empty cluster instead
              lowEnergyAdded = 0;                   // new empty cluster has E=0
              // declare an empty cluster for the next move statement
              clusters.addEmpty();                  // declare an empty cluster
              minCluster = clusters.getSize()-1;    // move to the new cluster
              cout << greenText  << "Adding an empty cluster at position " 
                   << minCluster << normText << endl;  // debugging
            } // end empty cluster case
/*
            cout << brownText << " b " << flush; // debugging

            // now perform regular move for either case
            cout << magentaText << "i before with " << clusters.getSize() 
                 <<  " total clusters" << normText << endl;  // debugging
            cout << brownText << "Step j = " << j << " with cluster i size = " 
                 << clusters[i].getSize() << " and minCluster size = " 
                 << clusters[minCluster].getSize() << normText << endl;  // debugging
            cout << redText  << "display i for i = " << i << " and j = " << j 
                 << normText << endl;  // debugging
            clusters[i].display();  // debugging
            clusters[minCluster].display();  // debugging

            cout << blueText << "\nnode move for node " 
                 << clusters[i].getCurrentNode() << " from " << i << " to " 
                 << minCluster << "\nThe node to cluster list is:  " 
                 << clusters.nodes << normText << endl;  // debugging
*/
            // update node to cluster list - this individual move is manually
            // updated due to the class structure - that is, we call the cluster
            // move function rather than a cluster list move
            clusters.nodes[clusters[i].getCurrentNode()] = minCluster;
            clusters[i].moveCurrent(clusters[minCluster],lowEnergyAdded,lowEnergyRemoved);
/*
            cout << blueText << "The node to cluster list is:  " << redText 
                 << clusters.nodes << normText << endl;  // debugging

            cout << magentaText << "i after  " << normText << endl;
            cout << brownText << "Step two j = " << j << " with cluster i size = " 
                 << clusters[i].getSize() << " and minCluster size = " 
                 << clusters[minCluster].getSize() << normText << endl;  // debugging
            cout << magentaText << "display now " << normText << endl;  // debugging
            if(clusters[i].getSize()>0)  clusters[i].display();  // debugging
            else  cout << "cluster i is empty with size " 
                       << clusters[i].getSize() << endl;  // debugging
            cout << magentaText << "done " << normText << endl;  // debugging
            cout << brownText << " c " << flush; // debugging
            clusters[minCluster].display();  // debugging
*/
            moved = 1;  // changes were made, so converged is false...
            // unless we have reached the zero energyChange move limit
            if(lowEnergy<0)        equalEnergyCount = 0;
            else if(lowEnergy==0)  equalEnergyCount++;
            //cout << brownText << " d " << flush; // debugging

            #ifdef DEBUG_MODE
            if(bClusterChecked.getSize()!=N) {                     // debugging
              cout << redText << "bClusterChecked size = "         // debugging
                   << bClusterChecked.getSize() << " and k = "     // debugging
                   << k << endl;                                   // debugging
              return 0;                                            // debugging
            } // end if                                            // debugging
           #endif


          } // end if minCluster

        } //end if wNode clean
  
        #ifdef DEBUG_MODE
        if(debugMode>0)  
          cout << brownText << "Step j=" << j << " with " << clusters.getSize() 
               << " total clusters " << normText << endl;  // debugging
        #endif

        clusters[i].next();  // increment to next node for user iteration

        #ifdef DEBUG_MODE
        if(debugMode>0)  cout << magentaText << "j after  " << normText << endl;
        #endif

      } // end for j - particle loop
      
      #ifdef DEBUG_MODE
      if(debugMode>0)  
        cout << magentaText << "Step i = " << i << " with nClusters = " 
             << clusters.getSize() << normText << endl;  // debugging
      #endif

      // and get rid the cluster now if it is empty - difficult to make 
      // this occur internally within the class methods naturally
      // if we do not delete a list then increment i as normal, otherwise do 
      // not increment i since that would skip a cluster in the list
      if(clusters[i].getSize()==0 && !bqConstrained) {

        #ifdef DEBUG_MODE
        if(debugMode>0) 
          cout << redText << "deleting empty cluster " << i 
               << normText << endl;  // debugging
        #endif

        clusters.erase(i);
      } // end delete empty cluster
      else   i++;  // if no erase, increment as normal

      #ifdef DEBUG_MODE
      if(debugMode>0)  cout << blueText << "i after " << i << normText << endl;
      #endif

    } // end while i - cluster loop

    
    #ifdef DEBUG_MODE
    if(debugMode>0)  cout << redText << "after completed i loop " << flush; // debugging
    #endif

    if(nStep%CheckInterval==0) {
      // try a run truncation based on number of convergence checks
      if(lastEnergy==thisEnergy)  constantEnergyCount++;
      else {
        constantEnergyCount = 0;
        lastEnergy = thisEnergy;
        thisEnergy = lowEnergy;
      } // end else

      #ifdef DEBUG_MODE
      if(debugMode>0)  cout << "In Check a... " << flush;  // debugging
      // output run status - current total energy
      if(debugMode>1)  cout << "r " << flush;  // debugging
      #endif

      if(verbosity>0)  
        cout << "  E = " << clusters.calcEnergy(Cp) << " with Q = "
             << clusters.getQ()      << "... " << endl;

      #ifdef DEBUG_MODE
      if(debugMode>0)  cout << "In Check b... " << flush;  // debugging
      #endif

    } // end if nStep

    #ifdef DEBUG_MODE
    if(debugMode>0)  cout << brownText << "after check b " << flush; // debugging
    #endif

    // some more stuff such as splitting suspect clusters...
    converged = !moved;
    // if last x# CheckInterval moves had a zero energy change, set converged
    // flag anyhow (don't want perpetual zero changes)
    if(lastNodeMoved==thisNodeMoved && lastCluster==thisCluster) {
      if(verbosity>0)  
        warningMsg("Flagging convergence due to oscillating node!");
      converged = 1;
    } // end if equalEnergyCount
    else if(useEnergy) { // convergence checks that rely strictly on energy
      if(equalEnergyCount>=CheckInterval*3*(int)sqrt((double)NMax)) {
        if(verbosity>0)  
          warningMsg("Flagging convergence due to maxed zero moves!");
        converged = 1;
      } // enf if equalEnergyCount
      else if(constantEnergyCount==TryHard) {
        if(verbosity>0)
          warningMsg("Flagging convergence due to maximum constant energy checks!");
        converged = 1;
      } // end if constantEnergyCount
    } // end if useEnergy

    #ifdef DEBUG_MODE
    if(debugMode>0) 
      cout << brownText << "d with nStep = " << nStep << endl; // debugging
    #endif
/*
    if(!moved) 
      if(!useZeroEnergy) {
        useZeroEnergy = 1; // turn on zero energy moves
        converged = 0;     // unset converged status
        //cout << redText << "Switching on zero energy moves" << normText << "\n";
      }
*/
  } // ************* end while nStep *****************************************
  if(verbosity>1)
    cout << "\n" << blueText << "Done after step " << nStep << normText << endl;

  return nStep;
// ---------------------------------------------------------------------------
} // end communityDetect (original N^2 and NZn versions)
// ---------------------------------------------------------------------------

int communitySolve(TClusterList &clusters, TCMatrix &Cp, TCMData a, TCMData b,   
       bool useEnergy = 1,      int NTrialsMax = 1, int NStepMax = 10, 
       int  CheckInterval = 10, int verbosity  = 0, int debugMode = 0) {
  // Takes a cluster list clusters and performs NTrialsMax sets of complete 
  // solutions and sets clusters to the best answer at the end (if applicable)
  // Each individual solution is limited to at most NStepMax iterations though
  // the convergence criteria generally end the run before that value.
  // CheckInterval - limits the number of high cost convergence checks
  //                 only applies to original communityDetect
  // -------------------------------------------------------------------------
  // begin algorithm
  // -------------------------------------------------------------------------
  int  NEdges = Cp.getNEdges(), N = Cp.getSize();
  int  NStepStartCheck = 10;    // step number to begin periodic checks
  int  nStep, nRuns = 0;
  int  nodeOffset = 0;          // offset for displaying the answer
  int  QMax = 20;

  TVectorInt   trialEnergies(NTrialsMax);      // store energies of the trials
  TVectorFloat trialModularities(NTrialsMax);  // modularities of the trials
  TClusterList bestClusters(N,N,NEdges);       // current best clusters list

  int          bestTrial = 1;     // which trial is best? (E or Q)
  int          bestStep = 1;      // what is number of steps for best trial
  trialEnergies.init(2 << 29);    // store energies of the trials
  trialModularities.init(-1.0);   // modularities of the trials
/*
  if(verbosity>0) {
    cout << green << "Starting solution using ";
    if(bDense) cout << brown << "dense";  else cout << brown << "sparse"; 
    cout << green << " cluster iteration:"  << normText << "\n";
  } // end if verbosity
*/  
  #ifdef MSL_TCMWEIGHTED
  #else
  //if(verbosity>0)  
  //  warningMsg("Using a dense connection matrix.  Assuming non-NZ algorithm!");
  #endif

  int  nTrials = 0;
  while(nTrials<NTrialsMax) {
    nTrials++;  // start trial loop with nTrials = 1
    
    if(verbosity>0 && NTrialsMax>1)  
      cout << green << "Starting trial " << brown << nTrials << ":  " 
           << normText << "\n";
    
    // -------------------------------------------------------------------------
    // initialize clusters
    // -------------------------------------------------------------------------
    //clusters.initRandom(N,QMax);         // initialize a random configuration
    //cout << green << "before init... " << normText << flush;
    clusters.initSymmetric(1);             // re-initialize nodes randomize 
    clusters.initClean();                  // initialize all to clean (unmoved)
    //cout << green << "after init... " << normText << flush;
    //cout << "Here inside solve B  " << endl; // debugging
    if(verbosity>1)  clusters.display(1,0,"Init state inside clusterSolve ");

    #ifdef DEBUG_MODE
    if(debugMode>0) {
      clusters.display(nodeOffset);           // display random initial state
      cout << magenta << "Initial energy = " << clusters.getEnergy()
           << " and mod = " << clusters.getQ() << normText << "\n";
    } // end if debugMode
    #endif

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    #ifdef DEBUG_MODE
    cout << normText << "In communitySolve " << blue << "before" << normText 
         << " community detect... " << endl; // debugging
    #endif

    // run the actual community detection algorithm for the unknown system
    // Note that if we use a dense connection matrix we assume the old algorithm
    // so that we can account for weighted negative edges.  This is reasonable
    // since the dense matrix is limited to 10,000 or so nodes anyhow.
    // If there is a need for an NZ algorithm that accounts for negative 
    // weighted edges, we can implement it at a later time. 
    #ifdef MSL_TCMDENSE
    // use old version
    nStep = communityDetect(clusters,Cp,NStepMax,CheckInterval,1,1,3,verbosity);
    #else
    nStep = communityDetectNZOpt(clusters,Cp,a,b,NStepMax,0,0,verbosity);
    #endif
    //cout << green << "before solve... " << normText << flush;
    //cout << green << "after solve" << normText << endl;
    clusters.calcEnergy(Cp);  // supposed to be a redundant calculation - not verified
    #ifdef DEBUG_MODE
    cout << normText << "In communitySolve " << magenta << "after" << normText 
         << " community detect... " << endl; // debugging
    #endif

    //bestClusters.display(0,0,"best ");    cout << endl; // debugging
    //clusters.display();    cout << endl; // debugging

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // keep track of best for several trial cases
    //cout << magenta << "useEnergy = " << useEnergy << endl;  // debugging
    if(NTrialsMax>1) { // only store if we have more than one trial attempted
      if(nTrials==1) {
        //bestClusters.display(nodeOffset,0,"Best first trial before ");
        //cout << "Before copy with " << clusters.getSize() << " clusters... " 
        //     << endl; // debugging
        bestClusters = clusters;  // start with first as default best
        //cout << "After copy... " << endl; // debugging
        bestTrial = 1;  bestStep = nStep;
        //bestClusters.display(nodeOffset,0,"Best first trial after ");
      } else {
        if(useEnergy) {
          if(clusters.getEnergy()<bestClusters.getEnergy()) {
            //errorMsg("Yup?  We are here?");  // debugging
            //bestClusters.display(nodeOffset,0,"Best trial before ");
            bestClusters = clusters;
            bestTrial = nTrials;  bestStep = nStep;
            //bestClusters.display(nodeOffset,0,"Best trial after ");
          } // end if on energy test
        } else { // use a modularity criterion instead
            if(clusters.getQ() > bestClusters.getQ()) {
              //errorMsg("Huh?  Why are we here?");  // debugging
              bestClusters = clusters;
              bestTrial = nTrials;  bestStep = nStep;
            } // end if on modularity test
          } // end else useEnergy
      } // end else nTrials
    } // end if NTrialsMax
    else  bestStep = nStep;

    //clusters.display();  // debugging
    //errorMsg("Done");    // debugging


    trialEnergies[nTrials-1]     = clusters.getEnergy();      // indexed from 0
    trialModularities[nTrials-1] = clusters.getQ();  // indexed from 0
  } // end for nTrials
  // end trials - output run information
  // -------------------------------------------------------------------------
  // these two 'best' values are based on the chosen criterion for convergence
  // if we use energy, modularity is not maximized necessarily
  int     bestEnergy     = bestClusters.getEnergy();
  double  bestModularity = bestClusters.getQ();

  // some verbose output
  if(verbosity>1) {
    bestClusters.display(nodeOffset,0,"Best trial ");
    // output all energies and modularities of the trials
  } // end if
  if(verbosity>0) {
    cout << green << "\nEnergies for the " << NTrialsMax << " trials were:  "
         << grey  << "(Best trial was " << green << bestTrial << grey << "):  ";
    TVectorOut(cout,trialEnergies,1,1,grey);  cout << "\n";
    cout << "Modularities for the " << NTrialsMax << " trials were:  ";
    TVectorOut(cout,trialModularities,1,1,grey);
    cout << magenta << "\nFinal energy = " << bestClusters.getEnergy() << "  ";
    cout << brown   << "Done on step "   << nStep << ".  There were " 
         << clusters.getSize() << " final clusters." << normText << "\n";
  } // end if

  // now assign the best answer if we do more than one trial and it is not 
  // already the last trial (would be a redundant copy)
  if(NTrialsMax>1 && bestTrial!=NTrialsMax)  clusters = bestClusters;  

  return bestStep; // return # of steps required for convergence of bestCluster
} // end communitySolve
//-----------------------------------------------------------------------------


inline int subcommunitySolve(TClusterList &a, TCluster &c, TCMatrix &Cp, 
              int NTrialsMax = 10, int NStepMax = 1000, int CheckInterval = 100, 
              bool useDense = 0) {
  // take a single cluster and solve it independently as a complete system
  // takes a cluster 'c' and returns a sub-divided clusterList 'a'
  a.erase();  // erase entire current list if it is not empty
  // initialize answer list with current cluster - one to one initialization
  a.add(c);
  a.initSymmetric(1);  // initSymmetric just distributes existing nodes evenly
  // now solve the system via the regular community solve routine
  int nStep;
  nStep = communitySolve(a,Cp,NTrialsMax,NStepMax,CheckInterval,useDense);
  return nStep;  // return # of steps required for convergence of best solution
} // end subcommunitySolve


/* ***************************************************************************
  dump CMatrix values
*************************************************************************** */
void coutCMatrix(TMatrixInt &CM, int bOutputLarge=0, bool bShowZeroes=1) {
  int i, j, cij;
  int N = CM.getCols();

  //cout << blueText << " CM " << CMatrix[1][0] << normText;
  if(N>96 || !bShowZeroes)  return;  // debugging
  //if(N>64 && !(bool)bOutputLarge)  return;  // debugging
  //if(bOutputLarge < 2) return;  // clearFlags()debugging

  if(bOutputLarge>0) { // scan for largest integer
    int maxInt = 0;
    for(i=0; i<N; i++)  for(j=0; j<N; j++) 
      if(abs(CM(i,j))>maxInt)  maxInt = abs(CM(i,j));
    if(maxInt==0) { 
      maxInt = 1; // dummy value for zero matrix
    //  warningMsg("Matrix is a zero matrix!  Exiting matrix output.");  
    //  return; 
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
        cij = CM(j,i);
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
      if(abs(CM(i,j))>maxInt)  maxInt = CM(i,j);
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
        cij = CM(j,i);
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
} // end coutCMatrix


