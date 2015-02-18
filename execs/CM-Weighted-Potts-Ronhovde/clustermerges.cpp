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
// ---------------------------------------------------------------------------
// community utility functions functions
// ---------------------------------------------------------------------------

// ------------------------- TMergeNode class --------------------------------
// TMergeNode used for the NZ merge routines to track and sort all merges on
// only a few NZ passes through the list.  Susequent passes resolve conflict 
// of merge clusters.
// Merges are performed in the order of the largest energy decrease, and the
// easiest way to track and sort the merges is to use a simple encapsulated 
// data structure to be used with the TList class.
// We define the corresponding boolean operators for the sort routine.
class TMergeNode {
 public:
  // overloaded operators
  inline TMergeNode& operator=(const TMergeNode& b);
  // the boolean operators are defined on the mergeE rather than the weight
  // because of the binary search in the reference TCMSparse::operator(i,j)
  inline bool operator==(const TMergeNode& b);
  inline bool operator> (const TMergeNode& b);
  inline bool operator< (const TMergeNode& b);
  inline bool operator>=(const TMergeNode& b);
  inline bool operator<=(const TMergeNode& b);
  // the operator functions are needed mostly to maintain basic compatibility
  // with TList operations such as display(...).  They are not well-defined 
  // for the TMergeNode.  They only change the energy value not the clusters.
  inline TMergeNode operator+(const TMergeNode& b);  // for energy only!
  inline TMergeNode operator+(int    i);             // for energy only!
  inline TMergeNode operator+(double d);             // for energy only!
  // Standard public member functions
  inline      TMergeNode(int cluster1=-1, int cluster2=-1, TCMData mergeVal=0);
  inline      TMergeNode(TMergeNode &b);                // copy constructor
  //         ~TMergeNode() {};                           // destructor
  // keep the data public for the implementation classes
         int      c1, c2;
         TCMData  mergeE;
}; // end class TMergeNode

// define the TMergeNode methods on functions
inline TMergeNode::TMergeNode(int cluster1, int cluster2, TCMData mergeVal) { 
  c1 = cluster1;  c2 = cluster2;  mergeE = mergeVal;};     // constructor
inline TMergeNode::TMergeNode(TMergeNode &b) { 
  c1 = b.c1;  c2 = b.c2;  mergeE = b.mergeE;};     // copy constructor
// define operators
inline bool TMergeNode::operator==(const TMergeNode& b) { 
  return mergeE==b.mergeE && 
         ( (c1==b.c1 && c2==b.c2) || (c1==b.c2 && c2==b.c1) ); }; 
inline bool TMergeNode::operator> (const TMergeNode& b) { 
  return ( mergeE>b.mergeE || ( mergeE==b.mergeE && c1>b.c1 ) ||
                              ( mergeE==b.mergeE && c1==b.c1 && c2>b.c2 )); };
inline bool TMergeNode::operator< (const TMergeNode& b) { 
  return ( mergeE<b.mergeE || ( mergeE==b.mergeE && c1<b.c1 ) ||
                              ( mergeE==b.mergeE && c1==b.c1 && c2<b.c2 )); };
inline bool TMergeNode::operator>=(const TMergeNode& b) { 
  return ( (*this)>b || (*this)==b ); }; // use already defined operators
inline bool TMergeNode::operator<=(const TMergeNode& b) { 
  return ( (*this)<b || (*this)==b ); }; // use already defined operators
inline TMergeNode& TMergeNode::operator=(const TMergeNode& b) { 
  c1 = b.c1;  c2 = b.c2;  mergeE = b.mergeE;  return *this; }; // op '='
// the operator functions are needed mostly to maintain basic compatibility
// with TList operations such as display(...).  They are not well-defined 
// for the TMergeNode.  They only change the energy value not the clusters.
inline TMergeNode TMergeNode::operator+(const TMergeNode& b) { 
  TMergeNode c = (*this);  c.mergeE += b.mergeE;    return c; };
inline TMergeNode TMergeNode::operator+(int i) { 
  TMergeNode c = (*this);  c.mergeE += (TCMData)i;  return c; };
inline TMergeNode TMergeNode::operator+(double i) { 
  TMergeNode c = (*this);  c.mergeE += (TCMData)i;  return c; };
// additional TMergeNode functions
inline ostream& operator<<(ostream& fout, TMergeNode& a) {
  // Outputs the weight of the TCMData node only
  fout << "(" << a.c1 << "," << a.c2 << "," << a.mergeE << ")";  return fout; };
// ------------------------ End TMergeNode class ------------------------------

// ---------------------- Merge detection functions ---------------------------
int  getAllMergesNZ(TClusterList &a, TCMatrix &Cp, 
       int &minCluster1, int &minCluster2, bool &bAllMerge, 
       TList<TMergeNode> &merges, TCMData gammaa = 1, TCMData gammab = -1) {
  // The "All" version tracks all possible merges on one O(NZ) pass through the 
  // neighbor list.  We store the lowest energy change for each cluster.  
  // These are then resolved in almost one pass, up to conflicts, by the
  // detectAndMerge"All" function below.  Conflicts occur when two clusters
  // would merge with another cluster.  We cannot merge both because the 
  // connections change; and due to memory constraints, the NZ operation only 
  // tracks connections for the current cluster (although we may be able to 
  // implement a secondary tracking on the respective number of edges).
  // The minimum merge is still manually returned (for the moment).
  // Assumes Cp is one of the sparse matrix representations.
  // We avoid O(q) scans of the energy and other vectors in order to avoid
  // an overall O(q^2) algorithm.  The other merge algorithm is already O(q^2).
  // The idea is that this version is used in 'very' sparse systems where q 
  // is more comparable to N than to Z or n.
  // The function returns the minimum energy detected over the scan (positive 
  // or negative), but does not perform the actual merge; therefore, q is 
  // constant over the operation of this function.
  
  #ifdef DEBUG_MODE
  //cout << "Starting getAllMergesNZ... " << flush;
  #endif

  if(a.getSize()==1) {
    minCluster1 = -1;  minCluster2 = -1;  // invalidate merged clusters
    return 0;
  } // end if size check

  int   iCluster, jCluster, iNode, jNode, jDegrees;
  int   N = a.getNNodes();  // q is constant over the operation of this function
  int   q = a.getSize();    // N is the total number of nodes
  //bool  bAllMerge = 1;
  int   minMergeEnergy = BigInteger;
  int  *pNeighborStart, *pNeighbor; // optimzing neighbor pointers
  int  *pjClustereSum;      // optimzing pointer for energy sum
  // track the minima for the current iCluster
  TMergeNode  m;  // the temporary node
  merges.erase();             // delete any existing data
  int nPossibleMerges = 0;

  TVectorInt  eSum(q,0);      // energy sum for each cluster for current cluster
  TVectorInt  bcChecked(q,0); // was cluster 'x' checked already?
  TVectorInt  edgeCount(q,0); // count number of nodes connected to a cluster

  //a.calcEnergy(Cp);  // just in case
  //a.initNodesList(); // the node to cluster identifier list

  // loop over particles and move each one in turn
  int i = 0;
  #ifdef DEBUG_MODE
  //cout << "Start cluster loop... " << flush;
  #endif
  while(i<q) {
    a[i].begin();  // begin user iteration
    iCluster = i;  // just to keep clear code

    m.mergeE = BigInteger; // reset temporary merge
    m.c1 = -1;  m.c2 = -1;  
    // ----------------------------------------------------------------
    for(int k=0; k<a[i].getSize(); k++) {
      //wNode = a[i].getCurrent();
      iNode = a[i].getCurrentNode();

      //cout << blue << k << ", " << normText << iNode << "| " << flush; // debugging

      // scan edge connection list for iNode and sum energy contributions
      jDegrees       = Cp.ki(iNode);
      pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7% speed boost
      pNeighbor      = pNeighborStart;
      for(int j=0; j<jDegrees; j++) {
        //jNode = Cp.ki(iNode,j);
        //jCluster = a.nodes[jNode];
        jNode = (*pNeighbor);
        jCluster = a.nodes[jNode];
        // count number of nodes jNode is connected to in jCluster
        edgeCount[jCluster] += 1;
        // sum (subtract) the energy for node iNode merge into jCluster for
        // the unweighted version, we can skip the energy sum and just use
        // the count with the stored and unstored values
        // I use pre-processor if's to avoid an extra inner loop if() call
        #ifdef MSL_TCMWEIGHTED
        eSum[jCluster] -= Cp(iNode,jNode);
        #else
        eSum[jCluster] -= gammaa;  // avoids an O(log Z) function call
        #endif

        pNeighbor++;  // increment neighbor pointer
      } // end for j - neighbor list search
      #ifdef DEBUG_MODE
      //cout << "Done with neighbor loop k... " << normText << flush;
      #endif

      a[i].next();  // increment to next node for user iteration
    } // end for k - particle loop
    #ifdef DEBUG_MODE
    //cout << "Done with node loop j... " << normText << flush;
    #endif
    // -----------------------------------------------------------------

    // now that we have looped over all nodes in cluster i, we now have to 
    // finish calculating the merge energies by reiterating over the neighbors
    a[i].begin();  // renew user iteration
    // ----------------------------------------------------------------
    for(int k=0; k<a[i].getSize(); k++) {
      iNode = a[i].getCurrentNode();

      // scan edge connection list and update energies for unconnected edges
      // first, set the energy change for removing iNode from iCluster
      // which is slightly different from the added case
      eSum[iCluster] = 0;       // cannot merge a cluster with itself, so 
      bcChecked[iCluster] = 1;  // skip iCluster in following energy test
      jDegrees       = Cp.ki(iNode);
      pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7% speed boost
      pNeighbor      = pNeighborStart;
      for(int j=0; j<jDegrees; j++) {
        //jNode = Cp.ki(iNode,j);
        jNode = (*pNeighbor);
        jCluster = a.nodes[jNode];
        pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
        if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
          // sum (subtract) the energy for node iNode merge into jCluster
          // note that the unstored value is a negative number
          (*pjClustereSum) += 
            (a[iCluster].getSize()*a[jCluster].getSize()
            -edgeCount[jCluster])*(-gammab);
          bcChecked[jCluster] = 1;
          
          // keep track of whether we merge all clusters in a
          // a problem is that we do not necessarily scan all cluster
          // combinations such as when disconnected cluster list exists!
          bAllMerge &= ( (*pjClustereSum)<0 );

          // include the lowest energy change but skip the current cluster
          // iCluster==jCluster is excluded by above if statement
          //if(iCluster!=jCluster && (*pjClustereSum)<minMergeEnergy) {
          // this iCluster minimum energy test
          if((*pjClustereSum)<0)  nPossibleMerges++;
          if((*pjClustereSum)<m.mergeE) {
            m.mergeE = (*pjClustereSum);
            m.c1 = min(iCluster,jCluster);  m.c2 = max(iCluster,jCluster);
          } // end if eSum
          // global minimum energy test
          if((*pjClustereSum)<minMergeEnergy) {
            minMergeEnergy = (*pjClustereSum);
            minCluster1 = iCluster;  minCluster2 = jCluster;
          } // end if eSum
        } // end if bcChecked

        pNeighbor++;  // increment neighbor pointer
      } // end for j - neighbor list search
      #ifdef DEBUG_MODE
      //cout << "Done with last node loop. " << normText << flush;
      #endif

      // scan edge connection list and clear changed values - assignments
      // can be redundant
      pNeighbor = pNeighborStart;  // optimizing neighbor node pointer
      for(int j=0; j<jDegrees; j++) {
        //jNode = Cp.ki(iNode,j);
        jNode = (*pNeighbor);
        jCluster = a.nodes[jNode];
        eSum[jCluster]      = 0;
        edgeCount[jCluster] = 0;
        bcChecked[jCluster] = 0;

        pNeighbor++;  // increment neighbor pointer
      } // end for j - neighbor list search

      a[i].next();  // increment to next node for user iteration
    } // end for k - particle loop
    // -----------------------------------------------------------------
    #ifdef DEBUG_MODE
    //cout << "Done with last node loop. " << normText << flush;
    #endif

    // now store the minimum energy for this iCluster
    if(m.mergeE<0)  merges.add(m);

    i++;
  } // end while i - cluster loop
  #ifdef DEBUG_MODE
  //cout << "Done with cluster loop. " << normText << flush;
  #endif

  // now sort the possible merges by energy (embedded in TMergeNode class)
  merges.sort();
  #ifdef DEBUG_MODE
  //cout << "Done with merge list sort... " << normText << flush;
  #endif

  // remove duplicate merges
  int k=0, eraseCount = 0;
  if(merges.getSize()>1) {
    #ifdef DEBUG_MODE
    //cout << "Starting " << flush;
    #endif
    merges.begin();
    TMergeNode *pn = &(merges.getStart());  // the previous data
    merges.next();
    #ifdef DEBUG_MODE
    //cout << "merge iteration on a size " << flush;
    //cout << merges.getSize() << " list... " << flush;
    #endif
    while(k<merges.getSize()-1) {
      //cout << "Merge nodes are: " << merges.getCurrent() << " and " 
      //     << (*pn) << endl;  // debugging
      if(merges.getCurrent() == (*pn)) {
        #ifdef DEBUG_MODE
        //cout << "Erasing current " << merges.getCurrent() << "... " 
        //     << endl;  // debugging
        #endif
        merges.eraseCurrent();
        eraseCount++;
      } // end if
      // store current as new previous before iterating unless we deleted 
      // current, then iPrevious is unchanged
      else  pn = &(merges.getCurrent());// update previous TMergeNode to current

      merges.next();
      k++;
    } // end while k
  } // end if merge size check
  #ifdef DEBUG_MODE
  //cout << "Done with duplicate merge removal... " << normText << flush;
  //if(((double)nPossibleMerges/((double)q*(double)(q-1)) )>0.75)
  //  cout << red << "Warning:  NPossibleMerges = " << nPossibleMerges 
  //       << " out of " << (q*(q-1)) << "... " << normText << flush; // debugging
  #endif

  return minMergeEnergy;
} // end getAllMergesNZ


int  detectAndMergeAllNZ(TClusterList &c, TCMatrix &Cp, 
                         TCMData a=1, TCMData b=-1) {
  // A Z^2 n^2 log Z version, we loop over the a and pick out the max 
  // energy decrease.  We then must repeat the entire process until no merges  
  // are detected.  Assumes Cp is one of the sparse matrix representations.
  
  // We cannot shortcut the process if all clusters in c need to be merged
  // as we do with the full cluster pair, O(q^2) or O(NZ log Z), version.
  // The neighbor search cannot detect a full merge because disconnected
  // clusters in c will not be sampled during a neighbor search (check if we
  // can add an O(NZ) pre-scan to try to detect disconnected c???)
  
  if(c.getSize()==1) {
    //cout << "Trivial merge attempt case.  Size is one.\n";  // debugging
    return 0;
  } // end if trivial case

  int   minMergeE = BigInteger, mergeE;
  int   nMerges = 0, nMergesNZ;  // number of merges total and for current loop

  int   maxq, minq, minC1, minC2;  // cluster references
  int   NMergesMax = c.getSize()-1;
  //int   N = c.getNNodes();  // q is constant over the operation of this function
  int   deltaE;
  
  bool  bAllMerge;
  TList<TMergeNode>  merges;
  TVectorInt         mergedClusters(c.getSize());

  //cout << red << "In mergeAllNZ\n" << normText << flush;         // debugging
  if(b>0)  errorMsg("b is positive?");  // debugging
  
  int mergeStep = 0;
  do {
    //cout << "Starting NZ pass " << mergeStep << "... " << flush;  // debugging

    bAllMerge = 0;                   // can the system fully merge?
    mergedClusters.init(0);          // reset the merged clusters list
    mergedClusters.resize(0,0,0,1);  // fake resize (no change in array length)
    //cout << "before get " << mergeStep << "... " << flush;   // debugging
    minMergeE = getAllMergesNZ(c,Cp,minC1,minC2,bAllMerge,merges,a,b);
    //cout << "after get " << mergeStep << "... " << flush;   // debugging
    //if(merges.getSize()>50)  // debugging
    //  cout << "Potential merges are: " << merges << endl;  // debugging
    //errorMsg("Done");        // debugging

    // we don't want to shortcut with a return statement since we update
    // cleared clusters at the end of the function
    //if(minMergeE>=0)  nMerges = 0;   // no possible merges
    if(bAllMerge) {                  // we have a complete merge
      for(int i=NMergesMax-1; i<0; i++)  c.move(i,0,1);  // move i to cluster 0
      // can't calculate energy with c.calcEnergy(.), or we have an O(N^2) 
      // calculation due to collapsed system
      // operate on some class specific variables - perhaps should be a method
      c.energyTotal = -c.cEdgeSum + c.uEdgeSum;
      c.modularity  = 0;
      // now update the cluster data
      c[0].energy   = c.energyTotal;  c[0].partialQ = 0;
      c[0].nNodes   = c.NNodes;       c[0].nEdges   = c.NEdges;
      c[0].cEdgeSum = c.cEdgeSum;     c[0].uEdgeSum = c.uEdgeSum;
      
      nMerges = NMergesMax;
    } // end if bAllMerge

    nMergesNZ = 0;
    int i = 0;
    // first merge is automatically minMergeE - also exits the loop if necessary
    mergeE = minMergeE;  
    if(merges.getSize()>0)  merges.begin();  // begin manual TList iteration
    while(i<merges.getSize() && mergeE<0) {
      //cout << "Starting iteration " << i << "... " << flush;       // debugging
      // at some point, we can optimize this for multiple merges by restricting 
      // erasures until the end of the function
      mergeE = merges.getCurrent().mergeE;  // get energy for this merge
      //localInitTime();
      //timeMsg("Merge time before ","...");  // debugging
      //cout << blue << "Attempting to merge clusters (list size " << c.getSize()
      //     << ") " << magenta << merges.getCurrent().c1 << blue << " and "
      //    << magenta << merges.getCurrent().c2 << blue << " at an energy of ";
      //if(mergeE<0) cout << red;  else cout << green;
      //cout << mergeE << normText << endl;  // debugging

      if(mergeE<0) {
        #ifdef DEBUG_MODE
        //cout << blue << "Attempting to merge clusters (size " << c.getSize()
        //     << ") " << magenta << merges.getCurrent().c1 << blue << " and "
        //     << magenta << merges.getCurrent().c2 << blue << " at an energy of ";
        //if(mergeE<0) cout << red;  else cout << normText;
        //cout << mergeE << normText << endl;  // debugging
        #endif

        // merge the a and erase maximum index cluster - use a fast move
        //c.display(0,0,"Pre-move "+itos(nMerges)+" ");            // debugging
        if(merges.getCurrent().c1 > merges.getCurrent().c2) {
          maxq = merges.getCurrent().c1;  minq = merges.getCurrent().c2; 
        } else {
          maxq = merges.getCurrent().c2;  minq = merges.getCurrent().c1; 
        } // end else

        // move list (and updates nodes list internally) if it is not in 
        // conflict with an already merged cluster.  That is, has either cluster
        // already been involved in a merge during this getALLMergesNZ call.
        // The current isMember() call in TVectorInt uses a linear search.
        if( !(mergedClusters.isMember(maxq) || mergedClusters.isMember(minq)) ){
          deltaE = c[maxq].getEnergy() + mergeE;
          // In this version, because we track the clusters that were merged, 
          // we cannot erase the empty lists until the end of the whole process.
          //c.move(maxq,minq,deltaE,1);
          c.updateNodesList(maxq,minq);  // update nodes list with known change
          c[maxq].move(c[minq],deltaE);
          // now record the clusters that were merged
          mergedClusters.add(maxq);  mergedClusters.add(minq);
          nMergesNZ++;
        } // end if member check
        //else cout << "Skipping conflicting merge" << endl;       // debugging
        //c.display(0,0,"Move "+itos(nMerges)+" ");                // debugging

      } // end if
      //cout << "Done with iteration " << flush;         // debugging
      //else
      //  if(minMergeEnergy==0 && minCluster1>=0 && minCluster2>=0)
      //    cout << blue << "a " << minCluster1 << " and " << minCluster2
      //         << " have a zero merge energy." << normText << "\n";// debugging
      
      //timeMsg("after ","\n");  // debugging
      merges.next();  // next in manual TList iteration
      //cout << i << endl;         // debugging
      i++;
    } // end while i loop

    nMerges += nMergesNZ;
    //cout << "Done with NZ pass " << mergeStep << endl;           // debugging
    mergeStep++;
    // we test for a max bound, a negative minimum energy and check whether
    // this iteration exhausts all possible merges (i.e. the were no conflicts)
    // Actually, the nMergesNZ check could be inconsistent since we pick the 
    // largest decrease for each node.  Need a final check to verify that all 
    // possible merges have occured.
  //} while(mergeStep<NMergesMax && minMergeE<0 && nMergesNZ<merges.getSize()); 
  } while(mergeStep<NMergesMax && minMergeE<0); 
  // end merge test do loop

  // now we can finally (and efficiently) erase all moved clusters  
  int nErased;  nErased = c.clearZeroClusters(1);
  #ifdef DEBUG_MODE
  //cout << magenta << "Merged " << nMerges << " clusters with " 
  //     << mergeStep << " iterations... " << normText << flush;  // debugging
  #endif
    
  return nMerges;
} // end detectAndMergeAllNZ


int  getBestMergeNZ(TClusterList &a, TCMatrix &Cp, 
       int &minCluster1, int &minCluster2, TCMData gammaa = 1, TCMData gammab = -1) {
  // A NZ log Z version, we loop over the a and pick out the max 
  // energy decrease.  Assumes Cp is one of the sparse matrix representations.
  // We avoid O(q) scans of the energy and other vectors in order to avoid
  // an overall O(q^2) algorithm.  The other merge algorithm is already O(q^2).
  // The idea is that this version is used in 'very' sparse systems where q 
  // is more comparable to N than to Z or n.
  // The function returns the minimum energy detected over the scan (positive 
  // or negative), but does not perform the actual merge; therefore, q is 
  // constant over the operation of this function.
  
  if(a.getSize()==1) {
    minCluster1 = -1;  minCluster2 = -1;  // invalidate merged clusters
    return 0;
  } // end if size check

  int   iCluster, jCluster, iNode, jNode, jDegrees;
  int   N = a.getNNodes();  // q is constant over the operation of this function
  int   q = a.getSize();    // N is the total number of nodes
  bool  bAllMerge = 1;
  int   minMergeEnergy = BigInteger;
  int  *pNeighborStart, *pNeighbor; // optimzing neighbor pointers
  int  *pjClustereSum;      // optimzing pointer for energy sum

  TVectorInt  eSum(q,0);      // energy sum for each cluster for current cluster
  TVectorInt  bcChecked(q,0); // was cluster 'x' checked already?
  TVectorInt  edgeCount(q,0); // count number of nodes connected to a cluster

  //a.calcEnergy(Cp);  // just in case
  //a.initNodesList(); // the node to cluster identifier list

    // loop over particles and move each one in turn
    int i = 0;
    while(i<q) {
      a[i].begin();  // begin user iteration
      iCluster = i;  // just to keep clear code

      // ----------------------------------------------------------------
      for(int k=0; k<a[i].getSize(); k++) {
        //wNode = a[i].getCurrent();
        iNode = a[i].getCurrentNode();

          // scan edge connection list for iNode and sum energy contributions
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            //jCluster = a.nodes[jNode];
            jNode = (*pNeighbor);
            jCluster = a.nodes[jNode];
            // count number of nodes jNode is connected to in jCluster
            edgeCount[jCluster] += 1;
            // sum (subtract) the energy for node iNode merge into jCluster for
            // the unweighted version, we can skip the energy sum and just use
            // the count with the stored and unstored values
            // I use pre-processor if's to avoid an extra inner loop if() call
            #ifdef MSL_TCMWEIGHTED
            eSum[jCluster] -= Cp(iNode,jNode);
            #else
            eSum[jCluster] -= gammaa;  // avoids an O(log Z) function call
            #endif

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

        a[i].next();  // increment to next node for user iteration
      } // end for k - particle loop
      // -----------------------------------------------------------------

      // now that we have looped over all nodes in cluster i, we now have to 
      // finish calculating the merge energies by reiterating over the neighbors
      a[i].begin();  // renew user iteration
      // ----------------------------------------------------------------
      for(int k=0; k<a[i].getSize(); k++) {
        iNode = a[i].getCurrentNode();

          // scan edge connection list and update energies for unconnected edges
          // first, set the energy change for removing iNode from iCluster
          // which is slightly different from the added case
          eSum[iCluster] = 0;       // cannot merge a cluster with itself, so 
          bcChecked[iCluster] = 1;  // skip iCluster in following energy test
          jDegrees       = Cp.ki(iNode);
          pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7% speed boost
          pNeighbor      = pNeighborStart;
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = a.nodes[jNode];
            pjClustereSum = &(eSum[jCluster]);  // optimizing pointer
            if(bcChecked[jCluster]==0) { // iCluster is omitted with bcChecked
              // sum (subtract) the energy for node iNode merge into jCluster
              // note that the unstored value is a negative number
              (*pjClustereSum) += 
                (a[iCluster].getSize()*a[jCluster].getSize()
                -edgeCount[jCluster])*(-gammab);
              bcChecked[jCluster] = 1;
              
              // keep track of whether we merge all clusters in a
              // a problem is that we do not necessarily scan all cluster
              // combinations such as when disconnected cluster list exists!
              bAllMerge &= ( (*pjClustereSum)<0 );

              // find the lowest energy change but skip the current cluster
              // iCluster==jCluster is excluded by above if statement
              //if(iCluster!=jCluster && (*pjClustereSum)<minMergeEnergy) {
              if((*pjClustereSum)<minMergeEnergy) {
                minMergeEnergy = (*pjClustereSum);
                minCluster1    = iCluster; // because we test over all a
                minCluster2    = jCluster;
              } // end if eSum
            } // end if bcChecked

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

          // scan edge connection list and clear changed values - assignments
          // can be redundant
          pNeighbor   = pNeighborStart;  // optimizing neighbor node pointer
          for(int j=0; j<jDegrees; j++) {
            //jNode = Cp.ki(iNode,j);
            jNode = (*pNeighbor);
            jCluster = a.nodes[jNode];
            eSum[jCluster]      = 0;
            edgeCount[jCluster] = 0;
            bcChecked[jCluster] = 0;

            pNeighbor++;  // increment neighbor pointer
          } // end for j - neighbor list search

        a[i].next();  // increment to next node for user iteration
      } // end for k - particle loop
      // -----------------------------------------------------------------

      i++;
    } // end while i - cluster loop

  //if(minMergeEnergy==BigInteger)  return 0; // remaining clusters are not connected
  return minMergeEnergy;
} // end getBestMergeNZ


int  detectAndMergeNZ(TClusterList &c, TCMatrix &Cp, 
                      TCMData gammaa=1, TCMData gammab=-1){
  // A Z^2 n^2 log Z version, we loop over the a and pick out the max 
  // energy decrease.  We then must repeat the entire process until no merges  
  // are detected.  Assumes Cp is one of the sparse matrix representations.
  
  // We cannot shortcut the process if all clusters in c need to be merged
  // as we do with the full cluster pair, O(q^2) or O(NZ log Z), version.
  // The neighbor search cannot detect a full merge because disconnected
  // clusters in c will not be sampled during a neighbor search (check if we
  // can add an O(NZ) pre-scan to try to detect disconnected c???)
  
  if(c.getSize()==1) {
    //cout << "Trivial merge attempt case.  Size is one.\n";  // debugging
    return 0;
  } // end if trivial case

  //bool  bAllMerge = 0;
  int   minMergeEnergy = BigInteger, nMerges = 0;
  int   maxq, minq, minCluster1, minCluster2;
  int   NStepMax = c.getSize();
  //int   N = c.getNNodes();  // q is constant over the operation of this function
  int   deltaE;

  //cout << red << "In mergeNZ\n" << normText << flush;         // debugging
  if(gammab>0)  errorMsg("b is positive?");  // debugging
  
  int i = 0;
  do {
    // at some point, we can optimize this for multiple merges by restricting 
    // erasures until the end of the function
    //localInitTime();
    minMergeEnergy = getBestMergeNZ(c,Cp,minCluster1,minCluster2,gammaa,gammab);
    //timeMsg("Merge time before ","...");  // debugging

    cout << blue << "Attempting to merge clusters (list size " << c.getSize()
         << ") " << magenta << minCluster1 << blue << " and "
         << magenta << minCluster2 << blue << " at an energy of ";
         if(minMergeEnergy<0) cout << red;  else cout << green;
    cout << minMergeEnergy << normText << endl;  // debugging

    if(minMergeEnergy<0) {
      #ifdef DEBUG_MODE
      cout << blue << "Attempting to merge clusters (list size " << c.getSize()
           << ") " << magenta << minCluster1 << blue << " and "
           << magenta << minCluster2 << blue << " at an energy of ";
           if(minMergeEnergy<0) cout << red;  else cout << normText;
      cout << minMergeEnergy << normText << endl;  // debugging
      #endif

      // merge the a and erase maximum index cluster - use a fast move

      //c.display(0,0,"Pre-move "+itos(nMerges)+" ");  // debugging
      if(minCluster1>minCluster2) { maxq = minCluster1;  minq = minCluster2; }
      else                        { maxq = minCluster2;  minq = minCluster1; }
      deltaE = c[maxq].getEnergy() + minMergeEnergy;
      //deltaE = c[minCluster2].getEnergy() + minMergeEnergy;

      // move list (and updates nodes list internally)
      c.move(maxq,minq,deltaE,1);

      //c.display(0,0,"Move "+itos(nMerges)+" ");  // debugging

      nMerges++;
    } // end if
    //else if(minMergeEnergy==0 && minCluster1>=0 && minCluster2>=0)
    //       cout << blue << "a " << minCluster1 << " and " << minCluster2
    //            << " have a zero merge energy." << normText << "\n";// debugging
    
    //timeMsg("after ","\n");  // debugging
    i++;
  } while(i<=NStepMax && minMergeEnergy<0); // end do loop
  
  //cout << magenta << "Merged " << nMerges << " total clusters " 
  //     << normText << "\n";    // debugging

  return nMerges;
} // end detectAndMergeNZ


int  getListMergesNZ(TClusterList &a, TCMatrix &Cp, TMatrixInt &mergeMatrix, 
                     TMatrixInt &countMatrix, bool bDisplayMerges=0) {
  // this routine assumes Cp is one of the sparse matrix representations
  int iCluster, jCluster, jNode, N = Cp.getSize(), q = a.getSize();
  bool bAllMerge = 0;
  int  nPossibleMerges = 0, mergeEnergy;

  mergeMatrix.resize(q,q);
  mergeMatrix.init(0);
  // use a counting matrix to correct for missing edges
  countMatrix.resize(q,q);
  countMatrix.init(0);
  //cout << red << "here B " << normText << endl;                  // debugging

  //a.calcEnergy(Cp);  // just in case
  //cout << red << "here C " << normText << endl;                  // debugging
  //a.initNodesList(); // the node to cluster identifier list
  //cout << red << "here D " << normText << endl;                  // debugging
  int uValue = Cp.getMin();  // the unstored value
  //cout << red << "unstored = " << uValue << normText << endl;  // debugging
  //cout << red << "here E with nClusters = "  << a.getSize()
  //     << " and N = " << a.getNNodes() << normText << endl;    // debugging

  //for(int i=0; i<a.getNNodes(); i++) {
  for(int i=0; i<N; i++) {
    iCluster = a.nodes[i];
    //cout << brown << Cp.kij[i] << blue << i << ": " << flush;    // debugging

    for(int j=0; j<Cp.ki(i); j++) {
      jNode    = Cp.ki(i,j);
      jCluster = a.nodes[jNode];
      //cout << "iCluster = " << iCluster << " " << flush;         // debugging
      //cout << "jCluster = " << jCluster << " " << flush;         // debugging

      // merge energy summation by clusters
      if(iCluster!=jCluster) {
        mergeEnergy = Cp(i,jNode);
        //cout << "Cp(" << i << "," << jNode << ") = " << mergeEnergy 
        //     << " for qi = " << iCluster   << "\n";         // debugging
        mergeMatrix[iCluster][jCluster] -= mergeEnergy;
        mergeMatrix[jCluster][iCluster] -= mergeEnergy;
        countMatrix[iCluster][jCluster] += 1;
        countMatrix[jCluster][iCluster] += 1;
      } // end if

    } // end for j
  } // end for i
  //if(bDisplayMerges)  coutCMatrix(mergeMatrix,1);  // debugging
  // we have double counted the energies, so we must correct this now
  mergeMatrix /= 2;  countMatrix /= 2;

  // correct energies for the missing edges and also check whether merges exist
  int eCorrection;
  for(int i=0; i<q; i++) {
    for(int j=i+1; j<q; j++) {
      // add energy for missing edges
      eCorrection = (countMatrix[i][j] - a[i].getn()*a[j].getn())*uValue;
      mergeMatrix[i][j] += eCorrection;
      mergeMatrix[j][i] += eCorrection;

      // check for and count possible merges
      mergeEnergy = mergeMatrix[i][j];
      if(mergeEnergy>0)       bAllMerge = 0;
      else if(mergeEnergy<0)  nPossibleMerges++; // count possible merges

    } // end for j
  } // end for i

  // crudely zero the bottom half of matrix for consistence with the merge code
  // do not know the problem except with the getMin on the symmetric matrix
  for(int i=0; i<q; i++)  for(int j=i+1; j<q; j++)  mergeMatrix[i][j] = 0;

  if(bDisplayMerges)  coutCMatrix(mergeMatrix,1);

  return nPossibleMerges;
} // end getListMergesNZ

inline int getListMergesNZ(TClusterList &a, TCMatrix &Cp, 
       TMatrixInt &mergeMatrix, bool bDisplayMerges=0) {  // a legacy version
  TMatrixInt  countMatrix;
  int nPossibleMerges;
  nPossibleMerges=getListMergesNZ(a,Cp,mergeMatrix,countMatrix,bDisplayMerges);
  return nPossibleMerges;
} // end getListMergesNZ


int  getListMerges(TClusterList &a, TCMatrix &Cp, TMatrixInt &mergeMatrix,
                   TMatrixInt &edgeMatrix, bool bDisplayMerges=0) {
  // generate a matrix of merge energies for all clusters in list a
  int  i = 0, j;
  int  mergeEnergy, minMergeEnergy = BigInteger, nEdgesAB;
  int  minMergeIndex = -1, nPossibleMerges = 0;
  bool bAllMerge = 0;
  int  q = a.getq();

  mergeMatrix.resize(q,q);  mergeMatrix.init(0);
  edgeMatrix.resize(q,q);   edgeMatrix.init(0);

  for(int i=0; i<q; i++) {
    minMergeEnergy = BigInteger;
    minMergeIndex  = -1;

    for(int j=i+1; j<q; j++) {
      mergeEnergy = a[i].getMergeEnergy(a[j],Cp,nEdgesAB);
      edgeMatrix[i][j] = nEdgesAB;
      edgeMatrix[j][i] = nEdgesAB;
      // crude merge checking
      if(mergeEnergy<minMergeEnergy) {
        minMergeIndex  = j;  minMergeEnergy = mergeEnergy;
      } // end if mergeEnergy
      mergeMatrix[j][i] = mergeEnergy;
      if(mergeEnergy>0)       bAllMerge = 0;
      else if(mergeEnergy<0)  nPossibleMerges++;  // count possible merges

    } // end for j
  } // end for i

  if(bDisplayMerges)  coutCMatrix(mergeMatrix,1);         // debugging

  return nPossibleMerges;
} // end getListMerges


int  detectAndMerge(TClusterList &a, TCMatrix &Cp, TMatrixInt &mergeMatrix,
                    bool bDisplayMerges=0) {
  int nPossibleMerges, nMerges = 0;
  int n = a.getSize();
  TMatrixInt edgeMatrix;

  //cout << magenta << "Entering detectAndMerge... " << flush;  // debugging
  #ifdef MSL_TCMWEIGHTED
  nPossibleMerges = getListMergesNZ(a,Cp,mergeMatrix,bDisplayMerges);
  #else
  nPossibleMerges = getListMerges(a,Cp,mergeMatrix,edgeMatrix,bDisplayMerges);
  #endif
  //cout << magenta << "back from get merge... " << flush;  // debugging

  if(nPossibleMerges==n*(n-1)/2) {
    //warningMsg("complete merge detected");  // debugging
    return nPossibleMerges;  // complete merge, so just exit
  }
  else if(nPossibleMerges==0) {
    //warningMsg("No merges detected");  // debugging
    return 0;        // no merges at all, so just exit
  }

  // now start merging clusters starting with the smallest (largest abs value)
  // energy change and (inefficiently) recalculate merge after each one.
  // If we make it to here, there is at least one merge.
  unsigned row, col, minIndex, maxIndex;
  int      minVal;

  // now that we have the merge matrix, we can manipulate the matrix itself so 
  // that we have an O(q^2) operation not an O(N^2) operation as above.
  minVal = mergeMatrix.getMin(col,row,1);  // using symmetric search
  //cout << magenta << "merging clusters " << row << " and " 
  //     << col << normal << endl;  // debugging
  while(minVal<0) {
    // pick the min and max values for the collapse function
    // it is an error if they are equal
    if(col>row) { minIndex = row;  maxIndex = col; }
    else        { minIndex = col;  maxIndex = row; }
    //cout << magenta  << "Merging clusters " << maxIndex << " and " << minIndex 
    //     << red << "\nThe minimum energy was " << minVal << normal << endl;  // debugging

    #ifdef DEBUG_MODE
    cout << green << "Attempting to merge clusters " << magenta << minIndex
         << green << " and " << magenta << maxIndex << blue 
         << " at an energy of " << red << minVal << normText << endl;
    #endif

    //a.display(0,0,"Pre-move "+itos(minVal)+" ");  // debugging
    // merge the clusters and erase maxIndex cluster
    a.move(maxIndex,minIndex,minVal+a[maxIndex].getEnergy(),0); 
    nMerges++;
    //a.display(0,0,"Move "+itos(minVal)+" ");  // debugging

    // omit merge collapse due to unknown bug   
    //mergeMatrix.collapse(minIndex,maxIndex); // resize and collapse energies
    //clear diagonal entry to stay consistent with the cluster implementation
    //mergeMatrix[minIndex][minIndex] = 0;
    //cout << brown << "done.\nNew cluster list is:" << endl;  // debuggging
    //a.display();  // debugging
    //cout << brown << "and the new merge matrix is:\n";  // debuggging
    //coutCMatrix(mergeMatrix,1);

    // there is a different result if we omit the redefinition of mergeMatrix,
    // so there is likely a bug somewhere in the q^2 matrix collapsing routine
    // need to debug this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // but for the moment, we just redo the operation for the matrix
    #ifdef MSL_TCMWEIGHTED
    //nPossibleMerges = getListMergesNZ(a,Cp,mergeMatrix,bDisplayMerges);
    nPossibleMerges = getListMergesNZ(a,Cp,mergeMatrix,0);  // debugging
    #else
    // an O(N^2) version for dense matrices with negative non-constant weights
    nPossibleMerges = getListMerges(a,Cp,mergeMatrix,edgeMatrix,bDisplayMerges);
    #endif

    if(bDisplayMerges && nMerges==0)  coutCMatrix(mergeMatrix,1);  // debugging
    minVal = mergeMatrix.getMin(col,row,1);  // using symmetric search
    //cout << red << "\nThe new minimum energy is " << minVal 
    //     << normal << endl;  // debugging
  } // end while

  return nMerges;
} // end detectAndMerge


/*
int  getMergeEnergy(TCluster &a, TCluster &b, TCMatrix &Cp) {
  // this function has been superceded by an efficient in-class implementation
  // takes two clusters and attempts to detect whether there is a possible 
  // energy barrier preventing a merger with node by node moves
  int         energyCut, energyCutA = 0;
  TVectorInt  edgeEnergiesA(a.getSize());
  //edgeEnergiesA.init(0);

  // these energy cuts are positive due to '-='
  for(int i=0; i<a.getSize(); i++) {
    energyCut   = b.getEnergyChangeAdded(a[i],Cp);
    energyCutA += energyCut;
    //edgeEnergiesA[i] = energyCut;
  } // end for i
  //cout << blueText << "The merge edges are:  " << edgeEnergiesA
  //     << normText << endl;                                      // debugging

  return energyCutA;
} // end getMergeEnergy
*/


bool  detectBarrier(TCluster &a, TCluster &b, TCMatrix &Cp, bool bFastCount=0) {
  // takes two clusters and attempts to detect whether there is a possible 
  // energy barrier preventing a merger with node by node moves
  // function has largely been superceeded by the more efficient NZ versions
  // and the efficient class-based merge detection routine
  int  energyCutA = 0, energyCutB = 0;
  int  maxEdge = 0,   edge;

  // these energy cuts are positive due to '-='
  for(int i=0; i<a.getSize(); i++)
    energyCutA -= b.getEnergyChangeAdded(a[i],Cp);
  for(int i=0; i<b.getSize(); i++)
    energyCutB -= a.getEnergyChangeAdded(b[i],Cp);

  for(int i=0; i<b.getSize(); i++)  for(int j=0; j<b.getSize(); j++) {
      edge = Cp(b.getNode(i),a.getNode(j));
      if(maxEdge>edge)  maxEdge = edge;
    } // end for j

  return (max(energyCutA, energyCutB) < maxEdge);
} // end detectBarrier


int  dynamicMerge(TClusterList &c, TCMatrix &Cp, bool bDisplayMerges = 0) {
  // take a cluster list description, generate the merge energy matrix for it,
  // construct a corresponding new system, and dynamically solve the new system
  // with the existing community detection routine.  Note that the collapsed
  // system will almost certainly be dense (even if only by missing edges); 
  // and thus we use the TCMDense connection matrix to represent
  // it without losing any information.
  // System dereferences mergeMatrix in case that is useful.
  int  nMerges = 0, nPossibleMerges, nSteps;
  int  deltaE, deltaCE, deltaUE, deltaL; // deltaL is the edge count sum
  int  q = c.getq();
  TMatrixInt  mergeMatrix, edgeMatrix;
  
  //cout << "Entering dynamicMerge with q = " << q << "... " << endl; // debugging

  //int  getListMerges(TClusterList &a, TCMatrix &Cp, TMatrixInt &mergeMatrix,
  //                 bool bDisplayMerges=0) {
  nPossibleMerges = getListMerges(c,Cp,mergeMatrix,edgeMatrix,bDisplayMerges);
  //coutCMatrix(mergeMatrix,1);  // debugging
  //coutCMatrix(edgeMatrix,1);   // debugging
  if(nPossibleMerges==0) {
    //cout << red << "No possible merges" << normText << endl;  // debugging
    return 0;
  } // end if nPossibleMerges
  //else  errorMsg("Done.");   // debugging

  // use mergeMatrix to construct a TCMDense object
  //cout << "Declare collapsed connection matrix... " << endl;  // debugging
  TCMDense  Cd(mergeMatrix,q,"Dense Cd constructed in dynamicMerge");
  //cout << "isValid = " << Cd.isValid() << endl;
  
  //cout << "Declare working cluster list" << flush;  // debugging
  // now create the working information for the communityDetect call
  TClusterList working;                 // the working cluster list
  //cout << "... " << endl;  // debugging
  working.createSymmetric(q,1); // create symmetric, randomize
  //working.display(0,0,"Working init ");
  // now dynamically collapse the system
  //cout << "Solving collapsed system... " << flush;  // debugging
  //cout << "isValid = " << Cd.isValid() << endl;
  nSteps = communityDetect(working,Cd); // returned is number of steps required
  //cout << "done with solve. " << endl;  // debugging
  //cout << "isValid = " << Cd.isValid() << endl;
  //cout << "Calculating energy... " << flush;  // debugging
  working.calcEnergy(Cd);
  //cout << "done with energy. " << endl;  // debugging
  
  
  //working.display(0,0,"Working after solve "); // debugging
  //Cd.display();               // debugging
  //coutCMatrix(edgeMatrix,1);  // debugging
  //cout << endl;               // debugging
  //errorMsg("Done.");        // debugging

  // now collapse the passed clusterlist c based on the dynamic results
  // if necessary, this collapsed system can be operated on by a collapse 
  // routine much more efficiently than an original very large system
  // We do not remove the empty lists until the end so that the original 
  // c cluster identities information are unchanged.
  // The working cluster node numbers are the passed cluster list c clusters.
  int iMerge, jMerge;
  int energySum = 0, cEnergySum = 0, uEnergySum = 0;
  //cout << "i = " << flush;
  for(int i=0; i<working.getSize(); i++) {
    iMerge = working[i][0].getNode();
    jMerge = -1; // an invalid starting value used below as an invalid flag
    deltaE = 0;  deltaCE = 0;  deltaUE = 0;  // clear sum variables
    deltaL = 0;
    
    //cout << blue << i << ": " << flush;
    for(int j=1; j<working[i].getSize(); j++) {
      //cout << grey << j << " " << flush;
      jMerge = working[i][j].getNode();
      c.updateNodesList(jMerge,iMerge); // update nodes list with known change
      // we don't know the deltaE here, so use the regular move function
      c[jMerge].move(c[iMerge]);
      nMerges++;
      deltaL  = c[jMerge].nEdges;
      deltaE  = c[jMerge].energy;
      deltaCE = c[jMerge].cEdgeSum;
      deltaUE = c[jMerge].uEdgeSum;
    } // end for j

    if(jMerge>=0) {  // only execute if we have an actual merge occuring
      // now update the energy
      // due to external edges, the only updated information in move() is nNodes
      deltaE  += working[i].energy;
      //cout << blue << "iMerge = " << iMerge 
      //             << " jMerge = " << jMerge << "  " << flush;
      //cout << blue << "energy change is " << deltaE << "  " << flush;
      //deltaCE += working[i].cEdgeSum;  // not exactly correct
      //deltaUE += working[i].uEdgeSum;  // not exactly correct
      deltaL  += edgeMatrix[iMerge][jMerge];
      //cout << blue << "edge change with is " << deltaL << endl;
      c[iMerge].energy   += deltaE;
      //c[iMerge].nNodes += working[i].nNodes;  // automatically updated in move
      c[iMerge].nEdges   += deltaL;
      energySum          += deltaE;
      // the edge sums use the edgeMatrix, so it is only valid if we are using
      // one of the sparse matrix representations
      cEnergySum         += deltaCE;// inconsistent - use edgeMatrix to correct
      uEnergySum         += deltaUE;// inconsistent - use edgeMatrix to correct
      c[iMerge].partialQ  = 0.0;    // inconsistent - no per move update tracked
      c[iMerge].partialQW = 0.0;    // inconsistent - no per move update tracked
    } // end if jMerge
  } // end for i
  c.energyTotal += energySum;
  c.cEdgeSum    += cEnergySum;
  c.uEdgeSum    += uEnergySum;

  // now we can finally (and efficiently) erase all moved clusters  
  int nErased;  nErased = c.clearZeroClusters(1);

  // finally we can check the case of a fully collapsed system (a common case
  // for dynamicMerge) and make the inconsistent variables consistent 
  if(c.getSize()==1) { 
   c[0].partialQ  = 0.0;  c.modularity  = 0.0; 
   c[0].partialQW = 0.0;  c.modularityW = 0.0; 
   //c[0].cEdgeSum  = 0;  c.cEdgeSum  = 0; 
   //c[0].uEdgeSum  = 0;  c.uEdgeSum  = 0; 
  } // end if c.getSize == 1
  // else modularity of c is inconsistent at this point
  // so... take the easy way out and just recalculate the overall energy
  //else c.calcEnergy(Cp);  // this is a costly O(N^2) calc if q is too small

  return nMerges;
} // end detectAndMergeAllNZ


bool detectDisjoint(TClusterList &clusters, TCMatrix &Cp, int verbosity = 0) {
  // Returns a boolean value specifying whether the passed cluster list is
  // disjoint.  That is, it returns true if it has *any* edges connecting any 
  // two communities.  It truncates the function as soon as one edge is found.
  // Note that the implementation uses optimizations that depend on the
  // specific structure of the connection data structures.  
  // Care should be exercised when changing the classes.

  int     i, k;
  int     N = Cp.getSize(), nChecked = 0;  // the number of nodes, " checked
  int     iNode, iCluster, jNode, jCluster, minCluster, jDegrees;
  int    *pNeighborStart,      *pNeighbor; // optimzing neighbor pointers
  bool    bDisjoint = 1;        // is the list disjoint?

  if(clusters.getSize()==1)  return 0;  // trivial undefined case

  // loop over particles and move each one in turn
  i = 0;
  while(i<clusters.getSize() && bDisjoint) {
    clusters[i].begin();  // begin user iteration
    //iCluster = i;       // just to keep clear code
    //cout << blue << i << ": " << normText << flush;  // debugging

    // ----------------------------------------------------------------
    k=0;
    while(k<clusters[i].getSize() && bDisjoint) {
      iNode = clusters[i].getCurrentNode();

      // scan edge connection list for iNode and sum energy contributions
      jDegrees       = Cp.ki(iNode);
      pNeighborStart = &(Cp.kij[iNode][0]); // gives ~7.5% speed boost
      pNeighbor      = pNeighborStart;
      for(int j=0; j<jDegrees; j++) {
        jNode = (*pNeighbor);
        jCluster = clusters.nodes[jNode];
        // check whether this edge is inside its own cluster
        bDisjoint &= (jCluster == i);
        //bDisjoint &= (clusters.nodes[(*pNeighbor)] == i);

        pNeighbor++;  // increment neighbor pointer
      } // end for j - neighbor list search

      nChecked++;          // how many nodes have been checked so far
      clusters[i].next();  // increment to next node for user iteration
      k++;
    } // end while k - particle loop
    // -----------------------------------------------------------------

    i++;  // increment cluster loop
  } // end while i - cluster loop

  if(verbosity>1) {
    cout << blue << "Done with disjoint check after " << nChecked << " nodes.  " 
         << " The partition ";
    if(bDisjoint)  cout << red << "was" << blue << " disjoint\n";
    else           cout << "was " << red << "not" << blue << " disjoint\n";
  } // end if verbosity

  return bDisjoint;
// ---------------------------------------------------------------------------
} // end detectDisjoint
// ---------------------------------------------------------------------------


