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
using namespace MSL;

// ---------------------------------------------------------------------------
// cluster list statistics and analysis
// ---------------------------------------------------------------------------
void nodeStats(TClusterList &c, TCMatrix &Cp, TCMData a, TCMData b, 
       string fRoot, int verbosity = 0,   int debugMode = 0) {
  // Takes a cluster list and a connection matrix and outputs several statistics
  // on both cluster and node levels.   Files are output in csv text format.
  // There is a general limit in OpenOffice of about 1024 columns, so the amount
  // of data is restricted by some multiple of q<1024 where applicable.  
  // The stats currently work with general weights but constant negative weights
  // Output is:
  // Node level statistics:
  //   1. Node ID No.
  //   1. Member of which cluster. (N vector)
  //   3. Degree of node.
  //   2. No. of edges of each *node* connected to nodes in each cluster
  //      (N vector for own cluster, Nxq for all clusters)
  //   3. % of edges of each *node* connected to nodes in each cluster
  //      (N vector for own cluster, Nxq for all clusters)
  //   4. Energy of each *node* for each cluster (This measure is our equivalent 
  //      of the "centrality" or "influence" of a node to each cluster.)
  //      (N vector for own cluster, Nxq for all clusters)
  //   File formats are:
  //      Node member cluster only:  "fRootNodesData.csv" 
  //        col1: (1), col2: (2), col3: (3), col4: (4) 
  //      Node member cluster only:  "fRootNodesDataAll.csv" 
  //        col1: (1), cols 2 - 2q+1: (2), 
  //        cols 2q+2 - 3q+1: (3), cols 3q+2 - 4q+1: (4) 
  //
  //   "Extra stuff at some point"
  //   5. Energy of member cluster on removal.  This measure partially indicates
  //      the "integrity" of the member cluster if the node is removed.  For
  //      unconstrained systems, a resulting positive energy will necessarily 
  //      split the cluster into at least two distinct clusters.  
  //      The integrity is more ambiguous if the cluster still has a negative 
  //      energy, but a quick check would be to check how many nodes are still 
  //      bound to the cluster after removal (partially prevents the need for 
  //      a sub-solution).
  //      a. Do we include an actual sub-solution of the clsuter if it splits?
  //   6. Perhaps implement a specific node information request at some point.
  // fRoot is the root file name to use for the various statistics files.
  int N = Cp.getSize(), L = Cp.getL(), q = c.getq();  // convenience variables
  b = abs(b);  // work with absolute values

  // Begin node level statistics
  string    nsOutput = fRoot + "NodesData.csv";
  if(verbosity>0) 
    cout << brown << "Exporting node statistics to file \"" << green << nsOutput 
         << normText << "\"\n";
  ofstream  nout(nsOutput.c_str());
  nout << "\"Node statistics for cluster list " << fRoot << ":\"\n"
       << "\"Each node's statistics for its own community:\"\n"
       << "\"No.\",\"q_i\",\"k_i\",\"l_qi\",\"p_qi\",\"w_qi\",\"E_qi\"\n";
  // calculate needed statistics - complete set of data
  TMatrixInt    li(N,q,0), wqi(N,q,0), eqi(N,q,0);
  TMatrixFloat  pqi(N,q,0.0);
  int jDegrees, iCluster, jCluster, jNode, kNode, kCluster, ni, nk;
  for(int i=0; i<q; i++) {
    c[i].begin();
    ni = c[i].getn();

    for(int j=0; j<ni; j++) {
      jNode = c[i].getCurrentNode();
      jCluster = c.nodes[jNode];
      jDegrees = Cp.ki(jNode);
      for(int k=0; k<jDegrees; k++) {
        kNode = Cp.ki(jNode,k);
        kCluster = c.nodes[kNode];
        li[jNode][kCluster]  += 1;              // summed number of edges
        wqi[jNode][kCluster] += Cp(jNode,kNode);// summed edge weight (positive)
      } // end for k

      // calculate other stats
      for(int k=0; k<q; k++) {
        // percentages using the completed degree data
        pqi[jNode][k] = (double)li[jNode][k]/(double)jDegrees;
        // energy for each node connected to each respective cluster
        nk = c[k].getn();
        if(k==jCluster)
          eqi[jNode][k] = -wqi[jNode][k] + (nk - li[jNode][k] - 1)*b;
        else
          eqi[jNode][k] = -wqi[jNode][k] + (nk - li[jNode][k]    )*b;
      } // end for k

      c[i].next();
    } // end for j
  } // end for i
  cout << "done with statistics... " << flush;

  // output basic node statistics data
  for(int i=0; i<N; i++) {
    iCluster = c.nodes[i];
    nout << i << "," << c.nodes[i] << "," << Cp.ki(i) << "," 
         << li[i][iCluster]  << "," << pqi[i][iCluster] << "," 
         << wqi[i][iCluster] << "," << eqi[i][iCluster] << "\n";
  } // end for i
  nout << "\n";
  cout << "done with basic data output... " << flush;


  // output all node statistics data
  string sNodeHeader;
  sNodeHeader = "\"Node\",\"q_i\",\"k_i\",\"q_j -> \"";
  for(int i=0; i<q; i++)  sNodeHeader += "," + itos(i);
  sNodeHeader += "\n";

  nout << "\"Each node's statistics for its all communities:\"\n"
       << "\"Edge count\"\n" << sNodeHeader;
  for(int i=0; i<N; i++) {
    nout << i << "," << c.nodes[i] << "," << Cp.ki(i) << ",";
    for(int j=0; j<q; j++)  nout << "," << li[i][j];
    nout << "\n";
  } // end for i
  nout << "\n\"Edge density\"\n" << sNodeHeader;
  for(int i=0; i<N; i++) {
    nout << i << "," << c.nodes[i] << "," << Cp.ki(i) << ",";
    for(int j=0; j<q; j++)  nout << "," << pqi[i][j]; 
    nout << "\n";
  } // end for i
  nout << "\n\"Weight (negative of energy contribution from edges)\"\n"
       << sNodeHeader;
  for(int i=0; i<N; i++) {
    nout << i << "," << c.nodes[i] << "," << Cp.ki(i) << ",";
    for(int j=0; j<q; j++)  nout << "," << wqi[i][j];
    nout << "\n";
  } // end for i
  nout << "\n\"Energy (including missing edge energy contributions)\"\n"
       << sNodeHeader;
  for(int i=0; i<N; i++) {
    nout << i << "," << c.nodes[i] << "," << Cp.ki(i) << ",";
    for(int j=0; j<q; j++)  nout << "," << eqi[i][j];
    nout << "\n";
  } // end for i
  nout << "\n";
  cout << "done with complete data output... " << flush;

  nout.close(); // done with node level statistics
  cout << "closed file. " << endl;
  return;
} // end nodeStats


void clusterStats(TClusterList &c, TCMatrix &Cp, TCMData a, TCMData b, 
       string fRoot, int verbosity = 0,   int debugMode = 0) {
  // Takes a cluster list and a connection matrix and outputs several statistics
  // on both cluster and node levels.   Files are output in csv text format.
  // There is a general limit in OpenOffice of about 1024 columns, so the amount
  // of data is restricted by some multiple of q<1024 where applicable.  
  // The stats currently work with general weights but constant negative weights
  // Output is:
  // Cluster level statistics:
  //   Information:  weighted, unweighted, ...
  //   1. Cluster ID No.
  //   1. Number of nodes in each cluster.
  //   2. Edges within each cluster.
  //   3. Total weight of each cluster's existing edges (no missing edges)
  //   4. Energy of each cluster
  //   5. Node members of each cluster (for small graphs)
  //   6. ?
  // Cluster to cluster statistics:
  //   1. No. of edges of each *cluster* connected to each cluster (qxq matrix)
  //   2. % of edges of each *cluster* connected to each cluster (qxq matrix)
  //   3. Energy merge cost for each cluster pair (qxq matrix)
  //   File formats are:
  //      Node member cluster only:  "fRootClustersData.csv" 
  //        col1: (1), col2: (2), col3: (3), col4: (4) 
  //      Node member cluster only:  "fRootClustersDataAll.csv" 
  //        cols 1 - q: (1), cols q+1 - 2q: (2), cols 2q+1 - 3q: (3)
  //   
  // fRoot is the root file name to use for the various statistics files.
  int N = Cp.getSize(), L = Cp.getL(), q = c.getq();  // convenience variables
  b = abs(b);  // work with absolute values
  
  // Begin node level statistics
  string nsOutput = fRoot + "ClustersData.csv";
  ofstream  nout(nsOutput.c_str());
  cout << "Calculating clusters stats... " << flush;
  // calculate needed statistics - complete set of data
  TMatrixInt    li(q,q,0), wqi(q,q,0), eqi(q,q,0);
  TMatrixFloat  pqi(q,q,0.0);
  int jDegrees, iCluster, jCluster, kCluster;
  int jNode, kNode, mNode, ni, nj, nk, cij;
  double  invnimax, invnm;
  for(int i=0; i<q; i++) {
    c[i].begin();
    ni = c[i].getn();
    invnimax = 2.0/( (double)ni*(double)(ni-1) );
    
    // loop over nodes in current cluster i
    for(int j=0; j<ni; j++) {
      jNode = c[i].getCurrentNode();
      for(int k=j+1; k<ni; k++) {
        kNode = c[i][k].node;
        cij = Cp(jNode,kNode);
        if(cij>0) {
          li[i][i]  += 1;               // summed number of edges
          wqi[i][i] += cij;             // summed edge weight (positive)
          pqi[i][i] += invnimax;        // summed percentage (cluster density)
        } // end if cij
        eqi[i][i]   += cij;             // summed energy
      } // end for k

      c[i].next();
    } // end for j
  } // end for i
  cout << "done with same cluster stats... " << flush;

  
  // sum over all other clusters
  for(int i=0; i<q; i++) {
    ni = c[i].getn();

    for(int j=i+1; j<q; j++) {
      nj = c[j].getn();
      invnm = 1.0/( (double)ni*(double)nj );

      c[i].begin();
      for(int k=0; k<ni; k++) {
        kNode = c[i].getCurrentNode();
        c[j].begin();
        for(int m=0; m<nj; m++) {
          mNode = c[j].getCurrentNode();
          cij = Cp(kNode,mNode);
          if(cij>0) {
            li[i][j]  += 1;             // summed number of edges
            wqi[i][j] += cij;           // summed edge weight (positive)
            pqi[i][j] += invnm;         // summed percentage (cluster density)
          } // end if cij
          eqi[i][j]   += cij;           // summed energy

          c[j].next();
        } // end for m

        c[i].next();
      } // end for k
    } // end for j
  } // end for i
  cout << "done with cluster stats... " << flush;


  // output all cluster statistics data
  nout << "\"Node statistics for cluster list " << fRoot << "\"\n"
       << "\"No.\",\"n\"\n";
  nout << "\"Number of edges connecting communities\"\n\"q_i\",\"n\",\"q_j -> \"";
  for(int i=0; i<q; i++)  nout << "," << i;
  nout << "\n";
  for(int i=0; i<q; i++) {
    nout << i << "," << c[i].getn() << ",";
    for(int j=0; j<i; j++)  nout << "," << li[j][i];
    for(int j=i; j<q; j++)  nout << "," << li[i][j];
    nout << "\n";
  } // end for i
  nout << "\n\"Density of edges connecting communities\"\n\"q_i\",\"n\",\"q_j -> \"";
  for(int i=0; i<q; i++)  nout << "," << i;
  nout << "\n";
  for(int i=0; i<q; i++) {
    nout << i << "," << c[i].getn() << ",";
    for(int j=0; j<i; j++)  nout << "," << pqi[j][i];
    for(int j=i; j<q; j++)  nout << "," << pqi[i][j];
    nout << "\n";
  } // end for i
  nout << "\n\"Weight of edges connecting communities\"\n\"q_i\",\"n\",\"q_j -> \"";
  for(int i=0; i<q; i++)  nout  << "," << i;
  nout << "\n";
  for(int i=0; i<q; i++) {
    nout << i << "," << c[i].getn() << ",";
    for(int j=0; j<i; j++)  nout << "," << wqi[j][i];
    for(int j=i; j<q; j++)  nout << "," << wqi[i][j];
    nout << "\n";
  } // end for i
  nout << "\n\"Merge energy required for paired communities\"\n\"q_i\",\"n\",\"q_j -> \"";
  for(int i=0; i<q; i++)  nout << "," << i;
  nout << "\n";
  for(int i=0; i<q; i++) {
    nout << i << "," << c[i].getn() << ",";
    // reverse sign on the energies so that it is an energy *cost* to merge
    for(int j=0; j<i; j++)  nout << "," << -eqi[j][i];
    for(int j=i; j<q; j++)  nout << "," << -eqi[i][j];
    nout << "\n";
  } // end for i
  nout << "\n";
  cout << "done with complete data output... " << flush;

  nout.close(); // done with node level statistics
  return;
} // end clusterStats


