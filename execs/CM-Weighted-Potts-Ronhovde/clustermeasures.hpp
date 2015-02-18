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
// equals tests and measures
// ---------------------------------------------------------------------------
double infoEntropy(TClusterList &c, int N = -1) {
  // calculate the information entropy (Shannon entropy) of the cluster list
  // solution using p_i = n_i/N
  if(c.getq()==1)  return 0.0;

  double  info = 0.0, p;
  // if N is not given, use the cluster list value (the cluster list node count
  // is not updated until a calcEnergy() call)
  if(N==-1)  N = c.getNNodes(); 

  for(int i=0; i<c.getSize(); i++) {
    #ifdef DEBUG_MODE
    if(c[i].getSize()<1)
      errorMsg("cluster "+itos(i)+" has size zero in infoEntropy?");
    #endif
    p = (double)( c[i].getSize() )/(double)N;
    //if(c[i].getSize()>1)  info -= p*log10(p)/log10(2.0);
    info -= p*log10(p)/log10(2.0);
  } // end for i

  //msg("done\n");
  return info;
}; // end infoEntropy


double  calcNMVI(TClusterList& a, TClusterList& b, int NNodes,
                 double &nmi, double &mi, double &vi) {
  // calculate Normalized Mutual Information (NMI), Mutual Information (mi),
  // and Variation of Information (vi) for two given cluster lists all in one
  // function since they are all very closely related
  // requires an existing valid cluster with compatible maximum sizes
  if(a.getNNodes()!=b.getNNodes())  
    errorMsg("Number of nodes do not match in calcNMVI!");

  double shannonA = 0.0, shannonB = 0.0, probX, probY;
  if(a.getSize()==1 || b.getSize()==1) {
    mi = 0.0;  nmi = 0.0;
    shannonA = infoEntropy(a,NNodes);
    shannonB = infoEntropy(b,NNodes);
    vi = shannonA + shannonB;
    return vi;
  } // end size 1's

  // confusion matrix N which is usually not symmetric
  TMatrixInt  N(b.getSize(),a.getSize(),0);
  // fill confusion matrix N
  //for(int i=0; i<NNodes; i++)  N[clusterVectorB[i]][clusterVectorA[i]] += 1;
  for(int i=0; i<NNodes; i++)  N[b.nodes[i]][a.nodes[i]] += 1;

  // now calculate Mutual Information value I(a,b)
  double  ISum = 0.0;
  int     Nij;
  for(int i=0; i<a.getSize(); i++) {
    for(int j=0; j<b.getSize(); j++) {
      Nij = N[j][i];
      // must catch zero entries since the sum term is only zero in the limit
      //if(Nij>0)  ISum += (double)Nij/(double)NNodes
      //  *log((double)( Nij*NNodes )/(double)( a[i].getSize()*b[j].getSize() ));
      if(Nij>0)  ISum += (double)Nij/(double)NNodes
                        *log((double)Nij*(double)NNodes
                        /((double)a[i].getSize()*(double)b[j].getSize()));
    } // end for j
  } // end for i
  ISum /= log(2.0);  // convert to bits

  // set reference Mutual Information value
  mi = ISum;  
  
  // calculate the shannon entropy of the two data sets to prepare for VI calc
  for(int i=0; i<a.getSize(); i++) {
    probX = (double)a[i].getSize()/(double)NNodes;
    shannonA -=  probX*log(probX);
  } // end for i
  for(int j=0; j<b.getSize(); j++) {
    probY = (double)b[j].getSize()/(double)NNodes;
    shannonB -=  probY*log(probY);
  } // end for j
  shannonA /= log(2.0);  // convert to bits
  shannonB /= log(2.0);  // convert to bits
 
  // finally calculate the variation in information
  vi = shannonA + shannonB - 2.0*mi; // set reference VI value

  // calculate the NMI value by the normalization term.  Cancelled '-' sign 
  // since shannonX is defined with a negative sign.
  // must catch zero case separately
  if(mi>0.0)  nmi = 2.0*mi/( shannonA + shannonB );
  else        nmi = 0.0;

  return vi;
} // end calcNMVI


// ---------------------------------------------------------------------------
// fast measures
// ---------------------------------------------------------------------------
double  calcNMVILL(TClusterList& a, TClusterList& b, int NNodes,
                   double &nmi, double &mi, double &vi, int &nMoved, 
                   double &shannonA, double &shannonB, bool bCheck = 0) {
  // calculate Normalized Mutual Information (NMI), Mutual Information (mi),
  // and Variation of Information (vi) for two given cluster lists all in one
  // function since they are all very closely related
  // requires an existing valid cluster with compatible maximum sizes
  //#ifdef DEBUG_MODE
  if(a.getq()==0 || b.getq()==0)  
    errorMsg("Attempting to calcNMVILL on a size 0 cluster!");
  if(a.getNNodes()!=b.getNNodes())  
    errorMsg("Number of nodes do not match in calcNMVILL (N_a = "
             +itos(a.getNNodes())+" and N_b = "+itos(b.getNNodes())+")!");
  //#endif

  nMoved = 0;
  // catch trivial cases
  if(a.getq()==1 || b.getq()==1) {
    mi = 0.0;  nmi = 0.0;  
    shannonA = infoEntropy(a,NNodes);
    shannonB = infoEntropy(b,NNodes);
    vi = shannonA + shannonB;
    return vi;
  } // end size 1's

  TCMDataNode w;
  // confusion matrix N which is usually not symmetric
  TListArray<TCMDataNode>  N(b.getSize());
  // fill confusion matrix N
  w.value = 1;
  TCMDataNode *pn = NULL;
  for(int i=0; i<NNodes; i++) {
    //N[b.nodes[i]][a.nodes[i]] += 1;
    w.index = a.nodes[i];
    pn = N[b.nodes[i]].getMember(w);
    if(pn==NULL)  N[b.nodes[i]].add(w);  // not found, so add it
    else          pn->value += 1;        // found it, so increment counter
  } // end for i

  //cout << N << endl;

  // now calculate Mutual Information value I(a,b)
  // this version only checks the non-zero entries
  TVectorInt  nSum(b.getq(),0), naMax(a.getq(),-1), nbMax(b.getq(),-1);
  double      ISum = 0.0, doubleN = (double)NNodes;
  int         Nij, jIndex;
  for(int i=0; i<b.getq(); i++) {
    if(N[i].getSize()>0) {  // stops a legacy error check on begin()
      N[i].begin();  // begin manual iterator to linked-list

      //for(int j=0; j<a.getSize(); j++) {
      for(int j=0; j<N[i].getSize(); j++) {
        //Nij = N[j][i];
        Nij    = N[i].getCurrent().value;
        jIndex = N[i].getCurrent().index;
        // track variables that find number of misplaced nodes
        nSum[i] += Nij;
        if(Nij>nbMax[i])       nbMax[i] = Nij;       // store row max, b
        if(Nij>naMax[jIndex])  naMax[jIndex] = Nij;  // store column max, a
        // must catch zero entries since the sum term is only zero in the limit
        //if(Nij>0)  ISum += (double)Nij/(double)NNodes
        //  *log((double)( Nij*NNodes )/(double)( a[i].getSize()*b[j].getSize() ));
        if(Nij>0)  ISum +=  (double)Nij/doubleN
                  *log( (double)Nij*doubleN/( a[jIndex].getnD()*b[i].getnD() ));
        else       errorMsg("Found zero entry in calcNMVILL()?");
        N[i].next();  // increment manual iterator to next linked-list element
      } // end for j
    } // end if size check
    else  errorMsg("Found zero cluster size in calcNMVILL()?");

    if(bCheck && b[i].getn()==1)  
      cout << i << " (" << b[i].getStartNode() << ")  " << flush;  // debugging
  } // end for i
  ISum /= log(2.0);  // convert to bits

  // set reference Mutual Information value
  mi = ISum;

  // find number of misplaced nodes - assumes no overlapping nodes
  // picking smallest q between cluster lists a and b makes the logic easier
  if(a.getq()<b.getq())
        for(int i=0; i<a.getq(); i++)  nMoved += a[i].getn() - naMax[i];
  else  for(int i=0; i<b.getq(); i++)  nMoved += b[i].getn() - nbMax[i];
  //cout << red << "No. moved nodes = " << nMoved << normText << endl; // debugging
  
  // calculate the shannon entropy of the two data sets to prepare for VI calc
  //cout << "N = " << NNodes << endl; // debugging
  shannonA = 0.0;  shannonB = 0.0;
  double  probX, probY;
  for(int i=0; i<a.getSize(); i++) {
    #ifdef DEBUG_MODE
    if(a[i].getSize()<1)
      errorMsg("cluster "+itos(i)+" has size zero in calcNMVILL?");
    #endif
    probX = (double)a[i].getSize()/(double)NNodes;
    shannonA -=  probX*log(probX);
    //cout << "i = " << i << ":  " << probX << " " << shannonA << endl; // debugging
  } // end for i
  for(int j=0; j<b.getSize(); j++) {
    #ifdef DEBUG_MODE
    if(b[j].getSize()<1)
      errorMsg("cluster "+itos(j)+" has size zero in calcNMVILL?");
    #endif
    probY = (double)b[j].getSize()/(double)NNodes;
    shannonB -=  probY*log(probY);
  } // end for j
  shannonA /= log(2.0);  // convert to bits
  shannonB /= log(2.0);  // convert to bits

  // finally calculate the variation in information
  vi = shannonA + shannonB - 2.0*mi; // set reference VI value
  // minimum possible VI is 2/N, so test and set to zero if below this threshold
  // we pick a very small epsilon since we run on almost O(10^8) nodes
  if(vi<(2.0/(double)NNodes-5.0e-12))  vi = 0.0;

  //cout << "H_A = " << shannonA << "  H_B = " << shannonB << endl; // debugging

  // calculate the NMI value by the normalization term.  Cancelled '-' sign 
  // since shannonX is defined with a negative sign.
  // must catch zero case separately
  if(mi>0.0)  nmi = 2.0*mi/( shannonA + shannonB );
  else        nmi = 0.0;

  return vi;
} // end calcNMVILL
inline double calcNMVILL(TClusterList& a, TClusterList& b, int NNodes,
            double &nmi, double &mi, double &vi, int &nMoved, bool bCheck = 0) {
  double dumA, dumB;
  calcNMVILL(a,b,NNodes,nmi,mi,vi,nMoved,dumA,dumB,bCheck);
  return vi;
} // end calcNMVILL without information references


// debugging
double calcMod(TClusterList &a, TCMatrix &Cp, double gamma=1.0) {
  // Assumes Cp is one of the sparse matrix representations.
  // The function returns the maximum RB modularity detected over the 
  // scan (positive or negative).
  if(a.getNNodes()==0)  errorMsg("Number of nodes = 0 in calcMod()?");
  if(a.getNNodes()!=Cp.getSize())  
    errorMsg("Number of nodes in list and matrix does not match in calcMod()?");

  int  N = Cp.getSize(), q = a.getq();  // total number of nodes and clusters
  int  iCluster, jCluster, iNode, jNode, jDegrees, nEdges = 0;
  int *pNeighbor; // optimzing neighbor pointer
  double modSum = 0.0, g4LD = 0.25*gamma/(double)Cp.getL();

  //TMatrixInt  edgeMatrix(q,q,0); // number of nodes connected to a cluster
  //TVectorInt  degrees(q,0);      // cluster degree sums
  int  degreeSum = 0, degreeTotal = 0, iClusterEdges;
  // loop over particles and move each one in turn
  for(iCluster=0; iCluster<q; iCluster++) {
    //a[iCluster].begin();  // begin user iteration
    iClusterEdges = 0;
    degreeSum = 0;

    // ----------------------------------------------------------------
    for(int k=0; k<a[iCluster].getn(); k++) {
      //iNode = a[iCluster].getCurrentNode();
      iNode = a[iCluster][k].node;
      // scan edge connection list for iNode and sum energy contributions
      jDegrees  = Cp.ki(iNode);
      pNeighbor = &(Cp.kij[iNode][0]);
      for(int j=0; j<jDegrees; j++) {
        jNode = (*pNeighbor);
        jCluster = a.nodes[jNode];
        // count number of nodes jNode is connected to in jCluster
        // double counts for iCluster
        nEdges++;
        if(iCluster==jCluster)  iClusterEdges++;

        pNeighbor++;  // increment neighbor pointer
      } // end for j - neighbor list search
      // now increment the sum for the total degree of cluster i
      degreeSum += Cp.ki(iNode);

      //a[iCluster].next();  // increment to next node for user iteration
    } // end for k - node loop
    // track total degree
    degreeTotal += degreeSum;

    // sum edge contribution to modularity, correct for double counting edges
    modSum += 0.5*(double)iClusterEdges - g4LD*pow((double)degreeSum,2);
  } // end for iCluster - cluster loop
  modSum /= (double)Cp.getL();  // correct for normalized scale
  // -------------------------------------------------------------------

  // error checks
  nEdges /= 2;  // correct for double counted edges
  int dCount = degreeTotal/2;
  if(nEdges!=Cp.getL() || dCount!=Cp.getL() || nEdges!=dCount)
    errorMsg("Number of edges in matrix "+itos(Cp.getL())
            +" does not match the counted value of "+itos(nEdges)
            +" or degree counted as "+itos(dCount)+" in calcMod()?");

  return modSum;
} // end calcMod
// debugging
