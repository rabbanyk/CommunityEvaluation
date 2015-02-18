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
#include "./msl/MSL_VectorBag_Template.h"
#endif
//#ifndef _ML_UTILS
//#include "ML_Utils.h"
//#endif
#include "./msl/MSL_Stats1D.h"
*/
using namespace MSL;

// ---------------------------------------------------------------------------
// main hierarchical community detection functions
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// get "best" solution sub-functions
// ---------------------------------------------------------------------------
int communityBestE(TCMatrix &Cp, TClusterList &clusters, TClusterList &best, 
         double &nIterations, double &nMerges, double &solveTime,
         int a = 1, int b = -1, TInitType initType = Symmetric, int q = -1,
         int NTrialsMax = 4, int NIterMax = 1000, 
         bool useMerges = 1, bool useZeroE = 0, bool useConstrainedq = 0, 
         int verbosity = 0) {
  // This is essentially a concise repeat of the communitySolve() function
  // nIterations and solveTime are both the average per trial (not totals)
  int     nIter = -1, nIterR, iterCount, iBest;
  //double  solveTime;
  int     minE = BigInteger;   // energy
  int     qNMVICutoff = 3000;  // rough limit on size of merge tests
  nMerges = 0.0;

  if(initType!=Symmetric && q<2)
    errorMsg("Must specify a valid q for a non-symmetric init. (q = "
             +itos(q)+")");
  else if(q>0)  warningMsg("Ignoring q for symmetric init. ("+itos(q)+")");
/*
  // debugging
  if(useMerges)        cout << brown << "Merges are on\n";
  else                 cout << red   << "Merges are off\n";
  if(useConstrainedq)  cout << red   << "q is fixed\n";
  else                 cout << brown << "q varies\n";
  if(useZeroE)         cout << brown << "Zero energy moves are on\n";
  else                 cout << red   << "Zero energy moves are off\n";
  cout << "NIterMax = " << NIterMax << "\n";
  cout << "NTrialsMax = " << NTrialsMax << "\n";
  cout << normText << "\n";
  // debugging
*/
  localInitTime();  // prep for solve time sum
  solveTime = 0.0;  iterCount = 0; // summing variables
  for(int iTrials=0; iTrials<NTrialsMax; iTrials++) {
    if(verbosity>2)  
      cout << "Starting trial " << iTrials << "... " << flush;
    // initialize the cluster we will solve
    //clusters.initSymmetric(1);  // use temporary init
    //if(initType==Symmetric)  
    //cout << "init... " << flush; // debugging
    #ifdef USE_POWER_INIT
    clusters.initPowerS(Cp.getSize(),8,24,-1,1);  // temporary
    #else
    //if(initType==Symmetric)  
      clusters.initSymmetric(1);
    #endif
    //clusters.display(0,0,"After init in old communityBestE()");  // debugging
    //else if(initType==Even)  clusters.initEvenS(q,1); // init to even state
    //else                     clusters.initRandomS(q); // init to random state
    //nIter = communityDetectNZRedo(clusters,Cp,a,b,NIterMax,useZeroE,0,0);
    //cout << "solve with Cp size of " << Cp.getSize() << "... " 
    //     << flush; // debugging
    if(verbosity>2)  cout << "Solving trial " << iTrials << "... " << flush;
    //clusters.display(0,0,"initialized clusters in solve "); // debugging
    //Cp.display();  // debugging
    nIter = communityDetectNZ(clusters,Cp,a,b,NIterMax,
                              useZeroE,useConstrainedq,verbosity-2,0);
    //cout << "energy... " << flush; // debugging
    clusters.calcEnergy(Cp);  // get starting energy
    //clusters.display(0,0,"After solve in old communityBestE()");  // debugging

    // use merge test
    //cout << "merge... " << flush; // debugging
    if(clusters.getq()<qNMVICutoff && useMerges && !useConstrainedq) {
      if(verbosity>2)  cout << red << "Merge " << flush;
      // this merge routine handles completely general merges (even those 
      // with non-constant negative weights)
      nMerges = 0;          // we store the total number of merges
      int nThisLoop, nMergeLoops = 0;  // track how many loops were made
      do {
        nThisLoop = -1;
        // use dense version - this catches fully merged systems better
        //nMerges += detectAndMergeAllNZ(clusters,Cp,a,b);
        //TMatrixInt  mergeMatrix(clusters.getq(),clusters.getq());
        //nThisLoop = detectAndMerge(clusters,Cp,mergeMatrix,0);
        //nMerges = dynamicMerge(clusters,Cp,0);  // debugging
        // use fast version
        nThisLoop = detectAndMergeAllNZ(clusters,Cp,a,b);
        //nThisLoop = dynamicMerge(clusters,Cp,0);  // debugging
        //cout << red << "post " << flush;  // debugging
        // perform a node-level refinement on the merged solution
        nIterR = 0;
        if(nThisLoop>0)
          //nIterR = communityDetectNZRedo(clusters,Cp,a,b,NIterMax,useZeroE,0,0);
          nIterR = communityDetectNZ(clusters,Cp,a,b,NIterMax,
                                     useZeroE,useConstrainedq,verbosity-2,0);
        if(verbosity>2)
          cout << "(nMerges = " << nThisLoop << ")... "<< normText << flush;
        if(nIterR>1 && verbosity>2)
          cout << red << "Performed " << nIterR << " iterations post "
               << "merge." << normText << " with " << nMergeLoops 
               << " loops\n";
        nMergeLoops++;  nMerges += nThisLoop;
      } while(nThisLoop>0 && nMergeLoops<NIterMax); // end while merge loop
      if(nMergeLoops>=NIterMax)
        warningMsg("Number of merge loops was maxed at NIterMax?");
    } // end if merge condition check 

    clusters.calcEnergy(Cp);  // perform final manual energy calculation
    if(clusters.getE()<minE) {
      minE = clusters.getE();  best = clusters;
    } // end if

    iterCount += nIter;
  } // end for iTrials
  // store the solution time for possible phase change analysis
  solveTime   = getLocalTime();
  solveTime  /= (double)NTrialsMax;  // convert to average time per trial
  nMerges    /= (double)NTrialsMax;  // convert to average merges per trial
  nIterations = (double)iterCount/(double)NTrialsMax;
  if(verbosity>2)  cout << green << "done with trials.\n"; // debugging

  return nIter;  // returns number of iterations on best trial (low E, high Q)
} // end communityBestE

inline double alpha(int a, int b) {
  // calculate the plot parameter alpha from the weights
  return ( -(double)(a+b)/(double)min(a,abs(b)) );
} // end alpha
inline double gamma(int a, int b) {
  // calculate the plot parameter alpha from the weights
  return ( -(double)b/(double)a );
} // end gamma


inline void gammaDensity(int &a, int &b, double density, int abScale = 100) {
  // calculate the approximate a and b for our solvers given a minimum density
  if(density>0.5) {a =  abScale;  b = -(int)((double)a*density/(1.0-density));}
  else            {b = -abScale;  a = -(int)((double)b*(1.0-density)/density);}
  return;
} // end gammaDensity


double gammaLogStepMRA(int &a, int &b, double &gammaAPM, double gammaLogStep, 
                       int verbosity=0) {
  // take the current integer weights a and b and take a geometric step size
  // specified by gammaLogStep.  Return value is the new gammaAPM.
  gammaAPM *= gammaLogStep;
  if(verbosity>3)
    cout << red << "Resetting a = " << a << " and b = " << b << ":   ";
  // now calculate the approximate a and b for our solvers
  if(gammaAPM>1.0) { a =  100;  b = -(int)( gammaAPM*100.0 + 0.5 ); }
  else             { b = -100;  a =  (int)( 100.0/gammaAPM + 0.5 ); }
  if(verbosity>3)  
    cout << red << ".  The new a = " << a << " and b = " << b << endl;

  // catch invalid cases for a or b
  if(b==0)  errorMsg("b = 0 in gammaLogStepMRA()");
  if(a==0)  errorMsg("b = 0 in gammaLogStepMRA()");
  
  return gammaAPM;
} // end gammaLogStepMRA


double gammaDensityStepMRA(int &a, int &b, double gammaAPM,
         double density, double densityStep, double densityScale, 
         double density_t = 0.1, double epsilon = 5.0e-8, int verbosity=0) {
  // take the current integer weights a and b and take a density step size.
  // Assumes an unweighted original graph (not current connection matrix).
  // Return value is the new gammaAPM.
  // density is the current density value
  // densityStep is the additive size of the step above density_t
  // densityScale is the multiplicative scaling below density_t
  if(verbosity>3)  
    cout << red << "Resetting a = " << a << " and b = " << b << ":   ";
  if(density>(density_t+epsilon))  density -= densityStep;
  // for the very sparse systems or end region of a run
  else                             density *= densityScale;
  gammaDensity(a,b,density);  // get the a and b values based on the density
  // a density of 0.5 is the transition between incrementing a or b
  //if(density>0.5) { a =  100;  b = -(int)((double)a*density/(1.0-density)); }
  //else            { b = -100;  a = -(int)((double)b*(1.0-density)/density); }
  if(verbosity>3) 
    cout << red << ".  The new a = " << a << " and b = " << b << endl;

  // catch invalid cases for a or b
  if(b==0)  errorMsg("b = 0 in gammaDensityStepMRA()");
  if(a==0)  errorMsg("b = 0 in gammaDensityStepMRA()");
  
  return gammaAPM;
} // end gammaLogStepMRA

/*
inline double infoMRA(TClusterListArray &ca, double &nmiAvg, double &viAvg, 
                    double &miAvg, double &infoAvg, int verbosity=0) {
  // calculate the replica information correlations
  // loop over all replica pairs (for NMI and VI calculations)
  int NReplicas = ca.getSize();
  TStats1D  nmiVals("NMI_infoMRA()"), viVals("VI_infoMRA()"), 
            miVals("I_infoMRA()"),    infoVals("H_infoMRA()");
  double infoA, infoB;

  if(verbosity>0) {
    localInitTime();
    msg("Starting correlations ("+itos(nReplicaPairs)+" max)... ",brown);
  } // end if verbosity
  for(int k=0; k<NReplicas-1; k++) { // skip last, already accounted for it
    if(verbosity>0)  cout << brown << k << " [ " << grey;

    for(int j=k+1; j<NReplicas; j++) {
      if(verbosity>0)  cout << j << " " << flush;
      calcNMVILL(stored[j],stored[k],N,nmiVal,miVal,viVal,nMoved,infoA,infoB); 
      nmiVals.store(nmiVal);  miVals.store(miVal);
      viVals.store(viVal);    movedVals.store(nMoved);
      if(nMoved>0 && verbosity>2)
        cout << "(nMoved = " << nMoved << ") " << flush; // debugging
    } // end for j

    infoVals.store(infoA);  // only store information for the j cluster list
    if(verbosity>0)  timeMsg("] "," ");
  } // end for k
  nmiAvg  = nmiVals.avg();  miAvg = miVals.avg();  viAvg  = viVals.avg();
  // add the last information variable for list k since we consider pairs of 
  // cluster lists above
  infoVals.store(infoB);    infoAvg = infoVals.avg();  
  if(verbosity>0)  cout << grey << " done." << normText << endl;

  return 0.0;
} // end replica information correlations
*/
inline double nmviMRA(TClusterListArray &ca, int N, int NReplicas, 
     TStats1D &nmis,  TStats1D &vis, TStats1D &mis,  TStats1D &infos, 
     TStats1D &moves, int verbosity=0) {
  // calculate the replica information correlations
  // This version uses reference TStats1D variables
  // loop over all replica pairs (for NMI and VI calculations)
  // crude parameter pass as a workaround.  TClusterListArray is not properly
  // tracking the number of lists stored.
  //int NReplicas = ca.getSize();  
  double nmiVal, viVal, miVal, infoA, infoB;
  int    nMoved;

  if(verbosity>0) {
    localInitTime();
    msg("Starting MRA information correlations... ",brown);
  } // end if verbosity
  for(int k=0; k<NReplicas-1; k++) { // skip last, already accounted for it
    if(verbosity>0)  cout << brown << k << " [ " << grey;

    for(int j=k+1; j<NReplicas; j++) {
      if(verbosity>0)  cout << j << " " << flush;
      calcNMVILL(ca[j],ca[k],N,nmiVal,miVal,viVal,nMoved,infoA,infoB); 
      nmis.store(nmiVal);  mis.store(miVal);
      vis.store(viVal);    moves.store(nMoved);
      if(nMoved>0 && verbosity>2)
        cout << "(nMoved = " << nMoved << ") " << flush; // debugging
    } // end for j

    infos.store(infoA);  // only store information for the j cluster list
    if(verbosity>0)  timeMsg("] "," ");
  } // end for k

  // add the last information variable for list k since we consider pairs of 
  // cluster lists above
  infos.store(infoB);
  if(verbosity>0)  cout << grey << " done." << normText << endl;

  return 0.0;
} // end replica information correlations using TStats1D variables


inline double pzStatsMRA(TCMatrix &Cp, TClusterListArray &ca, int N, int NReplicas, 
     TStats1D &zins,  TStats1D &zouts, TStats1D &zs,  TStats1D &pins, 
     TStats1D &pouts, TStats1D &ps, TStats1D &mus, int verbosity=0) {
  // calculate the replica degree statistics
  // This version uses reference TStats1D variables
  // loop over all replica pairs (for NMI and VI calculations)
  // crude parameter pass as a workaround.  TClusterListArray is not properly
  // tracking the number of lists stored.
  //int NReplicas = ca.getSize();  
  double zAvg, zinAvg, zoutAvg, pAvg, pinAvg, poutAvg;

  if(verbosity>0) {
    localInitTime();
    msg("Starting MRA information correlations... ",brown);
  } // end if verbosity
  for(int k=0; k<NReplicas-1; k++) { // skip last, already accounted for it
    ca[k].calcPZParams(Cp,zinAvg,zoutAvg,zAvg,pinAvg,poutAvg,pAvg);
    zs.store(zAvg);  zins.store(zinAvg);  zouts.store(zoutAvg);  
    ps.store(pAvg);  pins.store(pinAvg);  pouts.store(poutAvg);
    mus.store(zoutAvg/zAvg);  
  } // end for k

  return 0.0;
} // end replica information correlations using TStats1D variables
