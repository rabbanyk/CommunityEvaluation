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

void gammaToab(double &gammaAPM, int &a, int &b, int verbosity=0, 
               int    abScale=100) {
  // take the current gamma weight and set the correct a and b values
  if(verbosity>3)
    cout << red << "Resetting a = " << a << " and b = " << b << ":   ";
  // now calculate the approximate a and b for our solvers
  if(gammaAPM>1.0) { 
    a =  abScale;  b = -(int)( gammaAPM*(double)abScale + 0.5 );
  } else { 
    b = -abScale;  a =  (int)( (double)abScale/gammaAPM + 0.5 ); 
  } // end else
  if(verbosity>3)  
    cout << red << ".  The new a = " << a << " and b = " << b << endl;

  // catch invalid cases for a or b
  if(b==0)  errorMsg("b = 0 in gammaLogStepMRA()");
  if(a==0)  errorMsg("b = 0 in gammaLogStepMRA()");
  
  return;
} // end gammaLogStepMRA


int communityBestE(TCMatrix &Cp, TClusterList &best, 
         double &nIterAvg, double &nMerges, double &timeAvg,
         double gammaAPM, TInitClusterList &initMethod,
         int NTrialsMax = 4, int NIterMax = 1000, 
         bool useMerges = 1, bool useZeroE = 0, bool useConstrainedq = 0, 
         int verbosity = 0) {
  // This is essentially a concise repeat of the communityBestE() function
  // except that it is based on a passed gammaAPM (rather than a and b)
  // nIterAvg and timeAvg are both the average per trial (not totals)
  int     nIter = -1, nIterR, iterCount, iBest;
  int     minE = BigInteger, N = Cp.getSize();
  int     qNMVICutoff = 3000;  // rough limit on size of merge tests

  // error checks
  if(NTrialsMax<1)  errorMsg("Number of trials "+itos(NTrialsMax)
                            +" is invalid in general communityBestE()?");
  if(NTrialsMax<1)  errorMsg("Maximum number of iterations "+itos(NIterMax)
                            +" is invalid in general communityBestE()?");

  // get the integer weights for our integer valued energy and connection matrix
  int  a, b;
  gammaToab(gammaAPM,a,b);
/*
  // debugging
  cout << lightred << "gammaAPM = " << gammaAPM << " with a = " << a 
       << " and b = " << b << normText << endl;    // debugging
  //Cp.display();
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
  // declare a working cluster list
  TClusterList clusters, *pWorking;
  clusters.createSymmetric(N,1);  best.initSymmetric(1);  // debugging
  if(NTrialsMax>1)  pWorking = &clusters;
  // if NTrialsMax == 1, use a faster solve directly on best.  This saves
  // a redundant copy if we use a working cluster list.
  else              pWorking = &best;
  
  localInitTime();  // prep for solve time sum
  iterCount = 0;  nMerges = 0.0;  // summing variables
  for(int iTrials=0; iTrials<NTrialsMax; iTrials++) {
    if(verbosity>2)  cout << "Starting trial " << iTrials << "... " << flush;
    // initialize the cluster we will solve
    initMethod.init(*pWorking);  // call virtual init function
    //pWorking->display(0,0,"After init in communityBestE()");  // debugging

    if(verbosity>2)  cout << "Solving trial " << iTrials << "... " << flush;
    //clusters.display(0,0,"initialized clusters in solve "); // debugging
    //Cp.display();  // debugging
    nIter = communityDetectNZ(*pWorking,Cp,a,b,NIterMax,
                              useZeroE,useConstrainedq,verbosity-2,0);
    //cout << "energy... " << flush; // debugging
    pWorking->calcEnergy(Cp);  // get starting energy
    //pWorking->display(0,0,"After solve in communityBestE()");  // debugging

    // use merge test
    //cout << "merge... " << flush; // debugging
    if(pWorking->getq()<qNMVICutoff && useMerges && !useConstrainedq) {
      if(verbosity>2)  cout << red << "Merge " << flush;
      // this merge routine handles completely general merges (even those 
      // with non-constant negative weights)
      nMerges = 0;          // we store the total number of merges
      int nThisLoop, nMergeLoops = 0;  // track how many loops were made
      do {
        nThisLoop = -1;
        // dense version - this catches fully merged systems better
        //nMerges = dynamicMerge(clusters,Cp,0);  // debugging
        // fast version using neighbor search
        nThisLoop = detectAndMergeAllNZ(*pWorking,Cp,a,b);
        //nThisLoop = dynamicMerge(*pWorking,Cp,0);  // debugging

        // perform a node-level refinement on the merged solution
        nIterR = 0;
        if(nThisLoop>0)
          nIterR = communityDetectNZ(*pWorking,Cp,a,b,NIterMax,
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

    pWorking->calcEnergy(Cp);  // perform final manual energy calculation
    if(pWorking->getE()<minE && NTrialsMax>1) {
      minE = pWorking->getE();  best = *pWorking;
    } // end if

    iterCount += nIter;
  } // end for iTrials
  // store the solution time for possible phase change analysis
  timeAvg   = getLocalTime()/(double)NTrialsMax; // average time per trial
  nMerges  /= (double)NTrialsMax;  // convert to average merges per trial
  nIterAvg  = (double)iterCount/(double)NTrialsMax;
  if(verbosity>2)  cout << green << "done with trials.\n"; // debugging

  return nIter;  // returns number of iterations on best trial (low E, high Q)
} // end communityBestE


double rnMRA(TCMatrix &Cp, string outputFileName, TClusterList &best,
  TInitClusterList &initMethod,
  double gammaStart=0.001, double gammaStop = 25.0, int NDecadeSteps = 20,
  bool   useEnergy = 0,  bool useMerges = 1, bool useZeroEMoves = 0, 
  int    NReplicasMax=8, int NTrialsMax=4, int NStepsMax=20, int NIterMax=1000,
  int    nodeOffset=0,   int verbosity=0, string sRunDescription = "") {
  // a - the starting value for the connected edge weighting
  // b - the starting value for the unconnected edge weighting
  // clusters - the input clusters list
  // useEnergy = 1;         // optimize energy or RB modularity?
  // useZeroEMoves = 1;     // optimize energy or RB modularity?
  double  epsilon = 0.00000005;
  int     N = Cp.getSize(), NEdges = Cp.getNEdges(); // convenience definitions
  double  p = Cp.getDensity();
  int     a, b, aOld, bOld;  // store previous step values
  double  prevNMIVal = 0.0, nmiVal, miVal, viVal;

  // Which direction are we stepping gammaAPM (or gammaRB)?  If the user 
  // specifies gammaStop<gammaStart, we step down in gamma.  Additionally, 
  // the run is truncated if we detect a fully collapsed system.
  bool bgammaStepDown = ( gammaStart > gammaStop );

  const int NDisplayLimit = 1001; // max number of nodes for console output
  const int NFileLimit = 201;     // max number of items for file output
  int       NReplicaPairs = NReplicasMax*(NReplicasMax-1)/2;
  //int  nMerges, nThisLoop, nMergeLoops;
  // declare some statistics variables
  double  solveTime;           // temporary storage for the time elapsed to 
  double  nIterAvg, nMergesD;  // stat variables on community solve function
  double  gammaRBMultiplier = 1.0, gammaAPM, gammaRB;
  bool    bOutputData = 1;     // default state is to output MRA data file
  if(outputFileName=="") {
    bOutputData = 0;           // overide data output 
    if(verbosity>1)  
      cout << red << "Warning:  Output filename was a null string in rnMRA()."
           << "  Not writing data file, but the best cluster list is still set."
           << normText << "\n";
    //errorMsg("Invalid outputFileName \"" + outputFileName + "\" in rnMRA().");
  } // end if outputFileName
  else
    outputFileName += ".csv";

  // on the log step size, we step down in gamma so that the collapse and 
  // disjoint cutoff checks still work
  bool   useLogStep = 1;
  double gammaLogStepSize = pow(10.0, 1.0/(double)NDecadeSteps );

  // density parameters - are not used for the log step general version
  bool useDensityStep = 0;  
  double densityStart = 0.95, densityStep = 0.05, densityScale = 0.5;
  double density, densityThreshold = 0.1;

  // catch some errors
  if(useZeroEMoves && !useEnergy && verbosity>1)
    warningMsg("Zero energy moves are not allowed for optimizing RB energy in rnMRA().");

  // open the output file
  ofstream  fout;
  if(bOutputData)  fout.open(outputFileName.c_str(),ios::out | ios::app);
  // iteration parameters

  // -------------------------------------------------------------------------
  // Plateau detection variables
  // moving average vector - we store the NMI values in the slot modulo the 
  // and use a crude threshold to identify exclude spikes and emphasize stable
  // regions of the graph
  int nMovingAvgLength = 3;
  TVectorFloat movingNMIs(nMovingAvgLength);
  movingNMIs.init(0.0);
  // track best NMI or VI trials
  double  maxNMI = 0.0, minVI = log((double)N)/log(2.0),
          prevNMI, prevMI, prevInfo, SpikeFilterW = 0.8, bestGamma = -1.0;
  int     bestHCount = 0, currentHCount = 0, iBest = -1;
  // -------------------------------------------------------------------------

  //Cp.display();  // debugging
  cout << cyan << "Starting rnMRA() function with a maximum number of steps of "
       << NStepsMax << " and gamma = " << gammaStart << " " << lightcyan;
  if(useDensityStep)   cout << "using a gamma density step";
  else if(useLogStep)  
    cout << "using a gamma log step with " << NDecadeSteps << " per decade."
         << "  gammaStart = " << gammaStart << " and gammaStop = " << gammaStop;
  cout << normText << "\n";

  // output run information
  string outputHeader;
  outputHeader  = "\"MRA data:\",\"N = " + itos(N)   + "\","
                + "\"L = " + itos( Cp.getL()       ) + "\","
                + "\"p = " + ftos( Cp.getDensity() ) + "\",";
  //if(!useEnergy)  outputHeader += "\"rAFG = " + ftos(r) + "\",";
  if(!useEnergy)  outputHeader += "\"rAFG = 0.0\",";
  //outputHeader += "\"NGraphsMax = " + itos(NGraphsMax) + "\",";
  outputHeader += "\"NReplicas = " + itos(NReplicasMax)  + "\","
                + "\"NTrials = "   + itos(NTrialsMax)    + "\","
                + "\"NIterMax = "  + itos(NIterMax)      + "\"\n";
  if(Cp.description!="")  outputHeader += Cp.description;
  if(useEnergy) {
    outputHeader += "\"Optimizing energy ";
    if(useZeroEMoves)  outputHeader += "with zero energy moves ";
    else               outputHeader += "without zero energy moves ";
  } else               outputHeader += "\"Optimizing modularity ";
  if(useMerges)        outputHeader += "using merges with a ";
  else                 outputHeader += "not using merges with a";
  outputHeader += initMethod.initString(2) + "\"\n";
  outputHeader += csvDateTime() + csvSeedString();
  if(sRunDescription!="")  outputHeader += "\"" + sRunDescription + "\"\n";

  // output column data identifiers, note:  gamma = -b/a
  outputHeader += csvH("a") + csvH("b") + csvH("alpha");
  if(useEnergy)  outputHeader += csvH("gamma");     // use our gamma
  else           outputHeader += csvH("gamma_RB");  // use RB gamma
  //outputHeader += "\"pout\",\"poutr\",\"sigma_poutr\",";
  //outputHeader += "\"Z_in\",\"sigma_Zin\",\"Z_out\",\"sigma_Zout\",";
  //outputHeader += "\"Z\",\"sigma_Z\",\"nMoved\",\"sigma_nMoved\",";
  outputHeader += csvSH("NMI")   + csvSH("VI")    + csvSH("I")     + csvSH("H") 
                + csvSH("nDisp") + csvSH("q")     + csvSH("E")     + csvSH("Q")
                + csvSH("Z")     + csvSH("Zin")   + csvSH("Zout")
                + csvSH("p")     + csvSH("pin")   + csvSH("pout")  + csvSH("mu")
                + csvH("q_best") + csvH("E_best") + csvH("Q_bestE") 
                + csvSH("t")     + csvSH("nIter") + csvSH("nMerges")
                //+ csvH("E<E_best Count") + csvH("Q>Q_bestE Count")
                + csvHeol("cluster sizes...");

  // now write the header string to the output file
  if(bOutputData)  fout << outputHeader;  

  // allocate memory and initialize the subH starting point
  TClusterListArray  replicas(NReplicasMax); // replicas lists
  for(int i=0; i<NReplicasMax; i++)  replicas[i].createSymmetric(N,1);
  TClusterListArray  allBest(NStepsMax,0); // track the strongest partitions
  TClusterList       working(N);  // working cluster list for community solver
  working.createSymmetric(N,1);

  // track all NMI values
  TVectorFloat  storedNMI(NStepsMax);
  storedNMI.init(-1);         // init to invalid NMI values
  storedNMI.resize(0,0,0,1);  // start with size 0 and add later
  // declare statistics variables
  TStats1D  nmiVals("NMI"),     miVals("I"),     viVals("VI"),  infoVals("H"),
            movedVals("nDisp"), qVals("q"), 
            energyVals("E"),    modVals("Q"),    iterVals("iterations"),  
            tVals("t"),         mergeVals("nMerges"),           
            zVals("Z"),         zinVals("Zin"),  zoutVals("Zout"),
            pVals("p"),         pinVals("pin"),  poutVals("pout"), muVals("mu");
  // track optimal communities
  int   outputCount = 0;  // output values for aOut and bOut manual displays
  int   iMinE, minE, minTrialE, nMoved, nIter, nIterR;
  bool  bDisjoint = 0, bDisjointAll = 1, bDisjointAny = 0, bDisjointStop = 0;

  // **************************************************************************
  // begin main resolution loop
  int i = 0;
  // for possible density iteration
  //cout << normText << "a = " << a << ", b = " << b << "\n";  // debugging
  //if(useDensityStep) { density  = densityStart;  gammaDensity(a,b,density); }

  //gammaAPM = gamma(a,b);
  //gammaRB  = gammaAPM*gammaRBMultiplier;  // scale for RB modularity
  gammaAPM = gammaStart;
  gammaRB  = gammaAPM*gammaRBMultiplier;  // scale for RB modularity
  //cout << normText << "a = " << a << ", b = " << b << "\n";  // debugging

  // crude starting aOld and bOld to make things easier.  This assumes that the
  // original Cp connection matrix starts as unscaled.
  aOld = 1;  bOld = -1;
  gammaToab(gammaAPM,a,b,verbosity);

  do {
    if(verbosity>0) {
      cout << normText << "Beginning loop " << i << " with ";
      if(useEnergy)  cout << magenta << "gamma = "   << gammaAPM;
      else           cout << magenta << "gammaRB = " << gammaRB;
      cout << normText << endl;
      //debugPause("gammaAPM = "+ftos(gammaAPM)+"... ");  
    } // end if

    // scale the existing connection matrix
    Cp.scale(a,b,aOld,bOld);    //Cp.display(); // debugging

    if(i>0) { prevMI = miVals.avg();  prevNMI = nmiVals.avg(); }
    //prevInfo = infoRs.avg();  prevq = qRs.avg();  prevgamma = gamma(a,b);  

    // clear statistics for next set of replicas
    nmiVals.clear();  miVals.clear();     viVals.clear();   movedVals.clear();
    qVals.clear();    energyVals.clear(); modVals.clear();
    infoVals.clear(); mergeVals.clear();  iterVals.clear(); tVals.clear();
    zVals.clear();    zinVals.clear();    zoutVals.clear();
    pVals.clear();    pinVals.clear();    poutVals.clear(); muVals.clear();
    bDisjoint = 0;    bDisjointAll = 1;   bDisjointAny = 0;
    // ------------------------------------------------------------------------
    // calculate replica solutions and (some) data storage
    localInitTime();
    minE = BigInteger;  // track lowest energy replica
    for(int j=0; j<NReplicasMax; j++) {
      if(verbosity>0)  msg("Starting replica "+itos(j)+"... ");

      // ----------------------------------------------------------------------
      // solve the current system
      if(verbosity>2)  cout << brown << "Starting trials... " << flush;
      if(useEnergy) // use our Potts model energy optimization
        //communityBestE(Cp,working,replicas[j],nIterAvg,nMergesD,solveTime,a,b,
        //  Symmetric,-1,NTrialsMax,NIterMax,useMerges,useZeroEMoves,0,verbosity);
        communityBestE(Cp,replicas[j],nIterAvg,nMergesD,solveTime,gammaAPM,
          initMethod,NTrialsMax,NIterMax,useMerges,useZeroEMoves,0,verbosity);
      //else          // else use modularity, RB Potts model, or AFG optimization
      //  communityBestQ(Cp,working,replicas[j],nIterAvg,nMergesD,solveTime,
      //    gammaRB,0.0,Symmetric,-1,NTrialsMax,NIterMax,useMerges,verbosity);
      //replicas[j].calcEnergy(Cp);  // already done in communityBest...
      if(!useEnergy)   replicas[j].calcQ(Cp,gammaRB);
      if(verbosity>2)  cout << brown << "done with trials." << normText << endl;

      iterVals.store(nIterAvg);   tVals.store(solveTime);  
      mergeVals.store(nMergesD);
      //loopVals.store(nMergeLoops/(double)NTrialsMax);
      energyVals.store(replicas[j].getE());  modVals.store(replicas[j].getQ());
      if(replicas[j].getE()<minE) { minE = replicas[j].getE();  iMinE = j; }
      qVals.store(replicas[j].getq());
      //infoVals.store(infoEntropy(replicas[j],N)); // done below in correlations
      // ----------------------------------------------------------------------

      if(verbosity>1)  msg("disjoint check... ");// debugging
      if(!bDisjointStop) {
        bDisjoint     = detectDisjoint(replicas[j],Cp);
        bDisjointAll &= bDisjoint;  bDisjointAny |= bDisjoint;
      } // end if bDisjointStop
      if(verbosity>1)  msg("done\n",grey);// debugging
    } // end for j
    // end of replica solutions

    // If we have a complete set of disjoint clusters, we are done; 
    // but we allow a few extra iterations for presentation purposes.
    if(bDisjointAll && !bDisjointStop) {
      warningMsg("Disjoint lists detected.  Truncating run soon.");
      NStepsMax = i + 4;  bDisjointStop = 1;
    } // end if
    // ------------------------------------------------------------------------
/*
    if(a==100 && b==-100) {
      //clusters.display(nodeOffset,1,"sorted clusters ");
      nodeStats(clusters,Cp,a,b,"hier256p50");
      clusterStats(clusters,Cp,a,b,"hier256p50");
    } else if(a==566 && b==-100) {
      //clusters.display(nodeOffset,1,"sorted clusters ");
      nodeStats(clusters,Cp,a,b,"hier256p30");
      clusterStats(clusters,Cp,a,b,"hier256p30");
    }
*/
    // calculate the replica information correlations
    nmviMRA(replicas,N,NReplicasMax,nmiVals,viVals,miVals,infoVals,movedVals,verbosity);
    if(verbosity>0)  cout << grey << " done." << normText << endl;
    pzStatsMRA(Cp,replicas,N,NReplicasMax,zVals,zinVals,zoutVals,
               pVals,pinVals,poutVals,muVals);
    // ------------------------------------------------------------------------

    // ----------------------------------------------------------------
    // plateau detection
    // set the single best trial for this mu - used for accuracy comparison
    // check whether this is the best NMI resolution so far (for this mu)

    // store current NMI value in loop mod moving average length
    // we then use the moving average below to exclude spikes in NMI
    movingNMIs[i%nMovingAvgLength] = nmiVals.avg();

    // Check for machine precision errors when NMI is "exactly" equal to 1
    // on a perfect solution.
    if(nmiVals.avg()>1.0) {
      cout << cyan << "Warning:  NMI was detected within epsilon larger "
           << "than 1 (NMI = " << setprecision(16) << nmiVals.avg() << ")?";
      // check for a genuine problem if NMI > 1.0+epsilon
      if(nmiVals.avg()>(1.0+epsilon)) 
        errorMsg(" and it is greater than (1+epsilon)");
      cout << "\n\n" << normText;
    } // end if NMI upper range check

    // Test for the maximum NMI.  Otherwise, else checks for a plateau in H.
    //if(nmiVals.avg()>(maxNMI+epsilon)) { // original selection criterion
    if( nmiVals.avg()>(maxNMI+epsilon) && 
       // add a special condition for application to RB Potts model
       (qVals.avg()<0.9*(double)N || useEnergy) ) {
      // add a rough threshold to exclude spikes and focus on the stable 
      // region.  The selection criterion is still the max NMI, we just 
      // skip the assignment if the threshold is not exceeded.
      if(movingNMIs.getAvg()>SpikeFilterW*maxNMI || i<nMovingAvgLength) {
        // store the determined best list based on the best NMI (so far)
        // only copy best replica if this is not a NMI "spike"
        best = replicas[iMinE];             // copy lowest-E replica to best
        maxNMI = nmiVals.avg();             // reset best NMI obtained
        bestGamma = gammaAPM;
        //sNMIBest = sStep;                 // store run description string
        iBest = i;
      } // end if moving average check
      // reset plateau variables for either a spike or a new plateau
      currentHCount = 1;  bestHCount = 1;
    } else {
      // try to detect the best local "plateau" in average Info H(A)
      // using an "H-count" - That is, track how broad the plateau is 
      // by recording the number or width of the plateau in H.
      // If the distance between points is widely varied, we would need a 
      // different way to track the "best" plateau(s) than counting number 
      // of data points in the plateau. 
      // This works for the LFR test with the current density iteration.
      if(nmiVals.avg()>(maxNMI-epsilon) && floatEq(nmiVals.avg(),prevNMI) 
         // add a special condition for application to RB Potts model
         //&& (qVals.avg()<((double)N-epsilon) || useEnergy) 
         && (qVals.avg()<0.9*(double)N || useEnergy) 
        ) {
        if(floatEq(miVals.avg(),prevMI))
          currentHCount++;  // increment current plateau "length"
        else  currentHCount = 1; // else reset the plateau counter on new H

        // store the best list determined by secondary plateau criteria
        // at an equal maximum NMI.  
        // This copy is not efficient since it recopies best for almost 
        // every step of the plateau in H, but it is easier for the moment.
        if(currentHCount>bestHCount) {
          bestHCount = currentHCount;
          best = replicas[iMinE];
          bestGamma = gammaAPM;
          //sNMIBest = sStep;
          iBest = i;
        } // end if
      } // end if plateau check
      else  currentHCount = 1; // else reset the plateau counter on new NMI
    } // end else max NMI check
/*
    // debugging
    //if(verbosity>3)
      cout << setprecision(6) //<< setfill(' ') << setw(8)
           << red << "Moving avg (i = " << setfill(' ') << setw(2) << i 
           << ") = " << setfill(' ') << setw(8) << movingNMIs.getAvg() 
           << brown << " (threshold = " 
           << setfill(' ') << setw(8) << ( 0.8*movingNMIs.getAvg() )
           << ")" << cyan << " w/ max NMI = " 
           << setfill(' ') << setw(8) << maxNMI << blue 
           << " (plateau = " << cyan << setfill(' ') << setw(2) << bestHCount 
           << blue << ")" << " iBest = " 
           << cyan << setfill(' ') << setw(2) << iBest
           << normText << " w/ current NMI avg = " 
           << setfill(' ') << setw(8) << nmiVals.avg()
           << ", H avg = " << setfill(' ') << setw(8) << infoVals.avg()
           << ", I avg = " << setfill(' ') << setw(8) << miVals.avg() 
           << " and plateau = " << currentHCount << "\n";
    debugPause("checking the plateaus...");
    // debugging
*/
    // ----------------------------------------------------------------
    // Now for VI, store current best VI result although the LFR code 
    // specifically focuses on NMI as used above.  This would be used for a
    // separate output file based on the best VI results.
    if(viVals.avg()<(minVI-epsilon) && best.getq()>1) {
      //sVIBest = sStep;  
      minVI = viVals.avg(); 
    } // end if
    // ------------------------------------------------------------------


    // ------------------------------------------------------------------------
    // peak detection logic - display and/or store very high correlation values
    double  bestDisplayThreshold = 0.99, bestStoreThreshold = 0.99;
    if(nmiVals.avg()>bestDisplayThreshold && N<NDisplayLimit && verbosity>1)
      replicas[iMinE].display(nodeOffset,1,"Best replica answer ");
    // if we have an optimal assignment and it is not a duplicate, store list
    //if(nmiVals.avg()>(1.0-epsilon)) {
    if(nmiVals.avg()>bestStoreThreshold) {
      double nmiPast;
      // compare this best with the last one
      if(allBest.getSize()==0)  allBest.add(replicas[iMinE]);
      else {
        if(allBest.getSize()>0)
         calcNMVILL(replicas[iMinE],allBest[allBest.getSize()-1],N,nmiPast,
                    miVal,viVal,nMoved);
        if(nmiPast<(1.0-epsilon))  allBest.add(replicas[iMinE]);
        else  cout << "Skipped best since it is too close to previous\n";
      } // end else
    } // end best stored
    // end peak detection logic
    // ------------------------------------------------------------------------

    storedNMI.add(nmiVals.avg());  // append NMI average to the end of the list
    cout << cyan << "NMI = "     << red     << nmiVals.avg()
         << cyan << ", VI = "    << magenta << viVals.avg()
         << cyan << ", and Q = " << magenta << modVals.avg()
         << cyan << " for loop " << i << " with a = " << a << " and b = " << b 
         << " with average q = " << qVals.avg() << ", " 
         << " p_min = " << ((double)b/(double)(b-a)) 
         << " and " << " nMerges = " << mergeVals.avg() << normText << "\n";

    // output this trial's parameters to the output file
    if(bOutputData)  
      fout << csv( a ) << csv( b ) << csv( alpha(a,b) ) << csv( gamma(a,b) )
           << nmiVals.csv()   << viVals.csv() << miVals.csv()    << infoVals.csv()
           << movedVals.csv() << qVals.csv()  << energyVals.csv() << modVals.csv()
           << zVals.csv() << zinVals.csv() << zoutVals.csv() 
           << pVals.csv() << pinVals.csv() << poutVals.csv() << muVals.csv() 
           << csv( replicas[iMinE].getq() ) << csv( replicas[iMinE].getE() ) 
           << csv( replicas[iMinE].getQ() ) 
           << tVals.csv() << iterVals.csv() << mergeVals.csv(); 
    // output cluster sizes with a practical column limit (esp. in OpenOffice)
    replicas[iMinE].sortClusters();
    if(bOutputData) {
      if(replicas[iMinE].getq()<NFileLimit)  
        for(int k=0; k<replicas[iMinE].getq()-1; k++)
          fout << csv( replicas[iMinE][k].getn() );
      fout << csveol( replicas[iMinE][replicas[iMinE].getq()-1].getn() ) << flush;
    } // end if output data

    if(verbosity>3) {
      cout << brown << "\nReplicas:   " << normText;  nmiVals.output(2);
      cout << brown << "VIs:        "   << normText;  viVals.output(1);
      cout << brown << "infos:      "   << normText;  infoVals.output(1);
      cout << brown << "nClusters:  "   << normText;  qVals.output(2);
      cout << brown << "energies:   "   << normText;  energyVals.output(2);
    } // end if verbosity

    cout << "Done with iteration " << i << endl;
    
    // get next iterations correct weights
    // track old values for use with Cp scale at the beginning of the loop
    aOld = a;  bOld = b;  
    if(bgammaStepDown)  gammaAPM /= gammaLogStepSize;
    else                gammaAPM *= gammaLogStepSize;
    gammaToab(gammaAPM,a,b,verbosity);
    //if(useLogStep)  gammaLogStepMRA(a,b,gammaAPM,gammaLogStepSize);
    //else if(useDensityStep)  
    //  gammaDensityStepMRA(a,b,gammaAPM,density,densityStep,densityScale,densityThreshold);
    //else errorMsg("Invalid weight step in rnMRA()");
    // adjust scale for RB modularity
    gammaRB = gammaAPM*gammaRBMultiplier;

    //cout << "End step " << i << endl;  // debugging
    i++;
    // conditions to end multiresolution iteration steps
  } while( i<NStepsMax && (                               // failsafe step max
          // fully collapsed?         
          (gammaAPM>gammaStop && qVals.avg()>(1.0+epsilon) && bgammaStepDown) ||
          (gammaAPM<gammaStop && !bgammaStepDown) ));     // gamma constraints
  // ((useDensityStep && qVals.avg()>1) || (useLogStep && gammaAPM>gammaMin)));  
  // end do while i - main comparison loop

  if(bOutputData)  fout << "\n";  // end the run data with an empty line
  // end main resolution loop
  // **************************************************************************
  if(verbosity>1) {
    if(allBest.getSize()>0) {
      //best.display(nodeOffset,1,"Best list ");
      cout << "Best " << allBest.getSize() << " cluster lists " 
           << "(with perfect correlations) were:\n";
      for(int k=0; k<allBest.getSize(); k++)
        if(N<NDisplayLimit) allBest[k].display(nodeOffset,1,"Best "+itos(k)+" ");
    } // end if exists any "best" lists
    else  cout << grey << " (with no perfect matches)" << normText << "\n";
  } // end if verbosity

  cout << "Stored NMI values:  ";
  TVectorOut(cout,storedNMI,1,1,1);
  cout << endl;

  if(bOutputData)  fout.close();
  // return the value of gamma corresponding to the best answer
  return bestGamma;  
} // end new hierarchy by replica comparisons

