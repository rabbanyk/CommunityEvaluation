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
// ----------------------------------------------------------------------------
inline double kMinA2(double x, double xMean, double xf) {
  return ( log(xf) + xMean/xf - log(x) - xMean/x );
} // end kMinA2
double getkMinA2(double xnp1, double xn, double xMean, double xf,
                 int NMaxIterations = 100, double epsilon = 0.5E-10) {
  // from Wikipedia Secant Method
  double d;
  for(int n=1; n<=NMaxIterations; n++) {
    d = (xn - xnp1) / 
        (kMinA2(xn,xMean,xf) - kMinA2(xnp1,xMean,xf)) * kMinA2(xn,xMean,xf);
    if(fabs(d)<epsilon)  return xn;
    xnp1 = xn;
    xn  -= d;
  } // end for i
  return xn;
} // end getkMinA2


double noiseTestPowerRS(string outputFileNameRoot, int N, double kMax = 64.0, 
  int nMin = 8, int nMax = 24, double gammaLFR = -2.0,  double betaLFR = -1.0,
  double pin = 0.95, 
  bool useEnergy = 0, bool useZeroEMoves = 1, bool useMerges = 1,
  int  IterationType = 0,
  double r = 0.0, int a = 1, int b = -1, int NDecadeSteps = 20,
  double kStart=10, double kStop=50, double kStep=2, 
  int NTrialsMax=4, int NGraphsMax=100,  int NIterMax=100, 
  int NStepsMax=100, int verbosity=0, string sRunDescription = "") {
  // pout refers to the level of added noise at each iteration.  
  // a - the starting value for the connected edge weighting
  // b - the starting value for the unconnected edge weighting
  // answer - the input 'kernel' cluster list, best if it is disjoint
  double  epsilon = 5.0e-8;
  double  nmiVal, miVal, viVal, pgnVal;
  double  gammaRB = -(double)b/(double)a, gammaAPM = gammaRB;
  int     iBest, minE;
  double  maxQ, gammaRBiBest;  // track max modularity
  double  pout, pinActual, poutActual, zinAvg,zoutAvg,zAvg,pinAvg,poutAvg,pAvg;
  int     q, nMoved, nMovedVI, nEdgesOut, nEdgesIn, nIter;
  double  NEdgesInMax, NEdgesOutMax;  // use double because max size is large
  //bool    useZeroEMoves = 1, 
  bool    useEvenInit = 1;
  int     aOld = a, bOld = b, gammaStep = 1, aStepSize, abStepSize = 1, am1Log10;
  int     aStart = a, bStart = b, aStartLinear = 30*abs(b);
  bool    useLogStepSize = 1;  // for weight iteration
  //bool    useMerges = 1;
  double  gammaLogStep, gammaAPMStop;  // for weight iteration
  int     qNMVICutoff = 3000;
  
  // -------------------------------------------------------------------------
  // which iteration type are we using?
  //   0 - step only kMin
  //   1 - step kMin and kMax
  //   2 - step only kMean
  int  kStepType = IterationType;
  if(kStepType==0)       cout << green << "  Incrementing only kMin";
  else if(kStepType==1)  cout << green << "  Incrementing kMin and kMax";
  else if(kStepType==2)  cout << green << "  Incrementing only kMean";
  else if(kStepType==3)  cout << green << "  Incrementing kMean and kMax";
  else  errorMsg("Degree step type is not known in noiseTestPowerRS()?");
  cout << " by steps of size " << kStep << " stopping at kMin (or kMean) = " 
       << kStop << normText << "\n";
  //debugPause("Stop?");
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // allow the use of the multiresolution (MRA) algorithm
  bool    useMRA = 0;
  double  gammaStart, gammaStop, bestGamma;
  if(useEnergy) { gammaStart = gamma(a,b);  gammaStop = 0.001; } // APM solver
  else          { gammaStart = gamma(a,b);  gammaStop = 0.1;   } // RBPM solver
  int     NReplicasMax = 8;   // internal only number of replicas
  int     NStepsMaxMRA = 300; // limit the number of steps in MRA solver
  string  outputFileNameMRA = outputFileNameRoot + "_MRAData.csv";
  TInitSymmetric  initMethod(N,1);  // use a symmetric init in MRA solver
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // the integer implementation messes up the step size, an ad hoc fix is to 
  // scale the integer weights a and b by a multiplier only when using the 
  // RB Potts model
  double  gammaRBMultiplier = 1.0;
  // use abStepSize as the number of steps per decade instead
  //gammaLogStep = pow(10.0, 1.0/(double)abStepSize);
  gammaLogStep = pow( 10.0, 1.0/(double)NDecadeSteps );
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // define the number of pout densities to be tested
  //int     NpSteps = (int)( (poutStop-poutStart)/poutStep ) + 1;
  double  solveTime, nMergesD, nIterAvg;   
  // track how many graph solutions actually result in a *lower* energy or 
  // *higher* Q (based on the current solver) than the purported "answer"
  int     lowerECount = 0, higherQCount = 0;  
  int     iterCount = 0;  // track number of iterations
  // -------------------------------------------------------------------------

  // error detection check(s)
  if(kStart<1)  warningMsg("kStart is too small in noiseTestPowerRS()");

  // here are some variables to track the maxNMI values for each weight loop
  string    outputHeader, sNMIBest, sVIBest, fname, sStep;
  double    maxNMI, minVI;
  ofstream  fout, nmiout, viout;
  fname = outputFileNameRoot + ".csv";
  fout.open(fname.c_str(),ios::out | ios::app);
  if(NStepsMax>1) { // special case for multiple weight samples
    fname = outputFileNameRoot + "BestNMI.csv";
    nmiout.open(fname.c_str(),ios::out | ios::app);
    fname = outputFileNameRoot + "BestVI.csv";
    viout.open(fname.c_str(),ios::out | ios::app);
  } // end if

  //const int NDisplayLimit = 1001;  // max number of nodes for console output
  //const int NFileLimit = 901;      // max number of items for file output

  cout << cyan << "Starting power noise test with a NGraphsMax = "
       << NGraphsMax << " and NTrialsMax = " << NTrialsMax << " with " << green;
  if(useEnergy)  
    cout << "energy" << cyan << " optimization" << normText << endl;
  else          
    cout << "modularity" << cyan << " optimization with " << NDecadeSteps 
         << " steps per decade of gammaRB (step scale is " << gammaLogStep 
         << ")" << normText << endl;

  //TCMSparseW  Cp(N);
  TCMDense      Cp(N);

  double kMean, kMin; //, kMax;
  int    nMerges, nIterR;

  TVectorFloat  storedNMI(NStepsMax,-1.0), storedVI(NStepsMax,-1.0);

  // -------------------------------------------------------------------------
  // the Init lists are for debugging
  TClusterList  answer(N,N), clusters(N,N), best(N,N), 
                clustersInit(N,N), bestInit(N,N);
  answer.createSymmetric(N);
  clusters.createSymmetric(N);      best.createSymmetric(N);
  clustersInit.createSymmetric(N);  bestInit.createSymmetric(N);
  //cout << "Done with list declarations " << endl;
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  //cout << "Declaring stat variables " << flush;
  TStats1D nmiVals("NMI"),  miVals("I"),     viVals("VI"),      qVals("q"),  
           movedVals("nDisplaced"),          mergeVals("nMerged"),  
           zVals("z"),      zinVals("zin"),  zoutVals("zout"),  muVals("mu"), 
           infoVals("H"),   pVals("p"),  pinVals("pout"), poutVals("pout"), 
           energyVals("E"), eRBVals("E_RB"), modRBVals("Q_RB"), modVals("Q"),
           tVals("t"),      iterVals("iterations"),
           eAnsVals("E_ans"),  eRBAnsVals("Q"), modAnsVals("Q_ans"),
           modRBAnsVals("Q_RB"), 
           pAnsVals("p_ans"),  pinAnsVals("pout_ans"),  poutAnsVals("pout_ans"),
           zAnsVals("z_ans"),       zinAnsVals("zin_ans"),  qAnsVals("q_ans"),  
           zoutAnsVals("zout_ans"), muAnsVals("mu_ans");
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // output run information header information and column labels
  outputHeader += "\"Comparison test for increasing power law noise:\",";
  outputHeader += cssv("N = ",N) //+ cssv("kMax = ",kMax) 
                + cssv("alpha_LFR = ",gammaLFR) + cssv("beta_LFR = ",betaLFR) 
                + cssv("nMin = ",nMin)          + cssv("nMax = ",nMax)
                + cssv("pin = ",pin)            //+ cssv("r_AFG_Start = ",r)
                + cssv("NGraphs = ",NGraphsMax) + cssv("NTrials = ",NTrialsMax)
                + cssv("NIterMax = ",NIterMax);
  if(useEnergy) {
    outputHeader += "\n\"Optimizing energy ";
    if(useZeroEMoves)  outputHeader += "with zero energy moves\"\n";
    else               outputHeader += "without zero energy moves ";
  } else
    outputHeader += cssv("NDecadeSteps = ",NDecadeSteps)
                  + "\n\"Optimizing RB configuration Potts model ";
  if(useMerges)          outputHeader += "using merges with a ";
  else                   outputHeader += "not using merges with a";
  // the following is a hack addition so that I do not have to immediately
  // modify the communityBest E or Q functions
  #ifdef USE_POWER_INIT
  outputHeader += "power law inititialization (same parameters as above).\"\n";
  #else
  outputHeader += "symmetric initialization.\"\n";
  #endif
  if(kStepType==0)       outputHeader += "\"Incrementing only kMin\"\n";
  else if(kStepType==1)  outputHeader += "\"Incrementing kMin and kMax\"\n";
  else if(kStepType==2)  outputHeader += "\"Incrementing only kMean\"\n";
  else if(kStepType==3)  outputHeader += "\"Incrementing kMean and kMax\"\n";
  outputHeader += csvDateTime() + csvSeedString();
  if(sRunDescription!="")  outputHeader += "\"" + sRunDescription + "\"\n";
  outputHeader += Cp.description;

  // output column data identifiers, note:  gamma = -b/a
  outputHeader += "\"a\",\"b\",";
  if(useEnergy)  outputHeader += "\"alpha\",\"gamma\",";
  else           outputHeader += "\"r_AFG\",\"gamma_RB\",";
  if(kStepType==0)       outputHeader += "\"kMin_out_alg\",\"kMean_out_calc\",";
  else if(kStepType==1)  outputHeader += "\"kMin_out_alg\",\"kMean_out_calc\",";
  else if(kStepType==2)  outputHeader += "\"kMin_out_calc\",\"kMean_out_alg\",";
  else if(kStepType==3)  outputHeader += "\"kMin_out_calc\",\"kMean_out_alg\",";
  outputHeader += "\"kMax_out_alg\",";
  // output statistics headers

  outputHeader += "\"Z\",\"sigma_Z\",\"Z_in\",\"sigma_Zin\",";
  outputHeader += "\"Z_out\",\"sigma_Zout\",\"mu\",\"sigma_mu\",";
  outputHeader += "\"NMI\",\"sigma_NMI\",\"VI\",\"sigma_VI\",\"I\",\"sigma_I\",";
  outputHeader += "\"H\",\"sigma_H\",\"q\",\"sigma_q\",";
  outputHeader += "\"nDisp\",\"sigma_nDisp\",";
  outputHeader += "\"p\",\"sigma_p\",\"pin\",\"sigma_pin\",\"pout\",\"sigma_pout\",";
  outputHeader += "\"E\",\"sigma_E\",\"E_RB\",\"sigma_E_RB\",";
  outputHeader += "\"Q_RB\",\"sigma_Q_RB\",\"Q\",\"sigma_Q\",";
  outputHeader += "\"t_avg\",\"t_sigma\",\"iter_avg\",\"iter_sigma\",";
  outputHeader += "\"nMerges\",\"sigma_nMerges\",";
  outputHeader += "\"E_ans\",\"sigma_E_ans\",\"E_RB_ans\",\"sigma_E_RB_ans\",";
  outputHeader += "\"Q_RB_ans\",\"sigma_Q_RB_ans\",\"Q_ans\",\"sigma_Q_ans\",";
  outputHeader += "\"p_ans\",\"sigma_p_ans\",\"pin_ans\",\"sigma_pin_ans\",";
  outputHeader += "\"pout_ans\",\"sigma_pout_ans\",";
  outputHeader += "\"Z_ans\",\"sigma_Z_ans\",\"Z_in_ans\",\"sigma_Zin_ans\",";
  outputHeader += "\"Z_out_ans\",\"sigma_Zout_ans\",\"mu_ans\",\"sigma_mu_ans\",";
  outputHeader += "\"q_ans\",\"sigma_q_ans\",";
/*
  // output statistics columns headers
  outputHeader += zVals.csvH()     + zinVals.csvH()    + zoutVals.csvH() 
                + muVals.csvH()
                + nmiVals.csvH()   + viVals.csvH()     + miVals.csvH() 
                + infoVals.csvH()  + qVals.csvH()      + movedVals.csvH()
                + pVals.csvH()     + pinVals.csvH()    + poutVals.csvH()
                + energyVals.csvH() + eRBVals.csvH()
                + modRBVals.csvH() + modVals.csvH() 
                + tVals.csvH()     + iterVals.csvH()   + mergeVals.csvH() 
                + eAnsVals.csvH()  + eRBAnsVals.csvH() + modRBAnsVals.csvH() 
                + modAnsVals.csvH() 
                + pAnsVals.csvH()  + pinAnsVals.csvH() + poutAnsVals.csvH()
                + zAnsVals.csvH()  + zinAnsVals.csvH() + zoutAnsVals.csvH() 
                + muAnsVals.csvH() + qAnsVals.csvH();
*/
  outputHeader += "\"E<E_ans count\",\"Q_RB>Q_RB_ans count\"\n";
  fout << outputHeader;  // now send header string to output file
  if(NStepsMax>1) { // special case for RB model
    nmiout << outputHeader;  viout << outputHeader;
  } // end if
  // -------------------------------------------------------------------------

  // **************************************************************************
  // degree iteration is based on kMin or kMean
  if(kStepType==0 || kStepType==1)       kMin  = kStart; 
  else if(kStepType==2 || kStepType==3)  kMean = kStart;
  else errorMsg("Degree step type is not known in noiseTestPowerRS()?");
  do {
    if(kStepType==2 || kStepType==3) {
      // we use kMin and kMax below to determine the power law noise, so when 
      // stepping kMean, get kMin from kMean and kMax
      if(!floatEq(gammaLFR,-2.0))  
        errorMsg("We assume gamma_LFR = -2 (degree power law), but gamma_LFR = "
                +ftos(gammaLFR)+" in noiseTestPowerRS()?");
      else {
        kMin = getkMinA2(kMean/2.5,kMean/2.0,kMean,kMax);
        if(isnan(kMin) || kMin<1.0)  kMin = 1.0;
        if(kMin<(1.0-epsilon) || kMin>kMean || kMin>(kMax+epsilon))
          errorMsg("Invalid value of kMin = "+ftos(kMin)+" (with kMean = "+
                   ftos(kMean)+" and kMax = "+ftos(kMax)+")?");
      } // end else
      // output staring info
      cout << brown << "Beginning kMin " << kMin << " (kMin stop is " << kStop
           << ") and kMax = " << kMax << "\n" << normText;
    } else {
      // otherwise calculate kMean directly from the power law relation
      kMean = randomPowerMean(kMin,kMax,gammaLFR);
      if(isnan(kMean) || kMean<kMin || kMean>kMax)
        errorMsg("Invalid calculated value kMean = "+ftos(kMean)
                +" in noiseTestPowerRS()?");
      // output staring info
      cout << brown << "Beginning kMean " << kMean << " (kMin stop is " << kStop
           << ") and kMax = " << kMax << "\n" << normText;
    } // end else step type

    // mark trial in the large data file if we vary the weight for each kMean 
    // (usually only for RB Potts model)
    if(NStepsMax>1)  fout << "\"kMean " << kMean << "\"\n";  

    // begin main iteration loop for model weights
    int i = 0;
    a = aStart;  b = bStart;
    gammaAPM = gamma(a,b);  gammaRB = gammaAPM*gammaRBMultiplier;
    // string storage for identifying the max NMI (for RB model in noise test)
    maxNMI = -1.0;   minVI = log((double)N)/log(2.0)+epsilon;
    sNMIBest = "Init NMI String";  sVIBest = "Init VI String";  sStep = "";

    //CpStep = CpStart;
    do {
      if(verbosity>0) {
        cout << normText << "\nBeginning Step " << i 
             << " with a = " << a << ", b = " << b << " gamma = " << gammaAPM 
             << ", and gammaRB = " << gammaRB << "\n";
        cout << normText << "Previous a = " << aOld << " and b = " << bOld 
             << "\n"; // debugging
      } // end if

      // clear statistics for next set of graphs
      nmiVals.clear();     miVals.clear();       viVals.clear();   
      movedVals.clear();   infoVals.clear();     qVals.clear();    
      energyVals.clear();  modVals.clear();    
      eRBVals.clear();     modRBVals.clear();
      pVals.clear();       pinVals.clear();      poutVals.clear();
      zinVals.clear();     zVals.clear();        zoutVals.clear(); 
      muVals.clear(); 
      tVals.clear();       iterVals.clear();     mergeVals.clear();
      eAnsVals.clear();    modAnsVals.clear(); 
      eRBAnsVals.clear();  modRBAnsVals.clear();
      pAnsVals.clear();    pinAnsVals.clear();   poutAnsVals.clear();
      zinAnsVals.clear();  zAnsVals.clear();     qAnsVals.clear(); 
      zoutAnsVals.clear(); muAnsVals.clear(); 
      lowerECount = 0;     higherQCount = 0,   // summing variables

      // ----------------------------------------------------------------------
      // calculate replica solutions and (some) data storage
      localInitTime();
      for(int j=0; j<NGraphsMax; j++) {
        if(verbosity>2)  cout << "Starting graph " << j << "... " << flush;
        // output marker for separate graphs in data file
        //if(NStepsMax>1)  fout << "\"Graph " << j << "\"\n";  

        // this is somewhat crude, but we initialize a densly connected "answer"
        // and then add power-law distributed edges using the established
        // clusters as the floor for the degree of each node.  The power-law
        // variates are still generated by kMin however.
        if(verbosity>2)  cout << "Init answer... " << flush;   // debugging
        answer.initPowerS(N,nMin,nMax,betaLFR,1);

        // set some statistical and parameter variables that depend on answer
        NEdgesOutMax = 0.0, NEdgesInMax = 0.0;  // use as a summing variable
        nEdgesOut = 0;
        for(int i=0; i<answer.getq(); i++)
          NEdgesInMax += answer[i].getnD()*(answer[i].getnD()-1.0)/2.0;
        NEdgesOutMax = (double)N*(double)(N-1)/2.0 - NEdgesInMax;
        //debugPause();  // debugging

        // now define a system with a power-law degree distribution and pin
        // as the minimum average internal edge density
        if(verbosity>2)  cout << "Init Cp... " << flush;   // debugging
        // init with zero noise and add noise separately
        //Cp.initByClusters2(answer,N,pin,0.0,nEdgesIn,nEdgesOut,1,1);
        //nEdgesOut = Cp.addNoisePower(answer,kMin,kMax,gammaLFR);
        Cp.initByClustersPower(answer,N,pin,kMin,kMax,gammaLFR,nEdgesIn,nEdgesOut);
        if(useEnergy)  Cp.scale(a,b,1,-1);  // assumes original is unscaled!
        //cout << magenta << "L = " << Cp.getL() << "\n";  // debugging

        #if 0
        // debugging 
        TStats1D  zDisplay(N,0,"Display degrees");
        for(int i=0; i<N; i++)  zDisplay.store((double)Cp.ki(i)-epsilon);
        cout << magenta << "\nDegree stats:  "; 
        zDisplay.output(3);  //zDisplay.binTest();
        debugPause("Degree stats output...");
        // debugging
        #endif

        // calculate some parameters for the known answer
        answer.calcEnergy(Cp);  // calculate answer energy given new matrix
        //cout << "Answer Q = " << answer.getQ() << endl;  // debugging
        if(!useEnergy)  answer.calcQ(Cp,gammaRB,r);
        answer.calcPZParams(Cp,zinAvg,zoutAvg,zAvg,pinAvg,poutAvg,pAvg);
        pAnsVals.store(pAvg);       pinAnsVals.store(pinAvg);  
        poutAnsVals.store(poutAvg);
        zinAnsVals.store(zinAvg);   zoutAnsVals.store(zoutAvg);
        zAnsVals.store(zAvg);       muAnsVals.store(zoutAvg/zAvg);
        qAnsVals.store(answer.getq());
        if(verbosity>0)
          cout << "Answer zin = " << zinAvg << " and zout = " << zoutAvg 
               << " with q = " << qAnsVals.avg()   << endl;  // debugging

        if(verbosity>0)
          cout << green << "p = " << Cp.getDensity() 
               << ".  nEdgesIn = " << nEdgesIn 
               << " NEdgesInMax = " << NEdgesInMax
               << ".  nEdgesOut = " << nEdgesOut 
               << " NEdgesOutMax = " << NEdgesOutMax
               << ".  pk_in = " << pinActual << normText << "\n";

        // -------------------------------------------------------------------
        // get the current solution
        if(useMRA)
          bestGamma = rnMRA(Cp,outputFileNameMRA,best,initMethod,
            gammaStart,gammaStop,NDecadeSteps,useEnergy,useMerges,useZeroEMoves,
            //NReplicasMax,NTrialsMax,NStepsMax,NIterMax,nodeOffset,verbosity);
            NReplicasMax,NTrialsMax,NStepsMaxMRA,NIterMax,0,verbosity-1,
            "Graph "+itos(j)+" MRA data in noiseTestPowerRS()");
        else if(useEnergy)
          communityBestE(Cp,clusters,best,nIterAvg,nMergesD,solveTime,a,b,
              Symmetric,-1,NTrialsMax,NIterMax,useMerges,useZeroEMoves,
              0,verbosity);
        //else  // else use modularity, RB Potts model, or AFG optimization
        //  communityBestQ(Cp,clusters,best,nIterAvg,nMergesD,solveTime,
        //      gammaRB,0.0,Symmetric,-1,NTrialsMax,NIterMax,useMerges,verbosity);
        // -------------------------------------------------------------------

        best.calcPZParams(Cp,zinAvg,zoutAvg,zAvg,pinAvg,poutAvg,pAvg);
        pVals.store(pAvg);  pinVals.store(pinAvg);  poutVals.store(poutAvg);
        zVals.store(zAvg);  zinVals.store(zinAvg);  zoutVals.store(zoutAvg);
        muVals.store(zoutAvg/zAvg);

        calcNMVILL(best,answer,N,nmiVal,miVal,viVal,nMovedVI);
        movedVals.store(nMovedVI);  
        nmiVals.store(nmiVal);  viVals.store(viVal);  miVals.store(miVal);

        energyVals.store(best.getE());  modVals.store(calcMod(best,Cp));
        eAnsVals.store(answer.getE());  modAnsVals.store(calcMod(answer,Cp));
        eRBVals.store(-best.getQ()*Cp.getLD()); 
        modRBVals.store(best.getQ());   modRBAnsVals.store(answer.getQ());   
        eRBAnsVals.store(-answer.getQ()*Cp.getLD());  
        qVals.store(best.getq());       infoVals.store(infoEntropy(best,N));
        tVals.store(solveTime);     // avg solve time
        iterVals.store(nIterAvg);   // avg iterations
        mergeVals.store(nMergesD);  // avg merges
        // track number of potentially ill-defined "answers"
        if(best.getE()<answer.getE())   lowerECount++;
        if(best.getQ()>answer.getQ())   higherQCount++;

        if(verbosity>1)  msg("done\n",grey);// debugging
      } // end for j - graphs iteration

      if(verbosity>1 && N<5000) {
        best.display(0,1,"best ");  cout << "\n"; // debugging
        answer.display(0,1,"answer ");   cout << "\n"; // debugging
      } // end if
      // end of solutions
      // ----------------------------------------------------------------------

      if(verbosity>2)  cout << "Getting update... " << endl;
      storedNMI[i] = nmiVals.avg();  storedVI[i] = viVals.avg();
      if(verbosity>0) 
        cout << blue << "NMI = "  << red     << nmiVals.avg()
             << blue << ", VI = " << green   << viVals.avg()
             << blue << ", Q = "  << magenta << modVals.avg()
             << ", q = " << qVals.avg() << normText << "\n";

      // store parameters for current step as a line of csv format string
      if(verbosity>2)  cout << "Storing parameters... " << endl;
      if(useMRA) {
        // use MRA returned bestGamma to calculate the connection matrix weights
        int aGamma, bGamma;
        gammaToab(bestGamma,aGamma,bGamma);
        sStep = csv( aGamma ) + csv( bGamma );
        if(useEnergy)  sStep += csv( alpha(aGamma,bGamma) ) + csv( bestGamma );  
        else           sStep += csv( r )                    + csv( bestGamma );
      } else {
        // output data as normal with internal variable values
        sStep = csv( a ) + csv( b );
        if(useEnergy)  sStep += csv( alpha(a,b) ) + csv( gammaAPM );  
        else           sStep += csv( r )          + csv( gammaRB  );
      } // end else
      
      //if(verbosity>2)  cout << "a... " << flush;  // debugging

      sStep += csv(kMin)        + csv( kMean )     + csv( kMax )
             + zVals.csv()      + zinVals.csv()    + zoutVals.csv()
             + muVals.csv()     + nmiVals.csv()    + viVals.csv()
             + miVals.csv()     + infoVals.csv()   
             + qVals.csv()      + movedVals.csv()  
             + pVals.csv()      + pinVals.csv()    + poutVals.csv()
             + energyVals.csv() + eRBVals.csv()    
             + modRBVals.csv()  + modVals.csv()    
             + tVals.csv()      + iterVals.csv()   + mergeVals.csv()
             + eAnsVals.csv()   + eRBAnsVals.csv() + modRBAnsVals.csv()
             + modAnsVals.csv() 
             + pAnsVals.csv()   + pinAnsVals.csv() + poutAnsVals.csv()
             + zAnsVals.csv()   + zinAnsVals.csv() + zoutAnsVals.csv()
             + muAnsVals.csv()  + qAnsVals.csv()
             + csv( lowerECount ) + csveol( higherQCount );

      if(verbosity>2)  cout << "Writing basic data... " << endl;
      fout << sStep << flush; // output information for this step

      // store current best NMI results
      if(nmiVals.avg()>maxNMI) {
        sNMIBest = sStep;
        maxNMI   = nmiVals.avg();
        iBest    = i;
        gammaRBiBest = gammaRB;
      } // end if
      // store current best VI result
      if(viVals.avg()<minVI && best.getq()>1) {
        sVIBest = sStep;
        minVI   = viVals.avg();
      } // end if

      // ------------------------------------------------------------------
      // now adjust weight according to specify step type (constant or log)
      if(verbosity>1)  cout << "Adjusting weights\n" << endl;
      aOld = a;  bOld = b;       // store for scale operation

      if(useLogStepSize && NStepsMax>1) {
        gammaAPM *= gammaLogStep;
        if(verbosity>2)  cout << red << "Resetting and b:   ";
        // now calculate the approximate a and b for our solvers
        if(gammaAPM>1.0) { a =  100;  b = -(int)( gammaAPM*100.0 + 0.5 ); }
        else             { b = -100;  a =  (int)( 100.0/gammaAPM + 0.5 ); }
        if(verbosity>1)
          cout << red << ".  The new a = " << a << " and b = " << b << endl;

        gammaRB = gammaAPM*gammaRBMultiplier;
      } // end if useLogStepSize
      else if(NStepsMax>1) {
        gammaStep = abStepSize;  // use regular step size
        if(a>abs(bStart))  a -= gammaStep;
        else               b -= gammaStep;
        gammaRB = gamma(a,b)*gammaRBMultiplier;
      } // end if NStepsMax
      // catch invalid cases for a or b
      if(b==0 || a ==0) 
        errorMsg("Invalid zero value for a or b in noiseTestPowerRS()");
      // ------------------------------------------------------------------
      
      if(verbosity>1)  cout << "Done with iteration " << i << "\n" << endl;
      i++;
    } while(i<NStepsMax);

    cout << "Stored NMI values: "; TVectorOut(cout,storedNMI,1,1); cout << "\n";
    cout << "Stored VI values:  "; TVectorOut(cout,storedVI,1,1);  cout << "\n";
    cout << endl;

    if(NStepsMax>1) {
      fout << endl; // differentiate different k's in large file
      nmiout << sNMIBest << flush;  viout << sVIBest << flush;
    } // end special case for RB model

    if(kStepType==0)        kMin  += kStep;                 
    else if(kStepType==1) { kMin  += kStep;  kMax += kStep; }
    else if(kStepType==2)   kMean += kStep;
    else if(kStepType==3) { kMean += kStep;  kMax += kStep; }
    else errorMsg("Degree step type is not known in noiseTestPowerRS()?");
    
  } while( (kMin <(kStop+epsilon) && (kStepType==0 || kStepType==1)) ||
           (kMean<(kStop+epsilon) && (kStepType==2 || kStepType==3)) );
  //} while( pout<(poutStop+epsilon));
  // end do while i - main comparison loop

  // output the \gamma + \delta\gamma results for the RB model
  if(!useEnergy) {
    
  } // end if

  // if I end the output with an endl, I do not get recurring memory error?
  fout << endl;  // end the run data with an empty line
  fout.close();  nmiout.close();  viout.close();
  // end main loop
  // **************************************************************************

  return 0.0;
} // end noiseTestPowerRS
