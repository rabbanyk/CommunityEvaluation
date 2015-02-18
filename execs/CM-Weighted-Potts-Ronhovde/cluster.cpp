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
// standard includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <sys/times.h>

// local includes

#ifndef MSL_ARRAY_H
#include "MSL/MSL_Array.h"
#include "MSL/MSL_Array1D_Template.h"
#include "MSL/MSL_Array2D_Template.h"
#endif
#include "clusterclasses.h"
//#ifndef _ML_UTILS
#include "ML_Utils.h"
#include "MSL_Stats1D.h"
//#endif
#include "cluster.hpp"
#include "mytimefns.cpp"

using namespace std;
using namespace MSL;

// declare variables and program parameters
RandomLib::Random r;
int     NStepMax = 1000, NTrialsMax = 1, NRunsMax = 1;
int     NStepStartCheck = 100;    // step number to begin periodic checks
int     nStep;
int     CheckInterval = 100;      // step interval which to check convergence
int     debugMode = 0;            // debug level:  0 is none
unsigned long long ISeed = 0;     // random number seed
bool    flag = 0;
double  accSum = 0.0, acc2Sum = 0.0, accAvg;
//const   int BigInteger = 2 << 29; // constant large integer
int     verbosity = 0;            // how much non-debugging info is reported

// algorithm parameters
int     N = 64, NEdges = 0;       // starting number of spins and edges
double  p = 0.2, pin = 0.2, pout = 0.05; // connection probability
int     weightC = 1, weightU = -1;// weights of connected and unconnected edges
int     spinOffset = 0;           // offset for displaying the answer
double  epsilon = 0.5e-6;         // error tolerance on float equals tests
int     QMax = 20, q = 6;         // the number of communities
double  zout = -1.0;              // for fixed valency per spin, what is the 
                                  // avg number of spins outside the community?

bool   useFixedValency = 0;       // use an input data file?
//bool useRestart      = 0;       // use a restart file? - invalid right now
bool   useInputFile    = 0;       // do we use an input data file?
bool   exitwNoRun      = 0;       // exit just prior to main calculations?
bool   outputCMatrix   = 0;       // output CMatrix to file?
bool   onlyCMatrix     = 0;       // write CMatrix file and exit?
bool   showCMatrix     = 0;       // display version of CMatrix to console?
bool   useEnergy       = 1;       // use energy for convergence criteria?
                                  //   (otherwise use modularity for testing)
bool   useERGraph      = 0;       // use energy for convergence criteria?
                                  //   (otherwise use modularity for testing)
bool   useEqualsTest   = 0;       // use equals test for known systems

string inputFileName(""), outputFileName(""), outputFileRoot("");

// for debugging
int    outputSpinConnections = -1;// -1 flag value ignores at end, >0 is valid


// *** Begin MAIN ************************************************************
int main(int argc, char *argv[]) {
  cout << normText << boldText << "\ncluster.cpp code output\n";

  cout << normText << "Echo command line:  " << greyText << "\"" << argv[0];
  for(int i=1; i<argc; i++) { cout << " " << argv[i]; } // end command echo
  cout << "\"\n\n" << normText ;

  // initialize starting time variables
  initTimes(StartTime, theTime);
  // interpret command line parameters
#include "cluster_clparams.cpp"

  // -------------------------------------------------------------------------
  // make variable and parameter changes specified by user
  // -------------------------------------------------------------------------
  // setup system if a fixed outside valency is specified by user
  // At the moment we assume a total number of connect edges per spin of 16 
  // since that is the standard used in the comparison graphs in the literature
  if(useFixedValency) {
    warningMsg("  User specified to use a fixed valency.");
    warningMsg("    -pout and -pin flags are ignored for fixed valency.",greyText);
    double n = (double)N/(double)QMax;
    pin  = (16.0 - zout)/(n - 1.0);
    pout = (16.0 - (n-1.0)*pin)/((double)N-n);
    cout << normText << "    Fixed valency values pin = " << pin 
                     << " and pout = " << pout << "\n";
  } // end if useFixedValency

  // set random seed either by user value or by class initialization
  //srand(ISeed);
  cout << "  Maximum number of time steps is " << NStepMax << "\n";
  if(ISeed>0) {
    r.Reseed(ISeed);
    cout << greenText << "User initialized RandomLib seed to: "
         << brownText << r.SeedString() << normText << endl;
  } // end if Iseed
  else  cout << greenText << "RandomLib seed automatically set to: "
             << r.SeedString() << normText << endl;


  // -------------------------------------------------------------------------
  // declare standard connection matrix and answer clusters for initialization
  // -------------------------------------------------------------------------
  TMatrixInt   CM(N,N,0);    // resize CM by command line parameters
  TClusterList answer(QMax);
  // -------------------------------------------------------------------------
  // declare Cp ('Cij-prime') and other program parameters
  // -------------------------------------------------------------------------
  TMatrixInt    Cp(N,N);                      // weighted connection matrix
  //TClusterList clusters(QMax,NEdges);      // define the working cluster list
  //TClusterList bestClusters(QMax,NEdges);  // current best of working clusters
  TVectorInt    trialEnergies(NTrialsMax);    // store energies of the trials
  TVectorFloat  trialModularities(NTrialsMax);// modularities of the trials
  TVectorInt    k(N);                         // store spin degrees (# of conn)
  TMatrix       pij(N,N);                     // probability matrix
  int equalEnergyCount = 0; // number of runs where energy is equal or better
  TStats1D      accStats(NRunsMax);
  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  cout << redText << "here a " << normText << flush;       // debugging

  int nRuns = 0;
  while(nRuns<NRunsMax) {
  nRuns++;
  int          bestTrial = 1;                 // which trial is best? (E or Q)
  trialEnergies.init(2 << 29); // store energies of the trials
  trialModularities.init(-1.0);   // modularities of the trials
  if(nRuns>1) {
    for(int i=0; i<QMax; i++)  answer[i].erase();  // clear answer cluster list
    answer.nClusters = 0;                          // clear answer cluster list
  } // end if nRuns

  cout << redText << "here b " << normText << flush;       // debugging

  // initialize connection matrix
  // read graph from data file or generate from scratch if no file specified
  int    errorCode;
  bool   isBinary;
  string binOutputFile;
  if(useInputFile && NRunsMax==1) {  // only import external file once
    if(inputFileName[inputFileName.size()-2]=='.'  && 
       inputFileName[inputFileName.size()-1]=='b'  )  isBinary = 1;
    else                                              isBinary = 0;
    //errorCode = readCMatrix_dimacs(CM, NEdges, inputFileName, isBinary);
    errorCode = readCMatrix_gml(CM, NEdges, answer, inputFileName,spinOffset);
    //cout << "  Number of edges was set to be " << NEdges << "\n";// debugging
    if(errorCode<=0)  errorMsg("There was a problem reading the input file.");
  } // end if useInputFile
  else {
    // generate CMatrix from scratch since no input file was specified
    initCMatrix(CM, answer, NEdges, QMax, pin, pout, 1); // even communities
  } // end else useInputFile
  if(showCMatrix)  coutCMatrix(CM);

  // initialize degrees vector - could be made twice as efficient
  for(int i=0; i<N; i++) {
    int degreeCount = 0;
    for(int j=0; j<N; j++)  if(CM[i][j]>0)  degreeCount++;
    k[i] = degreeCount;
  } // end for i
  // now init pij
  p = (double)NEdges/(double)(N*(N-1)/2); // average connection probability
  cout << "  Graph has " << N << " vertices and " << NEdges << " edges "
       << "with an average density of " << p << "\n";
  if(useERGraph) {
    cout << "  Using Erdos-Renyi graph for modularity calc (constant p)\n";
    pij.init(p);
    for(int i=0; i<N; i++)  pij[i][i] = 0.0;  // now delete diagonal
  } // end if useERGraph
  // otherwise use degree model for the probability
  else {
    cout << "  Using pij = ki*kj/2L (standard) modularity calc\n";
    double  pijVal, pijSum = 0.0;
    for(int i=0; i<N; i++) {
      pij[i][i] = 0.0;  // delete diagonal - should not be called
      for(int j=i+1; j<N; j++) {
        // uses symmetric CM and thus results in symmetric pij
        pijVal = (double)k[i]*k[j];
        pij[i][j] = pijVal;
        pij[j][i] = pijVal;
        pijSum   += pijVal;  // double checking
      } // end for j
    } // end for i
    pij *= ( (double)NEdges/pijSum );  // normalize so that sum(i!=j)=2NEdges
    //pij *= ( 1.0/(1.0*(double)NEdges) );  // normalize so that sum(i!=j)=1.0

    //for(int i=0; i<N; i++)  for(int j=0; j<N; j++)  pij[i][j] *= ( (double)NEdges/pijSum );
    if(debugMode>0) {    // retest sum
      cout << "  Degrees vector is:  " << k << " and average is " 
           << k.getAvg() << "\n";
      pijSum = 0.0;
      for(int i=0; i<N; i++)  for(int j=0; j<N; j++)  pijSum += pij[i][j];
      cout << "  Sum of pij matrix is " << pijSum << "\n";
    } // end if debugMode
  } // end else


  if(exitwNoRun)  // user execution stop specified by -x parameter
    abortMsg("User directed program to stop prior to main calculation.");

  // -------------------------------------------------------------------------
  // declare other program parameters depending on CM initialization (NEdges)
  // -------------------------------------------------------------------------
  TClusterList clusters(QMax,NEdges);      // define the working cluster list
  TClusterList bestClusters(QMax,NEdges);  // current best of working clusters

/*
  cout << redText << "here 3 " << normText << flush;       // debugging
  string fname_dm("cmp0"  + itos((int)((pij+0.005)*100.0)) 
                          + "n" + itos(N) + ".dat");
  string fname_dmb("cmp0" + itos((int)((p+0.005)*100.0)) 
                          + "n" + itos(N) + ".b");
  if(outputCMatrix || onlyCMatrix) {
    if(outputFileName[outputFileName.size()-2]=='.'  && 
       outputFileName[outputFileName.size()-1]=='b'  )  isBinary = 1;
    else                                                isBinary = 0;
    // Select which file name to use
    if(outputFileName=="") {
      if(isBinary)  outputFileName = fname_dmb;
      else          outputFileName = fname_dm;
    } // end if outputFileName

    // output information for each nTrial iteration
    if(outputCMatrix || onlyCMatrix) {
      if(isBinary)   cout << "Writing CMatrix in binary format:  ";
      else           cout << "Writing CMatrix in ascii format:  ";
      outputFileName = outputFileRoot + "_" + itos(nTrials+1) + ".dat";
      binOutputFile  = outputFileName + ".b";
      outputCMatrix_dimacs(CM, outputFileName, isBinary);
      if(!isBinary) {
        string extCommand("./asc2bin " + outputFileName + " " + binOutputFile);
        cout << "  Executing external command: \"" << extCommand << "\"\n";
        system(extCommand.c_str());
      } // end if !isBinary
    } // end if outputCMatrix

  } // end if outputCMatrix
  if(onlyCMatrix)  abortMsg("User directed to only output CMatrix.");
*/

  // -------------------------------------------------------------------------
  // now fill Cp with weights (usually +1 or -1) but set diagonal values to 0
  // -------------------------------------------------------------------------
  cout << redText << "here c " << weightC << normText << endl;    // debugging
  int cijTemp;
  for(int i=0; i<N; i++) { for(int j=i+1; j<N; j++) {
    cijTemp = CM[i][j];  
    //if(debugMode>1)  cout << cijTemp << " for i = " << i << " and j = " << j << endl; // debugging
    if(cijTemp>0)       { Cp[j][i] = cijTemp*weightC;  Cp[i][j] = cijTemp*weightC; }
    else if(cijTemp==0) { Cp[j][i] = weightU;          Cp[i][j] = weightU; }
    else errorMsg("There was a problem with the input connection matrix!");
  }} // end for ij
  for(int i=0; i<N; i++)  Cp[i][i] = 0;  // zero the diagonal values
  if(showCMatrix)  coutCMatrix(Cp);                       // debugging
  cout << redText << "here d " << weightU << normText << endl;    // debugging

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------
  answer.calcEnergy(pij,Cp);       // pre-calculate the answer energy
  answer.display(spinOffset,1,"Original answer ");
  // -------------------------------------------------------------------------
  // pre-test the constructed answer with the community detection algorithm to 
  // check for an optimal initial system
  //communityDetect(answer,Cp,CM,N,pij,NStepMax,CheckInterval,useEnergy,3);
  // -------------------------------------------------------------------------
  answer.sortClusters();  
  answer.calcEnergy(pij,Cp);      // pre-calculate the answer energy
  answer.display(spinOffset,1,"Tested answer ");
 
  // -------------------------------------------------------------------------
  // begin algorithm
  // -------------------------------------------------------------------------
  int nTrials = 0;
  while(nTrials<NTrialsMax) {
  nTrials++;  // start trial loop with nTrials = 1
  
  //cout << redText << "  here 4\n" << normText << flush;  // debugging
/*
  // output information for each nTrial iteration
  if(outputCMatrix || onlyCMatrix) {
    if(isBinary)   cout << "Writing CMatrix in binary format:  ";
    else           cout << "Writing CMatrix in ascii format:  ";
    outputFileName = outputFileRoot + "_" + itos(nTrials+1) + ".dat";
    binOutputFile  = outputFileName + ".b";
    outputCMatrix_dimacs(CM, outputFileName, isBinary);
    if(!isBinary) {
      string extCommand("./asc2bin " + outputFileName + " " + binOutputFile);
      cout << "  Executing external command: \"" << extCommand << "\"\n";
      system(extCommand.c_str());
    } // end if !isBinary
  } // end if outputCMatrix
*/
  if(verbosity>0) {
    cout << greenText << "Starting iteration";
    if(NTrialsMax>1)  cout << " for trial " << nTrials;
    cout << ":\n";
  } // end if verbosity
  
  // -------------------------------------------------------------------------
  // initialize clusters
  // -------------------------------------------------------------------------
  if(nRuns>1) {
    for(int i=0; i<QMax; i++)  clusters[i].erase();  // clear clusters list
    clusters.nClusters = 0;                          // clear clusters list
  } // end if nRuns
  initRandomClusters(clusters,N,QMax); // initialize a random configuration
  clusters.calcEnergy(pij,Cp);// sum starting energy for all clusters
  clusters.initClean();                // initialize all to clean (unmoved)
  if(debugMode>0) {
    clusters.display(spinOffset);      // display random initial state
    cout << magentaText << "Initial energy = "   << clusters.getEnergy()
         << " and modularity = " << clusters.getModularity() << normText
         << "\n";
  } // end if debugMode

  //cout << "here a before community detection" << endl;  // debugging

  // -------------------------------------------------------------------------
  // run the actual community detection algorithm for the unknown system
  communityDetect(clusters,Cp,CM,N,pij,NStepMax,CheckInterval,useEnergy,3,debugMode);
  clusters.calcEnergy(pij,Cp);  // perform final manual energy calculation
                                //   should be redundant - debugging
  // -------------------------------------------------------------------------

  //cout << "here b after community detection with nStep = " << nStep << endl;  // debugging


  // keep track of best for several trial cases
  //cout << magentaText << clusters.getEnergy() << bestClusters.getEnergy() << endl;
  if(useEnergy) {
    if(nTrials==1) { 
      //cout << "here c " << bestClusters.nClusters << " " << clusters.nClusters << endl;  // debugging
      bestClusters = clusters;  
      //cout << "here d " << endl;  // debugging
      bestClusters.calcEnergy(pij,Cp); 
      //cout << "here e " << endl;  // debugging
    }
    else if(clusters.getEnergy() < bestClusters.getEnergy()) {
           //cout << redText << clusters.getEnergy() << bestClusters.getEnergy() << endl;
 
           bestClusters = clusters;  bestClusters.calcEnergy(pij,Cp);  bestTrial = nTrials;  }
  } // end if
  else {  // use modularity instead
    if(nTrials==1) { bestClusters = clusters;  bestClusters.calcEnergy(pij,Cp); }
    else if(clusters.getModularity() > bestClusters.getModularity()) {
           bestClusters = clusters;  bestTrial = nTrials;  }
  } // end else
  //cout << "here d after community detection" << endl;  // debugging

  trialEnergies[nTrials-1]     = clusters.getEnergy();       // indexed from 1
  trialModularities[nTrials-1] = clusters.getModularity();   // indexed from 1

  //cout << "here e after community detection" << endl;  // debugging

  if(debugMode>0) { // output trial by trial results
    // end program - output run information
    cout << "-------------------------------------------------------------------\n";  // debugging
    //cout << "Here A... " << endl;  // debugging
    clusters.display(spinOffset,1,"Trial " + itos(nTrials) + " ");
    //cout << "Here B... " << endl;  // debugging
    cout << "compared to constructed community energy of " << answer.getEnergy()
         << " and modularity of " << answer.getModularity() << "\n";
    cout << brownText << "There were " << clusters.getSize() << " final clusters.\n"
                      << "Trial " << nTrials << " done on step " << nStep << "." << normText << endl;
    cout << "-------------------------------------------------------------------\n";  // debugging
  } // end if debugMode for trial by trial results output
/*
  // try testing clustering "equality" for each trial
  clusters.sortClusters();
  int nMovedTotal = 0, nMoved;
  for(int i=0; i<clusters.getSize(); i++) {
    equalsTest(answer[i],clusters[i],nMoved);
    nMovedTotal += nMoved;
    // calculate percentage average of spins correct
    accSum += (double)nMoved;  // use as a running average
  } // end for i
  nMovedTotal /= 2;  // above double counts moved spins
*/
  } // end for nTrials
  // -------------------------------------------------------------------------

  // try testing clustering "equality" for the "best" clustering
  bestClusters.sortClusters();
  int     nMovedTotal = 0, nMoved;
  double  movedSum = 0.0, movedAvg;
/*
  for(int i=0; i<answer.getSize(); i++) {
    equalsTest(answer[i],bestClusters[i],nMoved);
    nMovedTotal += nMoved;
    // calculate percentage average of spins correct
    movedSum += (double)nMoved;  // use as a running average
  } // end for i
  nMovedTotal /= 2;  // above double counts moved spins
*/
  if(useEqualsTest)  equalsTest3(answer,bestClusters,N,nMoved);
  movedSum    = (double)nMoved;     // legacy assignments for below
  nMovedTotal = nMoved;             

  // now calculate accurace rate for standard tests
  movedAvg  = 1.0 - movedSum/(double)N;
  accSum += movedAvg;    acc2Sum += movedAvg*movedAvg;
  accStats.store(movedAvg);

  // end program - output run information
  //bestClusters.calcEnergy(pij,Cp);
  int     bestEnergy   = bestClusters.getEnergy();
  int     answerEnergy = answer.getEnergy();
  double  bestMod      = bestClusters.getModularity();
  double  answerMod    = answer.getModularity();
  cout << "Pre-defined communities are:" << normText << "\n";  // debugging
  //answer.calcEnergy(pij,Cp);  // calculated at beginning of trial
  cout << redText << "  here 4\n" << normText << flush;  // debugging
  answer.display(spinOffset,0,"Answer ");

  cout << "Here A... " << endl;  // debugging
  bestClusters.calcEnergy(pij,Cp); // redundant supposedly // debugging
  bestClusters.display(spinOffset,0,"Best trial ");
  cout << brownText;
  if(!useInputFile) {  // no answer set if an input file is used
    cout << redText   << "The accuracy rate for this run was " << movedAvg 
         << brownText << "\n";
    if(nMovedTotal>0)  cout << nMovedTotal << " spins changed position\n";
    else               cout << "Solution is exact compared to the answer.\n";
  } // end if

  cout << "Here B... " << endl;  // debugging
  //cout << magentaText << "Final energy was " << bestClusters.getEnergy() << "\n";
  if(bestEnergy>answerEnergy)  cout << redText;     // all trials failed
  else                         cout << magentaText; // ok
  if(bestEnergy<=answerEnergy) equalEnergyCount++;     // all trials failed
  cout << "Final energy was " << bestEnergy 
       << " compared to constructed community energy of " << answerEnergy;
  // flag color of trial text to indicate which solutions were good
  if(trialEnergies[0] <= answerEnergy)  cout << greenText; // 1st ok
  else if(bestEnergy>answerEnergy)      cout << redText;   // failed
  else                                  cout << greyText;  // both ok
  // output all energies and modularities of the trials
  if(NTrialsMax>1) { 
    cout << " (from trial " << bestTrial << greyText << ")\nEnergies for the " 
         << NTrialsMax << " trials were:      ";
    TVectorOut(cout,trialEnergies,1,1,greyText);  cout << "\n";
    cout << "Modularities for the " << NTrialsMax << " trials were:  "; 
    TVectorOut(cout,trialModularities,1,1,greyText);  cout << magentaText;
  } // end if NTrialsMax
  // output best modularity and ending information
  cout << "\nFinal modularity was " << bestClusters.getModularity() 
       << " compared to constructed modularity of " << answer.getModularity() << "\n";
  cout << brownText << "There were " << bestClusters.getSize() 
       << " final clusters.\nDone on step " << nStep << "." << normText << endl;
  if(outputSpinConnections>=0)  coutSpinConnections(CM,outputSpinConnections);


  } // end for nRuns

  double  accAvg   = accSum/(double)NRunsMax;
  double  accSigma = sqrt(acc2Sum/(double)NRunsMax - accAvg*accAvg);
  cout << redText << "The accuracy average was " << accAvg
       << " on a total of " << NRunsMax << " runs with a standard deviation of "
       << accSigma << normText << "\n";
  cout << magentaText << "The number of runs that has an equal or lower energy" 
       << " was " << equalEnergyCount << normText << "\n";
  accStats.output(3);

  // ------ Calculate and output total run time info -------------------------
  outputRuntime(StartTime, EndTime, CurrTime);

  // return memory
  //delete Today;  delete TodayEnd;
  //cout << brownText << " Here end 1\n" << normText << flush;
  return 1;
} // end MAIN
