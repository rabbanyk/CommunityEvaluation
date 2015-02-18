/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Copyright © 2009  Peter Ronhovde                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 *  This program is free software: you can redistribute it and/or modify       *
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
 *  Created by Peter Ronhovde on 10/15/09 (email: ronhovde[at]hbar.wustl.edu)  *
 *	Modified on 11/12/09                                                       *
 *  Location: Washington University in St. Louis, St. Louis, MO 63101          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

// standard includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
// Note:  Use the -DDEFINE_NAME for the g++ compile.  This is particularly 
// necessary for the weighted implementation (MSL_TCMWEIGHTED) of the source 
// code to be implemented across all files.

//# define DOS_OS
#define LINUX_OS
// local includes with verbose pre-processor guards
unsigned long long ISeed; // random number seed
#ifndef ML_UTILS_H
#include "./msl/ML_Utils.h"
#endif
#ifndef MSL_TSTATS1D_H
#include "./msl/MSL_Stats1D.h"
#endif
#ifndef MSL_TARRAY1D_H
#include "./msl/MSL_Array1D_Template.h"
#endif
#ifndef MSL_TARRAY2D_H
#include "./msl/MSL_Array2D_Template.h"
#endif
#ifndef MSL_TCLUSTERLIST_H
#include "clusterclasses.h"
#endif
#include <sys/times.h>
#ifndef MSL_TVECTORBAG_T
#include "./msl/MSL_VectorBag_Template.h"
#endif

#ifndef MSL_TCMSPARSE_UNWEIGHTED_H
// the edges only data structure
//#include "./msl/MSL_CMatrix_Sparse_Unweighted.h"
#endif
#ifndef MSL_TCMSPARSEW_H
// the edges and data data structure
#include "./msl/MSL_CMatrix_SparseW.h"
#endif
#ifndef MSL_TCMDENSE_H
// the complete matrix data structure
#include "./msl/MSL_CMatrix_Dense.h"
#endif


#include "cluster.hpp"
#include "clustermerges.cpp"
#include "clustermeasures.hpp"
#include "clusterstats.hpp"
#include "RNMain.hpp"
#include "RNMRA.hpp"
#include "noisetest.hpp"

//#include "mytimefns.cpp"
time_t theTime;
tms    StartTime;//  times(&StartTime);  // gcc specific time variables
tms    CurrTime, EndTime;                // gcc specific time variables

time_t localTime;
tms    localStartTime;//  times(&StartTime);  // gcc specific time variables
tms    localCurrTime, localEndTime;           // gcc specific time variables


using namespace std;
using namespace MSL;

// *** Begin MAIN ************************************************************
int main(int argc, char *argv[]) {
  cout << normText << boldText << "\ncluster.cpp code output\n";
  echoCMD(argc,argv);  // echo the command line (mostly for file record keeping)

  // interpret command line parameters
  #include "cluster_clparams.cpp"

  // -------------------------------------------------------------------------
  // make variable and parameter changes specified by user
  // -------------------------------------------------------------------------
  // set random seed either by user value or by class initialization
  initRNG(ISeed);  // if ISeed>0 we assume the user specified a seed

  // -------------------------------------------------------------------------
  // begin algorithm
  // -------------------------------------------------------------------------

  #ifdef LINUX_OS
  // track the initialization time of either solver
  double initTime, calcTime, timeSum;
  initTime = getRuntime(StartTime,EndTime,CurrTime);
  initTimes(StartTime, theTime);
  #endif

// ---------------------------------------------------------------------------
#if ( defined(MRA_SOLVER) || defined(APM_GREEDY_SOLVER) )
  // In either solver we need to input a connection matrix definition, and we 
  // may want to input an answer file if we are testing one of the solvers.
  // General comments:
  // (a) The solvers are implemented with a weighted sparse matrix using 
  //     weighted edges if the parameter MSL_TCMWEIGHTED is defined at compile 
  //     time.  However, most of the speed is recovered for unweighted systems 
  //     if the parameter is not defined.  If MSL_TCMWEIGHTED is not defined, 
  //     the connection matrix is treated as if it is *unweighted* even if it 
  //     contains weighted edges.
  // (b) The sparse data structure does not allow weighted *missing* links (all
  //     are assumed to have a constant weight of -1).  The dense array 
  //     structure will allow them, but the input routines do not currently 
  //     recognize them, so they will have to be manually set by 
  //     Cp.data[iNode][jNode] = value.
  // (c) Even though the APM naturally handles directed networks, none of the 
  //     solvers currently implement directed networks due to legacy code 
  //     implementations.

  //TCMDense  Cp(N);  // the dense matrix implementation
  TCMSparseW  Cp(N);  // the sparse matrix implementation
  if(useInputFile) {  // import external file
    int errorCode;
    // finish general descriptive information
    // multiple quotes are format for embedded quotes around "filename"
    Cp.description += "\"Input file name = \"\"" + inputFileName + "\"\"\"\n";
    cout << "Importing file in hier.cpp... " << flush;  // debugging
    errorCode = Cp.input(inputFileName,nodeOffset,useSymmetrizedCM);
    //if(errorCode<=0) errorMsg("There was a problem reading the input file ("
    //                          +itos(errorCode)+").");
  } // end if useInputFile
  else 
    errorMsg("User must specify a valid gml or ascii text array data file.");
  if(showCMatrix)  Cp.display();  // should only use for very small systems

  // parameters to control the maximum system size will display to the console
  int  maxDisplayq = 201, maxDisplayN = 10001;

  // input answer cluster list if directed by the user 
  // [-uansf:1 (default) parameter along with -infans:filename (required)]
  TClusterList answer;
  if(useAnswer) {
    if(inputAnswerFile=="")
      errorMsg("User specified to use an answer file, but there is no answer file name?");

    if(verbosity>0)  cout << "Reading answer cluster list... " << endl;
    answer.createSymmetric(N,1); // with random ordered nodes
    // the raw answer file is an undocumented text file with each "line" being
    // a cluster.  Each node in the cluster represented by a single integer
    // (indexed from 0 to N-1 generally) separated by spaces.
    if(verbosity>0)  cout << "Input answer list... " << endl;
    answer.inputRaw(inputAnswerFile,N);
    answer.calcEnergy(Cp);
    if(verbosity>1 && answer.getq()<maxDisplayq && 
       answer.getNNodes()<maxDisplayN)  
      answer.display(0,1,"answer ");  cout << "\n";
  } // end if useAnswer

  // define a "best" cluster list which is used for any solver
  TClusterList  best;
  best.createSymmetric(N,1); // creates the starting best cluster list
  double gammaBest;

  // What type of initialization do we use in either solver?
  // These initialization types are defined in cluster_clparams.cpp file.
  // N is the number of nodes in the system for all initialization types.
  // Command line parameters are:  N (-n:#), q (-q:#), 
  //     nMin (-nmin:#), nMin (-nmin:#), betaLFR (-betaLFR:negative#)
  // (1) initSymmetric(N,1);            // symmetric initialization with q=N
  // (2) initRandom(N,q,1);             // random assignments with q communities
  // (3) initEven(N,q,1);               // even sizes with q communities
  // (4) initPower(N,nMin,nMax,betaLFR,1);// power law init with beta negative
  pInitMethod = &initSymmetric;         // Symmetric is the implemented method.

  // user execution stop specified by -x parameter (useful to quickly check
  // passed command line parameters)
  if(exitwNoRun)
    abortMsg("User directed program to stop prior to main calculation.");
#else
  errorMsg("Please define either MRA_SOLVER or APM_GREEDY_SOLVER preprocessor definition");
#endif
// ---------------------------------------------------------------------------
  
  // track the solution time of either solver
  #ifdef LINUX_OS
  calcTime = getRuntime(StartTime,EndTime,CurrTime);
  timeSum  = initTime + calcTime;
  initTimes(StartTime, theTime);  // temporarily init times here
  #endif
  
// ---------------------------------------------------------------------------
#ifdef MRA_SOLVER
  // The general version of the MRA algorithm.  
  // (1) It returns by reference the "best" answer based on the highest NMI and 
  //     secondarily the most stable plateau in H.
  // (2) Output data files are the "best" community structure and a multi-
  //     resolution data file showing the behavior of the replica correlations
  //     over all tested resolutions.
  if(verbosity>0)
    cout << brown << "OK.  Starting MRA solver..." << normText << endl;

  // now apply the multiresolution algorithm using a log step size (-gsteps:#)
  // with gammaNSteps per decade of gamma.
  // The connection matrix Cp should not be scaled in any way (other than 
  // weighted edges) prior to the rnMRA function call.
  if(verbosity>0)  cout << "Solving system with MRA algorithm... " << endl;
  
  gammaBest = rnMRA(Cp,outputFileName,best,*pInitMethod,
        gammaStart,gammaStop,gammaNSteps,useEnergy,useMerges,useZeroE,
        NReplicasMax,NTrialsMax,NStepsMax,NIterMax,nodeOffset,verbosity);
#endif
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
#ifdef APM_GREEDY_SOLVER
  // The general version of a one use solver using the "Absolute Potts Model"
  // (APM) given a single value of gamma.
  // (1) It returns the "best" answer based on the lowest energy trial.
  // (2) Output data file is the same "best" community structure.
  if(verbosity>0)  
    cout << brown << "OK.  Starting greedy APM solver..." << normText << endl;

  #ifdef MSL_TCMWEIGHTED
  // get the integer a and b from the user gammaAPM (via the -g:# parameter)
  // since the Cp connection matrix is integer valued.
  int  a, b, CpIntegerScale = 100;
  gammaToab(gammaAPM,a,b,CpIntegerScale);
  if(verbosity>0)
    cout << lightred << "\ngammaAPM = " << gammaAPM << " with a = " << a 
         << " and b = " << b << normText << endl;
  // Note that this assumes that the original connection matrix is unscaled!
  // We only scale Cp for this solver (not the MRA solver above) since 
  // communityBestE performs no operations on Cp itself.  We do not undo the 
  // scale after the solution since Cp is no longer used in these solvers.
  Cp.scale(a,b,1,-1);  // assumes an unscaled connection matrix!
  #endif

  // Apply the greedy community detection algorithm using the APM
  double nIterAvg, nMergesAvg, solveTime;
  if(verbosity>0)  cout << "Solving with communityBestE()... " << endl;
  communityBestE(Cp,best,nIterAvg,nMergesAvg,solveTime,gammaAPM,*pInitMethod,
                 NTrialsMax,NIterMax,useMerges,useZeroE,0,verbosity);
  gammaBest = gammaAPM;  // for compatibility with general code below.
#endif
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Ending code for either solver
#if ( defined(MRA_SOLVER) || defined(APM_GREEDY_SOLVER) )
  // what do we do with the "best" answer with either solver?
  if(verbosity>0)  cout << "Calculating solution energy... " << endl;
  best.calcEnergy(Cp);
  // Use best value of gamma for the RB energy.  Note that we scale the RB 
  // energy by L (or W), the number of edges, so that there is a seemless 
  // equivalence to modularity at gamma = 1.  Also note, there is no equivalence
  // between the APM and the RB gammas.  Gamma is used strictly to have a 
  // corresponding value of RB energy.
  best.calcQ(Cp,gammaBest,0.0);

  // What do we do with the best solution?
  if(verbosity>1 && best.getq()<maxDisplayq && best.getNNodes()<maxDisplayN)  
    best.display(0,1,"best after solve ");
  if(outputFileName!="") {
    string filename = outputFileName + "Best.txt";
    best.outputRaw(filename);
  } else //if(verbosity>0)  
    warningMsg("No output file name for the best solution was specified (-outf:filename)");

  // allow a comparison to a known 'answer' if directed by the user
  if(useAnswer) {
    if(verbosity>0)  
      cout << brown << "Calculating information comparisons compared to known "
           << "\'answer\'... " << endl;
    double tNMI, tVI, tMI;
    int nMoves = -1;  // initialize to invalid value
    calcNMVILL(answer,best,N,tNMI,tMI,tVI,nMoves);
    if(verbosity>0)  
      cout << magenta << "NMI = " << tNMI << ", VI = " << tVI << ", MI = " 
           << tMI << " with " << nMoves << " nodes misplaced)\n";
  } // end i
#endif 
// ---------------------------------------------------------------------------

  #ifdef LINUX_OS
  outputTime(initTime,"initialization");  cout << "\n";
  outputTime(calcTime,"calculation");  cout << "\n";
  outputTime(timeSum,"run");  cout << "\n";
  double FinalizationTime = getRuntime(StartTime,EndTime,CurrTime);
  outputTime(FinalizationTime,"finalization");  cout << endl;
  #endif

  cout << brown << "Done with program." << normText << endl;
  return 1;
} // end MAIN
