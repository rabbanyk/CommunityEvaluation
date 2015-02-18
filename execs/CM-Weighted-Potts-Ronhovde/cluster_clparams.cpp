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
int     NStepsMax = 1, NTrialsMax = 4,   NGraphsMax = 10, 
        NRunsMax = 1,  NReplicasMax = 8, NIterMax = 1000;
        // step checking parameters - for ending trial runs
int     NStepStartCheck = 10, CheckInterval = 10, nStep;
int     wAnswer = 2;

int     N = 64,  NEdges = 0;       // starting number of nodes and edges
double  p = 0.2, pin = 0.75, pout = 0.05; // connection probability
int     weightC = 1, weightU = -1;// weights of connected and unconnected edges
int     weightStep = 1;           // weight increment for hierarchy loop
double  epsilon = 0.5e-10;        // error tolerance on float equals tests
int     QMax = 50;                // maximum possible number of communities (legacy)
int     q = 6;                    // the number of communities


//extern  unsigned long long ISeed = 0;     // random number seed
        ISeed = 0;     // random number seed
int     debugMode = 0;            // debug level:  0 is none
int     verbosity = 0;            // how much non-debugging info is reported - 
                                  //   is not implemented in all functions
string  sText = "";
double  temperature = 0.0;

double  rScale = 1.0;
double  vScale = 1.0;
double  vShift = 0.0;

double  gammaAPM    = 1.0;
double  gammaStart  = 19.0;        // starting weight for rnMRA solver
double  gammaStop   = 0.001;       // ending weight for rnMRA solver
int     gammaNSteps = 20;          // number of gamma steps per decade for rnMRA

bool   useFixedValency = 0;       // use an input data file?
//bool useRestart      = 0;       // use a restart file? - invalid right now
bool   useInputFile    = 0;       // do we use an input data file?
bool   exitwNoRun      = 0;       // exit just prior to main calculations?
bool   outputCMatrix   = 0;       // output CMatrix to file?
bool   onlyCMatrix     = 0;       // write CMatrix file and exit?
bool   showCMatrix     = 0;       // display version of CMatrix to console?
bool   useEnergy       = 1;       // use energy for convergence criteria? 
                                  //   (only case at the moment)
bool   useDensityStep  = 1;       // does hierarchy function step by density
                                  //   or by parameter weights?

bool   useMerges       = 1;       // use merges in some dynamic solutions?
bool   useNodeRefine   = 1;       // use node refinement after any merges?
bool   useAnswer       = 0;       // do we utilize an "answer" cluster list?

double densityStart    = 1.0;     // at what value do we start the iteration?
double densityScale    = 0.75;    // at what value do we scale the iteration?
double densityStep     = 0.025;   // at what value do we step the iteration?

double kStart          = 1.0;     // at what degree do we start the iteration?
double kStop           = 0.75;    // at what degree do we scale the iteration?
double kStep           = 0.025;   // at what degree do we step the iteration

bool   useDense        = 0;       // use dense method for cluster searching?
                                  //   - deprecated, only used in old search
bool   useERGraph      = 0;       // use energy for convergence criteria?
                                  //   (otherwise use modularity for testing)
bool   useHierarchy    = 0;       // use hierarchical initialization
bool   useEqualsTest   = 0;       // use equals test for known systems
bool   useConstrainedH = 0;       // use constrained hierarchical method
bool   useConstrainedq = 0;       // use constrained number of communities q?
bool   usePowerInit    = 0;       // do we use a power initialization of list?
bool   useEvenInit     = 0;       // do we use an even initialization of list?
bool   useRandomInit   = 0;       // do we use a random initialization of list?
bool   useSymmetricInit = 0;      // do we use a Symmetric initialization?
bool   useManualDisplay = 0;      // manually display certain cluster solutions?
bool   useExpNoise      = 0;      // do we use exponential or uniform noise?
bool   usePowerNoise    = 0;      // do we use power-law or uniform noise?
bool   useZeroE         = 0;      // do we allow zero energy moves?
bool   useSymmetrizedCM = 0;      // do we symmetrize a directed matrix

string inputFileName(""), outputFileName(""), outputFileRoot(""),
       inputAnswerFile("");

double  zout = 4.0;               // for fixed valency per node, what is the 
                                  // avg number of nodes outside the community?
double  zin  = 12.0;              // for fixed valency per node, what is the 
                                  // avg number of nodes outside the community?
double  zavg = 16.0;              // for fixed valency per node, what is the 
                                  // avg total links per node
double  zMin = 8.0;              // maximum possible degree per node (LFR test)
double  zMax = 50.0;              // maximum possible degree per node (LFR test)
int     nMax = -1;                // maximum size of communities (LFR test)
int     nMin = -1;                // minimum size of communities (LFR test)
int     IterType = 0;             // "iteration type" - generic value
double  rAFG = 0.0;               // AFG node self-weight
double  poutStop  = 0.35;         // ending density for noise test
double  poutStep  = 0.025;        // density step for noise test
double  poutStart = poutStep;     // starting density for noise test
double  muStart = 0.0;            // starting density for noise test
double  muStop  = 0.7;            // ending density for noise test
double  muStep  = 0.05;           // density step for noise test
double  gammaLFR = 2.0;           // exponent for degree distribution in LFR
double  betaLFR  = 1.0;           // exponent for cluster size distribution

int    nodeOffset = 0;            // offset for displaying the answer
int    outputNodeConnections = -1;// -1 flag value ignores at end, >0 is valid

  // Interpret command line parameters
  // Currently looks for command line parameter of form "a.out -xabc=#"
  // ex. command forms are:  clique.o -rseed=24357, clique.o /t:100.0
  char   cParam;
  string sParam, sFlag, sInputVal;
  int    eqPos, theCharSum;
  if(argc>1) {
    cout << greenText << "Interpreting command line parameters:" 
         << normText << "\n";
    for(int i=1; i<argc; i++) {
      //cout << "Testing with i = " << i << endl;  // debugging
      cParam = argv[i][0];
      sParam = argv[i];
      eqPos  = sParam.find_first_of("=",1); // defaults to last position
      // or maybe user used a ':' as the separation character
      if(eqPos==string::npos)  eqPos = sParam.find_first_of(":",1);
      // now find command line parameter string with above data
      if( ((cParam=='-')||(cParam=='/')) && (sParam.length()>1) ) {
        sFlag      = sParam.substr(1,eqPos-1);
        if(eqPos==string::npos)  sInputVal = "";  // no input value
        else sInputVal  = sParam.substr(eqPos+1,sParam.length());
        theCharSum = charSum(sFlag,1); // convert case to lower and sum chars
      } // end valid flag
      else theCharSum = 0; // invalid flag (null char)
      switch ( theCharSum ) {
        case 'a':
          cout << "  Changing which comparison answer for acccuracy check from " 
               << itos(wAnswer,1);
          wAnswer = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(wAnswer,1) << endl;
          break;
        case 'c':
          cout << "  Changing convergence interval check from " 
               << itos(CheckInterval,1);
          CheckInterval = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(CheckInterval,1) << endl;
          break;
        case ('c'+'s'+'t'): {
          cout << "  Changing convergence interval start check from " 
               << itos(NStepStartCheck,1);
          NStepStartCheck = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(NStepStartCheck,1) << endl;
          break;
        }
        case ('d'+'c'+'m'): {
          if(sInputVal=="")  showCMatrix = 1;
          else               showCMatrix = (bool)(argv[i][eqPos+1] - 48);
          cout << "  Display CMatrix set to " << showCMatrix << endl;
          break;
        }
        case ('d'+'e'+'s'+'c'+'r'): {
          //int k = 1;
          //while(argv[i][eqPos+k]!='\0') { sText += argv[i][eqPos+k];  k++; }
          //char *pc = &(argv[i][eqPos+1]);
          //while((*pc)!='\0') { sText += (*pc);  pc++; }
          sText = &(argv[i][eqPos+1]);
          cout << "  Description text set to " << sText << endl;
          break;
        }
        case ('d'+'m'): {
          cout << "  Changing debug mode level from " << itos(debugMode,1);
          debugMode = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(debugMode,1) << endl;
          break;
        }
        case ('d'+'i'): { 
          cout << "  Changing starting density from " << densityStart;
          densityStart = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << densityStart << endl;
          if(densityStart>1.0 || densityStart<0.0) { 
            //cerr << "Density value error!!!  "
            //     << "Invalid value is " << densityStart << ". Aborting.\n";
            errorMsg("Invalid start density!  Value is "+ftos(densityStart));
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('d'+'s'): { 
          cout << "  Changing density step from " << densityStep;
          densityStep = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << densityStep << endl;
          if(densityStep>1.0 || densityStep<0.0) { 
            //cerr << "Density value error!!!  "
            //     << "Invalid value is " << densityStart << ". Aborting.\n";
            errorMsg("Invalid density step.  Value is "+ftos(densityStep));
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('d'+'s'+'c'): { 
          cout << "  Changing density scale from " << densityScale;
          densityScale = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << densityScale << endl;
          if(densityScale<0.0 || densityScale>0.9999) { 
            errorMsg("Invalid density scale.  Value is "+ftos(densityScale));
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('r'+'s'): { 
          cout << "  Changing radius scale from " << rScale;
          rScale = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << rScale << endl;
          if(rScale<0.0) { 
            errorMsg("Invalid radius scale.  Value is "+ftos(rScale));
            return 0;
          } // end if radius scale error check
          break;
        }
        case ('v'+'s'+'c'): { 
          cout << "  Changing potential scale from " << vScale;
          vScale = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << vScale << endl;
          break;
        }
        case ('v'+'s'+'h'): { 
          cout << "  Changing potential shift from " << vShift;
          vShift = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << vShift << endl;
          break;
        }
        case ('e'+'r'+'i'+'n'+'i'+'t'): {
          if(sInputVal=="")  useERGraph = 1;
          else               useERGraph = (bool)(argv[i][eqPos+1] - 48);
          //cout << "Testing exit flag:  " << (argv[i][eqPos+1]) << endl;
          if(useERGraph) 
            cout << magenta << "  User specified to calculate modularity." 
                 << "based on an Erdos-renyi graph." << normText << endl;
          break;
        }
        case ('f'+'v'+'a'+'l'+'n'+'c'): { 
          if(sInputVal=="")  useFixedValency = 1;
          else               useFixedValency = (bool)(argv[i][eqPos+1] - 48);
          break;
        }
        case ('g'): { 
          cout << "  Changing degree exponent in LFR test from " << gammaAPM;
          gammaAPM = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << gammaAPM << endl;
          //if(gammaLFR<0.0) { 
          //  errorMsg("Invalid gamma exponent.  Value is "+ftos(gammaLFR));
          //  return 0;
          //} // end if
          break;
        }
        case ('g'+'s'+'t'+'a'+'r'+'t'): { 
          cout << "  Changing gamma Potts model weight from " << gammaStart;
          gammaStart = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << gammaStart << endl;
          if(gammaStart<0.0) { 
            errorMsg("Invalid gammaStart value "+ftos(gammaStart));
            return 0;
          } // end if
          break;
        }
        case ('g'+'s'+'t'+'o'+'p'): { 
          cout << "  Changing gamma Potts model weight from " << gammaStop;
          gammaStop = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << gammaStop << endl;
          if(gammaStop<0.0) { 
            errorMsg("Invalid gammaStop value "+ftos(gammaStop));
            return 0;
          } // end if
          break;
        }
        case ('g'+'s'+'t'+'e'+'p'+'s'): {
          cout << "  Changing number of gamma steps per decade from " 
               << itos(gammaNSteps,1);
          gammaNSteps = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(gammaNSteps,1) << endl;
          if(gammaNSteps<=0)  
            errorMsg("Invalid value for gammaNSteps = "+itos(gammaNSteps)+"?");
          break;
        }
        case ('g'+'l'+'f'+'r'): { 
          cout << "  Changing degree exponent in LFR test from " << gammaLFR;
          gammaLFR = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << gammaLFR << endl;
          //if(gammaLFR<0.0) { 
          //  errorMsg("Invalid gamma exponent.  Value is "+ftos(gammaLFR));
          //  return 0;
          //} // end if
          break;
        }
        case ('b'+'l'+'f'+'r'): { 
          cout << "  Changing cluster size exponent in LFR test from " << betaLFR;
          betaLFR = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << betaLFR << endl;
          //if(betaLFR<0.0) { 
          //  errorMsg("Invalid beta exponent.  Value is "+ftos(betaLFR));
          //  return 0;
          //} // end if
          break;
        }
        case ('i'+'n'+'f'): { 
          cout << "  Using CMatrix input data file \"";
          inputFileName = &argv[i][eqPos+1];
          cout << inputFileName << "\"\n";
          // now set the boolean variable to tell program to use the file
          useInputFile = 1;
          //cout << "Debugging:  Exiting input file assignment with argc = "
          //     << argc << endl;  // debugging 
          break;
        }
        case ('i'+'n'+'f'+'a'+'n'+'s'): { 
          cout << "  Using input \'answer\' data file \"";
          inputAnswerFile = &argv[i][eqPos+1];
          cout << inputAnswerFile << "\"\n";
          useAnswer = 1;  // turn on answer code
          break;
        }
        case ('i'+'t'+'y'+'p'+'e'): {
          cout << "  Changing iteration type from " 
               << itos(IterType,1);
          IterType = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(IterType,1) << endl;
          if(IterType<0 || IterType>3)  
            errorMsg("Invalid value of "+itos(IterType)+" for IterType?");
          break;
        }
        case ('n'): { 
          cout << "  Changing number of nodes from " << itos(N,1);
          N = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(N,1) << endl;
          break;
        }
/*
        case ('n'+'r'): {
          cout << "  Changing maximum number of runs from " 
               << itos(NRunsMax,1);
          NRunsMax = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(NRunsMax,1) << endl;
          break;
        }
*/
        case ('n'+'i'+'m'): {
          cout << "  Changing maximum community detection iterations from " 
               << itos(NIterMax,1);
          NIterMax = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(NIterMax,1) << endl;
          if(NIterMax<=0)  errorMsg("Invalid value for NIterMax?");
          break;
        }
        case ('n'+'r'+'m'): {
          cout << "  Changing maximum number of replicas from " 
               << itos(NReplicasMax,1);
          NReplicasMax = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(NReplicasMax,1) << endl;
          if(NReplicasMax<=0)  errorMsg("Invalid value for NReplicasMax?");
          break;
        }
        case ('n'+'g'+'m'): {
          cout << "  Changing maximum number of graphs from " 
               << itos(NGraphsMax,1);
          NGraphsMax = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(NGraphsMax,1) << endl;
          if(NGraphsMax<=0)  errorMsg("Invalid value for NGraphsMax?");
          break;
        }
        case ('n'+'s'+'m'): {
          cout << "  Changing maximum number of steps from " 
               << itos(NStepsMax,1);
          NStepsMax = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(NStepsMax,1) << endl;
          if(NStepsMax<=0)  errorMsg("Invalid value for NStepsMax?");
          break;
        }
        case ('n'+'t'+'m'): {
          cout << "  Changing maximum number of trials from " 
               << itos(NTrialsMax,1);
          NTrialsMax = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(NTrialsMax,1) << endl;
          if(NTrialsMax<=0)  errorMsg("Invalid value for NTrialsMax?");
          break;
        }
        case ('o'+'c'+'m'): { 
          cout << "  Specifying CMatrix output data file name of \"";
          outputFileRoot = &argv[i][eqPos+1];
          cout << outputFileRoot << "\"\n";
          // now set the boolean variable to tell program to use the file
          outputCMatrix = 1;
          break;
        }
        case ('o'+'c'+'m'+'x'): { 
          if(sInputVal=="")  onlyCMatrix = 1;
          else               onlyCMatrix = (bool)(argv[i][eqPos+1] - 48);
          cout << "  Only output CMatrix set to " << onlyCMatrix << endl;
          break;
        }
        case ('o'+'u'+'t'+'f'): { 
          cout << "  Using output data file \"" << green;
          outputFileName = &argv[i][eqPos+1];
          cout << outputFileName << normText << "\"\n";
          break;
        }
        case 'p': { 
          cout << "  Changing edge probability from " << p;
          p = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << p << endl;
          if(p>1.0 || p<0.0) { 
            cerr << "Edge connection probability error!!!  "
                 << "Invalid value is " << p << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('p'+'i'+'n'): { 
          cout << "  Changing edge probability inside a community from " << pin;
          pin = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << pin << endl;
          if(pin>1.0 || pin<0.0) { 
            cerr << "Edge connection probability error!!!  "
                 << "Invalid value is " << pin << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('p'+'o'+'u'+'t'): { 
          cout << "  Changing edge probability outside a community from " << pout;
          pout = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << pout << endl;
          if(pout>1.0 || pout<0.0) { 
            cerr << "Edge connection probability error!!!  "
                 << "Invalid value is " << pout << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('m'+'o'): { 
          cout << "  Changing mixing parameter start from " << muStart;
          muStart = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << muStart << endl;
          if(muStart<0.0) { 
            cerr << "Mixing parameter error!  "
                 << "muStart = " << muStart << "?  Aborting.\n";
            return 0;
          } // end if mixing parameter error check
          break;
        }
        case ('m'+'f'): { 
          cout << "  Changing mixing parameter stop from " << muStop;
          muStop = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << muStop << endl;
          if(muStop<0.0) { 
            cerr << "Mixing parameter error!  "
                 << "muStop = " << muStop << "?  Aborting.\n";
            return 0;
          } // end if mixing parameter error check
          break;
        }
        case ('m'+'s'+'t'+'e'+'p'): { 
          cout << "  Changing mixing parameter step size from " << muStep;
          muStep = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << muStep << endl;
          break;
        }
        case ('z'+'o'): { 
          cout << "  Changing degree start from " << kStart;
          kStart = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << kStart << endl;
          if(kStart<0.0) { 
            cerr << "k start parameter error!  "
                 << "kStart = " << kStart << "?  Aborting.\n";
            return 0;
          } // end if degree parameter error check
          break;
        }
        case ('z'+'f'): { 
          cout << "  Changing degree stop from " << kStop;
          kStop = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << kStop << endl;
          if(kStop<0.0) { 
            cerr << "k degree parameter error!  "
                 << "kStop = " << kStop << "?  Aborting.\n";
            return 0;
          } // end if degree parameter error check
          break;
        }
        case ('z'+'s'): { 
          cout << "  Changing degree step size from " << kStep;
          kStep = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << kStep << endl;
          break;
        }
        case ('p'+'o'): { 
          cout << "  Changing starting noise density from " << poutStart;
          poutStart = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << poutStart << endl;
          if(poutStart>1.0 || poutStart<0.0) { 
            cerr << "Edge density probability error!!!  "
                 << "poutStart is " << poutStart << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('p'+'f'): { 
          cout << "  Changing final noise density from " << poutStop;
          poutStop = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << poutStop << endl;
          if(poutStop>1.0 || poutStop<0.0) { 
            cerr << "Edge density probability error!!!  "
                 << "poutStop is " << poutStop << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('p'+'s'): { 
          cout << "  Changing noise density step size from " << poutStep;
          poutStep = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << poutStep << endl;
          if(poutStop>1.0 || poutStop<0.0) { 
            cerr << "Edge density probability error!!!  "
                 << "poutStep is " << poutStep << ". Aborting.\n";
            return 0;
          } // end if edgeProb error check
          break;
        }
        case ('q'): { 
          cout << "  Changing number of communities from " << itos(QMax,1);
          QMax = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(QMax,1) << endl;
          break;
        }
        case ('n'+'m'+'a'+'x'): { 
          cout << "  Changing maximum number of communities from " 
               << itos(nMax,1);
          nMax = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(nMax,1) << endl;
          break;
        }
        case ('n'+'m'+'i'+'n'): { 
          cout << "  Changing minimum number of communities from " 
               << itos(nMin,1);
          nMin = atoi(&argv[i][eqPos+1]);
          //NDim = NPar - 1;
          cout << " to " << itos(nMin,1) << endl;
          break;
        }
        case 'r': { 
          cout << "  Changing AFG self-weight from " << rAFG;
          rAFG = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << rAFG << endl;
          break;
        }
        case 't': { 
          cout << "  Changing temperature from " << temperature;
          temperature = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << temperature << endl;
          break;
        }
        case ('r'+'s'+'e'+'e'+'d'): { 
          cout << "  Changing random seed from " << itos(ISeed,1);
          //ISeed = atoi(&argv[i][eqPos+1]);
          ISeed = strtoul(&argv[i][eqPos+1],NULL,0);
          cout << " to " << ISeed << endl;
          break;
        }
        case ('s'+'o'+'f'+'f'+'s'+'e'+'t'): { 
          cout << "  Changing node display offset from " << nodeOffset;
          nodeOffset = atoi(&argv[i][eqPos+1]);
          cout << " to " << nodeOffset << endl;
          break;
        }
        case ('u'+'a'+'n'+'s'+'f'): { 
          if(sInputVal=="")  useAnswer = 1;
          else               useAnswer = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useAnswer)  
            cout << "  Using an answer cluster list (if available)\n";
          else 
            cout << "  Excluding an answer cluster list (if implemented)\n";
          cout << normText;
          break;
        }
        case ('u'+'h'+'c'): { 
          if(sInputVal=="")  useConstrainedH = 1;
          else               useConstrainedH = ((bool)(argv[i][eqPos+1] - 48));
          cout << brown << "  Using ";
          if(useConstrainedH)  cout << "constrained hierarchical method\n";
          else                 cout << "unconstrained hierarchical method\n";
          cout << normText;
          break;
        }
        case ('u'+'c'+'q'): { 
          if(sInputVal=="")  useConstrainedq = 1;
          else               useConstrainedq = ((bool)(argv[i][eqPos+1] - 48));
          if(useConstrainedq)  cout << brown << "  Using a constrained";
          else                 cout << brown << "  Using an unconstrained";
          cout << " number of communities q\n" << normText;
          break;
        }
        case ('u'+'d'): { 
          if(sInputVal=="")  useDensityStep = 1;
          else               useDensityStep = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useDensityStep) cout << "  Using density step in hierarchy iteration\n";
          else               cout << "  Using weight step in hierarchy iteration\n";
          cout << normText;
          break;
        }
        case ('u'+'m'+'e'+'r'): { 
          if(sInputVal=="")  useMerges = 1;
          else               useMerges = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useMerges)  cout << "  Using merges in certain dynamics\n";
          else           cout << "  Not using merges in certain dynamics\n";
          cout << normText;
          break;
        }
        case ('u'+'n'+'r'): { 
          if(sInputVal=="")  useNodeRefine = 1;
          else               useNodeRefine = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useNodeRefine)  cout << "  Using node refinements after any merges\n";
          else               cout << "  Not using node refinements after any merges\n";
          cout << normText;
          break;
        }
        case ('u'+'p'+'o'+'w'): { 
          if(sInputVal=="")  usePowerInit = 1;
          else               usePowerInit = ((bool)(argv[i][eqPos+1] - 48));
          usePowerNoise = usePowerInit;  // default to same for power noise
          cout << green;
          if(usePowerInit)  cout << "  Using power initialization\n";
          else              cout << "  Power initialization is turned off\n";
          if(usePowerNoise) cout << "  Using power noise (degree distribution)\n";
          cout << normText;
          break;
        }
        case ('u'+'e'+'v'+'e'+'n'): { 
          if(sInputVal=="")  useEvenInit = 1;
          else               useEvenInit = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useEvenInit)  cout << "  Using even initialization\n";
          else             cout << "  Using random n initialization\n";
          cout << normText;
          break;
        }
        case ('u'+'r'+'a'+'n'+'d'): { 
          if(sInputVal=="")  useRandomInit = 1;
          else               useRandomInit = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useRandomInit)  cout << "  Using random n initialization\n";
          else               cout << "  Using even initialization\n";
          cout << normText;
          break;
        }
        case ('u'+'s'+'y'+'m'+'m'): { 
          if(sInputVal=="")  useSymmetricInit = 1;
          else               useSymmetricInit = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useSymmetricInit)  cout << "  Using symmetric initialization\n";
          else                  cout << "  Turning off symmetric initialization\n";
          cout << normText;
          break;
        }
        case ('u'+'e'+'n'): { 
          if(sInputVal=="")  useEnergy = 1;
          else               useEnergy = ((bool)(argv[i][eqPos+1] - 48));
          cout << green << "  User specified to use ";
          if(useEnergy) cout << magenta << "energy";  
          else          cout << brown   << "modularity";
          cout << green << " as convergence test" << normText << "\n";
          break;
        }
        case ('u'+'e'+'q'): { 
          if(sInputVal=="")  useEqualsTest = 1;
          else  useEqualsTest = ((bool)(argv[i][eqPos+1] - 48));
          if(useEqualsTest)  cout << "  Turning on equals test\n";
          else               cout << "  Equals test is turned off\n";
          break;
        }
        case ('u'+'e'+'x'+'p'): { 
          if(sInputVal=="")  useExpNoise = 1;
          else               useExpNoise = ((bool)(argv[i][eqPos+1] - 48));
          cout << green;
          if(useExpNoise)  cout << "  Using exponential degree system noise\n";
          else             cout << "  Using uniform degree system noise\n";
          cout << normText;
          break;
        }
        case ('u'+'m'+'d'): { 
          if(sInputVal=="")  useManualDisplay = 1;
          else  useManualDisplay = ((bool)(argv[i][eqPos+1] - 48));
          if(useManualDisplay)  cout << "  Turning on manual display\n";
          else                  cout << "  Manual display is turned off\n";
          break;
        }
        case ('u'+'s'+'y'+'m'+'m'+'c'+'m'): { 
          if(sInputVal=="")  useSymmetrizedCM = 1;
          else  useSymmetrizedCM = ((bool)(argv[i][eqPos+1] - 48));
          if(useSymmetrizedCM)  warningMsg("  Using symmetrized matrix!\n");
          else                  cout << "  Not using matrix symmetrization\n";
          break;
        }
        case ('u'+'z'+'e'): { 
          if(sInputVal=="")  useZeroE = 1;
          else  useZeroE = ((bool)(argv[i][eqPos+1] - 48));
          if(useZeroE)  cout << "  Turning on zero energy moves (if available)\n";
          else          cout << "  Zero energy moves are off\n";
          break;
        }
        case 'v':
          cout << "  Changing verbosity (output quantity) from " 
               << itos(verbosity,1);
          verbosity = atoi(&argv[i][eqPos+1]);
          cout << " to " << itos(verbosity,1) << endl;
          break;
        case ('w'+'c'): { 
          cout << "  Changing connected edge weight from " << weightC;
          weightC = atoi(&argv[i][eqPos+1]);  
          cout << " to " << weightC << endl;
          break;
        }
        case ('w'+'s'): { 
          cout << "  Changing weight step from " << weightStep;
          weightStep = atoi(&argv[i][eqPos+1]);  
          cout << " to " << weightStep << endl;
          break;
        }
        case ('w'+'u'): { 
          cout << "  Changing unconnected edge weight from " << weightU;
          weightU = atoi(&argv[i][eqPos+1]);  
          cout << " to " << weightU << endl;
          break;
        }
        case ('x'): { 
          if(sInputVal=="")  exitwNoRun = 1;
          else               exitwNoRun = (bool)(argv[i][eqPos+1] - 48);
          //cout << "Testing exit flag:  " << (argv[i][eqPos+1]) << endl;
          if(exitwNoRun) 
            cout << redText  << "  User specified to exit before calculations." 
                 << normText << endl;
          else cout << "  Calculations will proceed as normal." << endl;
          break;
        }
        case ('z'+'o'+'u'+'t'): { 
          cout << "  Changing Zout from " << zout;
          zout = (double)atof(&argv[i][eqPos+1]);  
          cout << " to " << zout << endl;
          if(zout<0.0)
            errorMsg("Fixed Zout set to a negative value???");
          break;
        }
        case ('k'+'a'+'v'+'g'):
        case ('z'+'a'+'v'+'g'): {
          cout << "  Changing average degree from " << zavg;
          zavg = (double)atof(&argv[i][eqPos+1]);
          cout << " to " << zavg << endl;
          if(zavg<0.0)
            errorMsg("Average degree set to a negative value???");
          break;
        }
        case ('k'+'m'+'i'+'n'): {
          cout << "  Changing min degree from " << zMin;
          zMin = (double)atof(&argv[i][eqPos+1]);
          cout << " to " << zMin << endl;
          if(zMin<0.0)
            errorMsg("Minimum node degree was set to a negative value?");
          break;
        }
        case ('k'+'m'+'a'+'x'): 
        case ('z'+'m'+'a'+'x'): {
          cout << "  Changing max degree from " << zMax;
          zMax = (double)atof(&argv[i][eqPos+1]);
          cout << " to " << zMax << endl;
          if(zMax<0.0)
            errorMsg("Maximum node degree was set to a negative value?");
          break;
        }
        // for debugging
        case ('n'+'o'+'u'+'t'): {
          outputNodeConnections = atoi(&argv[i][eqPos+1]);
          break;
        }
        case '\0':   // invalid default case
        default  : { // do nothing but notify user of the invalid parameter
          cout << red << "  Invalid command line parameter:  \""
               << argv[i] << "\"... Ignoring." << normText << endl;
          //return 0;
          break;  // break is redundant after return
        }

      } // end command line parameter switch case check
      //cout << "Testing on exit with i = " << i << endl;  // debugging
    } // end command line parameter for loop
  } // end if command line parameter check
  //cout << endl;
/*
  if( !(usePowerInit || useEvenInit || useRandomInit || useSymmetricInit) &&
       verbosity>0)
    warningMsg("  User did not specify an initialization type?  Assuming Symmetric.");
*/
  q = QMax;
  TInitClusterList *pInitMethod;
  TInitSymmetric  initSymmetric(N,1);
  TInitRandom     initRandom(N,q,1);
  TInitEven       initEven(N,q,1);
  TInitPower      initPower(N,nMin,nMax,betaLFR,1);


