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
Peter Ronhovde
August 2005
ML - My Library - Generic utility functions not involving scientific or significant mathematical algorithms.  
Trivial functions such as max() and min() are included here however.

See MSL_Array.h for general comments and some documentation.

Generic Utilities Specific Notes:
*/
#ifndef ML_UTILS_H
#define ML_UTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <sys/times.h>

using namespace std;

#ifdef LINUX_OS
#include <values.h>          // unix version of limits.h
#else
#include <limits.h>
#endif

// double allows easily changing between double and float data types
#ifndef TFLOAT_T
typedef double TFloat; // float type definition
#define TFLOAT_T
#endif

// Define greek characters for WinXP dos window
//char cphi = 232;  char csigma = 229; // standard ascii versions
// In a dos windows, must set the Lucida font and then use "chcp 1253" to 
// change the codepage and then use the command "graftabl"
// Greek characters in a Linux terminal are the same character values as long
// as the terminal is changed to WINDOWS 1253 encoding
const char cphi = 246;   const char csigma = 243; const char cepsilon = 229;
const char comega = 249; const char crho   = 241;

extern  unsigned long long ISeed;     // random number seed

// define some common text ANSI sequences as cstrings (not at all complete)
/*
Black            \e[0;30m	Blue             \e[0;34m
Green            \e[0;32m	Cyan             \e[0;36m
Red              \e[0;31m	Purple           \e[0;35m
Brown            \e[0;33m	Gray             \e[0;37m
Dark Gray        \e[1;30m	Light Blue       \e[1;34m
Light Green      \e[1;32m	Light Cyan       \e[1;36m
Light Red        \e[1;31m	Light Purple     \e[1;35m
Yellow           \e[1;33m	White            \e[1;37m */
const char blackText[11]  = "\033[0;31m";  const char blueText[11]     = "\033[0;34m";
const char greenText[11]  = "\033[0;32m";  const char cyanText[11]     = "\033[0;36m";
const char redText[11]    = "\033[0;31m";  const char magentaText[11]  = "\033[0;35m";
const char brownText[11]  = "\033[0;33m";  const char greyText[11]     = "\033[0;37m";
const char dgreyText[11]  = "\033[1;30m";  const char lblueText[11]    = "\033[1;34m";
const char lgreenText[11] = "\033[1;32m";  const char lcyanText[11]    = "\033[1;36m";
const char lredText[11]   = "\033[1;31m";  const char lmagentaText[11] = "\033[1;35m";
const char yellowText[11] = "\033[1;33m";  const char whiteText[11]    = "\033[1;37m";
// text modification
const char normText[8]    = "\033[0m";     const char boldText[8]      = "\033[1m";
// shorter versions
const char black[11]      = "\033[0;31m";  const char blue[11]         = "\033[0;34m";
const char green[11]      = "\033[0;32m";  const char cyan[11]         = "\033[0;36m";
const char red[11]        = "\033[0;31m";  const char magenta[11]      = "\033[0;35m";
const char brown[11]      = "\033[0;33m";  const char grey[11]         = "\033[0;37m";
const char darkgrey[11]   = "\033[1;30m";  const char lightblue[11]    = "\033[1;34m";
const char lightgreen[11] = "\033[1;32m";  const char lightcyan[11]    = "\033[1;36m";
const char lightred[11]   = "\033[1;31m";  const char lightmagenta[11] = "\033[1;35m";
const char yellow[11]     = "\033[1;33m";  const char white[11]        = "\033[1;37m";
// text modification
const char normal[8]      = "\033[0m";     const char bold[8]          = "\033[1m";
/*
// If one wants to kill the color codes for the console output, delete the 
// above definitions and uncomment these data.
const char black[11]      = "";  const char blue[11]         = "";
const char green[11]      = "";  const char cyan[11]         = "";
const char red[11]        = "";  const char magenta[11]      = "";
const char brown[11]      = "";  const char grey[11]         = "";
const char darkgrey[11]   = "";  const char lightblue[11]    = "";
const char lightgreen[11] = "";  const char lightcyan[11]    = "";
const char lightred[11]   = "";  const char lightmagenta[11] = "";
const char yellow[11]     = "";  const char white[11]        = "";
// text modification
const char normal[8]      = "";     const char bold[8]       = "";
*/
// other useful char arrays
const char spc = ' ';

// Miscellaneous constants for reference
// 1 ly = 1 lightyear = 9.4605284 × 10^15 meters
// 1 pc = 1 parsec = 3.261636263 ly = 3.08568025 × 10^16 m
// 1 Mo = solar mass = 1.98892 × 10^30 kg
// G = 6.67300 × 10^-11 m^3/(kg*s^2) = 4.51737 × 10^-30 pc^3/(Mo*s^2)
// Phi scale 4.301179 × 10^-10 from [Mo/kpc] to [(100 km/s)^2]
// electron charge e = 1.60217646 × 10^-19 C
// speed of light c = 2.99792458 × 10^8 m/s
// 1 GeV/c^2/cm^3 = 2.63332455 × 10^7 Mo/kpc^3 conversion factor

const char Months[12][10] =
 {"January", "February", "March",     "April",   "May",      "June",
  "July",    "August",   "September", "October", "November", "December"};
const char MonthsShort[12][4] =
 {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
const char Days[7][10] =
 {"Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"};
const char DaysShort[7][4] =
 {"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};


const int BigInteger = 2 << 29; // constant large integer

// generic trivial math functions
inline double min(double x, double y) { if(x < y) return x; else return y; }
inline double max(double x, double y) { if(x > y) return x; else return y; }
inline int    min(int x,    int y)    { if(x < y) return x; else return y; }
inline int    max(int x,    int y)    { if(x > y) return x; else return y; }

//inline double abs(double x) { if(x > 0.0) return x; else return -x; }
inline int    abs(int x)    { if(x > 0.0) return x; else return -x; }


inline void echoCMD (int argc, char *argv[]) { // end echoCMDLine
  cout << normText << "Echo command line:  " << greyText << "\"" << argv[0];
  for(int i=1; i<argc; i++) { cout << " " << argv[i]; } // end command echo
  cout << "\"\n\n" << normText ;
} // end echoCMDLine

// function prototypes
string  toLowerCase(const string s);
string  toUpperCase(const string s);
int     charSum(const string s, int convertCase = 0);
string  itos(int I, bool UseThSep = 0, const char ThSep = ',',
            bool usePad = 0, unsigned int padDigits = 0, char padChar = '0');
// only basic operation is implemented in ftos() at the moment
string dtos(double f, short unsigned int NDigits = 5, 
  bool UseThSep = 0, bool UseSciNotation = 0, bool UseMathSciNotation = 0,
  bool UseDecSep = 0, short unsigned int DecSepInterval = 5,
  double UpperSciNotationStart = 9.9999999E9, 
  double LowerSciNotationStart = 1.0E-5, bool UseEuroNumbers = 0);
inline string ftos(float f, short unsigned int NDigits = 6, 
  bool UseThSep = 0, bool UseSciNotation = 0, bool UseMathSciNotation = 0,
  bool UseDecSep = 0, short unsigned int DecSepInterval = 5,
  double UpperSciNotationStart = 9.9999999E9, 
  double LowerSciNotationStart = 1.0E-5, bool UseEuroNumbers = 0) {
 return dtos((double)f,NDigits,UseThSep,UseSciNotation,UseMathSciNotation,
   UseDecSep,DecSepInterval,UpperSciNotationStart,LowerSciNotationStart,
   UseEuroNumbers); 
}
ostream& OutputLineSeparator(ostream& fo, int length = 78,
                  char separatorChar = '-', char commentChar = '%');
template <typename T> void selectSort(T data[], int Length);
template <typename T> void insertionSort(T data[], int Length);
bool  floatEq(double f1, double f2, double epsilon=5.0E-8);


inline void debugPause(string s = "") { 
  if(s!="")  cout << s << flush;
  string pauseString;  
  cin >> pauseString; 
  return;
} // end pause

inline double Phi(const double x1, const double x2) {
  // from opengroup
  return ( erf(x2*M_SQRT1_2) - erf(x1*M_SQRT1_2) ) / 2.0;
} // end Phi
inline double Norm(const double x) {
  //return ( exp(-x*x*0.5)*M_2_SQRTPI*0.5 );
  return ( exp(-x*x*0.5)*sqrt(M_PI*2.0) );
} // end Norm
inline double PhiAS(const double x) {
  // from Abramowitz & Stegun
  double  p = 0.33267, a1 = 0.4361836, a2 = -0.1201676, a3 = 0.9372980;
  double  t = 1.0/( 1.0 + p*x );
  //return ( 1.0 - Norm(x)*(a1*t + a2*t*t + a3*t*t*t) );
  return ( 1.0 - Norm(x)*t*(a1 + t*(a2 + a3*t)) ); // faster, 4 vs. 7 mult.
} // end PhiAS
inline double PhiInv(const double p) {
  // approximate PhiInv by John D'Errico based on Abramowitz & Stegun
  double  c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
  double  d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
  if(p<0.5){
    double t = sqrt( log(1.0/(p*p)) );
    return (-t + (c0 + t*(c1 + c2*t))/(1.0 + t*(d1 + t*(d2 + t*d3))));
  } else {
    double q = 1.0 - p;  // use symmetry for p>=0.5
    double t = sqrt( log(1.0/(q*q)) );
    // return negative of p<0.5 calculation
    return (t - ((c0 + t*(c1 + c2*t))/(1.0 + t*(d1 + t*(d2 + t*d3)))));
  } // end else
} // end PhiInv

// crude include after needed global constants
//#include "../mytimefns.cpp"
//*****************************************************************************
// # ifdef LINUX_OS
  // ------- timing variables -------------------------------------------------
  // single processor time calculation stuff
  //tm *Today;      //Today = new tm;
  //tm *TodayEnd;   //TodayEnd = new tm;
  //time_t theTime=time(NULL);  //Today = localtime(&theTime);
  extern time_t theTime;
  extern tms StartTime;//  times(&StartTime);  // gcc specific time variables
  extern tms CurrTime, EndTime;              // gcc specific time variables

  extern time_t localTime;
  extern tms    localStartTime;//  times(&StartTime);  // gcc specific time variables
  extern tms    localCurrTime, localEndTime;           // gcc specific time variables

  //clock_t UserRunTime, UserInitTime, SysRunTime, SysInitTime;
  //clock_t startClock, currClock, endClock;
  //startClock = clock();
  //sysconf(_SC_CLK_TCK);   // determine the clock ticks per second

  //double runTime  = 0.0;  // run time in seconds
  //double initTime = 0.0;  // time to do preliminary initializations in seconds
  //cout << "CLOCKS_PER_SEC = " << itos(CLOCKS_PER_SEC,1) << "\n";
  // -------------------------------------------------------------------------

inline void initTimes(tms &StartTime, time_t &theTime) {
  theTime=time(NULL);  //Today = localtime(&theTime);
  times(&StartTime);  // gcc specific time variables
  //clock_t startClock, currClock, endClock;
  //startClock = clock();
  //sysconf(_SC_CLK_TCK);   // determine the clock ticks per second

  return;
} // end initTimes

inline void localInitTime() {
  localTime=time(NULL);  //Today = localtime(&theTime);
  times(&localStartTime);  // gcc specific time variables
  //clock_t startClock, currClock, endClock;
  //startClock = clock();
  //sysconf(_SC_CLK_TCK);   // determine the clock ticks per second

  return;
} // end localInitTimes

inline void outputRuntime(tms &StartTime, tms &EndTime, tms &CurrTime) {
  double runTime  = 0.0;  // run time in seconds
  double initTime = 0.0;  // time to do preliminary initializations in seconds

  // ------ Calculate total run time info arrays -----------------------------
  times(&EndTime);
  clock_t UserRunTime, UserInitTime, SysRunTime, SysInitTime;

  // Dates and times here are gcc specific under Linux.
  UserInitTime = CurrTime.tms_utime - StartTime.tms_utime;
  SysInitTime  = CurrTime.tms_stime - StartTime.tms_stime;
  UserRunTime  = EndTime.tms_utime  - StartTime.tms_utime;
  SysRunTime   = EndTime.tms_stime  - StartTime.tms_stime;

  //endClock = clock(); // clock() is not working???

  // following seems to work, but I don't know why I need the factor of 10,000
  initTime = (double)(UserInitTime+SysInitTime)/(double)CLOCKS_PER_SEC*10000.0;
  runTime  = (double)UserRunTime/(double)CLOCKS_PER_SEC*10000.0;  // in secs
  //cout << "Total initalization time was:  " << initTime << " and ";
  cout << brownText << "Total run time was " << runTime << " secs" << normText;
  // integer runTime in seconds - helps with div and mod below
  double fracSec = runTime - floor(runTime);
  int intTime = (int)floor(runTime); 
  if(runTime>60.0) {
    cout << ", or ";
    if(intTime>86400){ cout << intTime/86400 << " days ";  intTime %= 86400; }
    if(intTime>3600) { cout << intTime/3600  << " hrs ";   intTime %= 3600;  }
    if(intTime>60)   { cout << intTime/60    << " mins ";  intTime %= 60;    }
    cout << ((double)intTime+fracSec) << " secs\n";
  } // end if runTime

  cout << endl;
  // -------------------------------------------------------------------------
  return;
} // end outputRuntime

inline double getRuntime(tms &StartTime, tms &EndTime, tms &CurrTime, bool bOutput=0) {
  double runTime  = 0.0;  // run time in seconds
  double initTime = 0.0;  // time to do preliminary initializations in seconds

  // ------ Calculate total run time info arrays -----------------------------
  times(&EndTime);
  clock_t UserRunTime, UserInitTime, SysRunTime, SysInitTime;

  // Dates and times here are gcc specific under Linux.
  UserInitTime = CurrTime.tms_utime - StartTime.tms_utime;
  SysInitTime  = CurrTime.tms_stime - StartTime.tms_stime;
  UserRunTime  = EndTime.tms_utime  - StartTime.tms_utime;
  SysRunTime   = EndTime.tms_stime  - StartTime.tms_stime;

  //endClock = clock(); // clock() is not working???

  // following seems to work, but I don't know why I need the factor of 10,000
  initTime = (double)(UserInitTime+SysInitTime)/(double)CLOCKS_PER_SEC*10000.0;
  runTime  = (double)UserRunTime/(double)CLOCKS_PER_SEC*10000.0;  // in secs
  if(bOutput) {
    //cout << "Total initalization time was:  " << initTime << " and ";
    cout << brownText << "Total run time was " << runTime << " secs" << normText;
    // integer runTime in seconds - helps with div and mod below
    double fracSec = runTime - floor(runTime);
    int intTime = (int)floor(runTime); 
    if(runTime>60.0) {
      cout << ", or ";
      if(intTime>86400){ cout << intTime/86400 << " days ";  intTime %= 86400; }
      if(intTime>3600) { cout << intTime/3600  << " hrs ";   intTime %= 3600;  }
      if(intTime>60)   { cout << intTime/60    << " mins ";  intTime %= 60;    }
      cout << ((double)intTime+fracSec) << " secs\n";
    } // end if runTime

    cout << endl;
  } // end if bOutput

  // -------------------------------------------------------------------------
  return runTime;
} // end outputRuntime


inline void outputTime(double timeInSecs, string sRunType = "run") {
  // integer runTime in seconds - helps with div and mod below
  double fracSec = timeInSecs - floor(timeInSecs);
  int intTime = (int)floor(timeInSecs); 
  cout << brown << "Total " << green << sRunType << brown << " time was "
       << magenta << timeInSecs << brown << " secs" << normText;
  if(timeInSecs>60.0) {
    cout << ", or ";
    if(intTime>86400){ cout << intTime/86400 << " days ";  intTime %= 86400; }
    if(intTime>3600) { cout << intTime/3600  << " hrs ";   intTime %= 3600;  }
    if(intTime>60)   { cout << intTime/60    << " mins ";  intTime %= 60;    }
    cout << ((double)intTime+fracSec) << " secs";
  } // end if runTime
} // end outputTime

inline double getLocalTime() {
  // returns time elapsed in sec since last local* initTimes or timeMsg call
  double tElapsed = getRuntime(localStartTime,localEndTime,localCurrTime,0);
  localInitTime();
  //if(bFlush)  cout << flush;
  return tElapsed; 
} // end msg

inline string dateString(bool bLong=1) {
  time_t t = time(0);
  struct tm *lt = localtime(&t);
  string s = "";
  if(bLong)  s += Months[lt->tm_mon];
  else       s += MonthsShort[lt->tm_mon];
  // add day and year
  s += " " + itos(lt->tm_mday) +", " + itos(lt->tm_year + 1900); 
/*
  cout << std::setfill('0');
  cout << std::setw(4) << lt->tm_year + 1900
       << std::setw(2) << lt->tm_mon + 1
       << std::setw(2) << lt->tm_mday
       << std::setw(2) << lt->tm_hour
       << std::setw(2) << lt->tm_min
       << std::setw(2) << lt->tm_sec
       << std::endl;
*/
  return s;
}

inline string timeString(bool bLong=1) {
  time_t t = time(0);
  struct tm *lt = localtime(&t);
  string s;

  s = itos(lt->tm_hour) + ":" + itos(lt->tm_min,0,0,1,2);
  if(bLong)  s += ":" + itos(lt->tm_sec,0,0,1,2); // add seconds
  return s;
}


inline string dateTimeString(bool bLong=1) {
  time_t t = time(0);
  struct tm *lt = localtime(&t);
  return (dateString(bLong) + "  " + timeString(bLong));
}
inline string csvDateTime(bool bLong=1, char endChar = '\n') {
  return ("\"Date:  "+dateString(bLong)+" at "+timeString(bLong)+"\""+endChar);
}

// # else
  // not sure what to do in windows
// # endif
// end time functions
// ****************************************************************************

//inline void errorMsg(const string s) {
inline void errorMsg(string s) {
  cout << redText << "Error:  " << s << " Aborting." << normText << endl;
  exit(1);  return;
} // end errorMsg
//inline void abortMsg(const string s) {

inline void abortMsg(string s) {
  cout << redText << s << " Aborting." << normText << endl;
  exit(1);  return; 
} // end abortMsg
//inline void warningMsg(const string s, const char colorCode[] = redText) {

inline void warningMsg(string s, const char colorCode[] = red) {
  cout << colorCode << s << normText << endl;  return; 
} // end warningMsg

inline void debugMsg(string s, const char color[] = green, int debugMode = 1) {
  // dm - debug mode level:  0 - none, 1 - basic level, 2 - higher...
  if(debugMode>0)  cout << color << s << normText << flush;
  return; 
} // end debugMsg

inline void msg(string s, const char colorCode[] = green) {
  cout << colorCode << s << normText << flush;
  //if(bFlush)  cout << flush;
  return;
} // end msg

inline double timeMsg(string sBefore, string sAfter = "", 
                      const char colorCode[] = brown) {
  // returns time elapsed in sec since last local* initTimes or timeMsg call
  double tElapsed = getRuntime(localStartTime,localEndTime,localCurrTime,0);
  cout << colorCode << sBefore << "(t=" << tElapsed << ")" 
       << sAfter << normText << flush;  
  localInitTime();
  //if(bFlush)  cout << flush;
  return tElapsed; 
} // end msg



// Initialize random generator (seeding is automatic)

#ifdef USE_RANDOMLIB
//class RandomLib;
// MersenneTwister random numbers class
//#include "/home/ronhovde/RandomLib/Random.hpp"
//# include "/home/peter/RandomLib/Random.hpp"


extern RandomLib::Random r;
inline double randomDouble() { // a floating point RNG
  return (double)r.Fixed(); }
inline double randomDouble(double scale) { // a floating point RNG
  return (double)r.Fixed()*scale; }
inline int randomInt(int iStart, int iEnd) { // an integer RNG
  return r.IntegerC(iStart,iEnd); }
inline void rngReseed(unsigned long long ISeed) {
  cout << magenta << "Reseeding the RandomLib seed" << endl; // debugging
  r.Reseed(ISeed);  return; }
// RandomLib uses a string reseed or an obscure (to me) vector value
inline void rngReseed(string ISeed) { 
  cout << magenta << "Reseeding the RandomLib seed" << endl; // debugging
  r.Reseed(ISeed);  return; }
inline string rngSeedString() { return r.SeedString(); }
inline void initRNG(unsigned long long ISeed) {
  if(ISeed>0) {
    r.Reseed(ISeed);
    cout << green << "  User initialized RandomLib seed to: "
         << brown << r.SeedString() << normText << endl;
  } // end if Iseed
  else  cout << green << "  RandomLib seed automatically set to: "
             << r.SeedString() << normText << endl;

  cout << grey << "  Here are two test numbers: " << r.Fixed() << " and " 
       << r.Fixed() << normText << endl;  // debugging
} // end initRNG()

#else 

// define standard random number generators using srand() and rand()
inline double randomDouble() {
  // a floating point random number generator
  return (double)rand()/((double)RAND_MAX+1.0); }
inline double randomDouble(double rScale) { // generic double rng
  return (double)rand()/((double)RAND_MAX+1.0)*rScale; }
inline double randomDouble(double low, double high) { // generic double rng
  return (double)rand()/((double)RAND_MAX+1.0)*(high-low) + low; 
}

inline int randomInt(int iStart, int iEnd) { // generic double rng with scale
  return (int)( (double)rand()/((double)RAND_MAX+1.0)
               *(double)(iEnd-iStart) + (double)iStart  + 0.5); }
inline double randomNormalCLT(double mean=0.0, double sigma=1.0) {
  // a crude normal generator using 12 average/shifted uniformly distributed
  // numbers and the Central Limit Theorem
  double  rSum = 0.0;
  for(int i=0; i<12; i++)  rSum += randomDouble();
  rSum -= 6.0;  // now this is approximately a N(0,1) random number
  // now scale according to mean and sigma
  return (mean + sigma*rSum);
} // end randomNormalCLT
inline double randomExp(double mean=1.0) {
  // base e log, exponential variate
  return ( -mean*log(randomDouble()) );
} // end randomExp
//inline double randomExp(double lambda=1.0) {
  // base e log, exponential variate
//  return ( -log(randomDouble())/lambda );
//} // end randomExp
inline double randomNormalPolar(double mean=0.0, double sigma=1.0) {
  // a normal generator using polar approximation
  // from finance.bi.no/~bernt/gcc_prog/recipes/recipes/node23.html
  double u1, u2, v1, v2, s=2.0, x1;
  while(s>=1.0) {
    u1 = randomDouble();  u2 = randomDouble();
    v1 = 2.0*u1 - 1.0;    v2 = 2.0*u2 - 1.0;
    s = v1*v1 + v2*v2;
  } // end while
  x1 = v1*sqrt( -2.0*log(s)/s );
  // scale according to mean and sigma
  return ( mean + sigma*x1 );
} // end randomNormalPolar
inline double randomNormal(double mean=0.0, double sigma=1.0) {
  // a normal generator using Abramowitz & Stegun approximation
  return ( mean + sigma*PhiInv(randomDouble()) );
  // scale according to mean and sigma
} // end randomNormalCLT
inline double randomPower(double min, double max, double alpha) {
  // x is the power-law variate
  double U = randomDouble(), x;
  if(floatEq(alpha,-1.0)) {
    // special case of 1/x power law
    x = min*exp( U*log(max/min) );
  } else {
    // else calculate a general power-law variate
    double  alphap1 = alpha+1.0, xp, minp, maxp;
    xp = pow(U,alphap1);  minp = pow(min,alphap1);  maxp = pow(max,alphap1);
    x = pow( U*(maxp - minp) + minp , 1.0/alphap1 );
  } // end else
  return x;  
} // end randomPower
inline int randomPowerInt(int min, int max, double alpha) {
  return ((int)( randomPower((double)min,(double)max,alpha) + 0.5 ));  
} // end randomPower
inline double randomPowerMean(double min, double max, double alpha) {
  double mean, Ck;
  if(floatEq(alpha,-2.0)) {
    Ck = (alpha+1)/( pow(max,alpha+1) - pow(min,alpha+1) );
    mean = Ck*log(max/min);
  } else if(floatEq(alpha,-1.0)) {
    Ck = 1.0/log(max/min);
    mean = Ck*( pow(max,alpha+2) - pow(min,alpha+2) )/ (alpha+2);
    errorMsg("power law mean has not been tested for alpha <> -2");
  } else {
    // all other cases
    Ck = (alpha+1)/( pow(max,alpha+1) - pow(min,alpha+1));
    mean = Ck*( pow(max,alpha+2) - pow(min,alpha+2) )/ (alpha+2);
    errorMsg("power law mean has not been tested for alpha <> -2");
  } // end else
  return mean;
} // end randomPowerMean


inline string csv(double value, int digits = 6) { 
  if(isnan(value))  errorMsg("csv() shows that the value is nan");
  if(isinf(value))  errorMsg("csv() shows that the value is nan");
  return ( ftos( value , digits ) + "," );
} // end double csv string
inline string csvH(string s)     { return ( "\""+s+"\","                  ); } 
inline string csvHeol(string s)  { return ( "\""+s+"\"\n"                 ); } 
inline string csvSH(string s)    { return ( "\""+s+"\",\"sigma_"+s+"\","  ); }
inline string csvSHeol(string s) { return ( "\""+s+"\",\"sigma_"+s+"\"\n" ); }
/*
inline string csvSH(vector<string &s, bool bEndEOL) { 
  // stats header in csv format including the value and sigma_value
  string h("");  
  for(int i=0; i<s.size()-1; i++)  h += "\""+s[i]+"\",\"sigma_"+s[i]+"\",";
  if(bEndEOL)  h += "\""+s[i]+"\",\"sigma_"+s[i]+"\"\n";
  else         h += "\""+s[i]+"\",\"sigma_"+s[i]+"\",";   // end as regular csv
  return h;
} // end double csv string
*/
inline string csv(int value) { 
  return ( itos( value ) + "," );
} // end double csv string
inline string csveol(double value, int digits = 6) { 
  return ( ftos( value , digits ) + "\n" );
} // end double csv string
inline string csveol(int value) { 
  return ( itos( value ) + "\n" );
} // end double csv string
inline string cssv(string d, int value) { 
  return ( "\"" + d + itos( value ) + "\"," );
} // end double css string (in quotes)
inline string cssv(string d, double value, int digits = 6) { 
  return ( "\"" + d + ftos( value , digits ) + "\"," );
} // end double css string (in quotes)
inline string cssveol(string d, int value) { 
  return ( "\"" + d + itos( value ) + "\"\n" );
} // end double css string (in quotes)
inline string cssveol(string d, double value, int digits = 6) { 
  return ( "\"" + d + ftos( value , digits ) + "\"\n" );
} // end double css string (in quotes)

inline void rngReseed(unsigned long long ISeed) { // generic integer rng
  //cout << magenta << "Reseeding the srand seed" << endl;  // debugging
  srand(ISeed);  return; }
inline unsigned long long getRNGSeed() { return ISeed; }
inline string rngSeedString() { 
  //char buffer[20];
  //string s = (string)ulltoa(ISeed,buffer);  // ???
  //return s;
  // use a crude double conversion for now
  string s = dtos( (double)ISeed );
  s.erase(s.find('.'));
  return s; 
}
inline string csvSeedString(char endChar = '\n') {
  return ( "\"RNG seed = " + rngSeedString() + " using srand()\"" + endChar );
}
inline void initRNG(unsigned long long ISeed) {
  if(ISeed>0) {
    srand(ISeed);
    cout << green << "  User initialized rand() seed to: " << brown 
         << itos(ISeed) << normText << endl;
  } // end if ISeed
  else { // initialize the seed automatically
/*
    streampos IStart = (int)fmod((double)time(NULL),(double)10000.0);
    cout << brown << "  RNG seed position is " << IStart << "\n";
    ifstream  fin("/home/peter/truerandoms_2007-12-26.bin",ios::in | ios::binary);
    cout << brown << "  Current file position is " << fin.tellg() << "\n";
    if(fin.good()) {
      fin.clear();
      cout << brown << "  Current file position is " << fin.tellg() << "\n";
      if(fin.eof())  cout << "  At end of file\n";
      else           cout << "  Not at end of file yet\n";
      long testInt;  fin >> testInt;  fin >> testInt; // fin >> testInt;
      if(fin.eof())  cout << "  At end of file\n";
      else           cout << "  Not at end of file yet\n";
      cout << brown << "  Current file position is " << fin.tellg() << "\n";
      fin.seekg(1,ios::beg);
      cout << brown << "  Current file position is " << fin.tellg() << "\n";
      //int testInt;   fin >> testInt;
      //fin >> testInt; 
      fin >> ISeed;  // now read the seed value at the 'random' location
      cout << brown << "  RNG test seed is " << testInt << "\n";
      srand(ISeed);
      cout << green << "  srand() seed automatically set to: "
           << itos(ISeed) << normText << endl;
      fin.close();
    } else cerr << red << "  Could not open RNG seed file" << normText << "\n";
*/
      srand(ISeed);
  } // end else
} // end initRNG()

#endif

#endif

