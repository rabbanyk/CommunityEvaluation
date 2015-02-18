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
#include <string>
#ifndef ML_UTILS_H
#include "./ML_Utils.h"
#endif
/*
Peter Ronhovde                                     May 2007
MSL - My (or Math) Scientific Library - TStats1D base class

Comments:
- All memory allocation is dynamic.
- Data type converts all input to floats (or doubles)
- The data array uses my MSL_Array1D array class.

- This is not production code, so it is assumed that the user will have some level 
of knowledge about how to handle cases where there is a potential "dummy" type
problem.  It is not sloppy code for the most part, just that there are a few places
where adding all of the necessary validation code would slow down the application.
It also reduces the requirements of coding the classes and reduces the class size 
somewhat. The user should be careful and knowledgable when manually making 
changes to member data elements.
- Since this package is developed with scientific applications in mind that tend
to use large arrays, value based operators [-(unary), + and -, * and /(scalar)]
should only be used with (very) small arrays.
- There is an issue with validation functions and validation checks taking up
*way* too much computational time.  If full class error checking is used 
(including validating array references and repeated callings of the isValid() 
function), the performance can be reduced by as much as a factor of 4!  
This is particularly a problem with validation on Array element referencing.  
The original code was constructed so that StrictErrorChecking could be set 
to 1 to eliminate the checks with the compiler removing the dead code; 
however, I had direct evidence that the older Borland compiler I had was 
*not* removing the trivially assessed dead/unreachable code.  Therefore, 
I had to manually eliminate them.  It is preferrable to reimplement this if a 
compiler is being used that with remove the unreachable code (and also 
the associated if statement check?).  Otherwise, I might have to implement 
the variation in the code with preprocessor commands which is much less 
preferrable.

General (but not strict) Coding Style:
- Constant values begin with a capital letter.
- Variables and functions usually begin with a lower case letter with the first
letter of subsequent whole words being capitalized i.e. isValid() or nParts.
However, this is sometimes violated for mathematical objects such as a matrix
equation Ax=b where the capital and lower case letters have other significance,
or where it may be otherwise confusing from a science/math point of view to 
follow the coding style.  The context of the program should make it clear when 
this is the case.
- Pointers begin with a single lower case 'p' as in pA or pPart1.
- A variable that stores a number (rows, objects...) begins with an 'n' or 'N'.
- iteration variables often begin with an 'i'.  Of course, standard iteration variables
are i, j, and k and sometimes 'l' (el).  Effort is made to keep i refering to the column,
j referring to the row.  The reason it was done this way was that developement 
was begun with a 2D Poisson's Equation solver in mind where i was used as the 
index for the r direction.  It then remained for the general implementation of the 
array classes.  k typically refers to the 3D level and 'l' (el) is a general index used 
to eliminate to eliminate possible confusion or when the algorithm specifically 
references it.
- string variables often begin with an 's' or 'S'.
*/

/* To-Do List for TStats1D):  *********************************************************
Items are not in order of importance.
1. Need to write a better sort routine and allow an ascending or descending parameter option.
2. Write a 'output stats function' so that we don't have to manually do so all the time.
3. Need to write the operator=() function.
4. Median needs to be checked and debugged.
5. Need to implement skew and kurtosis values and also add them to the output
function output.
6. It would be nice to implement an insertion operator in addition to the
output function.
7. Bin test is currently only valid for data sets between 0.0 and 1.0, need to fix.
8. Validation checks are not working and are therefore disabled.
9. Check for copy and constructor problems if data is of length data is zero 
such as when noDataStore is set to true.
10. Negated 'noDataStore' makes logic confusing - fix it
**************************************************************************************/

/* Future 'Wish' or 'Idea'-List:  **************************************************************************************
1. 
  
**************************************************************************************/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TSTATS1D_H
#define MSL_TSTATS1D_H

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
//#ifndef MSL_ARRAY1D_H
//#include "MSL/MSL_Array1D.h"
//#endif

namespace MSL {

const int MaxNDataElements = 1024*1024;

#ifndef MSL_double_T
// double allows easily changing between double and float data types
//typedef double TFloat;     // float type definition
#define MSL_double_T
#endif

/*
class TStats {
 public:
  // Data elements
  static  bool StrictErrorChecking;
  static  bool VerboseErrorReports;
  // User member functions

  // Standard member functions
  virtual bool isValid() const = 0;
  TStats() {};
  TStats(TStats& a) {};TStats1D
  ~TStats() {};

 protected:
}; // end virtual class TStats
//------------------------------------------------------------------------------
// initialize static variables
bool TStats::StrictErrorChecking = 1;
bool TStats::VerboseErrorReports = 0;
//------------------------------------------------------------------------------
*/

//------------------------------------------------------------------------------
class TStats1D {
 public:
  inline double& operator[](int i) const;     // standard C/C++ array reference
  // User member functions
         double  average();  
  inline double  avg();     // both are standard averages
         double  sigma();   // regular 'N-1' standard deviation
         double  sigmaN();  // sigma 'N' value rather than 'N-1' value
  inline double  skew();    // currently invalid
  inline double  kurt();    // currently invalid
         double  sum();     
         double  sumSq();   // if a user wants to work with stats manually
         double  max();     
         double  min();     
  inline int     getN() { return nData; };
  // csv() outputs a "comma separated value" list of important data
  //   wData = 0 outputs:  avg() + endChar
  //           1 outputs:  avg() + separator + sigma() + endChar
  inline string  csv(int wData = 1, char endChar = ',', char separator = ',');
  inline string  csvHeader(int wData = 1, char endChar = ',', char separator = ',');
  inline string  csvH(int wData = 1, char endChar = ',', char separator = ',');
  
  // less used statistics
  //double  mode()   { cout << "Mode is not currently implemented"   << endl; };
  double  median(bool inPlace);  // calculate median - force inPlace sort?

  // Other useful member functions
         void  store(double value);   // add a data element to the class
         void  clear();               // delete and clear all data
         void  calc();                // user forced stats (re)-calculation -
                                      //   normally not needed
         void  shift(double val);     // shift all data by val
  inline void  resize(int newSize, double fillValue=0.0); // resize nData elements
         bool  sort(int sortMethod=3);  // sort the current set of data
               // perform a bin test on the data
         void  binTest(int nBins=10, double xMin=0.0, double xMax=1.0); 
         void  output(int verbosity=0); // output stats in standard format
  

  // Standard member functions
  virtual inline bool isValid();  // overloading to make this a genuine class
                      TStats1D(string description);
                      TStats1D(int MaxSize=0, bool noDataStore=1, string description="");
                      TStats1D(TStats1D& a);
  virtual            ~TStats1D(); // compiler warning says that destructor should 
                                  // be virtual if there are virtual functions?

 string    desc;
 protected:
  int      nData;       // the number of data elements in the counter
  //Array1D  data;        // the actual data storage array in double data type
  double  *data;
  bool     validData;   // is the data set valid?
  bool     statsDirty;  // uses pre-calculated values when no change has occurred to data
  bool     unordered;   // Stats1D class does *not* change data order by itself
  bool     noDataStore; // does TStats1D object store all data or just running sums
                        //   (can only be set in initial constructor)
  // pre-calculated data variables to be used if called multiple times
  // (don't want to recalculate everything on every function call)
  double   *pAvg, *pSigma, *pSkew, *pKurt;
  double   *pSum, *pMax,  *pMin;
  long double  *pSumSq; // use long double only internal, return value as double
}; // end class TStats1D


inline TFloat  TStats1D::avg() { return average(); };
inline string  TStats1D::csv(int wData, char endChar, char separator) {
  //cout << desc << " average = " << flush << average() 
  //<< " and sigma = " << flush << sigma() << "\n" << flush;  // debugging
  if(isnan(average())) {
    warningMsg("TStats1D.csv() for "+desc+" shows that the average is nan?");
    return ( (string)"\"nan\"" + separator + "\"nan\"" + endChar );
  } // end if nan
  if(isinf(average())) {
    warningMsg("TStats1D.csv() for "+desc+" shows that the average is inf?");
    return ( (string)"\"inf\"" + separator + "\"inf\"" + endChar );
  } // end if nan
  //cout << average() << "  " << flush; // debugging
  if(wData==0)  return ( ftos(average()) + endChar );
  else if(wData==1) {
    if(nData>2) return ( ftos(average()) + separator + ftos(sigma()) + endChar);
    // catch case of no valid sigma()
    //else        return ( ftos(average()) + separator + "\"n/a\"" + endChar  );
    else        return ( ftos(average()) + ",\"n/a\"," );
  } else //if(wData>2)  return ( avg() + separator + sigma() + endChar )
    warningMsg("TStats1D::csv() output is not defined for input wData = "
               +itos(wData)+"\n");
  return ( "\"Error in TStats1D::csv() output\"," ) ;
}; // end csv string output
inline string  TStats1D::csvHeader(int wData, char endChar, char separator) {
  string header;
  if(wData==0)  header = "\"" + desc + "\"" + endChar;
  else if(wData==1)
    header = "\"" + desc + "\"" + separator + "\"sigma_" + desc + "\""+endChar;
  else //if(wData>2)  return ( avg() + separator + sigma() + endChar )
    warningMsg("TStats1D::csvHeader() output is not defined for input wData = "
               +itos(wData)+"\n");
  return header;
}; // end csv string output
inline string  TStats1D::csvH(int wData, char endChar, char separator) {
  return  csvHeader(wData,endChar,separator);  
};

inline double TStats1D::skew() { 
  warningMsg("Skewness is not currently implemented\n");  return 0.0; 
}; // end skew
inline double TStats1D::kurt() { 
  warningMsg("Kurtosis is not currently implemented\n");  return 0.0; 
}; // end kurtosis
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
class TStatsArray {
// Memory allocation is dynamic and uses standard C++ array indexing with 
// brackets as in arrayName[index] with zero indexing.

 public:
 // overloaded operators
           //TClusterList& operator=(const TClusterList& b);
  inline   TStats1D& operator[](int i) const { return *stats[i]; };
  // Public member functions
  inline   int    add(TStats1D &s);       // add a stats variable 
  inline   int    add(TStats1D *ps);      // add a stats variable by pointer
  inline   void   clear();  // clear all stats 
  inline   void   erase();  // erase stats array without affecting added stats
  //inline int    resize(int NewSize, bool bKeepData=1);
  //inline int    move(TClusterList &c);
  //inline void   move(int i, int j, TBool bFast=1);
  //inline int    erase(int i, TBool bFast=1);
  //inline void   swap(int i, int j);
  // data acquisition functions
  inline   int    getSize()    const { return nStats;    };
  inline   int    getMaxSize() const { return NMaxStats; };
  inline   string csvString();
  inline   string csvHeader();
 
  // Standard public member functions
  // the listSize = -1 is a flag that tell the constructor to assume that all
  // lists are valid and random access is expected
  TStatsArray(int MaxSize = 20) {
    NMaxStats = MaxSize;  stats = new TStats1D*[NMaxStats];  erase(); };
  ~TStatsArray() { delete[] stats;  nStats = 0; }; // does not destroy added stats!
  
 //protected:
  // Data elements
  int       nStats;          // the total number of cluster lists stored
  int       NMaxStats;       // the maximum total number of cluster lists

  // the cluster list data is stored as "2D" dynamic array to prevent a lot of 
  //   data movement during move operations
  TStats1D  **stats;         // an array of pointers
}; // end class TClusterList;

// --- inline functions ------------------------------------------------------
inline void TStatsArray::clear() {
  // reset all stats
  for(int i=0; i<nStats; i++)  stats[i]->clear();
  return; 
}; // end clear


inline void TStatsArray::erase() {
  // erase all stats
  for(int i=0; i<NMaxStats; i++)  stats[i] = NULL;
  nStats = 0;
  return; 
}; // end erase


inline int TStatsArray::add(TStats1D &s) {
  if(nStats==NMaxStats)  
    errorMsg("Cannot add another stats variable in TStatsArray::add()");
  stats[nStats] = &s;
  nStats++;
  return nStats; 
}; // end add
inline int TStatsArray::add(TStats1D *ps) {
  if(nStats==NMaxStats)
    errorMsg("Cannot add another stats variable in TStatsArray::add(*)");
  stats[nStats] = ps;
  nStats++;
  return nStats; 
}; // end add by pointer


inline string TStatsArray::csvHeader() {
  string s = "";
  for(int i=0; i<nStats-1; i++)  s += stats[i]->csvHeader();
  s += stats[nStats-1]->csvHeader(1,'\n');  // includes end of line
  return s; 
}; // end csvHeader


inline string TStatsArray::csvString() {
  string s = "";
  for(int i=0; i<nStats-1; i++)  s += stats[i]->csv();
  s += stats[nStats-1]->csv(1,'\n');
  return s; 
}; // end csvString

} // end namespace MSL
#endif
