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
- iteration variables often begin with an 'i'.  Of course, standard iteration 
variables are i, j, and k and sometimes 'l' (el).  Effort is made to keep i 
refering to the column, j referring to the row.  
The reason it was done this way was that developement was begun with a 2D 
Poisson's Equation solver in mind where i was used as the index for the r 
direction.  It then remained for the general implementation of the array classes.  
k typically refers to the 3D level and 'l' (el) is a general index used 
to eliminate to eliminate possible confusion or when the algorithm specifically 
references it.
- string variables often begin with an 's' or 'S'.
*/

/* To-Do List for TStats1D):  **************************************************
Items are not in order of importance.
1. Need to write a better sort routine and allow an ascending or descending 
   parameter option.
2. Write a 'output stats function' so that we don't have to manually do so all 
   the time.
3. Need to write the operator=() function.
4. Median needs to be checked and debugged.
5. Need to implement skew and kurtosis values and also add them to the output
function output.
6. It would be nice to implement an insertion operator in addition to the
output function.
7. Bin test is currently only valid for data sets between 0.0 and 1.0, needs fix
8. Validation checks are not working and are therefore disabled.
9. Check for copy and constructor problems if data is of length data is zero 
such as when noDataStore is set to true.
10. Negated 'noDataStore' makes logic confusing - fix it
*******************************************************************************/

/* Future 'Wish' or 'Idea'-List:  ********************************************************************************
1. 
  
*******************************************************************************/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
//#ifndef MSL_ARRAY1D_H
//#include "MSL/MSL_Array1D.h"
//#endif
#ifndef MSL_TSTATS1D_H
#include "MSL_Stats1D.h"
#endif

namespace MSL {
/*
//------------------------------------------------------------------------------
// initialize static variables
bool TStats::StrictErrorChecking = 1;
bool TStats::VerboseErrorReports = 0;
//------------------------------------------------------------------------------
*/


TFloat  TStats1D::average() { 
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    return *pAvg;  // return whether is was unchanged or newly calculated here
  } // end if isValid
  else cerr << "Error in average function for TStats1D object.\n";
  return 0.0;  // invalid if we make it to here
}; // end average


TFloat  TStats1D::median(bool inPlace) {
  // calculate median - force inPlace sort which destroys original data order
/*
  TFloat *pData; 
  if(inPlace) { 
    sort(3); 
    pData = &data[0]; 
  } // end if inPlace
  else { 
    // generate dummy array to sort
    TStats1D tempData = *this;  // crudely copies entire *this
    tempData.sort(3);
    pData = &(tempData.data[0]);
  } // end else inPlace


  // now continue with median calculation in either case with pData pointing
  // to the correct data set 
  if((nData%2)==1) return pData[nData/2 + 1];  // integer division in []'s
  else return ( pData[nData/2] + pData[nData/2 + 1] )/2.0;
*/
  return 0.0;
};  

TFloat  TStats1D::sigma() { 
  // Standard deviation of current Array1D
  //TFloat a2sum = 0.0, result;
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    return *pSigma;  // return whether is was unchanged or newly calculated here
  } // end if isValid
  else cerr << "Invalid sigma function for TStats1D object.\n";
  return -1.0;  // invalid if we make it to here
}; // end sigma


TFloat  TStats1D::sigmaN() { 
  // Standard deviation of current Array1D divided by 'N-1' rather than 'N'
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    // return whether is was unchanged or newly calculated here
    return ( (*pSigma)*(TFloat)((nData-1)*(nData-1))/(TFloat)(nData*nData) );  
  } // end if isValid
  else cerr << "Invalid sigmaN function for TStats1D object.\n";
  return -1.0;  // invalid if we make it to here
}; // end sigmaN


TFloat  TStats1D::sum() { 
  if(isValid()) {  // end if nData
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    return *pSum;  // return whether is was unchanged or newly calculated here
  } // end if isValid
  else cerr << "Invalid sum function for TStats1D object.\n";
  return 0.0;  // invalid if we make it to here
}; // end sum


double  TStats1D::sumSq() { // use long double only internally
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    // return whether is was unchanged or newly calculated here
    return (double)(*pSumSq);
  } // end if isValid
  else cerr << "Invalid sumSq function for TStats1D object.\n";
  return 0.0;  // invalid if we make it to here
}; // end sumSq


TFloat  TStats1D::max() { 
  // Standard deviation of current Array1D
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    return *pMax;  // return whether is was unchanged or newly calculated here
  } // end if isValid
  else cerr << "Invalid max function for TStats1D object.\n";
  return 0.0;  // invalid if we make it to here
}; // end sigma


TFloat  TStats1D::min() { 
  // Standard deviation of current Array1D
  if(isValid()) { 
    if(statsDirty)  calc();  // calculate all stats if data set is dirty
    return *pMin;  // return whether is was unchanged or newly calculated here
  } // end if isValid
  else cerr << "Invalid max function for TStats1D object.\n";
  return 0.0;  // invalid if we make it to here
}; // end sigma


void  TStats1D::calc() { 
  // Standard deviation of current Array1D
  TFloat asum, a2sum;
  //cout << redText << "In TStats calculation\n" << normText;  // debugging

  if(isValid() && nData>2) { 
    if(!noDataStore) { 
      // sums are NOT already calculated so we must calculate them here
      TFloat* pElement = &data[0];
      // calculate sums
      asum  = (*pElement);                // initialize asum
      a2sum = (*pElement)*(*pElement);    // initialize a2sum
      // initialize min and max (note this is NOT done in constructor)
      if(pMax==NULL)  pMax = new TFloat;
      if(pMin==NULL)  pMin = new TFloat;
      *pMax = (*pElement);
      *pMin = (*pElement);
      for(int i=1; i<nData; i++) { 
        pElement++;                       // assumes array[] data format
        asum  += (*pElement);  
        a2sum += (*pElement)*(*pElement);  
        // and some other stuff
        if((*pElement)>(*pMax))  *pMax = (*pElement);
        if((*pElement)<(*pMin))  *pMin = (*pElement);
      } // end for i
      *pSum   = asum;
      *pSumSq = a2sum;
    } // end if noDataStore

    //cout << redText << "In TStats calculation\n" << normText;    // debugging
    *pAvg   = (*pSum)/(TFloat)nData;
    // perform sigma calculation in long double to help avoid subtraction error
    // on squares of very large numbers
    *pSigma = (double)( (*pSumSq)/(long double)(nData-1) - 
                        (long double)(*pAvg)*(long double)(*pAvg)
                       *(long double)nData/(long double)(nData-1) );
    if(*pSigma<-0.00000005) {                                      // debugging
      warningMsg("sigma sum is negative? (variance = "+ftos(*pSigma)+" avg = "
                +ftos(*pAvg)+" with "+itos(nData)+" elements) on object "+desc);
      *pSigma = 0.0;  // manually force a valid value so that we can continue
    } // end if negative check
    *pSigma = sqrt(abs(*pSigma));  // debugging
    statsDirty = 0;  // all calculations are now complete - data is clean
    validData  = 1;  // if we made it to here, everything looks good
  } // end if isValid
  else { 
    if(!noDataStore) { 
      // sums are NOT already calculated so we must calculate them here
      TFloat* pElement = &data[0];
      // calculate sums
      asum  = (*pElement);                // initialize asum
      a2sum = (*pElement)*(*pElement);    // initialize a2sum
      // initialize min and max (note this is NOT done in constructor)
      if(pMax==NULL)  pMax = new TFloat; 
      if(pMin==NULL)  pMin = new TFloat;
      *pMax = (*pElement); 
      *pMin = (*pElement); 
      for(int i=1; i<nData; i++) { 
        pElement++;                       // assumes array[] data format
        asum  += (*pElement);  
        a2sum += (*pElement)*(*pElement);  
        // and some other stuff
        if((*pElement)>(*pMax))  *pMax = (*pElement);
        if((*pElement)<(*pMin))  *pMin = (*pElement);
      } // end for i
      *pSum   = asum;
      *pSumSq = a2sum;
    } // end if noDataStore

    //cout << redText << "In TStats calculation\n" << normText;  // debugging
    *pAvg   = (*pSum)/(TFloat)nData;
     statsDirty = 0;  // all calculations are now complete - data is clean
     validData  = 1;  // if we made it to here, everything looks good
    if(nData<=2) 
      cout << "Warning:  Number of data elements must be 3 or more for some "
           << "TStats1D object calculations for object " << desc << "\n";
    else         cerr << "Invalid calc function for TStats1D object.\n";
  } // end else
  return;  
}; // end calc


void TStats1D::store(TFloat val) {  
  // add another data element
  if(!noDataStore) {
    if(nData>=MaxNDataElements) {
      cerr << "Maximum number of TStats1D data elements reached." << endl;
      return;
    } // end if nData
  } // end if !noDataStore
  //cout << "Made it in the store function!\n";          // debugging
  //cout << "First element is " << data[0] << ".\n";     // debugging
  if(noDataStore) {
    *pSum   += val;
    *pSumSq += (long double)val*(long double)val;
    if(pMax==NULL)       { pMax = new TFloat;  *pMax = val; }
    else if(val>(*pMax))  *pMax = val;
    if(pMin==NULL)       { pMin = new TFloat;  *pMin = val; }
    else if(val<(*pMin))  *pMin = val;
    // have to modify for skew and kurtosis - by adding a3sum and a4sum???
  } // end if noDataStore
  else {
    data[nData] = val;  // taking advantage of 0 indexed array numbering
  } // end if else 
  nData++;
  statsDirty = 1;
  return;
};  // end store


void TStats1D::clear() { 
  // delete all data and zero *entire* data array for completeness 
  // (not just to size nData)
  //if(!noDataStore)  data.init(0.0);
  if(!noDataStore)  for(int i=0; i<nData; i++)  data[i] = 0.0;
  nData = 0;
  *pAvg  = 0.0;          *pSigma = 0.0;         *pSkew  = 0.0;  
  *pKurt = 0.0;          *pSum   = 0.0;         *pSumSq = 0.0;
  if(pMax!=NULL) { delete pMax;  pMax = NULL; }
  if(pMin!=NULL) { delete pMin;  pMin = NULL; }
  validData   = 0;  // invalid data indicator needed until data is added
  statsDirty  = 1;
  unordered   = 1;
}; // end clear


inline void TStats1D::resize(int newSize, TFloat fillVal) { 
  // resize the number of data elements up or down 
  // truncates data if newSize is less than current number of elements in data
  //if(newSize>0)  data.Resize(newSize, 1, fillVal);
  if(newSize>0) {
    if(data!=NULL)  delete[] data;
    data = new TFloat[newSize];
    for(int i=0; i<newSize; i++)  data[i] = fillVal;
  }
  else if(newSize==0) validData = 0;
  if(newSize<nData && newSize>0) { 
    nData = newSize;
    statsDirty = 1; 
  } // end if newSize
}; // end resize


void TStats1D::binTest(int nBins, TFloat xMin, TFloat xMax) { 
    // currently binTest is only valid for data between 0 and 1
    // variables needed to perform a bin test on the random numbers
    if(noDataStore) { 
      cout << "User specified no data storage, binTest is not available.\n";
      return;
    } // end if noDataStore
    
    TFloat binSize = fabs((xMax-xMin)/(TFloat)nBins);
    int i, wBin;
    //Array1D rBins(nBins,0.0);
    int *rBins;  rBins = new int[nBins];
    for(i=0; i<nBins; i++)  rBins[i] = 0;
    
    TFloat  rN = (TFloat)nData;

    cout << "Bin test info with a bin size of "  << binSize
         << " and a range of " << xMin << " to " << xMax 
         << " and number of data is " << nData << ":\n";
    for(i=0; i<nBins; i++)  cout << i << "\t" ;    cout << "\n";
    for(i=0; i<78; i+=4)    cout << "----";        cout << "\n";
    if(nData<1) { 
      warningMsg("No data elements?");  
      return;  // need a minimum ... 
    } // end if

    for(i=0; i<nData; i++) {  // make bin assignments
      wBin = (int)( (data[i]-xMin)/(xMax-xMin)*(TFloat)nBins );
      if(wBin==nBins)  wBin = nBins-1; // catch exact match case on end
      if(wBin<0 || wBin>nBins) {
        cerr << red << "There is a binTest problem!  wBin is " 
             << wBin << " but number of bins is " << nBins 
             << " and data[" << i << "] = " << data[i] << normText << endl;
        return;
      } // end if wBin
      rBins[wBin] += 1; // make bin increment
    } // end for i

    // output actual number distribution
    for(i=0; i<nBins; i++)  cout << rBins[i] << "\t" << flush;
    cout << "#'s out of " << itos(nData,1) << "\n" << flush;
    // output percentages
    string sBinPercent;
    for(i=0; i<nBins; i++) {
      sBinPercent = itos( (int)((TFloat)rBins[i]/rN*1000.0 + 0.5) );
      cout << sBinPercent.insert(sBinPercent.length()-1,".") << "\t" << flush;
    } // end for i
    cout << "%'s\n";
    // output number of sigmas for each bin
    for(i=0; i<nBins; i++) {
      cout << setprecision(3) << (xMin+binSize*(TFloat)i-(*pAvg))/(*pSigma)
           << "\t" << flush;
    } // end for i
    cout << "sigmas\n" << flush;

    delete[] rBins;
    return;
}; // end binTest


inline TFloat& TStats1D::operator[](int i) const {
  return data[i];
} // end standard array reference


inline bool TStats1D::isValid() {
#ifdef STATS_STRICT_ERROR_CHECKING_P
  if(StrictErrorChecking && validData) {
    // currently almost a dummy validation check
    if(nData>=0 && nData<=MaxNDataElements) return 1;  
    else  return 0;  // an invalid data set
  } // end if StrictErrorChecking
#else
  //return validData;
  return 1;  // debugging
#endif
};  // overloading isValid to make this a genuine class


void  TStats1D::shift(TFloat val) { 
  if(isValid() && nData>0) { 
    TFloat* pElement = &data[0];
    for(int i=0; i<nData; i++) { *pElement += val;  pElement++; }
    statsDirty = 1;
  } // end if isValid
  else cerr << "Invalid shift operation for TStats1D object.\n";
  return;  
}; // end shift


bool TStats1D::sort(int sortMethod) { 
  // return value is true if successful false if not successful
  int i, j;

  switch(sortMethod) {
    case 1: { // use selection sort
      // uses max and min selection sort the stats data in ascending order
      TFloat temp;  // swap variable
      TFloat *pMinN,  *pMaxN, *pCurN;  // pointer to respective numbers
      TFloat *pLastN, *pFirstN;  // pointer to current beginning and end of list

      pCurN = &data[0];
      pLastN = &data[nData-1];
      for(i=0; i<nData; i++) {
        pFirstN = pCurN;   // designate current first element for this loop
        pMinN   = pCurN;   // set current minimum as current first element in list
        pMaxN   = pLastN;  // set current maximum as current last element in list

        for(j=i+1; j<nData-i; j++) {
          pCurN++;         // increment current pointer to next number in list
          if ((*pCurN) < (*pMinN))       pMinN = pCurN;
          else if ((*pCurN) > (*pMaxN))  pMaxN = pCurN;
        } // end for j

        // now swap first value for minimum value
        temp = *pMinN;  *pFirstN = *pMinN;  *pMinN = temp;
        temp = *pMaxN;  *pLastN  = *pMaxN;  *pMaxN = temp;

        // now shrink part of list next loop looks at by shifting pointers
        pCurN++;   // increment current i index pointer
        pLastN--;  // decrement current last i index in the list
      } // end for i
      unordered = 0;  // now the original data order is permanently altered
      break;
    } // end if sortMethod==1

    case 2: { // use array[] based insertion sort
/*
      TFloat temp;  // temp variable
      for(i=1; i<nData; i++) {
        temp = data[i];
        j = i;
        while(j>0 && data[j-1]>temp) {
          //data[j]>data[j-1]; // there is a problem here
          j--;
        } // end while j
        data[j] = temp;
      } // end for i
      unordered = 0;  // now the original data order is permanently altered
*/
      break;
    } // end if sortMethod==2

    case 3: { // use pointer based insertion sort to eliminate some int mults
      TFloat  temp;      // temp variable
      TFloat *pJ, *pI;   // pointers to the i'th and j'th numbers
      TFloat *pPrevJ;    // pointer to previous j'th number

      pI = &data[0];   // dummy starting location is incremented in loop
      for(i=1; i<nData; i++) {
        pPrevJ = pI;   // we need j and j-1 pointer index locations stored
        pI++;          // increment to next number in list or 
                       //    start at i=1 in the beginning
        j = i;         // j and i counters start at same number index
        pJ = pI;       // j and i pointers point at same number to start
        temp = *pI;    // store next item in list to drop in correct place
        while(j>0 && (*pPrevJ)>temp) {
          *pJ = *pPrevJ;            // move data element up to next spot
          pJ  = pPrevJ;  pPrevJ--;  // decrement j'th and (j-1)'th pointers
          j--;
        } // end while j
        *pJ = temp;
      } // end for i
      unordered = 0;  // now the original data order is permanently altered
      break;
    } // end if sortMethod==2

    default: { return 0; } // invalid sort method - need to send error message
  } // end switch

  return 1;
}; // end sort


void TStats1D::output(int verbosity) { 
  // verbosity:  0 - outputs normal number stats and...
  //             1 - output more regular stats and...
  //             2 - additionally outputs binTest and...
  //             3 - additionally outputs all data
  if(isValid() && nData>2) { 
    if(statsDirty)  calc();
    // now output the stats to cout...
    if(verbosity>=1) { 
      cout << nData   << " elements with range of "
           << (*pMin) <<  " to " << (*pMax) << ", ";
    } // end verbosity 1
    // output standard stats
    cout << "avg = " << (*pAvg) << " and sigma = " << (*pSigma) << "\n";
    // need to add in skew and kurtosis values when they are implemented
    if(verbosity>=2 && !noDataStore)  binTest(10,*pMin,*pMax);
    if(verbosity>=4 && !noDataStore)  { // output entire data set
      cout << "List of data";
      if(unordered) cout << " in original order";
      cout << ":";
      for(int i=0; i<nData; i++)  cout << " " << data[i];
      cout << endl;
    } // end verbosity 4
  } // end if isValid
  else { 
    if(nData<=2) cerr << "Error:  Number of data elements must be 3 or more "
                      << "for TStats1D object output for object "
                      << desc << ".\n";
    else         cerr << "Invalid output function for TStats1D object.\n";
  } // end else
  return;
}; // end sort


// Constructors and destructor ---------------------------------------------
TStats1D::TStats1D(string description) {  
  // initialize an empty stats variable with maximum possible # of elements
  // (obviously can be wastefull of memory)
  int MaxElements = 0;  bool noStore = 1;
  nData = 0;
  if(MaxElements==0)  noStore = 1;      // override parameter if size is zero
  if(noStore==1)      MaxElements = 0;  // override parameter if size is zero
  if(MaxElements>0 && !noStore) {
    data = new TFloat[MaxElements];
    for(int i=0; i<MaxElements; i++)  data[i] = 0.0;  // initialize to zero
  } // end if MaxElements
  else data = NULL;

  //cout << "Made it in the constructor function with "  // debugging
  //     << MaxElements << " data storage capacity!\n";  // debugging
  //cout << "First element is " << data[0] << ".\n";     // debugging
  
  // initialize invalid values of stats data - original idea was to use 
  // a pointer to provide a genuine invalid result since skew, average, 
  // max, and min can take on any value, positive or negative
  // (currently only implemented for max and min however)
  //pAvg  = NULL;         pSigma = NULL;         pSkew  = NULL;  
  //pKurt = NULL;         pSum   = NULL;         pSumSq = NULL;
  pAvg   = new TFloat;    pSigma = new TFloat;   pSkew  = new TFloat;  
  pKurt  = new TFloat;    pSum   = new TFloat;   pSumSq = new long double;
  *pAvg  = 0.0;          *pSigma = 0.0;         *pSkew  = 0.0;  
  *pKurt = 0.0;          *pSum   = 0.0;         *pSumSq = 0.0;
  // max and min are NOT initilized since the range of data is unknown at start
  pMax = NULL;            pMin = NULL;
  validData   = 0;  // invalid data indicator needed until data is added
  statsDirty  = 1;
  unordered   = 1;
  noDataStore = noStore;  // do we store all data or just running sums?
  desc = description;
};
TStats1D::TStats1D(int MaxElements, bool noStore, string description) {  
  // initialize an empty stats variable with maximum possible # of elements
  // (obviously can be wastefull of memory)
  nData = 0;
  if(MaxElements==0)  noStore = 1;      // override parameter if size is zero
  if(noStore==1)      MaxElements = 0;  // override parameter if size is zero
  if(MaxElements>0 && !noStore) {
    data = new TFloat[MaxElements];
    for(int i=0; i<MaxElements; i++)  data[i] = 0.0;  // initialize to zero
  } // end if MaxElements
  else data = NULL;

  //cout << "Made it in the constructor function with "  // debugging
  //     << MaxElements << " data storage capacity!\n";  // debugging
  //cout << "First element is " << data[0] << ".\n";     // debugging
  
  // initialize invalid values of stats data - original idea was to use 
  // a pointer to provide a genuine invalid result since skew, average, 
  // max, and min can take on any value, positive or negative
  // (currently only implemented for max and min however)
  //pAvg  = NULL;         pSigma = NULL;         pSkew  = NULL;  
  //pKurt = NULL;         pSum   = NULL;         pSumSq = NULL;
  pAvg   = new TFloat;    pSigma = new TFloat;   pSkew  = new TFloat;  
  pKurt  = new TFloat;    pSum   = new TFloat;   pSumSq = new long double;
  *pAvg  = 0.0;          *pSigma = 0.0;         *pSkew  = 0.0;  
  *pKurt = 0.0;          *pSum   = 0.0;         *pSumSq = 0.0;
  // max and min are NOT initilized since the range of data is unknown at start
  pMax = NULL;            pMin = NULL;
  validData   = 0;  // invalid data indicator needed until data is added
  statsDirty  = 1;
  unordered   = 1;
  noDataStore = noStore;  // do we store all data or just running sums?
  desc        = description;
};  // end TStats1D default constructor


TStats1D::TStats1D(TStats1D& a) {  
  if(a.isValid()) { 
    nData = a.nData;
    if(data!=NULL)  delete[] data;
    if(nData>0) {  
      data = new TFloat[nData]; // resize data
      for(int i=0; i<nData; i++)  data[i] = a.data[i];
    }
    else data = NULL;
  
    pAvg   = new TFloat;    pSigma = new TFloat;     pSkew  = new TFloat;  
    pKurt  = new TFloat;    pSum   = new TFloat;     pSumSq = new long double;
    *pAvg  = *(a.pAvg);    *pSigma = *(a.pSigma);   *pSkew  = *(a.pSkew);  
    *pKurt = *(a.pKurt);   *pSum   = *(a.pSum);     *pSumSq = *(a.pSumSq);
    // now handle special cases of max and min with possible null pointers
    if(pMax!=NULL)   { delete pMax;          pMax = NULL; }
    if(pMin!=NULL)   { delete pMin;          pMin = NULL; }
    if(a.pMax!=NULL) { pMax = new TFloat;   *pMax  = *(a.pMax); } 
    if(a.pMin!=NULL) { pMin = new TFloat;   *pMin  = *(a.pMin); }
    validData   = a.validData;
    statsDirty  = a.statsDirty;
    unordered   = a.unordered;
    noDataStore = a.noDataStore;
  } // end if a.isValid
};  // end TStats1D copy constructor


TStats1D::~TStats1D() {  
  nData = 0;
  // initialize invalid values of stats data
  delete pAvg;    delete pSigma;   delete pSkew;
  delete pKurt;   delete pSum;     delete pSumSq;
  // now handle special cases of max and min with possible null pointers
  // the check is supposed to be redundant with C++ standards
  if(pMax!=NULL)  delete pMax;
  if(pMin!=NULL)  delete pMin;
  if(data!=NULL)  delete[] data;
};  // end TStats1D constructor

} // end namespace MSL

