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
#include <iostream>
#include <iomanip>
#include <cmath>
#ifndef ML_UTILS_H
#include "ML_Utils.h"
#endif
using namespace std;

// MersenneTwister random numbers class
//#include "/home/ronhovde/RandomLib/Random.hpp"
//#include "/users/peters/Desktop/Cluster/RandomLib/Random.hpp"
#ifdef LINUX_OS
#include <values.h>          // unix version of limits.h
  // Now, redefine *nix MAX* values that are used in the program code back to 
  // the dos names from limits.h.  This helps to avoid a lot of preprocessor
  // #ifdef DOS_OS... statements all throughout the code.
  //const int INT_MAX = MAXINT;  
#else
#include <limits.h>
#endif

// compiler won't let me use a define to rule out multiple definitions!!!???
//#ifndef double_T
// double allows easily changing between double and float data types
//typedef double TFloat; // float type definition
//#define double_T
//#endif


string toLowerCase(const string s) {
  char c;
  string t = s;
// My own lower case conversion routine since I am not interested in
// incorporating complicated generic string algorithms.
// Other characters are ignored.
  for(unsigned int i=0; i<t.length(); i++) { // unsigned eliminates a warning
    c = t[i];
    if(c>='A' && c<='Z') {
      c = (char)((int)c+32);
      t[i] = c;
    } // end if
  }
  return t;
} // end toLowerCase


string toUpperCase(const string s) {
  char c;
  string t = s;
// My own lower case conversion routine since I am not interested in
// incorporating complicated generic string algorithms.
// Other characters are ignored.
  for(unsigned int i=0; i<t.length(); i++) { // unsigned eliminates a warning
    c = t[i];
    if(c>='a' && c<='z') {
      c = (char)((int)c-32);
      t[i] = c;
    } // end if
  }
  return t;
} // end toUpperCase


int charSum(const string s, int convertCase) {
  int cSum = 0;
  string t=s;
// My own lower character sum routine.
// convertCase:  0 : no conversion performed on input string
// convertCase:  1 : converts alphabetical letters to lower case
// convertCase:  2 : converts alphabetical letters to upper case
// convertCase:  other : ignored
  if(convertCase==1)      t = toLowerCase(s);
  else if(convertCase==2) t = toUpperCase(s);
  for(unsigned int i=0; i<t.length(); i++) cSum += (int)s[i]; // unsigned eliminates a warning
  return cSum;
} // end charSum


string itos(int I, bool UseThSep, const char ThSep,
                   bool usePad, unsigned int padDigits, char padChar) {
/* For some reason (?!!), Borland doesn't supply an 'integer2string' coversion
method and defaulting to the c style itoa is clumsy.  This is a workaround
function which hopefully will be eliminated when I change compilers.
This function is only for decimal format.
Passed parameters:
  I - the integer
  UseThSep - use a (comma) separator for the thousands (and higher) places
  ThSep - the actual character to use as the thousands separator 
  pad is a boolean for whether the string is padded on the left with padChar
    up to a total of padDigits
*/
  //cout << "Entering itos\n";  // debugging
  string s = "", sISign = "", sTemp = " ", scTemp = "";
  string sThSep = "";  sThSep += ThSep;
  //if(I==0) { s = "";               }  // zero case
  if(I<0)  { sISign = "-";  I = -I; }  // store negative sign for negative integers
  int ILength;
  if(I==0)  ILength = 1;
  else      ILength = (int)floor(log10((double)I)) + 1;
  //cout << "After itos initializations\n";  // debugging
  for(int i=0; i<ILength; i++) {
    // This uses a quirky because (char) to (string) conversion since
    // a proper cast operation is not provided by Borland.
#ifdef DOS_OS
    sTemp.insert(i+1,(char)(I%10 + 48));
    //cout << "After itos char insert operation\n";  // debugging
#else
    scTemp += (char)(I%10 + 48);
    //cout << "After itos char to string cast\n";  // debugging
    sTemp.insert(i+1,scTemp);
    //cout << "After itos string insert operation\n";  // debugging
#endif
    //sTemp[i+1] = (char)(I%10 + 48);
    I = I/10;
    //cout << "After itos integer division\n";  // debugging
    scTemp = "";
    //cout << "After itos string reset\n";  // debugging
  }
  //cout << "After itos stringification\n";  // debugging
  for(unsigned int i=ILength; i>0; i--) { s += sTemp[i]; } // unsigned eliminates a warning
  //cout << "After itos copyfying\n";  // debugging
  //s += (string)(itoa(abs(I),charTemp,10));  // gcc doesn't recognize itoa()?

  // Add padChar characters on the left - some odd looking numbers can result
  //cout << s << " and its length is " << s.length() << endl; // debugging
  if(usePad && padDigits>s.length() && padDigits<=10) {
    string padString("");
    // unsigned eliminates a warning
    for(unsigned int i=0; i<(padDigits-s.length()); i++) padString += padChar;
    //cout << "padString is " << padString << endl;           // debugging
    s.insert(0,padString);
  } // end inclusion of thousands separator
  // Add appropriate thousands separators
  if(UseThSep) {
    if(s.length()>3)  s.insert(s.length()-3,sThSep);
    if(s.length()>7)  s.insert(s.length()-7,sThSep);
    if(s.length()>11) s.insert(s.length()-11,sThSep);
  } // end inclusion of thousands separator
  //cout << "After itos commification\n";  // debugging
  return (sISign + s);
} // end itos


string dtos(double f, short unsigned int NDigits, bool UseThSep, 
            bool UseSciNotation, bool UseMathSciNotation, bool UseDecSep, 
            short unsigned int DecSepInterval, double UpperSciNotationStart, 
            double LowerSciNotationStart, bool UseEuroNumbers) {
  // float version maxes out at 7 significant digits
  // NDigits is the number of decimal digits at the moment
  // this is a crude version to get things working at the moment and move on
  string s;
  int ePower;
  if(f<0.0) { s = "-";  f = -f; }  // set minus sign and take absolute value
  else        s = "";

  // override use of scientific notation for very large numbers
  if(f>2.0e9)  UseSciNotation = 1;  
  
  if(UseSciNotation) {
    if(log10(f)<0.0)  ePower = (int)log10(f)-1; // correct for negative powers
    else              ePower = (int)log10(f);
    f /= pow(10.0,ePower);
  } // end if UseSciNotation

  int wInt, fInt;
  // work with greater than zero part, temporarily drop fractional part
  // there is a catch if rounding changes the ones place, so we manually add
  // the rounding value to f first
  f += 0.5*pow(0.1,NDigits);
  if(f>1.0)  s += itos((int)f,UseThSep) + ".";
  else       s += "0.";

  // now fill the decimal part
  if(f>=1.0)  f -= (double)( (int)f );
  s += itos((int)( f*pow(10.0,NDigits) ),0,' ',1,NDigits);

  //cout << s[s.length()-1] << endl;
  // eliminate padded zeroes on the end
  int i = s.length()-1;
  while(i>1 && s[i]=='0') i--;  // find last non-zero character
  //if(!bPad)
  if(i<s.length()-1)  s.erase(i+1);

  if(UseSciNotation && ePower!=0)  s += "e" + itos(ePower);

  return s;
} // end dtos



ostream& OutputLineSeparator(ostream& fo, int length,
                  char separatorChar, char commentChar) {
  fo << commentChar;
  for (int i=0; i<length; i++) { fo << separatorChar; }  fo << endl;
  return fo;
} // end line separator


// Sort functions to-do's:  
// 1. At some point, I need to allow a parameter to sort in descending order.
template <typename T>
void selectSort(T data[], int Length) {
  // Uses an optimized selection (shaker) sort to sort the data in ascending
  // order working on both the maximum and minimum values on each loop pass.
  T   temp;                   // swap variable
  T  *pMinN, *pMaxN, *pCurN;  // pointer to respective numbers
  T  *pFirstN, *pLastN;       // pointer to current beginning and end of list

  pCurN = data;
  pLastN = &data[Length-1];
  for(int i=0; i<Length; i++) {
    pFirstN = pCurN;   // designate current first element for this loop
    pMinN   = pCurN;   // set current minimum as current first element in list
    pMaxN   = pLastN;  // set current maximum as current last element in list

    for(int j=i+1; j<Length-i; j++) {
      pCurN++;         // increment current pointer to next number in list
      if      ((*pCurN) < (*pMinN))  pMinN = pCurN;
      else if ((*pCurN) > (*pMaxN))  pMaxN = pCurN;
    } // end for j

    // now swap for min and max values to beginning and end of data
    temp = *pMinN;  *pFirstN = *pMinN;  *pMinN = temp;
    temp = *pMaxN;  *pLastN  = *pMaxN;  *pMaxN = temp;

    // now shrink the part of the list that is used for the next loop
    pCurN++;           // increment current i index pointer
    pLastN--;          // decrement current last i index in the list
  } // end for i

  return;
} // end selectSort

template <typename T>
void insertionSort(T data[], int Length) {
  // Uses an optimized insertion sort to sort a generic integer data array.
  int i, j; 
  T   temp;              // temporary value
  T  *pJ, *pI, *pPrevJ;  // pointers to i'th, j'th, and previous j'th numbers

  pI = data;             // dummy starting location is incremented in loop
  for(i=1; i<Length; i++) {
    pPrevJ = pI;         // we need j and j-1 pointer index locations stored
    pI++;                // increment to next number in list or start at i=1
    j = i;               // j and i counters start at same number index
    pJ = pI;             // j and i pointers point at same number to start
    temp = *pI;          // store next item in list to drop in correct place
    while(j>0 && (*pPrevJ)>temp) {
      *pJ = *pPrevJ;            // move data element up to next spot
      pJ  = pPrevJ;  pPrevJ--;  // decrement j'th and (j-1)'th pointers
      j--;
    } // end while j
    *pJ = temp;
  } // end for i

  return;
} // end insertionSort


bool  floatEq(double f1, double f2, double epsilon) {
  // return a boolean value indicating that f1 is within epsilon of f2 
  return ( (f1<f2+epsilon) && (f1>f2-epsilon) );
} // end floatEq

