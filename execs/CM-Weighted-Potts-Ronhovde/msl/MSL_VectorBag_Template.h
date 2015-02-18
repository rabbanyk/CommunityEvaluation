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
MSL - My (or Math) Scientific Library - TVectorBag<T>
See MSL_Array.h for general comments and some documentation.

TVectorBag<T> Specific Notes:
These notes handle only behavior specific to the TVectorBag<T> implementation.

1. If StrictErrorChecking != 1 then if the objects are of different sizes,
the assignment operator destroys the original object and recreates one with the
new 'correct' size to make the assignment.

2. The 'Assign' function is a (useful) function that equates the elements of the
TVectorBag<T> from Initial up to but not including the element arrayName[Length] (i.e.
it fills a total of (Length-Initial) elements of the array). The changed array
is returned by reference as well as having been passed by reference, so it can
be used either way.

3. The 'dot' product and 'mag' for TVectorBag<T> here is a convenience that just
treats the object as a single long vector performing the operations accordingly.
The motivation to include this, even though it doesn't match our idea of the
behavior of a "Matrix," is from (iterative) solutions of PDE's.  There, it can
be easier and more natural to treat the solution function as a 2D grid; but the
theory essentially treats the unknown vector as a 1D vector where at times dot
products and vector magnitudes are calculated.  Thus, having a method to treat
the 2D array as a 1D vector then does fit more naturally than it appears.

4. Operator[] is defined only for the T* return value since once we have a
T* (usually a double*), it is a intrinsic type.  Operator[] must be defined
as a member function which we cannot do for double* type.  Therefore, we are
forced to let the compiler handle out of bounds and error checking after this
point.  However, full error checking is done at the pointer return stage;
therefore, from there at least, it is assured that we have valid objects after
this point.
If needed, an alternate implementation at some point would be to create an
array of Array1D objects.  At this point though, I do not need any of the 
added functionality of such an implementation.  It might prove useful if I 
implemented a full TMatrix class.
*/
/* Known problems:
1. For some reason operator<< is trying to call the pure virtual isValid() function;
therefore, the error check is commented out for now.
*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TVECTORBAG_H
#define MSL_TVECTORBAG_H

#include <string>
#ifndef MSL_ARRAY_H
#include "MSL_Array.h"
#endif

namespace MSL {

template <typename T> class  TVectorBag;
template <typename T> TVectorBag<T>  operator*(T x, const TVectorBag<T>& a);
//template <typename T> istream&     operator>>(istream& fin, TVectorBag<T>& a);

#ifndef MSL_TVECTORBAG_T
// define some common vector types
typedef TVectorBag<TFloat> TVecBagFloat;// TFloat  array type definition
typedef TVectorBag<int>    TVecBagInt;  // integer array type definition
//typedef TVectorBag<bool> TSparseBool; // integer array type definition
#define MSL_TVECTORBAG_T
#endif


//------------------------------------------------------------------------------
template <typename T> 
class TVectorBag : public TMatrix {
// Memory allocation is dynamic, and this TVectorBag<T> uses the standard C++ array
// indexing with brackets as in arrayName[col][row] with zero indexing.
 //friend TVectorBag<T>  operator*<T>(T x, const TVectorBag<T>& a);
 //friend istream&    operator>><T>(istream& fin, TVectorBag<T>& a);
 public:
  // overloaded operators
         TVectorBag<T>&    operator=(const TVectorBag<T>& b);
         TVectorBag<T>&    operator*=(T x);          // *= by scalar
         TVectorBag<T>&    operator/=(T x);          // /= by scalar
         //TVectorBag<T>   operator-();              // unary '-'
         //TVectorBag<T>   operator*(T x);           // * by scalar
         //TVectorBag<T>   operator/(T x);           // / by scalar
  inline Array1D<T>&       operator[](int i);
         T&                operator()(int i, int j);
  // Public member functions
  inline int               getCols() const { return Cols; };
  inline T                 getMax(bool bAbs=0) const;
  inline T                 getMin(bool bAbs=0) const;
         T                 getMax(int &col, int &row, bool bAbs=0) const;
         T                 getMin(int &col, int &row, bool bAbs=0) const;
         int               getNElements() const;
  inline bool              isValid() const;               // valid array check
  inline bool              isValid() ;   // debugging     // valid array check
  // The following functions treat the 2D array as a single large 1D vector.
         void              resize(int Cols); // for now, it DELETES ALL DATA!
         void              init(T Value);    // Reinitialize with index changes
  // Standard public member functions
            TVectorBag(int Cols = 0);    // bare constructor
            //TVectorBag(int Cols, int Rows, T Value); // initializing constructor
            TVectorBag(const TVectorBag<T> &b);          // copy constructor
           ~TVectorBag();                              // destructor
 protected:
  // Data elements
  int Cols;//, Rows;         // number of columns and rows in the array
  Array1D<T> *data;   // actual data array as an array of TVector
}; // end class TVectorBag<T>;

//----------------------TVectorBag<T> Member Function Declarations--------------------

//Assign default static parameters
//TVectorBag<T>::StrictErrorChecking = 1;

template <typename T> 
TVectorBag<T>::TVectorBag(int Cs) {
  // Raw constructor th at skips the initialization
  // Of course, without the initialization, it takes the chance that the user
  // will access junk, but it is provided for "more efficient" applications
  // where the user wants to initialize himself.
  // As the default constructor, it generates a TVectorBag<T> object that is a zero
  // length real TVectorBag<T> object but officially an invalid one since it holds no
  // data.  This is useful for file io and assignment operations.
  Cols = Cs;
  if(Cols>0) {
    data = new Array1D<T>[Cols];   // Set up columns
    //for (int i=0; i<Cols; i++) { 
    //  data[i] = new Array1D<T>; 
      //if(Rows>0)  data[i].resize(Rows);
    //} // end for i
    ValidArray = Valid2DArrayFlag;
  } // end un-initialized array creation
  else if(Cols==0) {
    ValidArray = 0;  data = NULL;
    if(VerboseErrorReports!=0)
      cerr << "Attempting to create a zero length TVectorBag<T> object.  "
           << "Cols = " << Cols << endl;
  } // end zero length if
  else { // Cols<0 or Rows<0 case
    ValidArray = -1;  data = NULL;
    cerr << "Attempting to create a negative length TVectorBag<T> object.  "
         << "Cols = " << Cols << endl;
  } // end negative length if
}; // end TVectorBag<T> constructor
/*
// NEEDS SOME WORK!!!!
template <typename T> 
TVectorBag<T>::TVectorBag(int Cs, int Rs, Array1D<int> &v, T Value) { 
  // constructor with init vector and values
  Cols = Cs;  Rows = Rs;
  if(Cols>0 && Rows>0) {
    data = new Array1D<T>*[Cols];   // Set up columns
    for (int i=0; i<Cols; i++) { 
      data[i] = new Array1D<T>; 
      if(Rows>0)  data[i].resize(v[i]);
    } // end for i
    // Now, initialize array with 
    for (int j=0; j<Cols; j++)  for (int i=0; i<data[i].getSize(); i++)
      data[j][i]  = Value;
    ValidArray = Valid2DArrayFlag;
  } // end un-initialized array creation
  else { 
    cerr << "Error creating TVectorBag<T> object with initialization!  "
         << "Rows x Cols = " << Rows << " x " << Cols << endl;
    ValidArray = -1; data = NULL;  Cols = 0;  Rows = 0; 
  } // creation error 
}; // end TVectorBag<T> constructor
*/
template <typename T> 
TVectorBag<T>::TVectorBag(const TVectorBag<T> &b) { // copy constructor
  if(b.isValid()) { // check for valid array
    Cols = b.Cols;
    data = new Array1D<T>[Cols];   // Set up columns
    for (int i=0; i<Cols; i++) { 
      //data[i] = new Array1D<T>; 
      if(b[i].getSize()>0)  data[i].resize(b[i].getSize());
    } // end for i
    // Now, initialize array
    for(int i=0; i<Cols; i++)  for(int j=0; j<data[i].getSize(); j++) 
      data[i][j] = b.data[i][j]; // end ij
    ValidArray = Valid2DArrayFlag;
  } // end un-initialized array creation
  else if(Cols==0) { ValidArray = 0;  data = NULL; }
  else if(Cols<0)  { ValidArray = -1; data = NULL; }
}; // end TVectorBag<T> copy constructor

template <typename T> 
TVectorBag<T>::~TVectorBag() {
  cout << "Entering TVectorBag<T> destructor " << Cols << "... " << flush;
  if(isValid()) {
    //for (int i=0; i<Cols; i++)  delete data[i];  // Delete columns
    delete[] data;  data  = NULL;
  } // end valid array destruction
  else  cerr << "Error Destructing TVectorBag<T>.\n";
  Cols = 0;  ValidArray = 0;
  cout << "exiting TVectorBag<T> destructor" << endl;
}; // end TVectorBag<T> destructor

template <typename T> 
inline bool TVectorBag<T>::isValid()       { return 1; }; 
template <typename T> 
inline bool TVectorBag<T>::isValid() const { return 1; }; 
/*inline bool TVectorBag<T>::isValid() const { 
  if(StrictErrorChecking==1) {
    if((ValidArray==Valid2DArrayFlag) && Cols>0 && Rows>0) return 1;
    else { return 0;  cerr << "Invalid TVectorBag<T> error\n"; }
  } // if StrictErrorChecking
  else return 1;
}; // end TVectorBag<T>::isValid
*/
template <typename T> 
TVectorBag<T>& TVectorBag<T>::operator=(const TVectorBag<T>& b) { // assignment for 2 arrays
  if(b.isValid()) { // check for valid array
    //for (int i=0; i<Cols; i++)  delete data[i];  // Delete columns
    delete[] data;
    Cols  = b.Cols;
    fixed = b.fixed;
    data  = new Array1D<T>[Cols]; // Set up columns
    for (int i=0; i<Cols; i++) { 
      //data[i]  = new Array1D<T>; 
      data[i].resize(b[i].getSize());
    } // end for i
    // Now, initialize array
    for(int i=0; i<Cols; i++)  for(int j=0; j<data[i].getSize(); j++)
      data[i][j] = b.data[i][j];
    ValidArray = Valid2DArrayFlag;
  } // end un-initialized array creation
  else if(Cols==0) { ValidArray = 0;  data = NULL; }
  else if(Cols<0)  { ValidArray = -1; data = NULL; }
  return *this;
}; // end =
/*
template <typename T> 
TVectorBag<T>& TVectorBag<T>::operator+=(const TVectorBag<T>& b) {   // additive assignment
  // check for valid arrays first
  if(isValid() && b.isValid() && Cols==b.Cols && Rows==b.Rows) {
    for (int j=0; j<Cols; j++) { for (int i=0; i<Rows; i++) {
      data[j][i] += b.data[j][i]; }} // end ij
  } // end if
  return *this;
}; // end +=

template <typename T> 
TVectorBag<T>& TVectorBag<T>::operator-=(const TVectorBag<T>& b) {   // subtraction assignment
  // check for valid arrays first
  if(isValid() && b.isValid() && Cols==b.Cols && Rows==b.Rows) {
    for (int j=0; j<Cols; j++) { for (int i=0; i<Rows; i++) {
      data[j][i] -= b.data[j][i]; }} // end ij
  } // end if
  return *this;
}; // end -=
*/
template <typename T> 
TVectorBag<T>& TVectorBag<T>::operator*=(T x) {
  if(isValid()) {
    for (int i=0; i<Cols; i++)  for(int j=0; j<data[i].getSize(); j++)
      data[i][j]  *= x; 
  } // end if error checking
  return *this;
}; // end *=

template <typename T> 
TVectorBag<T>& TVectorBag<T>::operator/=(T x) {   // T division assignment
  if(isValid())  { // check for valid array first
    if(x!=0) {
      for (int j=0; j<Cols; j++)  for(int i=0; i<data[i].getSize(); i++) 
        data[j][i]  /= x;
    } // end if
    else cerr << "Division by zero in TVectorBag<T> scale operation!\n";
  } // end if error checking
  return *this;
}; // end /=



template <typename T> 
T& TVectorBag<T>::operator()(int i, int j) {
  if(i<0 || j<0)  errorMsg("Warning array index < 0 in op() for VectorBag!");
  else { if(j<Cols)  if(i<data[j].getSize())  return (data[i])[j]; }
  return 0;  // default case for sparse matrix value that was not found
}; // end array reference


template <typename T> 
inline Array1D<T>&  TVectorBag<T>::operator[](int i) {
  return data[i];  // default case for sparse matrix value not found
}; // end array reference


// Other useful array functions
template <typename T> 
void TVectorBag<T>::resize(int Cs) {
  // Crude resize function to "get things working" at the moment
  // At the moment it DELETES ALL DATA!!!
  // first delete existing array data
  if(Cs>0) {
    //for(int i=0; i<Cols; i++)  delete data[i];
    delete[] data;
    // now recreate in the new empty image (only empty column vectors declared)
    data = new Array1D<T>[Cs];
    //for(int i=0; i<Cols; i++)  data[i] = new Array1D<T>;
    Cols = Cs;
  } // end if Cs, Rs
  else  cerr << "TVectorBag resize() error!  Aborting." << endl;
  return;
};

template <typename T> 
void TVectorBag<T>::init(T Value) {
  // reinitialize the Array but not the indexes
  if(isValid()) {
    for(int i=0; i<Cols; i++)  
      for(int j=0; j<data[i].getSize(); j++)  
        (data[i])[j] = Value; 
  } // end if
  return;
};

template <typename T> 
int  TVectorBag<T>::getNElements() const {
  // returns the total number of elements in the Bag
  int rowSum = 0;
  //if(isValid())  
    for(int i=0; i<Cols; i++)  rowSum += data[i].getSize();
  return rowSum;
};

template <typename T> 
inline T  TVectorBag<T>::getMin(bool bAbs) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  int iDummy, jDummy;
  return getMin(iDummy,jDummy,bAbs);
}; // end getMin
template <typename T> 
T  TVectorBag<T>::getMin(int &col, int &row, bool bAbs) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  T minVal;
  if(bAbs) { minVal = abs(data[0][0] ); }
  else     { minVal = data[0][0] ;      }
  if(isValid()) {
    for(int i=0; i<Cols; i++) { for(int j=0; j<data[i].getSize(); j++) {
      if(bAbs) { 
        if(abs(data[j][i])<minVal) {
          minVal = abs(data[j][i] );  
          col = i;  row = j; 
        } // end if
      } else { 
        if(data[j][i]<minVal) {
          minVal = data[j][i] ;
          col = i;  row = j; 
        } // end if
      } // end else
    }} // end ij
  } // end if isValid
  return minVal;
}; // end getMin with row,col


template <typename T> 
inline T  TVectorBag<T>::getMax(bool bAbs) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  int iDummy, jDummy;
  return getMax(iDummy,jDummy,bAbs);
}; // end getMax
template <typename T> 
T  TVectorBag<T>::getMax(int &col, int &row, bool bAbs) const {
  // return the maximum value (or absolute value if bAbs==1) in the matrix
  T maxVal;
  if(bAbs) { maxVal = abs(data[0][0] ); }
  else     { maxVal = data[0][0] ;      }
  if(isValid()) {
    for(int i=0; i<Cols; i++) { for(int j=0; j<data[i].getSize(); j++) {
      if(bAbs) { 
        if(abs(data[j][i])>maxVal) {
          maxVal = abs(data[j][i] );
          col = i;  row = j;
        } // end if
      } else { 
        if(data[j][i]>maxVal) {
          maxVal = data[j][i] ;
          col = i;  row = j; 
        } // end if
      } // end else
    }} // end ij
  } // end if isValid
  return maxVal;
}; // end getMax with row,col


/*
template <typename T> 
TVectorBag<T>& TVectorBag<T>::Assign(TVectorBag<T>& a, TVectorBag<T>& b,
         int CIFinal, int RIFinal, int CInitial = 0, int RInitial = 0) {
// Assignment function that will equate the elements of the Array1D from Initial
// up to but not including the element arrayName[Length] (i.e. it fills a total
// of (IFinal-Initial) elements of the array). The changed array is returned by
// reference as well as having been passed by reference, so it can be used
// either way.
  if(a.isValid() && b.isValid() &&
     CIFinal<=b.Cols && RIFinal<=a.Rows && CInitial>=0 && RInitial>=0)
    for(int i=RInitial; i<RIFinal; i++) { for(int j=CInitial; j<CIFinal; j++) {
      a.data[j][i] = b.data[j][i]; }} // end ij
  return a;
}; // end partial assign

template <typename T> 
TVectorBag<T> TVectorBag<T>::operator+(const TVectorBag<T>& a) { // array addition
  TVectorBag<T> c(a.Cols,a.Rows);
  if(isValid()&&a.isValid() && Rows==a.Rows && Cols==a.Cols)
      for(int i=0; i<Rows; i++) { for(int j=0; j<Cols; j++) {
        c.data[j][i] = data[j][i] + a.data[j][i]; }} // end ij
  return c;
}; // end +

template <typename T> 
TVectorBag<T> TVectorBag<T>::operator-(const TVectorBag<T>& a) { // array subtraction
  TVectorBag<T> c(a.Cols,a.Rows);
  if(isValid()&&a.isValid() && Rows==a.Rows && Cols==a.Cols)
      for(int i=0; i<Rows; i++) { for(int j=0; j<Cols; j++) {
        c.data[j][i] = data[j][i] - a.data[j][i]; }} // end ij
  return c;
}; // end binary -

template <typename T> 
TVectorBag<T>  TVectorBag<T>::operator-() {             // unary negative
  TVectorBag<T> c(Cols,Rows);
  if(isValid()) { // check for valid array first
      for (int j=0; j<Cols; j++) { for(int i=0; i<Rows; i++) {
        c.data[j][i] = -data[j][i]; }} // end ij
  } // end array error checking
  return c;
}; // end unary -

template <typename T> 
TVectorBag<T> TVectorBag<T>::operator*(T x) {   // vector scaling vec*x
  TVectorBag<T> c(Cols,Rows);
  if(isValid())
    for(int i=0; i<Rows; i++) { for(int j=0; j<Cols; j++) {
      c.data[j][i] = data[j][i]*x; }} // end ij
  return c;
}; // end scalar multiplication

template <typename T> 
TVectorBag<T> TVectorBag<T>::operator/(T x) {   // scalar division
  TVectorBag<T> c(Cols,Rows);
  if(isValid() && x!=0.0)
    for(int i=0; i<Rows; i++) { for(int j=0; j<Cols; j++) {
      c.data[j][i] = data[j][i]/x; }} // end ij
  else if(x==0.0) cerr << "Scaling division by zero error on Array1d object.\n";
  return c;
}; // end scalar division
// x/TVectorBag<T> is not a valid operation, so it is not defined
*/
//------------------------------------------------------------------------------
// End TVectorBag<T> method definitions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// define related non-member functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Overloaded iostream operators
template <typename T> 
ostream& operator<<(ostream& fout, TVectorBag<T>& a) {
  // Outputs the Array1D in MatrixMarket format.  Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  //if(a.isValid()) { 
    fout << "NColumns:  " << a.getCols() << endl;
    for(int i=0; i<a.getCols(); i++)  
      fout << "Node " << i <<  ": \t" << a[i] << "\n"; 
  //} // end if
  return fout;
};
/*
template <typename T> 
istream& operator>>(istream& fin, TVectorBag<T>& a) {
  int M, N;  // cols and rows - here M = 1 since it is an Array1D
  int nComments=0;
  string sInput;
  bool inputErrorFlag = 0;
  // Inputs the Array1D from MatrixMarket format. Since the format specification
  // includes comment lines, the user  allowed to comment as needed, so this
  // only outputs the *data* and not the format header or any comments.
  if(a.isValid()) {
    // Skip all comment lines (which should all be at the top of the output file
    // All comment lines should begin with a '%'.  Comments cannot be on an
    // existing data line even if preceded by a '%' after the data.
    // 1024 is technically the maximum length in MM format, but we ignore a line
    // much larger than that just in case.

    // First, check first format line.  It must conform to MatrixMarket
    // specification by.  We only import the type.
    // "%%MatrixMarket matrix array real general"
    fin >> sInput; if(sInput!="%%MatrixMarket") { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="matrix")         { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="array")          { inputErrorFlag=1; }
    fin >> sInput; if(sInput!="real")           { inputErrorFlag=1; }
    std::getline(fin,sInput,'\n'); if(sInput!="general") { inputErrorFlag=1; }
    if(inputErrorFlag!=0) {
      cerr << "TVectorBag<T> MatrixMarket header incorrect.\n"; return fin; }
    // Now proceed to ignore all comment lines after the header.
    char checkCommentLine = fin.peek();
    while (checkCommentLine == '%') {
      if (checkCommentLine == '%') {
        fin.ignore(MaxLineLength,'\n');
        nComments++;
      } // end taking care of a comment line
      checkCommentLine = fin.peek();  // get next line's first char
    } // end comment line skipping
    // Import all actual data starting with the number of rows and columns.
    fin >> N >> M; // discard M here since it must be 1
    //cout << "Before destructor in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    a.~TVectorBag<T>();  // explicit destructor call since it is being redefined
    //cout << "After  destructor in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    a.Rows = N;  a.Cols = M;
    a.data = new T*[a.Cols];
    for (int i=0; i<a.Cols; i++) { a.data[i] = new T[a.Rows]; }
    a.ValidArray = Valid2DArrayFlag;
    for(int j=0; j<a.getCols(); j++) { for(int i=0; i<a.getRows(); i++) {
      fin >> a.data[j][i]; }} // end ij
    //cout << "After redefinition in fin " << a.Rows << " " << a.Cols << " " 
    //     << a.ValidArray << " " << a[1][1] << endl; // debugging
    if(a.VerboseErrorReports!=0) {
      cout << "Array Rows x Cols = " << N << " " << M << "  "
           << "Number of comment lines " << nComments << ".\n";
    } // end if VerboseErrorReports
  } // end if
  return fin;
};
/*
// Now, code the commuted x*vector operation as a NON-MEMBER function.
template <typename T> 
TVectorBag<T> operator*(T x, const TVectorBag<T>& a) {
  TVectorBag<T> c(a.Cols, a.Rows);
  if(a.isValid())
    for(int i=0; i<a.Rows; i++) { for(int j=0; j<a.Cols; j++) {
      c.data[j][i] = a.data[j][i]*x; }} // end ij
  return c;
};
*/
} // end namespace MSL

#endif

