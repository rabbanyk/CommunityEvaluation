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
MSL - My (or Math) Scientific Library - TArray base class

Comments:
- All memory allocation is dynamic.
- The arrays will use standard [col][row] zero indexed C/C++ multidimensional
arrays to prevent any confusion since it is a C++ library.
- I implemented my own classes rather than use an existing library in order to 
control the size of the classes as well as ensure that I am not taking any 
significant performance hit.
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
- The different dimension implementations of TArray implement a number of Vector
and Matrix like operations for obvious usefulness; but in order to reduce the class
size, the TArray#D class stops short of full TVector, TMatrix, or TTensor 
implementations for 1, 2, and 3 dimensional arrays respectively since they do not 
implement functions such as inverses, various decompositions, matrix 
multiplication, etc.   It makes some sense to separate an array from a matrix or 
vector since not all arrays are intended to be full matrices with appropriate available 
operations.  These various operations may be implemented in the future as 
non-member functions or as member functions of full TMatrix, etc. classes that 
inherit from TArray#D.
- Array '=' assignments cannot change the size of the object so that Arrays
with different sizes may not be assigned to each other.  There is a little debate 
(with myself) over which is the best way to implement the class.  Mathematically, 
one would expect that Arrays of different size would be incompatible for assignment
operations.  However, in actual practice, this turns out to be slightly problematic for
two reasons.  First, the user must know beforehand how big of an Array to create for
the assignment operation with functions that return an Array.  This is mostly just a 
convenience issue.  Second, the validation checks can add significant overhead 
for scientific applications if many assignments are performed.  (For example, some 
SOR iterative solutions for a 2D cylindrical system can require a million or more 
iterations since convergence is relatively slow.)
- A reimplementation of the Array2D and Array3D classes as internally large singly indexed vectors would allow
a strong optimization routine to at least the SOR algorithm.  See the GCG_SOR.h SORStep optimized 
algorithm function for an example of such an optimization.  This data structure would be particularly useful
in optimizing a 3D Poisson's Equation solver (most likely SOR initially).  After completion, it would also
be transparent to a general user.  The only possible drawback would be the requirement of always 
allocating large contiguous blocks of memory rather than many shorter blocks.  It might be a small
problem with large systems; however, since one would be restricted to 46000 x 46000 system which is 
not really much of a restriction on a single processor SOR or GCG type solver since systems this large
would take prohibitively long to solve (many years).
this case is really the target of the implementation, so if need be, allowances could be made to facilitate 
the optimization routine since it is very favorable to manually manage memory references (it eliminated
many integer multiplications be grid point that are required for array index referencing).
Note (Nov 2005) It does work for the Array2D class; however, it was found that (at least for 
the SOR method) there was only about a %1 performance gain over just manually handling memory 
management (SOR method in mind here) for the regular TFloat** version of Array2D.  Due to this, it 
may be questionable whether a full Array3D implementation would actually benefit the 3D SOR solver.
Note (Nov 2005).  This is difficult to do directly for Array3D since the reference [][][] operator does not
work since after a[] we have a TFloat** return value and the programer has no way to get the second
a[][] to properly reference the correct element of the single vector that is modelling the array.  

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

/* To-Do List (for all Array dimensions):  *********************************************************
Items are not in order of importance.
1. Make a "prettyPrint" function to output arrays in nice human readable format.
2. Resolve inconsistency in the implementation of the 1D and 2D constructors.
4. Check Array3D to ensure that different size sub Array2D's are not allowed.
7. Implement limited dot product and magnitude methods for Array2D.
8.  Change array definitions in IterationSolvers.h   M <-> N since col <-> row.
9. Add square Array2D and Array3D constructors.  Resolve related default constructor ambiguity.
11. Array2D (non LMB version) += operator appears to be broken.
12. The Array2D and Array3D classes should have a constructor for TFloat** or TFloat***
respectively along with a conversion routine back to a C pointer array.
(LOW PRIORITY ITEMS)
29. Add optional removal of a lot of error checking code by either pre-processor directives if needed.
24. Add sizeof() function for whole object, and individual columns of objects (may be done in 
by some compilers automatically?).
25. Write more memory optimized (manual pointer handling) assignment of Array3D loop routines
(actually makes very little difference though since the current Array3D is based on Array2D).
26. Output error messages for multidimensional arrays with only one dimension
specified or make them such that it assumes a square array.
21. Add additional error checking to stream io (end of file before data end, extra unused data, etc.)
22. Implement extra comment lines in Array#D output (implicitly allowed already).
23. Recheck validations especially with the assignments operators (made somewhat unneeded 
until some further major code writing since the explicit validations were taking up way too much computational time).
**************************************************************************************/

/* Future 'Wish' or 'Idea'-List:  *********************************************************
1. Implement full TMatrix and TVector classes with appropriate new methods.
  a. matrix multiplication
  b. matrix inversion (perhaps iterative)
  c. eigenvalue/eigenvector calculation
  d. various solvers:  LU, QR, SVD, etc. as needed
  e. other important decompositions: cholesky
2. Implement full TMatrix and TVector classes with appropriate new methods.
  a. row or column designator for TVector
  b. vector multiplication to inner or outer product
3. Other types of matrices (all real):  In all cases, I will have to re-implement all of the other standard 
operators and methods.
  A. triangular 
    1. upper or lower designator
    2. strict triangular also?
    3. triangular solver (back substitution, will probably be a function needed for LU decomposition anyhow)
    4. use different sized pointer arrays to appropriately reduce memory size
  B. symmetric or antisymmetric
    1. use different sized pointer arrays to appropriately reduce memory size
  C. tridiagonal
    1. use 3 Array1D's to reduce memory
    2. fast tridiagonal solver
  E. diagonal matrix?  (Is there as significant benefit as an explicit matrix?  easier for math operations, important?)
Less important:  (I haven't needed one of these yet unless I wanted to explicitly solve my Ax=b system)
  A. banded matrices?  
  B. sparse matrices?  
  C. Hessenburg matrices?
    1. implement as a triangular plus a vector
Just ideas (far out, if ever):
  A. block matrices?
  
**************************************************************************************/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef MSL_TARRAY_H
#define MSL_TARRAY_H

namespace MSL {

// define constants
const int Valid1DArrayFlag = 123450;
const int Valid2DArrayFlag = 234501;
const int Valid3DArrayFlag = 345012;
//const Array1D NullArray1D(0);

const char MMHeader[42] = "%%MatrixMarket matrix array real general\n";
const int  MaxLineLength = 1024*1024;

#ifndef MSL_TFLOAT_T
// TFloat allows easily changing between double and float data types
typedef double TFloat; // float type definition
#define MSL_TFLOAT_T
#endif

//------------------------------------------------------------------------------
// virtual Array base class
/*
Defines a base TArray class to hold basically only static common flags to unify Array
error checking and other operations for all instances of Arrays.
1. StrictErrorChecking=1 causes many more data validation checks to be made, and
prevents other things as well:
2. VerboseErrorReports is currently implemented sparingly not very "verbosely."  Extra
reporting will be added as needed.  See comment above about unremoved dead
code.

*/
class TArray {
 public:
  // Data elements
  static const bool StrictErrorChecking = 1;
  static const bool SloppyAssignments   = 0;
  static const bool VerboseErrorReports = 0;
  // Member functions
  virtual bool isValid() const = 0;
               TArray() {};
               TArray(TArray& a) {};
  virtual     ~TArray() {};
 protected:
  int ValidArray;
}; // end virtual class Array
//------------------------------------------------------------------------------


class TVector {
 public:
  // Data elements
  static const bool StrictErrorChecking = 1;
  static const bool SloppyAssignments   = 0;
  static const bool VerboseErrorReports = 0;
  // Member functions
  inline bool isValid(bool bEmptyOK = 0) const { return 0; };
        TVector(bool bEmptyOK = 0) {};
        TVector(TVector& a) {};
       ~TVector() {};
 protected:
  int ValidArray;
}; // end virtual class TVector


class TMatrix {
 public:
  // Data elements
  static const bool StrictErrorChecking = 1;
  static const bool SloppyAssignments   = 0;
  static const bool VerboseErrorReports = 0;
  // Member functions
  bool isValid(bool bEmptyOK = 0) const { return 0; };
  //virtual bool isValid(bool bEmptyOK = 0) = 0;
  //virtual bool isValid() = 0;
               TMatrix() {};
               TMatrix(TMatrix& a) {};
              ~TMatrix() {};
 protected:
  int ValidArray;
}; // end virtual class TMatrix


} // end namespace MSL
#endif
