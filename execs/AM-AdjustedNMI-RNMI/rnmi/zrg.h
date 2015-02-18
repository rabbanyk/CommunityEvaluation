/*   Random number generator of
 *   sbm version 1.1, release date 30/05/2012
 *   Copyright 2012 Aurelien Decelle, Florent Krzakala, Lenka Zdeborova and Pan Zhang
 *   sbm is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or 
 *   (at your option) any later version.

 *   sbm is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
*/
#ifndef __ZRG_H__
#define __ZRG_H__

#include <cmath>

using namespace std;
class ZRANDOMv3
{
 public:
  double rdflt( void );

  double poidev( double );

  double gasdev( void );
  void set_seed(long);

  ZRANDOMv3( long i=0 );
void ranseq(int*,int);
void ranseq(double*,int);
 private:
  static const long IM1=2147483563;
  static const long IM2=2147483399;
  static const long IMM1=( IM1-1 );
  static const int IA1=40014;
  static const int IA2=40692;
  static const int IQ1=53668;
  static const int IQ2=52774;
  static const int IR1=12211;
  static const int IR2=3791;
  static const int NTAB=32;
  static const int NDIV=( 1+IMM1/NTAB );
  //static const double EPS=1.2e-7;
  //static const double RNMX=9.9999988e-1; //=1.0-EPS;
  //static const float AM=( 1.0/IM1 );
  double EPS;
  double RNMX;
  float AM;
  int tmpint;
  int tmpdouble;
  long iv[NTAB];
  long idum;
  long idum2;
  long iy;
  double gas_gset;
  int    gas_iset;
  double gas_rsq;
  double gas_v1;
  double gas_v2;
  double gas_fac;

  double GammlnVal[1000];
};
#endif //ndef __ZRG_H__
