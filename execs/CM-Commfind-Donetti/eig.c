/*
** Copyright (c) 2005, Luca Donetti
**
** This file is part of commfind.
**
** Commfind is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**   
** Commfind is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with commfind; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
** Author contact information:
**   donetti@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "commfind.h"

struct { 
  int logfil, ndigit, mgetv0;
  int msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  int mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
  int mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
} F77NAME(debug);

// double precision symmetric routines.

void F77NAME(dsaupd)(int *ido, char *bmat, int *n, char *which,
		     int *nev, double *tol, double *resid,
		     int *ncv, double *V, int *ldv,
		     int *iparam, int *ipntr, double *workd,
		     double *workl, int *lworkl, int *info);

void F77NAME(dseupd)(int *rvec, char *HowMny, int *select,
		     double *d, double *Z, int *ldz,
		     double *sigma, char *bmat, int *n,
		     char *which, int *nev, double *tol,
		     double *resid, int *ncv, double *V,
		     int *ldv, int *iparam, int *ipntr,
		     double *workd, double *workl,
		     int *lworkl, int *info);

inline void debug(int dbg) {
  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = -3;
  F77NAME(debug).mgetv0 = 0;
  F77NAME(debug).msaupd = dbg;
  F77NAME(debug).msaup2 = 0;
  F77NAME(debug).msaitr = 0;
  F77NAME(debug).mseigt = 0;
  F77NAME(debug).msapps = 0;
  F77NAME(debug).msgets = 0;
  F77NAME(debug).mseupd = 0;
  F77NAME(debug).mnaupd = 0;
  F77NAME(debug).mnaup2 = 0;
  F77NAME(debug).mnaitr = 0;
  F77NAME(debug).mneigt = 0;
  F77NAME(debug).mnapps = 0;
  F77NAME(debug).mngets = 0;
  F77NAME(debug).mneupd = 0;
  F77NAME(debug).mcaupd = 0;
  F77NAME(debug).mcaup2 = 0;
  F77NAME(debug).mcaitr = 0;
  F77NAME(debug).mceigt = 0;
  F77NAME(debug).mcapps = 0;
  F77NAME(debug).mcgets = 0;
  F77NAME(debug).mceupd = 0;
}

/* y = L_g x */
void lap_mult(graph *g, double *x, double *y) {
  int k,j;
  node nd;
  
  for (k=0; k<g->n; k++) {
    nd = g->node[k];
    y[k] = x[k]*nd.deg;
    for (j=0; j<nd.deg; j++)
      y[k] -= x[nd.link[j]];
  }
}

void lap_norm_mult(graph *g, double *x, double *y) {
  int k,j,maxdeg, deg;
  static double *tmp=NULL;
  static double *sqrtdeg=NULL;
  
  
  if (g==NULL) {  /* end of calculation: reset static vars */
    free(sqrtdeg);
    sqrtdeg=NULL;
    free(tmp);
    tmp=NULL;
    return;
  }

  if (sqrtdeg==NULL) { /* start of calculation: initialize static vars */
    maxdeg = 0;
    for (k=0; k<g->n; k++)
      if (g->node[k].deg>maxdeg)
	maxdeg =g->node[k].deg;
    sqrtdeg = malloc((maxdeg+1)*sizeof(double));
    for (k=1; k<=maxdeg; ++k)
      sqrtdeg[k] = sqrt(k);
    tmp = realloc(tmp, g->n*sizeof(double));
  }

  /* tmp = D^(-1/2) x */
  for (k=0; k<g->n; k++)
    tmp[k] = x[k]/sqrtdeg[g->node[k].deg];
  /* y = -A tmp */
  for (k=0; k<g->n; k++) {
    y[k] = 0;
    deg = g->node[k].deg;
    for (j=0; j<deg; j++)
      y[k] -= tmp[ g->node[k].link[j]];
  }
  /* y = x + D^(-1/2) y */
  for (k=0; k<g->n; k++) {
    y[k] /= sqrtdeg[g->node[k].deg];
    y[k] += x[k];
  }
}

int gr_eig(graph *g, int *nev, int smallest, int rvec,
	   double *eval, double *evec,
	   int ncv, int maxitr, int dbg, int mode) {
  int ido = 0;
  int ldv = g->n;
  int iparam[11];
  int ipntr[11];
  int lworkl;
  int *select;
  int info = 0;
  int ierr;
  char bmat = 'I';
  char *which;
  double tol = 1.0e-4;
  double resid[g->n];
  double *v;
  double workd[3*(g->n)];
  double *workl;
  double sigma;
  void (*mult_func)(graph *, double *, double *) = NULL;
  int k,j;
  double fdeg;


  iparam[0] = 1;   /* ishfts */
  iparam[2] = maxitr; /* maxitr */
  iparam[6] = 1;   /* mode   */
  
  if (smallest)
    which = "SM";
  else 
    which = "LM";
  
  switch (mode) {
  case LAPLACIAN:
    mult_func = lap_mult;
    break;
  case NORM_LAPL:
    mult_func = lap_norm_mult;
    break;
  }

  if (ncv > g->n)
    ncv = g->n;

  if (ncv <= *nev)
    ncv = (*nev)+1;

/*   printf("%d\t%d\t%d\n", ncv, *nev, g->n); */
  lworkl = ncv*(ncv+8);
  workl = malloc(lworkl*sizeof(double));
  v = malloc((g->n)*ncv*sizeof(double));

  debug(dbg);
  do {
    F77NAME(dsaupd)(&ido, &bmat, &(g->n), which, nev, &tol, resid, &ncv,
		    v, &ldv, iparam, ipntr, workd, workl,
		    &lworkl, &info);

    if ((ido == -1) || (ido == 1)) 
      mult_func(g,workd+ipntr[0]-1, workd+ipntr[1]-1);
    else
      break;
  } while(1);

  if (info < 0) {
    printf("dsaupd error %d: %d\n", info, iparam[4]);
    free(workl);
    free(v);
    return -1;
  }

  *nev = iparam[4];  /* number of converged eigenvalues */
  select = malloc(ncv*sizeof(int));

  F77NAME(dseupd)(&rvec, "All", select, eval, v, &ldv, &sigma, 
		  &bmat, &(g->n), which, nev, &tol, resid, &ncv, v, &ldv, 
		  iparam, ipntr, workd, workl, &lworkl, &ierr);
  free(workl);
  free(select);

  if (ierr != 0) { 
    printf("dseupd error %d\n", ierr);
    free(v);
    return -1;
  }

  if (rvec) {
    memcpy(evec, v, (*nev)*(g->n)*sizeof(double));
    if (mode==NORM_LAPL)  /* adjust eigenvectors */
      for (k=0; k<g->n; ++k) {
	fdeg = 1.0/sqrt(g->node[k].deg);
	for (j=k; j<(*nev)*(g->n); j+=g->n)
	  evec[j] *= fdeg;
      }
  }

  if (mode==NORM_LAPL) /* free static arrays in lap_norm_mult */
    lap_norm_mult(NULL, NULL, NULL);
  free(v);
  return 0;
}
