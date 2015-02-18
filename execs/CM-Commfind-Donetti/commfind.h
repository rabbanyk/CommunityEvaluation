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

#define F77NAME(a) a ## _

#define LAPLACIAN 101
#define NORM_LAPL 102

typedef struct {
  int deg;             /* node degree */
  int *link;           /* array of neighbor indexes */
  int id;
} node;

typedef struct {
  node *node;          /* array of nodes */
  int n;               /* number of nodes */
  int nlink;           /* number of links */
} graph;

int gr_eig(graph *g, int *nev, int smallest, int rvec,
	   double *eval, double *evec,
	   int ncv, int maxitr, int dbg, int compl_mode);
