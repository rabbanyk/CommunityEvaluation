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

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <argp.h>
#include <math.h>
#include "commfind.h"

#define min(a,b) ((a>b)? b : a)
#define max(a,b) ((a<b)? b : a)

typedef struct {
  double val;
  int index;
} intdouble;

typedef struct {
  double d;
  int n1,n2;
} distlink;

typedef struct {
  int nmember;
  int *members;
  int nlink;
  int *linkto;
  int minlink;
  int *distk;
  double *dist;
  double mindist;
  double *eij;
  double eii;
  double ai;
} community;

FILE *logfile;

/*************************  argp stuff ******************************/     
const char *argp_program_version = "commfind 0.2";
     
/* Program documentation. */
static char doc[] =
"Find communities in a network described by edges in file FILE, computing NEV eigenvectors.";
     
/* A description of the arguments we accept. */
static char args_doc[] = "FILE NEV";

/* The options we understand. */
static struct argp_option options[] = {
  {"Laplacian",  'L',      0, 0, "Use Laplacian eigenvectors" },
  {"logfile",    'l', "FILE", 0, "Name of the log file (default out.log)" },
  { 0 }
};

/* Used by `main' to communicate with `parse_opt'. */
struct arguments {
  char *filename;
  char *logfile;
  unsigned int flag,nev;
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  struct arguments *args = state->input;
  
  switch (key) {
  case ARGP_KEY_INIT:
    /* Default values. */
    args->flag = NORM_LAPL;
    args->logfile = "out.log";
    break;
  case 'L':
    args->flag = LAPLACIAN;
    break;
  case 'l':
    args->logfile = arg;
    break;
  case ARGP_KEY_ARG:
    switch (state->arg_num) {
    case 0:  /* first arg: filename */
      args->filename = arg;
      break;
    case 1:  /* second arg: number of eigenvectors */
      args->nev = strtoul(arg,0,0);
      if (args->nev == 0) {
	fprintf(stderr,"Wrong format for number of eigenvectors\n");
	exit(-3);
      }
      break;
    default:  /* Too many arguments. */
      argp_usage (state);
    }
    break;
  case ARGP_KEY_END:
    /* chech if the argument (filename) is present) */
    if (state->arg_num < 2)
      argp_usage (state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

/* argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };
     
/******************* end of argp stuff ************************/     

int compare_distlinks (const void *a, const void *b) {
  const distlink *da = (const distlink *) a;
  const distlink *db = (const distlink *) b;
  
  return (da->d > db->d) - (da->d < db->d);
}
 
double angledist(double *c1, double *c2, int d) {
  int k;
  double n1=0.0;
  double n2=0.0;
  double ps=0.0;
  double x;

  for (k=0; k<d; k++) {
    n1 += c1[k]*c1[k];
    n2 += c2[k]*c2[k];
    ps += c1[k]*c2[k];
  }
  x = ps/sqrt(n1)/sqrt(n2);
  if (x > 1) 
    x = 1;
  if (x < -1)
    x = -1;
  return acos(x);
}

void comm_findmaxq(int n, graph *g, double *evec, int uev,
		   int *e1, int *e2, double *maxq, int *kmax) {
  int k,j, nl, dk, i1, i2, l1, l2, rl1, rl2, oldn, kk;
  double *c, q;
  double *coords[n];
  community *comm[n];
  community *c1, *c2, *cl;
  distlink *dlink;

  /* get coordinates from eigenvectors */
  c = malloc(n*uev*sizeof(double));
  for (k=0; k<n; k++) {
    coords[k] = c + uev*k;
    for (j=0; j<uev; j++)
      coords[k][j] = evec[(j+1)*n+k];
  }

  /* calculate distance for all links */
  nl = 0;
  dlink = malloc(g->nlink*sizeof(distlink));
  for (k=0; k<n; k++) {
    for (j=0; j<g->node[k].deg; j++)
      if (k < g->node[k].link[j]) {
	dlink[nl].n1 = k;
	dlink[nl].n2 = g->node[k].link[j];
	dlink[nl].d = angledist(coords[k], coords[dlink[nl].n2], uev);
	nl++;
      }
  }
  free(c);
  assert(g->nlink == nl);

  qsort(dlink, nl, sizeof(distlink), compare_distlinks);

  /* create single member communitities */
  q = 0;
  for (k=0; k<n; k++){
    comm[k] = malloc(sizeof(community));
    comm[k]->nmember = 1;
    comm[k]->members = malloc(sizeof(int));
    comm[k]->members[0] = k;
    dk = g->node[k].deg;
    comm[k]->nlink = dk;
    comm[k]->linkto = malloc(dk*sizeof(int));
    comm[k]->distk = malloc(dk*sizeof(int));
    comm[k]->eij = malloc(dk*sizeof(double));
    for (j=0; j<dk; j++) {
      comm[k]->linkto[j] = g->node[k].link[j];
      comm[k]->eij[j] = 0.5/nl;
    }
    comm[k]->ai = 0.5*dk/nl;
    q -= comm[k]->ai*comm[k]->ai;
  }
  /* fill distk values */
  for (k=0; k<nl; k++) {
    c1 = comm[dlink[k].n1];
    for (j=0; j<c1->nlink; j++)                    /* find link to n2 */
      if (c1->linkto[j] == dlink[k].n2)
	break;
    assert(j<c1->nlink);
    c1->distk[j] = k;
    /* do the same for reverse link */
    c2 = comm[dlink[k].n2];
    for (j=0; j<c2->nlink; j++)                    /* find link to n1 */
      if (c2->linkto[j] == dlink[k].n1)
	break;
    assert(j<c2->nlink);
    c2->distk[j] = k;
  }

  *maxq = q;
  *kmax = 0;
  kk = 0;
  for (k=0; k<nl; k++) {
    if (dlink[k].d < 0)
      continue;
    e1[kk] = dlink[k].n1;
    e2[kk] = dlink[k].n2;
    c1 = comm[e1[kk]];
    c2 = comm[e2[kk]];
    assert(c1 != c2); 
    i1 = c1->members[0];
    i2 = c2->members[0];
    
    /* find link from c1 to c2 */
    for (j=0; j<c1->nlink; j++)
      if (c1->linkto[j] == i2)
	break;
    assert(j<c1->nlink);     /* link found */
    q += 2*(c1->eij[j] - c1->ai*c2->ai);

    kk++;
    if (q>*maxq) {
      *maxq = q;
      *kmax = kk;
    }    

    /* add c1 and c2 data */
    c1->ai += c2->ai;
    for (l2=0; l2<c2->nlink; l2++) {                 /* adjusting links */
      cl = comm[c2->linkto[l2]];
      if (cl == c2)
	continue;
      for (rl2=0; rl2<cl->nlink; rl2++)            /* find reverse link */
	if (cl->linkto[rl2] == i2)
	  break;
      assert(rl2<cl->nlink);
      /* look for link to the same community in c1 */
      for (l1=0; l1<c1->nlink; l1++)
	if (c1->linkto[l1] == c2->linkto[l2])
	  break;
      if (l1<c1->nlink) {           /* found link to the same community */
	for (rl1=0; rl1<cl->nlink; rl1++)          /* find reverse link */
	  if (cl->linkto[rl1] == i1)
	    break;
	assert(rl1<cl->nlink);
	c1->eij[l1] += c2->eij[l2];
	/* find shortest link */
	if (dlink[c1->distk[l1]].d > dlink[c2->distk[l2]].d) 
	  dlink[c2->distk[l2]].d = -1;
	else {
	  dlink[c1->distk[l1]].d = -1;
	  c1->distk[l1] = c2->distk[l2];
	}
	cl->eij[rl1] += cl->eij[rl2];                      /* add eij's */
	cl->distk[rl1] = c1->distk[l1];
	cl->eij[rl2] = cl->eij[cl->nlink-1];       /* delete link to c2 */
	cl->linkto[rl2] = cl->linkto[cl->nlink-1];
	cl->distk[rl2] = cl->distk[cl->nlink-1];
	cl->nlink--;
	cl->eij = realloc(cl->eij, (cl->nlink)*sizeof(double));
	cl->linkto = realloc(cl->linkto, (cl->nlink)*sizeof(int));
	cl->distk = realloc(cl->distk, (cl->nlink)*sizeof(int));
      } else {                 /* link to the same community not found */
	c1->nlink++;
	c1->eij = realloc(c1->eij, (c1->nlink)*sizeof(double));
	c1->linkto = realloc(c1->linkto, (c1->nlink)*sizeof(int));
	c1->distk = realloc(c1->distk, (c1->nlink)*sizeof(int));
	c1->eij[c1->nlink-1] = c2->eij[l2];
	c1->linkto[c1->nlink-1] = c2->linkto[l2];
	c1->distk[c1->nlink-1] = c2->distk[l2];
	cl->linkto[rl2] = i1;
      }
    }
    /* find c1's link to c2 */
    for (l2=0; l2<c1->nlink; l2++)
      if (c1->linkto[l2] == i2)
	break;
    if (l2<c1->nlink) {                            /* if link exists */
      /* find link to c1 itself if it exists */
      for (l1=0; l1<c1->nlink; l1++)
	if (c1->linkto[l1] == i1)
	  break;
      if (l1<c1->nlink) {               /* previous link to c1 found */
	c1->eij[l1] += c1->eij[l2];                    /* add eij's  */
	c1->eij[l2] = c1->eij[c1->nlink-1];     /* delete link to c2 */
	c1->linkto[l2] = c1->linkto[c1->nlink-1];
	c1->distk[l2] = c1->distk[c1->nlink-1];
	c1->nlink--;
	c1->eij = realloc(c1->eij, (c1->nlink)*sizeof(double));
	c1->linkto = realloc(c1->linkto, (c1->nlink)*sizeof(int));
	c1->distk = realloc(c1->distk, (c1->nlink)*sizeof(int));
      } else                              /* no previous links to c1 */
	c1->linkto[l2] = i1;
    }
    /* add members */
    oldn = c1->nmember;
    c1->nmember += c2->nmember;
    c1->members = realloc(c1->members, c1->nmember*sizeof(int));
    for (j=0; j<c2->nmember; j++) {
      c1->members[oldn+j] = c2->members[j];
      comm[c2->members[j]] = c1;
    } 
    free(c2->members);
    free(c2->linkto);
    free(c2->distk);
    free(c2->eij);
    free(c2);
  }

  assert(n==comm[0]->nmember);
  free(comm[0]->members);
  free(comm[0]->linkto);
  free(comm[0]->distk);
  free(comm[0]->eij);
  free(comm[0]);
  free(dlink);
}

 void print_communities(graph *g, int *e1, int *e2, int kmax) { 
  int k, j, oldn, nc;
  int n = g->n;
  int ind[n];
  community *comm[n];
  community *c1, *c2, *cl;

  /* reset communities */
  for (k=0; k<n; k++){
    comm[k] = malloc(sizeof(community));
    comm[k]->nmember = 1;
    comm[k]->members = malloc(sizeof(int));
    comm[k]->members[0] = k;
  }
  
  for (k=0; k<kmax; k++) {
    c1 = comm[e1[k]];
    c2 = comm[e2[k]];
    assert(c1 != c2);
    /* swap communities so that the second one has a smaller # of members */
    if (c1->nmember < c2->nmember) {
      cl = c1;
      c1 = c2;
      c2 = cl;
    }
    /* add members */
    oldn = c1->nmember;
    c1->nmember += c2->nmember;
    c1->members = realloc(c1->members, c1->nmember*sizeof(int));
    for (j=0; j<c2->nmember; j++) {
      c1->members[oldn+j] = c2->members[j];
      comm[c2->members[j]] = c1;
    }
    free(c2->members);
    free(c2);
  }
  
  nc = 0;
  for (k=0; k<n; k++) {
    cl = comm[k];
    if (cl != NULL) {
      for (j=0; j<cl->nmember; j++) {
	/* assign index to all members of the community */
	ind[cl->members[j]] = nc;
	comm[cl->members[j]] = NULL;
      }
      nc++;
      free(cl->members);
      free(cl);
    }
  }

  /* print node indexes */
  for (k=0; k<n; k++) 
    printf("%d %d\n", k, ind[k]);
  fprintf(logfile,"%d communities found\n", nc);
 } 

void find_comm(graph *g, unsigned int nev, unsigned int evmode) {

  int k, uev, umax, ni, nf, nn;
  int n = g->n;
  int *kmax;
  double *maxq;
  double *eval;
  double *evec;
  int **e1,**e2;

  int evinfo = 0;
  int maxitr = 5000;
  int ncv = 30;

  eval = malloc(nev*sizeof(double));
  evec = malloc(n*nev*sizeof(double));

  /* calculating eigenvalues and eigenvectors */
  gr_eig(g, &nev, 1, 1, eval, evec, ncv, maxitr, evinfo, evmode);

  /* write EV in logfile */
  fprintf(logfile, "number of converged eigenvalues: %d\n", nev);
  for (k=0; k<nev; k++)
    fprintf(logfile, "%f\n",eval[k]);

  /* angular distance needs at least two coordinates */
  ni = 2;
  nf = nev-1;
  nn = nf-ni+1;

  kmax = malloc(nn*sizeof(int));
  maxq = malloc(nn*sizeof(double));
  e1 = malloc(nn*sizeof(int*));
  e2 = malloc(nn*sizeof(int*));
  for (k=0;k<nn; k++) {
    e1[k] = malloc((n-1)*sizeof(int));
    e2[k] = malloc((n-1)*sizeof(int));
  }
  
  /* loop on the number of used eigenvectors */
  umax = 0;
  for (uev=ni; uev<=nf; uev++) {
    comm_findmaxq(n, g, evec, uev, e1[uev-ni],e2[uev-ni], &(maxq[uev-ni]),
		  &(kmax[uev-ni]));
    fprintf(logfile, "D= %d\tQmax = %g\n", uev, maxq[uev-ni]);
    /* find max modularity */
    if (umax>0) {
      if (maxq[uev-ni] > maxq[umax-ni])
	umax = uev;
    } else
      umax = uev;
  }

  uev = umax;

  fprintf(logfile, "using %d eigenvalues: Qmax found = %f\n", 
	  uev, maxq[uev-ni]);
  print_communities(g, e1[uev-ni], e2[uev-ni], kmax[uev-ni]);

  for (k=0;k<nn; k++) {
    free(e1[k]);
    free(e2[k]);
  }
  free(e1);
  free(e2);
  free(kmax);
  free(maxq);
  free(eval);
  free(evec);
}

/* add an edge between two nodes in the graph */
void addlink(graph *g, int i1, int i2) {
  node *n1 = g->node+i1;
  node *n2 = g->node+i2;
  
  n1->deg++;
  n1->link = realloc(n1->link,n1->deg*sizeof(int));
  n1->link[n1->deg-1] = i2;
  n2->deg++;
  n2->link = realloc(n2->link,n2->deg*sizeof(int));
  n2->link[n2->deg-1] = i1;

  g->nlink++;
}

/* read edge list from file and build the graph */
graph *graph_from_file(const char *filename) {
  int k, minn, maxn;
  int ne = 0;
  int buflen = 0;
  unsigned int *n1 = NULL;
  unsigned int *n2 = NULL;
  int llength;
  char *line=NULL;
  FILE *fd = fopen(filename, "r");
  graph *gp = malloc(sizeof(graph));

  if (fd == NULL) {
    fprintf(stderr,"ERROR trying to open %s\n",filename);
    perror("in graph_from_file()");
    return NULL;
  }
  /* read edges */
  while (getline(&line,&llength,fd)>0) {
    if (ne > buflen-1) {
      buflen += 1024;
      n1 = realloc(n1, buflen*sizeof(unsigned int));
      n2 = realloc(n2, buflen*sizeof(unsigned int));
    }
    /* look for bad lines*/
    if (sscanf(line, "%u %u\n", &(n1[ne]), &(n2[ne])) < 2) {
      fclose(fd);
      return NULL;
    }
    ne++;
  }
  fclose(fd);

  /* look for graph size (and index of first node)*/
  minn = min(n1[0],n2[0]);
  maxn = max(n1[0],n2[0]);
  for (k=0; k<ne; k++) {
    if (n1[k] < minn)
      minn = n1[k];
    if (n2[k] < minn)
      minn = n2[k];
    if (n1[k] > maxn)
      maxn = n1[k];
    if (n2[k] > maxn)
      maxn = n2[k];
  }

  /* build graph */
  gp->n = maxn-minn+1;
  gp->node = calloc(gp->n,sizeof(node));
  gp->nlink = 0;
  for (k=0; k<ne; k++)
    addlink(gp, n1[k]-minn, n2[k]-minn);

  free(line);
  free(n1);
  free(n2);
  return gp;
}

/* check if the graph is connected */
int connected_graph(graph *g) {
  int k, llen, lnext, l, i;
  int n = g->n;
  int *mark = calloc(n, sizeof(int));
  int list[n];

  llen = 0;
  lnext = 0;
  list[llen++] = 1;
  mark[1] = 1;
  while (lnext<llen) {
    l = list[lnext++];
    for (k=0; k<g->node[l].deg; ++k) {
      i = g->node[l].link[k];
      if (!mark[i]) {
	mark[i] = 1;
	list[llen++] = i;
      }
    }
  }
  free(mark);
  return (llen==n);
}

void freegraph(graph *g) {
  int k;
  int n = g->n;

  for (k=0; k<n; k++)
    free(g->node[k].link);
  free(g->node);
  free(g);
}

int main(int argc, char **argv) {
  graph *g;
  struct arguments args;

  argp_parse (&argp, argc, argv, 0, 0, &args); /* parse command line */

  logfile = fopen(args.logfile,"w");
  if (logfile==NULL) {
    perror("could not open log file for writing");
    exit(-4);
  }
  fprintf(logfile, "file to read: %s\n", args.filename);
  fprintf(logfile, "number of eigenpairs to be computed: %u\n", args.nev);

  g = graph_from_file(args.filename);

  if (g == NULL) {
    fprintf(stderr, "could not build network from file %s\n", args.filename);
    exit(-2);
  }

  /* check if the graph is connected */
  if (!connected_graph(g)) {
    fprintf(stderr, "The network is not connected: exiting\n");
    exit(-5);
  }
  /* check nev */
  if (args.nev<3 || args.nev>=g->n) {
    fprintf(stderr, "number of eigenpairs out of range:\n");
    fprintf(stderr, "it must be larger than 2 and smaller than the number of nodes\n");
    exit(-1);
  }

  find_comm(g,args.nev,args.flag);

  freegraph(g);
  fclose(logfile);
  exit(0);
}
