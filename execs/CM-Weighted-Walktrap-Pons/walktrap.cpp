// File: walktrap.cpp
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons 
// Email    : pons@liafa.jussieu.fr
// Web page : http://www.liafa.jussieu.fr/~pons/
// Location : Paris, France
// Time	    : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "graph.h"
#include "communities.h"
#include <ctime>
#include <set>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

void info() {
  cerr << "WalkTrap v0.2 -- Finds community structure of networks using random walks." << endl;
  cerr << "Copyright (C) 2004 Pascal Pons." << endl << endl;
}

void help(char* prog_name) {
  info();
  cerr << "usage: " << prog_name << " [input_file] [-o output_file] [-i index_file] [options]" << endl << endl;
  cerr << "input_file: read the network from this file. stdin if not specified" << endl;
  cerr << "output_file: stdout if not specified" << endl;
  cerr << "index_file: index of the real name of the vertices" << endl << endl;
  cerr << "options:" << endl;
  cerr << "-s\t(silent) do not display the progress." << endl;
  cerr << "-tx\tset the length of random walks to x (default = 4)." << endl;
  cerr << "-dx\tset to x the detail level of the output (0 <= x <= 4, default = 1)." << endl;
  cerr << "-px\tprint the partition with x communities" << endl;
  cerr << "-b\tprint the best modularity partition" << endl;
  cerr << "-mx\tlimit the memory usage to x MB" << endl;
  cerr << "-h\tprint this help" << endl << endl;
  cerr << "see readme.txt for more details" << endl;
  exit(0);
}



int main(int argc, char** argv)
{  
  int length = 4;
  int details = 2;
  long max_memory = -1;
  bool silent = false;
  bool print_best_modularity = false;
  char* output_file = 0;
  char* input_file = 0;
  char* index_file = 0;
  set<int> print_partition;

  for (int i = 1; i < argc; i++) 
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
	case 't':
	  length = atoi(argv[i]+2);
	  if (length <= 0) help(argv[0]);
	  break;
	case 's':
	  if(argv[i][2] != 0) help(argv[0]);
	  silent = true;
	  break;
	case 'o':
	  if(argv[i][2] != 0) help(argv[0]);
	  if(++i < argc) 
	    if(!output_file) {
	      output_file = argv[i];
	      break;
	    }
	  help(argv[0]);
	case 'd':
	  details = atoi(argv[i]+2);
	  if((argv[i][2] != '0' || argv[i][3] != 0) && details <= 0) help(argv[0]);
	  break;
	case 'p':
	  if (atoi(argv[i]+2) <= 0) help(argv[0]);
	  print_partition.insert(atoi(argv[i]+2));
	  break;
	case 'b':
	  print_best_modularity = true;
	  break;	  
	case 'm':
	  max_memory = atol(argv[i]+2)*1024*1024;
	  if (max_memory <= 0) help(argv[0]);
	  break;
	case 'i':
	  if(argv[i][2] != 0) help(argv[0]);
	  if(++i < argc) 
	    if(!index_file) {
	      index_file = argv[i];
	      break;
	    }
	  help(argv[0]);
	default:
	  help(argv[0]);
      }
    }
    else {
      if(!input_file) input_file = argv[i];
      else help(argv[0]);
    }

  if(output_file) {
    ifstream ftmp;
    ftmp.open(output_file, ios::in);
    if(!ftmp.fail()) {cerr << "file " << output_file << " already exists" << endl; exit(0);}
    freopen(output_file, "w", stdout);
  }

  if(!silent) info();
  time_t begin = clock();
  
  Graph* G = new Graph;

  if(input_file) {
    ifstream f;
    f.open(input_file, ios::in);
    if(!f.is_open()) {
      cerr << "unable to open file " << input_file << endl; 
      exit(0);
    }
    if(!silent) cerr << "loading graph from file: \"" << input_file << "\""<<  endl;
    f >> *G;
    f.close();
  }
  else {
    if(!silent) cerr << "loading graph from standard input..." << endl;
    cin >> *G;
  }
  if(!silent) cerr << G->nb_vertices << " vertices and " << G->nb_edges << " edges sucessfuly loaded." << endl << endl;

  if(index_file) {
    if(!silent) cerr << "loading index from file " << index_file << endl;
    if(G->load_index(index_file)) 
      if(!silent) cerr << "index sucessfuly loaded" << endl << endl;
  }
  
  Communities C(G, length, silent, details, max_memory);  

  if(!silent) cerr << "merging the communities:";

  while(!C.H->is_empty()) {
    C.merge_nearest_communities();
  }

  if(!silent) cerr << endl;

  if(print_best_modularity)
    C.print_best_modularity_partition();
  
  for(set<int>::iterator it = print_partition.begin(); it != print_partition.end(); ++it)
    if(*it <= G->nb_vertices) C.print_partition(*it);
  
  if(!silent) cerr << endl << "computation successfully terminated in " << double(clock() - begin) / double(CLOCKS_PER_SEC) << "s" << endl;
  delete G;
  if(output_file) fclose(stdout);
  return 0;
}

