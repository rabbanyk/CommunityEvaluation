       #=======================================#
       #  WalkTrap v0.2                        #
       #  Copyright (C) 2004-2005 Pascal Pons  #
       #=======================================#


1) DESCRIPTION
   ===========

WalkTrap is a C++ program that finds community structure of a network.
It is based on the fact that a random walker tends to be trapped in dense
part of a network corresponding to communities.


2) AUTHOR & COPYRIGHT
   ==================

This program is copyright (C) 2004-2005 by Pascal Pons:
    Email:	pons@liafa.jussieu.fr
    Web page:	http://www.liafa.jussieu.fr/~pons/

Collaborator: Matthieu Latapy
    Email:	latapy@liafa.jussieu.fr
    Web page:	http://www.liafa.jussieu.fr/~latapy/

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (gpl.txt); if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


3) USAGE
   =====

compiling:
    To compile the program type the command "make".
    If you do have not the "make" utility try : 
    "g++ -O3 walktrap.cpp graph.cpp communities.cpp -o communities"
    This program requires the C++ Standard Template Library and has been
    tested with gcc 3.2.3

command line usage:
    walktrap [input_file] [-o output_file] [-i index_file] [options]

input_file:
    The file in which the network is stored.
    If this parameter is omitted, stdin is used.
    Be careful that the input format has changed since version 0.1

    INPUT FORMAT: (see example.net)
    The format is a list of undirected weighted edges. The vertices must be
    encoded as consecutive integers starting from 0. Each line contains two
    vertices and a weight (separated by spaces or tabulations) that define a
    weighted edge. The weight may be omitted, in this case a default weight
    equal to 1.0 is considered. The multi-edges will be considered as a single
    edge which weight is the sum of the corresponding weights. A comment line
    starts with "#"

    The format is less flexible than the format used in version 0.1 but it 
    allows considering weighted networks as well as unweighted networks.
    You can use the graph converter tool "Gconvert" available on my web page
    http://www.liafa.jussieu.fr/~pons/. It handles different formats and can 
    generate an index of the real name of the vertices that can be used in
    the output (see index_file section)

output_file:
    The file in which the community structure will be stored.
    If this parameter is omitted, stdout is used.

    OUTPUT FORMAT:
    according to the options the output may contain :
    - The list of the first communities containing a single vertex.
    - The successive merging of communities with
	- the modularity Q of the partition obtained at this step.
	- the value of delta_sigma that has been chosen for this merging.
	- the vertices that belong to the new community.
	- the description of the whole partition.
    - The partition with the best modularity (-b option)
    - More partitions asked by user with the options -p

index_file:
    If you wish to keep the real name of the vertices in the output, you
    can specify an index file.

    INDEX FORMAT:
    The index is a list of the real name of the vertices. Each line begins
    with a vertex number (as in the input file) followed by its real name
    that may be an arbitrary string. All the vertices must be defined once.

options:
    -s	: (silent) 

    -tx	: set the length of random walks to x. Default value is t = 4.

    -dx : set to x the detail level of the output (1 <= d <= 5).
	  d  = 1 nothing is written.
	  d >= 2 the successive mergings of the communities are written.
	  d >= 3 the modularity Q and the value of delta_sigma are written.
	  d >= 4 the new community is written at each step.
	  d >= 5 the whole partition is written at each step.
	  If ommited default value is d = 2.


    -b	: at the end of the process, print the partition that corresponds
	  to the best value of the modularity.

    -px	: at the end of the process, print the whole partition that 
	  corresponds to x communities. This option can be used several
	  times.

    -mx : limit the memory usage of the program to x MB. This option is
	  very useful for very large networks. But the algorithm runs 
	  slower with less memory.

    -h	: print the command line usage.

examples:
    "walktrap network.dat -o communities.txt -t5 -b"
    read network from file "network.dat" and writes the community structure
    found in file "communities.txt". The computation is done with random
    walks of length 5 and the best modularity partition is printed.

    "walktrap example.net -o test -d1 -p3"
    only the partition into 3 communities computed by the program is printed.

    "walktrap example.net -o test2 -m400" 
    the program will never use more that 400MB of memory.

    "walktrap example.net -d1 -s -b"
    only the best modularity partition is printed to the screen.

4) MORE INFORMATION
   ================

more information is available in the paper:
Matthieu Latapy and Pascal Pons, 'Computing communities in large networks
using random walks', submitted preprint to which this program actually is
associated.

It may be dowloaded, as well as the source code and material for WalkTrap
on: http://www.liafa.jussieu.fr/~pons/

5) COMMENTS & BUG REPORT
   =====================

If you find a bug, please send a bug report to pons@liafa.jussieu.fr
including the input file and the parameters that caused the bug.

You can also send me any comment or suggestion about the program.

6) HISTORY
   =======

v0.2 (June 2005)
new features :
    - support of weighted networks (input format has been modified)
    - an efficient memory manager has been implemented
    - many optimizations on the probability vectors storage and computation
    - some heuristic optimizations in the merging process
    - a heap structure has been implemented to store the distances 
      between communities.

v0.1 (November 2004)
first public version of Walktrap.

					June 7th, 2005. Pascal Pons.
