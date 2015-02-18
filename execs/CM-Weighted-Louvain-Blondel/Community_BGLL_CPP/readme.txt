-----------------------------------------------------------------------------

Community detection
Based on the article "Fast unfolding of community hierarchies in large networks"
Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre

This program or any part of it must not be distributed without prior agreement 
of the above mentionned authors.

-----------------------------------------------------------------------------

Author   : E. Lefebvre, adapted by J.-L. Guillaume
Email    : jean-loup.guillaume@lip6.fr
Location : Paris, France
Time	 : February 2008

-----------------------------------------------------------------------------

Disclaimer:
This is the first public version of this program, if you find a bug, please
send a bug report to jean-loup.guillaume@lip6.fr including if necessary the 
input file and the parameters that caused the bug.

Note that the program does not make much verifications about the arguments,
and is expecting a friendly use. Optimization have to be made directly in 
the code to work with very large graphs since memory usage is not so clean,
in particular the following two lines:
community.cpp:169:  g2.links    = (int *)malloc((long)100000000*8);
community.cpp:170:  g2.weights  = (int *)malloc((long)100000000*8);
The function neigh_comm in community.cpp can be optimized.

You can also send me any comment or suggestion about the program.

-----------------------------------------------------------------------------

This package offers a set of functions to use in order to compute 
communities on graphs weighted or unweighted. A typical sequence of 
actions is:

1. Conversion from a text format (each line contains a couple "src dest")
Note that nodes are renumbered to be consecutive AND starting from 0.
This can be confusing, so take care:
./convert -i graph.txt -o graph.bin


2. Computes communities and displays hierarchical tree:
./community graph.bin -l -1 > graph.tree

To ensure a faster computation (with a loss of quality), one can use
the -q option to specify that the program must stop if the increase of
modularity is below epsilon for a given iteration or pass:
./community graph.bin -l -1 -q 0.0001 > graph.tree

The program can deal with weighted networks using -w option:
./community graph.bin -l -1 -w > graph.tree
In this specific case, the convertion step must also use the -w option.


3. Displays information on the tree structure (number of hierarchical
levels and nodes per level):
./hierarchy graph.tree

Displays the belonging of nodes to communities for a given level of
the tree:
./hierarchy graph.tree -l 2 > graph_node2comm_level2

-----------------------------------------------------------------------------

Known bugs (to be corrected in the next release):
- if the number of links or the total weight in the case of weighted networks
is higher than 2^32, the modularity computation is meaningless and the outcome
of the algorithm is undefined.

