========================================================================
    Community detection by Community-Affiliation Graph Model 
========================================================================

The example implements community detection by the Cluster Affiliation Model 
for BIG networks (BIGCLAM).

This program detects network communities from a given network by fitting 
BIGCLAM, a probabilistic generative model for networks, to the given network 
by maximum likelihood estimation. User can specify how many communities 
she would detect, or let the program automatically determine the number 
of communities in the network based on the structure of the network.

This program supports multi-thread computing via OpenMP.

Fitting procedure and the cluster affiliation model for big networks are 
described in the following paper:
J. Yang and J. Leskovec, Overlapping Community Detection at Scale: 
A Nonnegative Matrix Factorization Approach, WSDM '13.
J. Yang and J. Leskovec, Community-Affiliation Graph Model for 
Overlapping Community Detection, ICDM '12.
J. Yang and J. Leskovec, Structure and Overlaps of Communities in 
Networks, ACM TIST '13.

The code works under Windows with Visual Studio or Cygwin with GCC,
Mac OS X, Linux and other Unix variants with GCC. Make sure that a
C++ compiler is installed on the system. Visual Studio project files
and makefiles are provided. For makefiles, compile the code with
"make opt".

////////////////////////////////////////////////////////////////////////
Parameters:
-nt: option for parallization. -nt:1 uses no parallelization -nt:N uses N threads
-i: input file name (edgelist by default.) ".ungraph" extension (SNAP binary file for PUNGraph) can be loaded as well
-o: output file name for detected community affiliation
-c: the number of communities to detect (-1:detect automatically by cross-validation)
-mc: minimum number of communities to try for cross-validation
-xc: maximum number of communities to try for cross-validation
-nc:how many numbers to try in cross-validation for the number of communities
-sa, -sb: Parameters for backtracking step size search. You can play with them if you want. 
sa, sb correspond to "Alpha" and "Beta" in the backtracking line search in the following book:
S. Boyd and L. Vandenberghe, Convex Optimization, Page 464.
////////////////////////////////////////////////////////////////////////
Usage:

-- To detect K communities from the network
./bigclam -c:K

-- To detect communities when you do not know how many communities you want. (It will search K from "L" different values from "m" to "M" based on Cross Validation likelihood)
./bigclam -c:-1 -mc:m -xc:M -nc:L

