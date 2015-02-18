/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Copyright � 2009  Peter Ronhovde                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 *  This program is free software: you can redistribute it and/or modify       *
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

The code is for a general purpose multiresolution algorithm (MRA) solver using the "absolute Potts model" (APM) presented in a paper by Peter Ronhovde and Zohar Nussinov (Phys. Rev. E vol. 80 article 016109).  I have trimmed the code, leaving the more essential code for general applications with our Potts model.

This algorithm differs slightly from the algorithm in the mentioned paper in that we use a geometric step size rather than a step size based on the predicted edge density (see below).

--------------------------------------------------------
Included files:
--------------------------------------------------------
README.txt ----> this file

RNMain.cpp ----> the main function that calls the MRA solver or the (single resolution) greedy solver.
RNMRA.hpp  ----> the MRA algorithm code
other *.h, *.hpp, *.cpp files ----> the rest of code base

Other files:
karate_answer_nooffset.txt ---->  the answer file corresponding to the known two cluster split of the well-known Zachary karate club example
karate_weighted_nooffset.gml ----> the corresponding network definition including weighted edges (initial gml file is by Mark Newman with the weighted edges from UCINet).  The file uses a zero node starting ID rather than one in most published representations of the network.
karateMRAExample.gif ----> an example MRA plot as applied to the karate club

example256.gml ----> a 256 node hierarchy that is analyzed in the PRE paper mentioned above.
example256_level2_raw.txt ----> the intermediate level 2 of the 256 node hierarchy with 5 clusters
example256_level3_raw.txt ----> the bottom level 3 of the 256 node hierarchy with 16 clusters

compileRN ----> the script file which will compile the program for the MRA solver and the single resolution greedy solver for both weighted and unweighted systems.  If the script is not executable, use the command "chmod 744 compileRN" to make it executable, and then execute it with "./compileRN"


--------------------------------------------------------
Compilation:
--------------------------------------------------------

I have not written a makefile, but the compilation time is relatively short even with optimization enabled.  One particular compile command is:
g++ RNMain.cpp msl/MSL_CMatrix.cpp msl/MSL_CMatrix_SparseW.cpp msl/MSL_CMatrix_Dense.cpp msl/ML_Utils.cpp clusterclasses.cpp msl/MSL_Stats1D.cpp -DMRA_SOLVER -DMSL_TCMWEIGHTED -o rnMRAw -O2
I have included a script file "compileRN" to compile the weighted and unweighted versions of the MRA solver and the (single resolution) greedy solver.  I have not tested the code in a MS Windows-based compiler.



--------------------------------------------------------
Program output:
--------------------------------------------------------

Console updates:
I use ANSI color codes to add color to the console output.  If the colors are causing problems (such as in vim), the color definitions are located in the file ./msl/ML_Utils.h where I have included a section with the color definitions being commented-out.

Output data:
I use OpenOffice to view the csv, comma separated value, output data files; but any spread sheet application should be able to view them.  (Colors are not used in the data file output.)  Briefly, two files output by the MRA solver are:
filename.csv  ----> The "complete" set of graph data for use if you want to analyze how the structure varies over all resolutions (only for the MRA solver).
filenameBest.txt ----> A text file what includes the "best" partition based on strong information correlations, a stable information plateau, and the lowest energy trial at that resolution (in that order of precedence).
In the interest of brevity, the main data columns that are needed are:  
Columns:
D ----> the value of gamma (our Potts model weight) for multiresolution analysis
E ----> the NMI average and standard deviation for the "best" trial (by average correlations over all pairs of replicas)
K ----> the value of partition Shannon entropy [H(A) average over all replicas]
O ----> the average cluster size q of each resolution (over all replicas)
Feel free to contact me at (ronhovde[at]hbar.wustl.edu) if the interpretation of the other data columns is unclear.
I have also included a sample plot in gif file format showing the application of the MRA to the Zachary karate club (weighted edges).


For simple systems with a single strongly defined resolution (range), the automated "best" partition will select a proper solution.  
For many problems, human interpretation of the MRA data is helpful if not necessary.


--------------------------------------------------------
Running the Program:
--------------------------------------------------------

I make frequent use of command line parameters.  The order of parameters does not matter.  

Example:  Executing the program(s)
A sample function call applied to the (weighted) Zachary karate club (included data file starts with a node ID of zero rather than one) is:
./rnMRAw -n:34 -nsm:300 -nim:100 -nrm:10 -ntm:4 -gstart:19 -gstop:0.001 -gsteps:20 -v:2 -rseed:8231593 -inf:karate_weighted_nooffset.gml -infans:karate_answer_nooffset.txt -outf:MRAKarateTestW

This command writes two files:
MRAKarateTestW.csv     ----> the comma separated value file with the full set of multiresolution statistics (see above description for the most important columns)
MRAKarateTestWBest.txt ----> the "best" answer (very high and most "stable" NMI correlations)

With this command, we see that the "best" (most stable) MRA solution is different than the known two cluster split.  From a plot of the MRA csv data file (gif file plot is included with zip file), we see another perfectly correlated two-cluster solution at \gamma ~= 0.1 (in addition to some other variations of these partition themes).  

We can use the single resolution solver to display this solution with the command:
./rnCDw -n:34 -nsm:300 -nim:1000 -ntm:10 -g:0.1 -v:2 -rseed:2436083 -inf:karate_weighted_nooffset.gml -infans:karate_answer_nooffset.txt -outf:APMCDTestW
Here we see that the answer agrees exactly with the actual known split of the karate club.



Other comments and explanation:
Of course, no "answer" is known on general problems, so simply omit the -infans:filename.txt parameter.  

To run the "unweighted" version, omit the 'w' at the end of the program name (if my included script was used to compile the program).  All edges will then be treated as unweighted.

All parameters have default values, but some of them are set with other tests in mind, so care should be used if relying on default values.  (Most of the command line parameters and their defaults are defined in cluster_clparams.cpp with a few defined at the top of RNMain.cpp.)

Directed networks are not currently implemenented in the code due to legacy function implementations that assume symmetric matrices.  To solve weighted networks, the compile parameter -DMSL_TCMWEIGHTED was added to the compile command (see the script compileRN).  This turns on the explicit connection matrix energy evaluation within the cluster solve and merge routines; and as a result, it will slow down the solver by a factor of approximately log k where k is the average node degree when using the sparse matrix structure.


The command line options used above are:
--------------------------------------------------------
-nrm:100 ---->  The number of replicas for the MRA solver [generally O(10)]
-nsm:300 ---->  The maximum number possible resolutions to examine (a failsafe condition)
-nim:100 ---->  The maximum number of iterations for one community solution attempt.  Each iteration is one complete pass through all nodes.  Generally, the number of iterations for the best resolution(s) is O(1) for small and medium sized systems.  Very large systems or "difficult" problems are generally O(10).
-ntm:10  ---->  The number of optimization trials to use for each replica [generally O(10) or less].

-n:34        ---->  The number of nodes.  This must always be specified.
-g:0.1       ---->  The floating point model weight \gamma for the single resolution solver
-gstart:19.0 ---->  The starting value of \gamma for the MRA solver
-gstop:0.001 ---->  The ending value of \gamma for the MRA solver.  If gstart > gstop, the solver steps down in gamma.  Additionally, the run will also terminate as soon as the system is found to fully collapse.
-gsteps:20  ---->  The number of geometric steps *per decade* of \gamma

-rseed:741023 ----> Set the random seed.  Unfortunately, automatic seed generation is on the "to do" list, so they must be specified manually.
-v:1          ----> the "verbosity" of the update information for the console output: -v=0 is the most terse and -v=4 is the generally the most verbose.


-inf:inputfilename.gml ----> The gml file that contains the connection matrix definition.  This file is required.

-outf:filename ----> The "root" filename used for the data output.  Note that it *does not* include the ".csv" filename suffix.  The two files output by the code are:
"filename.csv" is the "complete" set of graph data (only for MRA solver)
"filenameBest.txt" is the selected set of "best" answer.
This output root file name is required.  The first file is output in csv format, and a number of trial parameters are exported there.  The second file is a "raw" cluster data with each cluster on a separate line and each node separated by a space.

-infans:karate_answer_nooffset.txt ----> the file containing the "known" answer (for testing purposes).  Simply omit the parameter when no answer is known.

Peter Ronhovde
Updated October 15, 2009

