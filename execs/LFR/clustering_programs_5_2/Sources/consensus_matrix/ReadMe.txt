


This program computed the consensus matrix given  a number of partitions.


1. Compile:
   type ./compile.sh from terminal


2. Run on an example:
   ./pr_consensus pts 0.5


*** INPUT files ***

pts is a file which reports all the partitions from which we can compute the consensus matrix. In the example, there are three files, each of them stores a different partition.
The format of the partition file is such that every module is separated by a new line.

For instance, one.dat looks like:

1	2	3
7	9

which means that there are two modules, in the former there are nodes 1,2 and 3.
Nodes have to be integer numbers, but there is no need for them to be consecutive.

The second number is a threshold. 
0.5 is suggested. Links with lower values will be kept only to avoid nodes to become isolated.


*** OUTPUT files ***

[pts]_net.dat

contains the consensus matrix in the format

[node] [node] [weight]



*** HOW to use the program ***

You need to collect a number of partitions from your favorite clustering algorithm, saving them in certain files. This program will compute the consensus matrix from them. Then, you can use your clustering algorithm again on the consensus matrix.
I will include some popular clustering algorithm soon. Thanks for your patience.





