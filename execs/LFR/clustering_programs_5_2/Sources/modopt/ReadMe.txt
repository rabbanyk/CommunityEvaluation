
To compile:
g++ main.cpp -O3 -o modopt

To run the code just type:

./modopt [edgelist]

-------------------------------------------
The edgelist file has to be in the format:

node1 node2

or 

node1 node2 weight

-------------------------------------------

There are other optional arguments but if you want to use them make
sure to type them in the correct order.

Ex:
./modopt network.dat 23 1.5 5 1e-3 0.9999
23 is the random seed for the random number generator
1.5 is the resolution parameter
5 the number of runs
1e-3 the initial temperature
0.9999 the temperature step for cooling down