Extract the gzipped tar archive, run 'make' to compile and, for example, './infomap 345234 flow_undir.net 10' to run the code.

tar xzvf infomap_undir.tgz
cd infomap_undir
make
./infomap 345234 flow_undir.net 10

Here ./infomap is the name of executable, 345234 is a random seed (can be any positive integer value), flow_undir.net is the network to partition (in Pajek format), and 10 is the number of attempts to partition the network (can be any integer value equal or larger than 1). 

This code can handle undirected networks with or without weighted links (no need to be integer weights). If a link occurs more than once in the network file (same or opposite direction), the weights are aggregated. 

Here is an example of a small network with three undirected weighted links (for more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/):   

*Vertices 3                                                                 
1 "Name of first node"                                                      
2 "Name of second node"                                                     
3 "Name of third node"                                                      
*Edges 3                                                                    
1 2 1.0                                                                     
1 3 3.3                                                                     
2 3 2.2                                                                     

The output file has the extension .tree and corresponds to the best partition (shortest description length) of the attempts. The output format has the pattern
# Code length 3.15098 in 4 modules.
1:1 0.075 "Node 2"
1:2 0.075 "Node 1"
1:3 0.05 "Node 4"
1:4 0.05 "Node 3"
2:1 0.075 "Node 6"
2:2 0.075 "Node 5"
2:3 0.05 "Node 8"
2:4 0.05 "Node 7"
3:1 0.075 "Node 10"
3:2 0.075 "Node 9"
3:3 0.05 "Node 12"
3:4 0.05 "Node 11"
4:1 0.075 "Node 14"
4:2 0.075 "Node 13"
4:3 0.05 "Node 16"
4:4 0.05 "Node 15"

For each row, except the first one, which summarizes the result, the first number is the module assignment, the second number is the rank within the module, the decimal number is the steady state population of random walkers, and within quotation marks is the node name. 

Results are also written to three files in Pajek format. The partition file with extension .clu should be used together with the original network (import both files and use the command Draw->Draw-Partition in Pajek). The file with extension _map.net is a network with all nodes and links aggregated according to the modular map. Finally, the vector file with extension _map.vec gives the size of the modules and should be used together with the file _map.net (import both files and use the command Draw->Draw-Vector in Pajek).
