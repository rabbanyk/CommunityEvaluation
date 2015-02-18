
from igraph import *
import sys, getopt

''' 
http://igraph.sourceforge.net/doc/python/igraph.GraphBase-class.html
'''

directed=False 
loops=False
attribute=None
multiple=False
n=1000
p=.5 
m=5000
repeat = 1
simplify = True
fw_prob =.3
bw_factor=0.0
ambs=1
k =1
dim=3
size=n/3 
nei=5
exponent_out=2.5 
exponent_in=-1
type_dist = pref_matrix = None
fitness_out = None
fitness_in=None
outpref=False
power=1
zero_appeal=1
finite_size_correction=True

def generate(method, ):
    if (method ==  'Erdos_Renyi'):
        print 'Erdos_Renyi'
        return (Graph.Erdos_Renyi(n, p, m, directed, loops))#m - the number of edges. If given, p must be missing.
    elif (method ==  'Barabasi'):
        print 'Barabasi' 
        return (Graph.Barabasi(n, m, outpref, directed, power, zero_appeal, implementation="psumtree", start_from=None))
    elif (method ==  'Establishment'):
        print 'Establishment' 
        return (Graph.Establishment(n, k, type_dist, pref_matrix, directed))
    elif (method ==  'Forest_Fire'):
        print 'Forest_Fire' 
        return (Graph.Forest_Fire(n, fw_prob, bw_factor, ambs, directed))
    elif (method ==  'Preference'):
        print 'Preference' 
        return (Graph.Preference(n, type_dist, pref_matrix, attribute, directed, loops))
    elif (method ==  'Static_Power_Law'):
        print 'Static_Power_Law' 
        return (Graph.Static_Power_Law(n, m, exponent_out, exponent_in, loops, multiple, finite_size_correction))
    elif (method ==  'Static_Fitness'):
        print 'Static_Fitness' 
        return (Graph.Static_Fitness(m, fitness_out, fitness_in, loops, multiple))
    elif (method ==  'Watts_Strogatz'):
        print 'Watts_Strogatz' 
        return (Graph.Watts_Strogatz(dim, size, nei, p, loops, multiple))
#     elif (method ==  'Degree_Sequence'):
#         print 'Degree_Sequence' 
#         return (Graph.Degree_Sequence(out, in=None, method="simple"))


def generateGraph(method, outputPath = 'network.gml', format = 'gml'):
    for r in range(repeat):
        g = generate(method)
      
        if outputPath<>None:
            for c in comm:
                g.write_gml(outputPath+str(r)+'.gml')
        else:  
            print list(comm)
        if outputPath<>None:
            f.close()
 
 
 
 
# TODO: not implemenmted 
# finite_size_correction=True
# type_dist = pref_matrix = None
# fitness_out = None
# fitness_in=None
# outpref=False

def main(argv):
    inputfile = None
    outputfile = None
    method = None
    try:
        opts, args = getopt.getopt(argv,"hr:o:m:n:p:e:f:b:",
                                   ["repeat=","ofile=","method=",
                                    "nodes=","prob=","edges=","untouched",
                                    "directed","loops","multiple","attribute"
                                    , "fw=","bw=","ams=","k=",
                                    "dim=", "size=", "nei=", 
                                    "exponent_out=","exponent_in=",
                                    "power=", "zero_appeal="])
    except getopt.GetoptError:
        print 'graphGeneratorInterface.py  -o <outputfile> -m <method>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'graphGeneratorInterface.py  -o <outputfile> -m <method>'
            sys.exit()
        elif opt in ("-r", "--repeat"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--method"):
            method = arg
        elif opt in ("-n", "--nodes"):
            n = int(arg)
        elif opt in ("-p", "--prob"):
            p = float(arg)
        elif opt in ("-e", "--edges"):
            m = int(arg)
        elif opt in ("-f", "--fw"):
            fw_prob = float(arg)
        elif opt in ("-b", "--bw"):
            bw_factor = float(arg)
        elif opt in ("--ambs"):
            ambs=int(arg)
        elif opt in ("--k"):
            k =int(arg)
        elif opt in ("--dim"):
            dim=int(arg)
        elif opt in ("--size"):
            size=int(arg) 
        elif opt in ("--nei"):
            nei=int(arg)
        elif opt in ("--exponent_out"):
            exponent_out= float(arg)
        elif opt in ("--exponent_in"):
            exponent_in= float(arg)
        elif opt in ("--power"):
            power= float(arg)
        elif opt in ("--zero_appeal"):
            zero_appeal=1       
        
        elif opt in ("--untouched"):
            simplify = False
        elif opt in ( "--directed"):
            directed = False
        elif opt in ( "--loops"):
            loops = False
        elif opt in ( "--multiple"):
            multiple = False
        elif opt in ( "--attributes"):
            attributes = arg

       
            
    if method==None or inputfile ==None  :
        print 'usage: graphGeneratorInterface.py  -o <outputfile> -m <method>'
        sys.exit() 
    else:
        findCommunities(method, inputfile, outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])    
    
    
    
#g = Graph.Erdos_Renyi(50, .1,  directed=False, loops=False)
'''
Parameters:
n - the number of vertices
m - either the number of outgoing edges generated for each vertex or a list containing the number of outgoing edges 
    for each vertex explicitly.
outpref - True if the out-degree of a given vertex should also increase its citation probability (as well as its in-degree), but it defaults to False.
directed - True if the generated graph should be directed (default: False).
power - the power constant of the nonlinear model. It can be omitted, an d in this case the usual linear model will be used.
zero_appeal - the attractivity of vertices with degree zero.
implementation - the algorithm to use to generate the network. Possible values are:
"bag": the algorithm that was the default in igraph before 0.6. It works by putting the ids of the vertices into a bag (multiset) exactly as many times as their in-degree, plus once more. The required number of cited vertices are then drawn from the bag with replacement. It works only for power=1 and zero_appeal=1.
"psumtree": this algorithm uses a partial prefix-sum tree to generate the graph. It does not generate multiple edges and it works for any values of power and zero_appeal.
"psumtree_multiple": similar to "psumtree", but it will generate multiple edges as well. igraph before 0.6 used this algorithm for powers other than 1.
start_from - if given and not None, this must be another Graph object. igraph will use this graph as a starting point for the preferential attachment model.
Reference: Barabasi, A-L and Albert, R. 1999. Emergence of scaling in random networks. Science, 286 509-512.
'''
# g = Graph.Barabasi(500, power=0.6, m=10)
# plot(g)
# g.write_gml('barba.gml')
'''
Forest_Fire(n, fw_prob, bw_factor=0.0, ambs=1, directed=False)
source code 
Generates a graph based on the forest fire model

The forest fire model is a growin graph model. In every time step, a new vertex is added to the graph. The new vertex chooses an ambassador (or more than one if ambs>1) and starts a simulated forest fire at its ambassador(s). The fire spreads through the edges. The spreading probability along an edge is given by fw_prob. The fire may also spread backwards on an edge by probability fw_prob * bw_factor. When the fire ended, the newly added vertex connects to the vertices ``burned'' in the previous fire.

Parameters:
n - the number of vertices in the graph
fw_prob - forward burning probability
bw_factor - ratio of backward and forward burning probability
ambs - number of ambassadors chosen in each step
directed - whether the graph will be directed
'''




'''
Establishment(n, k, type_dist, pref_matrix, directed=False)
source code 
Generates a graph based on a simple growing model with vertex types.

A single vertex is added at each time step. This new vertex tries to connect to k vertices in the graph. The probability that such a connection is realized depends on the types of the vertices involved.

Parameters:
n - the number of vertices in the graph
k - the number of connections tried in each step
type_dist - list giving the distribution of vertex types
pref_matrix - matrix (list of lists) giving the connection probabilities for different vertex types
directed - whether to generate a directed graph.
'''
# import scipy as sp
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import powerlaw
# numargs = powerlaw.numargs
# [ a ] = [0.9,] * numargs
# rv = powerlaw(a)
# x = np.linspace(0, np.minimum(rv.dist.b, 3))
# h = plt.plot(x, rv.pdf(x))
def testing():
    N=1
    n = 500
    for n in [500,1000,10000,100000]:#range(N):
    #    gFF = Graph.Forest_Fire(n, fw_prob=.36,bw_factor = .32,ambs=1)
    #     dd = [d/5 for d in gFF.degree()]
    #     print dd
        g = Graph.Barabasi(n)#), power=0, m=dd)
    #     g = gFF
    #     Graph.simplify(g)#remove self-loop and multiple edges 
        g.write_gml('BBtesting'+str(n)+'.gml')
    
    #     print(g)
    
    # http://igraph.sourceforge.net/doc/python/igraph.GraphBase-class.html#transitivity_undirected
        mem = Graph.community_fastgreedy(g.simplify(multiple=True, loops=True, combine_edges=None), None).as_clustering()
        # print mem 
        print 'n:',g.vcount(),' e:', g.ecount(),' cc:', g.transitivity_undirected(), ' q:', g.modularity(membership = mem, weights=None),
        print  ' maxD:' , g.maxdegree(),' ass:', g.assortativity_degree()
    # g.vs['color']= [cm.hot(cccc) for m in cccc for cccc in mem]  
    
        colors = ["lightgray", "cyan", "magenta", "yellow", "blue", "green", "red"]
        # print g.vs['color']
        components = mem
        for i, component in enumerate(components):
        #     print component
            for n in component:
                color = colors[min(6, i)]
                g.vs[n]["color"] = color
          
        layout = g.layout("kk")
        plot(g, layout = layout)
