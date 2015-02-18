from igraph import *


''' 
,community_edge_betweenness           
,community_leading_eigenvector ---        
,community_spinglass ---
,community_fastgreedy                 
,community_leading_eigenvector_naive  
,community_walktrap
,community_infomap                    
,community_multilevel                 
,community_label_propagation          
,community_optimal_modularity     

usage:    
http://igraph.sourceforge.net/screenshots2.html#4
http://igraph.sourceforge.net/igraphbook/igraphbook-random.html
http://www.cs.rhul.ac.uk/home/tamas/development/igraph/tutorial/tutorial.html
http://igraph.sourceforge.net/doc-0.6/python/igraph.Graph-class.html#community_optimal_modularity
http://bommaritollc.com/2012/06/17/summary-community-detection-algorithms-igraph-0-6/
'''
clusters = None
directed=False 
weights=None
vertex_weights = None

def communities(g, method):
    if (method ==  'edge_betweenness'):
        print 'community_edge_betweenness'
        return (Graph.community_edge_betweenness(g,clusters ,  directed, weights))
    elif (method ==  'leading_eigenvector'):
        print 'community_leading_eigenvector' 
        return (Graph.community_leading_eigenvector(g,clusters))    
    elif (method ==  'spinglass'):
        print 'community_spinglass' 
        return (Graph.community_spinglass(g))
    elif (method ==  'walktrap'):
        print 'community_walktrap' 
        return (Graph.community_walktrap(g, weights, steps=4).as_clustering()) 
    elif (method ==  'fastgreedy'):
        print 'community_fastgreedy' 
        return (Graph.community_fastgreedy(g.simplify(multiple=True, loops=True, combine_edges=None), weights).as_clustering())   
    elif (method ==  'infomap'):
        print 'community_infomap' 
        return (Graph.community_infomap(g, weights, vertex_weights, trials=10))   
    elif (method ==  'multilevel'):
        print 'community_multilevel' 
        return (Graph.community_multilevel(g, weights, return_levels=False))                  
    elif (method ==  'label_propagation'):
        print 'community_label_propagation' 
        return (Graph.community_label_propagation(g, weights, initial=None, fixed=None))           
    elif (method ==  'optimal_modularity'):
        print 'community_optimal_modularity' 
        return (Graph.community_optimal_modularity(g))   

def findCommunities(method, inputPath, outputPath = 'res.com', format = 'gml'):
#     print inputPath, outputPath
    g = load(inputPath, format)

    f =None
    if outputPath<>None:
        f = open(outputPath, 'w')
    
    comps= Graph.clusters(g, mode="weak")
    print f
    
    for comp in comps:
        subg =  Graph.subgraph(g,comp)
        if len(comp)>10:
            comm = communities(subg, method) 
        else:
            comm = [range(len(comp))]
        if outputPath<>None:
            for c in comm:
                f.write(str( [comp[i] for i in c]) + "\n")
        else:  
            print list(comm)
    if outputPath<>None:
        f.close()

#method - the measure to use. "vi" or "meila" means the variation of information metric of Meila (2003), "nmi" or "danon" means the normalized mutual information as defined by Danon et al (2005), "split-join" means the split-join distance of van Dongen (2000), "rand" means the Rand index of Rand (1971), "adjusted_rand" means the adjusted Rand index of Hubert and Arabie (1985).
#compare_communities(comm1, comm2, method='vi', remove_none=False)

import sys, getopt

def main(argv):
    inputfile = None
    outputfile = None
    method = None
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:",["ifile=","ofile=","method"])
    except getopt.GetoptError:
        print 'communityMinerInterface.py -i <inputfile> -o <outputfile> -m <method>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'communityMinerInterface.py -i <inputfile> -o <outputfile> -m <method>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--method"):
            method = arg
    if method==None or inputfile ==None  :
        print 'usage: communityMinerInterface.py -i <inputfile> -o <outputfile> -m <method>'
        sys.exit() 
    else:
        findCommunities(method, inputfile, outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])

# findCommunities('../../../wkarate.net', 'community_infomap' )

# 
# def plotDegreeDist(g):
#     d =Graph.degree(g, mode="in")
#     dd =Graph.degree_distribution(g, mode="in")#, cumulative=True)
#     fpl =  statistics.power_law_fit(d, xmin=20, return_alpha_only=False)
#     print fpl    
#     plot(fpl)#, log="xy", xlab="degree", ylab="cumulative frequency", col=1")
# #     lines(10:500, 10*(10:500)^(-coef(alpha)+1))
#     
# 
# 
# plotDegreeDist(Graph.Barabasi(100000))




        