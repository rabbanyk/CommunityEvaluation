import numpy as np
import matplotlib.pyplot as plt
from igraph import *

# read network 
def read_graph(name):
    g = Graph.Read_Edgelist(name, directed=False)
    return g 
# read ground-truth communities
def read_communities(name):
#     g = Graph.Read_Edgelist(name, directed=True)
#     return g.es
    res = [] 
    with open(name) as f:
        for line in f:
            i, j = [int(x) for x in line.split()]
#             print i, j
            while j>=len(res):
                res.append([])
#             print res, len(res)
            res[j].append(i)
    return res[1:]
# plot degree distribution for each community


def plot_dists(g_name,c_name ):
  
    g = read_graph(g_name)
    c = read_communities(c_name)
    summary(g)
    g.write_gml('test.gml', creator=None, ids=None)
    degrees= g.degree()#g.vs.degree()
    # print degrees
    
    
    dis,_ = np.histogram(degrees,   max(degrees))
    # print dis
    # plt.show()
    # plt.loglog(range(0,len(dis)),dis,'ro')
    # plt.show()
    
    # layout = g.layout("kk")
    # plot(g, g_name+".pdf", layout = layout)
    
    rows = 4
    tot = len(c)//2 
    cols =  tot//rows+1
    f, axss = plt.subplots(rows, cols,  sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)
    # print axss.shape
    axss[0,0].loglog(range(0,len(dis)),dis,'ro')
    # vs = np.array(g.vs.degree())
    for i,com in enumerate(c):
    # #     print com
    #     degrees = vs[com]
    # #     print degrees
    # #     plt.hist(degrees,   max(degrees))
    # #     plt.show()
    #     subV = g.vs.select(com)#g.vs.select(lambda vertex: vertex.index in com)
    #     subE = g.es.select(_within=com)
    #     print len(com), len(subV), len(subE)
    #     t = t+len(subE)
        subG = g.induced_subgraph(com, implementation="auto")
        
    #     layout = subG.layout("kk")
    #     plot(subG, g_name+"__"+str(i)+".pdf", layout = layout)
        degrees = subG.degree()
        dist, _ = np.histogram(degrees,   max(degrees))
    #     print i, i/2, i/2//cols, i/2%cols
        ax = axss[(i+2)/2//cols,(i+2)/2%cols]
        ax.loglog(range(0,len(dist)),dist,'o')
        ax.set_xlim([0, 10**3]) 
        ax.set_ylim([0, 10**3]) 
        
    # print len(g.es), t
    # plt.show()
    plt.savefig(g_name+'__dist.png')


for i in range(1,9):
    g_name = '/home/reihaneh/projects/AttributeGraphJava/execs/LFR/synthetic1000/S3/data/network'+str(i)+'.dat'
    c_name = '/home/reihaneh/projects/AttributeGraphJava/execs/LFR/synthetic1000/S3/data/community'+str(i)+'.dat'
    plot_dists(g_name, c_name)
