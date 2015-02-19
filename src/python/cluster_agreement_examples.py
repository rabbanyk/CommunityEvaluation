import collections
from cluster_agreement import  *
import pprint
from pylab import *
from matplotlib import rc
from tabulate import tabulate

def log_matrix(name, M):
    logging.info(name+': {0}'.format(np.array(M)))
   # logging.info(name+'A='+latex(M,mode='inline'))

def log_matrix_and_overlap(name, M):
    if M is None:
        return
    log_matrix(name, M)
    log_matrix(name+name+'^T', M.dot(M.T))
    log_matrix(name+'^T'+name, M.T.dot(M))
    
def log_clustered_graph(A,V,U1,U2):
    log_matrix_and_overlap('A', A)
    log_matrix_and_overlap('V', V)
    log_matrix_and_overlap('U1', U1)
    log_matrix_and_overlap('U2', U2)
   
    log_matrix('OVU', V.T.dot(U1))
    log_matrix('VV^T-UU^T', V.dot(V.T) - U1.dot(U1.T))
    if U2 is None:
        return
    log_matrix('OVU_2', V.T.dot(U2))
    log_matrix('VV^T - U_2U_2^T  ', V.dot(V.T) - U2.dot(U2.T))

def experiment(A,V,U,U2=None, fix_pos=None):  
    logging.basicConfig( level=logging.ERROR)
    log_clustered_graph(A,V,U,U2)

    G = nx.from_numpy_matrix(A)
    plotTable = False
    nr = 2 if plotTable else 1
    nf = 2 #numFigs
    if not U2==None:
        nf =3
   
    if fix_pos==None:
        fix_pos = nx.spring_layout(G)#, nlist, dim, scale)spring_layout(G)  
        print '----- positions to repeat the layout:: \n ',fix_pos
    
    figure(figsize=(9, 9))
    subplot(nr,nf,1)#.set_title('$U$')
    draw_clustered_graph(G,V, fix_pos)
    subplot(nr,nf,2)#.set_title('$V1$')
    draw_clustered_graph(G,U, fix_pos)
    if not U2 == None:
        subplot(nr,nf,3)#.set_title('$V2$')
        draw_clustered_graph(G,U2, fix_pos)
    show(block = False)
    

    Res = get_all_agreement_variations(V, U)
    headrs = np.hstack( ([" - "], Res[0]))
    vals = [ tuple(["(V,U_1)"]+Res[1].tolist()) ] 
    if U2 is not None:
        vals = vals + [ tuple(["(V,U_2)"]+get_all_agreement_variations(V, U2)[1].tolist()) ]   
    
    vals = vals + [ tuple(["(N,V)"]+agreement_with_structure(G, V)[1].tolist()) ]   
    vals = vals + [ tuple(["(N,U_1)"]+agreement_with_structure(G, U)[1].tolist()) ]   
    if not U2==None:
        vals = vals + [ tuple(["(N,U_2)"]+agreement_with_structure(G, U2)[1].tolist()) ]   

    for m in ["trans","sum"]:
        vals = vals + [ tuple(["(V,U_1|G)"+m[0]]+agreement_based_on_structure(G, V,U, method=m)[1].tolist()) ]   
        if not U2==None:
            vals = vals + [ tuple(["(V,U_2|G)"+m[0]]+agreement_based_on_structure(G, V,U2, method=m)[1].tolist()) ]   
    
    vals = np.asarray(vals, dtype=dict(names = headrs, formats=["a18"]+["float32"]*len(Res[1]) ))
    print tabulate(vals, headers = "keys", tablefmt="latex",  floatfmt=".3f")
    
    if plotTable:
        args = { "tablefmt":"plain",  "floatfmt":".3f", "numalign":"right"}
        ax = plt.subplot2grid((2,nf), (1,0), colspan=nf)
        plt.axis('off')
        plt.text(0,0, tabulate(vals, headers = "keys", **args))#,backgroundcolor = 'y',alpha=1)
    show(block = True)
    
def test_agreements():
    U = np.matrix([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    V = np.matrix([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    Res= get_all_agreement_variations(U, V)
    print tabulate(Res[1:], headers=Res[0],  tablefmt="grid") #latex

def test_clustered_graph_draw():
    A = np.array([[0,1,1,0,0,0,0,0,0,0],
            [1,0,0,1,0,0,0,0,0,0],
            [1,0,0,1,0,0,0,1,0,0],
            [0,1,1,0,1,1,0,0,0,0],
            [0,0,0,1,0,1,1,0,0,0],
            [0,0,0,1,1,0,1,0,0,0],
            [0,0,0,0,1,1,0,0,0,0],
            [0,0,1,0,0,0,0,0,1,1],
            [0,0,0,0,0,0,0,1,0,1],
            [0,0,0,0,0,0,0,1,1,0]])
            
    G = nx.from_numpy_matrix(A)
    U = np.array([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    figure(figsize=(10, 4))
    draw_clustered_graph(G,U, nx.spring_layout(G))
    show()
 
def basic_example():
    A = np.array([[0,1,1,0,0,0,0,0,0,0],
            [1,0,0,1,0,0,0,0,0,0],
            [1,0,0,1,0,0,0,1,0,0],
            [0,1,1,0,1,1,0,0,0,0],
            [0,0,0,1,0,1,1,0,0,0],
            [0,0,0,1,1,0,1,0,0,0],
            [0,0,0,0,1,1,0,0,0,0],
            [0,0,1,0,0,0,0,0,1,1],
            [0,0,0,0,0,0,0,1,0,1],
            [0,0,0,0,0,0,0,1,1,0]])
    fix_pos_GraphExample = None
#     fix_pos_GraphExample =  {0: array([ 0.39094446,  0.        ]), 
#                              1: array([ 0.23930109,  0.07809667]), 
#                              2: array([ 0.51977985,  0.14737463]), 
#                              3: array([ 0.26160767,  0.28061077]), 
#                              4: array([ 0.15234891,  0.4785181 ]), 
#                              5: array([ 0.06896513,  0.40358507]), 
#                              6: array([ 0.        ,  0.56391659]), 
#                              7: array([ 0.8040627 ,  0.10953033]), 
#                              8: array([ 0.98430573,  0.02177705]), 
#                              9: array([ 1.        ,  0.14796979])}
    fix_pos_GraphExample = {0: array([ 0.05,  .35]), 
                            1: array([ 0.02,  0.62]), 
                            2: array([ 0.25,  0.4]), 
                            3: array([ 0.25,  0.73]), 
                            4: array([ .5,  0.9]), 
                            5: array([ 0,  0.9]), 
                            6: array([ 0.25,  1 ]), 
                            7: array([ 0.25,  0.1]), 
                            8: array([ 0,  0.   ]), 
                            9: array([ .5, 0])}

    return A, fix_pos_GraphExample
def disjoint_exp():
    print '__________________________disjoint_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    V = np.array([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    U = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    print agreement_from_omega(U, V  )
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)

def disjoint_exp_extreme():
    V = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]])
    U = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]])
    sklearn_measures(U,V)

def overlap2_exp():
    print '__________________________disjoint_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[1,0,0],[1,0,0],[1,0,0],[1,1,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    V = np.array([[1,0,0],[1,0,0],[1,0,0],[.6,.4,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    print agreement_from_omega(U, V  )
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)

def overlap3_exp():
    print '__________________________disjoint_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[1,0,0],[1,0,0],[1,0,0],[1,1,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    V = np.array([[1,0,0],[1,0,0],[1,0,0],[.6,.4,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    V2 = np.array([[2,0,0],[2,0,0],[1,0,0],[1,1,0],[0,2,0],[0,2,0],[0,3,0],[0,0,2],[0,0,2],[0,0,2]])
    experiment(A,U,V,V2, fix_pos=fix_pos_GraphExample)
        
def disjoint_generalized_exp():
    print '__________________________disjoint_generalized_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[2,0],[2,0],[2,0],[0,2],[0,2],[0,2],[0,2],[2,0],[2,0],[2,0]])
    V = np.array([[2,0,0],[2,0,0],[2,0,0],[2,0,0],[0,2,0],[0,2,0],[0,2,0],[0,0,2],[0,0,2],[0,0,2]])
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)
  
def disjoint_generalized_exp_2():
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[3,0],[3,0],[3,0],[0,3],[0,3],[0,3],[0,3],[3,0],[3,0],[3,0]])
    V = np.array([[2,0,0],[2,0,0],[2,0,0],[2,0,0],[0,2,0],[0,2,0],[0,2,0],[0,0,2],[0,0,2],[0,0,2]])
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)

def disjoint_fuzzy_exp():
    print '__________________________disjoint_fuzzy_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[0.2,0],[0.2,0],[0.2,0],[0,0.2],[0,0.2],[0,0.2],[0,0.2],[0.2,0],[0.2,0],[0.2,0]])
    V = np.array([[0.2,0,0],[0.2,0,0],[0.2,0,0],[0.2,0,0],[0,0.2,0],[0,0.2,0],[0,0.2,0],[0,0,0.2],[0,0,0.2],[0,0,0.2]])
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)
  
def overlap_exp():
    print '__________________________overlap_hard_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[1,0],[1,0],[1,0],[1,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    V = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)


def fuzzy_exp():
    print '__________________________overlap_fuzzy_exp________________________________'
    A, fix_pos_GraphExample = basic_example()
    U = np.array([[1,0],[1,0],[1,0],[.5,.5],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    V = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    experiment(A,U,V, fix_pos=fix_pos_GraphExample)

def graph_example():
    A = np.zeros((9,9),dtype='float')
    for e in [(0,1),(0,5),(0,6),
              (1,0),(1,2),(1,5),
              (2,1),(2,3),(2,5),
              (3,2),(3,5),(3,4),
              (4,3),(4,5),
              (5,0),(5,1),(5,2),(5,3),(5,4),(5,6),(5,8),
              (6,0),(6,5),(6,7),(6,8),
              (7,6),(7,8),
              (8,6),(8,7),(8,5),
              ]:
        A[e]=1

    fix_pos_GraphExample = {0: np.array([ 0.20385061-0.1,  0.23977043-0.1]), 
               1: np.array([ 0.43693314,  0.        ]), 
               2: np.array([ 0.75900616,  0.01462725]), 
               3: np.array([ 0.97478493,  0.24577521-0.1]), 
               4: np.array([ 1.        ,  0.53455287-0.1-0.1]), 
               5: np.array([ 0.56256428+0.05,  0.38383768 +0.1]), 
               6: np.array([ 0.16250501-0.1,  0.59711573]), 
               7: np.array([ 0.        ,  0.96007752]), 
               8: np.array([ 0.32060991,  0.8106598 ])}


    U = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1]])
    V1 = np.array([[0,1],[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1]])
    V2 = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1]])

    experiment(A,U,V1) #,V2,  fix_pos=fix_pos_GraphExample)
    print convert_cluster_to_labels(U)
    print convert_cluster_to_labels(V1)
    print convert_cluster_to_labels(V2)
    
def omega_example():
    A = np.zeros((5,5),dtype='float')
    for e in [(0,1),(0,4),
              (1,0),(1,2),
              (2,1),(2,3),
              (3,2),(3,4),
              (4,3),(4,0),
              ]:
        A[e]=1

    fix_pos_GraphExample =   {0: array([ 0.49682607,  0.95151713]), 
                              1: array([ 1.        ,  0.59183855]), 
                              2: array([  8.12052312e-01,   3.19389103e-05]), 
                              3: array([ 0.19258776,  0.        ]), 
                              4: array([ 0.        ,  0.58628989])}


    U = np.array([[1,0,0],[0,1,0],[0,1,1],[1,1,1],[1,0,1]])
    V1 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,0,0],[1,0,0]])
    V2 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,0,1],[1,0,0]])
    print 'V_1', agreement_from_omega(U, V1)
    print U.T.dot(V1)
    print 'V_2', agreement_from_omega(U, V2)
    print U.T.dot(V2)
#     return
    experiment(A,U,V1,V2,  fix_pos=fix_pos_GraphExample)

def convert_cluster_to_labels(U):
    return np.nonzero(U)[1]

def sklearn_measures(U, V):
    #     http://scikit-learn.org/stable/modules/classes.html#clustering-metrics
    import sklearn.metrics.cluster as sym
    U_labels = np.nonzero(U)[1]
    V_labels = np.nonzero(V)[1]
    print U_labels, V_labels
#     V2_labels = np.nonzero(V2)[1]
    print 'entro(U)=',sym.entropy(U_labels),'entro(V)=',sym.entropy(V_labels), 'entro(U,V)=',sym.mutual_info_score(U_labels, V_labels)
    res = [ ['ari', 'nmi', 'ami', 'vm' ], \
            [ sym.adjusted_rand_score(U_labels, V_labels),\
              sym.normalized_mutual_info_score(U_labels, V_labels),\
              sym.adjusted_mutual_info_score(U_labels, V_labels),\
              sym.v_measure_score(U_labels, V_labels)]]
    print res
    return res

def matching_exp():
    A = np.zeros((11,11),dtype='float')
    for e in [(0,1),(0,3),(0,4),
              (1,2),(1,3),(1,4),
              (2,3),
              (3,4),
              (5,6),(5,7),
              (6,7),
              (7,8),
              (8,9),(8,10),
              (9,10)
              ]:
        A[e]=A[(e[1],e[0])]=1

    fix_pos_GraphExample =   {0: array([ 0.5/2 ,  0.95]), 
                              1: array([ 1./2  ,  0.59]), 
                              2: array([ 0.81/2,  0.  ]), 
                              3: array([ 0.19/2,  0.  ]), 
                              4: array([ 0.  ,  0.59]), 
                              5: array([ 0.74,  1.  ]), 
                              6: array([ 1,  1.  ]), 
                              7: array([ 0.88,  0.69]), 
                              8: array([ 0.88,  0.35]), 
                              9: array([ 0.74,  0]), 
                              10: array([ 1,  0])}

    np.set_printoptions(precision=3,suppress=True)
    print fix_pos_GraphExample

    U = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],\
                   [0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])

    V1 = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],\
                   [0,1],[0,1],[0,1],[0,1],[0,1],[0,1]])
    
    V2 = np.array([[1,0,0],[0,0,1],[0,0,1],[1,0,0],[1,0,0],\
                   [0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    
    print U
    print V1
    print U.T.dot(V1)
    M1= sklearn_measures(U, V1)
    M2= sklearn_measures(U, V2)
    
    for i,m in enumerate(M1[0]):
        print m,'::', M1[1][i], ' ', M2[1][i] 
    
    print '---', nmi(U,V1)
    print '---', nmi(U,V2)
    
    experiment(A,U,V1,V2,  fix_pos=fix_pos_GraphExample)

#graph_example()

#

def wighted(directed =False):
    A = np.array([[ 0. , 0. , 9. , 0.],
                     [ 4.,  0.,  0.,  0.],
                     [ 0. , 1. , 0. , 16.],
                     [ 1.  ,0. , 0.,  0.]])
    if not directed: A = A+A.T
#     printlatex(A)
    print  '---- A: ---- \n',A
    
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph() if directed else None) 

#     pos =  nx.circular_layout(G)#spectral_layout(G))#hell_layout(G) )
    pos = {0: array([ 1. ,  0.5], dtype=float32), 
           1: array([ 0.49999997,  1.        ], dtype=float32), 
           2: array([ 0.        ,  0.49999997], dtype=float32), 
           3: array([ 0.5,  0. ], dtype=float32)}
    figure(figsize=(4, 4))
    draw_clustered_graph(G, None,pos)
    show(block=False)
    
 
    N = get_incident_matrix_from_adjecency(A,directed)
#     printlatex(N)
    print '---- N: ---- \n',N
    L2 = N.dot(abs(N.T))
#     printlatex(L2)
    print  '---- NNT: ---- \n',L2
    D = sum(A, axis=0)
    D = np.diag(D)
    print  '---- D: ---- \n',D
    print  '---- D+A: ---- \n',D+A
    show()
  
def wighted_directed():
    A = np.array([[ 0. , 0. , 9. , 0.],
                 [ 4.,  0.,  0.,  0.],
                 [ 0. , 1. , 0. , 16.],
                 [ 1.  ,0. , 0.,  0.]])
    
    printlatex(A)
    print  '---- A: ---- \n',A
#     Bout =np.array([[sqrt(2),0,0,1,0], [0,1,0,0,0],[0,0,sqrt(3),0,sqrt(4)],[0,0,0,0,0]])
#     Bin =np.array([[0,0,sqrt(3),0,0], [sqrt(2),0,0,0,0],[0,1,0,0,0],[0,0,0,1,sqrt(4)]])
    G = nx.from_numpy_matrix(A,create_using=nx.DiGraph())

#     pos =  nx.circular_layout(G)#spectral_layout(G))#hell_layout(G) )
    pos = {0: array([ 1. ,  0.5], dtype=float32), 
           1: array([ 0.49999997,  1.        ], dtype=float32), 
           2: array([ 0.        ,  0.49999997], dtype=float32), 
           3: array([ 0.5,  0. ], dtype=float32)}
    figure(figsize=(4, 4))
    draw_clustered_graph(G, None,pos)
    show(block=False)
    
    edges =  transpose(np.nonzero(A))
#     print edges
    Bout = np.zeros((A.shape[0], edges.shape[0]), dtype =np.float)
    Bin = np.zeros((A.shape[0], edges.shape[0]), dtype =np.float)
    
    for j,e in enumerate(edges):
#         print e[0] --> e[1], A[e[0],e[1]]
        Bout[e[0]][j] = sqrt(A[e[0],e[1]])
        Bin[e[1]][j] = sqrt(A[e[0],e[1]])
  
    print '---- Bout: ---- \n',Bout
    print '---- Bin: ---- \n',Bin 
    
    BoBi = Bout.dot(Bin.T) 
    print 'A = BoBi.T\n',BoBi
    #print nx.adjacency_matrix(G)
#     print G.edges(data=True)
    Di = Bin.dot(Bin.T) 
    Do = Bout.dot(Bout.T) 
    print 'In-degree matrix = BiBi.T\n',Di
    print 'Out-degree matrix = BoBo.T\n',Do
    D = Di+Do
    print 'In+Out-degree matrix = BiBi.T+BoBo.T = D\n',D

    M = Bout-Bin
    print 'M: incidence matrix = Bo-Bi \n' ,M
    L= M.dot(M.T)
    print 'L: laplacian matrix of undirected version= MM.T \n', L
    
    N = Bout+Bin
    print 'N: incidence matrix 2 = Bo+Bi \n' ,N
    L2 = N.dot(N.T)
    print 'NNT=\n',L2
    print 'NNT-D=\n',L2-D
   
    print 'NBi=\n',N.dot(Bin.T)
    print  '---- A: ---- \n',A
    print 'M|M|T=\n',M.dot(np.abs(M.T)).clip(min=0)
    print 'M|M|T=\n',M.dot(np.abs(M.T))+Di 
    print 'M|M|T=\n',M.dot(np.abs(M.T))+Di -Do

    show()
#     print 'adjacency of Line Graph, edges as nodes, nodes as egdes: \n',N.T.dot(N)


 
def wighted_directed_complex_try():
    A = np.array([[ 0. , 0. , 9. , 0.],
                 [ 4.,  0.,  0.,  0.],
                 [ 0. , 1. , 0. , 16.],
                 [ 1.  ,0. , 0.,  0.]])
    
    printlatex(A)
    print  '---- A: ---- \n',A
#     Bout =np.array([[sqrt(2),0,0,1,0], [0,1,0,0,0],[0,0,sqrt(3),0,sqrt(4)],[0,0,0,0,0]])
#     Bin =np.array([[0,0,sqrt(3),0,0], [sqrt(2),0,0,0,0],[0,1,0,0,0],[0,0,0,1,sqrt(4)]])
    G = nx.from_numpy_matrix(A,create_using=nx.DiGraph())

#     pos =  nx.circular_layout(G)#spectral_layout(G))#hell_layout(G) )
    pos = {0: array([ 1. ,  0.5], dtype=float32), 
           1: array([ 0.49999997,  1.        ], dtype=float32), 
           2: array([ 0.        ,  0.49999997], dtype=float32), 
           3: array([ 0.5,  0. ], dtype=float32)}
    figure(figsize=(4, 4))
    draw_clustered_graph(G, None,pos)
    show(block=False)
    
    edges =  transpose(np.nonzero(A))
#     print edges
    Bout = np.zeros((A.shape[0], edges.shape[0]), dtype =np.complex)
    Bin = np.zeros((A.shape[0], edges.shape[0]), dtype =np.complex)
    
    for j,e in enumerate(edges):
#         print e[0] --> e[1], A[e[0],e[1]]
        Bout[e[0]][j] = sqrt(A[e[0],e[1]]) +1j
        Bin[e[1]][j] = sqrt(A[e[0],e[1]])-1j
  
    print '---- Bout: ---- \n',Bout
    print '---- Bin: ---- \n',Bin 
    
    BoBi = Bout.dot(Bin.T) 
    print 'A = BoBi.T\n',BoBi
    #print nx.adjacency_matrix(G)
#     print G.edges(data=True)
    Di = Bin.dot(Bin.T) 
    Do = Bout.dot(Bout.T) 
    print 'In-degree matrix = BiBi.T\n',Di
    print 'Out-degree matrix = BoBo.T\n',Do
    D = Di+Do
    print 'In+Out-degree matrix = BiBi.T+BoBo.T = D\n',D

    M = Bout-Bin
    print 'M: incidence matrix = Bo-Bi \n' ,M
    L= M.dot(M.T)
    print 'L: laplacian matrix of undirected version= MM.T \n', L
    
    N = Bout+Bin
    print 'N: incidence matrix 2 = Bo+Bi \n' ,N
    L2 = N.dot(N.T)
    print 'NNT=\n',L2
    print 'NNT-D=\n',L2-D
   
    print 'NBi=\n',N.dot(Bin.T)
    print  '---- A: ---- \n',A
    print 'M|M|T=\n',M.dot(np.abs(M.T)).clip(min=0)
    print 'M|M|T=\n',M.dot(np.abs(M.T))+Di 
    show()
#     print 'adjacency of Line Graph, edges as nodes, nodes as egdes: \n',N.T.dot(N)



def transformation():
    A = np.zeros((9,9),dtype='float')
    edges = [(0,1),(0,5),(0,6),
              (1,0),(1,2),(1,5),
              (2,1),(2,3),(2,5),
              (3,2),(3,5),(3,4),
              (4,3),(4,5),
              (5,0),(5,1),(5,2),(5,3),(5,4),(5,6),(5,8),
              (6,0),(6,5),(6,7),(6,8),
              (7,6),(7,8),
              (8,6),(8,7),(8,5),
              ]
    for e in edges:
        A[e]=1

    fix_pos = {0: np.array([ 0.20385061-0.1,  0.23977043-0.1]), 
               1: np.array([ 0.43693314,  0.        ]), 
               2: np.array([ 0.75900616,  0.01462725]), 
               3: np.array([ 0.97478493,  0.24577521-0.1]), 
               4: np.array([ 1.        ,  0.53455287-0.1-0.1]), 
               5: np.array([ 0.56256428+0.05,  0.38383768 +0.1]), 
               6: np.array([ 0.16250501-0.1,  0.59711573]), 
               7: np.array([ 0.        ,  0.96007752]), 
               8: np.array([ 0.32060991,  0.8106598 ])}

    edges =  transpose(np.nonzero(np.tril(A)))
    edg_pos ={}
    for j,e in enumerate(edges):
        print j, e[0], e[1]
        edg_pos[j]= 0.5*(fix_pos[e[0]]+fix_pos[e[1]])
    print edg_pos
    U = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1]])
    V1 = np.array([[0,1],[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1]])
    V2 = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1]])
   
    N = get_incident_matrix(A.shape[0],edges)
    
    G=nx.from_numpy_matrix(A)#=N.dot(N.T)
    LG=nx.from_numpy_matrix(N.T.dot(N)) #Line graph
    
    print N.shape, U.shape, A.shape,( N.T.dot(N)).shape
    T = N.T.dot(U)
    TRA = T.dot(T.T)
    print TRA.shape, A.shape
#     nx.write_dot(LG,'graph.dot')
    
    subplot(2,4,1)
#     draw_clustered_graph(G ,None, pos= fix_pos)
    draw_clustered_graph(G ,U, pos= fix_pos,draw_graph = True)
    subplot(2,4,2)    
    draw_clustered_graph(G ,V1, pos= fix_pos,draw_graph = True)
    subplot(2,4,3)
    draw_clustered_graph(G ,V2, pos= fix_pos,draw_graph = True)
#     subplot(2,4,4)
    printlatex(N.T)
    printlatex(N.T.dot(U))
    printlatex(N.T.dot(V1))
    printlatex(N.T.dot(V2))


#     subplot(2,4,5)
#     draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(2,4,6)    
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(U)).dot(U.T.dot(N))) ,N.T.dot(U), pos= edg_pos,draw_graph = False)
    subplot(2,4,7)
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(V1)).dot(V1.T.dot(N))) ,N.T.dot(V1), pos= edg_pos,draw_graph = False)
    subplot(2,4,8)
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(V2)).dot(V2.T.dot(N))) ,N.T.dot(V2), pos= edg_pos,draw_graph = False)
    show()
    return 

    subplot(5,4,1)
    draw_clustered_graph(LG ,None, pos= edg_pos,draw_graph = False)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,2)    
    draw_clustered_graph(LG ,N.T.dot(U), pos= edg_pos,draw_graph = False)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,3)
    draw_clustered_graph(LG ,N.T.dot(V1), pos= edg_pos,draw_graph = False)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,4)
    draw_clustered_graph(LG ,N.T.dot(V2), pos= edg_pos,draw_graph = False)
    draw_clustered_graph(G ,None, pos= fix_pos)
    
    subplot(5,4,5)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,6)    
    draw_clustered_graph(nx.from_numpy_matrix(U.dot(U.T)) ,U, pos= fix_pos)
    subplot(5,4,7)
    draw_clustered_graph(nx.from_numpy_matrix(V1.dot(V1.T)) ,V1, pos= fix_pos)
    subplot(5,4,8)
    draw_clustered_graph(nx.from_numpy_matrix(V2.dot(V2.T)) ,V2, pos= fix_pos)
    
    subplot(5,4,9)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,10)    
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(U)).dot(U.T.dot(N))) ,N.T.dot(U), pos= edg_pos)
    subplot(5,4,11)
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(V1)).dot(V1.T.dot(N))) ,N.T.dot(V1), pos= edg_pos)
    subplot(5,4,12)
    draw_clustered_graph(nx.from_numpy_matrix((N.T.dot(V2)).dot(V2.T.dot(N))) ,N.T.dot(V2), pos= edg_pos)
   
    
    subplot(5,4,13)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,14)    
    draw_clustered_graph(nx.from_numpy_matrix((N.dot(N.T)).dot(U.dot(U.T))) ,U, pos= fix_pos)
    subplot(5,4,15)
    draw_clustered_graph(nx.from_numpy_matrix((N.dot(N.T)).dot(V1.dot(V1.T))) ,V1, pos= fix_pos)
    subplot(5,4,16)
    draw_clustered_graph(nx.from_numpy_matrix((N.dot(N.T)).dot(V2.dot(V2.T))) ,V2, pos= fix_pos)
   
    subplot(5,4,17)
    draw_clustered_graph(G ,None, pos= fix_pos)
    subplot(5,4,18)    
    draw_clustered_graph(nx.from_numpy_matrix((A.T.dot(U)).dot(U.T.dot(A))) ,A.T.dot(U), pos= fix_pos)
    subplot(5,4,19)
    draw_clustered_graph(nx.from_numpy_matrix((A.T.dot(V1)).dot(V1.T.dot(A))) ,A.T.dot(V1), pos= fix_pos)
    subplot(5,4,20)
    draw_clustered_graph(nx.from_numpy_matrix((A.T.dot(V2)).dot(V2.T.dot(A))) ,A.T.dot(V2), pos= fix_pos)
   
    show()

def main():
    overlap3_exp()
#     disjoint_exp()
#      disjoint_exp_extreme()

#     overlap2_exp()
#     matching_exp()
#     omega_example()
#     wighted()
#     wighted(directed=True)
#     transformation()
#     graph_example()
#     disjoint_generalized_exp()
#  
#     disjoint_generalized_exp_2()
#     disjoint_fuzzy_exp()
#     overlap_exp()
#     fuzzy_exp()


main()