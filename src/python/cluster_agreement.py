import math
import logging
import numpy as np
import numpy.linalg as LA
import networkx as nx
import igraph as ig
import networkx.linalg.graphmatrix as gm
from sympy.printing import latex
from sympy.matrices import *

FORMAT = "[%(lineno)s:%(funcName)20s()]\n %(message)s"
logging.basicConfig(format=FORMAT, level=logging.ERROR)
#logging.basicConfig(filename='run.log',level=logging.DEBUG)

def printlatex(A):
    print latex(Matrix(A),mode='inline')

def draw_clustered_graph(G,U, pos, draw_edge_wights = False, draw_graph = True):
    """
    Visualizes the given clustering of a graph.
    
    Parameters
    ----------
    G: A networkx graph 
    U: numpy matrix
        A nxk matrix representing a clustering, where n is the number of data points and k is the number of clusters in U, 
        so that U_{ik} shows if the ith data-point belongs to the kth cluster in U.   
    
    pos: A Python dictionary (optional)
        A dictionary of positions for graph nodes keyed by node

        
    Returns
    -------
    None
    
    Examples
    -------
    >>> import numpy as np
    >>> import networkx as nx
    >>> from pylab import *
    
   >>> np.matrix(  [[0,1,1,0,0,0,0,0,0,0],
                    [1,0,0,1,0,0,0,0,0,0],
                    [1,0,0,1,0,0,0,1,0,0],
                    [0,1,1,0,1,1,0,0,0,0],
                    [0,0,0,1,0,1,1,0,0,0],
                    [0,0,0,1,1,0,1,0,0,0],
                    [0,0,0,0,1,1,0,0,0,0],
                    [0,0,1,0,0,0,0,0,1,1],
                    [0,0,0,0,0,0,0,1,0,1],
                    [0,0,0,0,0,0,0,1,1,0]])
            
    >>> G = nx.from_numpy_matrix(A)
    >>> U = np.matrix([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    >>> figure(figsize=(10, 4))
    >>> draw_clustered_graph(G,U, nx.spring_layout(G))
    >>> show()
            
    """
    
    color_list = ['b','g','r','m','y','c']
    color_list =['#9999FF', '#FF9999', '#66FF99','#FF99FF', '#FFFF99', '#E0E0E0','#FF9933','#FF3333','#33FF99','#33FFFF','#3333FF']

    if pos==None:
        pos = nx.spring_layout(G)#, nlist, dim, scale)spring_layout(G)  
        print '----- positions to repeat the layout:: \n ',pos
    
    maxW = 0.0
    for (i,j, w) in G.edges(data=True):
        maxW = max(maxW,w['weight'] )
    if U is not None:
        UG = nx.from_numpy_matrix( U.dot(U.T))
        maxU = U.max()
        # draw clusters
        t =6.2
        for c in range(U.shape[1]):
            nsizes = [ int(t**4 * U[i,c]/maxU) for i in range(U.shape[0])]
            esizes = [ ((U[i,c]*U[j,c])*t**2.0/maxU**2) for (i,j, w) in UG.edges(data=True)]
            
            try:
                nx.draw(UG, pos = pos, node_size=nsizes, style ='solid' ,\
                        #labels ='.',\
                        node_color = color_list[c],edge_color =color_list[c],\
                        linewidths=0, width= esizes, alpha =.5)#, alpha =1
            except:
                pass
    # draw graph
#     if draw_graph:
    try:
        f= lambda w: np.log(w/maxW +1)*w*3/maxW +1  if draw_graph else 0 
        nx.draw(G, pos = pos, node_color = 'w', alpha =1 ,  linewidths=1, width= [f(w['weight']) for (i,j, w) in G.edges(data=True)])
    except:
        pass
    
    if draw_edge_wights:
        edge_lbls = dict([((u,v,),d['weight']) for u,v,d in G.edges(data=True)])
    #     print edge_lbls
        nx.draw_networkx_edge_labels(G, pos = pos, label_pos=0.5 ,edge_labels =edge_lbls)
   
    return pos
    
def agreement_from_overlaps(U,V, f = lambda x: x*(x-1)*0.5, exact=False):
    """ 
    Computes the agreement between two disjoint clusterings based on the overlap of their clusters.
    Parameters
    ----------
    U: numpy matrix
        A nxk matrix representing a clustering, 
        where n is the number of data points and k is the number of clusters in U, 
        so that U_{ik} shows if the ith data-point belongs to the kth cluster in U. 
    V: numpy matrix
     #    return { '0_RCM(N,U)':agreement_from_co_memberships(N,U), 
#             '1_A(N,U)': agreement_alternative_forms(N,U),
#             
#             'RCM(A,U)':agreement_from_co_memberships(A,U),
#             'A(A,U)': agreement_alternative_forms(A,U)}
       A nxr matrix representing a clustering similar to U
    f: function
        A Python function (non-linear) to apply on the overlap sizes to compute the divergence,
        f = lambda x: x*(x-1)*0.5 (default) will result in Rand Index and Adjusted Rand Index, 
        f = lambda x: 0 if x==0 else (x * math.log(x)) derives  normalized VI.
  
    Returns
    -------
        G: float
            agreement between input clusterings from 
            $$\GAM_{\varphi} = \frac{\mathbf{1} \varphi(N) \mathbf{1}^T  - \mathbf{1} \varphi(N \mathbf{1}^T )   } {\varphi(\mathbf{1} N \mathbf{1}^T) } + \frac{\mathbf{1} \varphi(N) \mathbf{1}^T - \varphi(\mathbf{1} N) \mathbf{1}^T } {\varphi(\mathbf{1} N \mathbf{1}^T) }  $$       
        AG: float
            agreement between input clusterings from 
            $$\AG_{\varphi} =\frac{\mathbf{1}  \varphi(N)  \mathbf{1}^T  -  E}{M-E} ,\quad M = \frac{ \mathbf{1} \varphi(N \mathbf{1}^T ) +  \varphi(\mathbf{1} N) \mathbf{1}^T  }{2 } ,\; E = \frac{ (\mathbf{1}  \varphi(N  \mathbf{1}^T ))  ( \varphi(\mathbf{1} N ) \mathbf{1}^T  )}{\varphi(\mathbf{1}  N \mathbf{1}^T) } $$

    Examples
    --------
    >>> U = Matrix([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    >>> V = Matrix([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    >>> agreement_from_overlaps(U,V)
    """
    return 1- distance_from_overlaps(U, V, f, exact)

def distance_from_overlaps(U,V, f = lambda x: x*(x-1)*0.5, exact=False):
    f = np.vectorize(f, otypes=[np.float])
    N = np.dot(U.T,V) # the overlaps, a.k.a contingency table
   
    logging.info('the overlap matrix: \n {0}'.format(np.matrix(N)))
    logging.debug('the overlap matrix transposed: \n {0}'.format(np.matrix(np.dot(V.T,U))))
  #  logging.info('$ N = $'+latex(N,mode='inline'))
    
    b2 = np.ones((N.shape[1],1),dtype='float')
    b1 = np.ones((1, N.shape[0]),dtype='float')
    
    m1=b1.dot(N)
    m2=N.dot(b2)
    m = b1.dot(N).dot(b2)
    s1 = (b1.dot(f(m2))) 
    s2 = (f(m1).dot(b2)) 
    n = f(m)
    I = (b1.dot(f(N)).dot(b2))
    
    # or 
    # s1 = f(N.sum(axis=1)).sum()
    # s2 = f(N.sum(axis=0)).sum()
    # n =  f(N.sum())
    # I =  f(N).sum()
    
    logging.debug('s1: {0}, s2: {1}, n: {2}, I: {3}'.format(s1,s2,n,I))
#     if (s1- I)/n>0.5 or (s2- I)/n >0.5 :  
#         print '---------------------------', (s1- I)/n, (s2- I)/n  
#     if (2*I) >n :  
#         print '--======================================', (2*I)/n    
    G =  float ((s1+s2- 2* I  )/( n ) )
#     if  float ((s1+s2)/( n ) <1 ):
#         print '--------------------------------------------------------------------------', G<1   

#     print m1.shape, m2.shape, m2.dot(m1)
#     print  b1.shape, b2.shape
    if exact:
        E = (b1.dot(f(m2).dot(f(m1))/f(m))).dot(b2)
    else:
        E = (b1.dot(f(m2.dot(m1)/m))).dot(b2)
#     EA = (s1*s2)/n
#     EA = (b1.dot(f(m2).dot(f(m1))/f(m))).dot(b2)
    
#     print '-----------', 1-G, s1+s2, I, E , EA , EM.sum(),'==',m
#     print '--', s1+s2-2*I, s1+s2-2*E, s1+s2-2*EA, ( s1+s2-2*I)/(s1+s2-2*E),( s1+s2-2*I)/(s1+s2-2*EA) 

    AG = float ((s1+s2- 2* I  )/( (s1+s2) - 2 *E ))
#     AG =  float ((s1+s2- 2* I  )/( (s1+s2) - 2*(s1*s2)/n ))
  #  AG = 1- float ((I  - (s1*s2)/n )/( (s1+s2)/2 - (s1*s2)/n ))
#     print '------ E = ',E, ', e/m-logm = ', E/m - math.log(m), '     ', 1-AG
#     print AG
    logging.debug( '-->'+ 'G= {0}, AG = {1}'.format( G,AG))
    return  np.array([G, AG])


def nmi(U,V):
    N = np.dot(U.T,V) # the overlaps, a.k.a contingency table
   
    b2 = np.ones((N.shape[1],1),dtype='float')
    b1 = np.ones((1, N.shape[0]),dtype='float')
    
    n = b1.dot(N).dot(b2)
    nv= b1.dot(N)
    nu= N.dot(b2)
    
    f = np.vectorize(lambda x: x*math.log(x) if x >0 else 0, otypes=[np.float])
    Hu = -(b1.dot(f(nu/n))) 
    Hv = -(f(nv/n).dot(b2)) 
    Huv = -(b1.dot(f(N/n)).dot(b2))
    Iuv = Hu+Hv-Huv
#     print 'H(U)=',Hu,'H(V)=',Hv,'I(U,V)=',Iuv 
    print '-------hu+hv-2huv = ', Hu+Hv-2*Huv, '----hu+hv = ', Hu+Hv, '-- huv =', Huv  

    vi = float((2*Huv-(Hu+Hv) )/math.log(n))
    nmi_sum = float(2*Iuv/(Hu+Hv)) # VM in sklearn
    nmi_sqrt = float(Iuv/math.sqrt(Hu*Hv)) # NMI in sklearn

#     G =  float ((Hu+Hv- 2* Huv  )/( n ) )
#     E = (b1.dot(f(nu.dot(nv)/n))).dot(b2)
#     AG = float ((Hu+Hv- 2* Huv  )/( (Hu+Hv) - 2 *E ))
    return  np.array([vi, nmi_sum, nmi_sqrt])


def agreement_from_co_memberships(U, V, same_nodes = True):
    return 1-distance_from_co_memberships(U, V, same_nodes)

def distance_from_co_memberships(U, V, same_nodes = True):

    """ 
    Computes the agreement between two clusterings based on the co_memberships of data points in their clusters.
    Parameters
    ----------
    U: Matrix
        A nxk matrix representing a clustering, 
        where upper is the number of data points and k is the number of clusters in U, 
        so that U_{ik} shows if the ith data-point belongs to the kth cluster in U. 
    V: Matrix
        A nxr matrix representing a clustering similar to U
 
    original: Boolean
        determines whether the formula should be exact same as the original RI formula or the approximate forms
        
    Returns
    -------
        G: float
            agreement between input clusterings from 
            $$\GAM_{\varphi} = 1 - \frac{|| UU^T - VV^T ||}{MAX}$$
            where MAX is $$MAX= \mathfrak{m} upper^2$$ and  $$ \mathfrak{m} = max(max(UU^T),max(VV^T))$$
        AG: float
            agreement between input clusterings from same formula as G but where MAX is 
            $$MAX=||UU^T||+||VV^T||-2 ||UU^T||||VV^T|| / (\mathfrak{m} upper^2)$$.
           
            
    Examples
    --------
    >>> U = Matrix([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1],[1,0],[1,0],[1,0]])
    >>> V = Matrix([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[0,0,1]])
    >>> agreement_from_co_memberships(U,V)
    """
    
    CU = U.dot(U.T).astype('float')
    CV = V.dot(V.T).astype('float')
    print CU
#     OU = U.T.dot(U).astype('float')
#     OV = V.T.dot(V).astype('float')
#     OUV =  N = U.T.dot(V).astype('float')
#     printlatex(U.astype('int'))
#     printlatex(V.astype('int'))
#     printlatex(CU.astype('int'))
#     printlatex(CV.astype('int'))
#     printlatex(((CU-CV)).astype('int'))
    sf = lambda A: np.multiply(A, A).sum() #squared frob norm
#     asum = lambda A: np.sum(abs(A)) #squared frob norm

    norm = sf
    
    """    
        #rank is the number of clusters?
    #     print LA.matrix_rank(U),  LA.matrix_rank(CU)
    #     print LA.matrix_rank(V), LA.matrix_rank(CV)
    #     print LA.matrix_rank(CU-CV), LA.matrix_rank(U.T.dot(V))
       
    #     print "||UU'||=", norm(CU), "=||U'U||=", norm(OU ), 
    #     print " and ||VV'||=", norm(CV), "=||V'V||=", norm(OV)
    #     print "=||UU' -VV'||=", norm(CU-CV),"=||UU'||+||VV'||-2||U'V|| =", norm(CU)+norm(CV) - 2* norm(N), 
    #     print "   where ||U'V||=", norm(OUV)
    #     print " || UU' -VV' || =<", norm(CU-CV), " < ||UU'||+||VV'|| =", norm(CU)+norm(CV) , "<||U||||U'||+||V||||V'||=", norm(U)**2+norm(V)**2 
    #     print "<||UU'||+||VV'|| = ", norm(CU)+norm(CV) ," ?< ||U||||V|| = ", norm(U) * norm(V) 
    
   
    D = CU
    print ">>>>> ",LA.norm(D, np.inf),LA.norm(D, - np.inf),LA.norm(D, 1),LA.norm(D, -1),LA.norm(D, 2),LA.norm(D, -2) 
    D = CV
    print ">>>>> ",LA.norm(D, np.inf),LA.norm(D, - np.inf),LA.norm(D, 1),LA.norm(D, -1),LA.norm(D, 2),LA.norm(D, -2) 
    D = CU+CV
    print ">>>>> ",LA.norm(D, np.inf),LA.norm(D, - np.inf),LA.norm(D, 1),LA.norm(D, -1),LA.norm(D, 2),LA.norm(D, -2) 
    D = CU-CV
    print ">>>>> ",LA.norm(D, np.inf),LA.norm(D, - np.inf),LA.norm(D, 1),LA.norm(D, -1),LA.norm(D, 2),LA.norm(D, -2) 
    """
#     upper = OUV.sum()**2 
    if not same_nodes:
        np.fill_diagonal(CU,0)
        np.fill_diagonal(CV,0)
#         upper = upper -  OUV.sum()
    print U
    print U.T

    print CU
    pairCount =  CU.shape[0]**2 - (0 if same_nodes else CU.shape[0])

    sUV = norm( CU-CV  )
    sU = norm(CU)#- (0 if same_nodes else CU.shape[0])
    sV = norm(CV)#- (0 if same_nodes else CU.shape[0])
    '''
    print ('M ' if same_nodes else ' ')+ '#pairs=', pairCount, '> Max= ', upper
    print "|| UU' -VV' || = ",sUV, " < ||UU'||+||VV'|| =", sU+sV , math.sqrt(sU*sV) #";; ||U'V||=" ,norm(OUV), 
    print " ;; ",norm(CU+CV),';;',  norm(np.maximum(CU,CV)) ,';;',  norm(np.multiply(CU,CV)) ,';;',  norm(np.dot(CU,CV)) 
 
    print "                           sum(U'V)**2=\t", OUV.sum()**2 ,'\t', OUV.sum()**2-OUV.sum() 
    print "              1/2 (||U||**2+||V||**2) =\t",  .5 * ( norm(U)**2 + norm(V)**2) ,'\t', .5 * ( norm(U)**2 + norm(V)**2) - .5*(norm(U)+norm(V))
    print "                           ||U||||V|| =\t ", norm(U) * norm(V) ,'\t', norm(U) * norm(V)- .5*( norm(U) + norm(V))
    print 'max= ',m,  'paircount = ', pairCount            
    print sUV ," == ", sU + sV - 2*np.sum(np.multiply(CU, CV))
    '''
    m = max(np.max(CU),np.max(CV))

    GI =  sUV /(pairCount*m**2)
    AGI = sUV / ((sU+sV) - 2*(np.sum(CU)*np.sum(CV))/(pairCount))
#     print "-----------", sUV, sU, sV
#     print pairCount, ((sU+sV) - 2*(np.sum(CU)*np.sum(CV))/(pairCount))
    logging.debug('-->'+('same_nodes' if same_nodes else 'exact') +' formulation::\upper G={0}, AG={1}'.format(GI, AGI))    
    return  np.array([GI, AGI])

#  {0,1,2,3},{4,5,6}
# {0,1,2},{3,4,5,6}
# {0},{1,2,3,4,5,6}
def agreement_from_omega(U, V, same_nodes = False):
    CU = U.dot(U.T).astype('float')
    CV = V.dot(V.T).astype('float')
    if not same_nodes:
        np.fill_diagonal(CU,0)
        np.fill_diagonal(CV,0)
    pairCount =  CU.shape[0]**2 - (0 if same_nodes else CU.shape[0])

    pair_counts = np.zeros((np.max(CU)+1,np.max(CV)+1))
    print '---------------------',pair_counts.shape
    print CU
    print CV
    for i in range(0, int(np.max(CU)+1)):
        for j in range(0, int(np.max(CV)+1)): 
            tmp = (CU==i)*(CV==j)
            tmp = np.tril(tmp) if same_nodes else np.tril(tmp,-1)
            pair_counts[i,j] = np.sum(tmp) 
    
    print pair_counts
    
#     d = np.nonzero(CU-CV)
#     print '---------',d, pairCount    
    w=  np.array( CU==CV, dtype='float')
    w = np.sum(w)- w.trace()
    
    print w
    w = w*1.0/ pairCount
#     print '----->w:: ',w
#     print np.diag(CU), U.shape[0]
#     print np.diag(CV), U.shape[0]
    tU = np.histogram(CU, bins=np.max(CU)+1)[0]
    tV = np.histogram(CV, bins=np.max(CV)+1)[0]
    print '> ', tU
    print '> ', tV
    tU[0]= tU[0] - (0 if same_nodes else CU.shape[0])
    tV[0]= tV[0] - (0 if same_nodes else CU.shape[0])

    print '> ', tU
    print '> ', tV
    m =  min(tU.shape[0], tV.shape[0])
    ew = [tU[j]*tV[j] for j in range(0 ,m)] 
    ew = np.sum(ew)*1.0 / (pairCount**2)
#     print ew
    aw = (w-ew)/(1-ew)
#     print '>>>>>> aw:: ', aw
    return  np.array([w, aw])



# U = np.matrix([[1,0],[1,0],[1,0],[1,0],[0,1],[0,1],[0,1]])
# V = np.matrix([[1,0],[1,0],[1,0],[0,1],[0,1],[0,1],[0,1]])
# print agreement_from_omega(U, V,False)

def agreement_from_matrix_norms(U,V, norm = lambda A: math.sqrt(np.multiply(A, A).sum()) ):
    return 1-distance_from_matrix_norms(U, V, norm)
  
def distance_from_matrix_norms(U,V, norm = lambda A: math.sqrt(np.multiply(A, A).sum()) ):
    '''
    http://en.wikipedia.org/wiki/Matrix_norm
    http://ac.els-cdn.com/0024379584900788/1-s2.0-0024379584900788-main.pdf?_tid=003c824c-e818-11e3-aff0-00000aab0f27&acdnat=1401467733_788f98b9c9636cc1be862364b08b3aff
    http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.norm.html#numpy.linalg.norm
    '''
#     L21_norm  = lambda A: np.sqrt(np.multiply(A, A).sum(axis=1)).sum()  
#     frobenius_norm = lambda A: math.sqrt(np.multiply(A, A).sum()) 
#     log_norm = lambda A: np.multiply(A, np.log(A.clip(1)).sum(axis=1)).sum()    
    print '--------------------', U

    CU = U.dot(U.T).astype('float')
    CV = V.dot(V.T).astype('float')
    

#     print  frobenius_norm(Diff), L21_norm(Diff), log_norm(Diff), LA.norm(Diff, 'fro')/(m*n)
#     print LA.norm(Diff, 'fro')

#     normalized_norm = lambda norm, M1, M2: norm(M1-M2)/(norm(M1)+norm(M2)) 
#     print 'nnnnnnn ',np.array([norm(CU-CV)/(norm(CU)+norm(CV)) ])
    print '--------------------', CU
    print '-----Norm---------------', norm(CU)

    return np.array([norm(CU-CV)/(norm(CU)+norm(CV)) ])
     
def distance_alternative_forms(U,V):
    return 1-agreement_alternative_forms(U, V)

def agreement_alternative_forms(U, V):
#     if log_mode:
#         print 'tr(U*U.T*V*V.T) = ', (U*U.T*V*V.T).trace(), ' =  eigens : ',(U*U.T*V*V.T).eigenvals()
#         print 'tr(U*U.T) = ', (U*U.T).trace(), ' =  eigens : ',(U*U.T).eigenvals()
#         print 'tr(U.T*U) = ', (U.T*U).trace(), ' =  eigens : ',(U.T*U).eigenvals()
#         print 'tr(V*V.T) = ', (V*V.T).trace(), ' =  eigens : ',(V*V.T).eigenvals()
    
    CU = U.dot(U.T).astype('float')
    CV = V.dot(V.T).astype('float')
    
    CU2 = CU.dot(CU)
    CV2 = CV.dot(CV)
    
#     CUP2 = np.multiply(CU,CU)
#     CVP2 = np.multiply(CV,CV)
#     
    Q = CU.dot(CV)
#     P = np.multiply(CU,CV) #
#     
    # Cauchy-Schwarz inequality ==> alt1<1
    alt1 =   float(Q.trace()) / math.sqrt( (CU2).trace() * (CV2).trace() )
#     alt3 =   float(Q.sum()) / math.sqrt( CU2.sum() * CV2.sum() )
    # can show that:: math.sqrt( sum(CU**2) ) <= sum(CU)  : also true for any other power, and trace if semitive definite
    # also that::  (su A)/n <= ((su A**m)/n )**1/m   : also true for trace, for sum should be symmetric
    # ===>  1 <(su A)/sqrt(su A**2) < sqrt(n) ==>  1 < (su A)**2/ su (A**2) < n 
    # su(A+B)**2 / su(A+B)  <= su(A**2)/su(A) + su(B**2)/su(B) : if symmetric also for trace 
    # Submultiplicativity ==> alt2<1
    
    alt2 =  float(Q.trace() / (CU.trace() * CV.trace()))
#     alt4 =  Q.sum() / (CU.sum() * CV.sum())
# 
#     alt5 =  P.sum() /math.sqrt(CUP2.sum()* CVP2.sum()) # 5==1
#     alt6 =  P.sum() /(np.multiply(U,U).sum()*np.multiply(V,V).sum()) # 6==2
#     
#     alt7 =  np.trace(CU.dot(np.log(CU.clip(1))) - CU.dot(np.log(CV.clip(1))) - CU + CV)
#     if log_mode: 
#     print alt1, alt2, alt3, alt4, alt5, alt6, alt7
#     print '---------- ',alt7
    return  np.array([alt1, alt2])#, alt3, alt4])#, alt7]) #, [alt5, alt6]])

def get_all_agreement_variations(U, V):
    dists=  get_all_distance_variations(U,V)
    return [dists[0], 1-dists[1]]

def get_all_distance_variations(U, V):
    return [\
#             'VI','AVIo',
#         "RIolog","ARIolog",
#         'RIo','ARIo',
#         "RI'o","ARI'o",
            'RI','ARI',
            "RI'","ARI'", 
            'nF',# 'I_log', 
            'sqtr','tr'#,'sqsu','su'#,'nv'
            ] ,\
        np.hstack((\
#             agreement_from_overlaps(U,V, f= lambda x: 0 if x==0 else (x * np.log(x))),\
#             distance_from_overlaps(U,V, f = lambda x: x*math.log(x) if x>0 else 0),\
#             #[0],\
#             distance_from_overlaps(U,V, f = lambda x: x*(x-1), exact=True),\
# #             agreement_from_overlaps(U,V),\
#             distance_from_overlaps(U,V, f = lambda x: x**2),\
            distance_from_co_memberships(U,V, False),\
            distance_from_co_memberships(U,V),\
            distance_from_matrix_norms(U,V) ,\
            distance_alternative_forms(U, V) 
                   )) 

def get_incident_matrix(n, edges, directed =False):
    N = np.zeros((n, len(edges)), dtype =np.float)
    for j,e in enumerate(edges):
        w = 1 if len(e)<3 else np.sqrt(e[2])
        N[e[0]][j] = w
        N[e[1]][j] = w * (-1 if directed else 1)
    return N

def get_incident_matrix_from_graph(G, directed =False):
    #N = gm.incidence_matrix(G)        
    edges = [(i,j, w['weight']) for (i,j, w) in G.edges(data=True)]
    return get_incident_matrix(G.number_of_nodes(), edges,directed)
    
def get_incident_matrix_from_adjecency(A, directed =False):
    edges =  np.nonzero(A if directed else np.tril(A))
#     print edges
    print  A[edges]
    weighted_edges = np.vstack((edges, A[edges]))
    return get_incident_matrix(A.shape[0], np.transpose(weighted_edges), directed)
#     N = np.zeros((A.shape[0], edges.shape[0]), dtype =np.float)
#     for j,e in enumerate(edges):
# #         print e[0] --> e[1], A[e[0],e[1]]
#         N[e[0]][j] = sqrt(A[e[0],e[1]])
#         N[e[1]][j] = sqrt(A[e[0],e[1]]) * (-1 if directed else 1)
#     return N

def distance_with_structure(G,U):
    N = get_incident_matrix_from_graph(G)
    return get_all_distance_variations(N, U)

def agreement_with_structure(G,U):
#     A = gm.adj_matrix(G)
    dists=  distance_with_structure(G,U)
    return [dists[0], 1-dists[1]]
#     N = get_incident_matrix_from_graph(G)
#     N = gm.incidence_matrix(G)
#     print '--A--\n', A
#     print '--N--\n', N
#     n = A.shape[0]
#     NNT = A 
#     Bi = A.sum(axis=0)
#     Bo = A.sum(axis=1).T
#     for i in range(0,n):
#         NNT[i,i] = 0.5*(Bi[0,i] +Bo[0,i])
#     print '==A+D===\n', NNT - N.dot(N.T)
#     print '==NNT===\n', N.dot(N.T)
#     return get_all_agreement_variations(N, U)
#    get_all_agreement_variations(A, U)
def agreement_based_on_structure(G, U, V, method = "trans"):
    dists=  distance_based_on_structure(G, U, V, method)
    return [dists[0], 1-dists[1]]

def distance_based_on_structure(G, U, V, method = "trans"):
#     A = gm.adj_matrix(G)    
    N = get_incident_matrix_from_graph(G)

#     N = gm.incidence_matrix(G)
    #D = np.diag(A.sum(axis=0)+A.sum(axis=1))    
    #NNT = A + D
    
    # avg dis
    if method=="trans":
#         print N.T.dot(U), N.T.dot(V)
#         return get_all_agreement_variations(N.T.dot(U), N.T.dot(V))
        return get_all_distance_variations(N.T.dot(U), N.T.dot(V))
    else:
        names, d = get_all_distance_variations(U,V)
        d1 = get_all_distance_variations(U,N)[1]
        d2 = get_all_distance_variations(V,N)[1]
        if method == "sum2":
            alpha =0.5
            return [d1[0], alpha*d+ (1-alpha)*np.absolute(d1-d2)/np.max([d1,d2],axis=0)]
        elif method == "sum":
            alpha =0.5
    #         return [d1[0], (1- 0.5*np.add(d1[1] ,d2[1]))]
            return [d1[0], alpha*d+ (1-alpha)*np.absolute(d1-d2) ]
        elif method == "prod":
            return [d1[0], 1./3*(np.sqrt(d*d1)+np.sqrt(d*d2)+np.sqrt(d1*d2))]
         
 #        'UAV: UUTAATVVT' :  float((U*U.T * A*A.T *V*V.T).trace()/ math.sqrt((U*U.T * A * A.T * U*U.T).trace() * (V*V.T*A*A.T *V*V.T).trace())),
 #        'UNV: UUTNNTVVT' :  float((U*U.T * N*N.T *V*V.T).trace()/ math.sqrt((U*U.T * N *N.T * U*U.T).trace() * (V*V.T*N*N.T *V*V.T).trace()))

  #       }
    return float( (U.dot(U.T).dot(N.dot(N.T)).dot(V.dot(V.T))).trace()  /\
                  math.sqrt(   (U.dot(U.T).dot(N.dot(N.T)).dot(U.dot(U.T))).trace() *\
                            (V.dot(V.T).dot(N.dot(N.T)).dot(V.dot(V.T))).trace()  )\
                 )


