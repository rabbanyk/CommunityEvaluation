# CommunityEvaluation
A Java framework for evaluating community mining algorithms. 
Which includes:

* Quality measures for a community structure (e.g. Q-modularity, ) presented in these paper:
  * [R Rabbany et al., *Relative Validity Criteria for Community Mining Algorithms*; ASONAM 2012](http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6425753&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D6425753)
  * [R Rabbany et al., *Communities validity: methodical evaluation of community mining algorithms*; SNAM 2013](http://link.springer.com/article/10.1007%2Fs13278-013-0132-x)
  * see `measure\criteria`
* Agreement indices to compare two given community structures (e.g. ), presented in following paper:
  * [R Rabbany et al., *Generalization of Clustering Agreements and Distances for Overlapping Clusters and Network Communities*; CORR 2014](http://arxiv.org/abs/1412.2601)
  * see `measure\cluster`
  * see a summary and examples in this [ipython notebook](http://nbviewer.ipython.org/github/rabbanyk/CommunityEvaluation/blob/master/Clustering%20Agreement.ipynb#)
* Implementation of the TopLeaders algorithm presented in:
  * [R Rabbany et al., *Top leaders community detection approach in information networks*, SNA-KDD 2010](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.419.1134)
  * see `algorithms\topleaders`
* Different graph based distances and centrality measures
  * see `measure\graph`  
* Java wrappers for a collection of community mining algorithms
  * see `algorithms\communityMining` and `execs\`  
  * most of the executable files would need to be re-build on the specific machine
