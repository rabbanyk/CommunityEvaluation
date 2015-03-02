package algorithms.communityMining;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import util.Statistics;
import algorithms.dev_topleaders.CommunityMiner;
import algorithms.dev_topleaders.Partitioning;
import algorithms.dev_topleaders.local.LocalTopLeader;
import edu.uci.ics.jung.graph.Graph;

/**
 * @author Reihaneh
 *
 * @param <V>
 * @param <E>
 */
/**
 * @author Reihaneh
 *
 * @param <V>
 * @param <E>
 */
public class LocalTopLeaders<V,E> extends CommunityMiner<V, E> {
 
	double o = -1, cc = -1, mc = -1, h = -1;
	
	Transformer<Graph<V, E>, Partitioning<V>>  topLeaders;
	
	public LocalTopLeaders() {
		super();
		o = -1;
		cc = -1;
		mc = -1;
		h = -1;
	}

	
	/**
	 * Finds communities in the given graph using local topleaders approach 
	 * @param graph: undirected, unweighted graph
	 * @return communities
	 */
	public Partitioning<V> findCommunities(Graph<V,E> graph){
		//centersClosenessThreshold = 600;
		double p[] = guessLocalParamsBasedOnGraph(graph);
	//	System.err.println(p[0]+" "+p[1]+" "+p[2]+" "+p[3]);
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		topLeaders = new LocalTopLeader<V, E>(p[0],p[1],p[2],p[3]);//(0.78,.2,6,.8);//
		return topLeaders.transform(graph);
	
//		int k = topLeaders.transform(graph).getNumberOfClusters();
//		return findCommunities(graph,k);
	}

	
	
	
	private double[] guessLocalParamsBasedOnGraph(Graph<V,E> graph){
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		double[] res = {0,0,0,0};
		//Effective Parameters: N, M, CC, Diameter, Average degree
		int n = graph.getVertexCount();
		System.err.println("_____________"+n);
		int m = graph.getEdgeCount();
//		double density =  m*1.0/(n*(n-1));
		
		if(mc==-1)
			res[2] = Math.round(Math.log(n))+1 ;//Clustering coefficient
		else res[2] = mc;
		
		
		
		double[] cluC = Statistics.avgCluteringCoefficient(graph), deg = Statistics.getAverageDegree(graph);
		System.err.println("density: " + 2.*m/(n*(n-1)) +", cc: "+cluC[0]  +" +- "+ cluC[1] +", d: "+ deg[0] +" +- "+ deg[1]);
		
	//	outlierThereshod = 1-cluC;
		
		//((arctan((4.3-6.8)*2)+pi/2 )/pi) * 600
		
//		outlierThereshod = Math.pow(outlierThereshod, 3);
//		outlierThereshod = Math.sin(Math.PI/2 * outlierThereshod);
//		outlierThereshod = Math.pow(outlierThereshod, 3)*2000;
		
		if(o==-1){
			res[0] = (Math.atan( (.6 -  (cluC[0]+cluC[1])) *10 )/Math.PI  + .5) * 15 ;
		//	res[0] *=  Math.atan(n-100)/Math.PI + .5; //Make it lower for small graphs
		}
		else res[0] = o * 4; 
		
		//cc = avgAgumentedCluteringCoefficient(graph);
		//System.err.println(">>>>>>>>>>>>>> "+ cc);
		
		if(cc==-1){
			res[3] = (deg[0] / deg[0] ) * 8;//+ deg[1] ;//+ Math.atan(deg[1]-1)*2 - res[0]*2;//(Math.atan(avgD + varD -3 )/Math.PI + .5) * 5;
			//res[3] *= Math.atan()/Math.PI + .5;
			
			//Make it lower for small graphs
			//res[3] =  (Math.atan((5 - (n*10/m))*10 )/Math.PI + .5)*4+5; //(Math.atan((5 - (n*10/m))*10 )/Math.PI + .5)*4+5; //Math.max(1800 * cc * cc  ,300);//(.44) 300 <-> (.56) 700 potentialLeaderThreshold * cc * 250;
		}else res[3] = cc * 4 + 3; 
		
		res[3]=6;
		
		if(h==-1)
			res[1] = 1;
		else res[1] = h ;
		
		res[1]= 0;
		
		return res;
	}
//	public double avgDegree(Graph<V,E> graph){
//	
//		return res;
//	}
	
	
	/**
	 * Sets the percentage/scale of noise/outliers that should be removed. Set this before any community detection attempt.
	 * If you don't explicitly set this, the algorithm would itself find the appropriate value based on the characteristics of the given network.
	 * @param o : a value between 0 and 1 that determines the scale of the noise that should be removed. <br\> 
	 * 0 (no noise should be detected) to 1 (be very hard with the noise) 
	 */
	public void setO(double o) {
		this.o = o;
	}

	/**
	 * Sets the percentage that centers of clusters could be placed closed together. Set this before any community detection attempt.
	 * If you don't explicitly set this, the algorithm would itself find the appropriate value based on the characteristics of the given network.
	 * @param c : a value between 0 and 1 that determines that scale. <br\> 
	 * 0 (shouldn't be close at all) to 1 (could be very close) 
	 */
	public void setCC(double c) {
		this.cc = c;
	}

	/**
	 * Sets a minimum for the size of communities that should be detected. Set this before any community detection attempt.
	 * If you don't explicitly set this, the algorithm would itself find the appropriate value based on the characteristics of the given network.
	 * @param mc : a value between 2 and n that limits the size of the smallest cluster, where n is the total number of nodes. <br\> 
	 */
	public void setMC(double mc) {
		this.mc = mc;
	}

	/**
	 * Sets the percentage of hubs/overlaps that should be considered. Set this before any community detection attempt.
	 * If you don't explicitly set this, the algorithm would itself find the appropriate value based on the characteristics of the given network.
	 * @param h : a value between 0 and 1 that determines that scale. <br\> 
	 * 0 (shouldn't be looking for overlap at all) to 1 (consider most of possible overlaps) 
	 */
	public void setH(double h) {
		this.h = h;
	}
	
}
