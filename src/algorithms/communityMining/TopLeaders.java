package algorithms.communityMining;

import java.util.ArrayList;

import org.apache.commons.collections15.Transformer;

import static io.Logger.*;
import algorithms.topleaders.CommunityMiner;
import algorithms.topleaders.Partitioning;
import algorithms.topleaders.global.GlobalTopLeaders;
import algorithms.topleaders.local.LocalTopLeader;
import edu.uci.ics.jung.graph.Graph;


/**
 * @author Reihaneh
 *
 * @param <V>
 * @param <E>
 */
public class TopLeaders<V,E> extends CommunityMiner<V, E> {
 
	double o = -1, cc = -1, mc = -1, h = -1;
	
	Transformer<Graph<V, E>, Partitioning<V>>  topLeaders;
	
	public TopLeaders() {
		super();
		o = -1;
		cc = -1;
		mc = -1;
		h = -1;
	}

	/**
	 * Finds k communities in the given graph using global topleaders approach 
	 * @param graph
	 * @param k: number of communities
	 * @return communities
	 */
	public Partitioning<V> findCommunities(Graph<V,E> graph, int k){
	//	System.err.print(" find communities ");
		
		long l = System.currentTimeMillis();
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		double p[] = {0,0,5,30};//guessLocalParamsBasedOnGraph(graph);
		
		logln(" params estimated in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();

		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		topLeaders = new GlobalTopLeaders<V, E>(k, p[0],p[1], p[3]);
		
		logln(" topleader initialized in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();

		Partitioning<V> res = topLeaders.transform(graph);
		
		logln(" topleader partitioned input in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();

		return res;
	}

	
	public Partitioning<V> findCommunities(Graph<V,E> graph, int k, ArrayList<V> initialCenters, int distanceMeasureType){
			double p[] = guessLocalParamsBasedOnGraph(graph);
			//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
			topLeaders = new GlobalTopLeaders<V, E>(k, p[0],p[1], p[3]);
			((GlobalTopLeaders)topLeaders).setHardCodedInitialCenters(initialCenters);
			((GlobalTopLeaders)topLeaders).setDistanceMeasureType(distanceMeasureType);
			return topLeaders.transform(graph);
		}

	
	/**
	 * Finds k communities in the given graph using global topleaders approach with this specific outlier threshold (in icloseness scale).<br/>
	 * This is not a good method to be called if you don't know the exact value of the outlier threshold, 
	 * instead of this method call <code> setO(double o)</code>, to set the percentage of the nodes that should be removed as outliers and then call the <code> findCommunities(Graph<V,E> graph, int k) </code>
	 * @param graph
	 * @param k: number of communities
	 * @param outlierThreshold: the specific outlier threshold
	 * @return communities
	 * @see algorithms.communityMining.TopLeaders#setO(double)
	 */
	private Partitioning<V> findCommunities(Graph<V,E> graph, int k, double outlierThreshold){
//		int centersClosenessThreshold = guessCenterClosenessThreshold(graph);;
		setO(outlierThreshold);
		return findCommunities(graph,k);
	}
	
	/**
	 * Finds k communities in the given graph using global topleaders approach with this specific outlier and centersCloseness threshold (in icloseness scale).<br/>
	 * This is not a good method to be called if you don't know the exact value of the outlier and/or centersCloseness threshold, 
	 * instead of this method call <code> setO(double o)</code>, to set the percentage of the nodes that should be removed as outliers and also call <code> setC(double c)</code>, to set the percentage to what the centers could be closed together. Then call the <code> findCommunities(Graph<V,E> graph, int k) </code>
	 * @param graph
	 * @param k : number of communities
	 * @param outlierThreshold : the specific outlier threshold ([0..10]: 0 is no outliers should be detected)
	 * @param centersClosenessThreshold : the specific centerCloseness threshold ([0..10]: 5)
	 * @return communities
	 * @see algorithms.communityMining.TopLeaders#setO(double)
	 * @see algorithms.communityMining.TopLeaders#setC(double)
	 */
	private Partitioning<V> findCommunities(Graph<V,E> graph, int k, double outlierThreshold, int centersClosenessThreshold){
		setCC(centersClosenessThreshold);
		return findCommunities(graph,k,outlierThreshold);
	}
	
	
	/**
	 * Finds communities in the given graph in two passes using local & global topleaders  
	 * @param graph: undirected, unweighted graph
	 * @return communities
	 */
	public Partitioning<V> findCommunities(Graph<V,E> graph){
		//centersClosenessThreshold = 600;
		double p[] = guessLocalParamsBasedOnGraph(graph);
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		topLeaders = new LocalTopLeader<V, E>(p[0],p[1],p[2],p[3]);
		//return topLeaders.transform(graph);
	
		int k = topLeaders.transform(graph).getNumberOfClusters();
		System.err.println("Rough clustering found " + k +  " clusters");
		return findCommunities(graph,k);
	}
	
	
	private double[] guessLocalParamsBasedOnGraph(Graph<V,E> graph){
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		double[] res = {0,0,0,0};
		//Effective Parameters: N, M, CC, Diameter, Average degree
		int n = graph.getVertexCount();
		//System.err.println("_____ "+ cc);
		int m = graph.getEdgeCount();
		
		if(mc<0)
			res[2] = Math.round(Math.log(n))+1 ;//Clustering coefficient
		else res[2] = mc;
		
		double[] cluC = avgCluteringCoefficient(graph), deg = getAverageDegree(graph);
	//	System.err.println("density: " + 2.*m/(n*(n-1)) +", cc: "+cluC[0]  +" +- "+ cluC[1] +", d: "+ deg[0] +" +- "+ deg[1]);
		
		if(o<0){
			res[0] = (Math.atan( (.6 -  (cluC[0]+cluC[1])) *10 )/Math.PI  + .5) * 18 ;
		}
		else if ( o >0 && o<=1)
			res[0] = o  * (deg[0] + deg[1] )* (deg[0] + deg[1] ); 
		else res[0] = o ;
		
		if(cc==-1){
			res[3] = 14 + (deg[0] + deg[1]) * (cluC[0]);//16; 
		}else if ( cc >0 && cc<=1)
			res[3] = cc * (deg[0] + deg[1])  * (deg[0] + deg[1]) * (cluC[0]) + 14  ; // 14-19 football, karate: 16-60, polboks:5-50 , stike: anything! 
		else
			res[3] = cc;
			
		//res[3]=16;
		
		if(h<0)
			res[1] = 1;
		else if(h>0 && h<=1)  
			res[1] = h*4 ;
		else res[1] = h;
		
		//res[1]= 0;
		
		return res;
	}
	
//	private int guessCenterClosenessThreshold(Graph<V,E> graph){
//		if(cc!=-1) return (int ) (cc * 10) ;
//		int centerClosenessThreshold = 5;
////		int n = graph.getVertexCount();
////		int m = graph.getEdgeCount();
////		centerClosenessThreshold = (int)  (Math.atan((5 - (n*10/m))*10 )/Math.PI + .5)*400+300; //Math.max(1800 * cc * cc  ,300);//(.44) 300 <-> (.56) 700 potentialLeaderThreshold * cc * 250;
//		return centerClosenessThreshold;
//	}
//	
	private double[] avgCluteringCoefficient(Graph<V,E> graph){
		double cc = 0, varCC=0;
		for (V v:graph.getVertices()){
			double cv = 0;
			for(V n1: graph.getNeighbors(v)){
				for(V n2: graph.getNeighbors(v)){
					if(n1!=n2 && graph.findEdge(n1, n2)!=null)
						cv++;
				}
			}
			if(cv!=0)
				cv/= graph.degree(v) * (graph.degree(v)-1);
			cc += cv;
			varCC += cv*cv;
		}
		cc/=graph.getVertexCount();
		varCC =Math.sqrt( varCC/graph.getVertexCount() - cc*cc);
	//	System.err.println(cc);
		return new double[]{cc,varCC};
	}

	private double[] getAverageDegree(Graph<V,E> graph){
		double avgD = 0 , varD=0;
		for (V v:graph.getVertices()){
			avgD += graph.degree(v);
			varD += graph.degree(v) * graph.degree(v);
		}
		avgD /= graph.getVertexCount();
		varD = Math.sqrt(varD/graph.getVertexCount() - avgD*avgD);
		return new double[]{avgD,varD};
	}
	
//	private double guessOutlierThreshold(Graph<V,E> graph){
//		if(o!=-1) return o * 7;
//		double outlierThereshod;
//	//	double cc = avgCluteringCoefficient(graph);
//		double[] cluC = avgCluteringCoefficient(graph), deg = getAverageDegree(graph);
//		
//		outlierThereshod = (Math.atan( (.6 -  (cluC[0]+cluC[1])) *10 )/Math.PI  + .5) * 15;
//		
////		outlierThereshod = (Math.atan((.36 - cc)*20 )/Math.PI + .5)*600 ;
////		outlierThereshod *=  Math.atan(graph.getVertexCount()-100)/Math.PI + .5; //Make it lower for small graphs
////		outlierThereshod *= 4.5/550;
//		return outlierThereshod;
//	}
	
	
//	private double avgCluteringCoefficient(Graph<V,E> graph){
//		double cc = 0;
//		for (V v:graph.getVertices()){
//			double cv = 0;
//			for(V n1: graph.getNeighbors(v)){
//				for(V n2: graph.getNeighbors(v)){
//					if(n1!=n2 && graph.findEdge(n1, n2)!=null)
//						cv++;
//				}
//			}
//			if(cv!=0)
//				cv/= graph.degree(v) * (graph.degree(v)-1);
//			cc += cv;
//		}
//		cc/=graph.getVertexCount();
//	//	System.err.println(cc);
//		return cc;
//	}

	
/*	private double[] guessLocalParamsBasedOnGraph(Graph<V,E> graph){
		//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
		double[] res = {0,0,0,0};
		//Effective Parameters: N, M, CC, Diameter, Average degree
		int n = graph.getVertexCount();
		int m = graph.getEdgeCount();
//		double density =  m*1.0/(n*(n-1));
		
		if(mc==-1)
			res[2] = Math.round(Math.log(n))+1 ;//Clustering coefficient
		else res[2] = mc;
		
		double cluC = avgCluteringCoefficient(graph);
		//System.err.println(cluC);
		
	//	outlierThereshod = 1-cluC;
		
		//((arctan((4.3-6.8)*2)+pi/2 )/pi) * 600
		
//		outlierThereshod = Math.pow(outlierThereshod, 3);
//		outlierThereshod = Math.sin(Math.PI/2 * outlierThereshod);
//		outlierThereshod = Math.pow(outlierThereshod, 3)*2000;
		
		if(o==-1){
			res[0] = (Math.atan((.36 - cluC)*20 )/Math.PI + .5)*600 ;
			res[0] *=  Math.atan(n-100)/Math.PI + .5; //Make it lower for small graphs
		}
		else res[0] = o * 600;
		//cc = avgAgumentedCluteringCoefficient(graph);
		//System.err.println(">>>>>>>>>>>>>> "+ cc);
		
		if(cc==-1)
			res[3] =   (Math.atan((5 - (n*10/m))*10 )/Math.PI + .5)*400+300; //Math.max(1800 * cc * cc  ,300);//(.44) 300 <-> (.56) 700 potentialLeaderThreshold * cc * 250;
		else res[3] = cc * 400 + 300; 
			
		if(h==-1)
			res[1] = 10;
		else res[1] = h;
		
		return res;
	}
	*/
	
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
	
	public String toString(){
		return "TopLeaders";
	}

	public String getName() {
		return "iTopLeaders";
	}

	public String getShortName() {
		return "iTL";
	}
	
}
