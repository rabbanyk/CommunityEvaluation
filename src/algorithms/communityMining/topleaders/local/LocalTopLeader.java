package algorithms.communityMining.topleaders.local;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.graph.centrality.ClusteringDegree;
import measure.graph.centrality.GraphVertexScorer;
import measure.graph.distance.ICloseness;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.topleaders.dev_.Partitioning;
import edu.uci.ics.jung.graph.Graph;

public class LocalTopLeader<V,E> implements Transformer<Graph<V, E>, Partitioning<V>>  {
	Graph<V, E> graph;
	ArrayList<V> sortedVertices;
	Set<V> seeds;
	//HashMap<V, Vector<Set<V>>> neighborhoods;
	private double outlierThereshod = 0,  hubThreshold = 0, potentialLeaderThreshold = 3, minCommunitySizeThreshold = 3, MAX_CLOSENESS = 10000; 
	
	double centersClosenessThreshold = 300;
	int nodeNeighborhoodThreshold = 3;//, sourceNeighborhoodThreshold = 2;
	Vector<V> powerfullButFree ;
	HashMap<V, HashMap<Community,Double>> memberships;
	
	ICloseness< V, E> iCloseness ;
	

	public LocalTopLeader() {
		super();
	}

	
	public LocalTopLeader(double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold) {
		super();
		this.outlierThereshod = outlierThereshod;
		this.hubThreshold = hubThreshold;
//		this.potentialLeaderThreshold = potentialLeaderThreshold;
		this.minCommunitySizeThreshold = minCommunitySizeThreshold;
		this.centersClosenessThreshold = centersClosenessThreshold;
//		this.nodeNeighborhoodThreshold = nodeNeighborhoodThreshold;
//		this.sourceNeighborhoodThreshold = sourceNeighborhoodThreshold;
	}

	public double avgCluteringCoefficient(Graph<V,E> graph){
		double cc = 0;
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
		}
		cc/=graph.getVertexCount();
	//	System.err.println(cc);
		return cc;
	}
	public double avgAgumentedCluteringCoefficient(Graph<V,E> graph){
		double cc = 0;
		int outees = 0;
		
		for (V v:graph.getVertices()){
			double cv = 0;
	
			for(V n1: graph.getNeighbors(v)){
				for(V n2: graph.getNeighbors(v)){
					if(n1!=n2 && graph.findEdge(n1, n2)!=null)
						cv++;
				}
			}
			if(graph.getNeighborCount(v)*100<= outlierThereshod) outees++;
			else {
				if(cv!=0)
					cv/= graph.degree(v) * (graph.degree(v)-1);
				cc += cv;
			}
		}
		System.err.println(outees);
		cc/= (graph.getVertexCount()-outees);
		return cc;
	}
	

	public double avgDegree(Graph<V,E> graph){
		double res = 0;
		for (V v:graph.getVertices()){
			res += graph.degree(v); 
		}
		res /= graph.getVertexCount();
		return res;
	}
	
	public Partitioning<V> transform(Graph<V, E> graph) {
		this.graph = graph;
		
		double avgD = 0 , varD=0;
		for (V v:graph.getVertices()){
			avgD += graph.degree(v);
			varD += graph.degree(v) * graph.degree(v);
		}
		avgD /= graph.getVertexCount();
		varD = Math.sqrt(varD/graph.getVertexCount() - avgD*avgD);
		
		potentialLeaderThreshold = avgD;//+varD;//+1; //+1 could be change to var
		initialize();
	
		System.err.println("Running Local with ["+"mc: "+ minCommunitySizeThreshold+", pl: "+potentialLeaderThreshold+", ot: "+outlierThereshod+ ", ht: "+hubThreshold+", cc: "+centersClosenessThreshold + ", nn: " +nodeNeighborhoodThreshold +" ] ");
//		System.err.println(iCloseness.getDistanceMap((V)new Integer(9)));
//		System.err.println(iCloseness.getDistanceMap((V)new Integer(10)));

		Partitioning<V> resultPartitioning = new Partitioning<V>();
		
		int c = 0;
		V seed = getNextSeed();
		while(seed != null){
			System.out.println("------------------------  "+seed);
			
			Set<V> followers = findFollowersOf(seed);
//			if(c++<3)
//			System.out.println(followers);
//			
			if(followers.size() >= minCommunitySizeThreshold){
				resultPartitioning.addCluster(followers);
				
				//Free previous commitments
				for (V v : followers){
//					System.err.println("--------------  "+v);
					Set<Community> toBeFreed = new HashSet<Community>();
					HashMap<Community,Double> commitments = memberships.get(v);
					for(Community belongness: commitments.keySet()){
						if(! belongness.contains(v)) {
							toBeFreed.add(belongness);
					//		System.out.println("- " + belongness);
						}
					}
					for(Community belongness: toBeFreed){
						commitments.remove(belongness);
					}
				//	System.out.println("- " + toBeFreed);
					//memberships.put(v, 	commitments);
					if(commitments.size() > 1){
						resultPartitioning.addHub(v);
						//System.out.println(commitments);
					}else if (resultPartitioning.isHub(v)){
						resultPartitioning.getHubs().remove(v);
					}
				}
				
				
//				System.out.println(seed);
//				System.out.println(followers);
			}else {
//				
//				if(seed.toString().equals("76")) {
//					System.err.println(memberships);
//					System.err.println("++++++++++++++++++++++++++++++++++");
//				}
//				
				seeds.remove(seed);
				
				Community selfBelongness = null;
				for (Community belongness : memberships.get(seed).keySet()) {
					if ( memberships.get(seed).get(belongness) == MAX_CLOSENESS)
						selfBelongness = belongness;
				} 
			
//				if(seed.toString().equals("5")) {
//					System.err.println(memberships);
//					System.err.println(resultPartitioning);
//				}
					
				//Return to your previous commitment cheaters
				for (V v : followers) {
					memberships.get(v).remove(selfBelongness);
					for (Community belongness : memberships.get(v).keySet()) {
						if (!belongness.contains(v))
							belongness.add(v);
					}
				}
				
				memberships.get(seed).remove(selfBelongness);
				if(memberships.get(seed).isEmpty())
					powerfullButFree.add(seed);
//
//				if(seed.toString().equals("76")) {
//					System.err.println(memberships);
//					System.err.println("-------------------------------------");
//				}
				
			}
			
			//System.out.println(memberships);
			seed = getNextSeed();
		}
		
		for (V v : powerfullButFree) {
			//System.out.println(v+" "+memberships.get(v));
			if(memberships.get(v)==null || memberships.get(v).isEmpty())
				resultPartitioning.addHub(v);
		}
		for(V v : sortedVertices){
			if(memberships.get(v) == null)
				resultPartitioning.addOutlier(v);
		}
		

		//if(cur.toString().equals("142"))

		//Non-overlapping
		for (V v :resultPartitioning.getHubs()){
			for (Set<V> cluster : resultPartitioning.getCommunities()) {
				if(cluster.contains(v))
					cluster.remove(v);
			}		
		}
	//	System.err.println(memberships.get(76));

		return resultPartitioning;
	}
	
	public void initialize(){
		seeds = new HashSet<V>();		
		sortedVertices = new ArrayList<V>(graph.getVertices());
		powerfullButFree = new Vector<V>();
		
		
		final GraphVertexScorer<V, E> centrality = new ClusteringDegree<V, E>();
		centrality.setGraph(graph);
		Collections.sort(sortedVertices,new Comparator<V>(){
			public int compare(V o1, V o2) {
			//	System.err.println(o1 +"("+centrality.getVertexScore(o1)+")" + (((int)(centrality.getVertexScore(o2) - centrality.getVertexScore(o1)))>0?"<":">") + o2 +"("+centrality.getVertexScore(o2)+")" );
				return (int)(centrality.getVertexScore(o2) - centrality.getVertexScore(o1));
			}});
		
//		final VertexScorer<V, ? extends Number> centrality =   new DegreeScorer<V>(graph);
//		//TODO: Sort vertexes based on their centrality(nlogn) or just take k most centrals(kn)? the current one is better if k << logn ??? 
//		Collections.sort(sortedVertices,new Comparator<V>(){
//			public int compare(V o1, V o2) {
//				return centrality.getVertexScore(o2).intValue() - centrality.getVertexScore(o1).intValue();
//			}});
		
	//	neighborhoods = new HashMap<V, Vector<Set<V>>>();
		iCloseness = new ICloseness<V, E>(nodeNeighborhoodThreshold,graph,null);
		memberships = new HashMap<V, HashMap<Community,Double>>();
	}

	// Returning next powerful leader
	public V getNextSeed(){
		V nextSeed = null;
		
		while(sortedVertices.size() > 0 && nextSeed == null){
			nextSeed = sortedVertices.get(0);
			
			//System.err.println(graph.degree(nextSeed));
			if(graph.degree(nextSeed) < potentialLeaderThreshold) 
				return null;
			sortedVertices.remove(0);
			
			powerfullButFree.add(nextSeed);
			for (V v : seeds) {
				if (iCloseness.getProximity(nextSeed, v).doubleValue() > centersClosenessThreshold){
					nextSeed = null;
					break;
				}
			}
			
		}
		
		if(nextSeed !=  null){
			seeds.add(nextSeed);
			powerfullButFree.remove(nextSeed);
		}
		
			
		return nextSeed;
	}
	
	public Set<V> findFollowersOf(V seed){
		//Set<V> followers = new HashSet<V>(); 
		Community followers = new Community();
		
		Vector<V> candidFollowers = new Vector<V>();
		Set<V> inQueue = new HashSet<V>();
		
		followers.add(seed);
		inQueue.add(seed);
		candidFollowers.add(seed);
		
		//candidFollowers.addAll(graph.getNeighbors(seed));
		while(!candidFollowers.isEmpty()){
			V cur = candidFollowers.get(0);
			
			candidFollowers.remove(0);
			
			if(branch(cur, seed, followers))
				for(V v: graph.getNeighbors(cur)){
					if((!inQueue.contains(v)) && (!seeds.contains(v))){
						inQueue.add(v);
						candidFollowers.add(v);
					}
				}
			
		}
		return followers.getMembers();
	}
	
	private void selfMembership (V seed, Community followers){
		//Seed should not belonging to any community
		if(memberships.get(seed)!=null){
			for (Community belongness: memberships.get(seed).keySet()){
				belongness.remove(seed);
			}
		}else
			memberships.put(seed, new HashMap<Community,Double>());
		
		memberships.get(seed).put(followers, MAX_CLOSENESS);
	 }
	 
	private boolean branch(V cur, V seed, Community followers){
		if (cur == seed){ 
			selfMembership(seed, followers);
			return true;
		}
		double belong = iCloseness.getProximity(cur, seed).doubleValue();
	//	System.err.println(cur+" "+seed + " : " + belong);

		if(belong < outlierThereshod) return false;
		
//		if(cur.toString().equals("76"))
		
		boolean add = false;

		if(memberships.containsKey(cur)){
			HashMap<Community,Double> commitments = memberships.get(cur);
			
			for(Community belongness: commitments.keySet()){
				if(Math.abs(belong - commitments.get(belongness)) <=  hubThreshold ){ // they are tied
					add = true;
				}else if(belong > commitments.get(belongness)) {  
					add = true;
//					if(cur.toString().equals("76"))
//						System.err.println("- "+commitments.get(belongness));
					belongness.remove(cur); //the other should be removed
				}
			}
		}else{
			memberships.put(cur, new HashMap<Community,Double>());
			add = true;
		}
		if (add){
			memberships.get(cur).put(followers, belong);
//			if(cur.toString().equals("76")){
//			System.err.println("+ "+cur+" > "+seed );
//			System.err.println(memberships.get(142));
//			}
			followers.add(cur);
		}
		
		return add;
	}
	
//	private double density(Set<V> nodes){	
//		if(nodes.size() <=1)	return 0;
//		
//		double den = 0;
//
//		for (V v : nodes) {
//			for (E e : graph.getIncidentEdges(v)) {
//				if(nodes.contains(graph.getOpposite(v, e))){
//					den += 1;//e.getNormalizedWeight();
//				}	
//			}
//		}
//		
//		return den/(2*nodes.size()*(nodes.size()-1));
//	}
	
//	public double iCloseness(V v1, V v2){
//		double res = 0, coun = 0;
//		double pow = 100;
//		
//		for (int i1 = 1, i2 = 1; i2 <= sourceNeighborhoodThreshold && i1<=nodeNeighborhoodThreshold; pow/=10, coun++) {
//			res+= pow*iCloseness(v1,i1, v2,i2);
//			if(coun%2==0 && i1<nodeNeighborhoodThreshold) i1++; else i2++;
//		}
//		
//		//res += iCloseness(v1,3, v2,2);
////		System.err.println(v1+" "+v2 + " : " + res);
//		return res;
//	}
	
//	public double iCloseness(V v1, int d1, V v2, int d2){
//		double res = 0;
//		
//		Set<V> tmp = new HashSet<V>();
//		//Computing intersection
//		for (V v : getNeighbors(v2).get(d2)) {
//			if(getNeighbors(v1).get(d1).contains(v)){
//				tmp.add(v);
//			//	res += 1;//getBelongness(v,v1) * getBelongness(v,v2);
//			}
//		}
//		
//		if (tmp.size()==0) return 0;
//		
//		res += tmp.size();
//		
////		res = 5 * Math.log(res+1);
//
//		double dens = density(tmp);
//		if(dens<0 || dens >1) System.err.println("Error in calculating density: " + dens);
//		res += dens;
//		
//		return res;
//	}
//	
//	private Vector<Set<V>> getNeighbors(V n) {
//
//		if (neighborhoods.containsKey(n))
//			return neighborhoods.get(n);
//
//		Vector<Set<V>> differentLevelNeighbors = new Vector<Set<V>>();
//
//		Set<V> neighbors = new HashSet<V>();
//
//		// length zero neighborhood
//		neighbors.add(n);
//		differentLevelNeighbors.add(new HashSet<V>(neighbors));
//
//		// length one to threshold neighborhoods
//		for (int i = 1; i <= sourceNeighborhoodThreshold|| i <= nodeNeighborhoodThreshold; i++) {
//			for (V v : differentLevelNeighbors.lastElement()){
//				neighbors.addAll(graph.getNeighbors(v));
//			}
//			differentLevelNeighbors.add(new HashSet<V>(neighbors));
//		}
//		neighborhoods.put(n, differentLevelNeighbors);
//		return differentLevelNeighbors;
//	}
//	
	
	private class Community {
		Set<V> nodes;
		
		public Community() {
			this.nodes = new HashSet<V>();
		}
		
		public Community(Set<V> nodes) {
			this.nodes = nodes;
		}
		public boolean contains(V v){
			return nodes.contains(v);
		}
		public void add(V v){
			nodes.add(v);
		}
		public void remove(V v){
			nodes.remove(v);
		}
		@Override
		public String toString() {
			return nodes.toString()+"\n";
		}
		
		public Set<V> getMembers(){
			return nodes;
		}
	}

	
}
