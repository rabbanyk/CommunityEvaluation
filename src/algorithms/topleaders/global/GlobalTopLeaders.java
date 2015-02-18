package algorithms.topleaders.global;

import static io.Logger.logln;
import io.Logger.DebugMode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
//import static io.Logger.*;


import measure.base.Proximity;
import measure.graph.GraphCentralityBasedMedoid;
import measure.graph.centrality.ClusteringDegree;
import measure.graph.centrality.Degree;
import measure.graph.centrality.GraphVertexScorer;
import measure.graph.distance.AdjacencyRelation;
import measure.graph.distance.ICloseness;
import measure.graph.distance.NeighborOverlapOld;
import measure.graph.distance.AdjacencyPearsonCorrelation;
import measure.graph.distance.ShortestPath;

import org.apache.commons.collections15.Transformer;

import algorithms.topleaders.Partitioning;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.Graph;

/**
 * An approach for community mining in social networks.
 * This Class detects given number (k) of communities in a given Network. 
 * The algorithm is based on k-means clustering method.
 * It consists of these two main steps
 * calculate centers of k clusters based on a centrality measure
 * determine clusters based on current k centers using iCloseness closeness measure
 * 
 */
public class GlobalTopLeaders<V , E> implements Transformer<Graph<V, E>, Partitioning<V>> {
    protected Random rand;
	 //TODO: there might not be k centers based on clustering method, caused null pointer, change every numofclusters to number of centers, or after init change it
    private Graph<V, E> graph;
   // private Initialization<V, E> initialization;
    Transformer<E, ? extends Number> weights;
    private HashSet<V> nodesToConsider;
    private HashSet<V> centerNodes;
    boolean checkSecond = true;
    private int iterated = 0;
    private int maxIteration = 10;
	private double outlierThereshod = 0, hubThreshold = 0, centersClossenessThreshold = 4;
//	private double MaxWeight = 0;
	
    private int numClusters;

    public int distanceMeasureType = 0;
    public final int ICLOSENESS=0,SP = 1, ADR=2, NOD=3, PCD=4;

	public void setDistanceMeasureType(int distanceMeasureType) {
		this.distanceMeasureType = distanceMeasureType;
	}

	private ArrayList<V> centers;
	int neighborhoodThreshold = 2;//3;
	HashMap<V, Set<V>> neighborhoods;
	Proximity< V> distanceMeasure ;
	private Partitioning<V> clustring;
	
	/**
	 * @param numClusters : the number of clusters to be detected in the graph
	 * @param initMethod : initialization method for choosing initial centers, "k_central_common_neighbor" is recommended from [random, k_centrals, ck_central_random, roulette_wheel, k_central_one_neighbor_farther, [k_central_common_neighbor]]
	 * @param neighborhoodThreshold : how far from each node we should search for the center (in terms of SP depth)
	 * @param sourceNeighborhoodThreshold : how far from a center we should consider the nodes 
	 * @param outlierThereshod : how close should be a node to a center to be considered close, otherwise it would be marked as an outlier 
	 * @param algorithm : use basic k-means or k-means refine, default is fine 
	 * @param centralityMethod : the centrality method to be used: recommended "degree" from [ distance, [degree], betweenness]
	 * @param maxIteration : maximum number of iteration for the algorithm
	 * @param centersClossenessThreshold : preferred farness of the initial centers, just used in "k_central_common_neighbor" initialization method
	 * @param initBuffer : the buffer size used in the initialization, it could be any number greater than 1
	 */
	public GlobalTopLeaders(int numClusters, double outlierThereshod, double hubThreshold, double centersClossenessThreshold) {
//		logFunction("kMeansCommunityDetector", DebugMode.brief);
//		logln("K: " + numClusters +", centralityMethodID: "+centralityMethod 
//				+",\n nodeNeighborhoodThresholdm: "+nodeNeighborhoodThresholdm	+", sourceNeighborhoodThreshold: ", DebugMode.brief);
		//this.maxIteration = maxIteration;
//		this.nodeNeighborhoodThreshold = nodeNeighborhoodThresholdm;
//		this.sourceNeighborhoodThreshold = sourceNeighborhoodThreshold;
	
//		System.err.println("Running Global with [ k:" +numClusters+", o: "+ outlierThereshod+", c: "+ centersClossenessThreshold+ ", h: " +hubThreshold+" ] ");
		this.hubThreshold = hubThreshold;
		this.outlierThereshod = outlierThereshod;
		this.numClusters = numClusters;
		this.centersClossenessThreshold = centersClossenessThreshold;
	}
	
	
	public Partitioning<V> transform(Graph<V, E> graph,Transformer<E, ? extends Number> weight){
		this.weights = weight;

//		//TODO: This is for normalizing the weight, input could be normalized for efficiency 
//		for (E e: graph.getEdges()){
//			if(weights.transform(e).doubleValue() > MaxWeight) MaxWeight = weights.transform(e).doubleValue();
//		}
//		
		return transform(graph);
	}
	
	public double getCloseness(V v1, V v2){ //TODO: check this i changed the distance
	//	System.err.println(v1 + " " + v2 + graph.containsVertex(v1) + " "+ graph.containsVertex(v2));
		double closeness =  distanceMeasure.getProximity(v1, v2).doubleValue();

		
		if(! ( distanceMeasure instanceof ICloseness<?, ?> )){//|| distanceMeasure instanceof OldICloseness<?, ?>  )){
			if (distanceMeasure instanceof NeighborOverlapOld<?, ?> || distanceMeasure instanceof AdjacencyPearsonCorrelation<?, ?>)
			closeness= 1-closeness;
			else closeness = 1/closeness;
		}
		//System.err.println(v1 + " " + v2 + " " + closeness);
		
		return closeness;
	}
	// K-means algorithm for finding clusters in in the given graph
	public Partitioning<V> transform(Graph<V, E> graph) {
//		logFunction("transform", DebugMode.normal);
		
		this.graph = graph;
		
//		logln("\nGraph:",DebugMode.detailed);
//		logln(graph,DebugMode.detailed);
//		
		//neighborhoods = new HashMap<V, Set<V>>();
		nodesToConsider = new HashSet<V>();
		
		
		if (numClusters < 0 || numClusters > graph.getVertexCount()) {
			System.err.println(numClusters);
			throw new IllegalArgumentException("Invalid number of clusters passed in.");
		}

		long l = System.currentTimeMillis();
		
		
//		System.err.println(" Distance Measure  ");
		
		switch(distanceMeasureType){
		case ICLOSENESS: distanceMeasure = new ICloseness<V, E>(graph,weights);
		break;
		case SP: distanceMeasure = new ShortestPath<V, E>(graph, weights);
		break;
		case ADR: distanceMeasure = new AdjacencyRelation<V, E>(graph,weights);
		break;
		case NOD: distanceMeasure = new NeighborOverlapOld<V, E>(graph,weights);
		break;
		case PCD: distanceMeasure = new AdjacencyPearsonCorrelation<V, E>(graph,weights);
		break;		
		}
		
	
		 logln(" distance measure initialized in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();
		
		//new DijkstraDistance<V, E>(graph);
		//new AdjacencyRelationDistance<V, E>(graph,weights);
					//new PearsonCorrelationDistance(graph,weights);
		//	new OldICloseness<V, E>(graph,weights);
		//new NeighborOverlapDistance<V, E>(graph,weights);
		
		
//		System.err.println(" Initialization ");
		
		centers = initializeCenters(numClusters, graph);//initialization.initializeCenters(numClusters, graph);
		centerNodes = new HashSet<V>();
		
		logln(" centers initialized in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();
		
//		logln("Initial Centers: " + centers, DebugMode.brief);
		
//		System.err.println(" Update NodeSet ");

		updateNodeSet();
		
		logln(" NodeSet initialized in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();
	
//		logln("neighbors figured out", DebugMode.brief);
//		logln("Number of nodes to consider: " + neighborhoods.keySet().size() +" out of "+graph.getVertexCount(), DebugMode.brief);
		
//		System.err.println("Number of nodes to consider: " + nodesToConsider.size()+" out of "+graph.getVertexCount());
		
		int itr = 0;
		boolean changed = true;
		iterated = 0;
		
//		System.err.println("Initial leaders: " + centers);
		
		while(itr++ < maxIteration && changed){
			changed = false;
		//	logFunction("itr: " + itr, DebugMode.brief);
			
			for (V c : centers) {
				centerNodes.add(c);
				//c.setCenter(true);
			}

	//		System.err.println(clustring);
			
			changed |= findClusters();
			
			logln("   findClusters in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
			 l = System.currentTimeMillis();
		
			iterated++;
			
			changed |= findCenters();
			
			logln("   findCenters in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
			 l = System.currentTimeMillis();
			
			updateNodeSet();
		
			logln("   updateNodeSet in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
			 l = System.currentTimeMillis();
			 
//			System.err.println("Iteration: " + itr +", leaders: "+centers);
			
//			System.err.println("Number of nodes to consider: " + nodesToConsider.size()+" out of "+graph.getVertexCount());
		}
		
		//System.out.println(clustring);
	//	clustring.setCenters(centers);
		return clustring;
	}
	//HashMap<Pair<V>, Double> belongness = new HashMap<Pair<V>, Double>(); 
		
//	private double getNormalizedWeight(E e){
//		if(weights != null)
//		return weights.transform(e).doubleValue()/MaxWeight;
//		else return 1;
//	}
//	
//	
	private Set<V> getNeighbors(V n){
		//iCloseness contains the neighbourhood
		if(distanceMeasure instanceof ICloseness<?, ?>){
			return ((ICloseness)(distanceMeasure)).getNeighbors(n).keySet();
		}else {
			if (neighborhoods == null){
				neighborhoods = new HashMap<V, Set<V>>();
			}
			if(neighborhoods.containsKey(n))
				return neighborhoods.get(n);
			
			HashSet<V> neighbors = new HashSet<V>();
			
			// length zero neighborhood
			neighbors.add(n);
			
			// length one to threshold neighborhoods 
			for (int i = 1;  i <= neighborhoodThreshold; i++) {
				Set<V> expanded_neighbours =  new HashSet<V>(neighbors);
				for (V v : neighbors){ 	
					expanded_neighbours.addAll(graph.getNeighbors(v));
				}
				neighbors.addAll(expanded_neighbours);
			}
			
			neighborhoods.put(n,neighbors);
			return neighbors;
		}
	}

	// returns the ID of the cluster that this node should be added to 
	private Vector<V> findCluster(V v){
		Vector<V> candidates = new Vector<V> ();

		//Find candidate centers
		for (V n1 : getNeighbors(v)){
				for (V n2 : getNeighbors(n1)) {
					if(!candidates.contains(n2) && centerNodes.contains(n2))
						candidates.add(n2);
				}
			}
		
		//Find the most probable leader
		Vector<V> ties = new Vector<V>();
		double closenest = -1, closeness;
			
		for (int i = 0; i < candidates.size(); i++) {
		
			closeness = getCloseness(v, candidates.get(i));
			
		//	System.err.print(candidateClusters.get(i) + " " + closeness);
			
			if(closeness >= outlierThereshod){
				if(Math.abs(closeness - closenest) <=  hubThreshold ){//tie
					ties.add(candidates.get(i));
				//	System.err.println(node+ "  :   " +ties);
				}else if(closeness>closenest) {
					ties.clear();
					ties.add(candidates.get(i));
					closenest = closeness;
				}
			}
		}

	//	System.err.println(" > " + ties);
		return ties;
	}
	
	private boolean findClusters(){
		boolean changed = true;
		
		long l;
		 l = System.currentTimeMillis();

		
		Partitioning<V> new_clustring = new Partitioning<V>();

		for (int i = 0; i < centers.size(); i++) {
			new_clustring.addCluster(new HashSet<V>());
			new_clustring.getCommunities().lastElement().add(centers.get(i));
		}
		
		int counter = 0;
		
		logln("      findClusters setup in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();
		
		for (V v : nodesToConsider)//graph.getVertices()) 
		if(! centerNodes.contains(v)){
			counter++;
			//Find the most probable cluster for node V 
			Vector<V> tmp = findCluster(v);
			
			if(tmp == null || tmp.size()==0){
				if(graph.getNeighborCount(v)>5){
					new_clustring.addHub(v);//Powerful but free nodes
				}
				else
					new_clustring.addOutlier(v);
			}
			else if(tmp.size() == 1)
				new_clustring.getCommunities().elementAt(centers.indexOf(tmp.firstElement())).add(v);
			else 
				new_clustring.addHub(v);
			
		}
		
		logln("      findClusters main in:  " + (System.currentTimeMillis() - l ), DebugMode.normal );
		 l = System.currentTimeMillis();
		
		if(clustring != null && new_clustring.equals(clustring))
			changed = false;
		
		clustring = new_clustring;
		
		return changed;
	}

	private V centerOfCluster(Set<V> clusterNodes){
	
		GraphCentralityBasedMedoid<V, E> centrality = new GraphCentralityBasedMedoid<V, E>(new Degree<V,E>() ,  graph);
		//centrality.setGraph();
		//FilterUtils.createInducedSubgraph(clusterNodes
		return centrality.findCentroid(clusterNodes);
		
//		V result = null;
//		double maxCentrality = 0, tmp;
//		for (V v2 : clusterNodes) {
//			tmp = 0;//centrality.getVertexScore(v2);
//			for (E e : graph.getIncidentEdges(v2)) {
//				V v1 = graph.getEndpoints(e).getFirst();
//				if(v1.equals(v2))v1 = graph.getEndpoints(e).getSecond();
//				
//				if(clusterNodes.contains(v1)) tmp += getNormalizedWeight(e); 
//			}
//			
//			if (result == null || tmp > maxCentrality){
//				maxCentrality = tmp;
//				result = v2;
//			}
//		}
//		
//		return result;
	}


	private boolean findCenters(){
		boolean changed = true;
		
		for (V c : centers) {
			centerNodes.remove(c);
		}
				
		ArrayList<V> new_centers = new ArrayList<V>(); 
		for (Set<V> cluster : clustring.getCommunities()) {
			new_centers.add(centerOfCluster(cluster));
		}
		
		if (new HashSet<V>(new_centers).equals(new HashSet<V>(centers)))
			changed = false;
		
		
		centers = new_centers;

		return changed;
	}
	
	
	private void updateNodeSet(){
		nodesToConsider.clear();
		for (V center : centers) {
			for (V v : getNeighbors(center)){//sourceNeighborhoodThreshold)) {
				nodesToConsider.addAll(getNeighbors(v));
			}
		}
		
	}

	public int getIterated() {
		return iterated;
	}
	
	
	ArrayList<V> hardCodedInitialCenters;
	
	public void setHardCodedInitialCenters( ArrayList<V> centers){
		
		hardCodedInitialCenters = centers;
	}
	public ArrayList<V> initializeCenters(int numberOfCenters, Graph<V, E> graph) {
		ArrayList<V> centers = new ArrayList<V>();
		
		if (hardCodedInitialCenters != null) {
			for (V v : hardCodedInitialCenters) {
				V mapped = null;
				for (V v2 : graph.getVertices()) {
					if(v2.equals(v))
						mapped = v2;
				}
				if(mapped!=null)
				centers.add(mapped);
				else System.err.println("Invalid inital center");
			}		
			return centers;
		}
		
		int maxCandidates = 1; 

		ArrayList<ArrayList<V> > candidateCenters = new ArrayList<ArrayList<V>>();
		ArrayList<ArrayList<V> > remainingOfCandidateCenters = new ArrayList<ArrayList<V>>();

		final GraphVertexScorer<V, E> centrality = new ClusteringDegree<V, E>();
		centrality.setGraph(graph);
		//TODO: Sort vertexes based on their centrality(nlogn) or just take k most centrals(kn)? the current one is better if k << logn ??? 
		ArrayList<V> ver = new ArrayList<V>(graph.getVertices());

		Collections.sort(ver,new Comparator<V>(){
			public int compare(V o1, V o2) {
			//	System.err.println(o1 +"("+centrality.getVertexScore(o1)+")" + (((int)(centrality.getVertexScore(o2) - centrality.getVertexScore(o1)))>0?"<":">") + o2 +"("+centrality.getVertexScore(o2)+")" );
				return (int)(centrality.getVertexScore(o2) - centrality.getVertexScore(o1));
			}});
		
		
		boolean found = false;
		while (!found) {

			if(ver.size()==0) break;
			
			V result = ver.get(0);
	
			ver.remove(result);

				if(candidateCenters.size() == 0) {
					candidateCenters.add(new ArrayList<V>());
					remainingOfCandidateCenters.add(new ArrayList<V>());
					candidateCenters.get(candidateCenters.size()-1).add(result);
					continue;
				}

				Vector<ArrayList<V>> probs = new Vector<ArrayList<V>>();

				for (ArrayList<V> cc : candidateCenters) {
					probs.add(new ArrayList<V>());
			//		boolean add = true;
					for (V v : cc) {
//						Set<V> tmp = new HashSet<V>(graph.getNeighbors(v));
//						tmp.retainAll(graph.getNeighbors(result));
//						if (tmp.size() >= centersClossenessThreshold)// TODO: Could be improved by having depth
//							probs.lastElement().add(v);
						
							double tmp = getCloseness(v, result);
							if (tmp >= centersClossenessThreshold)// TODO: Could be improved by having depth
								probs.lastElement().add(v);
					}
				}
			
				ArrayList<ArrayList<V>> newCandidateCenters = new ArrayList<ArrayList<V>>();
				ArrayList<ArrayList<V> > newRemainingOfCandidateCenters = new ArrayList<ArrayList<V>>();
				
				for (int i = 0; i < candidateCenters.size(); i++) {
					if(probs.get(i).size()==0){
						newCandidateCenters.add(new ArrayList<V>(candidateCenters.get(i)));
						newCandidateCenters.get(i).add(result);
						newRemainingOfCandidateCenters.add(new ArrayList<V>(remainingOfCandidateCenters.get(i)));
						if(newCandidateCenters.get(newCandidateCenters.size()-1).size() == numberOfCenters){
							found = true;
							centers = newCandidateCenters.get(newCandidateCenters.size()-1);
						}
					}else{
						newCandidateCenters.add(new ArrayList<V>(candidateCenters.get(i)));
						newRemainingOfCandidateCenters.add(new ArrayList<V>(remainingOfCandidateCenters.get(i)));
						newRemainingOfCandidateCenters.get(newRemainingOfCandidateCenters.size()-1).add(result);
						
						if(newCandidateCenters.get(newCandidateCenters.size()-1).size() == numberOfCenters){
							found = true;
							centers = newCandidateCenters.get(newCandidateCenters.size()-1);
						}
					}
				}
				
				for (int i = 0; i < candidateCenters.size(); i++) {
					if(probs.get(i).size()!=0){
						newCandidateCenters.add(new ArrayList<V>(candidateCenters.get(i)));
						newCandidateCenters.get(newCandidateCenters.size()-1).add(result);
						newRemainingOfCandidateCenters.add(new ArrayList<V>(remainingOfCandidateCenters.get(i)));
						for (V prob : probs.get(i)) {
							newCandidateCenters.get(newCandidateCenters.size()-1).remove(prob);
							newRemainingOfCandidateCenters.get(newRemainingOfCandidateCenters.size()-1).add(prob);
						}
						
						if(newCandidateCenters.get(newCandidateCenters.size()-1).size() == numberOfCenters){
							found = true;
							centers = newCandidateCenters.get(newCandidateCenters.size()-1);
						}
					}
				}
				
				if(newCandidateCenters.size()>maxCandidates){
					candidateCenters = new ArrayList<ArrayList<V>>(newCandidateCenters.subList(0, maxCandidates));
					remainingOfCandidateCenters = new ArrayList<ArrayList<V>>(newRemainingOfCandidateCenters.subList(0, maxCandidates));
				}
				else {
					candidateCenters = newCandidateCenters;
					remainingOfCandidateCenters = newRemainingOfCandidateCenters;
				}
		}
		
			while(!found){
				for (int i = 0; i < candidateCenters.size(); i++) {
					if(remainingOfCandidateCenters.get(i).size()>0){
						candidateCenters.get(i).add(remainingOfCandidateCenters.get(i).get(0));
					}else found = true;
					
					//candidateCenters.get(i).add(remainingOfCandidateCenters.get(i).get(0));
				//	remainingOfCandidateCenters.get(i).get(0);
					if(candidateCenters.get(i).size() == numberOfCenters){
						centers = candidateCenters.get(i);
						found = true;
						break;
					}
				}
			
		}

		return centers;
	}
	
}

