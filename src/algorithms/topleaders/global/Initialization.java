package algorithms.topleaders.global;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.Vector;

import measure.graph.GraphCentralityBasedMedoid;
import measure.graph.centrality.GraphVertexScorer;

import edu.uci.ics.jung.algorithms.scoring.BetweennessCentrality;
import edu.uci.ics.jung.algorithms.scoring.DegreeScorer;
import edu.uci.ics.jung.algorithms.scoring.DistanceCentralityScorer;
import edu.uci.ics.jung.algorithms.scoring.VertexScorer;
import edu.uci.ics.jung.graph.Graph;

public class Initialization<V, E> {
	Random rand;
	public enum InitializationMethod {random, k_centrals, ck_central_random, roulette_wheel, k_central_one_neighbor_farther, component_based, personalized_PageRank, k_central_common_neighbor }
	private InitializationMethod initMethod = InitializationMethod.k_central_common_neighbor;
	int maxCandidates = 1; 
	
	GraphVertexScorer<V, E> centrality;
	
	public Initialization(GraphVertexScorer<V, E> centralityMethod) {
		rand = new Random();		
		this.centrality = centralityMethod;
	}

//	public ArrayList<V> initializeCenters(int numberOfCenters, Graph<V, E> graph){
//		centrality.setGraph(graph);
//
//		if(initMethod == InitializationMethod.random) return naiveInit(numberOfCenters, graph);
//		else if(initMethod == InitializationMethod.ck_central_random) return notSoNaiveInit(numberOfCenters, graph);
//		else if(initMethod == InitializationMethod.roulette_wheel) return rouletteWheel(numberOfCenters, graph);
//		//else if(initMethod == InitializationMethod.component_based) return componentBased(numberOfCenters, graph);
//		else return centralityInit(numberOfCenters, graph);
//	}
	 
	/*
	 * Choose k initial centers randomly
	 */
	 public ArrayList<V> naiveInit(int numberOfCenters, Graph<V, E> graph) {
		 
		ArrayList<V> centers = new ArrayList<V>();
		Vector<V> ver = new Vector<V>(graph.getVertices());
		int tmp;
		
		while (centers.size() < numberOfCenters) {
			tmp = rand.nextInt(ver.size());
			if (!centers.contains(ver.elementAt(tmp)))
				centers.add(ver.elementAt(tmp));
		}
		
		return centers;
	}
	 
     /*
	 * Choose k initial centers randomly from ck most central ones
	 */
	 public ArrayList<V> notSoNaiveInit(int numberOfCenters, Graph<V, E> graph) {
		 
		ArrayList<V> centers = new ArrayList<V>();

		int c = 4; //TODO: set c as an input?
		
		//initMethod = InitializationMethod.k_central_one_neighbor_farther;
		ArrayList<V> ver = centralityInit(c * numberOfCenters, graph);
		
		int tmp;
		
		while (centers.size() < numberOfCenters) {
			tmp = rand.nextInt(ver.size());
			if (!centers.contains(ver.get(tmp)))
				centers.add(ver.get(tmp));
		}
		
		return centers;
	}

 /*
	 * Choose k initial centers randomly
	 */
	 public ArrayList<V> rouletteWheel(int numberOfCenters, Graph<V, E> graph) {
		 
		ArrayList<V> centers = new ArrayList<V>();

	//	VertexScorer<V, ? extends Number> centrality = null;
		
//		if(centralityMethod == Centrality.betweenness)
//			centrality = new BetweennessCentrality<V, E>(graph);
//		else if (centralityMethod == Centrality.distance)
//			centrality = new DistanceCentralityScorer<V, E>(graph,true,true,true);
//		else if (centralityMethod == Centrality.degree)
//			centrality = new DegreeScorer<V>(graph);

		Vector<V> ver = new Vector<V>(graph.getVertices());
		Vector<Double> probs = new Vector<Double>(); 
		double sum = 0;
		
		for (V v : ver) {
			probs.add(centrality.getVertexScore(v).doubleValue());
			sum += centrality.getVertexScore(v).doubleValue();
		}
		
//		for (int i = 0; i < probs.size(); i++) {
//			probs.set(i, probs.get(i)/sum);
//		}
//		
		
		while (centers.size() < numberOfCenters) {
			double random = rand.nextDouble()*sum;
			
			for (int i =0; i< ver.size(); i++) {
				random -= probs.get(i);
				if(random<=0) {
					centers.add(ver.get(i));
					sum -= probs.get(i);
					probs.remove(i);
					ver.remove(i);
					break;
				}
			}
		}
			
		
		return centers;
	}
		HashMap<V, Vector<Set<V>>> neighborhoods;
	
	private Vector<Set<V>> getNeighbors(V n, Graph<V, E> graph) {
		if (neighborhoods.containsKey(n))
			return neighborhoods.get(n);

		Vector<Set<V>> differentLevelNeighbors = new Vector<Set<V>>();

		HashSet<V> neighbors = new HashSet<V>();

		// length zero neighborhood
		neighbors.add(n);
		differentLevelNeighbors.add(new HashSet<V>(neighbors));

		// length one to threshold neighborhoods
		for (int i = 1; i <= initSourceNeighborhoodThreshold
				|| i <= initNodeNeighborhoodThreshold; i++) {
			for (V v : differentLevelNeighbors.lastElement())
				neighbors.addAll(graph.getNeighbors(v));

			differentLevelNeighbors.add(new HashSet<V>(neighbors));
		}

		neighborhoods.put(n, differentLevelNeighbors);
		return differentLevelNeighbors;
	}
	 
	 int initNodeNeighborhoodThreshold = 1, initSourceNeighborhoodThreshold = 1;
		double centersClosenessThreshold = 2;
	 
	 
	/*
	 * Choose k initial centers from the most central nodes
	 */
	public void setcentersClosenessThreshold(double threshold) {
		this.centersClosenessThreshold = threshold;
	}
	
	
	
	public void setMaxCandidates(int maxCandidates) {
		this.maxCandidates = maxCandidates;
	}

	public ArrayList<V> centralityInit(int numberOfCenters, Graph<V, E> graph) {

		ArrayList<V> centers = new ArrayList<V>();
		neighborhoods = new HashMap<V, Vector<Set<V>>>();

//		ArrayList<V> reservedCenters = new ArrayList<V>();
		
		ArrayList<ArrayList<V> > candidateCenters = new ArrayList<ArrayList<V>>();
		ArrayList<ArrayList<V> > remainingOfCandidateCenters = new ArrayList<ArrayList<V>>();

		
//		final VertexScorer<V, ? extends Number> centrality = 
//			(centralityMethod == Centrality.betweenness)? new BetweennessCentrality<V, E>(graph): 
//				(	 (centralityMethod == Centrality.distance)? new DistanceCentralityScorer<V, E>(graph,true,true,true): new DegreeScorer<V>(graph));
//		else if (centralityMethod == Centrality.degree)

		
		// Sort vertexes based on their centrality(nlogn) or just take k most centrals(kn)? the current one is better if k << logn ??? 
		ArrayList<V> ver = new ArrayList<V>(graph.getVertices());
		Collections.sort(ver,new Comparator<V>(){
			public int compare(V o1, V o2) {
				return centrality.getVertexScore(o2).intValue() - centrality.getVertexScore(o1).intValue();
			}});
		
//		Vector<V> ver = new Vector<V>(graph.getVertices());

		boolean found = false;
		while (!found) {

			if(ver.size()==0) break;
			
			V result = ver.get(0);

			//Found the most central node
//			double maxCentrality = 0, tmpCen = 0;
//			for (V v2 : ver) {
//				tmpCen = centrality.getVertexScore(v2).doubleValue();
//				if (result == null || tmpCen > maxCentrality){
//					maxCentrality = tmpCen;
//					result = v2;
//				}
//			}
//			if(result == null) break;
	
			ver.remove(result);

			//Add result to centers based on the given initialization method 
			if(initMethod == InitializationMethod.k_centrals){
				centers.add(result);
				if (centers.size() == numberOfCenters) found = true;
			}else {
				if(candidateCenters.size() == 0) {
					candidateCenters.add(new ArrayList<V>());
					remainingOfCandidateCenters.add(new ArrayList<V>());
					candidateCenters.get(candidateCenters.size()-1).add(result);
					continue;
				}

				Vector<ArrayList<V>> probs = new Vector<ArrayList<V>>();

				for (ArrayList<V> cc : candidateCenters) {
					probs.add(new ArrayList<V>());
					for (V v : cc) {
						if (initMethod == InitializationMethod.k_central_one_neighbor_farther) {
							if (graph.getNeighbors(v).contains(result))
								probs.lastElement().add(v);
						} else if (initMethod == InitializationMethod.k_central_common_neighbor) {
							Set<V> tmp = new HashSet<V>(getNeighbors(v, graph).get(1));
							tmp.retainAll(getNeighbors(result, graph).get(1));
							//double tmp = 
							if (tmp.size() >= centersClosenessThreshold)// TODO: Could be improved by having depth
								probs.lastElement().add(v);
						}
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

		}
		
//		if(initMethod == InitializationMethod.k_central_common_neighbor){
//		for (V v : reservedCenters) {
//			boolean toAdd = true;
//			for (V c : centers) {
//				if(graph.getNeighbors(c).contains(v))
//					toAdd = false;
//			} 
//			if(toAdd){
//				centers.add(v);
//		//		reservedCenters.remove(v);
//			}
//			
//		}
//		}
		
		
//		for (int i = 0; i < candidateCenters.size(); i++) {
//			System.err.print(candidateCenters.get(i));
//			System.err.println(" + "+remainingOfCandidateCenters.get(i));
//		}
//		
		if(initMethod == InitializationMethod.k_central_common_neighbor || initMethod == InitializationMethod.k_central_one_neighbor_farther){
			while(!found){
				for (int i = 0; i < candidateCenters.size(); i++) {
					candidateCenters.get(i).add(remainingOfCandidateCenters.get(i).get(0));
					remainingOfCandidateCenters.get(i).get(0);
					if(candidateCenters.get(i).size() == numberOfCenters){
						centers = candidateCenters.get(i);
						found = true;
						break;
					}
				}
			}
			
		}
//		
		
//		while (centers.size() < numberOfCenters) {
//			if(!centers.contains(reservedCenters.get(0)))
//			centers.add(reservedCenters.get(0));
//			reservedCenters.remove(0);
//		}
		
//		System.err.println("FOUND: " + found);
		return centers;
	}
	
	/*
	 * This function could be used to assign the desired random_seed
	 * */
	protected void setRandomSeed(int random_seed){
		rand = new Random(random_seed);
	}
	
}
