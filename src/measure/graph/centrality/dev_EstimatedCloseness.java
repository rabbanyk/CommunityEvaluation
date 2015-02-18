package measure.graph.centrality;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;


import measure.graph.GraphCentralityBasedMedoid;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
import edu.uci.ics.jung.graph.Graph;

/**
 * An estimated version for closeness centrality based on the following paper:
 * David Eppstein and Joseph Wang, Fast approximation of centrality, In
 * Proceedings of the twelfth annual ACM-SIAM symposium on Discrete algorithms
 * (SODA '01).
 * 
 * In the estimated version only k sample vertices is selected and
 * Single-source-shortest-path is calculated only using these selected vertices.
 * The k samples is selected uniformly at random and for k=Theta(logn/e^2) the
 * error is e*Delta where Delta is the diameter of the graph.
 **/
//TODO: Justin's implementation, haven't check it, modify it like the original closeness to handle unconnected graphs
public class dev_EstimatedCloseness<V, E>  extends GraphVertexScorer<V,E>{

	public Distance<V> distance;
	// the number of sample vertices
	public int sampleSize;
	public ArrayList<Map<V, Number>> v_distances;
	public ArrayList<V> sampleVertices;

	// 0<epsion<=1 is the error bound
	public double epsilon = 1;

	protected void clear(){
		super.clear();
		distance = null;
		// the number of sample vertices
		sampleSize = 0;
		v_distances= null;
		sampleVertices = null;

		// 0<epsion<=1 is the error bound
		epsilon = 1;
	}
	
	protected void initializeAllScores() {
		vertexScores = new HashMap<V, Double>();		
    	Distance<V> distance = weighted()?( new DijkstraDistance<V,E>(graph, weights)):( new UnweightedShortestPath<V,E>(graph));
		
		sampleSize = (int) Math
				.floor((Math.log10(graph.getVertexCount()) / Math.pow(epsilon,
						2)));
		v_distances = new ArrayList<Map<V, Number>>(sampleSize);
		sampleVertices = new ArrayList<V>(sampleSize);

		V[] vertices = (V[]) graph.getVertices().toArray();

		Random random = new Random();

		while (sampleVertices.size() < sampleSize) {
			V sample = vertices[random.nextInt(graph.getVertexCount())];
			if (sampleVertices.contains(sample))
				continue;
			else {
				sampleVertices.add(sample);
				v_distances.add(distance.getDistanceMap(sample));
			}
		}

		for (V v : this.graph.getVertices()) {
			double score = getVertexScore(v);
			this.vertexScores.put(v, score);
		}
	
	
	}

	/**
	 * Calculates the score for the specified vertex. Here I assume the graph is connected.
	 */
	public Double getVertexScore(V v) {
		if(vertexScores == null) initializeAllScores();
		
		Double sum = 0.0;
		for (int w=0; w<sampleVertices.size(); w++) {
			if (sampleVertices.get(w).equals(v))
				continue;
			
			Number w_distance = v_distances.get(w).get(v);
			// TODO what can I do with do disconnected node?
			if (w_distance == null)
				continue;

			else
				sum += w_distance.doubleValue();
		}
		// TODO what if it is not connected?
		double score = sum == 0 ? Double.POSITIVE_INFINITY : ((graph.getVertexCount()-1)*sampleSize) / (sum*graph.getVertexCount());
			
		return score;
	}
	
	public String toString() {
		return "EstimatedCloseness";
	}

}
