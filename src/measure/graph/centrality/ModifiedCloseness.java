package measure.graph.centrality;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.collections15.Transformer;

import measure.graph.GraphCentralityBasedMedoid;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
import edu.uci.ics.jung.graph.Graph;

// Modified closeness centrality to handle unconnected graphs, based on page 185 Newman's Network Book 
public class ModifiedCloseness<V,E> extends GraphVertexScorer<V,E>{
/*
 * Computing betweenness centrality requires computing shortest path between all pair of nodes in the graph, therefore it is wiser to compute all the scores at one round
 */

	
    protected  void initializeAllScores() {
    	vertexScores = new HashMap<V, Double>();		
    	Distance<V> distance = weighted()?( new DijkstraDistance<V,E>(graph, weights)):( new UnweightedShortestPath<V,E>(graph));
		
		for (V v : graph.getVertices()) {
			Map<V, Number> v_distances = new HashMap<V, Number>(distance.getDistanceMap(v));
			v_distances.remove(v);
		
			Double score = 0.0;
			for (Number d : v_distances.values()){
				score += 1/d.doubleValue();
			}
			
			score*= 1/(graph.getVertexCount()-1);
			
			vertexScores.put(v,score);
		}
	}

	public Double getVertexScore(V v) {
		if(vertexScores == null) initializeAllScores();
		return vertexScores.get(v);
	}

	public String toString(){
		return "Closeness";
	}

}
