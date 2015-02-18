package measure.graph.centrality;

import java.util.HashMap;
import java.util.Map;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;

public class DistanceVariance<V,E>  extends GraphVertexScorer<V,E>{

	
	public  void initializeAllScores() {
		//DijkstraDistance<V, E> dd = new DijkstraDistance<V, E>(graph);
		vertexScores = new HashMap<V, Double>();
		
		Distance<V> distance = weighted()?( new DijkstraDistance<V,E>(graph, weights)):( new UnweightedShortestPath<V,E>(graph));
		
		for (V v : graph.getVertices()) {
			Map<V, Number> v_distances = new HashMap<V, Number>(distance.getDistanceMap(v));
			v_distances.remove(v);
			
			Double score = 0.0, mean = 0.0;

			for (Number d : v_distances.values()){
				score += d.doubleValue();
			}
			
			mean = score / (double) v_distances.values().size();

			score = 0.0;
			
			for (Number d : v_distances.values()) { //Note the justin's version was d/mean!
				score += (d.doubleValue() - mean)* (d.doubleValue() - mean);
			}
			
			score = v_distances.size()/score; //reverse variance as large variance are not good
			
			this.vertexScores.put(v, score);	
		}
	}
	
	public Double getVertexScore(V v) {
		if(vertexScores == null) initializeAllScores();
		return vertexScores.get(v);
	}

	public String toString(){
		return "DistanceVariance";
	}


}
