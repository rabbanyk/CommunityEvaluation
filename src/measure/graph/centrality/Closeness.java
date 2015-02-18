package measure.graph.centrality;

import java.util.HashMap;
import edu.uci.ics.jung.algorithms.scoring.ClosenessCentrality;

// Not good for unconnected graphs, 
public class Closeness<V,E>  extends GraphVertexScorer<V,E>{

	public void initializeAllScores() {
		// Using Jung Closeness Centrality, 
		ClosenessCentrality<V,E> cc = weighted()?(new ClosenessCentrality<V, E>(this.graph, weights)):(new ClosenessCentrality<V, E>(graph));
		
		//to normalized the betweenness centrality 
			
		vertexScores = new HashMap<V, Double>();
		for(V v : this.graph.getVertices()){
			double score = cc.getVertexScore(v);
			this.vertexScores.put(v, score );// / divisor);	
		}
				
	}
	
/*
 * Computing betweenness centrality requires computing shortest path between all pair of nodes in the graph, therefore it is wiser to compute all the scores at one round
 */
	public Double getVertexScore(V v) {
		if(vertexScores == null) initializeAllScores();
		return vertexScores.get(v);
	}

	public String toString(){
		return "Closeness";
	}

}
