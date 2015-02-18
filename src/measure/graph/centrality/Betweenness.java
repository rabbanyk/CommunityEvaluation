package measure.graph.centrality;

import java.util.HashMap;

import edu.uci.ics.jung.algorithms.scoring.BetweennessCentrality;


public class Betweenness<V,E> extends GraphVertexScorer<V,E>{

	protected void computeScores(){
		// Using Jung Betweenness Centrality, TODO: weighted version might have some problems as indicated in their documentations
		BetweennessCentrality<V,E> bc = weighted()?(new BetweennessCentrality<V, E>(graph, weights)):(new BetweennessCentrality<V, E>(graph));
		
		//to normalized the betweenness centrality 
		//double n = graph.getVertexCount();
		//TODO: do we really need this? double divisor = (n-1)*(n-2)/2; 
			
		vertexScores = new HashMap<V, Double>();
		for(V v : graph.getVertices()){
			double score = bc.getVertexScore(v);
			this.vertexScores.put(v, score );// / divisor);	
		}
				
	}
	
/*
 * Computing betweenness centrality requires computing shortest path between all pair of nodes in the graph, therefore it is wiser to compute all the scores at one round
 */
	public Double getVertexScore(V v) {
		if(vertexScores == null) computeScores();
		return vertexScores.get(v);
	}

	public String toString(){
		return "Betweenness";
	}


}
