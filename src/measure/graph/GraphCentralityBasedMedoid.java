package measure.graph;

import java.util.Set;
import measure.base.Centroid;
import measure.graph.centrality.GraphVertexScorer;
import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.Graph;

/*
 * The base centrality class that should be extended by a specific centrality metric, which scores nodes based on their centrality within their network
 * It also caches the computed scores locally, each metric should have a getVertexScore where it computes and returns the score of that particular node or initializes all the score depending on which is more efficient 
TODO: add other centrality measures from Networks and see if we need to do normalization or not
 */

public class GraphCentralityBasedMedoid <V,E>  implements Centroid<V> {

	//The Graph
	protected Graph<V, E> graph;
	//Weight transformer
	protected Transformer<E, ? extends Number> weights;
	
	GraphVertexScorer<V, E> vertexScorer;
	//VertexScorerFactory<V, E, S> scorerFactory;
	
	public GraphCentralityBasedMedoid(GraphVertexScorer<V, E> vertexScorer, Graph<V, E> graph) {
		this(vertexScorer, graph, null);
	}

	public GraphCentralityBasedMedoid(GraphVertexScorer<V, E> vertexScorer,Graph<V, E> graph, Transformer<E, ? extends Number> weightTransformer) {
		this.graph = graph;
		this.weights = weightTransformer;
		this.vertexScorer = vertexScorer;
	}
	
	public V findCentroid(Set<V> nodes) {
		return findCentroid(FilterUtils.createInducedSubgraph(nodes, graph));
	}
	
	/*
	 * Returns the node with the highest centrality score
	 */
	protected V findCentroid(Graph<V, E> subGraph){
		vertexScorer.setGraph(subGraph,weights);
		
		V centroid = null;
		Double centroidScore = null; 
		
		for(V v : subGraph.getVertices())
			if(centroidScore == null || vertexScorer.getVertexScore(v) > centroidScore){
				centroid = v;
				centroidScore = vertexScorer.getVertexScore(v);
			}
		return centroid;
	}
	
	
//	protected void computeScores(Graph<V, E> subGraph){
//	vertexScores = new HashMap<V, Double>();
//	for (V v : subGraph.getVertices()){
//		vertexScores.put(v, getVertexScore(v));
//	}
//}

	
//	public V findCentroid(Graph<V,E> graph){
//		if(this.subGraph == null || graph != this.subGraph){
//			this.subGraph = graph;
//			initializeAllScores();
//		}
//		return findCentroid();
//	}
	
	public String toString(){
		return vertexScorer.toString();
	}
	
}
