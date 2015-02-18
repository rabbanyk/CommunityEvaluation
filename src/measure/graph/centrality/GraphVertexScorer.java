package measure.graph.centrality;

import java.util.HashMap;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.algorithms.scoring.VertexScorer;
import edu.uci.ics.jung.graph.Graph;

public abstract class GraphVertexScorer<V,E> implements VertexScorer<V, Double>{

	//The Graph in which scores must be calculated, could be a subgraph of the bigger graph
	protected Graph<V, E> graph;
	//Weight transformer
	protected Transformer<E, ? extends Number> weights;
	//Cache for computed scores
	protected HashMap<V, Double> vertexScores;
	

	public void setGraph(Graph<V, E> graph) {
		setGraph(graph,null);
	}

	public void setGraph(Graph<V, E> graph, Transformer<E, ? extends Number> weights) {
		this.graph = graph;
		this.weights = weights;
	}

	public abstract Double getVertexScore(V v);
	
	protected boolean weighted(){
		return weights!=null;
	}

	protected void clear(){
		vertexScores = null;
	}
	
}
