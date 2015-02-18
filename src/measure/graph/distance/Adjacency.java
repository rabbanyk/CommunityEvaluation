package measure.graph.distance;

import measure.graph.GraphProximity;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;

public class Adjacency <V, E> extends GraphProximity<V,E> {
	{
		type = Type.SIMILARITY;
	}
	
	public Adjacency(Graph<V, E> g) {
		this(g,null);
	}
	public Adjacency(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		this(g, weights,false);
	}
	public Adjacency(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		this(graph, weights, selfWeight,false);
	}
	public Adjacency(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	
	public Number computeMeasure(V source, V target) {
		return W(source, target);
	}
	
	public Number getDistance(V source, V target){
		return MaxWeight  - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return "Adjacency"+ type + ((selfWeight)?"Aug":"");
	}
	
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{A}$";
		else return "A";
	}

}
