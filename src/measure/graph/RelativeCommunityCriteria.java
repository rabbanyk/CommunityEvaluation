package measure.graph;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import measure.base.RelativeCriteria;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.Graph;

public abstract class RelativeCommunityCriteria<V,E> implements RelativeCriteria<V>{
	protected Graph<V, E> graph;
	protected Transformer<E, ? extends Number> weights;
	//protected Vector<Graph<V,E>> communities;
	protected double MaxWeight, SumWeight;
	
	protected boolean selfWeight = false; // (A(i,i) == MaxWeight )? 
	protected boolean normalizedWeight = false; // (A(i,i) == (A(i,i)/MaxWeight )?
	
	protected boolean weighted(){
		return weights!=null;
	}
	
	public RelativeCommunityCriteria (Graph<V, E> graph){
		this(graph,null);
	}
	
	public RelativeCommunityCriteria (Graph<V, E> graph,Transformer<E, ? extends Number> weights){
		this(graph,weights, false, false);
	}
	
	public RelativeCommunityCriteria(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super();
		this.graph = graph;
		this.weights = weights;
		this.selfWeight = selfWeight;
		this.normalizedWeight = normalized;
		
		if (weights==null){
			MaxWeight = 1;
			SumWeight = graph.getEdgeCount();
		}else{
			MaxWeight = 0;
			SumWeight = 0;
			for (E e : graph.getEdges()) {
				double w = weights.transform(e).doubleValue();
				if( w > MaxWeight) MaxWeight = w;
				SumWeight += w;
			}
		}
	}
	
	protected double W(E e){
		if(e==null) return 0; 
		return weighted()? weights.transform(e).doubleValue() : 1;
	}
	
	protected double W(V source, V target){
		E e = graph.findEdge(source, target);

		if(e==null) return (selfWeight && (target == source))? MaxWeight :  0; 
		return weighted()? weights.transform(e).doubleValue()/(normalizedWeight?MaxWeight:1) : 1;
	}
	
	
	protected  Set<V> getUnionNeighbours (V source, V target) {
		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
		sNeigh.addAll(tNeigh);
		sNeigh.add(source);
		sNeigh.add(target);
		return sNeigh;
	}
	
	protected abstract double evaluateCommunities(Vector<Graph<V,E>> communities);

	public double evaluate(Vector<Set<V>> clusters) {
		return evaluateCommunities(new Vector<Graph<V,E>>(FilterUtils.createAllInducedSubgraphs(clusters, graph)));
	}
	
}
