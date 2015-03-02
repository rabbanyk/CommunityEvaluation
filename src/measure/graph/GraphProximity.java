package measure.graph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import measure.base.Proximity;
import measure.graph.distance.AdjacencyPearsonCorrelation;

import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;

public abstract class GraphProximity <V, E> extends Proximity<V> {
	//TODO: Should we cache data here? 
	protected Graph<V, E> graph;
	protected Transformer<E, ? extends Number> weights;
	protected double MaxWeight, SumWeight;
	
	protected boolean selfWeight = false; // (A(i,i) == MaxWeight )? 
	protected boolean normalizedWeight = false; // (A(i,i) == (A(i,i)/MaxWeight )?
	

	public GraphProximity(Graph<V, E> graph) {
		this(graph, null);
	}
	public GraphProximity(Graph<V, E> graph,
			 boolean selfWeight) {
		this(graph, null, selfWeight, false);
	}
	public GraphProximity(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		this(graph, weights, selfWeight, false);
	}

	public GraphProximity(Graph<V, E> graph, Transformer<E, ? extends Number> weights) {
		this(graph, weights, false, false);
	}
	
	public GraphProximity(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super();
		this.graph = graph;
		this.weights = weights;
		this.selfWeight = selfWeight;
		this.normalizedWeight = normalized;
		cache = new HashMap<Pair<V>, Number>();
		
		if(graph!=null)
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
	protected boolean weighted(){
		return weights!=null;
	}
	
	protected  Set<V> getUnionNeighbours (V source, V target) {
		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
		sNeigh.addAll(tNeigh);
		sNeigh.add(source);
		sNeigh.add(target);
		return sNeigh;
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
	
//	protected double nW(E e){
//		if(e==null) return 0; 
//		return weighted()? weights.transform(e).doubleValue()/MaxWeight : 1;
//	}
	
	protected void checkInput(V v) {
		if (graph.containsVertex(v) == false) throw new IllegalArgumentException("Specified target vertex " + v + " is not part of graph " + graph);
	}
	
	HashMap<Pair<V>, Number> cache;
	
	public Number getProximity(V source, V target) {//TODO: only if Zero?
		Pair<V> chacheKey = new Pair<V>(source,target);
		if(cache.containsKey(chacheKey)) return cache.get(chacheKey); 
		
		checkInput(target);
		checkInput(source);
		
		Number proximity = computeMeasure(source, target); 
		if (proximity == null) {
				System.err.println("Distance Null Exception "+ this);
			//	Thread.currentThread().dumpStack();
			return (type==Type.DISTANCE)? 1./epsilon : epsilon;
		//	throw new DistanceFailedExeption(this);
		} 
		double res = proximity.doubleValue() + epsilon;
		
		if(res < 0) if(!(this.getClass().isAssignableFrom(AdjacencyPearsonCorrelation.class))){	
					System.err.println(" Distance Negative Exception: "+ res +" from : "+this);
				//	Thread.currentThread().dumpStack();
				}
		
		cache.put(chacheKey, res);
		cache.put(new Pair<V>(target,source), res); //if symetric
		return res;
	}
	
	public abstract Number computeMeasure(V source, V target);

}
