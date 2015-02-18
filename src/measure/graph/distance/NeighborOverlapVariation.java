package measure.graph.distance;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class NeighborOverlapVariation<V, E> extends GraphProximity<V,E> {
	{
		type = Type.SIMILARITY;
	}
	public NeighborOverlapVariation(Graph<V, E> g) {
		super(g);
	}
	public NeighborOverlapVariation(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public NeighborOverlapVariation(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	public NeighborOverlapVariation(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		super(graph, weights, selfWeight);
	}

	public Number computeMeasure(V source, V target) {
		Set<V> sNeigh = getUnionNeighbours(source, target);
		
		double dist1 = 0, dist2 = 0, sk, tk; 
		for (V k : sNeigh) {
			sk =  W(source, k);//(k == source)? W(source, target): W(source, k);
			tk = W(target, k); //(k == target)? W(source, target): W(target, k);

			dist1 += (sk+tk)*(sk+tk);
			dist2 += (sk-tk)*(sk-tk);
		}
		
		return  (dist1-dist2)/(dist1+dist2);
	}

	//TODO: or division?
	public Number getDistance(V source, V target){
		return 1 - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return "NeighbourOverlapVariation"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{NOV}$";
		else return "NOV";
	}

}
