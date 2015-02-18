package measure.graph.distance;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class NeighborOverlap<V, E> extends GraphProximity<V,E> {
	{
		type = Type.SIMILARITY;
	}
	public NeighborOverlap(Graph<V, E> g) {
		super(g);
	}
	public NeighborOverlap(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public NeighborOverlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	public NeighborOverlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		super(graph, weights, selfWeight);
	}
	public Number computeMeasure(V source, V target) {
		Set<V> sNeigh = getUnionNeighbours(source, target);
		
		double  sk, tk,	sumST = 0, 
		sum2S=0,  sum2T=0;
		
		for (V k : sNeigh) {
			sk = W(source, k);
			tk = W(target, k);
			
			sum2S += sk*sk;
			sum2T += tk*tk;
			sumST+= sk*tk;
		}
		
		return (sumST / (sum2S + sum2T - sumST));
	}

	public Number getDistance(V source, V target){
		return 1  - getProximity(source, target).doubleValue();
	}
	public String toString(){
		return "NeighbourOverlap"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{NO}$";
		else return "NO";
	}
	
//	public Number computeMeasureAlt(V source, V target) {
//		
//		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
//		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
//		
//		Set<V> union = new HashSet<V>(sNeigh), common = new HashSet<V>(sNeigh);
//		union.addAll(tNeigh);
//		common.retainAll(tNeigh);
//		
//		return (1.*common.size())/union.size();
//		}
}
