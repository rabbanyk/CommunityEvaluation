package measure.graph.distance;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class Overlap<V, E> extends GraphProximity<V,E> {
	{
		type = Type.SIMILARITY;
	}
	public Overlap(Graph<V, E> g) {
		super(g);
	}
	public Overlap(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public Overlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	public Overlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		super(graph, weights, selfWeight);
	}
	public Number computeMeasure(V source, V target) {
		Set<V> sNeigh = getUnionNeighbours(source, target);

		double ds = 0, dt = 0, sk, tk, sim = 0; 
		for (V k : sNeigh) {
			sk = W(source, k);
			ds += sk*sk;
			
			tk = W(target, k);
			dt += tk*tk;
			
			sim += sk*tk;
		}
		
		sim +=  W(source, target)*W(source, target);

		return (sim);///(Math.min(ds, dt));
	}

	public Number getDistance(V source, V target){
		return 1./ getProximity(source, target).doubleValue();
	}
	public String toString(){
		return "Overlap"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{O}$";
		else return "O";
	}

}
