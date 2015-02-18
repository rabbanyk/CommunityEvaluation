package measure.graph.distance;

import java.util.Set;
import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class TopologicalOverlap<V, E> extends GraphProximity<V,E> {
	{
		type = Type.SIMILARITY;
	}
	public TopologicalOverlap(Graph<V, E> g) {
		super(g);
	}
	public TopologicalOverlap(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public TopologicalOverlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	public TopologicalOverlap(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		super(graph, weights, selfWeight);
	}
	public Number computeMeasure(V source, V target) {
		Set<V> cNeigh = getUnionNeighbours(source, target);
		
		double ds = 0, dt = 0, sk, tk, nom = 0; 
		for (V k : cNeigh) {
			sk = W(source, k);
			ds += sk*sk;
			
			tk = W(target, k);
			dt += tk*tk;
			
			nom += sk*tk;
		}
		
		if(!selfWeight)
			nom +=  W(source, target)*W(source, target);

		double dom =(Math.min(ds, dt));
	
		double sim = nom/dom;
		
		if(sim >2 || sim < 0 ) System.err.println(sim+ " ( " +source+" , "+target+" ) " + " : "+nom+" / "+ dom + "  " + cNeigh.size());
		return sim;
	}

	public Number getDistance(V source, V target){
		return 1  - getProximity(source, target).doubleValue();
	}

	public String toString(){
		return "TopologicalOverlap"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{TO}$";
		else return "TO";
	}

}
