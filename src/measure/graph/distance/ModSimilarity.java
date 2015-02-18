package measure.graph.distance;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class ModSimilarity<V, E> extends GraphProximity<V,E> {
	public enum Norm{MINUS,DIVIDE};
	Norm norm = Norm.DIVIDE;
	{
		type = Type.SIMILARITY;
	}

	public ModSimilarity(Graph<V, E> g) {
		super(g);
	}
	public ModSimilarity(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public ModSimilarity(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	public ModSimilarity(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		this(graph, weights, selfWeight,Norm.DIVIDE);
	}
	public ModSimilarity(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight, Norm norm) {
		super(graph, weights, selfWeight);
		this.norm = norm;
	}

	public Number computeMeasure(V source, V target) {
		double sim = 0;
		sim +=  W(source, target);
		Set<V> sNeigh = getUnionNeighbours(source, target);
//		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
//		sNeigh.add(source);
//		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
//		tNeigh.add(target);

		double sk=0, tk=0;
		for (V k : sNeigh) sk += W(source, k);
		for (V k : sNeigh) tk += W(target, k);
		
		double tmp = sk*tk/((SumWeight+(selfWeight?2*MaxWeight:0))*4); 
	//	System.err.println(tmp + "   " + SumWeight +"  " + sk*tk/(SumWeight*4) );
		if(norm==Norm.DIVIDE) 
			sim = sim / tmp;//sim*(sk+tk)/(sk*tk) 
		else 
			sim = (sim - tmp+1)/2; 
		
		return sim;
	}
	public Number getDistance(V source, V target){
		return 1 - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return "Mod"+type+ norm+((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{M"+(norm==Norm.DIVIDE?"D":"M")+"}$";
		else return "M"+(norm==Norm.DIVIDE?"D":"");
	}

}
