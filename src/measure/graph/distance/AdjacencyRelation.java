package measure.graph.distance;

import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

/**
 * 
 * @author rabbanyk
 * Adjacency Relation Distance: the structural equivalence of two nodes; 
 * \[	d_{ij}^{A} = \sqrt{\sum_{k \neq{j,i}}{ (A_{ik} - A_{jk})^2}}\]
 */
public class AdjacencyRelation<V, E> extends GraphProximity<V,E> {
	{
		type = Type.DISTANCE;
	}
	
	public AdjacencyRelation(Graph<V, E> g) {
		super(g);
	}
	public AdjacencyRelation(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight) {
		super(graph, weights, selfWeight);
	}
	public AdjacencyRelation(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}
	public AdjacencyRelation(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, boolean selfWeight,
			boolean normalized) {
		super(graph, weights, selfWeight, normalized);
	}
	
	public Number computeMeasure(V source, V target) {
		Set<V> sNeigh = getUnionNeighbours(source, target);
		if(!selfWeight){
			sNeigh.remove(source);
			sNeigh.remove(target);
		}
		
		double dist = 0, sk, tk; 
		for (V k : sNeigh) {
			sk = W(source, k);
			tk = W(target, k);
			
			dist += (sk-tk)*(sk-tk);
		}
		return Math.sqrt(dist);
	}

	
	public Number getSimilarity(V source, V target){
		return graph.getVertexCount() * MaxWeight  - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return "AdjacencyRelation"+ type + ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{AR}$";
		else return "AR";
	}

}
