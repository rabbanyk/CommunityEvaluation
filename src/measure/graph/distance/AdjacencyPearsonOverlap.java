package measure.graph.distance;

import java.util.Set;
import measure.graph.GraphProximity;

import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;

public class AdjacencyPearsonOverlap<V, E> extends GraphProximity<V,E> {
	boolean normalize = true;
	{
		type = Type.SIMILARITY;
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph) {
		this(graph,true);
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph, Transformer<E, ? extends Number> weights) {
		this(graph, weights, true);
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph, boolean normalize) {
		this(graph,null,normalize);
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean normalize) {
		this(graph, weights, true, normalize);
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean selfWeight,boolean normalize) {
		this(graph, weights, selfWeight,true , normalize);
	}
	public AdjacencyPearsonOverlap(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean selfWeight,boolean normalizedWeight,boolean normalize) {
		super(graph, weights, selfWeight,normalizedWeight);
		this.normalize = normalize;
	}
	
	
	public Number computeMeasure(V source, V target) {
		//sNeigh.addAll(graph.getVertices());//if we iterate over all, this is same as PC
		Set<V> sNeigh = getUnionNeighbours(source, target);
		
		double  sk, tk,  n = graph.getVertexCount()+(selfWeight?1:0), sumS= 0, sumT=0,  a=0,b=0,c=0;
		double m  = n*n;
		
		for (V k : sNeigh) {
			sk =  W(source, k);
			tk =  W(target, k);
			sumS+=sk;
			sumT+=tk;
		}
		
		for (V k : sNeigh) {
			sk =  W(source, k);
			tk =  W(target, k);
			
			a+= sk*tk - sumS*sumT/m;
			b+= sk*sk - sumS*sumS/m;
			c+= tk*tk - sumT*sumT/m;
		}
		
		double res = a;
		res /= Math.sqrt(b*c);
		
		if(normalize) res = (res+1)/2;
		
		return res;
	}
	
	public Number getDistance(V source, V target){
		return 1  - getProximity(source, target).doubleValue();
		//TODO: or 1./ ?
	}
	
	public String toString(){
		return (normalize?"Normalized":"")+"PearsonOverlap"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{"+(normalize?"N":"")+"PO}$";
		else return (normalize?"N":"")+ "PO";
	}



}
