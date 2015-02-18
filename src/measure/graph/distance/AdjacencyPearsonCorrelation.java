package measure.graph.distance;

import java.util.Set;
import measure.graph.GraphProximity;
import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;

public class AdjacencyPearsonCorrelation<V, E> extends GraphProximity<V,E> {
	boolean normalize = false;
	{
		type = Type.SIMILARITY;
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph) {
		this(graph,false);
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph, Transformer<E, ? extends Number> weights) {
		this(graph, weights, false);
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph, boolean normalize) {
		this(graph,null,normalize);
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean normalize) {
		this(graph, weights, false, normalize);
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean selfWeight,boolean normalize) {
		this(graph, weights, selfWeight,false , normalize);
	}
	public AdjacencyPearsonCorrelation(Graph<V, E> graph, Transformer<E, ? extends Number> weights, boolean selfWeight,boolean normalizedWeight,boolean normalize) {
		super(graph, weights, selfWeight,normalizedWeight);
		this.normalize = normalize;
	}
	
	public Number computeMeasure(V source, V target) {
		Set<V> sNeigh = getUnionNeighbours(source, target);
		
		double  sk, tk, n = graph.getVertexCount()+(selfWeight?1:0),
		sumST = 0, sumS= 0, sumT=0,  sum2S=0,  sum2T=0;
		
		for (V k : sNeigh) {
			sk = W(source, k);
			tk = W(target, k);
			
			sumS+=sk; sum2S += sk*sk;
			sumT+=tk; sum2T += tk*tk;
			sumST+= sk*tk;
		}
		
		double res = sumST - (sumS*sumT)/n;
		res /= Math.sqrt((sum2S - (sumS*sumS/n) ) * (sum2T - (sumT*sumT/n)));
		
		if(normalize) res = (res+1)/2;
		return res;
	}
	
	public Number getDistance(V source, V target){
		return 1  - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return (normalize?"Normalized":"")+"PearsonCorrelation"+type+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{"+(normalize?"N":"")+"PC}$";
		else return (normalize?"N":"")+ "PC";
	}

	
//	public HashMap<V, Double> means; //cache for means
//	public HashMap<V, Double> vars; //cache for variances
//	public Number computeMeasureAlt(V source, V target) {
//		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
//		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
//		sNeigh.addAll(tNeigh);
////		sNeigh.remove(source);
//	//	sNeigh.remove(target);
//		
//		double dist = 0, sk, tk,
//		muS= getAverage(source), muT=getAverage(target);
//		
//		for (V k : sNeigh) {
//			sk = W(source, k);
//			tk = W(target, k);
//			
//			dist+= (sk*tk - sk*muT - tk*muS);
//		}
//		dist+= graph.getVertexCount()*muS*muT;
//			
//		dist /= (graph.getVertexCount()*getVariance(source)*getVariance(target));
//		return (dist+1)/2;
//	}
//
//	private double getAverage(V v) {
//		if (means == null) means = new HashMap<V, Double>();
//		
//		if(!means.containsKey(v)){
//			double m = 0;
//			for (E e : graph.getIncidentEdges(v)) {
//				m += W(e);
//			}
//			means.put(v, m / graph.getVertexCount());
//		}
//		return means.get(v);
//	}
//
//	private double getVariance(V v) {
//		if (vars == null) vars = new HashMap<V, Double>();
//		
//		if(!vars.containsKey(v)){
//			double mu = getAverage(v), m=0, w;
//			
//			for (E e : graph.getIncidentEdges(v)) {
//				w = W(e) ;
//				m += ((w*w) - 2*w*mu);
//			}
//			//WTF???  m+= m*m * (graph.getEdgeCount()-graph.getIncidentEdges(v).size());
//			vars.put(v, Math.sqrt(mu*mu + m/ graph.getVertexCount()));
//		}
//		return vars.get(v);
//		
//	}



}
