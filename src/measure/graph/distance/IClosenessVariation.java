package measure.graph.distance;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import measure.graph.GraphProximity;

import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;

public class IClosenessVariation<V, E> extends GraphProximity<V,E> {
	int neighborhoodThreshold = 2;//, sourceNeighborhoodThreshold = 1;
	HashMap<V, HashMap<V,Double>> neighborhoods;
	{
		type = Type.SIMILARITY;		
	}
	//TODO: normalize it !!!
	public IClosenessVariation(Graph<V, E> g) {
		this(g, null);
	}
	
	public IClosenessVariation(Graph<V, E> g, Transformer<E, ? extends Number> weight) {	
		this( 2, g, weight);
	}
	
	public IClosenessVariation(int neighborhoodThreshold, Graph<V, E> g, Transformer<E, ? extends Number> weight) {
		super(g, weight);		
		this.neighborhoods = new HashMap<V, HashMap<V,Double>>();
		this.neighborhoodThreshold = neighborhoodThreshold;
	}

	public HashMap<V,Double> getNeighbors(V n) {
		if (neighborhoods.containsKey(n))
			return neighborhoods.get(n);
		HashMap<V, Double> neighbors = new HashMap<V,Double>(), newly_neighbours;
		neighbors.put(n, graph.degree(n) * 1.0);//TODO:shouldn't this be weighted sum?

		HashSet<E> explored = new HashSet<E>();
		HashSet<E> newly_explored;// = new HashSet<E>();

		newly_neighbours = new HashMap<V, Double>();
		newly_explored = new HashSet<E>();

		for (int i = 1; i <= neighborhoodThreshold ; i++) {			
			for (V v : neighbors.keySet()){
				double weightedDegree = 0;
				for (E e : graph.getIncidentEdges(v))
					weightedDegree += W(e);//weights.transform(e).doubleValue();
				
				for (E e : graph.getIncidentEdges(v))
					if(!explored.contains(e)) {
						newly_explored.add(e);
						for (V m :  graph.getIncidentVertices(e)) if(m!=v){
								if(newly_neighbours.containsKey(m))
									newly_neighbours.put(m, newly_neighbours.get(m) + neighbors.get(v) * W(e)/ weightedDegree);//(graph.degree(v)));
								else
									newly_neighbours.put(m, neighbors.get(v)* W(e) /weightedDegree);//(graph.degree(v)));
						}
					}
			}
			neighbors.putAll(newly_neighbours);
			explored.addAll(newly_explored);
		}
	
		neighbors.remove(n);
		neighborhoods.put(n, neighbors);
		return neighbors;
	}
	
	private double ns(HashMap<V, Double> scoredNeighbourhood, V v){
		Double score = scoredNeighbourhood.get(v);
		return score==null?0:score.doubleValue();
	}
	
	public double iCloseness(V v1, V v2){
		HashMap<V, Double> Nv1 = getNeighbors(v1) , Nv2 = getNeighbors(v2);
		Set<V> union = new HashSet<V>(Nv1.keySet());
		union.addAll(Nv2.keySet());
		union.add(v1);union.add(v2);
		
		double res = 0 , dom = 0;
		for (V v: union) {
				res += (ns(Nv1,v)-ns(Nv2,v))*(ns(Nv1,v)-ns(Nv2,v));
				dom += (ns(Nv1,v)+ns(Nv2,v))*(ns(Nv1,v)+ns(Nv2,v));;
		}
//		res +=  ns(Nv1,v2)*ns(Nv2,v1);
//		dom += ns(Nv1,v2)*ns(Nv1,v2) + ns(Nv2,v1)*ns(Nv2,v1);

		return ((dom - res) / (dom + res));
	}
	
	public Double computeMeasure(V source, V target) {
		double ic= iCloseness(source, target);
		return ic ; 
	}
	
	public Number getDistance(V source, V target){
		return 1 - getProximity(source, target).doubleValue();
	}
	
	public String toString(){
		return "iClosenessVariation"+type+neighborhoodThreshold;
	}
	@Override
	public String getName() {
		return "ICV"+neighborhoodThreshold;
	}

}
