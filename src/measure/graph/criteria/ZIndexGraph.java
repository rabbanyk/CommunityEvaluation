package measure.graph.criteria;

import edu.uci.ics.jung.graph.Graph;
import measure.criteria.RelativeClusterCriteria;

public class ZIndexGraph<V,E>  extends RelativeClusterCriteria<V> {
	{
		//order = Order.MINIMIZER; // check the formula again...
		withinSetMethod = Within.SUM;
	}
	protected Graph<V, E> graph;
	
	
	public boolean isMaximizer() {
		return (order==order.MAXIMIZER);
	}

	public ZIndexGraph(Graph<V, E> graph){
		this(graph,true);
	}
	public ZIndexGraph(Graph<V, E> graph, boolean closnessCompatible) {
		super();
		this.closnessCompatible = closnessCompatible;
		this.graph = graph;
	}
	
	public double evaluate(){
		double wSum = 0, E = 0, V = 0, d, m = 0;
		
		for (int i = 0; i < k; i++) {
			V v[] = (V[])clusters.get(i).toArray();
			for (int v1=0; v1<v.length;v1++ ) 
				for (int v2=v1+1; v2<v.length;v2++ ){
					wSum +=getDistance(v[v1],v[v2]) ;
					m++;
				}
//			m +=  N(i) * (N(i) -1) /2;
		}
				
//		wSum/=2;
		
		double n = 0;
//		for(E e: graph.getEdges()) //if A, we can only iterate over edges
		for (V v1 : allPoints) {
			for (V v2 : allPoints) {
				d = getDistance(v1, v2);
				E += d;
				n++;
			}	
		}
		
		E = (E)/(n);//(N*N);
		
		V = 0.;
		
		for (int i = 0; i < k; i++) {
			V v[] = (V[])clusters.get(i).toArray();
			for (int v1=0; v1<v.length;v1++ ) 
				for (int v2=v1+1; v2<v.length;v2++ ){
					V += Math.pow(getDistance(v[v1], v[v2]) - E ,2);
				}
		}
		
		E*=m;
		
		
		double res = (isSimInUse()?(wSum - E):( E-wSum))/(Math.sqrt(V));
		return res;
	}

	public String toString(){
		return "ZIndexG " +closnessCompatible+ super.toString();
	}
	public String getName() {
		return "ZIndexG"+(closnessCompatible?"' ":" ") + proximityMeasure.getName();
	}
}

