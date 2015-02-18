package measure.graph.distance;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

/**
 * Node neighbourhood distance: the ratio of common neighbours between two nodes;
 * \[	d_{ij}^{\aleph} = 1 - \frac{|\aleph_i \cap \aleph_j|}{|\aleph_i \cup \aleph_j|} \] 
 * where $\aleph_i$ is the set of nodes directly connected to $i$.
 * 
 * since $|\aleph_i|=\sum_{k\neq{i}}{A_{ik}}$ formula could be re-written as:
 * \[d_{ij}^{\aleph} = 1 - \frac{a-b}{a+b} \]
 * where a = \sqrt{\sum_{k \neq{j,i}}{ (A_{ik} + A_{jk})^2}} and b = \sqrt{\sum_{k \neq{j,i}}{ (A_{ik} - A_{jk})^2}}
 * which could be then generalized for weighted networks by simply substituting A with W
 */
public class NeighborOverlapOld<V, E> extends GraphProximity<V,E> {
	{
		type = Type.DISTANCE;
	}
	public NeighborOverlapOld(Graph<V, E> g) {
		super(g);
	}
	public NeighborOverlapOld(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		super(g, weights);
	}

	public Number computeMeasure(V source, V target) {
	
		Set<V> sNeigh = new HashSet<V>(graph.getNeighbors(source));
		Set<V> tNeigh = new HashSet<V>(graph.getNeighbors(target));
		sNeigh.addAll(tNeigh);
//		selfWeight = true;
		double dist1 = 0, dist2 = 0, sk, tk; 
		for (V k : sNeigh) {
			if(k!=source)
				sk = W(source, k);
			else sk = W(source,target);
			if(k!=target)
				tk = W(target, k);
			else tk = W(target, source);
//			sk = W(source, k);
//			tk = W(target, k);
			
			dist1 += (sk+tk)*(sk+tk);
			dist2 += (sk-tk)*(sk-tk);
		}
		
		dist1 = Math.sqrt(dist1);
		dist2 = Math.sqrt(dist2);

		return 1-(dist1-dist2)/(dist1+dist2);
	}
	
	public Number getSimilarity(V source, V target){
		return 1 - getProximity(source, target).doubleValue();
	}

	public String toString(){
		return "NeighbourOverlapDistanceOLD";
	}
	@Override
	public String getName() {
		return "NODOLD";
	}

}
