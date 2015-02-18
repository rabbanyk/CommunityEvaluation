package measure.base;

import java.util.Set;


/**
 * @author rabbanyk
 *
 * @param <V>
 * 
 * <pre> argmin_m d(x,m) </pre>
 */

public class Medoid<V> implements Centroid<V>{
	Proximity<V>  proximity;
	
	public Medoid(Proximity<V>  distance){
		this.proximity = distance; 
	}

	
	public V findCentroid(Set<V> X) {
		V center = null;
		double sum = 0, sumd; 
		
		for (V x : X) {
			
			sumd = 0;
			for (V y : X) if(x!=y){
				sumd += proximity.getProximity(y, x).doubleValue();
			}
			
			if(center == null || (proximity.isSimilarity()?(sum < sumd):(sum > sumd))){
				center = x;
				sum = sumd;
			}
		
		}
		
		return center;
	}
	
	public String toString(){
		return "MedoidBasedOn "+proximity;
	}

}
