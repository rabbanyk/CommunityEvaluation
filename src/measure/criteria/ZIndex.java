package measure.criteria;


public class ZIndex<V>  extends RelativeClusterCriteria<V> {
	{
		order = Order.MINIMIZER;
		withinSetMethod = Within.SUM;
	}

	public ZIndex() {
		this(true);
	}
	public ZIndex(boolean closnessCompatible) {
		super();
		this.closnessCompatible = closnessCompatible;
	}
	
	public double evaluate(){
		double wSum = 0, E = 0, V = 0, d, m = 0;
			
		for (int i = 0; i < k; i++) {
			wSum += getWithin(i);
			m +=  N(i) * (N(i) -1);
		}
				
//		wSum/=2;
		
		for (V v1 : allPoints) {
			for (V v2 : allPoints) {
				d = getDistance(v1, v2);
				E += d;
			}	
		}
		
		E = (E)/(N*N);
		
		V = 0.;
		
		for (int i = 0; i < k; i++) 
		//	for (int j = 0; j < k; j++) if(i<j)
				for (V v1: clusters.get(i)) 
					for (V v2: clusters.get(i)) if(v1!=v2){
						V += Math.pow(getDistance(v1, v2) - E ,2);
					}
		
		E*=m;
		
		
		double res = (wSum - E)/(Math.sqrt(V));
		return res;
	}

	public String toString(){
		return "ZIndex " +closnessCompatible+ super.toString();
	}
	public String getName() {
		return "ZIndex"+(closnessCompatible?"' ":" ") + proximityMeasure.getName();
	}
}

