package measure.criteria;

public class ZIndexOld<V>  extends RelativeClusterCriteria<V> {
	{
		//order = Order.MINIMIZER; // check the formula again...
		withinSetMethod = Within.SUM;
	}
	public boolean isMaximizer() {
		return (order==order.MAXIMIZER);
	}

	public ZIndexOld() {
		this(true);
	}
	public ZIndexOld(boolean closnessCompatible) {
		super();
		this.closnessCompatible = closnessCompatible;
	}
	
	public double evaluate(){
		double wSum = 0, E = 0, V = 0, d, a1 = 0, a2 = 0, a3 = 0 , tmp;
	
		for (int i = 0; i < k; i++) {
			wSum += getWithin(i);
		}
				
		wSum/=2;
		
		V[] v= (V[])allPoints.toArray();
		
		for (int i=0; i<v.length;i++) {
			tmp =0;
			for (int j=i+1; j<v.length;j++) {
				d = getDistance(v[i], v[j]);
				E += d;
				tmp += d;
				a3 += d*d;
			}	
			a2+= tmp * tmp;
		}
		
		a1 = E*E;
		
		E/= N;
		
		V = a1 /(N*(N-1) )- (2*a2 )/(N*(N)) - a1/(N*N) + a3/N; //Different from the paper!!! a2/N*(N-1) should be!!
		
		double res = (isSimInUse()?(wSum - E):( E-wSum))/(Math.sqrt(V));
		
//		if((res+"").equals("NaN")){
//			System.err.println(res+" -- a:   "+a1 +"  "+a2 +" "+a3);
//			System.err.println(res+" -----     "+wSum +"  "+E +" "+V);
//		}

		return res;
	}

	public String toString(){
		return "ZIndexOld " +closnessCompatible+ super.toString();
	}
	public String getName() {
		return "ZIndexOld"+(closnessCompatible?"' ":" ") + proximityMeasure.getName();
	}
}

