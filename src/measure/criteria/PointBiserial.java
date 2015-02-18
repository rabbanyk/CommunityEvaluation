package measure.criteria;


public class PointBiserial<V> extends RelativeClusterCriteria<V> {
	{
		withinSetMethod = Within.SUM;
		betweenSetMethod = Between.SUM;
		order =  Order.MINIMIZER;
	}
	
	public PointBiserial() {
		this(true);
	}	
	public PointBiserial(boolean closnessCompatible) {
		super();
		this.closnessCompatible = closnessCompatible;
	}

	// Trying to maximize, this currently does not produce accurate results for OVERLAP or OUTLIERS.
	public double evaluate(){
		//  (M1_M0)/S * sqrt (m1m0/m^2)
		double M1=0, M0=0, S=0, m1=0, m0=0, m=0 , avg= 0, d;
		
		m = N * (N-1) /2;
		
		for (V v1 : allPoints) {
			for (V v2 : allPoints) {
				avg += getDistance(v1, v2);
			}
		}
		avg /= (m*2);
		
		for (V v1 : allPoints) {
			for (V v2 : allPoints) {
				d = getDistance(v1, v2) - avg;
				S += d*d ;
			}
		}
		S /= (m*2);
				
		for (int i = 0; i < k; i++) {
			m1 += N(i) * (N(i)-1) / 2;
			m0 += N(i) * (N - N(i)) / 2;
		}		
		
		for (int i = 0; i < k; i++) {
			M1 += getWithin(i);
			for (int j = 0; j < k; j++) if(i!=j){
				M0 += getBetween(i,j);
			}
		}
		
		M1/=2;
		M0/=2;
		
//		return (isSimInUse()?(M0 - M1):(M1 - M0))/S * Math.sqrt((m1 * m0)/(m*m));
		return (M1 - M0)/S * Math.sqrt((m1 * m0)/(m*m));

	}
	
	public String toString(){
		return "PointBiserial "+closnessCompatible + super.toString();
	}
	public String getName() {
		return "PB" +(closnessCompatible?"' ":" ")+ proximityMeasure.getName();
	}
}
