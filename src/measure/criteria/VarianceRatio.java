package measure.criteria;

/**
 * 
 * @author rabbanyk
 *
 * Variance Ratio criterion adopted for graphs
 * 
 */
public class VarianceRatio<V> extends RelativeClusterCriteria<V> {
	{
		order = Order.MAXIMIZER;
		withinSetMethod = Within.SUM2C;
	}
	public VarianceRatio() {
		this(true);
	}
	public VarianceRatio(boolean closenessCompatible) {
		super();
		this.closnessCompatible = closenessCompatible;
	}
	public boolean isMaximizer() {
		return (order==order.MAXIMIZER);
	}
	// Higher value is better.
	protected double evaluate(){		
		double wSum = 0;
		double bSum = 0;
		
		V graphCentroid = getCentroid();
		
		for (int i = 0; i < k; i++) {
			wSum += getWithin(i);
			bSum += N(i) * getDistance(getCentroid(i), graphCentroid);	
		}
	
		double res = (bSum / wSum) ;
	//	if(wSum == 0) throw new CriteriaFailedException(this);
	//	return ((isSimInUse()? (1/res):res)* ((N - k)/(k - 1)));
		//return (res * ((isSimInUse()?((k-1)*1./(N - k)):((N - k)*1./(k - 1)))));
		return ((isSimInUse()? (1/res):res) * (N - k)*1. / (k - 1));
	}
	
	public String toString(){
		return "VarianceRatio "+closnessCompatible +" "+ super.toString();
	}
	public String getName() {
		return "VR" +(closnessCompatible?"' ":" ")+ proximityMeasure.getName();
	}
}
