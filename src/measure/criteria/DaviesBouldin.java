package measure.criteria;


public class DaviesBouldin<V> extends RelativeClusterCriteria<V> {
	{
		withinSetMethod = Within.AVG2C;
		betweenSetMethod =  Between.CENTROID;
		order = Order.MINIMIZER;
	}
	public DaviesBouldin() {
		this(false);
	}
	public DaviesBouldin(boolean closenessCompatible) {
		super();
		this.closnessCompatible = closenessCompatible;
	}
	
//	public boolean isMaximizer() {
//		return (order==order.MAXIMIZER);
//	}
	public double evaluate(){
		double dbSum = 0;  double d;
	
		for (int i = 0; i < k; i++) {
			Double worst = null;
			for (int j = 0; j < k; j++) 
				if (i!=j){
					
			//		d = isSimInUse()? (getBetween(i, j)/(getWithin(i) + getWithin(j) )):((getWithin(i) + getWithin(j) ) / getBetween(i, j));
					//	if(worst ==null || (d > worst)) worst = d;

					d = ((getWithin(i) + getWithin(j) ) / getBetween(i, j));
					if(worst ==null || (isSimInUse()?(d < worst):(d > worst))) worst = d;
				}	
			dbSum += worst;
		}	
		
		return dbSum / k ;
		
		/*
		 * change minimization criteria to maximization
		 */
	//	return flipValues(value); 1-s should not work, what else?
	}

	public String toString(){
		return "DaviesBouldin "+closnessCompatible  + super.toString();
	}
	public String getName() {
		return "DB" +(closnessCompatible?"' ":" ")+ proximityMeasure.getName();
	}

}
