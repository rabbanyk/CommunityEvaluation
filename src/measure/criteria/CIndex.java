package measure.criteria;

import java.util.Vector;

public class CIndex<V> extends RelativeClusterCriteria<V> {
	{
		order = Order.MINIMIZER;
		withinSetMethod = Within.SUM;
	}
	public CIndex(){
		this(true);
	}
	public CIndex(boolean closnessCompatible ) {
		super();
		this.closnessCompatible = closnessCompatible;
	}

	// Trying to minimize, this implementation DOES NOT handle overlap or outliers. 
	public double evaluate(){
		double wSum = 0, min = 0, max = 0, m = 0;
	
		for (int i = 0; i < k; i++) {
			wSum += getWithin(i);
			m +=  N(i) * (N(i) -1);
		}
		
		Vector<Double> mSmallest = new Vector<Double>() , mLargest = new Vector<Double>();
		
		double d;
		for (V v1 : allPoints) {
			for (V v2 : allPoints) {
				d = getDistance(v1, v2);
				
				int i=0;
				while(i < mSmallest.size() && mSmallest.get(i) <= d) i++;
				if(i<m){ // d is in the m smallest distances
					mSmallest.add(i, d);
					if(mSmallest.size() > m) mSmallest.removeElementAt(mSmallest.size()-1);
				} 
				
				i=0;
				while(i < mLargest.size() && mLargest.get(i) >= d) i++;
				if(i<m){ // d is in the m largest distances
					mLargest.add(i, d);
					if(mLargest.size() > m) mLargest.removeElementAt(mLargest.size()-1);
				} 
			}	
		}

		for (Double sd : mSmallest) min += sd.doubleValue();
		for (Double ld : mLargest) 	max += ld.doubleValue();
		
		double res = (wSum - min)/(max - min);
		return res;
	}

	public String toString(){
		return "CIndex" +closnessCompatible+" "+ proximityMeasure.toString();
	}
	
	public String getName() {
		return "CIndex" +(closnessCompatible?"' ":" ") + proximityMeasure.getName();
	}

	
}
