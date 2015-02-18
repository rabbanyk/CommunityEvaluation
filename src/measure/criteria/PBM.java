package measure.criteria;


public class PBM<V> extends RelativeClusterCriteria<V> {
	{
		withinSetMethod = Within.SUM2C;
		betweenSetMethod = Between.CENTROID;
	}
	public PBM() {
		this(false);
	}
	public PBM(boolean closnessCompatible) {
		super();
		this.closnessCompatible = closnessCompatible;
	}
	// Want to maximize
	public double evaluate(){
		double wSum = 0;
		Double maxC = null;
		//V graphCentroid = getCentroid();
		for (int i = 0; i < k; i++) {
			wSum += getWithin(i);
//			for (V v : getCluster(i)){
//				bSum += getDistance(v, graphCentroid);	
//			}
		}
		for (int i = 0; i < k; i++) {
			for (int j = i+1; j < k; j++) {
				double tmp = getBetween(i, j);
				if((maxC==null) || (isSimInUse()?(tmp<maxC): (tmp>maxC))) maxC = tmp;
			}
		}
		return  ( maxC / wSum) / k ;
	}
	
	public String toString(){
		return "PBM " +closnessCompatible + super.toString();
	}
	
	public String getName() {
		return "PBM"+(closnessCompatible?"' ":" ") + proximityMeasure.getName();
	}
}
