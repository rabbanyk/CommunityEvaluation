package measure.criteria;


// the dissimilarity measure should be of ratio scale 
public class Silhouette<V> extends RelativeClusterCriteria<V> {

	boolean alternative ;
	

	public Silhouette(){
		this(false);
	}
	public Silhouette(boolean closenessCompatible){
		this(Between.AVERAGE,closenessCompatible);
	}
	
	// node2SetMethod = AVG is the original Silhouette while = AVG2C is the simplified Silhouette, 
	public Silhouette(Between node2SetMethod,boolean closenessCompatible) {
		this(node2SetMethod,false, closenessCompatible);
	}
	
	public Silhouette(Between node2SetMethod, boolean alternative,boolean closenessCompatible) {
		super();
		this.betweenSetMethod = node2SetMethod;
		this.alternative = alternative;
		this.closnessCompatible = closenessCompatible;
	}
	
	// Thus a high Silhouette metric is desired. 
	protected double evaluate(){
		double res = 0;
		
		for (int i = 0; i < k; i++) 
			if(N(i)<=1) continue; //it is in the original definition that if it is a single node cluster s is zero
			else {
			for (V v : getCluster(i)){
				double a = getDistance(v, i) , d;
				
				Double b = null;
				for (int j = 0; j < k; j++) if(i!=j){
					d = getDistance(v, j);
					if(b==null || (isSimInUse()? (d > b) :(d < b)) ) b = d;
				}
				
				if(alternative)
					res += b/a ;
				else 
					res += (b - a) /Math.max(a, b);// (isSimInUse()?Math.min(a, b):Math.max(a, b));
			}
		}
		
		return res / N;
		
	}
	
	public String toString(){
		return (alternative?"Altarnative Silhouette ":""+" Silhouette ")+closnessCompatible  + betweenSetMethod +  " " + super.toString();
	}
	public String getName() {
		return (alternative?"ASWC":""+"SWC") + (betweenSetMethod.ordinal()-Between.AVERAGE.ordinal())+(closnessCompatible?"'":"")+ " "+ proximityMeasure.getName();

	}
	
}
