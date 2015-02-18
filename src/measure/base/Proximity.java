package measure.base;


public abstract class Proximity <V> {
	public enum Type {DISTANCE, SIMILARITY};
	protected Type type ;//= Type.DISTANCE;
	public static double epsilon = 1.0E-9; //.0000001 ;

	
	public abstract Number getProximity(V v1, V v2);
	public abstract String getName();
	
	public boolean isSimilarity(){
		return type==Type.SIMILARITY;
	}
	public Number getDistance(V source, V target){
		if(type==Type.DISTANCE) return getProximity(source, target);
		else return 1/getProximity(source, target).doubleValue();
	}
	public Number getSimilarity(V source, V target){
		if(type==Type.DISTANCE) return 1/getProximity(source, target).doubleValue();
		else return getProximity(source, target); 
	}
	
}
