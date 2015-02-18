package measure.correlation;

public abstract class  Correlation<V> {
	public abstract Number compute (V v1, V v2);
	
	public String toString(){
		return this.getClass().getCanonicalName();
	}
}
