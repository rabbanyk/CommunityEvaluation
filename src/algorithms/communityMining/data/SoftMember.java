package algorithms.communityMining.data;

public class SoftMember<V> {
	V v;
	double strenght;
	
	public SoftMember(V v, double strenght) {
		super();
		this.v = v;
		this.strenght = strenght;
	}
	public SoftMember(V v) {
		this(v,1);
	}
	public V getV() {
		return v;
	}
	public void setV(V v) {
		this.v = v;
	}
	public double getStrenght() {
		return strenght;
	}
	public void setStrenght(double strenght) {
		this.strenght = strenght;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof SoftMember)
			return v.equals(((SoftMember) obj).v);
		else return v.equals(obj);
	}
	
	@Override
	public int hashCode() {
		return v.hashCode();
	}
	
	
}
