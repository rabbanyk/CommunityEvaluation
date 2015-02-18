package measure.cluster.agreement.partitioning.generalization;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.util.Pair;

public class GAMA<V> extends PartiotioningAgreement<V>{
	
	Transformer<Double, Double> phi ;
	Transformer< Pair<Set<V>>, Double> etha;
	Transformer<Double, Double> linear ;
	
	public static enum Type { RI, VI ,ELSE};
	Type type;
	
	public GAMA(){
		phi = getRIPhi(); //Default
		etha = getIntersectionEtha(); //Default
		linear = getRILinear(); //Default
	}
	
	public GAMA(Type type) {
		this.type = type;
		
		etha = getIntersectionEtha(); 

		switch (type) {
		case RI:
			phi = getRIPhi(); 
			linear = getRILinear();
			break;
		case VI:
			phi = getVIPhi();
			linear = getVILinear();
			break;
		}
	}

	public GAMA(Transformer<Double, Double> phi) {
		this.phi = phi;
		if(phi == null) phi = getRIPhi(); //Default
		etha = getIntersectionEtha(); //Default
		linear = getUniformLinear(); //Default
	}

	public GAMA(Transformer<Double, Double> phi, Transformer<Pair<Set<V>>, Double> etha) {
		this.phi = phi;
		this.etha = etha;
		if(phi == null) phi = getRIPhi(); //Default
		if (etha == null) etha = getIntersectionEtha(); //Default
		linear = getUniformLinear(); //Default
	}

	
	protected double r(Vector<Set<V>> U, Vector<Set<V>> V) {
		double n = 0, d = 0;

		for (int j = 0; j < V.size(); j++) {
			Set<V> v = V.get(j);
			double sum = 0, sumPhi = 0, tmp;
			for (int i = 0; i < U.size(); i++) {
				Set<V> u = U.get(i);
				
				tmp = etha.transform(new Pair<Set<V>>(u,v));
				
				d+=tmp;
				sum += tmp; 
				sumPhi += phi.transform(tmp);
			}
			
			n += sumPhi - phi.transform(sum);
			
		}
		d = phi.transform(d);
		if(d==0)
			System.err.println("--------->> GAMA shall not be zero, error" );
		if(d> 1 || d<-1)
			System.err.println("--------->> GAMA out of range!, error " + d );
		return d==0?0:n/d;
	}
	
	
	public double getAgreement(Vector<Set<V>> U, Vector<Set<V>> V) {
		return linear.transform(r(U,V) + r(V,U));
	}

	public static Transformer<Double, Double> getRIPhi(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return x*(x-1);}
		};
	}
	public static Transformer<Double, Double> getRILinear(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return x+1;}
		};
	}
	public static Transformer<Double, Double> getVIPhi(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return (x==0 ? 0 : x * Math.log(x));}
		};
	}
	public static Transformer<Double, Double> getVILinear(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return x*-1;}
		};
	}
	
	public static Transformer<Double, Double> getUniformLinear(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return x;}
		};
	}
	public static<V>  Transformer<Pair<Set<V>>, Double> getIntersectionEtha(){
		return  new Transformer<Pair<Set<V>>, Double>() {
			public Double transform(Pair<Set<V>> uv) {
				Set<V> z = new HashSet<V>(uv.getFirst());
				z.retainAll(uv.getSecond());
				return z.size()*1.;
			}
		};
	}
	
	public static<V> Transformer<Pair<Set<V>>, Double> getAdjustedIntersectionEtha(){
		return  new Transformer<Pair<Set<V>>, Double>() {
			public Double transform(Pair<Set<V>> uv) {
				Set<V> z = new HashSet<V>(uv.getFirst());
				z.retainAll(uv.getSecond());
				double m = 17;
				double tmp = z.size()*1.;
				tmp /= (uv.getFirst().size() * uv.getSecond().size() *1.)/(m);
				//if(tmp<0) tmp *=-1;
//				tmp*=tmp;
//				tmp /= (uv.getFirst().size()*1./m) * ( uv.getSecond().size()*1./m); 				
				return tmp;
			}
		};
	}
	public String toString(){
		return "GAM_" + type;
	}
}
