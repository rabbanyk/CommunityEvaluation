package measure.cluster.agreement.partitioning.generalization;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.util.Pair;

public class GAM<V> extends PartiotioningAgreement<V>{
	
	Transformer<Double, Double> phi ;
	Transformer< Pair<Set<V>>, Double> etha;
	Transformer<Double, Double> linear ;
	double Epsilon = 0.000000000001;
	public static enum Type { RI, VI , X2, NVI, ELSE};
	Type type;
	
	public GAM(){
		phi = getRIPhi(); //Default
		etha = getIntersectionEtha(); //Default
		linear = getRILinear(); //Default
	}
	
	public GAM(Type type) {
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

	public GAM(Transformer<Double, Double> phi) {
		this.phi = phi;
		if(phi == null) phi = getRIPhi(); //Default
		etha = getIntersectionEtha(); //Default
//		linear = getUniformLinear(); //Default
		linear = getRILinear(); //Default

	}

	public GAM(Transformer<Double, Double> phi, Transformer<Pair<Set<V>>, Double> etha) {
		this.phi = phi;
		this.etha = etha;
		if(phi == null) phi = getRIPhi(); //Default
		if (etha == null) etha = getIntersectionEtha(); //Default
//		linear = getUniformLinear(); //Original \in [-2 .. 0]
		linear = getRILinear(); //Default \in [-1 .. 1]

	}

	
	protected double r(Vector<Set<V>> U, Vector<Set<V>> V) {
		double n = 0, d = 0;
		
//		//TODO: testing remove these
//		int m=0;
//		for (int i = 0; i < U.size(); i++) {
//			Set<V> u = U.get(i);
//			m+=u.size();
//		}//System.err.println(m);
//		double[][] e = new double[U.size()][V.size()];
//		for (int i = 0; i < U.size(); i++) {
//			Set<V> u = U.get(i);
//			for (int j = 0; j < V.size(); j++) {
//				Set<V> v = V.get(j);
//				e[i][j] = (u.size() * v.size() *1.) / m;
//			}
//		} 
//		//TODO: testing remove till here
//		for (Set<V> v : V) {
		for (int j = 0; j < V.size(); j++) {
			Set<V> v = V.get(j);
			double sum = 0, sumPhi = 0, tmp;
			//for (Set<V> u : U) {
			for (int i = 0; i < U.size(); i++) {
				Set<V> u = U.get(i);
				
				tmp = etha.transform(new Pair<Set<V>>(u,v));
					
//				tmp -= e[i][j];
//				if(tmp<0) tmp *=-1;
////				tmp /= Math.sqrt(e[i][j]);
//				tmp /= (1-u.size()*1./m) * (1-v.size()*1./m); 
				
				d+=tmp;
				sum += tmp; 
				sumPhi += phi.transform(tmp);
			}
			
			n += sumPhi - phi.transform(sum);
			
		}
		d = phi.transform(d);
		
		double res = n!=0?n/d:0;
//		if (n!=0)
//			n/=d;
		
		if(d==0 && n!=0)
			System.err.println("-->> GAM_"+" ERROR:: is zero, n not zero, error" );
		if(res> Epsilon || res<-1-Epsilon){
			System.err.println(">> GAM_"+" ERROR::out of range! " +res+" n: "+ n + " d: " + d   );
			
		}
		return res ;//d==0?-1:n/d;
	}
	
	
	public double getAgreement(Vector<Set<V>> U, Vector<Set<V>> V) {
		return linear.transform(r(U,V) + r(V,U));
	}
	public static Transformer<Double, Double> getX2(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return x*(x);}
		};
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
			public Double transform(Double x) {return (x==0 ? 0 : ( x * Math.log(x)));}
		};
	}
	public static Transformer<Double, Double> getNegativeXlogx(){
		return new Transformer<Double, Double>() {
			public Double transform(Double x) {return (x==0 ? 0 :( (x+1)* Math.log(x+1)));}
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
	
//	public static<V> Transformer<Pair<Set<V>>, Double> getAdjustedIntersectionEtha(){
//		return  new Transformer<Pair<Set<V>>, Double>() {
//			public Double transform(Pair<Set<V>> uv) {
//				Set<V> z = new HashSet<V>(uv.getFirst());
//				z.retainAll(uv.getSecond());
//				double m = 17;
//				double tmp = z.size()*1.;
//				tmp /= (uv.getFirst().size() * uv.getSecond().size() *1.)/(m);
//				//if(tmp<0) tmp *=-1;
////				tmp*=tmp;
////				tmp /= (uv.getFirst().size()*1./m) * ( uv.getSecond().size()*1./m); 				
//				return tmp;
//			}
//		};
//	}
	public String toString(){
		return "GAM_" + type;
	}
}
