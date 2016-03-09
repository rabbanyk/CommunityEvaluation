package measure.cluster.overlapping;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;
import java.util.function.Function;
import org.apache.commons.collections15.Transformer;

import data.Pair;
import measure.cluster.agreement.ClusteringAgreement;

public class OCI<d > extends ClusteringAgreement<d> {

	public static Function<Double, Double> x2 = x -> x*x;
	public static Function<Double, Double> xlogx = x -> (x==0 ? 0 : ( x * Math.log(x)) );
	Function<Double, Double> phi = x2;

	Transformer< Pair<Set<d>,Set<d>>, Double> etha =   new Transformer<Pair<Set<d>,Set<d>>, Double>() {
		public Double transform(Pair<Set<d>,Set<d>> uv) {
			Set<d> u =uv.first , v=uv.second; 
			if (uv.first.size()>uv.second.size()) {u=uv.second ; v=uv.first; } 
			double res =0;
			for (d i: u)
				if(v.contains(i))
					res += 1;
			return res;
		}
	};
	Transformer<Set<d>, Double> metha = new Transformer<Set<d>, Double>() {
		public Double transform(Set<d> u) {
			double res =0;
			for (d i: u) res += 1;
			return res;
		}
	};;


	public OCI(Function<Double, Double> phi){
		this.phi = phi;
	}
	@Override
	public double getAgreement(Vector<Set<d>> U, Vector<Set<d>> V) {
		Set<d> all=new HashSet<>();
		for (Set<d> u : U) 	all.addAll(u);
		double n = all.size();
	
		return measure(n, U,V);
	}

	protected double O(Vector<Set<d>> U, Vector<Set<d>> V){
		double res =0;
		for (Set<d> u :U) for (Set<d> v :V){
			res+=phi.apply(o(u,v));
		}
		return res;
	}
	protected double E(double n, Vector<Set<d>> U, Vector<Set<d>> V){
		double res =0, ou=0;
		double[] oV = new double[V.size()];
		for (int i = 0; i < oV.length; i++) oV[i]=o(V.get(i));

		for (int i = 0; i < U.size(); i++){
			ou = o(U.get(i));
			for (int j = 0; j < oV.length; j++) res+=phi.apply(ou*oV[j]/n);
		}
		return res;
	}
	
	protected double measure(double n, Vector<Set<d>> U, Vector<Set<d>> V){
		double Euv =  E(n, U,V);
//		System.err.println("O(U, V): "+O(U, V)+"  O(U, U) " +  O(U, U) + " O(V, V) "+O(V, V) +"  Euv: "+Euv);
		return ( O(U, V) - Euv)/(0.5* ( O(U, U)+O(V, V)) - Euv);
	}
	
	protected double o(Set<d> u){
		return metha.transform(u);
	} 
	protected double o(Set<d> u, Set<d> v){
		return etha.transform(new Pair<Set<d>,Set<d>>(u,v));
	} 
	

	
	public String toString() {
		if(phi==x2) return "ORI";
		if(phi==xlogx) return "OMI";
		return "OCI";
	}
	
	public static void main(String[] args){
		System.err.println();
	}
	
	
}
