package measure.cluster.agreement.partitioning.generalization;

import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.generalization.GAM.Type;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.util.Pair;

public class AGAMA<V> extends GAM<V>{

	
	public AGAMA(Type type) {
		super(type);
	}
	public AGAMA(Transformer<Double, Double> phi, Transformer<Pair<Set<V>>, Double> etha) {
		super(phi,etha);
	}

	public double getAgreement(Vector<Set<V>> U, Vector<Set<V>> V) {
		
		double uvp=0, upv=0, vpu=0, puv=0;
		double  tmp;
		for (Set<V> v : V) {
			double  su=0;
			for (Set<V> u : U) {
				tmp = etha.transform(new Pair<Set<V>>(v,u));
				
				su += tmp; 
				puv += tmp;
				uvp += phi.transform(tmp);
			}
			vpu += phi.transform(su);
		}
		
		for (Set<V> u : U) {
			double  sv=0;
			for (Set<V> v : V) {
				tmp = etha.transform(new Pair<Set<V>>(u,v)); 
				sv += tmp;
			}
			upv += phi.transform(sv);
		}
		
		puv =  phi.transform(puv);
		
		//double E = vpu*upv;
		//if(E!=0) E/=puv;
		double res = (uvp );//- E);
		//double dom = (upv + vpu)/2;
//		if(res!=0 ) if(dom==0) res = 0; else res/=( (upv + vpu)/2 );//- E );
		return res;
	}
	public String toString(){
		return super.toString();
	}
}
