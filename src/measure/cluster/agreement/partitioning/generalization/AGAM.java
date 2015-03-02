package measure.cluster.agreement.partitioning.generalization;

import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.generalization.GAM.Type;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.util.Pair;

public class AGAM<V> extends GAM<V>{

	
	public AGAM(Type type) {
		super(type);
	}
	public AGAM(Transformer<Double, Double> phi, Transformer<Pair<Set<V>>, Double> etha) {
		super(phi,etha);
//		System.err.println(phi);
	}
	
	public double getAgreement(Vector<Set<V>> U, Vector<Set<V>> V) {
		double  tmp, I=0, E=0,n=0, sU=0, sV=0;
		double[] margV = new double[V.size()];
		double[] margU = new double[U.size()];
		
		for (int i = 0; i < V.size(); i++) {
			for (int j = 0; j < U.size(); j++) {
				tmp = etha.transform(new Pair<Set<V>>(V.elementAt(i),U.elementAt(j)));
				margV[i]+=tmp;
				margU[j]+=tmp;
				n += tmp;
				I += phi.transform(tmp);
			}
		}
		
		for (int i = 0; i < V.size(); i++) {
			for (int j = 0; j < U.size(); j++) {
				E+= phi.transform((margV[i]*margU[j])/n);
			}
		}
		for (int i = 0; i < V.size(); i++) {
			sV+= phi.transform(margV[i]);
		}
		for (int j = 0; j < U.size(); j++) {
			sU+= phi.transform(margU[j]);
		}
		
		double res = sU+sV-2*I;
		double dom = sU+sV-2*E;
		
//		if(this.type == Type.VI)
//			System.err.println( "------------------- sU:"+sU+" sV: "+sV+" I: "+I + " E: " +E);
//		if(res>Epsilon || res<-1*Epsilon){
//			if(dom <Epsilon && dom >-1*Epsilon){
//				System.err.println("--- sU+sV-2*I != 0 but sU+sV-2*E==0  :::: sU+sV-2*I="+ res+" / dom ="+ 
//							dom+ " sU:"+sU+" sV: "+sV+" I: "+I + " E: " +E);
//			}else 
//				res/=dom;
//		}else{
//			res =0;
//		}
		
		res = 1-res/dom;
		if(res> 1+Epsilon|| res<-1-Epsilon)
			System.err.println("--->> AGAM_"+type+" out of range!, error " + res );
		
		return res;
	}
	public double getAgreementOld(Vector<Set<V>> U, Vector<Set<V>> V) {
		
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
		
		double E = vpu*upv;
		
//		if(E==0)
//			System.err.println("--->> AGAM_"+type+" Expectaion is zero" );	
//		else if(puv==0)
//			System.err.println("--->> AGAM_"+type+" PUV shall not be zero, error" );
		
		
		if(E!=0){//{>Epsilon){
			if(puv<Epsilon){
				System.err.println("E!=0 but puv=0  :: E:"+ E
						+" puv:"+puv+" upv: "+upv+" vpu: "+vpu + " uvp: " +uvp);
			}else
				E/=puv;
		}
		double res = (uvp - E);
		
		double dom = (upv + vpu)/2.0 - E;

		if(res>Epsilon || res<-1*Epsilon){
			if(dom <Epsilon && dom >-1*Epsilon){
				System.err.println("--- (uvp-E) != 0 but 1/2(upv+vpu)-E==0  :::: (uvp-E)="+ res+" / dom ="+ 
							dom+ "=((upv + vpu)/2.0 :"+(upv + vpu)/2.0 +" -E:"+ E
						+") puv:"+puv+" upv: "+upv+" vpu: "+vpu + " uvp: " +uvp);
			}else 
				res/=dom;
		}else{
			res =0;
		}
	
		if(res> 1+Epsilon|| res<-1-Epsilon)
			System.err.println("--->> AGAM_"+type+" out of range!, error " + res );
		
		return res;
	}
	public String toString(){
		return "A"+super.toString();
	}
}
