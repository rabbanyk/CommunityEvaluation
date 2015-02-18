package measure.cluster.agreement.partitioning.classics;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

public class VI <V> extends PartiotioningAgreement<V>{
	//TODO: why returns more that 1? 
	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double  n = Math.max(n(partitioning), n(groundTruth));
		if(n(partitioning)!=n(groundTruth))
			System.err.println("Should have equal size ------------- Erooooor"+n(partitioning)+"  "+n(groundTruth));
		if ( Math.min(n(partitioning), n(groundTruth))==0) {
			System.err.println("n is zeroooooooo---- VI erroor");
			return 1;
		}
		//Transformed to be normalized and larger better (1-VI/log(n)), differs form original VI
//		System.err.println("n: " +n+"  "+H(partitioning,n)+" "+H(groundTruth,n)+" "+I(partitioning, groundTruth));
		return 1- (H(partitioning,n)+H(groundTruth,n)-2*I(partitioning, groundTruth) )/ Math.log(n);
	}
	
	
	public double getAgreementAlt(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		double n=0,nx,ny,nxy;
		double VI=0;//, Hx=0, Hy=0;
		
		for(Set<V> y : groundTruth ){
			ny = y.size();
			n+=ny;
		}
		
		for(Set<V> x : partitioning ){
			nx = x.size();

			for(Set<V> y : groundTruth ){
				ny = y.size();
				
				Set<V> z = new HashSet<V>(x);
				z.retainAll(y);
				
				nxy = z.size();
				
				VI += nxy==0? 0 : (nxy/n * Math.log(nx * ny / (nxy * nxy))); //Not sure what the base should be
			}
		}
//		for(Set<V> x : partitioning ){
//			nx = x.size();
//			denominator += nx * (nx==0? 0 : Math.log(nx / n )); //Not sure what the base should be
//		}
//		
//		for(Set<V> y : groundTruth ){
//			ny = y.size();
//			denominator += ny *(ny==0? 0 :  Math.log(ny / n )); //Not sure what the base should be
//		}
//		
		return VI / Math.log(n);
	}
	
	
	public double getAgreementAlt2(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		double n=0,nx,ny,nxy;
		double H1=0,H2=0;//, Hx=0, Hy=0;
		
		for(Set<V> y : groundTruth ){
			ny = y.size();
			n+=ny;
		}
		
		for(Set<V> x : partitioning ){
			nx = x.size();

			for(Set<V> y : groundTruth ){
				ny = y.size();
				
				Set<V> z = new HashSet<V>(x);
				z.retainAll(y);
				
				nxy = z.size();
				
				H1 += nxy==0? 0 : ((nxy/n) * Math.log( nx/nxy)); //Not sure what the base should be
				H2 += nxy==0? 0 : ((nxy/n) * Math.log( ny/nxy)); //Not sure what the base should be
			}
		}
//		for(Set<V> x : partitioning ){
//			nx = x.size();
//			denominator += nx * (nx==0? 0 : Math.log(nx / n )); //Not sure what the base should be
//		}
//		
//		for(Set<V> y : groundTruth ){
//			ny = y.size();
//			denominator += ny *(ny==0? 0 :  Math.log(ny / n )); //Not sure what the base should be
//		}
//		
//		System.err.println(H1/ Math.log(n) +"  "+ H2/ Math.log(n));
		return ((H1+H2) / Math.log(n) ) ;
	}
	
	public String toString(){
		return "VI";
	}
}