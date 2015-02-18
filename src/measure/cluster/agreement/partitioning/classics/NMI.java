package measure.cluster.agreement.partitioning.classics;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on Apr 19, 2012
 *<p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
*   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same dataset based on the Normalised Mutual Information presented in: <br/>  
 * Kvalseth, T. O. Entropy and Correlation: Some Comments Systems, Man and Cybernetics, IEEE Transactions on, 1987, 17, 517 -519 <br/>
 *    
 * <img src='http://mathurl.com/dx7yazl.png' alt='formula' title='formula'  /> <br/>
 * </p>
*   
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V> type of the data object 
 */
public class NMI <V> extends PartiotioningAgreement<V>{
	//TODO: why returns more that 1? 
	public enum Mode{SUM, MIN, MAX, SQRT};
	Mode mode = Mode.SUM;
	public NMI() {
		this(Mode.SUM);
	}
	public NMI(Mode mode) {
		super();
		this.mode = mode;
	}
	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		//TODO: this formulation would return zero if both clustering put all the datapoints in one cluster
				//sklearn has a if, checks if it is the case and returns zero, 
				//we did that too. However the formula (both versions) gives zero. Shall we keep it?
		if(partitioning.size()==1 && partitioning.equals(groundTruth))	return 1;

		switch (mode) {
		case SUM:
			return getAgreementSum(partitioning, groundTruth);
		case SQRT:
			return getAgreementSqrt(partitioning, groundTruth);
		case MIN:
			return getAgreementMin(partitioning, groundTruth);
		case MAX:
			return getAgreementMax(partitioning, groundTruth);
		default:
			return getAgreementSum(partitioning, groundTruth);
		}
	}
	public double getAgreementSum(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = n(groundTruth);
		double res = 2*I(partitioning, groundTruth);
		if(res!=0) res/= (H(partitioning,n)+H(groundTruth,n));
		return res;
//		return 2*I(partitioning, groundTruth) / (H(partitioning,n)+H(groundTruth,n));
	}
	
	public double getAgreementMin(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = n(groundTruth);
		double res = I(partitioning, groundTruth);
		if(res!=0) res/= Math.min(H(partitioning,n),H(groundTruth,n));
		return res;
//		return I(partitioning, groundTruth) / Math.min(H(partitioning,n),H(groundTruth,n));
	}
	
	public double getAgreementMax(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = n(groundTruth);
		double res = I(partitioning, groundTruth);
		if(res!=0) res/=  Math.max(H(partitioning,n),H(groundTruth,n));
		return res;
//		return I(partitioning, groundTruth) /  Math.max(H(partitioning,n),H(groundTruth,n));
	}
	
	public double getAgreementSqrt(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = n(groundTruth);
		double res = I(partitioning, groundTruth);
		if(res!=0) res/= Math.sqrt(H(partitioning,n)*H(groundTruth,n));
		return res;
	}
	
	public double getAgreementJoint(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		return I(partitioning, groundTruth) / H(partitioning,groundTruth);
	}
	
	
	public double getAgreementStandalone(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		double n=0,nx,ny,nxy;
		double numerator=0,denominator=0;
		
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
				
				numerator += nxy==0? 0 : nxy * Math.log(nxy * n / (nx * ny)); //Not sure what the base should be
			}
		}
		
		for(Set<V> x : partitioning ){
			nx = x.size();
			denominator += nx * (nx==0? 0 : Math.log(nx / n )); //Not sure what the base should be
		}
		for(Set<V> y : groundTruth ){
			ny = y.size();
			denominator += ny *(ny==0? 0 :  Math.log(ny / n )); //Not sure what the base should be
		}
		
		return numerator==0? numerator: (-2*numerator/denominator);
	}

	
	public String toString(){
//		return "NMI";
		return "NMI"+(mode==Mode.MIN?("min"):(mode==Mode.MAX?("max"):(mode==Mode.SQRT?("\u221A"):("\u208A"))));

	}
	public String toLatexString(){
		return "NMI"+(mode==Mode.MIN?("_{min}"):(mode==Mode.MAX?("_{max}"):(mode==Mode.SQRT?("_\\sqrt{.}"):("_\\Sigma"))));

	}

}
