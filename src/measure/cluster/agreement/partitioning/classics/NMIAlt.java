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
public class NMIAlt <V> extends PartiotioningAgreement<V>{
	
	
	public double getAgreement(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		double n=0,nx,ny,nxy;
		double I=0,denominator=0;
		
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
				
				I += nxy==0? 0 : nxy * Math.log(nxy * n / (nx * ny)); //Not sure what the base should be
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
		
		return -2*I/denominator;
	}

	
	public String toString(){
		return "NMIAlt";
	}
	

}
