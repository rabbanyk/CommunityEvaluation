package measure.cluster.agreement.partitioning.classics;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on Apr 18, 2012
 *
 * <p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
 *   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same data-set based on the Adjusted Normalized Mutual Information presented in: <br/>  
 *   Vinh, N. X.; Epps, J. & Bailey, J. Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance Journal of Machine Learning Research, 2010, 11, 2837-2854 <br/>
 *   More specifically, this is implementation of AMI_{sum} variation mentioned in the paper.
 *   </p>
 *   
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a> 
 *
 * @param <V> type of the data object 
 */
public class AMI <V> extends PartiotioningAgreement<V>{
	
	/**
	 * 
	 * This function calculates multiplication of all the numbers in <tt>num</tt> divided by all the numbers in <tt>den</tt>. <br/> 
	 * It doses this calculation by first crossing out all the same numbers in the numerator and the denominator, which is very useful if we are computing fractions of several factorials.  
	 * 
	 * @param num
	 * @param den
	 * @return    multiplication of all the numbers in <tt>num</tt> divided by all the numbers in <tt>den</tt>. <br/>
	 */ //this returns more than one sometimes
	private double fac(Vector<Integer> num, Vector<Integer> den){
		double res = 1;

		Vector<Integer> crossout = new Vector<Integer>();
		for (Integer integer : den) {
			if(num.contains(integer)){
				num.remove(integer);
				crossout.add(integer);
			}
		}
		
		for (Integer integer : crossout) {
			den.remove(integer);
		}
		
		for (Integer integer : num) {
			res*=integer;
		}
		for (Integer integer : den) {
			res/=integer;
		}
		
		return res;
	}
	
	private Vector<Integer> range(double n){
		Vector<Integer> res = new Vector<Integer>();
		for (int k = 1; k <= n; k++) {
			res.add(k);
		}
		return res;
	}
	
	/**
	 * 
	 * Refer to the paper for detailed description
	 * 
	 * @param U
	 * @param V
	 * @return expected mutual information between two clusterings <tt>U</tt> and <tt>V</tt>
	 */
	private double E(Vector<Set<V>> U, Vector<Set<V>> V){
		double n=0,nx=0,ny=0,nxy=0,res = 0;

		for(Set<V> y : V ){
			ny = y.size();
			n+=ny;
		}
		
		for(Set<V> x : U ){
			nx = x.size();
			for(Set<V> y : V ){
				ny = y.size();
				
				for(nxy = Math.max(1, nx+ny-n); nxy < Math.min(nx, ny); nxy++){ //TODO: shouldn't be less equal?
					double d;
					
					d = nxy/n;
					d *= Math.log((nxy * n )/ (nx * ny)); 
					
					Vector<Integer> nom = new Vector<Integer>();
					nom.addAll(range(nx));
					nom.addAll(range(ny));
					nom.addAll(range(n-nx));
					nom.addAll(range(n-ny));
					
					Vector<Integer> dom = new Vector<Integer>();
					dom.addAll(range(n));
					dom.addAll(range(nxy));
					dom.addAll(range(nx-nxy));
					dom.addAll(range(ny-nxy));
					dom.addAll(range(n-nx-ny+nxy));
					
					d*= fac(nom,dom);
					
					res+= d;
				}
			}
		}
		return res;
	}

	@Override
	public double getAgreement(Vector<Set<V>> partitioing, Vector<Set<V>> groundtruth) {
		double i= I(partitioing,groundtruth),
		n = n(groundtruth),
		max = ( H(groundtruth,n) + H(partitioing,n))/2 , 
		exp = E(groundtruth, partitioing);
		
		return (i - exp )/(max - exp);  
	}

	
	public String toString(){
		return "AMI";
	}
	
//	
//	public static void main(String[] params){
//		AMI<Integer> ami = new AMI<Integer>();
//
//		double n=34,nx=10,ny=5,nxy=5;//,res = 0;
//
//		Vector<Integer> nom = new Vector<Integer>();
//		nom.addAll(ami.range(nx));
//		nom.addAll(ami.range(n-nx));
//		nom.addAll(ami.range(ny));
//		nom.addAll(ami.range(n-ny));
//		
//		Vector<Integer> dom = new Vector<Integer>();
//		dom.addAll(ami.range(n));
//		dom.addAll(ami.range(nxy));
//		dom.addAll(ami.range(nx-nxy));
//		dom.addAll(ami.range(ny-nxy));
//		dom.addAll(ami.range(n-nx-ny+nxy));
//		
//		System.err.println(ami.fac(nom,dom));
//		
//		
//	}

}
