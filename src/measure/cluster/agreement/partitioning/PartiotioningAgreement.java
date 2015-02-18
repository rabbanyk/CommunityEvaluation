package measure.cluster.agreement.partitioning;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.ClusteringAgreement;

/**
 * Created on Apr 19, 2012
 
 *<p> This class is the base class for different <tt>Partitioning Agreement</tt> indices. <br/> 
*   It returns the agreement between two given <b><i>disjoint</i></b>  partitionings/groupings/clusterings.<br/>  
*   </p>
 * 
 * TODO: extension to overlapping clusters and fuzzy clusters e.g. "A fuzzy extension of the Rand index and other related indexes for clustering and classification assessment R.J.G.B. Campello" <br/> 
 * 
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V> type of the data object
 *  
 */

// TODO: also it should be checked that how should we handle clusterings that do not cover each other? say the given clustering dosen't include some data points in its results,
// Or if looking at the general picture that the second partitioning is not necessarily the ground-truth, there might be some data points that are in the second clustering but not in the first one.  
public abstract class PartiotioningAgreement<V> extends ClusteringAgreement<V> {
	
	
	protected double n(Vector<Set<V>> U){
		double n = 0;
		for (Set<V> x : U) 	n += x.size();
		return n;
	}
	
	protected double n(Vector<Set<V>> U,Vector<Set<V>> V){
		double  sumxy=0;

		for (Set<V> x : V) {
		for (Set<V> y : U) {
			Set<V> z = new HashSet<V>(x);
			z.retainAll(y);
			sumxy += z.size();
			}
		}
		return sumxy;
	}
	protected double n2(Vector<Set<V>> U){
		double  nx,sumx=0;
		for (Set<V> x : U) {
			nx = x.size();
			sumx += nx*nx;
		}
		return sumx;
	}
	protected double n2(Vector<Set<V>> U,	Vector<Set<V>> V){
		double nxy, sumxy=0;
		for (Set<V> x : V) {
		for (Set<V> y : U) {
			Set<V> z = new HashSet<V>(x);
			z.retainAll(y);
			nxy = z.size();
			sumxy += nxy *nxy ;
			}
		}
		return sumxy;
	}
	
	private double nlog(Vector<Set<V>> U,	Vector<Set<V>> V){
		double sumxy=0;
		for (Set<V> x : V) {
		for (Set<V> y : U) {
			Set<V> z = new HashSet<V>(x);
			z.retainAll(y);
			sumxy += z.size()==0 ? 0 : z.size() * Math.log(z.size()) ;
			}
		}
		return sumxy;
	}
		
	private double nlog(Vector<Set<V>> U){
		double sumx=0;
		for (Set<V> x : U) {
			sumx += x.size()==0 ? 0 : (x.size() * Math.log(x.size()));
		}
		return sumx;
	}
	
	/**
	 * @param U
	 * @return entropy of clustering <tt>U</tt> on a dataset consisting of n elements
	 */
	protected double H(Vector<Set<V>> U, double n){
		n=n(U);
		if(n==0) {
			System.err.println("----------------------in H n==0");
			return 0;
		}
		return Math.log(n)*n(U)/n - nlog(U)/n;
	}
	
	/**
	 * @param U
	 * @return entropy of clustering <tt>U</tt> on a dataset consisting of n elements
	 */
	protected double HAlt(Vector<Set<V>> U, double n){
		double res = 0;
		double nx=0;
		
		//probability that an data points belongs to cluster x is nx/n 
		for(Set<V> x : U ){
			nx = x.size();
			res += nx/n * (nx==0?0:Math.log(nx / n )); 
		}
		return res * -1;
	}
	
	/**
	 * @param U
	 * @param V
	 * @return joint entropy of clustering <tt>U</tt> and <tt>V</tt>
	 */
	//works only if n(V) = n(U) otherwise returns non-valid results like 1.4
	protected double H(Vector<Set<V>> U, Vector<Set<V>> V){
		double n = n(V);
		return Math.log(n)*n(U)/n + nlog(U,V)/n;
	}
	
	//works only if n(V) = n(U) otherwise returns non-valid results like 1.4
	protected double Ialt(Vector<Set<V>> U, Vector<Set<V>> V){
		double n = n(V);
		return  Math.log(n)+ nlog(U,V)/n - nlog(V)/n- nlog(U)/n;
	}
	
	/**
	 * @param U
	 * @param V
	 * @return mutual information of clustering U and V
	 */
	protected double I(Vector<Set<V>> U, Vector<Set<V>> V){
		double n=0,nx=0,ny=0,nxy=0,res = 0;
//		for(Set<V> y : V ) n += y.size();
		
		n = Math.max(n(U), n(V));
		
		if (n==0){
			System.err.println("----------------------in I n==0");
			return 0;
		}
		for(Set<V> x : U ){
			nx = x.size();
			
			for(Set<V> y : V ){
				ny = y.size();
				Set<V> z = new HashSet<V>(x);
				z.retainAll(y);
				nxy = z.size();
				res += (nxy/n) * (nxy==0?0:(Math.log((nxy * n) / (nx * ny)))); 
			}
		}
		return res;
	}
	
		
	
}
