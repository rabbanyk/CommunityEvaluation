package measure.cluster.agreement.partitioning.classics;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Created on Apr 19, 2012
 * 
 * <p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
 *   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same data-set based on the Jaccard Index, 
 *   look up the definition in  <a href="http://en.wikipedia.org/wiki/Jaccard_index"> wikipedia </a>.<br/>
 *  
 *  <b>Jaccard</b> similarity coefficient measures similarity of two sets as the fraction of their intersection to their union. 
 *  If we consider co-membership of data points in the same or different clusters as a binary variable, we could Jaccard agreement between clustering <tt>U</tt> and <tt>V</tt> as follows:
 *  <br/>
 *  <img src='http://mathurl.com/crazzyg.png' alt='Rand' title='Rand'  /> <br/>
 *  </p>
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V> type of the data object
 */
public class Jaccard<V> extends PartiotioningAgreement<V>{

	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = 0,sumxy=0,sumx=0,sumy=0;
		
		n = n(groundTruth);
		sumxy = n2(partitioning,groundTruth);
		sumx = n2(partitioning);
		sumy = n2(groundTruth);
		
		double Jaccard = (n - sumxy ) / (n + sumxy - sumx - sumy);
		return Jaccard;
	}


	
	public String toString(){
		return "Jaccard";
	}
	
	/**
	 * A different implementation for {@link #getAgreement(Vector, Vector) getAgreement} method
	 * 
	 * @param partitioning
	 * @param groundTruth
	 * @return same value as {@link #getAgreement(Vector, Vector) getAgreement} method
	 */
	public double getAgreementAlt2(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = 0, nx, ny, nxy, sumxy=0,sumx=0,sumy=0;
		
		for (Set<V> x : groundTruth) {
			nx = x.size();
			for (Set<V> y : partitioning) {
				ny = y.size();

				Set<V> z = new HashSet<V>(x);
				z.retainAll(y);
				nxy = z.size();

				sumxy += nxy *nxy ;
			}
		}
		
		for (Set<V> x : groundTruth) {
			nx = x.size();
			n += nx;
			sumx += nx* nx;
		}
		
		for (Set<V> y : partitioning) {
			ny = y.size();
			sumy += ny*ny ;
		}
		
		double Jaccard = (n - sumxy ) / (n+ sumxy - sumx - sumy);
		return Jaccard;
	}
	
	
	
	private Set<Pair<V>> computeSameCommunitySet(Vector<Set<V>> part){		
		Set<Pair<V>> sameCommunitySet = new HashSet<Pair<V>>();
		
		for(Set<V> community : part){
			for(V vertex : community){
				for(V vertex2 : community){
					if (vertex == vertex2) 
						continue;					
					
					sameCommunitySet.add(new Pair<V>(vertex, vertex2));					
				}
			}
		}
		
//		Set<V> hubs = part.getHubs();
//		for(V vertex : hubs){
//			for(V vertex2 : hubs){
//				if (vertex == vertex2) 
//					continue;					
//					
//				sameCommunitySet.add(new Pair<V>(vertex, vertex2));					
//			}
//		}
//		
//		Set<V> outliers = part.getOutliers();
//		for(V vertex : outliers){
//			for(V vertex2 : outliers){
//				if (vertex == vertex2) 
//					continue;					
//					
//				sameCommunitySet.add(new Pair<V>(vertex, vertex2));					
//			}
//		}
		return sameCommunitySet;
	}

	/**
	 * Another different (Justin's) implementation for {@link #getAgreement(Vector, Vector) getAgreement} method
	 * 
	 * @param partitioning
	 * @param groundTruth
	 * @return same value as {@link #getAgreement(Vector, Vector) getAgreement} method
	 */
	public double getAgreementAlt(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		Set<Pair<V>> partSet = computeSameCommunitySet(partitioning);
		Set<Pair<V>> gtSet = computeSameCommunitySet(groundTruth);
		Set<Pair<V>> union = new HashSet<Pair<V>>(partSet);
		union.addAll(gtSet);
		
		Set<Pair<V>> intersection = new HashSet<Pair<V>>(partSet);
		intersection.retainAll(gtSet);
		
		return (double)intersection.size() / (double)union.size();	
	}
	
}
