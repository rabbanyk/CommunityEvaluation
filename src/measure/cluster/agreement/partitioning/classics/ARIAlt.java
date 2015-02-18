package measure.cluster.agreement.partitioning.classics;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on Apr 18, 2012
 *
 *<p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
*   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same dataset based on the Adjusted Rand Index presented in: <br/>  
*   Hubert, L. & Arabie, P. Comparing partitions Journal of Classification, 1985, 2, 193-218
*   </p>
*   
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V> type of the data object 
 */
public class ARIAlt<V> extends PartiotioningAgreement<V>{
	
	public double getAgreement(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		double ARI = 0;

		double n = 0, nx, ny, nxy, sumxy=0,sumx=0,sumy=0;
		
		for (Set<V> x : groundTruth) {
			nx = x.size();
			for (Set<V> y : partitioning) {
				ny = y.size();

				Set<V> z = new HashSet<V>(x);
				z.retainAll(y);
				nxy = z.size();

				sumxy += nxy * (nxy - 1) / 2;
			}
		}
		
		for (Set<V> x : groundTruth) {
			nx = x.size();
			n += nx;
			sumx += nx*(nx-1)/2;
		}
		
		for (Set<V> y : partitioning) {
			ny = y.size();
			sumy += ny*(ny-1)/2;
		}
		
		double mulxy = sumx * sumy / (n*(n-1)/2); 
			
		ARI = (sumxy - mulxy)/( (sumx+sumy)/2 - mulxy);
		return ARI;
	}
	
	
//	Alternative implementation:
	/**
	 * 
	 * An alternative implementation of the ARI by formula in Jorge M. Santos and Mark Embrechts. On the use of the adjusted rand index as a metric for evaluating supervised classification. In ICANN (2), pages 175–184, 2009.
	 *
	 * @param nodes
	 * @param partitioning
	 * @param groundTruth
	 * @return same value as the getAgreement, 
	 */
	public double getAgreement(Vector<V> nodes, Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		if(nodes == null){
			nodes = new Vector<V>();
			for (Set<V> set : groundTruth) {
				nodes.addAll(set);
			}
		}
		
		double ARI = 0;

		int a = 0;// # of pairs that are in same community in both clusters and
					// groundTruth
		int b = 0;// # of pairs that are in same community in groundTruth but
					// not in clusters
		int c = 0;// # of pairs that are in same community in clusters but not
					// in groundTruth
		int d = 0;// # of pairs that are in different community in both clusters
					// and groundTruth

		Set<V> cc1, cg1, cc2, cg2;

		double np = 0;

		for (int i =0; i<  nodes.size(); i++) {
			V v1 = nodes.get(i);
			cc1 = getCommunity(v1, partitioning);

			cg1 = getCommunity(v1, groundTruth);

			for (int j =i+1; j<  nodes.size(); j++) {
				V v2 = nodes.get(j);

				if (!v1.equals(v2)) {
					np++;

					cc2 = getCommunity(v2, partitioning); 
					cg2 = getCommunity(v2, groundTruth);

					if(cc1 == null || cc2 == null) {
						if(cg2.equals(cg1)) b++; else d++;
						continue;
					}
					
					if (cc2.equals(cc1) && (cg2.equals(cg1)))
						a++;
					else if (!cc2.equals(cc1) && cg2.equals(cg1))
						b++;
					else if (cc2.equals(cc1) && !cg2.equals(cg1))
						c++;
					else if (!cc1.equals(cc2) && !cg1.equals(cg2))
						d++;
				}
			}
		}
		double tmp = (a + b) * (a + c) + (c + d) * (b + d);
		ARI = np * (a + d) - tmp;
		ARI /= np * np - tmp;

		return ARI;
	}
	
	
	private Set<V> getCommunity(V e, Vector<Set<V>> communities){
		for (Set<V> cluster : communities) {
			if(cluster.contains(e))
				return cluster;
		}
//		if (hubs.contains(e)) return hubs;
//		return outliers;
		return null;
	}
	
	public String toString(){
		return "ARIAlt";
	}

}
