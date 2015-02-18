package measure.cluster.agreement.partitioning.classics;

import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on May 7, 2012
*<p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
 *   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same data-set based on the Rand Index.
 *   <b>Rand Index</b> is similar to Jaccard index but also prizes true negatives. 
 *   
 *  <br/>
 *  <img src='http://mathurl.com/75w8567.png' alt='formula' title='formula'  /> <br/>
 *  </p>
 *  
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V>
 */
public class Rand<V> extends PartiotioningAgreement<V>{

	
	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = 0,sumxy=0,sumx=0,sumy=0;
		
		n = n(groundTruth);
		sumxy = n2(partitioning,groundTruth)-n(partitioning,groundTruth);
		sumx = n2(partitioning)-n(partitioning);
		sumy = n2(groundTruth)-n(groundTruth);
		
		double Rand = 1+ 1/(n*(n-1)) * (2* sumxy -  (sumx  + sumy));
		return Rand;
	}

	public double getAgreementStandalone(Vector<V> nodes, Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		if(nodes == null){
			nodes = new Vector<V>();
			for (Set<V> set : groundTruth) {
				nodes.addAll(set);
			}
		}
		
		double RI = 0;

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
		RI = np * (a + d) - tmp;
		RI /= np * np - tmp;

		return RI;
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
		return "Rand";
	}
	

}
