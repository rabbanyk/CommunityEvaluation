package measure.cluster.agreement.partitioning.classics;

import java.awt.font.TextAttribute;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on Apr 18, 2012
 *
 <p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
 *   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same dataset based on the Fmeasure presented in: <br/>  
 *   <a href="http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html"> Christopher D. Manning, Prabhakar Raghavan and Hinrich Schütze, Introduction to Information Retrieval, Cambridge University Press. 2008.</a> <br/>
 *   <br/>
 *  <img src='http://mathurl.com/86nrgmw.png' alt='formula' title='formula'  /> <br/>
 *  </p>
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V> type of the data object 
 */

public class FMeasure<V>  extends PartiotioningAgreement<V>{

	double betha ;//= .5; // how much to penalize false negatives more strongly than false positives by setting this value higher than 1, thus giving more weight to recall.
	private double b ;//= betha*betha;
	
	public FMeasure() {
		this(1);
	}
	public FMeasure(double betha) {
		super();
		this.betha = betha;
		this.b = betha*betha;
	}

	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double  sumxy=0, sumx=0, sumy=0;
		
//		n = n(groundTruth);
//		sumxy = n2(partitioning,groundTruth) - n; // n^2 -> n(n-1)
//		sumx = n2(partitioning) - n;
//		sumy = n2(groundTruth) - n;
	
		sumxy = n2(partitioning,groundTruth)-n(partitioning,groundTruth);
		sumx = n2(partitioning)-n(partitioning);
		sumy = n2(groundTruth)-n(groundTruth);
		
//		if(sumx == 0||sumy==0 ||sumxy==0) {
//			System.err.println("-------------- "+sumx+" "+sumy+" "+sumxy+" ");
//			
//			System.err.println("-------------- "+partitioning);
//			System.err.println("-------------- "+groundTruth);
//		}
	
		double p = sumxy==0? 0: (sumxy / sumx);//precision
		double r = sumxy==0? 0: (sumxy / sumy);//recall
		
		double res = (( b + 1 ) *p*r);
		res = (res==0)? 0: (res/(b * p + r));
		
		return res;//(( b + 1 ) *p*r)/(b * p + r);
	}

	public double getAgreementAlt(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		Vector<V>nodes = new Vector<V>();
		for (Set<V> set : groundTruth) {
			nodes.addAll(set);
		}
		int tp = 0;// # of pairs that are in same community in both clusters and  groundTruth
		int fn = 0;// # of pairs that are in same community in groundTruth but not in clusters
		int fp = 0;// # of pairs that are in same community in clusters but not in groundTruth
		int tn = 0;// # of pairs that are in different community in both clusters and groundTruth
		Set<V> cc1, cg1, cc2, cg2;

		for (int i =0; i<  nodes.size(); i++) {
			V v1 = nodes.get(i);
			cc1 = getCommunity(v1, partitioning);
			cg1 = getCommunity(v1, groundTruth);

			for (int j =i+1; j<  nodes.size(); j++) {
				V v2 = nodes.get(j);
				cc2 = getCommunity(v2, partitioning); 
				cg2 = getCommunity(v2, groundTruth);

				if(cc1 == null || cc2 == null) {
					if(cg2.equals(cg1)) fn++; else tn++;
					continue;
				}
				
				if (cc2.equals(cc1) && (cg2.equals(cg1)))
					tp++; // TruePositive
				else if (!cc2.equals(cc1) && cg2.equals(cg1))
					fn++; // FalseNegative
				else if (cc2.equals(cc1) && !cg2.equals(cg1))
					fp++; // FalsePositive
				else if (!cc1.equals(cc2) && !cg1.equals(cg2))
					tn++; // TrueNegative
			}
		}
		
		double p = tp, r = tp;//precision and recall
		p/=(tp+fp);
		r/=(tp+fn); 
		
		return (( b + 1 ) *p*r)/(b * p + r);
		}
		
		
		private Set<V> getCommunity(V e, Vector<Set<V>> communities){
			for (Set<V> cluster : communities) {
				if(cluster.contains(e))
					return cluster;
			}
//			if (hubs.contains(e)) return hubs;
//			return outliers;
			return null;
		}
		//http://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts
		public String toLatexString(){
			return "F_{\\beta="+betha+"}";
		}
	
		public String toString(){
			return ("F\u1D66"+"\u208C"+(betha==2?"\u2082":(betha==.5?("\u2080"+"."+"\u2085"):betha)));//.addAttribute(TextAttribute.SUPERSCRIPT_SUB, TextAttribute.SUPERSCRIPT_SUB, 13, 14);
		}
	
		
	
}
