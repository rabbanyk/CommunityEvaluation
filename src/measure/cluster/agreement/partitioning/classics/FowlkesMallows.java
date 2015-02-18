package measure.cluster.agreement.partitioning.classics;

import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

/**
 * Created on May 7, 2012
 * <p> This class is an extension of the abstract class <tt>Partitioning Agreement</tt>. <br/> 
 *   It returns the agreement between two given <b><i>disjoint</i></b> partitionings/groupings/clusterings of a same data-set based on the Fowlkes and Mallows Index.
 *   <b>Fowlkes and Mallows</b> is geometric mean of precision and recall defined in F-measure which are also known as the two nonsymmetric <b>Wallace</b> measures of partition correspondence. 
 *   See [Fowlkes and Mallows(1983)] Fowlkes EB, Mallows CL (1983) A method for comparing two hierarchical clusterings. Journal of the American Statistical Association 78(383):pp. 553–569
 *  <br/>
 *  <img src='http://mathurl.com/cl9lc6l.png' alt='formula' title='formula'  /> <br/>
 *  </p>
 *  
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 * @param <V>
 */
public class FowlkesMallows<V> extends PartiotioningAgreement<V>{

	public double getAgreement(Vector<Set<V>> partitioning,	Vector<Set<V>> groundTruth) {
		double n = 0,sumxy=0,sumx=0,sumy=0;
		
		n = n(groundTruth);
		sumxy = n2(partitioning,groundTruth)-n;
		sumx = n2(partitioning)-n;
		sumy = n2(groundTruth)-n;
		
		double FM = sumxy / Math.sqrt(sumx * sumy);
		return FM;
	}

	
	public String toString(){
		return "FowlkesMallows";
	}
	
}
