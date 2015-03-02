package measure.cluster.agreement;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

public abstract class ClusteringAgreement<d> {
	public  abstract double getAgreement(Vector<Set<d>> U, Vector<Set<d>> V);

	public  double getAgreementCovering(Collection<d> datapoints, Vector<Set<d>> partitioning, Vector<Set<d>> groundTruth){
		return getAgreement(coveringClustering(datapoints,partitioning), coveringClustering(datapoints, groundTruth));
	}
	
	public <d> Vector<Set<d>> coveringClustering(Collection<d> datapoints, Vector<Set<d>> C ){
		Set<d> unAssignedNodes = new HashSet<d>(datapoints);
		for (Set<d> c: C){
			c.retainAll(datapoints);
			unAssignedNodes.removeAll(c);
		}
		
		if (unAssignedNodes.size()>0)
			C.add(new HashSet<>(unAssignedNodes));
		
//		for (Set<d> c: C){
//			for (d datapoint : c)
//				if (! datapoints.contains(datapoint))
//					c.remove(datapoint);
//		}
		return C;
	}
	
	
	public String toString(){
		return this.getClass().getCanonicalName();
	}
	public String toLatexString(){
		return this.toString();
	}
}
