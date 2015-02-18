package measure.cluster.agreement;

public class NMI_Overlapping_LFR<V> extends ClusteringAgreementExecWrapper<V> {
	{
		exePath = "./execs/AM-Overlapping-NMI-Lancich/mutual";
		infoMarker = "mutual3:	";
	}

	public String toString(){
		return "NMI'";
	}
}
