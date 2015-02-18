package measure.cluster.agreement;

public class NMI_Overlapping_MacDaid<V> extends ClusteringAgreementExecWrapper<V> {
	{
		exePath = "./execs/AM-Overlapping-NMI-McDaid/onmi";
		infoMarker = "NMI<Max>:	";
	}
	public String toString(){
		return "NMI''";
	}
	// execs/LFR/MultiplexNMI-Benchmark/MultiplexNMI/nmi$

		// NMI<Max>: 0.53781
		// Other measures:
		// lfkNMI: 0.610344
		// NMI<Sum>: 0.617992
		// reihaneh@Rey:~/projects/AttributeGraphJava/execs/AM-Overlapping-NMI-McDaid$ ./onmi tv.txt u2.txt
		// Warning: two consecutive tabs, or tab at the start of a line. Ignoring empty fields like this
		// NMI<Max>: 0.615842
		// Other measures:
		// lfkNMI: 0.628119
		// NMI<Sum>: 0.615842

}
