package io.graph.dev_plosOne;

import io.graph.pajek.PajekGraphReader;


public class dev_PlosGraphReader<V,E> extends PajekGraphReader<V, E> {

	{
		commentIndicator = "#";
		tokenizationPattern = "[;\\s]+";
		edgeListStartIndicator = "Edges";
		labeled = false;
	}
	/*
	 * # Vertices
0;-6.921557263198199;-0.17932736736429467;2
1;-1.955129009290761;-0.39765585887995547;0
#
# Edges
0;89
0;24
0;18
0;90
	 * */
}
