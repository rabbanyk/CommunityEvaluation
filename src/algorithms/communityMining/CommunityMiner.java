package algorithms.communityMining;


import java.util.Map;

import algorithms.communityMining.data.Grouping;
import edu.uci.ics.jung.graph.Graph;

/*
 * List: http://complexinfo.wordpress.com/
 * WebService: http://hub.iis.sinica.edu.tw/spotlight/Help/main.htm
 * GraphSteream in Java: http://www.graphstream-project.org/download/
 * iGraph in R: http://igraph.sourceforge.net/download.html
 * LP in iGraph: http://rgm2.lab.nig.ac.jp/RGM2/func.php?rd_id=igraph:label.propagation.community
 * http://rgm2.lab.nig.ac.jp/RGM2/func.php?rd_id=igraph:label.propagation.community
 * 

*/
public abstract class CommunityMiner<V,E> {
	
	public  Grouping<V> findCommunities(Graph<V,E> graph){
			return findCommunities(graph, null);
	}
	public abstract Grouping<V> findCommunities(Graph<V,E> graph,  Map<E, ? extends Number> weights);
	
	
	public abstract String getName();
	public abstract String getShortName();
	public String toString(){
		return getName();
	}
}
