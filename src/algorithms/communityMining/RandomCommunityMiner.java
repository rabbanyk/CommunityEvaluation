package algorithms.communityMining;


import java.util.Map;
import java.util.Set;
import java.util.Vector;

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
public  class RandomCommunityMiner<V,E> extends CommunityMiner<V, E>{
	
	public  Grouping<V> findCommunities(Graph<V,E> graph){
			return findCommunities(graph, null);
	}
	public Grouping<V> findCommunities(Graph<V,E> graph,  Map<E, ? extends Number> weights){
		return findCommunities(graph, weights, null);
	}
	public Grouping<V> findCommunities(Graph<V,E> graph,  Map<E, ? extends Number> weights, Vector<Set<V>> groundTruth){
		
		Grouping<V> res = null;
		
		
		
		
		return res;
	}
	
	
	public  String getName(){return "random";};
	public  String getShortName(){return "random";};
	public String toString(){
		return getName();
	}
}
