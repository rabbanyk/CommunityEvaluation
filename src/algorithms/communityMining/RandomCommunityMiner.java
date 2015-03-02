package algorithms.communityMining;


import java.util.Map;
import java.util.Set;
import java.util.Vector;

import algorithms.communityMining.data.Grouping;
import edu.uci.ics.jung.graph.Graph;

public  class RandomCommunityMiner<V,E> extends CommunityMiner<V, E>{
	
	public  Grouping<V> findCommunities(Graph<V,E> graph){
			return findCommunities(graph, null);
	}
	public Grouping<V> findCommunities(Graph<V,E> graph,  Map<E, ? extends Number> weights){
		return findCommunities(graph, weights, null);
	}
	public Grouping<V> findCommunities(Graph<V,E> graph,  Map<E, ? extends Number> weights, Vector<Set<V>> groundTruth){
		
		Grouping<V> res = null;
		
		//TODO: Bring from relative experiment
		
		
		return res;
	}
	
	
	public  String getName(){return "random";};
	public  String getShortName(){return "random";};
	public String toString(){
		return getName();
	}
}
