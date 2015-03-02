package algorithms.dev_topleaders;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import org.apache.commons.collections15.Predicate;
import edu.uci.ics.jung.algorithms.filters.VertexPredicateFilter;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;


public abstract class CommunityMiner<V,E> {
	
	HashMap<Integer, V> idVertHash;
	
	public abstract Partitioning<V> findCommunities(Graph<V,E> graph);
	
	//public abstract Partitioning<V> findCommunities(Transformer<String, V> vertexTransformer, String network, String path);
	
	public String toString(){
		return getClass().getSimpleName();
	}
	
	
	public String getName(){
		return this.toString();
	};
	
	public String getShortName(){
		return getName();
	}

}
