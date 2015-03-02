package algorithms.communityMining.topleaders;

import java.util.ArrayList;
import java.util.Map;

import org.apache.commons.collections15.Transformer;

import static io.Logger.*;
import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.dev_.Partitioning;
import algorithms.communityMining.topleaders.global.GlobalTopLeaders;
import edu.uci.ics.jung.graph.Graph;


/**
 * @author Reihaneh
 *
 * @param <V>
 * @param <E>
 */
public class TopLeaders<V,E> extends CommunityMiner<V, E> {
 
	double outlierThereshod = 0, centersClosenessThreshold = .1, 
			minCommunitySizeThreshold = 4, hubThreshold = 0;
	
	public TopLeaders() {
		super();
	}
	public Grouping<V> findCommunities(int k, Graph<V, E> graph, Map<E, ? extends Number> weights) {
		Transformer<Graph<V, E>, Partitioning<V>>  topLeaders=  
				new GlobalTopLeaders<V, E>(k, outlierThereshod, hubThreshold, centersClosenessThreshold);
		//TODO: add weights
		return new Grouping<V>(topLeaders.transform(graph).getCommunities());
	
	}
	@Override
	public Grouping<V> findCommunities(Graph<V, E> graph, Map<E, ? extends Number> weights) {
		return findCommunities(11,graph,weights);
	}
	
//	public Grouping<V> findCommunities(Graph<V,E> graph, int k, ArrayList<V> initialCenters, int distanceMeasureType){
//			double p[] = guessLocalParamsBasedOnGraph(graph);
//			//double outlierThereshod, double hubThreshold, double minCommunitySizeThreshold,double centersClosenessThreshold
//			topLeaders = new GlobalTopLeaders<V, E>(k, p[0],p[1], p[3]);
//			((GlobalTopLeaders)topLeaders).setHardCodedInitialCenters(initialCenters);
//			((GlobalTopLeaders)topLeaders).setDistanceMeasureType(distanceMeasureType);
//			return new Grouping(topLeaders.transform(graph).getCommunities());
//		}

	
	public String toString(){
		return "TopLeaders";
	}

	public String getName() {
		return "TopLeaders";
	}

	public String getShortName() {
		return "TL";
	}

	
}
