package measure.graph.distance;

import measure.graph.GraphProximity;
import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;
import edu.uci.ics.jung.graph.Graph;

public class ShortestPath<V,E> extends GraphProximity<V,E> {
		{
			type = Type.DISTANCE;
		}
		Distance<V> distance;
		
		public ShortestPath(Graph<V, E> g) {
			super(g);
			distance =new UnweightedShortestPath<V,E>(graph);
		}
		public ShortestPath(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
			super(g, weights);
			if(weights!=null) distance = new DijkstraDistance<V,E>(graph, weights);
			else distance =new UnweightedShortestPath<V,E>(graph);

		}

		public Number computeMeasure(V source, V target) {
			return distance.getDistance(source, target);
		}
		
		public Number getSimilarity(V source, V target){
			return graph.getVertexCount()*MaxWeight  - getProximity(source, target).doubleValue();
		}
		
		public String toString(){
			return "ShortestPath"+type;
		}
		@Override
		public String getName() {
			return "SP";
		}

}
