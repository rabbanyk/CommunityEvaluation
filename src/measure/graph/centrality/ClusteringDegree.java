package measure.graph.centrality;

import java.util.Collection;
import java.util.HashMap;

/**
 * @author rabbanyk
 *
 * A Centrality measure that combines degree centrality and local clustering coefficient.
 * This helps lowering the centrality of hubs as they have low clustering coefficient -- their neighbors are not tied together.
 * It is computed by counting the number of triangles attached to a node, which is same as 2x(degree)x(degree-1)x(clustering_coefficient)  
 * weighted version is defined arbitrary, could be changed. 
 * One possible alternative is to implement it to match Clustering in weighted networks Tore Opsahl, however his definition ignores the weight on the 3rd edge
 */
public class ClusteringDegree<V,E>  extends GraphVertexScorer<V,E>{

	
	public Double getVertexScore(V v){
		if(vertexScores == null)
			vertexScores = new HashMap<V, Double>();
		
		if(!vertexScores.containsKey(v)){
			Collection<E> incidentEdges = graph.getIncidentEdges(v);
			double score = 0;
			
			if(weighted()){ 
				//Check weighted with Clustering in weighted networks Tore Opsahl, 
				for(E e1: incidentEdges){
					for(E e2: incidentEdges){
						if(e1!=e2){
							E e = graph.findEdge(graph.getOpposite(v, e1), graph.getOpposite(v, e2));
							if(e!=null) {
								score += Math.min(weights.transform(e).doubleValue(), weights.transform(e1).doubleValue()); 
								//TODO: could be anything: arithmetic mean, geometry mean, min, max, product ...
							}
						}
					}
				}
			}
			else{
				for(V n1: graph.getNeighbors(v)){
					for(V n2: graph.getNeighbors(v)){
						if(n1!=n2 && graph.findEdge(n1, n2)!=null) score++;
					}
				}
			}
			
			vertexScores.put(v , score);
		}
		
		return vertexScores.get(v);	
	}

	public String toString(){
		return "ClusteringDegree";
	}
}
