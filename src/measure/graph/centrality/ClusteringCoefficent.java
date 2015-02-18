package measure.graph.centrality;

import java.util.Collection;
import java.util.HashMap;


// local clustering coefficient = close triplets/total triplets
public class ClusteringCoefficent<V, E>  extends GraphVertexScorer<V,E>{
//TODO: do we need a weighted an unweighted version? weighted is same as unweighted
	
	public Double getVertexScore(V v){
		if(vertexScores == null)
			vertexScores = new HashMap<V, Double>();
		
		if(!vertexScores.containsKey(v)){
			Collection<E> incidentEdges = graph.getIncidentEdges(v);
			double score = 0;
		
			if(weighted()){ //Check weighted with Clustering in weighted networks Tore Opsahl, 
				double dom = 0;
				for(E e1: graph.getIncidentEdges(v)){
					for(E e2: graph.getIncidentEdges(v)){
						if(e1!=e2){
							dom +=  Math.min(weights.transform(e1).doubleValue(),weights.transform(e2).doubleValue());
							E e = graph.findEdge(graph.getOpposite(v, e1), graph.getOpposite(v, e2));
							if(e!=null){
								score += Math.min(weights.transform(e).doubleValue(),weights.transform(e1).doubleValue()); 
								//TODO: could be anything: arithmetic mean, geometry mean, min, max, product ..., Tore only considered e1 and e2
							}
						}
					}
				}
				score /= dom;
			}else{
				for(V n1: graph.getNeighbors(v)){
					for(V n2: graph.getNeighbors(v)){
						if(n1!=n2 && graph.findEdge(n1, n2)!=null) 
							score++;
					}
				}
				
				double deg = incidentEdges.size();
				score /= deg * (deg-1) ; 
			}
			
			vertexScores.put(v , score);
		}
		
		return vertexScores.get(v);	
	}
	
	public String toString(){
		return "ClusteringCoefficent";
	}
}
