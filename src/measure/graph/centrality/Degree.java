package measure.graph.centrality;

import java.util.Collection;
import java.util.HashMap;

public class Degree<V,E>  extends GraphVertexScorer<V,E>{

	
	public Double getVertexScore(V v){
		if(vertexScores == null)
			vertexScores = new HashMap<V, Double>();
		
		if(!vertexScores.containsKey(v)){
			Collection<E> incidentEdges = graph.getIncidentEdges(v);
			double deg = 0;
			if(weighted()){
				for (E e : incidentEdges) {
					deg += weights.transform(e).doubleValue(); //TODO: check this with wighted degree centrality definition also check if normalization is needed?
				}
			}else{
				deg = incidentEdges.size();//graph.degree(v)
			}
			vertexScores.put(v , deg);
		}
		
		return vertexScores.get(v);	
	}
	
	public String toString(){
		return "Degree";
	}

}
