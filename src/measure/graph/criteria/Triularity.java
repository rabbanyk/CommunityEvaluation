package measure.graph.criteria;


import java.util.Vector;

import measure.graph.RelativeCommunityCriteria;
import edu.uci.ics.jung.graph.Graph;

/**
 * <pre> Q(undirected)= 1/2m Sum(i,j){Aij - (K(i)-K(j))/2m } Delta(i,j) </pre> 
 * where m is number of links, K(i) is degree of node i: # of links for node i  
 * <p> Delta(i,j) is one if i and j are in the same cluster, otherwise 0
 * <p> Aij is 1 for un-weighted networks: can be any number for weighted networks
 */

public class Triularity<V, E> extends RelativeCommunityCriteria<V,E> {

	public Triularity(Graph<V, E> graph) {
		super(graph);
	}

	//TODO: make it weighted !!!! edge count = sum W, size = sum W, 1 = W 
	public double evaluateCommunities(Vector<Graph<V,E>> communities){
		
		double  triularity = 0;
		double  MaxT = 0, ET = 0 , T =0; 

		double m = graph.getEdgeCount();//, n = graph.getVertexCount();
		MaxT = (2*m)*(2*m);//+(2*m);
		
		for (Graph<V,E> cluster : communities) {
		
			for (V v1 : cluster.getVertices()) {
				for (V v2 : cluster.getVertices()) {
					
					double di=0,dj=0,dk=0;//, dij=0, eij;
					
					for (V k : graph.getNeighbors(v1)) di++;
					for (V k : graph.getNeighbors(v2)) dj++;
					
//					if(cluster.getNeighbors(v1).contains(v2))
//						T += 1 ;
//					ET += di*dj/(2*m) ;
					
					for (V k : cluster.getNeighbors(v1)){
						dk=0;
						for (V kk : graph.getNeighbors(k)) dk++;
					
						if (cluster.getNeighbors(v2).contains(k))//	if (cluster.containsVertex(k))
								T += 1;							
						ET += (di*dk/(2*m))*(dj*dk/(2*m));
					}
				}
			}
		}
		triularity = (T-ET)/(MaxT-ET);
		return triularity;
	}
	
	public String toString(){
		return "Triularity";
	}

	public String getName() {
		return "T" ;
	}

}
