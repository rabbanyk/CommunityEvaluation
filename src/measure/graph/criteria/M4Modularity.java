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

public class M4Modularity<V, E> extends RelativeCommunityCriteria<V,E> {

	public M4Modularity(Graph<V, E> graph) {
		super(graph);
	}

	public double evaluateCommunities(Vector<Graph<V,E>> communities){
		
		double modularity = 0;
		double Max , E=0; 
		
		Max =  graph.getVertexCount(); //Original formula

		for (Graph<V,E> cluster : communities) {
		
			for (V v1 : cluster.getVertices()) {
				for (V v2 : cluster.getVertices()) if(v1!=v2){
					double pij = 0, pi=0,pj=0;
					
					for (V k : graph.getNeighbors(v1)) pi++;	pi/=Max;
					for (V k : graph.getNeighbors(v2)) pj++;	pj/=Max;
				
					if (cluster.getNeighbors(v1).contains(v2))  pij += 1;
					pij/=Max;
						
					modularity += pij;
					E += pi*pj /Max;				
				}
			}
			
		}
		modularity  = (modularity - E ) / (1-E);// (modularity - E );  Original formula
		return modularity;
	}
	
	public String toString(){
		return "M4Q";
	}

	public String getName() {
		return "M4Q" ;
	}

}
