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

public class M1Modularity<V, E> extends RelativeCommunityCriteria<V,E> {

	public M1Modularity(Graph<V, E> graph) {
		super(graph);
	}

	//TODO: make it weighted !!!! edge count = sum W, size = sum W, 1 = W 
	public double evaluateCommunities(Vector<Graph<V,E>> communities){
		
		double modularity = 0;
		double Max , E=0; 
		
		Max = 2* graph.getEdgeCount();

		for (Graph<V,E> cluster : communities) {
		
			for (V v1 : cluster.getVertices()) {
				for (V v2 : cluster.getVertices()) {
					double pi=0,pj=0;
					
					for (V k : graph.getNeighbors(v1)) pi++;
					for (V k : graph.getNeighbors(v2)) pj++;
				
					if (cluster.getNeighbors(v1).contains(v2))
						modularity += 1;

					E += pi*pj / Max ;				
				}
			}
			
		}
		modularity = (modularity - E ) / (Max - E);//	modularity = (modularity - E ) / (Max); Original formula
		//if(modularity>1 || modularity<0  ) System.err.println("------------------------------------ >>>>>>  "+modularity);
		return modularity;
	}
	
	public String toString(){
		return "M1Q";
	}

	public String getName() {
		return "M1Q" ;
	}

}
