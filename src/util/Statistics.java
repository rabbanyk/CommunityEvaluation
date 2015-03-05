package util;

import edu.uci.ics.jung.graph.Graph;

public class Statistics {

	public static<V,E> double[] avgCluteringCoefficient(Graph<V,E> graph){
		double cc = 0, varCC=0;
		for (V v:graph.getVertices()){
			double cv = 0;
			for(V n1: graph.getNeighbors(v)){
				for(V n2: graph.getNeighbors(v)){
					if(n1!=n2 && graph.findEdge(n1, n2)!=null)
						cv++;
				}
			}
			if(cv!=0)
				cv/= graph.degree(v) * (graph.degree(v)-1);
			cc += cv;
			varCC += cv*cv;
		}
		cc/=graph.getVertexCount();
		varCC =Math.sqrt( varCC/graph.getVertexCount() - cc*cc);
	//	System.err.println(cc);
		return new double[]{cc,varCC};
	}
	public static <V,E> double[] getAverageDegree(Graph<V,E> graph){
		double avgD = 0 , varD=0;
		for (V v:graph.getVertices()){
			avgD += graph.degree(v);
			varD += graph.degree(v) * graph.degree(v);
		}
		avgD /= graph.getVertexCount();
		varD = Math.sqrt(varD/graph.getVertexCount() - avgD*avgD);
		return new double[]{avgD,varD};
	}
}
