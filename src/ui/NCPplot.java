package ui;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import algorithms.topleaders.Partitioning;
import edu.uci.ics.jung.graph.Graph;

public class NCPplot <V,E>{
	public Vector<Vector<Double>> ncpPlot(Partitioning<V> partitioning, Graph<V, E> graph){
		Vector<Vector<Double>> result = new Vector<Vector<Double>>();
		
		double tmp =0;
		for(Set<V> community: partitioning.getCommunities()){
			while(community.size()>= result.size()){
				result.add(new Vector<Double>());
			}
			
			tmp = conductance(community,graph);
			int i = 0;
			while(i<result.get(community.size()).size() && result.get(community.size()).get(i)>tmp){
				i++;
			}
			result.get(community.size()).add(i, tmp);
		}
		
		return result;
	}
	
	public double conductance(Set<V> set, Graph<V, E> graph){
		double res = 0, div = 0;
		
		//TODO: change it
		for(V v:set){
			Set<V> neighs = new HashSet<V>(graph.getNeighbors(v));
			div += neighs.size();
			neighs.retainAll(set);
			res += neighs.size();
		}
		
		div = Math.max(2*graph.getEdgeCount()-div, div);
		
		return res/div;
	}
	
	
	public Vector<Double> ncpPlotLine(Partitioning<V> partitioning, Graph<V, E> graph){
		Vector<Vector<Double>> tmp = ncpPlot(partitioning, graph);
		Vector<Double>  res = new Vector<Double>();
		for(int i=0;i<tmp.size();i++){
			if(tmp.get(i).size()>0)
				res.add(tmp.get(i).get(0));
			else res.add(0.);
		}
		return res;
	}

	public void printNcpPlot(Partitioning<V> partitioning, Graph<V, E> graph){
		Vector<Vector<Double>> tmp = ncpPlot(partitioning, graph);
		for(int i = 0; i< tmp.size(); i++){
			for (int j = 0; j< tmp.get(i).size(); j++){
				System.out.println(i +"\t"+tmp.get(i).get(j));
			}
		}
	}
}

