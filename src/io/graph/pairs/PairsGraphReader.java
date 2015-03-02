package io.graph.pairs;

import io.graph.GraphInputStream;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.SparseMultigraph;


public  class PairsGraphReader<V,E> extends GraphInputStream<V, E>{
	
	public Graph<V, E> readGraph(InputStream inputStream, Factory<V> nodeFactory,	Factory<E> edgeFactory) throws IOException {
		Graph<V, E> graph = new SparseMultigraph<V, E>();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
		
		V v1 = null ,v2 = null;
		int edgeIds = 1;
		
		String tmp =  br.readLine() , l1,l2, w;
		while( tmp!=null && tmp.length()>0){
			tmp = tmp.trim();
			w = null;
			
			int ind = (tmp.indexOf('\t'));
			if(ind == -1) ind = tmp.indexOf(' ');
			
			l1 = tmp.substring(0, ind).trim();

			tmp = tmp.substring(ind+1).trim();

			ind = (tmp.indexOf('\t'));
			if(ind == -1) ind = tmp.indexOf(' ');
			
			l2 = ((ind!=-1)?tmp.substring(0, ind):tmp).trim();
			
			if(ind!=-1 && tmp.length() > ind+1){
				tmp = tmp.substring(ind+1).trim();
			
				if(tmp.length()>0){ //weighted
					ind = (tmp.indexOf('\t'));
					if(ind == -1) ind = tmp.indexOf(' ');
					
					w = ((ind!=-1)?tmp.substring(0, ind):tmp).trim();
					//if(weights == null) weights = new HashedMap<E, Double>();
				}
			}
//			System.err.println(l1+" "+l2);
			
			v1 = labels_vertices.get(l1);
			v2 = labels_vertices.get(l2);
			
			if(v1 == null){
				if(nodeFactory != null)
					v1 = nodeFactory.create();
				else {
					try{
						v1 = (V)(new Integer(Integer.parseInt(l1)));
					}catch(NumberFormatException n){
						v1 = (V)(l1);
					}
					}
				graph.addVertex(v1);
				vertex_labels.put(v1,l1);
				labels_vertices.put(l1, v1);
			}
			if(v2 == null){
				if(nodeFactory != null)
					v2 = nodeFactory.create();
				else {
					try{
						v2 = (V)(new Integer(Integer.parseInt(l2)));
					}catch(NumberFormatException n){
						v2 = (V)(l2);
					}
					}
				graph.addVertex(v2);
				vertex_labels.put(v2,l2);
				labels_vertices.put(l2, v2);
			}
			
			// Prohibit repeated edges, LFR benchmarks have 2 lines per each edge
			// TODO: should be change for the case of multigraphs
			if(graph.getNeighbors(v1) != null && !graph.getNeighbors(v1).contains(v2)){
				E e = (edgeFactory != null)?edgeFactory.create():(E)(new Integer(edgeIds++));
				graph.addEdge(e,v1,v2);
				if(w!= null){
					if (weights==null)
						weights = new HashMap<E, Double>(); 
					weights.put(e, ( new Double(w))  );
				}
				
			}
			
			tmp = br.readLine();
		}
		
		return graph;
	}

	

}
