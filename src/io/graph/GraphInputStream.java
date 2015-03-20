package io.graph;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Factory;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;

public abstract class GraphInputStream <V,E> {
	 /**
     * Mapping of the Object created for this vertex and its string representation in the input file.
     */
	protected String commentIndicator = "*";
	protected String tokenizationPattern = "[\\s]+";
	
	protected int edgeIds = 1, nodeIds = 1;
	protected Graph<V, E> graph;
	protected Factory<V> nodeFactory;
	protected Factory<E> edgeFactory;
	protected BufferedReader br;
	
	protected HashMap<V, HashMap<Object, Vector<Object>>> nodeAttributes = new HashMap<V, HashMap<Object,Vector<Object>>>();

	protected Map<V, String> vertex_labels = new HashMap<V,String>();
	protected Map<E, Double> weights;

	public HashMap<V, HashMap<Object, Vector<Object>>> getNodeAttributes() {
		return nodeAttributes;
	}
	public Map<E, Double> getWeights() {
		return weights;
	}

	/**
     * Mapping of string representation of vertices in the input file and their objects: reverse of vertex_labels.
     */
	protected Map<String, V> labels_vertices = new HashMap<String, V>();
	public V transform(String vertex){
		return labels_vertices.get(vertex);
	}
	
	
	public Map<String, V> getLabels_vertices() {
		return labels_vertices;
	}


	/**
	 * @param inputStream
	 * @param nodeFactory
	 * @param edgeFactory
	 * @return graph from the input stream, if the nodeFactory or edgeFactory are null, it would assumes the nodes and edges are Integer and 
	 * It would constructs the vertex IDs equal to the exact string that was parsed in and appeared in the file
	 * @throws IOException
	 */
	public Graph<V,E> readGraph(BufferedReader br, Factory<V> nodeFactory ,Factory<E> edgeFactory) throws IOException{
		this.nodeFactory = nodeFactory;
		this.edgeFactory = edgeFactory;
		this.br = br;
		
		init();
		String tmp ;
		while( (tmp=br.readLine())!=null ){
			if(tmp.length()>0)
				parse(tmp);
		}
		return graph;	
	}
	
	protected void init(){
		graph = new SparseMultigraph<V, E>();
		edgeIds = nodeIds = 1;
	}

	protected abstract void parse(String tmp) throws IOException;
	
	public Graph<V,E> readGraph(InputStream inputStream, Factory<V> nodeFactory ,Factory<E> edgeFactory) throws IOException{
		return  readGraph( new BufferedReader(new InputStreamReader(inputStream)), nodeFactory, edgeFactory);
	}
		
	public  Graph<V,E> readGraph(InputStream inputStream) throws IOException{
		return  readGraph( inputStream, null, null);
	}
	public  Graph<V,E> readGraph(String path) throws IOException{
		return  readGraph( new BufferedReader(new InputStreamReader(new FileInputStream(new File(path)))), null, null);
	}
	
	
	protected Object parseValue(String v){
		try{Integer d = Integer.parseInt(v);return d;}catch(Exception exception){}
		try{Double d = Double.parseDouble(v);return d;}catch(Exception exception){}
		if(v.startsWith("\"")) return v.substring(1, v.length()-1);
		return v;
	}

	protected V getAddVertex(String label){
		V v = labels_vertices.get(label);
		if (v!=null) return v;
		
		if(nodeFactory != null)
			v = nodeFactory.create();
		else {
			if(label!=null)
				try{v = (V)parseValue(label);}catch(Exception exception){}
			if(v==null)
				v = (V) new Integer(nodeIds++);
		}
		graph.addVertex(v);
		vertex_labels.put(v,label);
		labels_vertices.put(label, v);
		return v;
	}
	protected E getAddEdge(V v1, V v2, Object w){
		return getAddEdge(v1, v2, w, null);
	}
	@SuppressWarnings("unchecked")
	protected E getAddEdge(V v1, V v2, Object w, Object label){
		E e = null;
		if(graph.getNeighbors(v1) != null && graph.getNeighbors(v1).contains(v2))
			return graph.findEdge(v1, v2); 
		
		if(edgeFactory != null)
			e = edgeFactory.create();
		else {
			if(label!=null)
				try{e = (E)parseValue(label.toString());}catch(Exception exception){}
			if(e==null)
				e = (E) new Integer(edgeIds++);
		}
		graph.addEdge(e,v1,v2);
		if(w!=null){
			if (weights==null)
				weights = new HashMap<E, Double>(); 
			weights.put(e, new Double(w.toString()));
		}
		return e;
	}
	

	
}
