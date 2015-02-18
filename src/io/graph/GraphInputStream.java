package io.graph;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.collections15.Factory;

import edu.uci.ics.jung.graph.Graph;

public abstract class GraphInputStream <V,E> {
	 /**
     * Mapping of the Object created for this vertex and its string representation in the input file.
     */
	protected HashMap<V, HashMap<Object, Object>> nodeAttributes = new HashMap<V, HashMap<Object,Object>>();

	protected Map<V, String> vertex_labels = new HashMap<V,String>();
	protected Map<E, Double> weights;

	public String transform(V vertex){
		return vertex_labels.get(vertex);
	}
	public HashMap<V, HashMap<Object, Object>> getNodeAttributes() {
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
	
	
	public Map<V, String> getVertex_labels() {
		return vertex_labels;
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
	public abstract Graph<V,E> readGraph(InputStream inputStream, Factory<V> nodeFactory ,Factory<E> edgeFactory) throws IOException;
	
	public  Graph<V,E> readGraph(InputStream inputStream) throws IOException{
		return  readGraph( inputStream, null, null);
	}
	public  Graph<V,E> readGraph(String path) throws IOException{
		return  readGraph( new FileInputStream(new File(path)), null, null);
	}
}
