package util;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.collections15.Factory;
import org.xml.sax.SAXException;

import data.GraphDataSet;
import algorithms.communityMining.data.Grouping;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;
import edu.uci.ics.jung.io.GraphMLReader;
import io.graph.GraphInputStream;
import io.graph.GraphOutputStream;
import io.graph.gml.GMLGraphReader;
import io.graph.gml.GMLGraphWriter;
import io.graph.pairs.PairsGraphReader;
import io.graph.pairs.PairsGraphWriter;
import io.graph.pajek.PajekGraphReader;
import io.graph.pajek.PajekGraphWriter;
import io.group.GroupingReader;
import io.group.GroupingWriter;
import io.group.ListGroupingWriter;
import io.group.PairGrouingReader;

public class IOUtils  {
	
	enum Type {pairs, edges, gml, graphml, net, dat, pajek, plosone};
	
	public static <V,E> Graph<V,E> load (String path) throws FileNotFoundException, IOException{
		return IOUtils.<V,E>getReader(path).readGraph(new FileInputStream(path));
	}
	public static <V,E> Graph<V,E> load (Type type, String path) throws FileNotFoundException, IOException{
		return IOUtils.<V,E>getReader(type).readGraph(new FileInputStream(path));
	}
	
	protected static <V,E> boolean isKnownFile(String extension){
		for(Type t : Type.values())
			if(t.toString().equals(extension))
				return true;
		return false;
	}
	public static <V,E> GraphInputStream<V, E> getReader (String file){
		String extension = file.substring(file.lastIndexOf('.')+1);
//		System.err.println(extension +"   "+ isKnownFile(extension));
		if (!isKnownFile(extension)) return null;
		return getReader(Type.valueOf(extension));
	}
	
	
	public static <V,E> GraphDataSet<V,E> loadGraphDataset (String path) throws FileNotFoundException, IOException{
		GraphDataSet< V,  E> res = new GraphDataSet<V, E>(path);
		GraphInputStream<V, E> reader = IOUtils.<V,E>getReader(path);
		res.graph = reader.readGraph(new FileInputStream(path));
		res.weights = reader.getWeights();
		res.attributes = reader.getNodeAttributes();
		res.labels_vertices = reader.getLabels_vertices();
//		reader.vertex_labels;
		return res;
	} 
	
	public static <V,E> void  write (String path, Graph<V,E> graph) throws FileNotFoundException, IOException{
		write(path, graph, null);
	}
	public static <V,E> void  write (String path, GraphDataSet<V,E> graph) throws FileNotFoundException, IOException{
		write(path, graph, null);
	}
	
	public static <V,E> void  write (String path, Graph<V,E> graph, Map<V, String> vertex_Ids) throws FileNotFoundException, IOException{
		IOUtils.<V,E>getWriter(path).writeGraph(path, graph, vertex_Ids);
	}
	public static <V,E> void  write (String path, GraphDataSet<V,E> graph, Map<V, String> vertex_Ids) throws FileNotFoundException, IOException{
		IOUtils.<V,E>getWriter(path).writeGraph(path, graph, vertex_Ids);
	}
	
	public static <V,E> GraphOutputStream<V, E> getWriter (Type type){
		switch (type) {
		case pairs: case dat: case edges: return new PairsGraphWriter<V, E>();
		case gml: return new GMLGraphWriter<V, E>();
		case net: case pajek: return new PajekGraphWriter<V, E>();
		case plosone: return new PajekGraphWriter<>(";", false,"#");
		}
		return null;
	}
	
	public static <V,E> GraphOutputStream<V, E> getWriter (String file){
		String extension = file.substring(file.lastIndexOf('.')+1);
		if (!isKnownFile(extension)) return null;
		return getWriter(Type.valueOf(extension));
	}
	
	@SuppressWarnings("unchecked")
	public static <V,E> GraphInputStream<V, E> getReader (Type type){
		switch (type) {
		case pairs: case dat: case edges: return new PairsGraphReader<V, E>();
		case gml: return new GMLGraphReader<V, E>();
		case net: case pajek: return new PajekGraphReader<V, E>();
		case plosone: return new PajekGraphReader<>(false, "[;\\s]+","#");
		case graphml: try {
				return new GraphInputStream<V, E>() {
					//TODO: do we need to give it the vertex and edge factory for consistency?
					GraphMLReader jungReader  = new GraphMLReader();
					@Override
					public Graph<V, E> readGraph(InputStream inputStream,
							Factory<V> nodeFactory, Factory<E> edgeFactory)
							throws IOException {
						Graph<V, E> graph = new SparseMultigraph<V, E>();
						jungReader.load(new InputStreamReader(inputStream), graph);
						//this.weights = jungReader.getEdgeDescriptions(); TODO: weights?
						return graph;
					}
					@Override
					protected void parse(String tmp) throws IOException {
					}
				};
			} catch (ParserConfigurationException | SAXException e) {
				e.printStackTrace();
			}
		default:
			try {
//				System.err.println(type.toString());
				return (GraphInputStream< V, E> ) Class.forName("io.graph.*."+type.toString()).newInstance();
			} catch (Exception e) {
				//e.printStackTrace();
			}		}
		return null;
	}
	
	
	public static<V,E> void writeGML(String path, Graph<V, E> graph,
			final Map<V, String> vertex_Ids, final  Map<E, ? extends Number> weights,	 
			String weightsOutputFormat, Map<V, HashMap<Object,Vector<Object>>> vertex_attributes) throws IOException {
		
		GMLGraphWriter<V,E> graphWriter = new GMLGraphWriter<V, E>();
		graphWriter.writeGraph(path, graph ,vertex_Ids, 
				weights, weightsOutputFormat ,vertex_attributes);
	}
	
	public static<V> void writeListGrouping(String path , Grouping<V> U, Map<V, String> vertex_labels) throws IOException {
		GroupingWriter<V> groupingWriter = new ListGroupingWriter<>();
		groupingWriter.writeGrouping(path , U, vertex_labels); 
	}
	
	
	
	
//	public static <V> GroupingReader<V> getGTReader (String file){
//		String filename = file.substring(0,file.lastIndexOf('.'));
////		System.err.println(extension +"   "+ isKnownFile(extension));
//		if (! new File(filename+".gt").exists()) return null;
//		return new PairGrouingReader<V>(true);
//	}
	
	
	
}
