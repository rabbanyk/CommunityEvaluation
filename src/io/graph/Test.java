package io.graph;

import io.graph.gml.GMLGraphReader;
import io.graph.gml.GMLGraphWriter;
import io.graph.pairs.PairsGraphReader;
import io.graph.pairs.PairsGraphWriter;
import io.graph.pajek.PajekGraphReader;
import io.graph.pajek.PajekGraphWriter;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.TransformerUtils;

import data.GraphDataSet;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

public class Test {

	public static void testWriters(){
		GraphDataSet<String, String>  dataset = new GraphDataSet<String, String>("test");
		dataset.graph = new SparseGraph<String, String>() ;
		dataset.graph.addEdge("ab","a", "b");
		dataset.graph.addEdge("bc","b", "c");
		dataset.graph.addEdge("ac","a", "c");
		dataset.graph.addEdge("ad", "a", "d");
		dataset.graph.addEdge("de","d", "e");

		dataset.setWeight("ab", 1.);
		dataset.setWeight("bc", 2.);
		dataset.setWeight("ac", 3.);
		dataset.setWeight("ad", 4.);
		dataset.setWeight("de", 5.);
		
		dataset.addAttribute("a", "label" ,"a");
		dataset.addAttribute("b", "label" ,"b");
		dataset.addAttribute("c", "label" ,"c");
		dataset.addAttribute("d", "label" ,"d");
		dataset.addAttribute("e", "label" ,"e");
		
		dataset.addAttribute("a", "att1" ,1);
		dataset.addAttribute("b", "att1" ,2);
		dataset.addAttribute("a", "att1" ,3);
		dataset.addAttribute("d", "att1" ,4);
		dataset.addAttribute("e", "att1" ,5);
		
		Map<String, String> vertex_Ids = new  HashMap<String,String>();
		vertex_Ids.put("a", ""+1);
		vertex_Ids.put("b", ""+2);
		vertex_Ids.put("c", ""+3);
		vertex_Ids.put("d", ""+4);
		vertex_Ids.put("e", ""+5);
		for(String v : vertex_Ids.keySet()){
			dataset.addAttribute(v,"id", vertex_Ids.get(v));
		}

		System.err.println(dataset.attributes);

		try {
			GraphOutputStream<String, String> graphWriter = new PairsGraphWriter<String, String>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.wpairs", dataset.graph ,  dataset.getAttMap("label") , dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/testGraphF.wpairs", dataset);
//			graphWriter.writeGraph("dev_temp/testGraphF.wpairs", dataset, vertex_Ids);

			graphWriter = new GMLGraphWriter<String, String>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset.graph ,  dataset.getAttMap("label") ,  dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset);
//			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset, vertex_Ids);
//			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset, dataset.getAttMap("label"));



			graphWriter = new PajekGraphWriter<String, String>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.net", dataset.graph ,   dataset.getAttMap("label") ,  dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/testGraphF.net", dataset);
//			graphWriter.writeGraph("dev_temp/testGraphF.net", dataset, vertex_Ids);
//			graphWriter.writeGraph("dev_temp/testGraphF.net", dataset, dataset.getAttMap("label"));

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static<V,E> void testReaders(){
		try {
			Graph< V, E> g ;
			GraphDataSet<V,E> dataset=null;
		
			
			GraphInputStream<V,E> graphReader = new PairsGraphReader<V,E>()	;
			g = graphReader.readGraph("dev_temp/testGraphF.wpairs");
			System.err.println(graphReader.getNodeAttributes());
			
			GraphOutputStream<V,E> graphWriter = new PairsGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testGraphF_rew.wpairs", g , graphReader.vertex_labels ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());
		
			graphReader = new GMLGraphReader<V,E>();
			g = graphReader.readGraph("dev_temp/testGraphF.gml");
			System.err.println(graphReader.getNodeAttributes());

			graphWriter = new GMLGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testGraphF_rew.gml", g , graphReader.vertex_labels  ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());

			graphReader = new PajekGraphReader<V,E>();
			g = graphReader.readGraph("dev_temp/testGraphF.net");
			System.err.println(graphReader.getNodeAttributes());

			graphWriter = new PajekGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testGraphF_rew.net", g , graphReader.vertex_labels  ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	//Testing
	public static void main(String[] args){
		//TODO: debug
		testWriters();
		testReaders();
	
	}
}
