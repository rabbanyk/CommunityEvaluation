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

import util.IOUtils;
import data.GraphDataSet;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

public class Test {

	private static GraphDataSet<String, String> getTestDataset(){
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
		
		return dataset;
	}


	public static <V,E> void test(GraphDataSet<V,E>  dataset){
		dataset.print();
		Map<V, String> vertex_Ids = null;
//				new  HashMap<V,String>();
//		int counter = 1;
//		for (V v : dataset.graph.getVertices())
//			vertex_Ids.put(v, ""+(counter++));
	
		try {
			IOUtils.write("dev_temp/testg.pairs", dataset, vertex_Ids);
			IOUtils.write("dev_temp/testg.net", dataset, vertex_Ids);
			IOUtils.write("dev_temp/testg.gml", dataset, vertex_Ids);
			IOUtils.write("dev_temp/testg.plosone", dataset, vertex_Ids);

			dataset = IOUtils.loadGraphDataset("dev_temp/testg.pairs");
			IOUtils.write("dev_temp/testg_n.pairs", dataset, dataset.getLabels());
			dataset.print();

			dataset = IOUtils.loadGraphDataset("dev_temp/testg.net");
			IOUtils.write("dev_temp/testg_n.net", dataset, dataset.getLabels());
			dataset = IOUtils.loadGraphDataset("dev_temp/testg.gml");
			IOUtils.write("dev_temp/testg_n.gml", dataset, dataset.getLabels());
			dataset = IOUtils.loadGraphDataset("dev_temp/testg.plosone");
			dataset.print();
			IOUtils.write("dev_temp/testg_n.plosone", dataset, dataset.getLabels());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <V,E> void testWriters(GraphDataSet<V,E>  dataset){
		Map<V, String> vertex_Ids = new  HashMap<V,String>();
		int counter = 1;
		for (V v : dataset.graph.getVertices())
			vertex_Ids.put(v, ""+(counter++));
	
//		for(V v : vertex_Ids.keySet()){
//			dataset.addAttribute(v,"id", vertex_Ids.get(v));
//		}
//		System.err.println(dataset.attributes);

		try {
			GraphOutputStream<V,E> graphWriter = new PairsGraphWriter<V,E>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.wpairs", dataset.graph ,  dataset.getAttMap("label") , dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/test.wpairs", dataset);
			graphWriter.writeGraph("dev_temp/testIded.wpairs", dataset, vertex_Ids);

			graphWriter = new GMLGraphWriter<V,E>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset.graph ,  dataset.getAttMap("label") ,  dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/test.gml", dataset);
			graphWriter.writeGraph("dev_temp/testIded.gml", dataset, vertex_Ids);
//			graphWriter.writeGraph("dev_temp/testGraphF.gml", dataset, dataset.getAttMap("label"));

			graphWriter = new PajekGraphWriter<V,E>()	;
//			graphWriter.writeGraph("dev_temp/testGraphF.net", dataset.graph ,   dataset.getAttMap("label") ,  dataset.weights, "%.0f", dataset.attributes);
			graphWriter.writeGraph("dev_temp/test.net", dataset);
			graphWriter.writeGraph("dev_temp/testIded.net", dataset, vertex_Ids);
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
			g = graphReader.readGraph("dev_temp/testg.pairs");
			System.err.println(graphReader.getNodeAttributes());
			
			GraphOutputStream<V,E> graphWriter = new PairsGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testg_n.pairs", g , graphReader.vertex_labels ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());
		
			graphReader = new GMLGraphReader<V,E>();
			g = graphReader.readGraph("dev_temp/testg.gml");
			System.err.println(graphReader.getNodeAttributes());

			graphWriter = new GMLGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testg_n.gml", g , graphReader.vertex_labels  ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());

			graphReader = new PajekGraphReader<V,E>();
			g = graphReader.readGraph("dev_temp/testg.net");
			System.err.println(graphReader.getNodeAttributes());

			graphWriter = new PajekGraphWriter<V,E>()	;
			graphWriter.writeGraph("dev_temp/testg_n.net", g , graphReader.vertex_labels  ,
					graphReader.getWeights(), "%.0f",
					graphReader.getNodeAttributes());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	//Testing
	public static void main(String[] args){
		//TODO: debug
		test(getTestDataset());
	}
}
