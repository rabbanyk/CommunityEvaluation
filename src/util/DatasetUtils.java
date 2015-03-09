package util;
import io.Logger;
import io.graph.GraphInputStream;
import io.graph.gml.GMLGraphReader;
import io.graph.gml.GMLGraphWriter;
import io.graph.pairs.PairsGraphReader;
import io.group.ListGrouingReader;
import io.group.PairGrouingReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.data.Grouping;
import data.GraphDataSet;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

public class DatasetUtils{
	
	static String datasetLocation = "./Datasets/";
	public static enum ClassicDataset {karate, wkarate, football, polbooks, polblogs};//, strike};
	public static enum DummyDataset {Structure, Omega, Weights, NMIexample};
	
	public static GraphDataSet<Integer, Integer> loadClassic(ClassicDataset datasetName){
		GraphDataSet<Integer, Integer> dataset=null;
		try {
			GMLGraphReader<Integer, Integer> reader = new GMLGraphReader<Integer, Integer>();
			dataset = new GraphDataSet<Integer, Integer>(datasetName.toString());
			dataset.graph = (Graph<Integer, Integer>) reader.readGraph(datasetLocation + "classics/"+datasetName.toString()+".gml");
			dataset.weights = reader.getWeights();
			dataset.attributes = reader.getNodeAttributes();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return dataset;
	}
	
	public static GraphDataSet<Integer, Integer> loadDummy(DummyDataset datasetName){
		
		GraphDataSet<Integer, Integer> dataset=new GraphDataSet<Integer, Integer>(datasetName.toString());
		dataset.attributes = new HashMap<Integer, HashMap<Object,Vector<Object>>>();
		Graph<Integer, Integer> g =  new SparseGraph<>() ;
		dataset.graph = g;
		int c=0;
		switch (datasetName) {
		case Structure:
			g.addEdge(c++,0,1);g.addEdge(c++,0,5);g.addEdge(c++,0,6);
			g.addEdge(c++,1,2);g.addEdge(c++,1,5);
			g.addEdge(c++,2,3);g.addEdge(c++,2,5);
			g.addEdge(c++,3,5);g.addEdge(c++,3,4);
			g.addEdge(c++,4,5);
			g.addEdge(c++,5,6);g.addEdge(c++,5,8);
			g.addEdge(c++,6,7);g.addEdge(c++,6,8);
			g.addEdge(c++,7,8);
			dataset.addPartitioingAttribute("V",createClusteringFromArray(new int[][]{ {0,1,2,3,4,5},{6,7,8}}));
			dataset.addPartitioingAttribute("U1",createClusteringFromArray(new int[][]{ {1,2,3,4,5},{0,6,7,8}}));
			dataset.addPartitioingAttribute("U2",createClusteringFromArray(new int[][]{ {0,1,2,3,4},{5,6,7,8}}));
			break;
		case Omega:
			g.addEdge(c++,0,1);g.addEdge(c++,0,4);
			g.addEdge(c++,1,2);
			g.addEdge(c++,2,3);
			g.addEdge(c++,3,4);
			dataset.addPartitioingAttribute("V", createClusteringFromArray(new int[][]{ {0,3,4},{1,2,3},{2,3,4}}));
			dataset.addPartitioingAttribute("U1", createClusteringFromArray(new int[][]{ {0,3,4},{1},{2}}));
			dataset.addPartitioingAttribute("U2", createClusteringFromArray(new int[][]{ {0,3,4},{1},{2,3}}));
			break;
		case Weights:
			dataset.weights = new HashMap<Integer, Double>();
			dataset.weights.put(c, 1.0);g.addEdge(c++,0,1);dataset.weights.put(c, 1.0);g.addEdge(c++,0,5);dataset.weights.put(c, 15.0);g.addEdge(c++,0,6);
			dataset.weights.put(c, 1.0);g.addEdge(c++,1,2);dataset.weights.put(c, 1.0);g.addEdge(c++,1,5);
			dataset.weights.put(c, 1.0);g.addEdge(c++,2,3);dataset.weights.put(c, 1.0);g.addEdge(c++,2,5);
			dataset.weights.put(c, 1.0);g.addEdge(c++,3,5);dataset.weights.put(c, 1.0);g.addEdge(c++,3,4);
			dataset.weights.put(c, 1.0);g.addEdge(c++,4,5);
			dataset.weights.put(c, 5.0);g.addEdge(c++,5,6);
			dataset.weights.put(c, 1.0);g.addEdge(c++,6,7);dataset.weights.put(c, 1.0);g.addEdge(c++,6,8);
			dataset.weights.put(c, 1.0);g.addEdge(c++,7,8);
			break;
		case NMIexample:
			g.addEdge(c++,1,2);g.addEdge(c++,1,4);g.addEdge(c++,1,5);
			g.addEdge(c++,2,3);g.addEdge(c++,2,4);g.addEdge(c++,2,5);
			g.addEdge(c++,3,4);g.addEdge(c++,4,5);
			g.addEdge(c++,6,7);g.addEdge(c++,6,8);g.addEdge(c++,7,8);
			g.addEdge(c++,8,9);
			g.addEdge(c++,9,10);g.addEdge(c++,9,11);g.addEdge(c++,10,11);
			dataset.addPartitioingAttribute("V", createClusteringFromArray(new int [][]{{1, 2, 3 ,4 ,5 },{6 ,7, 8},{9 ,10 ,11}}));
			dataset.addPartitioingAttribute("U1", createClusteringFromArray(new int [][]{{1 ,2 ,3 ,4, 5}, {6, 7 ,8, 9, 10, 11}}));
			dataset.addPartitioingAttribute("U2", createClusteringFromArray(new int [][]{{1 ,2 ,3},{6, 7 ,8},{4, 5, 9, 10, 11}}));
			break;
		default:
			break;
		}
		
//		try {
//			IOUtils.<Integer,Integer>writeGML("examples/"+dataset.name+"_attributed.gml",
//					dataset.graph, null, dataset.weights, null, dataset.attributes);
//			IOUtils.<Integer,Integer>writeGML("examples/"+dataset.name+".gml",
//					dataset.graph, null, dataset.weights, null, null);
//			
//			IOUtils.<Integer>writeListGrouping("examples/"+dataset.name+"_V.list", dataset.getGrouping("V"), dataset.getLabels());
//			IOUtils.<Integer>writeListGrouping("examples/"+dataset.name+"_U1.list", dataset.getGrouping("U1"), dataset.getLabels());
//			IOUtils.<Integer>writeListGrouping("examples/"+dataset.name+"_U2.list", dataset.getGrouping("U2"), dataset.getLabels());
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		return dataset;

	}
	public static Vector< Set<Integer>> createClusteringFromArray(int[][] aa){
		Vector<Set<Integer>> clus = new Vector<Set<Integer>>(); 

		for (int[] a : aa) {
			Set<Integer> c = new HashSet<Integer>();
			for (int i : a) {
				c.add(i);
			}
			clus.add(c);
		}
		return clus;
	}
	
	public static Vector<GraphDataSet<Integer, Integer>> loadAllClassics(){
		Vector<GraphDataSet<Integer, Integer>> datasets = new Vector<GraphDataSet<Integer, Integer>>();
		for(ClassicDataset dataset: ClassicDataset.values()){
			datasets.add(loadClassic(dataset));
		}
		return datasets;
	}
	public static <V,E> GraphDataSet<V, E> load(String network){
		return load(new File(network));
	}
	public static <V,E> GraphDataSet<V, E> load(File network){
		try {
			FileInputStream fgraph = new FileInputStream(network);
			GraphInputStream<V, E> graphReader =  IOUtils.getReader(network.getName());
			if(graphReader != null){
				System.err.println("loading " + network.getName()+"...");
//				System.err.println(network.getName());
//				System.err.println(graphReader);
				final GraphDataSet<V, E> dataSet = new GraphDataSet<V, E>(network.getName());
				dataSet.graph = graphReader.readGraph(fgraph);
				dataSet.weights= graphReader.getWeights();
				dataSet.attributes = graphReader.getNodeAttributes();
				dataSet.labels_vertices = graphReader.getLabels_vertices();
				Transformer<String, V> vertexTransformer = new Transformer<String, V>() {
					public V transform(String arg0) {
					return dataSet.labels_vertices.get(arg0);
				}};
				//Read GT if exists
				{
					String filename = network.getPath();
					filename = filename.substring(0,filename.lastIndexOf('.'));
//					System.err.println(filename);

					File gt = new File(filename+".gt");
//					System.err.println(gt.getName());
//					final Vector<V> vertexes = new Vector<V>(dataSet.graph.getVertices());
//					for(int i = 0; i< vertexes.size(); i++){
//						vertexLabels.put(vertexes.elementAt(i), i+startIndexId);
//					}
//					new Transformer
					if (gt.exists()){
//						System.err.println("  is here  ");
//						final Map<String, V> Labels_vertices = graphReader.getLabels_vertices();
//						Transformer<String, V> vertexTransformer = new Transformer<String, V>() {
//							public V transform(String arg0) {
//								return Labels_vertices.get(arg0);
//							}
//						};
						Grouping<V> groundT = (new PairGrouingReader<V>(false)).
								readPartitioning(new FileInputStream(gt),vertexTransformer);
						System.err.println("> found ground-truth with "+groundT.getNumberOfGroups()+" clusters for " + filename);
						dataSet.addPartitioingAttribute("value", groundT.getGroups());
					}else{
					 gt = new File(filename+".lgt");
					if (gt.exists()){
//						System.err.println("  is here  ");
//						final Map<String, V> Labels_vertices = graphReader.getLabels_vertices();
//						Transformer<String, V> vertexTransformer = new Transformer<String, V>() {
//							public V transform(String arg0) {
//								return Labels_vertices.get(arg0);
//							}
//						};
						Grouping<V> groundT = (new ListGrouingReader<V>()).
								readPartitioning(new FileInputStream(gt),vertexTransformer);
						System.err.println("> found ground-truth with "+groundT.getNumberOfGroups()+" clusters for " + filename);
						dataSet.addPartitioingAttribute("value", groundT.getGroups());
					}
					}
			
				}
				return dataSet;
			}
			fgraph.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	final static int MAX_EDGE_LIMIT= 1000000;

	
	
	public static <V,E> Vector<GraphDataSet<V, E>> loadAllDataSets(String path){
		return loadAllDataSets( new File(path));
	}
	public static  Vector<File> loadAllNames(String path){
		return loadAllNames( new File(path));
	}
	public static<V,E>   Iterator<GraphDataSet<V, E>> loadAll(String path){
		final Vector<File> datasets = loadAllNames(path);
		return  new Iterator<GraphDataSet<V,E>>() {
			int index = 0;
			@Override				

			public boolean hasNext() {
				return index<datasets.size();
			}
			@Override
			public GraphDataSet<V, E> next() {
				GraphDataSet<V, E> dataset = null;
				while (hasNext() && (dataset==null || dataset.graph.getEdgeCount()>MAX_EDGE_LIMIT )){
					try{
						dataset = DatasetUtils.<V,E>load(datasets.get(index));
						index++;
					}catch(Exception e){
						System.err.println("Could not read " + datasets.get(index-1));
//						e.printStackTrace();
					}
				}
				return dataset;
			}
			@Override
			public void remove() {
			}
		};
	}
	
	public static <V,E> Vector<GraphDataSet<V, E>> loadAllDataSets(File directory){
		Vector<GraphDataSet<V, E>> datasets = new Vector<GraphDataSet<V, E>>();
		for (File file : loadAllNames(directory)) {
			GraphDataSet<V, E> dataSet=null;
			try{
				dataSet = DatasetUtils.<V,E>load(file);
			}catch(Exception e){
				System.err.println("Could not read " + file);
				e.printStackTrace();
			}
			if (dataSet!=null ) datasets.add(dataSet);
		}
		return datasets;
	}
	public static  Vector<File> loadAllNames(File directory){
		Vector<File> datasets = new Vector<File>();
		if (!directory.canRead()) {
			System.err.println(" Can not read " + directory );
			return datasets;
		}
		for (File network : directory.listFiles())  
			if(network.isFile()){
					datasets.add(network);
			}else if(network.isDirectory()){
				datasets.addAll(DatasetUtils.loadAllNames(network));
			}
		Collections.sort(datasets);
//		, new Comparator<File>(){
//			@Override
//			public int compare(File o1, File o2) {
//				return o1.getName().compareTo(o2.getName());
//			}});
		return datasets;
	}
	
	private static void convertKaratewithCommunitiestoGML(){
		// Convert wkarate.pairs to gml
	
		GMLGraphWriter<Integer, Integer> gmlGraphWriter = new GMLGraphWriter<>();
		PairsGraphReader<Integer, Integer> reader = new PairsGraphReader<Integer, Integer>();
		GraphDataSet< Integer, Integer> dataset = new GraphDataSet<>("wkarate");
	
		try {
			dataset.graph = reader.readGraph("wkarate.pairs");

			System.err.println(dataset.graph.getEdges());
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		System.err.println(dataset.name);
		int[][] coms = {{1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 17, 18, 20, 22} , 
				{9, 10, 15, 16, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34}};
		dataset.attributes = new HashMap<Integer, HashMap<Object,Vector<Object>>>();
	
		for (Integer v :dataset.graph.getVertices()) {
			dataset.addAttribute(v, "label", ""+v);
		}
		for (int c = 0; c < coms.length; c++) {
			for (int v : coms[c]){
				dataset.addAttribute(v, "value", c+1);
			}
		}
	
		try {
			System.err.println(dataset.attributes);
			gmlGraphWriter.writeGraph(	dataset.name+".gml", dataset.graph, null, reader.getWeights(),"%.0f", dataset.attributes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args) {
//		convertKaratewithCommunitiestoGML();
		for (GraphDataSet<Integer, Integer> dataset : DatasetUtils.<Integer,Integer>loadAllDataSets("./Datasets/classics")) {
			dataset.printStats();;
		}
	}

}








