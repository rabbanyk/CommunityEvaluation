package dev_Experiments;

import io.Logger;
import io.graph.GraphInputStream;
import io.graph.pairs.PairsGraphReader;
import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;
import io.group.PairGrouingReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;

import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.data.Partitioning;
import measure.base.RelativeCriteria;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.correlation.Correlation;
import measure.criteria.RelativeClusterCriteria;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class BasicExperiment<V,E> {
	/**
	 * @return collection of data-sets containing ground-truth
	 */
	public Vector<NetworkWithGroundTruth> loadNetworks(String path){

		Vector<NetworkWithGroundTruth> truths = new Vector<NetworkWithGroundTruth>();
		
		File directory = new File(path);//

		Logger.logln(directory);
		
		FilenameFilter LFRnetFilter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.contains("network");
			}
		};
		
		for (File network : directory.listFiles(LFRnetFilter))  {
			try {
				Logger.logln(network.getName());
				System.err.println(network.getName());
				//Loading the Graph
				FileInputStream fgraph = new FileInputStream(network);
				PairsGraphReader<V, E> graphReader = new PairsGraphReader();
				Graph<V,E> graph = graphReader.readGraph(fgraph, null, null);
				Map<E, ? extends Number> weights = graphReader.getWeights();
				
				fgraph.close();
				
							
				//Loading Ground truth
				FileInputStream fgt = new FileInputStream(path  + network.getName().replaceAll("network", "community"));
				PairGrouingReader<V> gtReader = new PairGrouingReader<>();
				Grouping<V>  groundTruth = gtReader.readPartitioning(fgt, graphReader.getLabels_vertices());
				fgt.close();

				
				//Adding to data set
			//	if(isConnected(graph))
					truths.add(new NetworkWithGroundTruth(graph ,( weights==null? null:TransformerUtils.mapTransformer(weights)), groundTruth.getGroups(), network.getName().substring(0,network.getName().indexOf('.'))));
			
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		FilenameFilter pairsFilter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.contains(".pairs");
			}
		};
		
		for (File network : directory.listFiles(pairsFilter))  {
			try {
				Logger.logln(network.getName());
				
				System.err.println(network.getName());
				//Loading the Graph
				FileInputStream fgraph = new FileInputStream(network);
				PairsGraphReader< V, E> graphReader = new PairsGraphReader< V, E>();
				Graph<V,E> graph = graphReader.readGraph(fgraph, null, null);
				Map<E, ? extends Number> weights = graphReader.getWeights();
				fgraph.close();
				
							
				//Loading Ground truth
				FileInputStream fgt = new FileInputStream(path  + network.getName().replaceAll(".pairs", ".gt"));
				ListGrouingReader<V> gtReader = new ListGrouingReader<V>();
				Grouping<V>  groundTruth = gtReader.readPartitioning(fgt, graphReader.getLabels_vertices());
				fgt.close();
				
				//Adding to data set
				truths.add(new NetworkWithGroundTruth(graph, ( weights==null? null:TransformerUtils.mapTransformer(weights)), groundTruth.getGroups(),network.getName().substring(0,network.getName().indexOf('.'))));
			
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		return truths; // synthetic (change mu), real
		
//		int sep = network.getName().indexOf('.');
//		String netname = network.getName().substring(0,sep);//, type = network.getName().substring(sep+1);
//		int netnumber = Integer.parseInt(netname.substring(7));
		
//		if(type.equals(".net") || type.equals(".pajek")){
//		graphReader = new PajekGraphReader<V,E>();
//	}else if(type.equals(".pairs")){
//		graphReader = new PairsGraphReader<V,E>();
//	}else if(type.equals(".gml")){
//		graphReader = new GMLGraphReader<V, E>();
//	groundTruth = ((GMLGraphReader<V,E>)graphReader).groundTruth;//The GML format contains true clustering
//	}
	}

	public Number computeCorrelation(Correlation<Vector<? extends Number>> correlation, Vector<? extends Number> a, Vector<? extends Number> b){
		return correlation.compute(a, b);
	}

	DecimalFormat df = new DecimalFormat("#.###");
	DecimalFormat df2 = new DecimalFormat("#.##");
	DecimalFormat df1 = new DecimalFormat("#.#");

	/**
	 * @param agreementMethod
	 * @param groundTruth
	 * @param results
	 * @return Compute the agreement of each of the <tt>clustering</tt> results against the <tt>groundtruth</tt> based on the given {@link #agreementMethod} method
	 */
	public Vector<Number> externalEvaluation( PartiotioningAgreement<V> agreementMethod,Vector<Set<V>> groundTruth, Vector<Vector<Set<V>>> results){
		Vector<Number> res = new Vector<Number>();
		
		for (Vector<Set<V>> cluster : results) {
			Number agreement = agreementMethod.getAgreement(cluster, groundTruth);
			res.add(agreement);
		}
		
//		System.err.println(  "  AVG("+agreementMethod+")= " + avgK +" +- "+ varK + " in [ " + minK + " , " + maxK + " ] ");
		return res;
	}
	
	/**
	 * @param criteria
	 * @param results
	 * @return computes the score of given criteria on each of the <tt>clustering</tt> results 
	 */
	@SuppressWarnings("unchecked")
	public Vector<Number> internalEvaluation(RelativeCriteria< V> criteria, Vector<Vector<Set<V>>> results){
		Vector<Number> res = new Vector<Number>();
		
		boolean maximizer = true;
		
		if(criteria instanceof RelativeClusterCriteria) 
			maximizer = ((RelativeClusterCriteria)criteria).isMaximizer(); 
		
		double avg = 0;
		for (Vector<Set<V>> cluster : results) {
			res.add(criteria.evaluate(cluster));
			//if((criteria.evaluate(cluster)+"").equals("NaN")) System.err.println("***"+criteria);

			if(!maximizer) avg += res.lastElement().doubleValue();
		}
//		if(maximizer)
//		System.err.println(criteria + " is  MAXIMIZER");

		//flip values around the mean
		if(!maximizer){
//			System.err.println(criteria + " is MINIMIZER");
			avg /= res.size();
			avg*=2;
			
			for (int i = 0; i < res.size(); i++) {
				res.set(i, avg - res.get(i).doubleValue());
			} 
		}
		
		return res;
	}

	boolean isConnected (Graph<V,E> graph){
		//if it is disconnected
		WeakComponentClusterer<V, E> clusterer= new WeakComponentClusterer<V, E>();
		if(clusterer.transform(graph).size()>1) return false;
		return true;
	}
	
	boolean validate (Graph<V,E> graph, Vector<Set<V>> part){
		if (part.size()<2 ) return false;

		//check if it is disjoint
		Vector<V> partitionedVertexes = new Vector<V>();
		for (Set<V> set : part) {
			for (V v : set) {
				if(!partitionedVertexes.contains(v)) partitionedVertexes.add(v);
				else return false;
			}
		}

		//check if it covers all
		if(partitionedVertexes.size()!=graph.getVertexCount())
			return false;
		
		return true;
	}
	
	protected class NetworkWithGroundTruth{
		protected Graph<V,E> graph;
		protected Transformer<E, ? extends Number> weights;  
		protected Vector<Set<V>> groundTruth;
		protected String name;
		
		public NetworkWithGroundTruth(Graph<V, E> graph, Transformer<E, ? extends Number> weights, Vector<Set<V>> groundTruth, String name) {
			super();
			this.graph = graph;
			this.weights = weights;
			this.groundTruth = groundTruth;
			this.name = name;
		}

		public NetworkWithGroundTruth(Graph<V,E> graph, Transformer<E, ? extends Number> weights, Vector<Set<V>> groundTruth) {
			super();
			this.graph = graph;
			this.weights = weights;
			this.groundTruth = groundTruth;
		}

		public Graph<V,E> getGraph() {
			return graph;
		}

		public Vector<Set<V>> getGroundTruth() {
			return groundTruth;
		}
		
		public Transformer<E, ? extends Number> getWeights() {
			return weights;
		}

		public void setWeights(Transformer<E, ? extends Number> weights) {
			this.weights = weights;
		}

		public void setName(String name){
			this.name= name;
		}
		
		public String toString() {
			return name;
		}
		
	}
	
	
}
