package dev_Experiments;

import io.Logger;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import algorithms.*;
import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.external_methods.FastModularity;
import algorithms.communityMining.external_methods.Infomap;
import algorithms.communityMining.topleaders.TopLeaders;
import measure.MeasuresUtil;
import measure.cluster.agreement.*;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.cluster.agreement.partitioning.generalization.GAM.Type;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.AdjustionMethod;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.ExternalOverlap;
import measure.correlation.*;
import edu.uci.ics.jung.graph.Graph;
import io.Logger.DebugMode;

public class RandomExpreriment<V,E> extends BasicExperiment<V, E>{
	public String expPath ;//= "exp1/";
	Vector<PartiotioningAgreement<V>> agreementMethods;// = new ARI<V>();
	Vector<Correlation<Vector<? extends Number>>> correlationMethods;// = new PearsonCorrelation(); 

//	public boolean addAgreement(PartiotioningAgreement<V> e) {
//		return agreementMethods.add(e);
//	}

	public boolean addCorrelation(Correlation<Vector<? extends Number>> e) {
		return correlationMethods.add(e);
	}
	
	public RandomExpreriment( String experimentPath) {
		super();
		expPath = experimentPath;
		agreementMethods = new Vector<PartiotioningAgreement<V>>();
		correlationMethods = new Vector<Correlation<Vector<? extends Number>>>();
	}
	public RandomExpreriment() {
		super();
		agreementMethods = new Vector<PartiotioningAgreement<V>>();
		correlationMethods = new Vector<Correlation<Vector<? extends Number>>>();
	}
	public void  logExternalCorrelations(Vector<Vector<Vector<Number>>> extRess){
		//System.err.println(extRess);
//		---------------------------------- How much External Indexes Agree? ----------------------------
		Logger.logFunction("Correlation of external Criteria");
		
		double exCorr[][][] = new double[agreementMethods.size()][agreementMethods.size()][correlationMethods.size()];
		double exVarTmp[][][] = new double[agreementMethods.size()][agreementMethods.size()][correlationMethods.size()];
		for (int d = 0; d< extRess.size(); d++) {
			for (int a1 = 0; a1 < agreementMethods.size(); a1++) {
				for (int a2 = 0; a2 < agreementMethods.size(); a2++) if(a1!=a2){
					for (int cm =0; cm< correlationMethods.size(); cm ++ ){
						double corrScore = computeCorrelation(correlationMethods.get(cm), extRess.get(d).get(a1),extRess.get(d).get(a2)).doubleValue();
						exCorr[a1][a2][cm]+=corrScore;
						exVarTmp[a1][a2][cm]+=corrScore*corrScore;
						Logger.logln("% " + agreementMethods.get(a1)+ " & " + agreementMethods.get(a2) +" agree " + df.format(corrScore) + " on " + correlationMethods.get(cm));
					}
				}
			}
		}
		
		for (int cm =0; cm< correlationMethods.size(); cm ++ ){	
			Logger.log("\\begin{table}\n\\centering\n\\begin{tabular}{|l|", DebugMode.result);
			for (int a1 = 0; a1 < agreementMethods.size(); a1++) Logger.log("l |", DebugMode.result);
			Logger.logln("}\n\\hline", DebugMode.result);
			Logger.log(" Index ", DebugMode.result);
			for (int a1 = 0; a1 < agreementMethods.size(); a1++) Logger.log( " & " + agreementMethods.get(a1) +" " , DebugMode.result);
			Logger.logln(" \\\\ \n\\hline" , DebugMode.result);
	
			for (int a1 = 0; a1 < agreementMethods.size(); a1++) {
				Logger.log(agreementMethods.get(a1)  , DebugMode.result);
				for (int a2 = 0; a2 < agreementMethods.size(); a2++) {
					if(a1!=a2){
						exCorr[a1][a2][cm] /= extRess.size();
						exVarTmp[a1][a2][cm] = Math.sqrt((exVarTmp[a1][a2][cm]/extRess.size() - exCorr[a1][a2][cm]*exCorr[a1][a2][cm]));
					}else {
						exCorr[a1][a2][cm] =1;
						exVarTmp[a1][a2][cm] =0;
					}
					Logger.log(" & " + df.format(exCorr[a1][a2][cm]) + (extRess.size()>1?("$\\pm$" +   df.format(exVarTmp[a1][a2][cm])):"")+"" , DebugMode.result);
//					Logger.log(" & " + exCorr[a1][a2][cm] + (("$\\pm$" +   df.format(exVarTmp[a1][a2][cm])))+"" , DebugMode.result);

				}
				Logger.logln(" \\\\ " , DebugMode.result);
			}
			Logger.logln("\\hline  ", DebugMode.result);
			Logger.logln("\\end{tabular}", DebugMode.result);
			Logger.logln("\\caption{Correlation between external indexes on "+""+new File(expPath).getAbsolutePath()+" dataset, " + " based on " +correlationMethods.get(cm) +"}",DebugMode.result);
			Logger.logln("\\label{table:res1}", DebugMode.result);
			Logger.logln("\\end{table}", DebugMode.result);
		}		
	}
	
	public Vector<PartiotioningAgreement<V>> getExternalAlternatives(Graph<V, E> graph,  Transformer<E, ? extends Number> weights){
		Vector< PartiotioningAgreement<V>> res = new Vector< PartiotioningAgreement<V>>();
		
		MeasuresUtil.getAgreementAlternatives(graph, ( Transformer<E, Double>)weights);
		//res.add( new ARI<V>());//same //	comparer.add( new AGAM <Integer>(GAM.Type.RI));
		//same as ARI//comparer.add( new AGAM <Integer>(gam.getRIPhi(),gam.getIntersectionEtha()));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.Nodes ));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.AdjustedNodes ,AdjustionMethod.NORMALIZE));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.AdjustedNodes ,AdjustionMethod.NORMALIZE));
//
//		//	res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.AdjustedNodes ,false));
//	//	res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Nodes ));
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.Edges ));
//	//	res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Edges ));
////
//		Transformer phi =  new Transformer<Double, Double>() {
//				public Double transform(Double x) {return x==0 ? 0 : Math.log(x);}
//			};
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.ELSE,ExternalOverlap.Nodes,phi ));
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.ELSE,ExternalOverlap.NodesWeightedByDegree,phi ));
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.ELSE,ExternalOverlap.Edges,phi ));
//		
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.NodesWeightedByDegree ));
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.NodesWeightedByTRIANGLES ));
////		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.NodesWeightedByClusteringCoefficient));
//		//new GAM <Integer>(GAM.Type.VI);
////		res.add(new Rand<V>());//res.add(new  GAM <V>(GAM.Type.RI));
////		res.add(new NMI<V>(NMI.Mode.SUM));//res.add(new  NMIAlt <V>());
//		res.add(new NMI<V>(NMI.Mode.SQRT));
//	//	res.add(new VI<V>());
//	//	res.add(new  GAM <V>(GAM.Type.VI));
////	res.add(new  GAM <V>(phi));
////		res.add(new Jaccard<V>());
//		res.add(new AMI<V>()); //res.add(new  AGAM <V>(GAM.Type.VI));
////		res.add(new FMeasure<V>(.5));
////		res.add(new FMeasure<V>(2));
////		res.add(new FMeasure<V>(1));

		return res;
	}
	public void addCorrelationAlternatives(){
		addCorrelation( new SpearmanCorrelation());
		addCorrelation( new PearsonCorrelation());
	}

	public Vector<Vector<Vector<Number>>>  externalEvaluation(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings){
		//----------------------------------External Evaluation of the results----------------------------
		Logger.logFunction("Computing External Agreements " );
		Vector<Vector<Vector<Number>>> extRess = new Vector<Vector<Vector<Number>>>();
		for (int d = 0; d< dataset.size(); d++) {
			System.err.println(dataset.get(d).name+"--------------");
			agreementMethods = getExternalAlternatives(dataset.get(d).graph, dataset.get(d).weights);
			addCorrelationAlternatives();
			
			Vector<Vector<Number>> tmp = new Vector<Vector<Number>>();
			for (int am = 0; am < agreementMethods.size(); am++) {
				Vector<Number> extRes = externalEvaluation(agreementMethods.get(am), dataset.get(d).groundTruth, partitionings.get(d));
				tmp.add(extRes);
			}
			extRess.add(tmp);
		}
		logExternalCorrelations(extRess);
		return extRess;
	}
	
//	public Vector<Vector<Vector<Number>>>  externalEvaluation(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings){
//		//----------------------------------External Evaluation of the results----------------------------
//		Logger.logFunction("Computing External Agreements " );
//		Vector<Vector<Vector<Number>>> extRess = new Vector<Vector<Vector<Number>>>();
//		for (int d = 0; d< dataset.size(); d++) {
//			System.err.println(dataset.get(d).name+"--------------");
//			Vector<Vector<Number>> tmp = new Vector<Vector<Number>>();
//			for (int am = 0; am < agreementMethods.size(); am++) {
//				Vector<Number> extRes = externalEvaluation(agreementMethods.get(am), dataset.get(d).groundTruth, partitionings.get(d));
//				tmp.add(extRes);
//			}
//			extRess.add(tmp);
//		}
//		logExternalCorrelations(extRess);
//		return extRess;
//	}
	// returns k*(k-1) /2 variations
	   @SuppressWarnings("unchecked")
	protected Vector<Set<V>> generateMergedVariations(Vector<Set<V>> groundTruth, int mergeChancePercentage){
		Vector<Set<V>> result = null;
		
		if(groundTruth.size()<=2) return result;
	
		//1 level merge
		for (int i = 0; i < groundTruth.size(); i++) {
			for (int j = i+1 ; j < groundTruth.size(); j++) 
				if(random.nextInt(100) < mergeChancePercentage){
					if(result == null ) result =  (Vector<Set<V>> ) groundTruth.clone();					
					if(result.remove(groundTruth.get(i)) &&	result.remove(groundTruth.get(j))){
						Set<V> unioun= new HashSet<V>();
						unioun.addAll(groundTruth.get(i));
						unioun.addAll(groundTruth.get(j));
						result.add(unioun);
					}
			}
			
		}
		return result;
	}
	
	//returns k variations
	@SuppressWarnings("unchecked")
	protected Vector<Set<V>> generateSplittedVariations(Vector<Set<V>> groundTruth, int splitChancePercentage){
		Vector<Set<V>> result = null;

		//1 level split
		for (Set<V> c1 : groundTruth) if(random.nextInt(100) < splitChancePercentage){
            if(result == null ) result =  (Vector<Set<V>> ) groundTruth.clone();
            result.remove(c1);
			Set<V> splitted1 = new HashSet<V>();
			Set<V> splitted2 = new HashSet<V>();
			
			for (V v : c1) {
				if(random.nextBoolean())
					splitted1.add(v);
				else splitted2.add(v);
			}
			
//			HierarchicalBetweeness<V, E> betweeness = new HierarchicalBetweeness<V, E>();
//			tmp.addAll(betweeness.findCommunities(FilterUtils.createInducedSubgraph(c1,graph), 2));
			if(splitted1.size()>1 && splitted2.size()>1){
				result.add(splitted1);
				result.add(splitted2);
			}else result.add(c1);
		}
		return result;
	}
	protected Vector<Set<V>> moveNodesArround (Vector<Set<V>> vector, int moveChancePercentage){ //0..100
		Vector<Set<V>> result = new Vector<Set<V>>();
		
		for (Set<V> set : vector) {
			result.add(new HashSet<V>(set));
		}
		
		class TMP {
			V v;
			Set<V> from;
			Set<V> to;
			public TMP(V v, Set<V> from, Set<V> to) {
				super();
				this.v = v;
				this.from = from;
				this.to = to;
			}
		}
		
		Vector<TMP> toMove = new Vector<TMP>();
		
		//randomly moving nodes between clusters
		for (Set<V> set : result) {
			for (V v : set) {
				//change it : D
				if(random.nextInt(100) < moveChancePercentage){
					toMove.add(new TMP(v, set, result.get(random.nextInt(result.size()))));
				}
			}
		}
		
		for (TMP tmp : toMove) {
			if(tmp.from.size()>1){ // preventing emergence of empty node clusters
				tmp.from.remove(tmp.v);
				tmp.to.add(tmp.v);
			}
		}
		
		return result;
	}
	
	int maxResult = 25;
	int mergeChance = 30, swapChance = 20 , splitChance = 20;
	
	Random random = new Random();
	protected Vector<Vector<Set<V>>> generateDifferentVariationFromGroundTruth(Graph<V,E> graph, Vector<Set<V>> groundTruth){
		Vector<Vector<Set<V>>> results = new Vector<Vector<Set<V>>>();
		results.add(groundTruth);
		
//		System.err.println("---------------------------------  K: " + groundTruth.size());
//		System.err.println(" generateSplittedVariations");
		 
		maxResult/=2;
		while(results.size() < maxResult) {
			Vector<Vector<Set<V>>> vars = new Vector<Vector<Set<V>>>();
			for (Vector<Set<V>> vector : results) {
				Vector<Set<V>> tmp =  generateSplittedVariations(vector,splitChance);
				if(tmp!=null && tmp.size()>=2 && validate(graph, vector)) vars.add(tmp);
				if(results.size() + vars.size() > maxResult) break;
			}
//			System.err.println("   "+vars.size()+":"+results.size());
			
			for (Vector<Set<V>> vector : results) {
				Vector<Set<V>> tmp =   generateMergedVariations(vector,mergeChance);
				if(tmp!=null && tmp.size()>=2 && validate(graph, vector)) vars.add(tmp);

				if(results.size() + vars.size() > maxResult) break;
			}
//			System.err.println("        "+vars.size()+":"+results.size());

			for (Vector<Set<V>> vector : vars) {
				if(results.size() < maxResult && validate(graph, vector))
					results.add( vector);
				else break;
			}
		}
//		System.err.println(" Randomize");

		maxResult*=2;
		while (results.size() < maxResult) {
			Vector<Vector<Set<V>>> vars = new Vector<Vector<Set<V>>>();
			for (Vector<Set<V>> vector : results){
				 if(validate(graph, vector)) vars.add(moveNodesArround(vector, swapChance));
				if (results.size() + vars.size() > maxResult)
					break;
			}

			for (Vector<Set<V>> vector : vars) {
				if (results.size() < maxResult)
					 if(validate(graph, vector)) results.add(vector);
				else
					break;
			}
		}
//	
//		while(results.size() < maxRandomizeResult) {
//			Vector<Vector<Set<V>>> vars = new Vector<Vector<Set<V>>>();
//			for (Vector<Set<V>> vector : results) {
//		
//			}
////			System.err.println( " ransomized into :" + vars.size()+" variations, ");
//			
//			for (Vector<Set<V>> vector : vars) {
//				if(results.size() < maxRandomizeResult)
//					results.add( vector);
//				else break;
//			}
//		}
//		System.err.println();
	//	System.err.println("--------->" + results.size());
	
		return results;
	}
	//Some Statistics about the current data set
	protected void logDatasetStat(Vector<NetworkWithGroundTruth> dataset, Vector<Vector<Vector<Set<V>>>> partitionings , Vector<Vector<Vector<Number>>> extRess, int am){
//		for (int am = 0; am < agreementMethods.size(); am++) {
				Logger.logFunction(agreementMethods.get(am).toString() , DebugMode.result);
				Logger.logln("\\begin{table}", DebugMode.result);
				Logger.logln("\\begin{tabular}{| l |l | l| l|l|l|l|}", DebugMode.result);
				Logger.logln("\\hline" , DebugMode.result);
				  
				Logger.log("Dataset & $K^*$ & $\\#$  & $\\overline{K}$  & $\\overline{"+agreementMethods.get(am)+"}$ \\\\ "  , DebugMode.result);
				Logger.logln("\\hline" , DebugMode.result);
				
				for (int d = 0; d< partitionings.size(); d++) {
					double avgAM=0, varAM=0, minAM=10000000, maxAM=0, eind , k , avgK = 0, varK = 0, minK = 10000000, maxK = 0, N = partitionings.get(d).size();
					for (int i = 0; i< partitionings.get(d).size(); i++){
						eind = extRess.get(d).get(am).get(i).doubleValue();
						k = partitionings.get(d).get(i).size();
						
						avgK += k;
						varK += k * k;
						if (k < minK) minK = k;
						if (k > maxK) maxK = k;
						
						avgAM += eind;
						varAM += eind* eind;
						if(eind<minAM) minAM =eind;
						if(eind>maxAM) maxAM =eind;
					}
					avgK /= N;
					varK /= N;
					varK -= avgK * avgK;
					varK = Math.sqrt(varK);
					
					avgAM /= N;
					varAM /=N;
					varAM -= avgAM*avgAM;
					varAM = Math.sqrt(varAM);
					
//					football &  12 & 60 & 10.1$\pm$4.90$\in$[3,19]&  0.73$\pm$0.14$\in$[0.38,1] \\
					Logger.log(dataset.get(d).name + " & "	+ dataset.get(d).groundTruth.size() + " & ", DebugMode.result);
					Logger.log( partitionings.get(d).size() + " & ", DebugMode.result);
					Logger.log(df2.format(avgK) + "$\\pm$" + df2.format(varK) + "$\\in$"+ "[" +  df2.format(minK) + "," +  df2.format(maxK) + "]"+ " & ", DebugMode.result);
					Logger.log(df2.format(avgAM) + "$\\pm$" + df2.format(varAM) + "$\\in$"+ "[" +  df2.format(minAM) + "," +  df2.format(maxAM) + "]", DebugMode.result);

					Logger.logln(" \\\\ " , DebugMode.result);
				}
//				for (int j = 0; j < winners.get(am).get(cm).size() && j<TOPWINNERS; j++) {
//					Logger.log((j+1) + " & " , DebugMode.result);
//					Logger.log(winners.get(am).get(cm).get(j) , DebugMode.result);
//				}
				Logger.logln("\\hline  ", DebugMode.result);
				Logger.logln("\\end{tabular}", DebugMode.result);
				Logger.logln("\\caption{"+new File(expPath).getPath()+" Dataset with sampling parameter: "+"s="+splitChance+"\\%, m="+ mergeChance +"\\%, r="+swapChance+ "\\%, "   +"}",DebugMode.result);
				Logger.logln("\\label{table:res1}", DebugMode.result);
				Logger.logln("\\end{table}", DebugMode.result);
				
//			}
		
	}
	protected Vector<Vector<Set<V>>> generateDifferentKmeans(Graph<V,E> graph, int maxK){
		Vector<Vector<Set<V>>> vector = new Vector<Vector<Set<V>>>();
		for (int i = 2; i < maxK; i++) {
			vector.add((new TopLeaders<V, E>(i)).findCommunities(graph, null).getGroups());
		}
		
		return vector;
	}
	
//	protected Vector<Vector<Set<V>>> generateDifferentHierarical(Graph<V,E> graph, int maxK){
//		//Vector<Vector<Set<V>>> vector = new Vector<Vector<Set<V>>>();
//		//TODO:  why this much bad results? debug the hierarical or change it to dist based sth
//	//	HierarchicalBetweeness<V, E> hierarchicalBetweeness = new HierarchicalBetweeness<V, E>();
//		
//		return hierarchicalBetweeness.findCommunities(graph, 2, maxK);
//	}
	
	protected Vector<Vector<Set<V>>> generateDifferentCommunityMiningMethods(Graph<V,E> graph){
		Vector<CommunityMiner<V, E>> communityMiners =  new Vector<CommunityMiner<V,E>>();
		
//		communityMiners.add(new CliquePercolation<V, E>());
		communityMiners.add(new FastModularity<V, E>());
//		communityMiners.add(new MaxMinModularity<V, E>());
//		communityMiners.add(new DensityMethod<V, E>());
//		communityMiners.add(new Local<V, E>());
		//communityMiners.add(new LocalMetric<V, E>(criteria));
//		communityMiners.add(new SCAN<V, E>());
		communityMiners.add(new Infomap<V, E>());
		//communityMiners.add(new Infomod<V, E>());
		
		Vector<Vector<Set<V>>> vector = new Vector<Vector<Set<V>>>();
		for (CommunityMiner<V, E> communityMiner : communityMiners) {
			Logger.logln(communityMiner);
			vector.add(communityMiner.findCommunities(graph).getGroups());
		}
		
		return vector;
	}
	
	
	boolean addCMs =false, addkmeans=false, addHierar=false, addGTalt = true;

	/**
	 * @param graph
	 * @return a collection of possible partitioning/results for the given data set 
	 */
	protected Vector<Vector<Set<V>>> generateDifferentPartitionings(Graph<V,E> graph , Vector<Set<V>> groundTruth){
		Vector<Vector<Set<V>>> res = new Vector<Vector<Set<V>>>();
		
		int maxK = (int) Math.sqrt(graph.getVertexCount()) ;
		
		Vector<Vector<Set<V>>> parts;
		if(addCMs){
			parts = generateDifferentCommunityMiningMethods(graph);
			for (Vector<Set<V>> vector : parts) {
				if(validate(graph, vector)) {
					res.add(vector);
					Logger.logln(vector,DebugMode.detailed);
				}
			}
		}
		if(addkmeans){
			parts =  generateDifferentKmeans(graph, 10);
			for (Vector<Set<V>> vector : parts) {
				if(validate(graph, vector)) {
					res.add(vector);
					Logger.logln(vector,DebugMode.detailed);
				}
			}
		}
//		if(addHierar){
//			parts = generateDifferentHierarical(graph, maxK);
//			for (Vector<Set<V>> vector : parts) {
//				if(validate(graph, vector)) {
//					res.add(vector);
//					Logger.logln(vector,DebugMode.detailed);//System.err.println(vector);
//				}
//			}
//		}
		if(addGTalt){
		//	System.err.println(" to call generateDifferentVariationFromGroundTruth");
			parts = generateDifferentVariationFromGroundTruth(graph, groundTruth);
			for (Vector<Set<V>> vector : parts) {
				if(validate(graph, vector)) {
					res.add(vector);
					Logger.logln(vector,DebugMode.detailed);//System.err.println(vector);
				}
			}
		}
		Logger.logFunction("DatasetGenerationCompleted");

		return res;
			
	}
	
	protected Vector<Vector<Set<V>>> generateRandomClusters(Collection<V> nodes, int minK, int maxK, int nSamples){
		Vector<Vector<Set<V>>> res = new Vector<Vector<Set<V>>>();
		
		for (int k = minK; k < maxK; k++) {
			for (int sample = 0; sample < nSamples;) {
				Vector<V> nodesRemain = new Vector<V>(nodes);
//				System.err.println(nodes.size() + " " + maxK);
				
				Vector<Set<V>> randomCluster = new Vector<Set<V>>();
				for (int c = 0; c < k; c++) {
					randomCluster.add(new HashSet<V>());
					//we don't want an empty cluster
					int ind = (nodesRemain.size()>0)?random.nextInt(nodesRemain.size()):0;
					randomCluster.lastElement().add(nodesRemain.get(ind));
					nodesRemain.remove(ind);
				}
//				System.err.println("---- "+nodesRemain.size() + " " + k);
			
				for(V v: nodesRemain){// (int i = 0; i < nodes.size(); i++) {
					randomCluster.get(random.nextInt(k)).add(v);//nodes.get(i));
				}

				if (res.indexOf(randomCluster) == -1) {
					res.add(randomCluster);
					 sample++;
				}
			}
		}
//		System.err.println(res);
		Logger.logFunction("DatasetGenerationCompleted");

		return res;
			
	}
	
	

}
