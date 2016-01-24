package dev_Experiments;

import io.Logger;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import measure.base.*;
import measure.cluster.agreement.*;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.criteria.*;
import measure.criteria.RelativeClusterCriteria.Between;
import measure.criteria.RelativeClusterCriteria.Within;
import measure.graph.criteria.*;
import measure.graph.distance.*;
import measure.graph.distance.ModSimilarity.Norm;
import edu.uci.ics.jung.graph.Graph;
import io.Logger.DebugMode;

public class CriteriaComparer<V,E> extends RandomExpreriment<V, E>{
	public CriteriaComparer(String experimentPath) {
		super( experimentPath);
	}
	
	//TODO: add other criteria
	public Vector<RelativeCriteria<V>> getCriteiaAlternatives(Graph<V, E> graph,  Transformer<E, ? extends Number> weights){
		Vector<RelativeCriteria<V>> criterias =  new Vector<RelativeCriteria<V>>();

		// Adding Graph based Criteria
		criterias.add(new Modularity<V, E>(graph));
//		criterias.add(new AltModularity<V, E>(graph));
//		criterias.add(new M1Modularity<V, E>(graph));
//		//criterias.add(new M1AltModularity<V, E>(graph));
//		criterias.add(new M2Modularity<V, E>(graph));
//
//		criterias.add(new M3Modularity<V, E>(graph));
//		criterias.add(new M31Modularity<V, E>(graph));
//		criterias.add(new M32Modularity<V, E>(graph));
//		criterias.add(new M33Modularity<V, E>(graph));
//		criterias.add(new M34Modularity<V, E>(graph));
//
//		criterias.add(new M4Modularity<V, E>(graph));
//		criterias.add(new M5Modularity<V, E>(graph));
//		criterias.add(new M6Modularity<V, E>(graph));
//		criterias.add(new Triularity<V, E>(graph));
////		
		// Adding ClusterBased Criteria 
		Vector<Proximity<V>> distances = getGraphDistanceAlternatives(graph, weights);
//
		for (Proximity<V> distanceMethod: distances ) {
			Centroid<V> centroidMethod = new Medoid<V>(distanceMethod);
				
			criterias.add(new VarianceRatio<V>(false).setMetrics(distanceMethod, centroidMethod));
			if(distanceMethod.isSimilarity())
				criterias.add(new VarianceRatio<V>(true).setMetrics(distanceMethod, centroidMethod));
			
			if(distanceMethod.isSimilarity())
				criterias.add(new DaviesBouldin<V>(true).setMetrics(distanceMethod, centroidMethod));
			criterias.add(new DaviesBouldin<V>(false).setMetrics(distanceMethod, centroidMethod));
//////				
//				// Different variation of Dunn
			if(!distanceMethod.isSimilarity())
			for (Between between : Between.values()) {
				for (Within within : Within.values()) {
					if(between!=Between.SUM && within!=Within.SUM &&  within!=Within.SUM2C ){ // ignoring unnecessary ones
//						//if(!distanceMethod.isSimilarity())
						criterias.add(new Dunn<V>(between, within,false).setMetrics(distanceMethod, centroidMethod));
//						//if(distanceMethod.isSimilarity()) //Dunn is not sim friendly!!!
//							//criterias.add(new Dunn<V>(between, within,true).setMetrics(distanceMethod, centroidMethod));
					}
				}
			}
//			
			if(distanceMethod.isSimilarity()) 		
				criterias.add(new PBM<V>(true).setMetrics(distanceMethod, centroidMethod));
			criterias.add(new PBM<V>(false).setMetrics(distanceMethod, centroidMethod));
////				
////////	//Different Variation of Silhouette
			criterias.add(new Silhouette<V>(Between.AVERAGE,false,false).setMetrics(distanceMethod, centroidMethod)); //SWC
			criterias.add(new Silhouette<V>(Between.CENTROID,false,false).setMetrics(distanceMethod, centroidMethod)); //SSWC
			if(distanceMethod.isSimilarity()){
				criterias.add(new Silhouette<V>(Between.AVERAGE,false,true).setMetrics(distanceMethod, centroidMethod)); //SWC
				criterias.add(new Silhouette<V>(Between.CENTROID,false,true).setMetrics(distanceMethod, centroidMethod)); //SSWC
			}
			criterias.add(new Silhouette<V>(Between.AVERAGE,true,false).setMetrics(distanceMethod, centroidMethod));//ASWC
			criterias.add(new Silhouette<V>(Between.CENTROID,true,false).setMetrics(distanceMethod, centroidMethod));//ASSWC
			if(distanceMethod.isSimilarity()){
				criterias.add(new Silhouette<V>(Between.AVERAGE,true,true).setMetrics(distanceMethod, centroidMethod));//ASWC
				criterias.add(new Silhouette<V>(Between.CENTROID,true,true).setMetrics(distanceMethod, centroidMethod));//ASSWC
			}
//			
			if(distanceMethod.isSimilarity()) 					
				criterias.add(new PointBiserial<V>(true).setMetrics(distanceMethod, centroidMethod));
			else 
				criterias.add(new PointBiserial<V>(false).setMetrics(distanceMethod, centroidMethod));
			
			if(distanceMethod.isSimilarity())
				criterias.add(new CIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
			else 
				criterias.add(new CIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
//			
			if(distanceMethod.isSimilarity()) 
				criterias.add(new ZIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
			else
				criterias.add(new ZIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
			
//			if(distanceMethod.isSimilarity()) //This is similar to Zindex no need to include
//				criterias.add(new ZIndexGraph<V,E>(graph,true).setMetrics(distanceMethod, centroidMethod));
		//	criterias.add(new ZIndexGraph<V,E>(graph,false).setMetrics(distanceMethod, centroidMethod));
				
		}
		
		return criterias;
	}
	
	/**
	 * @param graph
	 * @param weights
	 * @return possible alternative distance measures that could compute distance between two nodes in the given graph
	 */
	public Vector<Proximity<V>>  getGraphDistanceAlternatives(Graph<V, E> graph,  Transformer<E, ? extends Number> weights){
		Vector<Proximity<V>> distances = new Vector<Proximity<V>>();
		//TODO: 1) getDist, use it instead of division for sim to dis conversion. 2) Icloseness check 
		distances.add(new Adjacency<V, E>(graph ,weights,false));
		distances.add(new Adjacency<V, E>(graph ,weights,true));

		distances.add(new AdjacencyRelation<V, E>(graph ,weights,false));
		distances.add(new AdjacencyRelation<V, E>(graph ,weights,true));
		
////		distances.add(new NeighborOverlapOld<V, E>(graph ,weights));		//Not documented
////		distances.add(new Overlap<V, E>(graph ,weights,false)); //unnormalized 		//Not documented
////		distances.add(new Overlap<V, E>(graph ,weights,true)); //unnormalized 		//Not documented
		distances.add(new NeighborOverlap<V, E>(graph ,weights,false));
		distances.add(new NeighborOverlap<V, E>(graph ,weights,true));
		distances.add(new NeighborOverlapVariation<V, E>(graph ,weights,false));		//Not documented
		distances.add(new NeighborOverlapVariation<V, E>(graph ,weights,true));		//Same as overlap; Not documented
		distances.add(new TopologicalOverlap<V, E>(graph ,weights,false));
		distances.add(new TopologicalOverlap<V, E>(graph ,weights,true));
//
		distances.add(new AdjacencyPearsonCorrelation<V, E>(graph ,weights,false,false));
		distances.add(new AdjacencyPearsonCorrelation<V, E>(graph ,weights,true,false)); //Normalized
		distances.add(new AdjacencyPearsonCorrelation<V, E>(graph ,weights,true,true)); //Normalized Augmented
		distances.add(new AdjacencyPearsonOverlap<V, E>(graph ,weights)); //Normalized Augmented

		distances.add(new ShortestPath<V, E>(graph ,weights));

		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.UNI,2));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.UNI,3));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.LIN,2));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.LIN,3));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.EXP,2));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.EXP,3));
		
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.UNI,2,true));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.UNI,3,true));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.LIN,2,true));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.LIN,3,true));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.EXP,2,true));
		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.EXP,3,true));
////		
		distances.add(new ModSimilarity<V, E>(graph ,weights,false,Norm.DIVIDE));		
		distances.add(new ModSimilarity<V, E>(graph ,weights,true,Norm.DIVIDE));		
		distances.add(new ModSimilarity<V, E>(graph ,weights,false,Norm.MINUS));
		distances.add(new ModSimilarity<V, E>(graph ,weights,true,Norm.MINUS));	
		distances.add(new ICloseness<V, E>(graph ,weights)); 
		distances.add(new ICloseness<V, E>(3, graph ,weights)); 
		distances.add(new ICloseness<V, E>(1, graph ,weights)); 
		distances.add(new IClosenessVariation<V, E>(graph ,weights)); 
		distances.add(new IClosenessVariation<V, E>(3, graph ,weights)); 
		distances.add(new IClosenessVariation<V, E>(1, graph ,weights)); 
//		
		System.err.println("We have "+distances.size()+" distance variations ");
		return distances;
	}
	
	protected class Winner{
		RelativeCriteria<V> c;
		DecimalFormat df = new DecimalFormat("#.###");
		double avgScore,varScore;
		
		public Winner(RelativeCriteria<V> c, double avg, double var) {
			super();
			this.c = c;
			this.avgScore = avg;
			this.varScore = var;
		}
		
		public String toString(){
			return c.getName() + " & " + df.format(avgScore) + "$\\pm$" + df.format(varScore);
		}
		
		@SuppressWarnings("unchecked")
		public boolean equals(Object o){
			if(! (o.getClass().equals(this.getClass()))) return false;
			else return ((Winner)o).c.equals(c);
		}
	}
	
	int TOPWINNERS = 50;
	public Vector<Winner> rank(Vector<RelativeCriteria<V>> criterias , Vector<double[]> correlations){ 
		//TODO: also consider statistical significance?
		
		final double[] avgs = new double[criterias.size()];
		final double[] vars = new double[criterias.size()];
		
		for (int i = 0; i < criterias.size(); i++) {
	//		Logger.log(criterias.get(i) + " : ");
			for (int j = 0; j < correlations.get(i).length; j++) {
				avgs[i] += correlations.get(i)[j];
				vars[i] += correlations.get(i)[j] * correlations.get(i)[j] ;
			
	//			Logger.log(correlations.get(i)[j] + " , ");
			}
			avgs[i] /= correlations.get(i).length;
			vars[i] /= correlations.get(i).length;
			vars[i] -= avgs[i]*avgs[i];
			vars[i] = Math.sqrt(vars[i]);
			
	//		Logger.logln("");
		}
		
		
		Vector<Integer> ranks = new Vector<Integer>();
		for (int i = 0; i < criterias.size(); i++) {
			ranks.add(i);
		}
		
		Collections.sort(ranks, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				return new Double(avgs[o2]).compareTo(avgs[o1]) ;
			}
		});
		
		Vector<Winner> winners= new Vector<Winner>();
//		Logger.logFunction("########  AVG Results  #########");
		for (int i = 0; i < ranks.size(); i++) { //TOP TEN
//			Logger.log(criterias.get(ranks.get(i)) + " :: " + avgs[ranks.get(i)] + " +- " + vars[ranks.get(i)] + " :::: [ ");
//			for (int j = 0; j < correlations.get(ranks.get(i)).length; j++) {
//				Logger.log(correlations.get(ranks.get(i))[j] + " , ");
//			}
//			Logger.logln("]");
			winners.add(new Winner(criterias.get(ranks.get(i)), avgs[ranks.get(i)] , vars[ranks.get(i)]));
		}
		return winners;
	
	}


	Transformer<V, Double> nodeWeights;

	public void compareCriteria(){
		//Preparing Data set
		String  DATASET = "/data/";
		Vector<NetworkWithGroundTruth> dataset =  loadNetworks(expPath + DATASET);
		Vector<Vector<Vector<Set<V>>>> partitionings = new Vector<Vector<Vector<Set<V>>>>(); // d x r x Vector<Set<V>>
		for (NetworkWithGroundTruth networkWithGroundTruth : dataset) {
			final Graph<V,E> curGraph = networkWithGroundTruth.graph;
			
			nodeWeights = new Transformer<V, Double>() {
				public Double transform(V arg0) {
					double res =0;
					for (V v1: curGraph.getNeighbors(arg0)) 
						for (V v2: curGraph.getNeighbors(arg0)) 
							if(curGraph.getNeighbors(v1).contains(v2)){
								res+=1;
							}
//					res/=2;
					//res/=(curGraph.degree(arg0)*(curGraph.degree(arg0)-1));
					return res;//new Double(curGraph.degree(arg0));
				}
			};
			partitionings.add(generateDifferentPartitionings(networkWithGroundTruth.graph, networkWithGroundTruth.groundTruth));
			System.err.println(partitionings.lastElement().size() + " partitions generated for "+ networkWithGroundTruth.name );
		}
		Vector<Vector<Vector<Number>>> extRess = externalEvaluation( dataset , partitionings); //d x am x p
		for (int rm = 0; rm < agreementMethods.size(); rm++) {
			logDatasetStat(dataset, partitionings , extRess,rm);
		}
		Vector<Vector<Vector<Number>>> criteriaResultss = internalEvaluation( dataset , partitionings); //d x cm x p
		plotAllCorrelationFigures(dataset,partitionings,extRess,criteriaResultss);		
		
		//overall
	//	System.err.println(criterias.size());
	
		//divide partitioings of each dataset based on this agreementMethods.get(am)			
		for (int rm = 0; rm < agreementMethods.size(); rm++) {
			Vector<Vector<Vector<Vector<Winner>>>> winners = new Vector<Vector<Vector<Vector<Winner>>>>();
			
			//l x d x am x p
			Vector<Vector<Vector<Vector<Number>>>> leveledExtRes = new Vector<Vector<Vector<Vector<Number>>>>();
			//l x d x cm x p
			Vector<Vector<Vector<Vector<Number>>>> leveledCriteriaResults = new Vector<Vector<Vector<Vector<Number>>>>();
			//l x d x r x Vector<Set<V>>
			Vector<Vector<Vector<Vector<Set<V>>>>> leveledPartitionings = new Vector<Vector<Vector<Vector<Set<V>>>>> ();
			
			for (int l = 0; l < 3; l++) {//three level 
				leveledExtRes.add(new Vector<Vector<Vector<Number>>>());
				leveledCriteriaResults.add(new Vector<Vector<Vector<Number>>>());
				leveledPartitionings.add(new Vector<Vector<Vector<Set<V>>>>());
			}
			for (int d = 0; d < dataset.size(); d++) {
				for (int l = 0; l < 3; l++) {//each a vector of datasets
					leveledCriteriaResults.get(l).add(new Vector<Vector<Number>>());
					leveledExtRes.get(l).add(new Vector<Vector<Number>>());
					leveledPartitionings.get(l).add(new Vector<Vector<Set<V>>>());
				}	
				
				for (int l = 0; l < 3; l++) {
					for (int am = 0; am < agreementMethods.size(); am++) {
						leveledExtRes.get(l).get(d).add(new Vector<Number>());//each cm
					}
					for (int crit = 0; crit < criterias.size(); crit++) {
						leveledCriteriaResults.get(l).get(d).add(new Vector<Number>());//each am
					}
				}
				
				double i1 = .3 ,i2 = .6;
				//TODO: should be divided to 3 equal part? set i1 i2 respectively   
				double b = 10000,e = -10000 ;
				for (int i = 0; i < extRess.get(d).get(rm).size(); i++) {
					double ei = extRess.get(d).get(rm).get(i).doubleValue();
					if(ei<b) b = ei;
					if(ei>e) e = ei;
				}
				i1 = (e-b)/3 + b;
				i2 = (e-b)*2/3 + b;
				
				for (int i = 0; i < extRess.get(d).get(rm).size(); i++) {
					double ei = extRess.get(d).get(rm).get(i).doubleValue();
					int level =0;
					if(ei<i1){
						level = 0;
					}else if(ei < i2){
						level = 1;
					}else {
						level = 2;
					}
					
					for (int am = 0; am < agreementMethods.size(); am++) {
						leveledExtRes.get(level).get(d).get(am).add(extRess.get(d).get(am).get(i));
					}
					for (int crit = 0; crit < criterias.size(); crit++) {
						leveledCriteriaResults.get(level).get(d).get(crit).add(criteriaResultss.get(d).get(crit).get(i));
					}
//					System.err.println( .." : "+leveledPartitionings.get(level).get(d).size() + " "+ am + " "+partitionings.get(d).size());
					leveledPartitionings.get(level).get(d).add(partitionings.get(d).get(i));
				}
				
			}
			
			winners.add(compareCorrelationsAndRank(dataset, partitionings, extRess, criteriaResultss));
			//far far results
			winners.add(compareCorrelationsAndRank(dataset, leveledPartitionings.get(0), leveledExtRes.get(0), leveledCriteriaResults.get(0)));
			//medium far results
			winners.add(compareCorrelationsAndRank(dataset, leveledPartitionings.get(1), leveledExtRes.get(1), leveledCriteriaResults.get(1)));
			//near optimal results: // d x r x Vector<Set<V>>
			winners.add(compareCorrelationsAndRank(dataset, leveledPartitionings.get(2), leveledExtRes.get(2), leveledCriteriaResults.get(2)));
		
			Logger.logFunction("");
			Logger.logFunction("");
			Logger.logFunction("Difficulty analysis, sampled partitioned by + " +agreementMethods.get(rm) );
			
			Logger.logFunction("Overall Dataset Statistics");
			logDatasetStat(dataset, partitionings , extRess,rm);
			logWinners(winners.get(0),"Overall Results",rm);
			
			Logger.logFunction("Near Optimal Dataset Statistics");
			logDatasetStat(dataset, leveledPartitionings.get(2) , leveledExtRes.get(2),rm);
			logWinners(winners.get(3),"Near Optimal Results",rm);
		
			Logger.logFunction("Medium Far Dataset Statistics");
			logDatasetStat(dataset, leveledPartitionings.get(1) , leveledExtRes.get(1),rm);
			logWinners(winners.get(2),"Medium Far Results" ,rm);
			
			Logger.logFunction("Far Far Dataset Statistics");
			logDatasetStat(dataset, leveledPartitionings.get(0) , leveledExtRes.get(0),rm);
			logWinners(winners.get(1),"Far Far Results",rm);
			
		}
		
	}

//	public Vector<Vector<Vector<Number>>>  externalEvaluation(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings){
//		//----------------------------------External Evaluation of the results----------------------------
//		Logger.logFunction("Computing External Agreements " );
//		Vector<Vector<Vector<Number>>> extRess = new Vector<Vector<Vector<Number>>>();
//		for (int d = 0; d< dataset.size(); d++) {
//			Vector<Vector<Number>> tmp = new Vector<Vector<Number>>();
//			for (int am = 0; am < agreementMethods.size(); am++) {
//				Vector<Number> extRes = externalEvaluation(agreementMethods.get(am), dataset.get(d).groundTruth, partitionings.get(d));
//				tmp.add(extRes);
//			}
//			extRess.add(tmp);
//		}
////		---------------------------------- How much External Indexes Agree? ----------------------------
//		Logger.logFunction("Correlation of external Criteria");
//		
//		double exCorr[][][] = new double[agreementMethods.size()][agreementMethods.size()][correlationMethods.size()];
//		double exVarTmp[][][] = new double[agreementMethods.size()][agreementMethods.size()][correlationMethods.size()];
//		for (int d = 0; d< dataset.size(); d++) {
//			for (int a1 = 0; a1 < agreementMethods.size(); a1++) {
//				for (int a2 = 0; a2 < agreementMethods.size(); a2++) if(a1!=a2){
//					for (int cm =0; cm< correlationMethods.size(); cm ++ ){
//						double corrScore = computeCorrelation(correlationMethods.get(cm), extRess.get(d).get(a1),extRess.get(d).get(a2)).doubleValue();
//						exCorr[a1][a2][cm]+=corrScore;
//						exVarTmp[a1][a2][cm]+=corrScore*corrScore;
//						Logger.logln("% " + agreementMethods.get(a1)+ " & " + agreementMethods.get(a2) +" agree " + df.format(corrScore) + " on " + correlationMethods.get(cm));
//					}
//				}
//			}
//		}
//		
//		for (int cm =0; cm< correlationMethods.size(); cm ++ ){	
//			Logger.log("\\begin{table}\n\\centering\n\\begin{tabular}{|l|", DebugMode.result);
//			for (int a1 = 0; a1 < agreementMethods.size(); a1++) Logger.log("l |", DebugMode.result);
//			Logger.logln("}\n\\hline", DebugMode.result);
//			Logger.log(" Index ", DebugMode.result);
//			for (int a1 = 0; a1 < agreementMethods.size(); a1++) Logger.log( " & " + agreementMethods.get(a1) +" " , DebugMode.result);
//			Logger.logln(" \\\\ \n\\hline" , DebugMode.result);
//	
//			for (int a1 = 0; a1 < agreementMethods.size(); a1++) {
//				Logger.log(agreementMethods.get(a1)  , DebugMode.result);
//				for (int a2 = 0; a2 < agreementMethods.size(); a2++) {
//					if(a1!=a2){
//						exCorr[a1][a2][cm] /= dataset.size();
//						exVarTmp[a1][a2][cm] = Math.sqrt((exVarTmp[a1][a2][cm]/dataset.size() - exCorr[a1][a2][cm]*exCorr[a1][a2][cm]));
//					}else {
//						exCorr[a1][a2][cm] =1;
//						exVarTmp[a1][a2][cm] =0;
//					}
//					Logger.log(" & " + df.format(exCorr[a1][a2][cm]) + (dataset.size()>1?("$\\pm$" +   df.format(exVarTmp[a1][a2][cm])):"")+"" , DebugMode.result);
//				}
//				Logger.logln(" \\\\ " , DebugMode.result);
//			}
//			Logger.logln("\\hline  ", DebugMode.result);
//			Logger.logln("\\end{tabular}", DebugMode.result);
//			Logger.logln("\\caption{Correlation between external indexes on "+""+new File(expPath).getAbsolutePath()+" dataset, " + " based on " +correlationMethods.get(cm) +"}",DebugMode.result);
//			Logger.logln("\\label{table:res1}", DebugMode.result);
//			Logger.logln("\\end{table}", DebugMode.result);
//		}		
//		return extRess;
//	}
//	
	Vector<RelativeCriteria<V>> criterias = null; 
	public Vector<Vector<Vector<Number>>>  internalEvaluation(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings){
		Vector<Vector<Vector<Number>>> criteriaResultss = new Vector<Vector<Vector<Number>>>();
		for (int d = 0; d< dataset.size(); d++) {
			Logger.logln( "Dataset " +  dataset.get(d).name , DebugMode.result);
	//		logDatasetStat(partitionings.get(d));
		
			criterias = getCriteiaAlternatives(dataset.get(d).graph, null);
			Vector<Vector<Number>> criteriaResultD = new Vector<Vector<Number>>(); //criteria x partitioning
			Logger.logFunction("Computing Criteria Values for " + criterias.size()  + " criterias");
			System.err.println("Computing Criteria Values for " + criterias.size()  + " criterias");
			//For each criterion
			for(int i=0; i< criterias.size();i++){ 
//				Logger.log(criterias.get(i).getName()+", ");
				if(criterias.size()>10 && (i% (criterias.size()/10) )==0)
				Logger.log("("+((i*100)/criterias.size())+"% of "+(d*100/dataset.size())+")");

				Logger.log(".");
				
				//Internal evaluation: Compute the score of all the partitioning using relative criterion indexes  
				Vector<Number> intRes = internalEvaluation(criterias.get(i), partitionings.get(d));
				criteriaResultD.add(intRes);
				
				boolean invalid = false;
				for (Number number : intRes) if( ( number+"").equals("NaN") || ( number+"").equals("Infinity") ) invalid = true;
				if (invalid) {
						Logger.logln("\n Error in computing:  " + criterias.get(i) ,DebugMode.result);//dataset.get(d)  + 
						Logger.logln( intRes ,DebugMode.result);
				}
			}
			criteriaResultss.add(criteriaResultD);
			Logger.logln(" ");
		}
		return criteriaResultss;

	}
	
	void plotAllCorrelationFigures(Vector<NetworkWithGroundTruth> dataset,  Vector<Vector<Vector<Set<V>>>> partitionings, Vector<Vector<Vector<Number>>> extRess, Vector<Vector<Vector<Number>>> criteriaResultss ){
		if(!plot) return;
		for (int d = 0; d< dataset.size(); d++) {
			for (int am = 0; am < agreementMethods.size(); am++) {
				for (int crit = 0; crit < criteriaResultss.get(d).size(); crit++) {
					plot(dataset.get(d), agreementMethods.get(am), criterias.get(crit), extRess.get(d).get(am), criteriaResultss.get(d).get(crit), partitionings.get(d) );
				}
			}
		}
	}
	public Vector<Vector<Vector<Winner>>> compareCorrelationsAndRank(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings, Vector<Vector<Vector<Number>>> extRess, Vector<Vector<Vector<Number>>> criteriaResultss){
		//----------------------------------Internal Evaluation of the results----------------------------
//	System.err.println(dataset.size());
//	System.err.println(partitionings.size());
//	System.err.println(extRess.size());
//	System.err.println(criteriaResultss.size());
	
	for (int d = 0; d< dataset.size(); d++) {
//		System.err.println(dataset.get(d).name);
//		System.err.println(partitionings.get(d).size());
//		System.err.println(extRess.get(d).size());
//		System.err.println(criteriaResultss.get(d).size());
		
//		for (int am = 0; am < agreementMethods.size(); am++) {
//			System.err.println("a+ " +extRess.get(d).get(am).size());
//		}
//		for (int cm =0; cm< criterias.size(); cm ++ ){
//			System.err.println("c+ "+ criteriaResultss.get(d).get(cm).size());
//		}
	}
	
	// Vector Initialization! 
		Vector<Vector<Vector<double[]>>> correlations = new Vector<Vector<Vector<double[]>>>(); // agreementMethod x correlationmethod x criteria x dataset
		for (int am = 0; am < agreementMethods.size(); am++) {
			correlations.add(new Vector<Vector<double[]>>());
			for (int cm =0; cm< correlationMethods.size(); cm ++ ){
				correlations.get(am).add(new Vector<double[]>());
			}
		}
		//-----------------------------------  Computing correlation of internal and external ------------------------------------------
		for (int d = 0; d< dataset.size(); d++) {
			Logger.logFunction("Computing Correlations");
			for (int am = 0; am < agreementMethods.size(); am++) {
				for (int cm =0; cm< correlationMethods.size(); cm ++ ){
					if(correlations.get(am).get(cm).size() == 0){ //initializing vector!
						for (int crit = 0; crit < criteriaResultss.get(d).size(); crit++) {
							correlations.get(am).get(cm).add(new double[dataset.size()]);
						}
					}
					for (int crit = 0; crit < criteriaResultss.get(d).size(); crit++) {
						double corrScore = computeCorrelation(correlationMethods.get(cm),  extRess.get(d).get(am), criteriaResultss.get(d).get(crit)).doubleValue();
						correlations.get(am).get(cm).get(crit)[d] = corrScore;
					}
				}
			}

		}		
		
		//-----------------------------------  Ranking criteria based on the correlations ------------------------------------------
		Logger.logFunction("Ranking");
		Vector<Vector<Vector<Winner>>> winners = new Vector<Vector<Vector<Winner>>> ();
		//Compare correlations of different criterion and rank them
		for (int am = 0; am < agreementMethods.size(); am++) {
			winners.add(new Vector<Vector<Winner>>());
			for (int cm =0; cm< correlationMethods.size(); cm ++ ){
				Logger.logFunction(agreementMethods.get(am));
				Logger.logFunction(correlationMethods.get(cm));
				
				winners.get(am).add(rank(criterias, correlations.get(am).get(cm) ));
			}
		}
		return winners;
	}
	
	public void logWinners(Vector<Vector<Vector<Winner>>> winners, String caption, int agreementMethod){
		//-----------------------------------  Print out the result in Latex table format ------------------------------------------
		Logger.logln("");
		Logger.logFunction("Winners for " + caption, DebugMode.result);
		//Print out winners
//		for (int am = 0; am < agreementMethods.size(); am++) {
			//winners.add(new Vector<Vector<Winner>>());
			for (int cm =0; cm< correlationMethods.size(); cm ++ ){
				Logger.logFunction(agreementMethods.get(agreementMethod) + " : " +correlationMethods.get(cm) , DebugMode.result);
				Logger.logln("\\begin{table}", DebugMode.result);
				Logger.logln("\\begin{tabular}{| l |l | l| l|l| l|}", DebugMode.result);
				Logger.logln("\\hline" , DebugMode.result);
				Logger.log("Rank & Criterion & " + agreementMethods.get(agreementMethod) + " " , DebugMode.result);
				for (int am2 = 0; am2 < agreementMethods.size(); am2++) if(agreementMethod!=am2){
					Logger.log( " & " + agreementMethods.get(am2) +" " , DebugMode.result);
				}
				Logger.logln(" \\\\ " , DebugMode.result);
				Logger.logln("\\hline" , DebugMode.result);
				for (int j = 0; j < winners.get(agreementMethod).get(cm).size() && j<TOPWINNERS; j++) {
					Logger.log((j+1) + " & " , DebugMode.result);
					Logger.log(winners.get(agreementMethod).get(cm).get(j) , DebugMode.result);
					
					for (int am2 = 0; am2 < agreementMethods.size(); am2++) if(agreementMethod!=am2){
						Logger.log( " & " + (winners.get(am2).get(cm).indexOf(winners.get(agreementMethod).get(cm).get(j)) + 1) , DebugMode.result);
					}
					Logger.logln(" \\\\ " , DebugMode.result);
				}
				
				Logger.logln("\\hline  ", DebugMode.result);
				Logger.logln("\\end{tabular}", DebugMode.result);
				Logger.logln("\\caption{"+caption+" "+new File(expPath).getAbsolutePath()+" Dataset: "+ agreementMethods.get(agreementMethod) + " : " +correlationMethods.get(cm) +"}",DebugMode.result);
				Logger.logln("\\label{table:res1}", DebugMode.result);
				Logger.logln("\\end{table}", DebugMode.result);
			}
//		}
	}
	
	boolean plot = false;
	void plot( NetworkWithGroundTruth dataset, PartiotioningAgreement<V> am, RelativeCriteria<V> crit, final Vector<? extends Number> a, final Vector<? extends Number> b,Vector<Vector<Set<V>>> partitionings){
		XYSeriesCollection collection = new XYSeriesCollection();
		if(a.size()!=b.size()) return;
		
		Logger.logln(dataset.name+"_"+am.toString()+"_"+crit.toString(), DebugMode.detailed);
		XYSeries seriesA = new XYSeries(am.toString());
		XYSeries seriesB = new XYSeries(crit.getName());
	
	//sort them first based on ARI	
		Vector<Integer> ranks = new Vector<Integer>();
		for (int i = 0; i < a.size(); i++) {
			ranks.add(i);
		}
		
		Collections.sort(ranks, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				return (new Double(a.get(o2).doubleValue())).compareTo(a.get(o1).doubleValue()) ;
			}
		});
		
		
		double mina=Double.MAX_VALUE, maxa=Double.MIN_VALUE;
		for (int i = 0; i< a.size(); i++) {
			if(a.get(ranks.get(i)).doubleValue() < mina) mina = a.get(ranks.get(i)).doubleValue();
			if(a.get(ranks.get(i)).doubleValue() > maxa) maxa = a.get(ranks.get(i)).doubleValue();
		}
		
		Logger.log("{", DebugMode.detailed);
		for (int i = 0; i< a.size(); i++) {
			seriesA.add(i,a.get(ranks.get(i)));
			Logger.log(a.get(ranks.get(i))+",", DebugMode.detailed);
		}
		Logger.logln("}", DebugMode.detailed);
		
		double minb=Double.MAX_VALUE, maxb=Double.MIN_VALUE;
		for (int i = 0; i< b.size(); i++) {
			if(b.get(ranks.get(i)).doubleValue() < minb) minb = b.get(ranks.get(i)).doubleValue();
			if(b.get(ranks.get(i)).doubleValue() > maxb) maxb = b.get(ranks.get(i)).doubleValue();
		}
		
		Logger.log("{", DebugMode.detailed);
		for (int i = 0; i< b.size(); i++) {
			seriesB.add(i,((b.get(ranks.get(i)).doubleValue() -minb)  / (maxb-minb) )*(maxa-mina) + mina);
			Logger.log( ((b.get(ranks.get(i)).doubleValue() -minb) / (maxb-minb))*(maxa-mina) + mina +",", DebugMode.detailed);
		}
		Logger.logln("}", DebugMode.detailed);

		collection.addSeries(seriesB);
		collection.addSeries(seriesA);
		JFreeChart chart = ChartFactory.createXYLineChart("", "Partitionings", "Value", collection, PlotOrientation.VERTICAL, true, true, false);
//		JFreeChart chart = ChartFactory.createXYLineChart(dataset.name.substring(0,dataset.name.indexOf('.'))+" "+am.toString()+" "+crit.toString(), "Sample Partitionings", "Value", collection, PlotOrientation.VERTICAL, true, true, false);
	//	chart.setBackgroundPaint(new Color(0, 0, 0, 0));
		chart.getPlot().setBackgroundPaint(Color.white);
		chart.getXYPlot().getRangeAxis().setLabelFont(new Font("sansserif", Font.BOLD, 20));
		
		chart.getXYPlot().getDomainAxis().setLabelFont(new Font("sansserif", Font.BOLD, 20));
		chart.getXYPlot().getRangeAxis().setTickLabelFont(new Font("serif", Font.BOLD, 8));
		chart.getXYPlot().getDomainAxis().setTickLabelFont(new Font("serif", Font.BOLD, 10));
		
		chart.getLegend().setItemFont(new Font("sansserif", Font.BOLD, 24));
		chart.getXYPlot().getRenderer().setSeriesStroke(   0,  new BasicStroke(4.0f));
		chart.getXYPlot().getRenderer().setSeriesStroke(   1,    new BasicStroke(4.0f));
		chart.getXYPlot().getRenderer().setSeriesPaint(0, new Color(0xe05b5b));
		chart.getXYPlot().getRenderer().setSeriesPaint(1, new Color(0x093276));

				//new BasicStroke(       4.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,      1.0f, new float[] {12.0f, 20.0f}, 0.0f    ));
	
		chart.setBackgroundPaint(Color.white);
	//	chart.addSubtitle(new TextTitle( "s="+splitChance+"%, m="+ mergeChance +"%, r="+swapChance+ "%, " + logDatasetStat(partitionings)  ));
	//	ChartUtilities.writeChartAsJPEG(new FileOutputStream("cct2"), 1, chart, 500, 300, info);
		try {
			new File(expPath+"/img/").mkdir();
			File path = new File(expPath+"/img/"+am.toString());
			path.mkdir();
			ChartUtilities.saveChartAsPNG(new File(path.getPath() +"/"+crit.getName()+""+dataset.name+".png"), chart, 600, 400);
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}
		
	/**
	 * @param args
	 */
	public static void main(String[] params) {
		String args[] = params ;
		args = new String[] {"exps/karate", "+G"};
		
		while(args.length<=0 || args[0].contains("help")){
			System.out.println("---- You have following options: ");
			System.out.println("+G: consider randomized versions from the groundtruth");
			System.out.println("+plot: plot all correlation images");
			System.out.println("N=: the number of random clusterings that should be generated");
			System.out.println("s=: split chance for the clustering randomizer");
			System.out.println("m=: merge chance for the clustering randomizer");
			System.out.println("r=: swap chance for the clustering randomizer");
			System.out.println("W=: the number of top/winner criteria to be printed out");
			
			
			System.out.println("Please provide path of the experiments followed by parameters; or type help.");
			Scanner in = new Scanner(System.in);
			args = in.nextLine().split(" ");
//			System.err.println(args.length);
			for (int i = 0; i < args.length; i++) {
//				System.err.println(args[i]);
			}
		    in.close();       
		}
		
		final CriteriaComparer<Integer, Integer> comparer = new CriteriaComparer<Integer, Integer>(args[0]);
		Logger.setLog_level(DebugMode.normal);
		new Logger(comparer.expPath); 
		
		Logger.logFunction(comparer.expPath);
				
//		comparer.addAgreement( new ARI<Integer>());//same //comparer.add( new AGAM <Integer>(GAM.Type.RI));
//		GAM<Integer> gam = new GAM<Integer>();
//		comparer.addAgreement( new AGAM <Integer>(gam.getRIPhi(), new Transformer<Pair<Set<Integer>>, Double>() {
//			public Double transform(Pair<Set<Integer>> uv) {
//				double res = 0;
//				Set<Integer> z = new HashSet<Integer>(uv.getFirst());
//				z.retainAll(uv.getSecond());
//				for (Integer i: z) {
//					res += comparer.nodeWeights.transform(i);
//				}
//				return res;
//			}}));
//
//		
//		comparer.addAgreement(new Rand<Integer>());//same //comparer.add( new GAM <Integer>(GAM.Type.RI));	
//	//	comparer.addAgreement(new FMeasure<Integer>());
//		comparer.addAgreement( new Jaccard<Integer>());
//	
//		comparer.addAgreement( new NMI<Integer>());//same //comparer.add( new NMIAlt<Integer>());
////		comparer.add(new VI<Integer>()); //same//comparer.add(new GAM <Integer>(GAM.Type.VI));
////	//	comparer.add( new AGAM<Integer>(GAM.Type.VI));
//		comparer.addAgreement( new AMI<Integer>());
		
//		comparer.addCorrelation( new PearsonCorrelation());
//		comparer.addCorrelation( new SpearmanCorrelation());
		
		comparer.splitChance = 20;
		comparer.mergeChance = 30;
		comparer.swapChance = 25;
		comparer.maxResult= 60;//180
		comparer.TOPWINNERS = 10000;//20;
		
		
		//exps/test/ +G S<20 s=50 M<30

		for (String string : args) {
			if(string.contains("+C"))
				comparer.addCMs = true;
			if(string.contains("+H"))
				comparer.addHierar = true;
			if(string.contains("+K"))
				comparer.addkmeans = true;
			if(string.contains("+G"))
				comparer.addGTalt = true;
			if(string.contains("+plot"))
				comparer.plot = true;
			if(string.contains("N="))
				comparer.maxResult = Integer.parseInt(string.substring(2));
			if(string.contains("s="))
				comparer.splitChance = Integer.parseInt(string.substring(2));
			if(string.contains("m="))
				comparer.mergeChance = Integer.parseInt(string.substring(2));
			if(string.contains("r="))
				comparer.swapChance = Integer.parseInt(string.substring(2));
			if(string.contains("L="))
				Logger.setLog_level(DebugMode.values()[  Integer.parseInt(string.substring(2)) ] );
			if(string.contains("W="))
				comparer.TOPWINNERS = Integer.parseInt(string.substring(2));
			if(string.contains("help")){
				System.out.println("You have following options: ");
				System.out.println("+G: consider randomized versions from the groundtruth");
				System.out.println("+plot: plot all correlation images");
				System.out.println("N=: the number of random clusterings that should be generated");
				System.out.println("s=: split chance for the clustering randomizer");
				System.out.println("m=: merge chance for the clustering randomizer");
				System.out.println("r=: swap chance for the clustering randomizer");
				System.out.println("W=: the number of top/winner criteria to be printed out");
				return;
			}

		}
		Logger.logln("N="+ comparer.maxResult +", s="+ comparer.splitChance  +", m="+ comparer.mergeChance +", r="+ comparer.swapChance );
		comparer.compareCriteria();
	}

}
