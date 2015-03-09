package measure.cluster.agreement;

import io.group.GroupingWriter;
import io.group.ListGroupingWriter;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import util.DatasetUtils;
import util.DatasetUtils.DummyDataset;
import util.IOUtils;
import data.GraphDataSet;
import data.Pair;
import measure.MeasuresUtil;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.cluster.agreement.partitioning.classics.AMI;
import measure.cluster.agreement.partitioning.classics.ARI;
import measure.cluster.agreement.partitioning.classics.NMI;
import measure.cluster.agreement.partitioning.classics.NMIAlt;
import measure.cluster.agreement.partitioning.classics.Rand;
import measure.cluster.agreement.partitioning.classics.VI;
import measure.cluster.agreement.partitioning.generalization.AGAM;
import measure.cluster.agreement.partitioning.generalization.AGAMA;
import measure.cluster.agreement.partitioning.generalization.GAM;
import measure.cluster.distance.AlgebricClusteringAgreement;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;

/**
 * Created on Apr 19, 2012
 *
 * A simple test for different partitioning agreement measures
 * 
 * @author Reihaneh Rabbany, see the <a href="{@docRoot}/license.txt">copyright</a>
 *
 */
public class test {

	public static  Set<Integer> createDummy(int[] a){
		Set<Integer> c = new HashSet<Integer>();
		for (int i : a) {
			c.add(i);
		}
		return c;
	}
	
	
	public static void testAltenativeImplementations(){
		Vector<Set<Integer>> p1 = new Vector<Set<Integer>>(); 
		Vector<Set<Integer>> p2 = new Vector<Set<Integer>>();
		Vector<Set<Integer>> groundTruth = new Vector<Set<Integer>>();
	
		p1.add(createDummy(new int[]{1,2,3,4,5,6}));
		p1.add(createDummy(new int[]{7,8,9,10,11,12}));
		p1.add(createDummy(new int[]{13,14,15,16,17}));

		p2.add(createDummy(new int[]{1,4,7,10,13,16}));
		p2.add(createDummy(new int[]{2,5,8,14,11,17}));
		p2.add(createDummy(new int[]{3,6,9,11,15}));
		
		groundTruth.add(createDummy(new int[]{1,3,4,5,6,7,13,15}));
		groundTruth.add(createDummy(new int[]{2,8,9,10,12}));
		groundTruth.add(createDummy(new int[]{11,14,16,17}));

		PartiotioningAgreement<Integer> meas = new AMI<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth));
		meas = new AGAM <Integer>(GAM.Type.VI);
		System.err.println(meas + ": " + (meas.getAgreement(p1, groundTruth)));
		System.err.println();
		
		meas = new ARI<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth)+" == "+((ARI<Integer>)meas).getAgreementStandalone(p1, groundTruth)+" == "+((ARI<Integer>)meas).getAgreementStandaloneAlt(null,p1, groundTruth));
	
		meas = new AGAM <Integer>(GAM.Type.RI);
		System.err.println(meas + ": " + (meas.getAgreement(p1, groundTruth)));
		System.err.println();
		
		
		meas = new AGAMA <Integer>(GAM.Type.RI);
		System.err.println(meas + ": " + (meas.getAgreement(p1, groundTruth)));
		System.err.println();
		//		meas = new FMeasure<Integer>();
//		System.err.println(meas + ": " +meas.getAgreement(partitioning, groundTruth)+" == "+((FMeasure<Integer>)meas).getAgreementAlt(partitioning, groundTruth));
//		meas = new Jaccard<Integer>();
//		System.err.println(meas + ": " +meas.getAgreement(partitioning, groundTruth)+" == "+((Jaccard<Integer>)meas).getAgreementAlt2(partitioning, groundTruth) +" == "+((Jaccard<Integer>)meas).getAgreementAlt(partitioning, groundTruth));
		meas = new NMI<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth)+ ": " +((NMI<Integer>)meas).getAgreementStandalone(p1, groundTruth));
	
		meas = new NMIAlt<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth));
		System.err.println();

		meas = new Rand<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth));
		
		meas = new GAM <Integer>(GAM.Type.RI);
		System.err.println(meas + ": " + (meas.getAgreement(p1, groundTruth)));
	//	System.err.println(meas + ": " + (meas.getAgreement(p2, groundTruth)));
		System.err.println();

		meas = new VI<Integer>();
		System.err.println(meas + ": " +meas.getAgreement(p1, groundTruth)+ ": " +((VI<Integer>)meas).getAgreementAlt(p1, groundTruth)+ ": " +((VI<Integer>)meas).getAgreementAlt2(p1, groundTruth));
		
		meas = new GAM <Integer>(GAM.Type.VI);
		System.err.println(meas + ": " + (meas.getAgreement(p1, groundTruth)));
		//System.err.println(meas + ": " + (meas.getAgreement(partitioning, groundTruth)));
		
	}
	
	public static void graphExample(){
		System.err.println("Graph Example");
		GraphDataSet<Integer, Integer> dataset = DatasetUtils.loadDummy(DummyDataset.Structure);
		Graph<Integer, Integer> g =  dataset.graph;
		Vector<Set<Integer>> groundTruth = dataset.getGrouping("V").getGroups();//DatasetUtils.createClusteringFromArray(new int[][]{ {0,1,2,3,4,5},{6,7,8}});
		Vector<Set<Integer>> p1 = dataset.getGrouping("U1").getGroups();//DatasetUtils.createClusteringFromArray(new int[][]{ {1,2,3,4,5},{0,6,7,8}});
		Vector<Set<Integer>> p2 = dataset.getGrouping("U2").getGroups();//DatasetUtils.createClusteringFromArray(new int[][]{ {0,1,2,3,4},{5,6,7,8}});
		
		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer,Integer>getAgreements(dataset,true, null)) {
			System.err.println(meas.toLatexString());
			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
		}
		
		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer,Integer>getAgreementAlternatives(g, null)) {
			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
		}
		
		Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(g,null, groundTruth, groundTruth);
		Pair<Vector<String>, Vector<Double>> tmp1 = AlgebricClusteringAgreement.getAllAgreements(g,null, p1, groundTruth);
		Pair<Vector<String>, Vector<Double>> tmp2 = AlgebricClusteringAgreement.getAllAgreements(g,null, p2, groundTruth);

		for (int i = 0; i < tmp.first.size(); i++) {
			System.err.println(tmp.first.get(i) + "(U,V) =" +tmp.second.get(i));
			System.err.println(tmp1.first.get(i) + "(U1,V) =" +tmp1.second.get(i));
			System.err.println(tmp2.first.get(i) + "(U2,V) =" +tmp2.second.get(i));
		}
	
	}
	
	public static void omegaExample(){
	GraphDataSet<Integer, Integer> dataset = DatasetUtils.loadDummy(DummyDataset.Omega);
	Graph<Integer, Integer> g =  dataset.graph;
	Vector<Set<Integer>> groundTruth = dataset.getGrouping("V").getGroups();
	System.err.println(groundTruth);
	Vector<Set<Integer>> p1 = dataset.getGrouping("U1").getGroups();
	Vector<Set<Integer>> p2 = dataset.getGrouping("U2").getGroups();
		
		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer,Integer>getAgreementAlternatives(g, null)) {
			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
		}
//		Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(g,null, groundTruth, groundTruth);
//		Pair<Vector<String>, Vector<Double>> tmp1 = AlgebricClusteringAgreement.getAllAgreements(g,null, p1, groundTruth);
//		Pair<Vector<String>, Vector<Double>> tmp2 = AlgebricClusteringAgreement.getAllAgreements(g,null, p2, groundTruth);

//		System.err.println("-------------------------");
//		for (int i = 0; i < tmp2.first.size(); i++) {
//			System.err.println(tmp.first.get(i) + "(V,V) =" +tmp.second.get(i));
//			System.err.println(tmp1.first.get(i) + "(U1,V) =" +tmp1.second.get(i));
//			System.err.println(tmp2.first.get(i) + "(U2,V) =" +tmp2.second.get(i));
//
//		}
		
	
	}
	
	public static void testExamples(){
		//Example for set matching
//		Vector<Set<Integer>> groundTruth = createClusteringFromArray(new int[][]{ {1,2,3,4,5},{6,7,8},{9,10,11}});
//		Vector<Set<Integer>> p1 = createClusteringFromArray(new int[][]{ {1,2,3,4,5},{6,7,8,9,10,11}});
//		Vector<Set<Integer>> p2 = createClusteringFromArray(new int[][]{ {1,2,3},{6,7,8},{4,5,9,10,11}});
		Vector<Set<Integer>> groundTruth = DatasetUtils.createClusteringFromArray(new int[][]{ {0,1,2,3,4,5},{6,7,8}});
		Vector<Set<Integer>> p1 = DatasetUtils.createClusteringFromArray(new int[][]{ {1,2,3,4,5},{0,6,7,8}});
		Vector<Set<Integer>> p2 = DatasetUtils.createClusteringFromArray(new int[][]{ {1,2,3,4},{5,6,7,8}});
//		Vector<Set<Integer>> t = createClusteringFromArray(new int[][]{ {1,2,3,4,5},{4,5,6,7,8}});
		
		// Melia example
//		groundTruth = createClusteringFromArray(new int[][]{{1 ,2, 3 ,4 ,5, 6 ,7 ,8 ,9 ,10, 11, 12},
//		{13, 14 ,15, 16, 17 ,18, 19, 20, 21 ,22 ,23 ,24},
//		{25, 26, 27, 28, 29, 30, 31, 32, 33, 34 ,35 ,36}});
//
//		p1 = createClusteringFromArray(new int[][]{{13 ,14, 15 ,16, 5 ,6 ,7 ,8 ,9 ,10 ,11 ,12},
//			{25 ,26,27,28 ,17 ,18, 19 ,20 ,21, 22, 23, 24},
//			{1 ,2,3,4, 29, 30, 31 ,32 ,33 ,34, 35 ,36}});
//		p2 = createClusteringFromArray(new int[][]{{13, 14 ,25 ,26, 5, 6, 7 ,8 ,9 ,10 ,11, 12},
//		{1, 2, 27 ,28 ,17 ,18 ,19 ,20, 21 ,22,23 ,24},
//		{3, 4 ,15 ,16, 29, 30 ,31 ,32 ,33, 34 ,35 ,36}});

		// Melia simpler example
//		groundTruth = createClusteringFromArray(new int[][]{{1 ,2, 3 ,4 ,5, 6 ,7 ,8 },
//		{13, 14 ,15, 16, 17 ,18, 19, 20 },
//		{25, 26, 27, 28, 29, 30, 31, 32}});
//
//		p1 = createClusteringFromArray(new int[][]{{13 ,14, 15 ,16, 5 ,6 ,7 ,8 },
//			{25 ,26,27,28 ,17 ,18, 19 ,20 },
//			{1 ,2,3,4, 29, 30, 31 ,32 }});
//		p2 = createClusteringFromArray(new int[][]{{13, 14 ,25 ,26, 5, 6, 7 ,8},
//		{1, 2, 27 ,28 ,17 ,18 ,19 ,20},
//		{3, 4 ,15 ,16, 29, 30 ,31 ,32 }});

		
		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer>getAllAgreementImplementations()) {
//			System.err.println(meas + "(t,t) =" +meas.getAgreement(t, t));
//			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
//			System.err.println();
		}
		
		Graph<Integer, Integer> g =  new SparseGraph<>() ;
		int c=0;
		g.addEdge(c++,0,1);g.addEdge(c++,0,5);g.addEdge(c++,0,6);
		g.addEdge(c++,1,0);g.addEdge(c++,1,2);g.addEdge(c++,1,5);
		g.addEdge(c++,2,1);g.addEdge(c++,2,3);g.addEdge(c++,2,5);
		g.addEdge(c++,3,2);g.addEdge(c++,3,5);g.addEdge(c++,3,4);
		g.addEdge(c++,4,3);g.addEdge(c++,4,5);
		g.addEdge(c++,5,0);g.addEdge(c++,5,1);g.addEdge(c++,5,2);g.addEdge(c++,5,3);g.addEdge(c++,5,4);g.addEdge(c++,5,6);g.addEdge(c++,5,8);
		g.addEdge(c++,6,0);g.addEdge(c++,6,5);g.addEdge(c++,6,7);g.addEdge(c++,6,8);
		g.addEdge(c++,7,6);g.addEdge(c++,7,8);
		g.addEdge(c++,8,6);g.addEdge(c++,8,7);g.addEdge(c++,8,5);
		
		
		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer,Integer>getAgreementAlternatives(g, null)) {
//			System.err.println(meas + "(t,t) =" +meas.getAgreement(t, t));
//			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
//			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
//			System.err.println();
		}
		Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(g,null, p1, groundTruth);
		for (int i = 0; i < tmp.first.size(); i++) {
			System.err.println(tmp.first.get(i) + "(U1,V) =" +tmp.second.get(i));
		}
		
//		final Map<Integer, Double> weights = new HashMap<Integer, Double>();
//		weights.put(0, 1.);
		
		
		//Example for omega
//		Vector<Set<Integer>> groundTruth = createClusteringFromArray(new int[][]{ {1,2,3,4,5},{2,3,4},{3,4,5}});
//		Vector<Set<Integer>> p1 = createClusteringFromArray(new int[][]{ {1,4,5},{2},{3}});
//		Vector<Set<Integer>> p2 = createClusteringFromArray(new int[][]{ {1,4,5},{2},{3,4}});
//		
//		for (PartiotioningAgreement<Integer> meas : MeasuresUtil.<Integer>getAgreementAlternatives()) {
//			
//			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
//			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
//			System.err.println();
//		}
	}
	
	public static void setmatching_example(){
		GraphDataSet<Integer, Integer> dataset = DatasetUtils.loadDummy(DummyDataset.NMIexample);
		Graph<Integer, Integer> g =  dataset.graph;
		
		
		Vector<Set<Integer>> V = dataset.getGrouping("V").getGroups();
		Vector<Set<Integer>> U1 = dataset.getGrouping("U1").getGroups();
		Vector<Set<Integer>> U2 = dataset.getGrouping("U2").getGroups();
			
		ClusteringAgreement<Integer> onmi = new NMI_Overlapping_MacDaid<>();
		System.err.println("-- "+onmi.getAgreement(U1, V));
		System.err.println("-- "+onmi.getAgreement(U2, V));
////		ClusteringAgreement<Integer> onmi2 = new NMI_Overlapping_LFR<>();
////		System.err.println("-- "+onmi2.getAgreement(U1, V));
////		System.err.println("-- "+onmi2.getAgreement(U2, V));
//		
		
		
	}
	
	
	/**
	 *
	 * Books example:  <a href="http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html"> Christopher D. Manning, Prabhakar Raghavan and Hinrich Sch��tze, Introduction to Information Retrieval, Cambridge University Press. 2008.</a>
	 *
	 * @param params
	 */
	public static void main(String[] params){
	
//		testExamples();
//		testAltenativeImplementations();
		graphExample();
		omegaExample();
		setmatching_example();
	}


}
