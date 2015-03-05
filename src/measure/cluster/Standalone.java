package measure.cluster;

import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import measure.MeasuresUtil;
import measure.MeasuresUtil.Implementation;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.distance.AlgebricClusteringAgreement;

import org.apache.commons.collections15.Transformer;

import data.GraphDataSet;
import util.DatasetUtils;
import util.IOUtils;
import util.DatasetUtils.DummyDataset;
import algorithms.communityMining.data.Grouping;

public class Standalone {
	public <V,E> void compare(String netPath,  String g1Path,  String g2Path,  
								boolean doExacts , boolean doRI, boolean doNorm , boolean doTrace ,
								boolean doVI , boolean doNMIsum , boolean doNMIsqrt ,boolean  doAMI,
								boolean doJacc ,boolean  doFM,
								boolean doQ  , boolean doNMIvars , boolean doOmega ,
								boolean doAll , boolean doDegreeWeightedStructureBased  , 
								boolean doEdgeStructureBased  , boolean doTransStructureBased  , boolean doSumStructureBased  ,
								boolean doWithStructureTrans  , boolean doWithStructurePlus  ,MeasuresUtil.Implementation implementation){
		
		
		GraphDataSet<V,E> dataset= new GraphDataSet<V, E>(netPath);
	
		Transformer<String , V> vertexTransformer ;
		
		if(netPath!=null) {
			dataset = DatasetUtils.<V,E>load(netPath);
			final Map<String, V> labels_vertices = dataset.labels_vertices;
			System.err.println(labels_vertices);
			vertexTransformer = new Transformer<String, V>() {
				@Override
				public V transform(String input) {
					return labels_vertices.get(input);
				}
			};
		}else {
			vertexTransformer = new Transformer<String, V>() {
				@Override
				public V transform(String input) {
					return (V)input;
				}
			};
		}
		Grouping<V> grouping1=null, grouping2 =null; 
		ListGrouingReader<V > grouingReader = new ListGrouingReader<>();
		try {
			grouping1 = grouingReader.readPartitioning(new FileInputStream(g1Path), vertexTransformer);
			grouping2 = grouingReader.readPartitioning(new FileInputStream(g2Path), vertexTransformer);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		compare(dataset, grouping1, grouping2, doExacts, doRI, doNorm, doTrace, doVI,
				doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doQ, doNMIvars, doOmega, 
				doAll, doDegreeWeightedStructureBased, doEdgeStructureBased, 
				doTransStructureBased, doSumStructureBased, doWithStructureTrans, 
				doWithStructurePlus, implementation);
	}
	public <V,E> void compare(GraphDataSet<V,E> dataset, Grouping<V> grouping1,  Grouping<V> grouping2,  
				boolean doExacts , boolean doRI, boolean doNorm , boolean doTrace ,
				boolean doVI , boolean doNMIsum , boolean doNMIsqrt ,boolean  doAMI,
				boolean doJacc ,boolean  doFM,
				boolean doQ  , boolean doNMIvars , boolean doOmega ,
				boolean doAll , boolean doDegreeWeightedStructureBased  , 
				boolean doEdgeStructureBased  , boolean doTransStructureBased  , boolean doSumStructureBased  ,
				boolean doWithStructureTrans  , boolean doWithStructurePlus  ,MeasuresUtil.Implementation implementation){
		Collection<V> datapoints= null;
		if (dataset!=null && dataset.graph!=null)
			datapoints = dataset.graph.getVertices();
		if(datapoints==null) 
			datapoints = AlgebricClusteringAgreement.getAllDatapoints( 
					grouping1.getGroups(), grouping2.getGroups());
		
		if (doAll){
			 doRI= doNorm= doTrace=	doVI= doNMIsum= doNMIsqrt= doAMI=doJacc = doFM= doQ= doNMIvars= doOmega= doAll=true;
		}
		if (dataset.graph==null && (doDegreeWeightedStructureBased || doEdgeStructureBased || doTransStructureBased 
				||doSumStructureBased||doWithStructureTrans||doWithStructurePlus)){
			System.err.println("Structure based measures only applicaple if the structure, i.e. graph, is given as input. Type -h for help.");
		} else {
			doDegreeWeightedStructureBased= doEdgeStructureBased= doTransStructureBased=
					doSumStructureBased= doWithStructureTrans= doWithStructurePlus = true;
		}
//		Vector<ClusteringAgreement<V>> measures = MeasuresUtil.getAgreementAlternatives(dataset.graph, null);
		Vector<ClusteringAgreement<V>> measures = MeasuresUtil.getAgreements(dataset, doExacts, doRI, doNorm,
				doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega, 
				 doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased, 
				 implementation);

		for (ClusteringAgreement<V> clusteringAgreement : measures) {
			System.out.print(clusteringAgreement.toLatexString()+ " , " );
		}
		System.out.println();
		for (ClusteringAgreement<V> clusteringAgreement : measures) {
			System.out.print(clusteringAgreement.getAgreementCovering(datapoints ,grouping1.getGroups() , grouping2.getGroups()) + " , ");
		}
		System.out.println();
		
//		TODO: Could make it faster this way
//		if(doAll){
//			measures.addAll( MeasuresUtil.getAgreementAlternatives(
//					dataset!=null?dataset.graph:null, dataset!=null?dataset.getWeights():null));
//			if(dataset!=null){
//				Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(dataset.graph, dataset.getWeights(), grouping1.getGroups(), grouping2.getGroups(),false, doTransStructureBased);
//			}
//		}

		
		
	}
	
	public static void main (String[] args){
		// ./compare grouping1.pairs grouping2.pairs [-g network.pairs] -
		// Results in :
		// ARI' 
		// Other Measures: RI', RI, ARI
		
		String netPath=null, g1Path, g2Path;
		boolean doExacts = false, doRI=false, doNorm = false, doTrace=false;
		boolean doVI=false, doNMIsum=false, doNMUsqrt=false, doAMI = false, doJacc =false, doFM =false;
			
		boolean doQ = false, doNMIvars= false, doOmega=false;
		boolean doAll = false;
		boolean doDegreeWeightedStructureBased = false, doEdgeStructureBased = false, 
				doTransStructureBased = false, doSumStructureBased = false;
		boolean doWithStructureTrans = false, doWithStructurePlus = false;
		MeasuresUtil.Implementation implementation = null;
		
		//To test
		args ="NMIexample_v NMIexample_u1 -g ./NMIexample_net.gml -exact".split("[,\\s]+");
//		for (int i = 0; i < args.length; i++) {
//			System.err.print(i+": "+args[i]+", ");
//		}System.err.println();
		
		if(args.length>0 && args[0].equals("-h")){
			System.out.println("Usage: ./compare grouping1 grouping2 [-g network.pairs]\n"
					+ " groupings format: each line one group, lists nodes in that group"
					+ "-all : compute all measures \n"
					+ "-exact : compute exact variation, i.e. do not consider pairs of same-nodes \n"
					+ "+ri : also compute the RI \n"
					+ "+norm : also compute the norm variation\n"
					+ "+trace : also compute the trace variation \n"
					
					+ "+vi : also compute VI \n"
					+ "+nmi : also compute NMI_sum \n"
					+ "+nmiSq : also compute NMI_sqrt \n"
					+ "+ami : also compute AMI \n"
					+ "+Jacc: also compute Jaccard"
					+ "+F: also compute FMeasure"
					
					+ "+sbwd : also compute structure based weighted degree overlap  \n"
					+ "+sbwe : also compute structure based edge overlap  \n"
					+ "+sbt : also compute structure based transposed \n"
					+ "+sbs : also compute  structure based sum variation\n"
					
//					+ "+wst : also compare with structure  transposed \n"
//					+ "+wsp : also compare with structure  sum variation\n"
//					
					+ "+omega : also compute omega index \n"
					+ "+nmiVars : also compute  overlapping nmi variations\n"
					
					+ "-originals: use original implementations if available \n"
					+ "-ethaBased: use implementations of generalized formula based on eta if available \n"
					+ "-algCont: use algebric implementations based on contingency table if available \n"
					+ "-algDelta [default]: use algebric implementations based on delta if available \n"

//					+ "+q: also compute q modularities measures  \n"
					);
			return;
		}
		if(args.length<2 ){
			System.err.println("Takes at least two arguments: name of the two files to compare. Type -h for help.");
			return;
		}
		g1Path = args[0];
		g2Path = args[1];
		
		for (int i = 2; i < args.length; i++) {
			switch (args[i]) {
			case "-g":
				if(i+1<args.length)
					netPath = args[i+1];
				else {
					System.err.println("Should provide the path to the network file. Type -h for help.");
					return;
				}
			case "-all":
				doAll = true;
				break;
			case "-exact":
				doExacts = true;
				break;
			case "+ri":
				doRI = true;
				break;
			case "+norm":
				doNorm = true;
				break;
			case "+trace":
				doTrace = true;
				break;
			case "+vi":
				doVI = true;
				break;
			case "+nmi":
				doNMIsum = true;
				break;
			case "+nmiSq":
				doNMUsqrt = true;
				break;
			case "+ami":
				doAMI = true;
				break;
			case "+sbwd":
				doDegreeWeightedStructureBased = true;
				break;
			case "+sbwe":
				doEdgeStructureBased = true;
				break;
			case "+sbt":
				doTransStructureBased = true;
				break;
			case "+sbs":
				doSumStructureBased = true;
				break;
//			case "+wst":
//				doWithStructureTrans = true;
//				break;
//			case "+wsp":
//				doWithStructurePlus = true;
//				break;
			case "+omega":
				doOmega = true;
				break;
			case "+nmiVars":
				doNMIvars = true;
				break;
			case "-originals":
				if (implementation!=null)
					System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.ORIGINAL;
			case "-ethaBased":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.GAM;
			case "-algCont":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.ALGEBRIC_OVERLAP;
			case "-algDelta":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.ALGEBRIC_DELTA;
			case "+q":
				doQ = true;
				break;
			default:
				break;
			}
		}
		
		Standalone  standalone = new Standalone();
		standalone.compare(netPath, g1Path, g2Path, doExacts, doRI, doNorm, doTrace,
				doVI, doNMIsum, doNMUsqrt, doAMI,doJacc , doFM, 
				doQ, doNMIvars, doOmega, doAll, 
				doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased,
				doSumStructureBased, doWithStructureTrans, doWithStructurePlus,implementation);
		
		//Writing Examples
//		for(DummyDataset dummyDataset: DummyDataset.values()){
//			GraphDataSet<Integer, Integer> dataset = DatasetUtils.loadDummy(dummyDataset);
//			try {
//				IOUtils.writeGML(dummyDataset+"_net.gml", dataset.graph, null, dataset.weights, null, null);
//				IOUtils.writeListGrouping(dummyDataset+"_v", dataset.getGrouping("V"),dataset.getLabels());
//				IOUtils.writeListGrouping(dummyDataset+"_u1", dataset.getGrouping("U1"),dataset.getLabels());
//				IOUtils.writeListGrouping(dummyDataset+"_u2", dataset.getGrouping("U2"),dataset.getLabels());
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
//		for (boolean doExact :new boolean []{true,false})
//		for(Implementation imp:Implementation.values()){
////			System.out.println("-----------------"+imp+"-----------------");
//			standalone.compare(dataset, dataset.getGrouping("V"), dataset.getGrouping("U1"), doExact, 
//				doRI, doNorm, doTrace,
//				doVI, doNMIsum, doNMUsqrt, doAMI,doJacc , doFM, 
//				doQ, doNMIvars, doOmega, doAll, 
//				doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased,
//				doSumStructureBased, doWithStructureTrans, doWithStructurePlus,imp);
//		}
	
	}
	
}
