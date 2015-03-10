package measure.cluster;

import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane.MaximizeAction;

import measure.MeasuresUtil;
import measure.MeasuresUtil.Implementation;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.distance.AlgebricClusteringAgreement;

import org.apache.commons.collections15.Transformer;
import org.jfree.data.statistics.MeanAndStandardDeviation;

import data.GraphDataSet;
import util.DatasetUtils;
import util.IOUtils;
import util.DatasetUtils.DummyDataset;
import algorithms.communityMining.data.Grouping;

public class Standalone {
	DecimalFormat decimalFormat = new DecimalFormat("#.###");

	public <V,E> void compare(String netPath,  String g1Path,  String g2Path,  boolean overlapping,
								boolean doExacts , boolean doRI, boolean doNorm , boolean doTrace ,
								boolean doVI , boolean doNMIsum , boolean doNMIsqrt ,boolean  doAMI,
								boolean doJacc ,boolean  doFM,
								boolean doQ  , boolean doNMIvars , boolean doOmega ,
								boolean doAll , boolean doDegreeWeightedStructureBased  , 
								boolean doEdgeStructureBased  , boolean doTransStructureBased  , boolean doSumStructureBased  ,
								boolean doWithStructureTrans  , boolean doWithStructurePlus  ,MeasuresUtil.Implementation implementation, boolean printDetailed){
		
		
		GraphDataSet<V,E> dataset=null;
		Transformer<String , V> vertexTransformer ;
		
		if(netPath!=null) {
			dataset = DatasetUtils.<V,E>load(netPath);
			final Map<String, V> labels_vertices = dataset.labels_vertices;
//			System.err.println(labels_vertices);
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
		
		compare(dataset, grouping1, grouping2,overlapping, doExacts, doRI, doNorm, doTrace, doVI,
				doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doQ, doNMIvars, doOmega, 
				doAll, doDegreeWeightedStructureBased, doEdgeStructureBased, 
				doTransStructureBased, doSumStructureBased, doWithStructureTrans, 
				doWithStructurePlus, implementation, printDetailed);
	}
	public <V,E> void compare(GraphDataSet<V,E> dataset, Grouping<V> grouping1,  Grouping<V> grouping2,  
				boolean overlapping,boolean doExacts , boolean doRI, boolean doNorm , boolean doTrace ,
				boolean doVI , boolean doNMIsum , boolean doNMIsqrt ,boolean  doAMI,
				boolean doJacc ,boolean  doFM,
				boolean doQ  , boolean doNMIvars , boolean doOmega ,
				boolean doAll , boolean doDegreeWeightedStructureBased  , 
				boolean doEdgeStructureBased  , boolean doTransStructureBased  , boolean doSumStructureBased  ,
				boolean doWithStructureTrans  , boolean doWithStructurePlus  ,MeasuresUtil.Implementation implementation, boolean printDetailed){
		Collection<V> datapoints= null;

		if (dataset!=null && dataset.graph!=null)
			datapoints = dataset.graph.getVertices();
		if(datapoints==null) 
			datapoints = AlgebricClusteringAgreement.getAllDatapoints(grouping1.getGroups(), grouping2.getGroups());

		if (doAll){
			if (overlapping)
				doRI=doNorm=doTrace=doOmega=true;
			else if (implementation==Implementation.ALGEBRIC_DELTA)
				doRI=doNorm=doTrace=doVI=doNMIsum=doJacc=doFM= true;
			else 
				doRI=doVI= doNMIsum=doJacc=doFM= true;
			
		}		

		if ((dataset!=null && dataset.graph!=null) || doDegreeWeightedStructureBased || doEdgeStructureBased || doTransStructureBased 
				||doSumStructureBased||doWithStructureTrans||doWithStructurePlus){
			if (dataset.graph==null)
				System.err.println("Structure based measures only applicaple if the structure, i.e. graph, is given as input. Type -h for help.");
			else {//if (doAll){
				if (overlapping){
					doTransStructureBased=doSumStructureBased= doWithStructureTrans= doWithStructurePlus = true;		
				}else {
					doDegreeWeightedStructureBased= doEdgeStructureBased= true;
				}
				if (implementation==Implementation.ALGEBRIC_DELTA) {
					doTransStructureBased=doSumStructureBased= doWithStructureTrans= doWithStructurePlus = true;		
				}
			}	
		}
//		Vector<ClusteringAgreement<V>> measures = MeasuresUtil.getAgreementAlternatives(dataset.graph, null);
		Vector<ClusteringAgreement<V>> measures = MeasuresUtil.getAgreements(dataset, overlapping, doExacts, doRI, doNorm,
				doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega, 
				 doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased, 
				 implementation);

		Vector<String> formats = new Vector<String>();
		Vector<String> names = new Vector<String>();
		Vector<String> values = new Vector<String>();
		
		String  format ;int detailedNameStartsFrom;
		double agreement ;
		for (int i=0;i<measures.size();i++){
			names.add(measures.get(i).toLatexString());
			detailedNameStartsFrom = names.get(i).indexOf(':');
			if(!printDetailed && detailedNameStartsFrom!=-1)
				names.set(i, names.get(i).substring(0, detailedNameStartsFrom));
			
			agreement = measures.get(i).getAgreementCovering(datapoints ,grouping1.getGroups() , grouping2.getGroups());
			values.add((printDetailed?agreement: decimalFormat.format(agreement))+"");
			
			formats.add("%-"+(Math.max(names.lastElement().length(), values.lastElement().length())+4)+"s");
		}
		for (int i=0;i<measures.size();i++){
			System.out.format(formats.get(i),   (i!=0?", ":"")+  names.get(i));
		}
		System.out.println();
		for (int i=0;i<measures.size();i++){
			System.out.format(formats.get(i),   (i!=0?", ":"")+  values.get(i));
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
	
	public void Compare(String[] args){
		
		// ./compare grouping1.pairs grouping2.pairs [-g network.pairs] -
		// Results in :
		// ARI' 
		// Other Measures: RI', RI, ARI
		
		String netPath=null, g1Path, g2Path;
		boolean doExacts = false, doRI=false, doNorm = false, doTrace=false;
		boolean doVI=false, doNMIsum=false, doNMUsqrt=false, doAMI = false, doJacc =false, doFM =false;
		boolean overlapping = false;
		boolean doQ = false, doNMIvars= false, doOmega=false;
		boolean doAll = false;
		boolean doDegreeWeightedStructureBased = false, doEdgeStructureBased = false, 
				doTransStructureBased = false, doSumStructureBased = false;
		boolean doWithStructureTrans = false, doWithStructurePlus = false;
		MeasuresUtil.Implementation implementation = null;
		boolean printDetailed=false;

		//To test
		//		for (int i = 0; i < args.length; i++) {
//			System.err.print(i+": "+args[i]+", ");
//		}System.err.println();
		
		if(args.length>0 && ( args[0].equals("-h")||args[0].equals("--help"))){
			System.out.println("Usage: ./compare grouping1 grouping2 [-g network]\n"
					+ "\t-all    : compute all measures \n"
					+ "\t-o      : compute overlapping variations\n"
					+ "\t-exact  : compute exact variation, i.e. do not consider pairs of same-nodes \n"
					+ "\t+ri     : also compute the RI \n"
					+ "\t+norm   : also compute the norm variation\n"
					+ "\t+trace  : also compute the trace variation \n"
					+ "\t+details: print detailed names for measures \n"
					
					+ "\t+vi     : also compute VI \n"
					+ "\t+nmi    : also compute NMI_sum \n"
					+ "\t+nmiSq  : also compute NMI_sqrt \n"
					+ "\t+ami    : also compute AMI \n"
					+ "\t+Jacc   : also compute Jaccard\n"
					+ "\t+F      : also compute FMeasure\n"
					
					+ "\t+sbwd   : also compute structure based weighted degree overlap  \n"
					+ "\t+sbwe   : also compute structure based edge overlap  \n"
					+ "\t+sbt    : also compute structure based transposed \n"
					+ "\t+sbs    : also compute  structure based sum variation\n"
					
//					+ "\t+wst : \t also compare with structure  transposed \n"
//					+ "\t+wsp : \t also compare with structure  sum variation\n"
//					
					+ "\t+omega  : also compute omega index \n"
					+ "\t+nmiVars: also compute  overlapping nmi variations\n"
					
					+ "\t-originals : use original implementations \n"
					+ "\t-ethaBased : [default for disjoint] use implementations of generalized formula based on eta \n"
					+ "\t-algCont   : use algebric implementations based on contingency table\n"
					+ "\t-algDelta  : [default for overlapping] use algebric implementations based on delta \n"

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
				break;
			case "-all":
				doAll = true;
				break;
			case "-exact":
				doExacts = true;
				break;
			case "-o":
				overlapping = true;
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
				
			case "+details":
				printDetailed = true;
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
				break;

			case "-ethaBased":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.GAM;
				break;

			case "-algCont":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.ALGEBRIC_OVERLAP;
				break;

			case "-algDelta":				
				if (implementation!=null)
				System.err.println("Only one implementation mode should be chosen, last one would be effective!");
				implementation = Implementation.ALGEBRIC_DELTA;
				break;

			case "+q":
				doQ = true;
				break;
			default:
				break;
			}
		}
		compare(netPath, g1Path, g2Path,overlapping, doExacts, doRI, doNorm, doTrace,
				doVI, doNMIsum, doNMUsqrt, doAMI,doJacc , doFM, 
				doQ, doNMIvars, doOmega, doAll, 
				doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased,
				doSumStructureBased, doWithStructureTrans, doWithStructurePlus,implementation, printDetailed);
		
	
	}
	public void Compare(String command){
		System.err.println("./cag.jar "+command);
		Compare(command.split("[\\s]+"));
	}
	public void test (){
		String comm="";
		comm+="-all  ";
		comm+="-o ";
//		comm+="-algDelta ";
		comm+="+nmiVars ";

		for (String dataset : new String[]{"NMIexample","Omega","Structure"}){
			System.out.println(dataset);
			Compare("examples/"+dataset+"_V.list "+"examples/"+dataset+"_U1.list "+comm);
			Compare("examples/"+dataset+"_V.list "+"examples/"+dataset+"_U2.list "+comm);
			Compare("examples/"+dataset+"_V.list "+"examples/"+dataset+"_U1.list "+"-g examples/"+dataset+".gml "+comm);
			Compare("examples/"+dataset+"_V.list "+"examples/"+dataset+"_U2.list "+"-g examples/"+dataset+".gml "+ comm);

		}
		
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
	
	public static void main (String[] args){
		Standalone  standalone = new Standalone();
//		standalone.test();
		standalone.Compare(args);
	}
	
}
