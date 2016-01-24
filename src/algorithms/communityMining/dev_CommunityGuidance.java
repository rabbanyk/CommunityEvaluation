package algorithms.communityMining;

import java.awt.Color;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JSplitPane;

import org.apache.commons.collections15.Transformer;

import measure.MeasuresUtil;
import measure.base.RelativeCriteria;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM;
import measure.cluster.agreement.partitioning.generalization.GAM.Type;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.AdjustionMethod;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.ExternalOverlap;
import measure.cluster.distance.AlgebricClusteringAgreement;
import measure.graph.criteria.Modularity;
import ui.PartitionView;
import util.DatasetUtils;
import util.DatasetUtils.ClassicDataset;
import util.IOUtils;
import io.graph.gml.GMLGraphWriter;
import algorithms.AlgorithmUtils;
import algorithms.AlgorithmUtils.Method;
import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.external_methods.Donetti;
import algorithms.communityMining.external_methods.FastModularity;
import algorithms.communityMining.external_methods.Infomap;
import algorithms.communityMining.external_methods.Louvain;
import algorithms.communityMining.external_methods.PottsModel;
import algorithms.communityMining.external_methods.WalkTrap;
import algorithms.communityMining.topleaders.TopLeaders;
import algorithms.communityMining.topleaders.data.Partitioning;
import data.GraphDataSet;
import data.Pair;
import dev_Experiments.ExternalIndexComparer;
import dev_Experiments.ExternalIndexComparer.Mode;
import edu.uci.ics.jung.graph.Graph;

public class dev_CommunityGuidance {
	DecimalFormat df = new DecimalFormat("#.###");
	DateFormat formatter = new SimpleDateFormat("mm:ss:SSS");  
	int missingValueIndicator = 0;
	boolean doAlgebric = false;
	boolean doMissingValueAsCluster=false,doMissingValueRemove=false,doMissingVlauesUntreated=true;
	static boolean doCommunities = true;
	
	public static void main(String[] args){
		String EXPNAME=".";
		boolean isOverlapping = false, doTrans = false, doRelative = false,
				doAtt = true, doAlphas=false, writeCommRes=true;
		
		if (args.length ==0)
			args = new String[] {"/home/reihaneh/projects/exps/fbtest/data/" , "-a"};

		if(args.length>0) EXPNAME = args[0];
//		deleteDir(new File(EXPNAME));
		new File(EXPNAME).mkdir();
		if(doCommunities)
			createCleanDir(new File(EXPNAME+"com"));


		for (String arg : args) {
			switch (arg) {
			case "-o":
				isOverlapping = true;
				break;
			case "-t":
				doTrans = true;
				break;
			case "-a":
				doAtt = true;
				break;
			case "-q":
				doRelative = true;
				break;
			case "-alpha":
				doAlphas = true;
				break;
			case "-w":
				writeCommRes = true;
				break;
			case "-h":
				System.err.println("./compare  datasetPath (should have a folder called data) \n "
								+  "-o: use overlapping methods and measures \n "
								+  "-a: compute agreement with all attributes \n"
								+  "-q: compute modularity q \n"
								+  "-t: compute transformed variations (time consuming) \n "
								+  "-w: write resulted communities on file \n"
								+  "-alpha: compute the +ARI agreement with different alphas \n");
			default:
				break;
			}
		}
		Iterator<GraphDataSet<Integer, Integer>> datasets = DatasetUtils.<Integer,Integer>loadAll(EXPNAME);//+"data");

		dev_CommunityGuidance test = new dev_CommunityGuidance();
		try {
//			test.compareMethodsWithGT(datasets, EXPNAME, isOverlapping, doQ, doTrans);
//			test.compareMethodsAllAtt(datasets, EXPNAME, isOverlapping, doQ, doTrans);
			test.compareMethods(datasets, EXPNAME, isOverlapping, doAtt, doRelative, doTrans, doAlphas, writeCommRes);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
				
	public <V,E> void compareMethodsWithGT(Iterator<GraphDataSet<V,E>> datasets, String runname, 
			boolean isOverlapping,boolean doRelative, boolean  dotrans) throws IOException{
		compareMethods(datasets, runname, isOverlapping, false, doRelative, dotrans, false, false);
	}
	public <V,E> void compareMethodsAllAtt(Iterator<GraphDataSet<V,E>> datasets, String runname, 
			boolean isOverlapping,boolean doRelative, boolean  dotrans) throws IOException{
		compareMethods(datasets, runname, isOverlapping, true, doRelative, dotrans, false, false);
	}
	
	
	public <V,E> void compareMethods(Iterator<GraphDataSet<V,E>> datasets, String runname,
			boolean isOverlapping, boolean doAtt, boolean doRelative, boolean  dotrans, boolean doalpha, boolean writeResultedCommunities) throws IOException{
		FileOutputStream outResultsE = new FileOutputStream(runname+"resE"+".csv");		
		FileOutputStream outResultsR = new FileOutputStream(runname+"resR"+".csv");			

		long startTime;
		HashMap<String,Long> community_runtimes;
		Vector<String> community_names= new Vector<String>();
		
		//------------------------Write the header of file, measure names, etc.
		outResultsR.write( ("dataset,weights,#nodes,#edges,CM_alg,t,k").getBytes());
		if(doRelative)
			for ( RelativeCriteria<V> relative: MeasuresUtil.<V,E>getRelativeAlternatives(null, null)){
				outResultsR.write( ("," + relative.getName()	).getBytes());
			}
		outResultsR.write(("\n").getBytes());

		outResultsE.write( ("dataset,weights,#nodes,#edges,CM_alg,k"+ (doAtt?",att":"")).getBytes());
		for (ClusteringAgreement<V> clusterAgreement :  MeasuresUtil.<V,E>getDisjointAgreements(new GraphDataSet<V, E>("dummy"))) {
			if(doMissingVlauesUntreated)
				outResultsE.write(( "," + clusterAgreement.toLatexString() ).getBytes());
			if(doMissingValueRemove)
				outResultsE.write(( ",Mre:" + clusterAgreement.toLatexString() ).getBytes());
			if( doMissingValueAsCluster)
				outResultsE.write(( ",MaC:" + clusterAgreement.toLatexString() ).getBytes());

		}
		if (doAlgebric)
			for(String name: AlgebricClusteringAgreement.getAllAgreements(null, null,null, null,doalpha, dotrans).first){
				outResultsE.write(( "," + name ).getBytes());
			}
		outResultsE.write(("\n").getBytes());
		
		//-------------------------------------------------------------------------
		while(datasets.hasNext()) {
			GraphDataSet<V, E> dataset = datasets.next();
			if (dataset==null) break;
			dataset.printStats();
			int n  = dataset.graph.getVertexCount();
			community_runtimes = new HashMap<>();
			
			//			if(att.equals("value")) groundth = attClustering;
			Vector<String> atts = new Vector<String>();
			if(dataset.attributes!=null && dataset.attributes.values().size()>0){
				for (Object att : dataset.attributes.values().iterator().next().keySet() ){
					String attName = att.toString().toLowerCase();
					if(attName.equals("id") || attName.equals("label")) continue;
					atts.add(attName);
				}
			}
			
			//Find All Communities
			if (doCommunities){
				//TOPLEADER
//				ClusteringAgreement<V> ARI = new GraphAGAM<V,E>(dataset.graph,dataset.getWeights(),Type.X2,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED);

				startTime = System.currentTimeMillis();
				int upper =(int) Math.sqrt(dataset.graph.getVertexCount());
				int lower = 2;
//				double agg=0, max_agg = 0, downHillCounter =0, aggPrev = 0;
//				double precision = 0.001;
//				int step = 1;
				Grouping<V> grouping =null; 
				for (int k=lower; k<=upper; k++){
					System.err.println(k+" ...  " );
//					aggPrev = agg;
					Grouping<V> res = (new TopLeaders(k)).findCommunities( dataset.graph, dataset.weights);
				
//					for(String att: atts) if ( (doAtt || att.equals("value") ) && !community_names.contains(att))
//						if (att.equals("year") || att.equals("dorm/house")){				
//							Grouping< V> groupingAtt = dataset.getGrouping(att,missingValueIndicator);
//							System.err.println("----------> Attribute: "+att + " with "  +groupingAtt.getNumberOfGroups());
//
//					if (groupingAtt.getNumberOfGroups()<upper) upper = groupingAtt.getNumberOfGroups();
					System.err.println( " topLeaders will find at most " + upper);
//						agg = ARI.getAgreement(res.getGroups(), groupingAtt.getGroups());
//						System.err.println(" ...  " + agg + " from" +ARI.toString());
//						if (agg>max_agg) {
//							grouping = res;
//							downHillCounter =0;
//						}
//						else {
//							if(aggPrev-agg > precision)
//								downHillCounter++;
//							if (downHillCounter>3) k = upper+1;
//						}
//					}
					String comName = "TopLeaders_"+k;
					community_names.add(comName.toLowerCase());
					community_runtimes.put(comName,System.currentTimeMillis() - startTime);
					System.err.println(">>  --  finished in: " +community_runtimes.get(comName) +" milisecond");
					if(grouping==null || grouping.getNumberOfGroups()<1){
						System.err.println("--------> ERROR :: "+comName+" failed on dataset "+ dataset.name + " resulted in " +	grouping);
						continue;
					}
					System.err.print(">>  -- resulted in " + grouping.getNumberOfGroups()+	" clusters " );
					dataset.addPartitioingAttribute(comName.toLowerCase(), grouping.getGroups());	
				}

				for (CommunityMiner<V, E> communityMiner : AlgorithmUtils.<V, E>getSelectedCommunititMiners(isOverlapping)) {
					System.err.println(">>  --  Method: " + communityMiner.getName());
					if (atts.contains( communityMiner.getName().toLowerCase())) continue;
	//				if(!community_names.contains(communityMiner.getName()))community_names.add(communityMiner.getName());
					
					startTime = System.currentTimeMillis();
					grouping =  communityMiner.findCommunities(dataset.graph, dataset.weights );
					
					community_names.add(communityMiner.getName().toLowerCase());
					community_runtimes.put(communityMiner.getName(),System.currentTimeMillis() - startTime);
					System.err.println(">>  --  finished in: " +community_runtimes.get(communityMiner.getName()) +" milisecond");
	
					if(grouping==null || grouping.getNumberOfGroups()<1){
						System.err.println("--------> ERROR :: "
								+communityMiner+" failed on dataset "+ dataset.name + " resulted in " +
								grouping);
						continue;
					}
					System.err.print(">>  -- resulted in " + grouping.getNumberOfGroups()+	" clusters " );
					dataset.addPartitioingAttribute(communityMiner.getName().toLowerCase(), grouping.getGroups());				
				}
				if (writeResultedCommunities){
					IOUtils.<V,E>writeGML(runname+"com/"+dataset.name+".gml",
							dataset.graph, null, dataset.weights, null, dataset.attributes);
				}
			}else {
				for (CommunityMiner<V, E> communityMiner : AlgorithmUtils.<V, E>getSelectedCommunititMiners(isOverlapping)) {
					if(!community_names.contains(communityMiner.getName()))
						community_names.add(communityMiner.getName().toLowerCase());
				}
			}
			//Do Evaluation


			//=================
			for(String att: atts) {
				System.err.println("----------> Attribute: "+att);
				Grouping< V> grouping =dataset.getGrouping(att,missingValueIndicator) ;
				System.err.println("In  "+ dataset.name+"  "+ att+" has "+ grouping.getNumberOfGroups()+
						" clusters "+grouping.getNumberOfGroupsOfSizeAtLeast(1)+" >0 "+
						grouping.getNumberOfGroupsOfSizeAtLeast(5)+" >4 ");
				
				outResultsR.write((dataset.name+( (dataset.weights== null)?",n":",y") +","+n+","+dataset.graph.getEdgeCount()+","+ 
				att+","+ (community_runtimes.containsKey(att)?community_runtimes.get(att):0)*0.001	+"," + grouping.getNumberOfGroups()).getBytes());
				if (doRelative){
					for ( RelativeCriteria<V> relative: MeasuresUtil.<V,E>getRelativeAlternatives(dataset.graph, dataset.getWeightsTransformer())){
						startTime = System.currentTimeMillis();
						System.out.print(relative.getName()+" for "+att);
						double val = relative.evaluate(grouping.getGroups());
						System.out.println(" : "+val + " (computed in "+(System.currentTimeMillis() - startTime)+" miliseconds )");
						outResultsR.write( ("," + val	).getBytes());
					}
				}
				outResultsR.write(("\n").getBytes());
			}  
			
			//=================
			for (String communityMiner : community_names) {
				Grouping<V> groupingComm = dataset.getGrouping(communityMiner);
				if (groupingComm.getNumberOfGroups()<=0) continue;
				System.err.println(communityMiner+" with "+groupingComm.getNumberOfGroups());
				for(String att: atts) if ( (doAtt || att.equals("value") ) && !community_names.contains(att)){
					Grouping< V> groupingAtt = dataset.getGrouping(att,missingValueIndicator);
					System.err.println("----------> Attribute: "+att + " with "  +groupingAtt.getNumberOfGroups());
						
					outResultsE.write((dataset.name+( (dataset.weights== null)?",n":",y") +","+n+","+dataset.graph.getEdgeCount()
								+","+ communityMiner	+ "," + groupingComm.getNumberOfGroups()+ (doAtt?"," + att:"")	).getBytes());
						startTime = System.currentTimeMillis();
						for (ClusteringAgreement<V> clusteringAgreement : MeasuresUtil.<V,E>getDisjointAgreements(dataset)) {
							double accu;
							//MissingValues handeled with measures directly
							if(doMissingVlauesUntreated){
								accu= clusteringAgreement.getAgreement(
										groupingComm.getGroups(),groupingAtt.getGroups());
								outResultsE.write(( "," + accu ).getBytes());
								if(((accu+"").equals("NaN")))System.err.println("ERROR:"+dataset.name+":"+clusteringAgreement+" = " + accu + " c: " + groupingComm.getStatistics() + " with g: "+groupingAtt.getStatistics());
							}
							//MissingValues removed
							if(doMissingValueRemove){
								accu= clusteringAgreement.getAgreementCovering(
										groupingAtt.getDataPoints(), groupingComm.getGroups(),groupingAtt.getGroups());
								outResultsE.write(( "," + accu ).getBytes());
								if(((accu+"").equals("NaN")))System.err.println("ERROR:"+dataset.name+":"+clusteringAgreement+" = " + accu + " c: " + groupingComm.getStatistics() + " with g: "+groupingAtt.getStatistics());
							}
							//MissingValues will be added as a cluster 
							if(doMissingValueAsCluster){
								accu= clusteringAgreement.getAgreementCovering(
										dataset.graph.getVertices(), groupingComm.getGroups(),groupingAtt.getGroups());
								outResultsE.write(( "," + accu ).getBytes());
								if(((accu+"").equals("NaN")))System.err.println("ERROR:"+dataset.name+":"+clusteringAgreement+" = " + accu + " c: " + groupingComm.getStatistics() + " with g: "+groupingAtt.getStatistics());
							}
						}						
						System.err.println(" agreements computed in "+(System.currentTimeMillis() - startTime)+" miliseconds )");
						
						if(doAlgebric){
							startTime = System.currentTimeMillis();
							Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(dataset.graph, dataset.getWeightsTransformer(), groupingComm.getGroups(), groupingAtt.getGroups(),doalpha, dotrans);
							for (int i = 0; i < tmp.first.size(); i++) {
								double accu= tmp.second.get(i);
								outResultsE.write(( "," + accu ).getBytes());
								if(((accu+"").equals("NaN")))
									System.err.println("ERROR:"+dataset.name+":"+tmp.first.get(i)+" = " + accu + " c: " + groupingComm.getStatistics() + " with g: "+groupingAtt.getStatistics());
	//								evalMetrics.put((tmp.first.get(i)+(att.equals("value")?"":("."+att))), df.format(accu));
							}
							System.err.println(" algebrics ---- computed in "+(System.currentTimeMillis() - startTime)+" miliseconds)");
						}
						outResultsE.write(("\n").getBytes());
					}
			}
		}
		outResultsE.close();
		outResultsR.close();

	}
		    
	public static boolean createCleanDir(File dir) {
		    if (dir.isDirectory()) {
		        String[] children = dir.list();
		        for (int i=0; i<children.length; i++) {
		            boolean success = createCleanDir(new File(dir, children[i]));
		            if (!success) {
		                return false;
		            }
		        }
		    }
		    return dir.mkdir();
	//		    return dir.delete();
	}
	  

	
	
	
	
	
	
	public <V,E> void __compareMethodsWithGTOld(Iterator<GraphDataSet<V,E>> datasets, String runname, boolean overlapping) throws IOException{
			
			FileOutputStream outResults = new FileOutputStream(runname+"res"+".csv");			
			outResults.write( ("dataset,weights,#nodes,#edges,CM_alg,q,t,k").getBytes());

			for (ClusteringAgreement<V> clusterAgreement :  MeasuresUtil.<V,E>getAgreementAlternatives(null, null, overlapping)) {
				outResults.write(( "," + clusterAgreement.toLatexString() ).getBytes());
			}
			for(String name: AlgebricClusteringAgreement.getAllAgreements(null, null,null, null).first){
				outResults.write(( "," + name ).getBytes());
			}

			outResults.write(("\n").getBytes());
			
//			System.err.println("Comparing methods on "+datasets.size()+"   datasets ");
//			for (GraphDataSet<V, E> dataset : datasets) {
			while(datasets.hasNext()) {
				GraphDataSet<V, E> dataset = datasets.next();
				dataset.printStats();
//				System.out.print("------>  Graph "+dataset.name+" : ");
//				System.out.print(dataset.graph.getVertexCount()+" nodes and " + dataset.graph.getEdgeCount() + " Edges");
				
				Modularity<V, E> modularity = new Modularity<V, E>(dataset.graph,dataset.getWeightsTransformer());
				Vector<CommunityMiner<V, E>> communityMiners =  AlgorithmUtils.<V, E>getSelectedCommunititMiners(overlapping);
				initializePlot(dataset.name, communityMiners.size());
				for (CommunityMiner<V, E> communityMiner : communityMiners) {
					Map<String , String> evalMetrics = new  HashMap<String , String>(); 
					System.err.println(">>  --  Method: " + communityMiner.getName());
				
					long startTime = System.currentTimeMillis();
					Grouping<V> grouping =  communityMiner.findCommunities(dataset.graph, dataset.weights );
					long duration = (System.currentTimeMillis() - startTime); //from nano to milisecond
					System.err.println(">>  --  finished in: " + duration +" milisecond");

					if (grouping == null) continue;
					dataset.addPartitioingAttribute(communityMiner.getName(), grouping.getGroups());				
					IOUtils.<V,E>writeGML(runname+"com/"+dataset.name+"."+communityMiner+".gml",
							dataset.graph, null, dataset.weights, null, dataset.attributes);
					
					System.err.print(">>  -- resulted in " + grouping.getNumberOfGroups()+
							" clusters with Q = " );
					startTime = System.currentTimeMillis();
					if(grouping==null || grouping.getNumberOfGroups()<1){
//						System.in.read();
						System.err.println("erooor :: "+communityMiner+" failed on dataset "+ dataset.name);
					}
					double q = modularity.evaluate(grouping.getGroups());
					System.out.println(" : "+q + 
							" (Q computed in "+(System.currentTimeMillis() - startTime)+" miliseconds )");

					int n  = dataset.graph.getVertexCount();
					outResults.write((dataset.name+( (dataset.weights== null)?",n":",y") +","+n+","+dataset.graph.getEdgeCount()
							+","+ communityMiner.getName()+
							"," + df.format(q) + "," + duration*0.001 + "," +
								+ grouping.getNumberOfGroups()).getBytes());
					
					evalMetrics.put("time", formatter.format(new Date(duration)));
					evalMetrics.put("Q", df.format(q));
					
					Vector<String> atts = new Vector<String>();
					if(dataset.attributes!=null && dataset.attributes.values().size()>0){
						for (Object att : dataset.attributes.values().iterator().next().keySet() ){
							String attName = att.toString().toLowerCase();
							if(attName.equals("id") || attName.equals("label")) continue;
							atts.add(attName);
						}
					}
					Map<V,String> labels = null;
					Grouping<V> groundth = null;
					for(String att: atts) if (att.equals("value")||att.equals("label")){
						Grouping< V> attClustering =null ;
						if(att.equals("label")) labels = dataset.getAttMap("label");
						else {
							attClustering = dataset.getGrouping(att, new Transformer<Object, Integer>() {
							Vector<Object> clusters = new Vector<>();
							@Override
							public Integer transform(Object obj) {
								if (!clusters.contains(obj)) clusters.add(obj);
								return clusters.indexOf(obj);
							}
							});
						}
						if(att.equals("value")) groundth = attClustering;
						startTime = System.currentTimeMillis();

						Vector< ClusteringAgreement<V>> res = MeasuresUtil.<V,E>getAgreementAlternatives(dataset.graph, dataset.getWeightsTransformer(), overlapping);
						for (ClusteringAgreement<V> clusteringAgreement : res) {
							double accu=0;
							if (att.equals("value")){
								accu = clusteringAgreement.getAgreementCovering(dataset.graph.getVertices(), grouping.getGroups(),attClustering.getGroups());
								outResults.write(( "," + accu ).getBytes());
								
								if(((accu+"").equals("NaN")))
								System.err.println(dataset.name+":"+clusteringAgreement+" = " + accu 
										+ " c: " + grouping.getStatistics() 
										+ " with g: "+attClustering.getStatistics()
										);
							}
//							System.err.print("-----"+att+":"+accu+":"+"."+clusteringAgreement);
//							if (att.equals("value"))

							evalMetrics.put((clusteringAgreement.toString()+(att.equals("value")?"":("."+att))), df.format(accu));
						}
						System.err.println(" classics ---- computed in "+(System.currentTimeMillis() - startTime)+" miliseconds )");
						startTime = System.currentTimeMillis();

						Pair<Vector<String>, Vector<Double>> tmp = AlgebricClusteringAgreement.getAllAgreements(dataset.graph, dataset.getWeightsTransformer(), grouping.getGroups(), attClustering.getGroups());
						for (int i = 0; i < tmp.first.size(); i++) {
							double accu=0;
							if (att.equals("value")){
								accu = tmp.second.get(i);
								outResults.write(( "," + accu ).getBytes());
								
								if(((accu+"").equals("NaN")))
								System.err.println(dataset.name+":"+tmp.first.get(i)+" = " + accu 
										+ " c: " + grouping.getStatistics() 
										+ " with g: "+attClustering.getStatistics()
										);
							}
							evalMetrics.put((tmp.first.get(i)+(att.equals("value")?"":("."+att))), df.format(accu));
				
						}
						System.err.println(" algebrics ---- computed in "+(System.currentTimeMillis() - startTime)+" miliseconds)");

					}
					outResults.write(("\n").getBytes());
					
					evalMetrics.put("k",  grouping.getNumberOfGroups()+"");
				
					
					plotCommunities(runname, dataset.name //.substring(0,dataset.name.indexOf('.'))
							+" "+ communityMiner.getName() + " k="+ evalMetrics.get("k") + 
							" Q="+ evalMetrics.get("Q") + " NMI="+ evalMetrics.get("NMI") // toString()
							, dataset.graph,dataset.weights, labels, groundth, grouping);
				}
			}
			outResults.close();
		}
		  
	  public <V,E> void __compareMethodsSynth(Vector<GraphDataSet<V, E>> datasets, String runname) throws IOException{
			FileOutputStream outResults = new FileOutputStream(runname+"res"+".csv");			
			outResults.write( ("dataset,weights,#nodes,#edges,CM_alg,q,t,k").getBytes());

			for (ClusteringAgreement<V> clusteringAgreement :  MeasuresUtil.<V,E>getAgreementAlternatives(null, null)) {
				outResults.write(( "," + clusteringAgreement.toLatexString() ).getBytes());
			}
			outResults.write(("\n").getBytes());
			
			System.err.println("Comparing methods on "+datasets.size()+"   datasets ");
			for (GraphDataSet<V, E> dataset : datasets) {
				dataset.printStats();
//				System.out.print("------>  Graph "+dataset.name+" : ");
//				System.out.print(dataset.graph.getVertexCount()+" nodes and " + dataset.graph.getEdgeCount() + " Edges");
				
				Modularity<V, E> modularity = new Modularity<V, E>(dataset.graph,dataset.getWeightsTransformer());

				
				Vector<String> atts = new Vector<String>();
				if(dataset.attributes!=null && dataset.attributes.values().size()>0){
					for (Object att : dataset.attributes.values().iterator().next().keySet() ){
						String attName = att.toString().toLowerCase();
						if(attName.equals("id") || attName.equals("label")) continue;
						atts.add(attName);
					}
				}
				System.err.println(atts);
				Map<V,String> labels = null;
				Grouping<V> groundth = null;
				for(String att: atts){
					Grouping< V> attClustering = dataset.getGrouping(att, new Transformer<Object, Integer>() {
						Vector<Object> clusters = new Vector<>();
						@Override
						public Integer transform(Object obj) {
							if (!clusters.contains(obj)) clusters.add(obj);
							return clusters.indexOf(obj);
						}
					});
					if(att.equals("label")) labels = dataset.getAttMap("label");
					if(att.equals("value")) groundth = attClustering;
				}
				
				Vector<CommunityMiner<V, E>> communityMiners =  AlgorithmUtils.<V, E>getSelectedCommunititMiners();
//				initializePlot(dataset.name, communityMiners.size());
				
				for (CommunityMiner<V, E> communityMiner : communityMiners) {
//					Map<String , String> evalMetrics = new  HashMap<String , String>(); 

					System.err.println("--Method: " + communityMiner.getName());
				
					long startTime = System.nanoTime();
					Grouping<V> grouping =  communityMiner.findCommunities(dataset.graph, dataset.weights );
					long endTime = System.nanoTime();
					long duration = (endTime - startTime)/1000000; //from nano to milisecond

					if (grouping == null) continue;
//					GMLGraphWriter<V, E> gmlGraphWriter = new GMLGraphWriter<V, E>();
//					gmlGraphWriter.writeGraph(runname+"com/"+dataset.name+"."+communityMiner+".gml",
//							dataset.graph, null, dataset.weights, null, dataset.attributes);
//					dataset.addPartitioingAttribute(communityMiner.getName(), grouping.getGroups());				
					
					System.err.print("--done! in " + duration+ "   q = ");
					double q = modularity.evaluate(grouping.getGroups());
					System.out.println("--"+q);

					System.err.println("mu = "+dataset.name.substring(1,2));
					int n  = dataset.graph.getVertexCount();
					outResults.write((dataset.name+( (dataset.weights== null)?",n":",y") +","+n+","+dataset.graph.getEdgeCount()
							+","+ communityMiner.getName()+
							"," + df.format(q) + "," + duration*0.001 + "," +
							dataset.name.substring(1,2) //mu
								//+ grouping.getNumberOfGroups()
								).getBytes());
					
//					evalMetrics.put("time", formatter.format(new Date(duration)));
//					evalMetrics.put("Q", df.format(q));
					
//					Vector<String> atts = new Vector<String>();
//					if(dataset.attributes!=null && dataset.attributes.values().size()>0){
//						for (Object att : dataset.attributes.values().iterator().next().keySet() ){
//							String attName = att.toString().toLowerCase();
//							if(attName.equals("id") || attName.equals("label")) continue;
//							atts.add(attName);
//						}
//					}
//					System.err.println(atts);
//					Map<V,String> labels = null;
//					Grouping<V> groundth = null;
//					for(String att: atts){
//						Grouping< V> attClustering = dataset.getGrouping(att, new Transformer<Object, Integer>() {
//							Vector<Object> clusters = new Vector<>();
//							@Override
//							public Integer transform(Object obj) {
//								if (!clusters.contains(obj)) clusters.add(obj);
//								return clusters.indexOf(obj);
//							}
//						});
//						if(att.equals("label")) labels = dataset.getAttMap("label");
//						if(att.equals("value")) groundth = attClustering;

						Vector< ClusteringAgreement<V>> res = MeasuresUtil.getAgreementAlternatives(dataset.graph, dataset.getWeightsTransformer());
						for (ClusteringAgreement<V> clusteringAgreement : res) {
							System.err.println("-----"+clusteringAgreement.toLatexString());
//							System.err.println(attClustering.getGroups());
							double accu=0;
//							if (att.equals("value")){
								accu = clusteringAgreement.getAgreementCovering(dataset.graph.getVertices(), grouping.getGroups(),groundth.getGroups());
//							}
//							System.err.print("-----"+att+":"+accu+":"+"."+clusteringAgreement);
//							if (att.equals("value"))
								outResults.write(( "," + accu ).getBytes());

//							evalMetrics.put((clusteringAgreement.toString()+(att.equals("value")?"":("."+att))), df.format(accu));
						}
//
//						
//					}
					outResults.write(("\n").getBytes());
					
//					evalMetrics.put("k",  grouping.getNumberOfGroups()+"");
				
					
//					plotCommunities(runname, dataset.name //.substring(0,dataset.name.indexOf('.'))
//							+" "+ communityMiner.getName() + " k="+ evalMetrics.get("k") + 
//							" Q="+ evalMetrics.get("Q") + " NMI="+ evalMetrics.get("NMI") // toString()
//							, dataset.graph,dataset.weights, labels, groundth, grouping);
				}
			}
			outResults.close();
		}
	  
	  public <V,E> void __compareExternals(Vector<GraphDataSet<V, E>> datasets, String runname) throws IOException{
			FileOutputStream outResults = new FileOutputStream(runname+"resExt"+".csv");			
			outResults.write( ("dataset,weights,#nodes,#edges,CM_alg,q,t,k").getBytes());

			for (ClusteringAgreement clusteringAgreement :  MeasuresUtil.getAgreementAlternatives(null, null)) {
				outResults.write(( "," + clusteringAgreement.toLatexString() ).getBytes());
			}
			outResults.write(("\n").getBytes());
			
			System.err.println("Comparing methods on "+datasets.size()+"   datasets ");
			for (GraphDataSet<V, E> dataset : datasets) {
				dataset.printStats();
				
				Vector<String> atts = new Vector<String>();
				if(dataset.attributes!=null && dataset.attributes.values().size()>0){
					for (Object att : dataset.attributes.values().iterator().next().keySet() ){
						String attName = att.toString().toLowerCase();
						if(attName.equals("id") || attName.equals("label")) continue;
						atts.add(attName);
					}
				}
				Map<V,String> labels = null;
				Grouping<V> groundth = null;
				for(String att: atts){
					Grouping< V> attClustering = dataset.getGrouping(att, new Transformer<Object, Integer>() {
						Vector<Object> clusters = new Vector<>();
						@Override
						public Integer transform(Object obj) {
							if (!clusters.contains(obj)) clusters.add(obj);
							return clusters.indexOf(obj);
						}
					});
					if(att.equals("label")) labels = dataset.getAttMap("label");
					if(att.equals("value")) groundth = attClustering;
				}
				
				Modularity<V, E> modularity = new Modularity<V, E>(dataset.graph,dataset.getWeightsTransformer());

				ExternalIndexComparer<V, E> comparer = new ExternalIndexComparer<>(Mode.RAND2G);//FRAGKNEE);//RAND2G);
				Vector<Vector<Set<V>>> randomClusters = comparer.getParts(dataset.graph, groundth.getGroups());
				for (Vector<Set<V>> grouping : randomClusters) {
					Map<String , String> evalMetrics = new  HashMap<String , String>(); 

					if (grouping == null) continue;
					double q = modularity.evaluate(grouping);
//					System.out.print("--"+q);

					int n  = dataset.graph.getVertexCount();
					outResults.write((dataset.name+( (dataset.weights== null)?",n":",y") +","+n+","+dataset.graph.getEdgeCount()
							+","+ comparer.mode +
							"," + df.format(q) + "," + 0*0.001 + "," +
								+ grouping.size()).getBytes());
					

					Vector< ClusteringAgreement<V>> res = MeasuresUtil.getAgreementAlternatives(dataset.graph, dataset.getWeightsTransformer());
					for (ClusteringAgreement<V> clusteringAgreement : res) {
//							System.err.println("-----"+clusteringAgreement.toLatexString());
//							System.err.println(attClustering.getGroups());
						double accu= clusteringAgreement.getAgreementCovering(dataset.graph.getVertices(), grouping, groundth.getGroups());
							outResults.write(( "," + accu ).getBytes());
					}
					outResults.write(("\n").getBytes());
				}					
			}
			outResults.close();
		}
	  
	boolean plotGraphs = false;
	JFrame jf;
	private void initializePlot(String title, int n){
			if (!plotGraphs) return;
			jf = new JFrame(title);
		//	jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			if(n>4)
		    jf.setLayout( new GridLayout( n/2 +n%2,2)); 
			else 		    jf.setLayout( new GridLayout( n,1)); 

		    jf.setBackground(Color.white);
	}
	public <V,E> void plotCommunities(String runname,String title,Graph<V,E> graph, 	
			Map<E, Double> edgeWeights ,Map<V,String> labels ,Grouping<V> colourGroup, Grouping<V> shapeGroup){
		if (!plotGraphs) return;

		if (graph.getVertexCount()>500) return;
		PartitionView< V, E> partitionView = new PartitionView< V, E>(graph,
				edgeWeights,labels,
				colourGroup,
				shapeGroup);
		partitionView.setSize(200, 200);
		JSplitPane pane = new JSplitPane(0, new JLabel(title,JLabel.CENTER) ,partitionView);
		pane.setBackground(Color.white);
		pane.setDividerSize(0);
		jf.getContentPane().add(pane);
		jf.pack();
		jf.setVisible(true);
		partitionView.saveAsFigure(runname +"/img/"+title);
		
		//Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		//jf.setLocation((int)screenSize.getWidth()/6, (int)screenSize.getHeight()/6);
		//jf.setSize((int)screenSize.getWidth()*2/3,(int)screenSize.getHeight()*2/3);
	}
	

			
		
}
	


