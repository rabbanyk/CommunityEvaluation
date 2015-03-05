package dev_Experiments;

import java.awt.Color;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.sql.ClientInfoStatus;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JSplitPane;

import org.apache.commons.collections15.Transformer;

import measure.MeasuresUtil;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.graph.criteria.Modularity;
import ui.PartitionView;
import util.DatasetUtils;
import util.DatasetUtils.ClassicDataset;
import io.graph.gml.GMLGraphWriter;
import algorithms.AlgorithmUtils;
import algorithms.AlgorithmUtils.Method;
import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.data.Partitioning;
import data.GraphDataSet;
import edu.uci.ics.jung.algorithms.metrics.Metrics;
import edu.uci.ics.jung.graph.Graph;

public class ExternalExperiment {
	
	  DecimalFormat df = new DecimalFormat("#.###");
	  DateFormat formatter = new SimpleDateFormat("mm:ss:SSS");
	  
	public <V,E> void compareMethods(Vector<GraphDataSet<V, E>> datasets, String runname) throws IOException{
//		communityMiners.addAll(minerUtils.getAllCommunititMiners());
		//communityMiners.add(new communityMinerIGraph<V,E>());																																																																																							
//		communityMiners.addAll( minerUtils.getVariant(Method.BP));

//		FileOutputStream out = new FileOutputStream("comparRes"+runName+".csv");
		FileOutputStream outTime = new FileOutputStream(runname+"res"+".csv");
					
		System.err.println("Comparing methods on "+datasets.size()+"   datasets ");
	
		for (GraphDataSet<V, E> dataset : datasets) {
//			dataset.printStats();
			System.out.print("------>  Graph "+dataset.name+" : ");
			System.out.print(dataset.graph.getVertexCount()+" nodes and " + dataset.graph.getEdgeCount() + " Edges");
			
			Vector<CommunityMiner<V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
//			communityMiners.add(minerUtils.getCommunititMiner(Method.FastModularity));
//			communityMiners.addAll(minerUtils.getModularityBasedCommunititMiners());

			communityMiners.add( AlgorithmUtils.<V, E>getCommunititMiner(Method.Louvain));
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.fastgreedy));		
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.spinglass));		
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.leading_eigenvector));		
//
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.infomap));		
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.walktrap));
//			communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.label_propagation));

			
			initializePlot(dataset.name, communityMiners.size());
			
			Modularity<V, E> modularity = new Modularity<V, E>(dataset.graph,dataset.getWeights());
			for (CommunityMiner<V, E> communityMiner : communityMiners) {
				Map<String , String> evalMetrics = new  HashMap<String , String>(); 

				System.err.println("Method: " + communityMiner.getName());
				long startTime = System.nanoTime();
				Grouping<V> grouping =  communityMiner.findCommunities(dataset.graph, dataset.weights );
				System.err.println(grouping);
			
				long endTime = System.nanoTime();
				long duration = (endTime - startTime)/1000000; //from nano to milisecond
				
//				FileOutputStream outComm = new FileOutputStream("res/"+dataset.name+"."+communityMiner+"."+runname+".com");
//				outComm.write(grouping.getGroups().toString().getBytes());
//				outComm.close();

				dataset.addPartitioingAttribute(communityMiner.getName(), grouping.getGroups());				
				GMLGraphWriter<V, E> gmlGraphWriter = new GMLGraphWriter<V, E>();
				gmlGraphWriter.writeGraph("res/"+dataset.name+"."+communityMiner+"."+runname+".gml",
						dataset.graph, null, dataset.weights, null, dataset.attributes);
				
				System.err.print("done! in " + duration+ "   q = ");
				double q = modularity.evaluate(grouping.getGroups());
				System.out.println(q);
				//w,n, e, m , q, t 
				int n  = dataset.graph.getVertexCount();
				
				outTime.write((dataset.name+( (dataset.weights== null)?" , n":" , y") +" , "+n+" , "+dataset.graph.getEdgeCount()
						+" , "+ communityMiner.getName()+
						"  ,  " + df.format(q) + " , " + duration*0.001 + "  ,  " +
						" , "+ grouping.getNumberOfGroups()).getBytes());
				
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
//				String [] atts = {"value"};
//				String [] atts ={"label","major",
//			    "second major/minor (if applicable)",
//			    "dorm/house","gender","a student/faculty status flag" ,"high school", "year"};
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

					Vector< ClusteringAgreement<V>> res = MeasuresUtil.<V,E>getAgreementAlternatives(dataset.graph, dataset.getWeights());
					for (ClusteringAgreement<V> partiotioningAgreement : res) {
						System.err.println(partiotioningAgreement);
						double accu = partiotioningAgreement.getAgreement(grouping.getGroups(), attClustering.getGroups());
						
						System.err.println(att+":"+accu+":"+"."+partiotioningAgreement);
						outTime.write(( "  ,  " + accu ).getBytes());

						evalMetrics.put((partiotioningAgreement.toString()+(att.equals("value")?"":("."+att))), df.format(accu));
					}

					
				}
				outTime.write(("\n").getBytes());
				
				evalMetrics.put("k",  grouping.getNumberOfGroups()+"");
			
				plotCommunities(dataset.name //.substring(0,dataset.name.indexOf('.'))
						+" "+ communityMiner.getName() + " k="+ evalMetrics.get("k") + 
						" Q="+ evalMetrics.get("Q") + " NMI="+ evalMetrics.get("NMI") // toString()
						, dataset.graph,dataset.weights, labels, groundth, grouping);
			}
		}
//		out.close();
		outTime.close();
	}
	
	JFrame jf;
	private void initializePlot(String title, int n){
			jf = new JFrame(title);
		//	jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			if(n>4)
		    jf.setLayout( new GridLayout( n/2 +n%2,2)); 
			else 		    jf.setLayout( new GridLayout( n,1)); 

		    jf.setBackground(Color.white);
	}
	
	
	public <V,E> void plotCommunities(String title,Graph<V,E> graph, 	
			Map<E, Double> edgeWeights ,Map<V,String> labels ,Grouping<V> colourGroup, Grouping<V> shapeGroup){
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
		partitionView.saveAsFigure("./img/"+title);
		
		//Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		//jf.setLocation((int)screenSize.getWidth()/6, (int)screenSize.getHeight()/6);
		//jf.setSize((int)screenSize.getWidth()*2/3,(int)screenSize.getHeight()*2/3);
	}
	
	
	public static void main(String[] args){
		new File("res").mkdir();
		new File("img").mkdir();

		Vector<GraphDataSet<Integer, Integer>> datasets = new Vector<GraphDataSet<Integer,Integer>>();
// ---------------  TEST -------------------------
//		datasets.add(DatasetUtils.<Integer,Integer>loadClassic(ClassicDataset.wkarate));

// ---------------  EXP 1 -------------------------
//		for (GraphDataSet<Integer, Integer> dataset: DatasetUtils.<Integer,Integer>loadAll("../Datasets/classics/")){
//			if (dataset.graph.getEdgeCount()<20000){
//				dataset.printStats();
//				datasets.add(dataset);
//			}
//		}
// ---------------  EXP 2 -------------------------
		for (GraphDataSet<Integer, Integer> dataset: DatasetUtils.<Integer,Integer>loadAllDataSets("../Datasets/SNAP/facebook/")){
			if (dataset.graph.getEdgeCount()<20000){
				dataset.printStats();
				datasets.add(dataset);
			}
		}
//		
//		for (GraphDataSet<Integer, Integer> dataset: DatasetUtils.<Integer,Integer>loadAll("../Datasets/Biological/linkGroup/")){
//				{dataset.printStats();
//				datasets.add(dataset);}
//		}
		
//		for (GraphDataSet<Integer, Integer> dataset: DatasetUtils.<Integer,Integer>loadAll("src/models/generative/iGraph/")){
//			if (dataset.graph.getEdgeCount()<20000){
//				dataset.printStats();
//				datasets.add(dataset);
//			}
//		}
//		datasets.addAll(DatasetUtils.<Integer,Integer>loadAll("../Datasets/classics/power/")) ;
//		datasets.add(DatasetUtils.<Integer,Integer>load("ff.gml")) ;

		
		//	datasets.addAll(DatasetUtils.loadAllClassics());
//		datasets.add(DatasetUtils.<Integer,Integer>load(new File("Rice31.gml")));/home/reihaneh/projects/Datasets/SNAP/facebook
		//datasets.addAll(ClassicDatasetLoader.loadLFR(false));
	//	datasets.add(ClassicDatasetLoader.loadFacebook(false));
		
		ExternalExperiment test = new ExternalExperiment();
		try {
			test.compareMethods(datasets, "ForestFireTestTest");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
			
			
		
		
	}
	

