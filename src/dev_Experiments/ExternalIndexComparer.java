package dev_Experiments;

import io.Logger;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Scanner;
import java.util.Set;
import java.util.Vector;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.renderer.xy.DeviationRenderer;
import org.jfree.data.xy.YIntervalSeries;
import org.jfree.data.xy.YIntervalSeriesCollection;

import measure.cluster.agreement.*;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import edu.uci.ics.jung.graph.Graph;
import io.Logger.DebugMode;

public class ExternalIndexComparer<V,E> extends RandomExpreriment<V, E>{
	public ExternalIndexComparer(String experimentPath) {
		super( experimentPath);
	}
	public ExternalIndexComparer(Mode mode) {
		this.mode = mode;
	}
	static DateFormat formatter = new SimpleDateFormat("MM:dd-HH:mm:ss");

	public Vector<Vector<Set<V>>> generateDifferentVariationFromGroundTruth(Graph<V,E> graph, Vector<Set<V>> groundTruth){
		if(! (mode == Mode.FRAGKNEE)) return super.generateDifferentVariationFromGroundTruth(graph, groundTruth);
//		System.err.println("HERE");
		Vector<Vector<Set<V>>> mergeResults = new Vector<Vector<Set<V>>>();
		Vector<Vector<Set<V>>> splitResults = new Vector<Vector<Set<V>>>();

		mergeResults.add(groundTruth);
		splitResults.add(groundTruth);
		
		//maxResult/=2;
		while(mergeResults.size()+splitResults.size() < maxResult) {
			Vector<Vector<Set<V>>> splitVars = new Vector<Vector<Set<V>>>();
			Vector<Vector<Set<V>>> mergeVars = new Vector<Vector<Set<V>>>();

			for (Vector<Set<V>> vector : splitResults) {
				Vector<Set<V>> tmp =  generateSplittedVariations(vector,splitChance);
				if(tmp!=null && tmp.size()>=2 && validate(graph, vector)) splitVars.add(tmp);
				if(mergeResults.size()+splitResults.size() + splitVars.size() > maxResult) break;
			}
//			System.err.println("   "+vars.size()+":"+results.size());
			for (int i = 0; i < 100; i++) 
			for (Vector<Set<V>> vector : mergeResults) {
				Vector<Set<V>> tmp =   generateMergedVariations(vector,mergeChance);
				if(tmp!=null && tmp.size()>=2 && validate(graph, vector)) mergeVars.add(tmp);

				if(mergeResults.size()+splitResults.size() + mergeVars.size() > maxResult) break;
			}
//			System.err.println("        "+vars.size()+":"+results.size());

			for (Vector<Set<V>> vector : mergeVars) {
				if(mergeResults.size()+splitResults.size() < maxResult && validate(graph, vector))
					mergeResults.add( vector);
				else break;
			}
			for (Vector<Set<V>> vector : splitVars) {
				if(mergeResults.size()+splitResults.size() < maxResult && validate(graph, vector))
					splitResults.add( vector);
				else break;
			}
		}
		mergeResults.addAll(splitResults);
////		System.err.println(" Randomize");

		for (int i = 0; i < mergeResults.size(); i++) {
			// if(validate(graph, mergeResults.get(i)))
			Vector<Set<V>> randomized = moveNodesArround(mergeResults.get(i), swapChance);
			if(validate(graph, randomized)) mergeResults.set(i, randomized );
		}
	
	//	maxResult*=2;
//		while (mergeResults.size() < maxResult) {
//			Vector<Vector<Set<V>>> vars = new Vector<Vector<Set<V>>>();
//			for (Vector<Set<V>> vector : mergeResults){
//				 if(validate(graph, vector)) vars.add(moveNodesArround(vector, swapChance));
//				if (mergeResults.size() + vars.size() > maxResult)
//					break;
//			}
//
//			for (Vector<Set<V>> vector : vars) {
//				if (mergeResults.size() < maxResult)
//					 if(validate(graph, vector)) mergeResults.add(vector);
//				else
//					break;
//			}
//		}
		return mergeResults;
	}
	
	
	public Vector<Vector<Set<V>>> getParts (Graph<V,E> graph , Vector<Set<V>> groundTruth){
		if(MAXK>graph.getVertexCount())  
			MAXK =graph.getVertexCount();///3;
		Vector<Vector<Set<V>>> partitionings = new Vector<Vector<Set<V>>>();
		
		Vector<Vector<Set<V>>> parsToBeAdded = new Vector<Vector<Set<V>>>();
		Vector<Double> nk = new Vector<Double>(); // k 
		int n = 0;
		for(int run =0; run<MaxRun; run++) { //repeat 10 times
			if((parsToBeAdded.size()>MAX) ||
					((mode!=Mode.KNEE) && (parsToBeAdded.size() >= BUCK*(MAXK-MINK)))) break;
			Vector<Vector<Set<V>>> tmp = (mode == Mode.KNEE||mode == Mode.FRAGKNEE)?
					generateDifferentPartitionings(graph,groundTruth):
					generateRandomClusters(graph.getVertices(), MINK, MAXK, NR);

			System.err.print(tmp.size());
			Vector<Vector<Set<V>>> uniquePa = new Vector<Vector<Set<V>>>();
			for (int i = 0; i < tmp.size(); i++) {
				int k = tmp.get(i).size();
				if(parsToBeAdded.indexOf(tmp.get(i))==-1){ 
					for(int kk = nk.size(); kk <= k; kk++) nk.add(new Double(0));
					if(nk.get(k) < BUCK && k<MAXK){
						uniquePa.add(tmp.get(i));
						nk.set(k,  nk.get(k).doubleValue() + 1);
						n++;
					}
				}
			}
			parsToBeAdded.addAll(uniquePa);
			System.err.print("->"+uniquePa.size()+":"+parsToBeAdded.size()+" ");
			if(run%50==0)System.err.println();
		}
		
//		partitionings.add(new Vector<Vector<Set<V>>>());
		for (Vector<Set<V>> randPart : parsToBeAdded) {
			if((nk.get(randPart.size()) == BUCK))//||(nk.get(randPart.size())>=1 && mode == Mode.FRAGKNEE))
				partitionings.add(randPart);
		}
		
		System.err.println("\n"+partitionings.size() + " partitions generated ");//+ dataset.get(d).name );
		System.err.println(nk);
	
		return partitionings;
	}
	
	public void compareExternalIndexes(){
		//Preparing Data set
		String  DATASET = "/data/";
	
		Vector<NetworkWithGroundTruth> dataset =  loadNetworks(expPath + DATASET);
		
		Vector<Vector<Vector<Set<V>>>> partitionings = new Vector<Vector<Vector<Set<V>>>>(); // d x r x Vector<Set<V>>
		for (int d = 0; d< dataset.size();d++) {
			partitionings.add(getParts( dataset.get(d).graph,  dataset.get(d).groundTruth));
		}
		
		Vector<Vector<Vector<Number>>> extRess = externalEvaluation(dataset, partitionings); //d x am x p
		for (int rm = 0; rm < agreementMethods.size(); rm++) {
			logDatasetStat(dataset, partitionings , extRess,rm);
		}
		plotAllExternalFigures(dataset,partitionings,extRess);		

	}
	
//	public Vector<Vector<Vector<Number>>>  externalEvaluation(Vector<NetworkWithGroundTruth> dataset , Vector<Vector<Vector<Set<V>>>> partitionings){
//		//----------------------------------External Evaluation of the results----------------------------
//		Logger.logFunction("Computing External Agreements " );
//		Vector<Vector<Vector<Number>>> extRess = new Vector<Vector<Vector<Number>>>();
//		for (int d = 0; d< dataset.size(); d++) {
//			System.err.println(dataset.get(d).name+"--------------");
//			agreementMethods = getExternalAlternatives(dataset.get(d).graph, dataset.get(d).weights);
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
	
	public Vector<Number> externalEvaluation( PartiotioningAgreement<V> agreementMethod,Vector<Set<V>> groundTruth, Vector<Vector<Set<V>>> results){
//		System.err.println("GGGGGGGGGGGGGGGGGG "+groundTruth);
//		System.err.println("RRRRRRRRRRRRRRRRRR "+results);
		if(! (mode==Mode.RAND2RAND)) return super.externalEvaluation(agreementMethod, groundTruth, results);
		Vector<Number> res = new Vector<Number>();
		
		for (Vector<Set<V>> cluster : results) {
			Number agreement = agreementMethod.getAgreement(cluster, results.get(random.nextInt(results.size())));
			res.add(agreement);
		}
		return res;
	}
//	
    void plotAllExternalFigures(Vector<NetworkWithGroundTruth> dataset,  Vector<Vector<Vector<Set<V>>>> partitionings, Vector<Vector<Vector<Number>>> extRess){
		if(!plot) return;
		for (int d = 0; d< dataset.size(); d++) {
				//System.err.println(extRess.get(d));
				plot(dataset.get(d), extRess.get(d), partitionings.get(d) );
		}
	}
	
	boolean plot = true;
	
	void plot( NetworkWithGroundTruth dataset, Vector<Vector<Number>> extres, final Vector<Vector<Set<V>>> partitionings){
		System.err.println(dataset.name+"_ Comparing "+ agreementMethods.size()+"_ agreement methods on " +partitionings.size() + " partitionings");
		
		Vector<Double> nk = new Vector<Double>(); // k 
		Vector<Vector<Double>> s = new Vector<Vector<Double>>(); //  am x k
		Vector<Vector<Double>> s2 = new Vector<Vector<Double>>(); //
		
		int maxK = 0;
		for (int i = 0; i < partitionings.size(); i++) {
			int k = partitionings.get(i).size();
			if(k>maxK) maxK = k;
		}
		
		for (int i = 0; i <= maxK; i++) {
			nk.add(new Double(0));
		}
		
		for (int am = 0; am < agreementMethods.size(); am++) {
			s.add(new Vector<Double>());
			s2.add(new Vector<Double>());
			for (int i = 0; i <= maxK; i++) {
				s.get(am).add(new Double(0));
				s2.get(am).add(new Double(0));
			}
		}		
		for (int i = 0; i < partitionings.size(); i++) {
			int k = partitionings.get(i).size();
			nk.set(k,  nk.get(k).doubleValue() + 1);
			
			for (int am = 0; am < agreementMethods.size(); am++) {
				double er =  extres.get(am).get(i).doubleValue();
			
				s.get(am).set(k, s.get(am).get(k).doubleValue() + er);
				s2.get(am).set(k, s2.get(am).get(k).doubleValue() +  er*er);
			}
		}
		System.err.println(nk);
		
//		Vector<XYSeries> amSeries =  new Vector<XYSeries>(); 
		Vector<YIntervalSeries> amSeriesVar =  new Vector<YIntervalSeries>(); 

		for (int am = 0; am < agreementMethods.size(); am++) {
//			amSeries.add(new XYSeries(agreementMethods.get(am).toString()));
			amSeriesVar.add(new YIntervalSeries(agreementMethods.get(am).toString()));
		}
		
		
		for (int am = 0; am < agreementMethods.size(); am++) {
			for (int i = 2; i < s.get(am).size(); i++) {
				double avg , var;
				avg = s.get(am).get(i)/nk.get(i);
//				amSeries.get(am).add(i, avg );
				var = (s2.get(am).get(i)/nk.get(i) > avg*avg)?Math.sqrt(s2.get(am).get(i)/nk.get(i) - avg*avg ):0;
				//System.err.println(i+" "+avg+"+-"+var);
				amSeriesVar.get(am).add(i, avg, avg-var, avg+var );
			}
		}
	
		//Normalize???
//		for (int i = 0; i< b.size(); i++) {
//			seriesB.add(i,((b.get(ranks.get(i)).doubleValue() -minb)  / (maxb-minb) )*(maxa-mina) + mina);
//			Logger.log( ((b.get(ranks.get(i)).doubleValue() -minb) / (maxb-minb))*(maxa-mina) + mina +",", DebugMode.detailed);
//		}
	//	XYSeriesCollection collection = new XYSeriesCollection();
		YIntervalSeriesCollection collection = new YIntervalSeriesCollection();
		
		 for (int i=0; i < amSeriesVar.size(); i++) {
	//		 collection.addSeries(amSeries.get(i));
			 collection.addSeries(amSeriesVar.get(i));
		 }
	
		JFreeChart chart = ChartFactory.createXYLineChart("", "Random Partitionings", "Index", collection, PlotOrientation.VERTICAL, true, true, false);
//		JFreeChart chart = ChartFactory.createXYLineChart(dataset.name.substring(0,dataset.name.indexOf('.'))+" "+am.toString()+" "+crit.toString(), "Sample Partitionings", "Value", collection, PlotOrientation.VERTICAL, true, true, false);
	//	chart.setBackgroundPaint(new Color(0, 0, 0, 0));
		chart.getPlot().setBackgroundPaint(Color.white);
		chart.getXYPlot().getRangeAxis().setLabelFont(new Font("sansserif", Font.BOLD, 20));
		chart.getXYPlot().getDomainAxis().setLabelFont(new Font("sansserif", Font.BOLD, 20));
		chart.getXYPlot().getRangeAxis().setTickLabelFont(new Font("serif", Font.BOLD, 8));
		chart.getXYPlot().getDomainAxis().setTickLabelFont(new Font("serif", Font.BOLD, 10));
		 
//		Font font = new Font("SansSerif", 0, 9);
//		 
//		 XYTextAnnotation xytextannotation = null;
//	        xytextannotation = new XYTextAnnotation("3rd", 2D, 2D);
//	        xytextannotation.setFont(font);
//	        xytextannotation.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
//	        chart.getXYPlot().addAnnotation(xytextannotation);
	        
	//	chart.getLegend().setPosition(RectangleEdge.RIGHT); 
		if(mode!=Mode.RAND2RAND)
			chart.getXYPlot().addDomainMarker(new ValueMarker(dataset.groundTruth.size(),Color.darkGray,new BasicStroke(2.0f)));
		
		chart.getXYPlot().setDomainGridlinesVisible(true);
		chart.getXYPlot().setDomainGridlinePaint(Color.lightGray);
//		chart.getXYPlot().setRangeGridlinePaint(Color.lightGray);
        
		chart.getXYPlot().getDomainAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		chart.getLegend().setItemFont(new Font("sansserif", Font.BOLD, 24));
		chart.setBackgroundPaint(Color.white);
		//chart.getXYPlot().getRangeAxis().setUpperBound(1);//Range(new Range(0, 1));
		chart.getXYPlot().getDomainAxis().setLowerBound(1);
		//http://www.rapidtables.com/web/color/RGB_Color.htm
//		Color[] colors = {new Color(0xCC0000), new Color(0x004C99), new Color(0x009900), new Color(0xFFFF00) ,new Color(0xFF8000), new Color(0x60606060),new Color(0xCC00CC),new Color(0xCC0066),  Color.RED, Color.YELLOW, Color.BLUE};
//		Color[] Tcolors = { new Color(0xFFCCCC), new Color(0xCCE5FF), new Color(0xCCFFCC), new Color(0xFFFFCC) , new Color(0xFFE5CC), new Color(0xE0E0E0E0),new Color(0xFFCCCC),new Color(0xFFCCE5), Color.RED, Color.YELLOW, Color.BLUE};
		Color[] colors = {new Color(0xb70408), new Color(0x86a24b), new Color(0x2880b7) ,new Color(0xe8e013), new Color(0x70578e),new Color(0xf18e3d),new Color(0x363e4d),  Color.RED, Color.YELLOW, Color.BLUE};
		Color[] Tcolors = { new Color(0xeab6b8), new Color(0xd8e8b7), new Color(0xd3e7f4) , new Color(0xe8e482), new Color(0xb9a5d1),new Color(0xf1cbac),new Color(0xc9d5ec), Color.RED, Color.YELLOW, Color.BLUE};

		DeviationRenderer renderer = new DeviationRenderer(true, false);//set shapes on/off
		for (int i=0; i <amSeriesVar.size(); i ++) {
			renderer.setSeriesStroke(   i,  new BasicStroke(4.0f,1,1, 1.0f, new float[] {i*4+0.f, i+2.0f}, 0.0f));
			renderer.setSeriesPaint(i, colors[i]);
			renderer.setSeriesFillPaint(i,  Tcolors[i]);
			renderer.setSeriesItemLabelsVisible(i, true);
			//SeriesShape(i, ShapeUtilities.createDiagonalCross(3, 1));
		 }
		
		//renderer.setAutoPopulateSeriesPaint(true);//ulateSeriesShape(true);
		
		{
        chart.getXYPlot().setRenderer(renderer);

	        // change the auto tick unit selection to integer units only...
	        //NumberAxis yAxis = (NumberAxis)  chart.getXYPlot().getRangeAxis();
	        //yAxis.setAutoRangeIncludesZero(false);
	    //    yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		}
		
		//new BasicStroke(       4.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,      1.0f, new float[] {12.0f, 20.0f}, 0.0f    ));
	//	chart.addSubtitle(new TextTitle( "s="+splitChance+"%, m="+ mergeChance +"%, r="+swapChance+ "%, " + logDatasetStat(partitionings)  ));
	//	ChartUtilities.writeChartAsJPEG(new FileOutputStream("cct2"), 1, chart, 500, 300, info);
		try {
			new File(expPath+"/img/").mkdir();
			File path = new File(expPath+"/img/ext");
			path.mkdir();
			ChartUtilities.saveChartAsPNG(new File(path.getPath() +"/"+""+dataset.name+"_"+mode+"_"+BUCK+"__"+swapChance+"%_"+(formatter.format(System.currentTimeMillis()))+".png"), chart, 800, 600);
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}
	
	
	
	int MINK = 2, MAXK=24, NR = 500;//maxk=24 for strike 
	int BUCK = 20
//			100
			, MAX = 35000, MaxRun=1000;//1000;//100, 5000
	{
		splitChance = 30;//20;
		mergeChance = 5;//30; //was 20 for the rand2G, 
		swapChance = 35;//25;
		maxResult = 500;//1000;//180;
	}
	public enum Mode{KNEE,FRAGKNEE,RAND2G,RAND2RAND};
	public Mode mode;
	{
		mode = Mode.RAND2G;
		mode = Mode.FRAGKNEE;
	}
	public static void main(String[] params) {
		String args[] = params ;
	//	args = new String[] {"exps/GN", "+G"};
	//	args = new String[] {"exps/wkarate", "+G"};

//		args = new String[] {"exps/strike", "+G"};
		args = new String[] {"exps/real", "+G"};

		//	args = new String[] {"exps/foot", "+G"};
	//	args = new String[] {"exps/karate", "+G"};

		
		while(args.length<=0 || args[0].contains("help")){
			System.out.println("---- You have following options: ");
			System.out.println("+G: consider randomized versions from the groundtruth");
			System.out.println("+plot: plot all correlation images");
			System.out.println("N=: the number of random clusterings that should be generated in each run");
			System.out.println("s=: split chance for the clustering randomizer");
			System.out.println("m=: merge chance for the clustering randomizer");
			System.out.println("r=: swap chance for the clustering randomizer");
			System.out.println("M=: the total number of random clusterings that should be generated");
			System.out.println("R=: the maximum number of runs for random clusterings that should be generated");
			System.out.println("B=: the total number of random clusterings with specific number of clusters that should be generated");
			System.out.println("Please provide path of the experiments followed by parameters; or type help.");
			Scanner in = new Scanner(System.in);
			args = in.nextLine().split(" ");
//			System.err.println(args.length);
			for (int i = 0; i < args.length; i++) {
//				System.err.println(args[i]);
			}
		    in.close();       
		}
		
		final ExternalIndexComparer<Integer, Integer> comparer = new ExternalIndexComparer<Integer, Integer>(args[0]);
		Logger.setLog_level(DebugMode.normal);
		new Logger(comparer.expPath); 
		
		Logger.logFunction(comparer.expPath);
				

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
			if(string.contains("M="))
				comparer.MAX = Integer.parseInt(string.substring(2));
			if(string.contains("B="))
				comparer.BUCK = Integer.parseInt(string.substring(2));
			if(string.contains("R="))
				comparer.MaxRun = Integer.parseInt(string.substring(2));
			if(string.contains("L="))
				Logger.setLog_level(DebugMode.values()[  Integer.parseInt(string.substring(2)) ] );
			if(string.contains("help")){
				System.out.println("You have following options: ");
				System.out.println("+G: consider randomized versions from the groundtruth");
				System.out.println("+plot: plot all correlation images");
				System.out.println("N=: the number of random clusterings that should be generated in each run");
				System.out.println("s=: split chance for the clustering randomizer");
				System.out.println("m=: merge chance for the clustering randomizer");
				System.out.println("r=: swap chance for the clustering randomizer");
				System.out.println("R=: the maximum number of runs for random clusterings that should be generated");
				System.out.println("M=: the total number of random clusterings that should be generated");
				System.out.println("B=: the total number of random clusterings with specific number of clusters that should be generated");
				System.out.println("Please provide path of the experiments followed by parameters; or type help.");
				return;
			}

		}
		Logger.logln("N="+ comparer.maxResult +", s="+ comparer.splitChance  +", m="+ comparer.mergeChance +", r="+ comparer.swapChance );
		comparer.compareExternalIndexes();
	}

}
