package algorithms.communityMining.external_methods;


import io.graph.GraphOutputStream;
import io.graph.pairs.PairsGraphWriter;
import io.group.GroupingReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Vector;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.data.Partitioning;
import edu.uci.ics.jung.algorithms.filters.VertexPredicateFilter;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/*
 
LP: http://perso.uclouvain.be/vincent.blondel/research/louvain.html
Bondel: https://sites.google.com/site/findcommunities/
RN: http://physics.wustl.edu/zohar/communitydetection/
WT: http://www-rp.lip6.fr/~latapy/PP/walktrap.html
Infomap,mod: http://www.tp.umu.se/~rosvall/code.html
DM: http://wdb.ugr.es/~donetti/

Hunter: http://hub.iis.sinica.edu.tw/Hunter/
Netcarto (also role analysis!): http://etseq.urv.cat/seeslab/downloads/
MCL: http://www.micans.org/mcl/

*/
public abstract class CommunityMinerExecutableWrapper<V,E> extends CommunityMiner<V, E>{
	protected GraphOutputStream<V, E> graphWriter ;
	protected GroupingReader< V> groupingReader;
	protected String weightFormat = "%s", extention, extentionWeighted ;
	protected int startIndexId = 0;
	protected int n;
	
	protected boolean verbousMode = false;
	/**
	 * 
	 * Runs an executable comment for running the approaches implemented by original authors in other languages
	 * See http://alvinalexander.com/java/edu/pj/pj010016
	 * 
	 * @param command
	 * @throws IOException
	 * @throws InterruptedException
	 */
	protected synchronized void runExecCommand(String command) throws IOException, InterruptedException{
		runExecCommand(command, null);
	}
	
	// runs the command and redirects output to given path if provided 
	protected synchronized void runExecCommand(String command, String outputPath) throws IOException, InterruptedException{
		System.out.println(command);
		Process child = Runtime.getRuntime().exec(command);
		Scanner sc = new Scanner(child.getInputStream());    		
		FileOutputStream outputStream = null;
		String output ;
		if(outputPath!=null)
			outputStream = new FileOutputStream(outputPath);
		while ( sc.hasNext()) {
			output =sc.nextLine();
			if(verbousMode) System.out.println(output);
			if(outputPath!=null)
				outputStream.write((output+"\n").getBytes());
		}
		if(outputPath!=null){
			outputStream.flush();
			outputStream.close();
		}
		if(!child.waitFor(1, TimeUnit.HOURS)) {
			child.destroy(); 
		}
		sc.close();
		sc = new Scanner(child.getErrorStream());    		
		while (sc.hasNext()){
			output =sc.nextLine();
			if(verbousMode) System.out.println(output);
		}
		sc.close();

		if(verbousMode){
			int exitVal = child.waitFor();
			System.out.println("Exited with error code " + exitVal);
		}
		
//		TODO: might be more efficient if calling C++ code directly, or using pipe instead of writing the actual file
//		http://www.javaworld.com/article/2077513/learn-java/java-tip-17--integrating-java-with-c--.html
		
	}
	
	
	
	public  abstract FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)  throws IOException, InterruptedException;
	public  FileInputStream runCommunityMiningExecutable(Graph<V,E> graph, Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)  throws IOException, InterruptedException{
		return runCommunityMiningExecutable(vertexTransformer, network, path, isWeighted);
	}
	//public abstract Grouping<V> readPartitioning(InputStream is,	 Transformer<String,V> vertexTransformer) throws IOException ;
	public String tmpFilename="tmpNet";
	public synchronized Grouping<V> findCommunities(Graph<V,E> graph, Map<E, ? extends Number> weights){
		n = graph.getVertexCount();
		String network = getName()+tmpFilename, path = "./temp";
		new File(path).mkdir();
		
		Map<V, Integer> vertexLabels = new HashMap<V, Integer>();
		final Vector<V> vertexes = new Vector<V>(graph.getVertices());
		for(int i = 0; i< vertexes.size(); i++){
			vertexLabels.put(vertexes.elementAt(i), i+startIndexId);
		}
		
		Grouping<V> detectedCommunities = null ;
		try {
			graphWriter.writeGraph(path +"/"+ network + (weights==null?extention:extentionWeighted), graph, vertexLabels, weights, weightFormat); 
		
			Transformer<String, V> vertexTransformer = new Transformer<String, V>() {
				public V transform(String arg0) {
					return vertexes.get(Integer.parseInt(arg0)- startIndexId);
				}
			};
			
			detectedCommunities =  groupingReader.readPartitioning(runCommunityMiningExecutable(graph, vertexTransformer,  network,  path,(weights==null?false:true)), vertexTransformer); 
		
			detectedCommunities.removeEmptyGroups();
			if (!verbousMode)
			runExecCommand("rm -rf " + path);
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return detectedCommunities;

	}

}
