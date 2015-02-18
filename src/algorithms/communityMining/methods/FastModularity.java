package algorithms.communityMining.methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.GroupingReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.CommunityMinerExecutableWrapper;
import algorithms.communityMining.data.Grouping;


/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 *
 *	Wrapper for FastModularity Community mining approach provided here:
 *	http://www.cs.unm.edu/~aaron/research/fastmodularity.htm
 *	Makefile and c++ sources slightly changed to be compatible with newer G++ complier
 *
 *	weights: OK
 * 	overlapping: No
 * 
 * @param <V> Type of Vertices
 * @param <E> Type of Edges
 */
public class FastModularity<V,E> extends CommunityMinerExecutableWrapper<V,E> {
	String exePath = "./execs/CM-Weighted-FastModularity-Clauset/FastCommunity_GPL_v1.0.3/FastCommunityMH"; 
	String exePathWeighted =  "./execs/CM-Weighted-FastModularity-Clauset/FastCommunity_w_GPL_v1.0.1/FastCommunity_wMH";
	{
		graphWriter = new PairsGraphWriter< V, E>();
		weightFormat = "%.0f"; //only works with integer weights
		extention = ".pairs";
		extentionWeighted = ".wpairs";
		startIndexId = 0;
	}
	public FastModularity(){
		
	}
	
	/*
	 * ./FastCommunityMH -f community500.pairs -l firstRun
	 * ./FastCommunityMH -f community500.pairs -l secondRun -c 494
	 * 
	 *  .wpairs for weighted, weights should be integer
	 * Based on original authors description: Input file requirements
	 *  SINGLE COMPONENT, NO SELF-LOOPS or MULTI-EDGES 
	 *  MINIMUM NODE ID = 0 , MAXIMUM NODE ID can be anything, but the program will use less memory if nodes are labeled sequentially
	 * */
	
	public synchronized FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted) throws IOException, InterruptedException {
			// This algorithm requires pairs format and the network need .pairs extension
			//./FastCommunityMH -f community500.pairs -l firstRun
			runExecCommand((isWeighted?exePathWeighted:exePath)+ " -f " + path + "/" + network + (isWeighted?".wpairs":".pairs")+"  -l firstRun");

			BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(path + "/" + network +"-fc_firstRun.info")));			  
			 String line = null;
			 int step = 0;
			 String STEP = "STEP------:";
			 while ((line = input.readLine()) != null) {
				 if (line.startsWith(STEP))
					 step = Integer.parseInt(line.substring(STEP.length()).trim()); 
			 }
			 input.close();

			 //./FastCommunityMH -f community500.pairs -l secondRun -c 494
			runExecCommand((isWeighted?exePathWeighted:exePath)+" -f " + path + "/" + network +   (isWeighted?".wpairs":".pairs")+"  -l secondRun -c "+ step);
		
			return  new FileInputStream(path + "/" + network + "-fc_secondRun.groups");

	}
	
	{
		groupingReader = new GroupingReader<V>() {
	/*GROUP[ 3 ][ 49 ]
	3
	10
	15
	16
	GROUP[ 3 ][ 49 ]
	55
	19
	12
	*/
	public Grouping<V> readPartitioning(InputStream is,	 Transformer<String,V> vertexTransformer) throws IOException {
		Grouping<V> partitioning =  new Grouping<V>();
		BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));

		String tmp;
		while ((tmp = bufferedInputStream.readLine()) != null) if(tmp.length()>0){
			if (tmp.startsWith("GROUP[", 0)) {
				partitioning.addGroup(new HashSet<V>());
			} else 
				partitioning.getLastGroup().add(vertexTransformer.transform(tmp));
		}
		bufferedInputStream.close();
		return partitioning;
	}
		};
	
		}
	
	public String getName() {
		return "FastModularity";
	}

	public String getShortName() {
		return "FM";
	}
	
	
}
