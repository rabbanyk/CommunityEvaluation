package algorithms.communityMining.methods;

import io.graph.GraphOutputStream;
import io.graph.pairs.PairsGraphWriter;
import io.group.GroupingReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.CommunityMinerExecutableWrapper;
import algorithms.communityMining.data.Grouping;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)

 * http://www-rp.lip6.fr/~latapy/PP/walktrap.html

 *	Wrapper for authors original implementation. 
 *	minor changes by adding include cstdlib, cstring and algorithm in graph.cpp to be able to compile
	
 * @param <V>
 * @param <E>
 */
public class WalkTrap<V,E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath = "./execs/CM-Weighted-Walktrap-Pons/walktrap"; 

	{
		graphWriter = new PairsGraphWriter< V, E>();
		// Dummy to test if not considering weights, walktrap performance would dropcp 
//		graphWriter = new GraphOutputStream<V, E>(){
//			@Override
//			protected String formatEdge(int v1, int v2, String weight) {
//				return ((v1) + "\t" + (v2) + "\n");
//			}};
		weightFormat = "%s"; 
		extention = ".net";
		extentionWeighted = ".net";
		startIndexId = 0; // integers starting from 0
	}
	

	
	@Override
	public FileInputStream runCommunityMiningExecutable(
			Transformer<String, V> vertexTransformer, String network,
			String path, boolean isWeighted) throws IOException,
			InterruptedException {
		//walktrap [input_file] [-o output_file] [-i index_file] [options] // see read me for more options 
		//TODO: can configure options: # of clusters, next best partition, length of random walk, etc...
		//TODO: can improve by directly reading input and output stream, no need for actual hard access
	    //"walktrap network.dat -o communities.txt -t5 -b"
	
		String args = " -d1 -b ",
				inputPath = path + "/" + network + extention,
				outputPath =  path + "/" + network +".com ";
		runExecCommand(exePath+ " "+ inputPath + " -o " + outputPath + args );

		return  new FileInputStream(path + "/" + network +".com");	
		}
	{
		groupingReader = new GroupingReader<V>() {
			@Override
			public Grouping<V> readPartitioning(InputStream is,
					Transformer<String, V> vertexTransformer) throws IOException {
				
				/*	Maximal modularity Q = 0.440181 for partition :
					community 55 = {27, 31, 2, 26, 29}
					community 63 = {28, 33, 11, 32, 30, 12, 22, 17, 20, 24, 19, 0, 4}
					community 61 = {7, 13, 18, 8, 9}
					community 62 = {14, 25, 1, 21, 15, 5, 16, 23, 6, 3, 10}
				*/

				Grouping<V> partitioning =  new Grouping<V>();
				BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));

				String tmp;
				Set<V> cc = new HashSet<V>();
				while ((tmp = bufferedInputStream.readLine()) != null) if(tmp.length()>0){
					if (tmp.startsWith("community", 0)) {
						cc = new HashSet<V>();
						partitioning.addGroup(cc);
						tmp = tmp.substring(tmp.indexOf('{')+1, tmp.indexOf('}'));
					
						for( String v : tmp.split("[,\\s]+"))
							cc.add(vertexTransformer.transform(v));
					}
				}
				bufferedInputStream.close();
				
				return partitioning;
			}
		};
	}
	

	@Override
	public String getName() {
		return "WalkTrap";
	}

	@Override
	public String getShortName() {
		return "WT";
	}

}
