package algorithms.communityMining.methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.PairGrouingReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.collections15.Transformer;

import util.DatasetUtils;
import util.DatasetUtils.ClassicDataset;
import data.GraphDataSet;
import algorithms.AlgorithmUtils;
import algorithms.AlgorithmUtils.Method;
import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.CommunityMinerExecutableWrapper;
import algorithms.communityMining.data.Grouping;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 * 
 *         https://sites.google.com/site/findcommunities/ in NetworkX http://perso.uclouvain.be/vincent.blondel/research/louvain.html
 * 			weights: OK
 * 			overlapping: No
 * 			hierarchy: Yes
 */
public class Louvain<V, E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath = "./execs/CM-Weighted-Louvain-Blondel/Community_latest/";
	{
		graphWriter = new PairsGraphWriter<V, E>();
		weightFormat = "%s";
		extention = ".pairs";
		extentionWeighted = ".pairs";
		startIndexId = 0; // integers starting from 0
		groupingReader = new PairGrouingReader<V>();
	}

	int level = -1;
	//TODO: we can find the best using agreement with a particular attribute instead of maximizing the modularity

	@Override
	public synchronized FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)
			throws IOException, InterruptedException {
		// ./convert -i input_file -o outfile [-r] [-w outfile_weight]
		// ./community input_file [-w weight_file] [-p part_file] [-q epsilon] [-l display_level] [-v] [-h]
		// ./hierarchy input_file [options]

		runExecCommand(exePath + "convert -i " + path + "/" + network + extention + " -o " + path + "/" + network + ".bin "	+ (isWeighted ? " -w "+ path + "/" +network+".weights" : ""));
		runExecCommand(exePath + "community " + path + "/" + network + ".bin " + (isWeighted ? " -w "+ path + "/" +network+".weights" : "") + " -l -1 ", path + "/" + network + ".tree");
		
		//To print out how many levels are in tree
		runExecCommand(exePath + "hierarchy " + path + "/" + network + ".tree" , path + "/" + network + ".cominfo");
		
		//Read the number of levels
		BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(path + "/" + network +".cominfo")));			  
		 String line = null;
		 String InfoMarker = "Number of levels: ";
		 while ((line = input.readLine()) != null) {
			 if (line.startsWith(InfoMarker))
				 level = Integer.parseInt(line.substring(InfoMarker.length()).trim()); 
		 }
		 input.close();
		runExecCommand(exePath + "hierarchy " + path + "/" + network + ".tree" + " -l "+(level-1)+" ", path + "/" + network + ".com");
		return new FileInputStream(path + "/" + network + ".com");
	}

	

	/*
	 * 
	 * From Readme File:
	 * 
	 * 1. Conversion from a text format (each line contains a couple "src dest") ./convert -i graph.txt -o graph.bin This program can also be used to convert
	 * weighted graphs (each line contain a triple "src dest w") using -w option: ./convert -i graph.txt -o graph.bin -w graph.weights Finally, nodes can be
	 * renumbered from 0 to nb_nodes - 1 using -r option (less space wasted in some cases): ./convert -i graph.txt -o graph.bin -r
	 * 
	 * 
	 * 2. Computes communities and displays hierarchical tree: ./community graph.bin -l -1 -v > graph.tree
	 * 
	 * To ensure a faster computation (with a loss of quality), one can use the -q option to specify that the program must stop if the increase of modularity is
	 * below epsilon for a given iteration or pass: ./community graph.bin -l -1 -q 0.0001 > graph.tree
	 * 
	 * The program can deal with weighted networks using -w option: ./community graph.bin -l -1 -w graph.weights > graph.tree In this specific case, the
	 * convertion step must also use the -w option.
	 * 
	 * The program can also start with any given partition using -p option ./community graph.bin -p graph.part -v
	 * 
	 * 
	 * 3. Displays information on the tree structure (number of hierarchical levels and nodes per level): ./hierarchy graph.tree
	 * 
	 * Displays the belonging of nodes to communities for a given level of the tree: ./hierarchy graph.tree -l 2 > graph_node2comm_level2
	 */
	
	@Override
	public String getName() {
		return "Louvain";
	}

	@Override
	public String getShortName() {
		return "LV";
	}

}
