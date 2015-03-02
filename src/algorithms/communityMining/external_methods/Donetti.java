package algorithms.communityMining.external_methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.PairGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;

import org.apache.commons.collections15.Transformer;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 * 
 * DM: http://wdb.ugr.es/~donetti/
 *	Wrapper for authors original implementation. 
 *  TODO: Performs worse than expected, should be a bug, or the number of eignvalue should not be constant of 3!
 *  
 * @param <V>
 * @param <E>
 */
public class Donetti<V,E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath = "./execs/CM-Commfind-Donetti/commfind"; 
	{
		graphWriter = new PairsGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".net";
		extentionWeighted = ".net";
		startIndexId = 0; // integers starting from 0
		groupingReader = new PairGrouingReader<V>();
	}
	
	int numberOfEignValues = 3;
	//TODO: we can find the best say the one that maximizes the modularity

	@Override
	public synchronized FileInputStream runCommunityMiningExecutable(
			Transformer<String, V> vertexTransformer, String network,
			String path, boolean isWeighted) throws IOException,
			InterruptedException {
//			execs/commfind-0.2$ ./commfind example.net
		String args =" "+ numberOfEignValues+ " ", //number of eign values
				inputPath = path + "/" + network + extention,
				outputPath =  path + "/" + network +".com ";
		runExecCommand(exePath+ " "+ inputPath  + args , outputPath);

		return  new FileInputStream(outputPath);	
		}
	

	@Override
	public String getName() {
		return "Donetti";
	}

	@Override
	public String getShortName() {
		return "DN";
	}

}
