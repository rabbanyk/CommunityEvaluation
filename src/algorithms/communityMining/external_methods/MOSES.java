package algorithms.communityMining.external_methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;

import org.apache.commons.collections15.Transformer;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 * Wrapper for authors Implementations 
 * 
 * https://sites.google.com/site/aaronmcdaid/moses
 * 
 * Unweighted!
 * 
 * McDaid, Aaron, and Neil Hurley. 
 * "Detecting highly overlapping communities with model-based overlapping seed expansion." 
 * Advances in Social Networks Analysis and Mining (ASONAM), 2010 International Conference on. IEEE, 2010.
 * 
 * In Makefile, the booslibrary is updated to the installed library in /usr/share/
 * It fails if there is an empty line at the end of the input file
 * Also minor change to remove an unused parameter ID, otherwise gives a warning which will be treated as error with the make parameters, -Werror
 * 
 * @param <V>
 * @param <E>
 */
public class MOSES <V, E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath ="./execs/CM-Overlapping-MOSES-McDaid/moses-2011-01-26/";

	{
		graphWriter = new PairsGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".pairs";
		extentionWeighted = ".pairs";
		startIndexId = 0;
		groupingReader = new ListGrouingReader<>();
	}
	
	@Override
	public FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)
			throws IOException, InterruptedException {

		runExecCommand(exePath+"moses " +  path + "/" + network + extention +
				" "+ path + "/" + network + ".res" ) ; //weights?

		return new FileInputStream(path + "/" +  network + ".res");
	}

	@Override
	public String getName() {
		return "MOSES";
	}

	@Override
	public String getShortName() {
		return "MOSES";
	}

}
