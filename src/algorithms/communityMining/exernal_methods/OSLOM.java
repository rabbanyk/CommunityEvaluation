package algorithms.communityMining.exernal_methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.CommunityMinerExecutableWrapper;
//http://www.oslom.org/
//https://sites.google.com/site/andrealancichinetti/software

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 * Wrapper for authors Implementations 
 *	
 *	gives hierarchy not used currently
 * @param <V>
 * @param <E>
 */
public class OSLOM <V, E> extends CommunityMinerExecutableWrapper<V, E> {

	
	String exePath ="./execs/CM-Overlapping-OSLOM-Lanci/OSLOM2/";

	{
		graphWriter = new PairsGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".dat";
		extentionWeighted = ".dat";
		startIndexId = 1;
		groupingReader = new ListGrouingReader<>("#");
	}
	
	@Override
	public FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)
			throws IOException, InterruptedException {
		
		// in case of directed "oslom_dir" 
		runExecCommand(exePath + "oslom_undir -f "+ path + "/" + network + extention 
				+(isWeighted ? " -w ": " -uw ") + " -hr 0  " ) ; 
		//-singlet to detect outliers, as single modules, otherwise they'll be assigned to a module
		//-hr 0 to detect no higher hierarchy levels
		return new FileInputStream(path + "/" + network + extention+ "_oslo_files/tp");
	}
	
	//TODO: weighetd version give strange results, check it! increasing hr won't help, all nodes in one comm for wkarate
	@Override
	public String getName() {
		return "OSLOM";
	}

	@Override
	public String getShortName() {
		return "OSLOM";
	}

}
