package algorithms.communityMining.methods;

import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.CommunityMinerExecutableWrapper;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 * Wrapper for authors Implementations 
 * 
 * Supports bipartite networks, not used currently
 * @param <V>
 * @param <E>
 */
public class COPRA<V, E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath ="./execs/CM-Overlapping-COPRA-Gregory/";

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
		//-prop p : Limits the maximum number of iterations to p
		//-v the maximum number of communities per vertex; v=1 is default.
		int maxNumberOFcommunitiesperV  =2;
		int repeat= 10;
		runExecCommand("java -cp "+exePath+"copra.jar COPRA " +  path + "/" + network + extention 
				+(isWeighted ? " -w ": "") + " -v "+maxNumberOFcommunitiesperV 
				+(repeat==0?"": " -mo -repeat "+repeat)) ; //
		// writes results in the current directory!

		runExecCommand("mv clusters-" + network + extention +" "+  path +"/" ) ;
		if (repeat!=0)		runExecCommand("mv best-clusters-" + network + extention +" "+  path +"/" ) ;

		return new FileInputStream(path + "/" +(repeat==0? "clusters-":"best-clusters-")+  network + extention);//best-clusters-
		// -bi bipartite
		
//		To compute the overlap modularity of a cover with respect to a network, use:
//			java -cp copra.jar ModOverlap karate.txt clusters-karate.txt
	}

	@Override
	public String getName() {
		return "COPRA";
	}

	@Override
	public String getShortName() {
		return "COPRA";
	}

}
