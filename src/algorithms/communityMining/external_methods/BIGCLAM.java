package algorithms.communityMining.external_methods;

import io.graph.GraphOutputStream;
import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;

public class BIGCLAM <V, E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath ="./execs/CM-Overlapping-AGM-SNAP/bigclam/";

	{
		graphWriter = new PairsGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".edgelist";
		extentionWeighted = ".edgelist";
		startIndexId = 0;
		groupingReader = new ListGrouingReader<>();
//		verbousMode = true;
	}
	
	@Override
	public FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)
			throws IOException, InterruptedException {
		
		int c = -1, mc = 2, xc = (int) Math.log(n), nc = 10;
		runExecCommand(exePath+"bigclam "
				+ " -i:" +  path + "/" + network + extention
//				+ " -l:" +  path + "/" + network + extention+".labels"  
				+ " -o:" +  path + "/" + network +".res" 
				+ " -c:" + c 
				+ " -mc:" + mc
				+ " -xc:" + xc
				+ " -nc:" + nc
				);
//		agmfitmain -i:football.edgelist -l:football.labels -c:12 -e:0.1
		return new FileInputStream(path + "/" +  network + ".res");
	}
	
	@Override
	public String getName() {
		return "BIGCLAM";
	}
	
	@Override
	public String getShortName() {
		return "BIGCLAM";
	}
	
}
