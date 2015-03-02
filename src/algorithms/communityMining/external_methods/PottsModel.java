package algorithms.communityMining.external_methods;

import io.graph.gml.GMLGraphWriter;
import io.group.GroupingReader;
import io.group.ListGrouingReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import algorithms.communityMining.data.Grouping;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)

 *  http://physics.wustl.edu/zohar/communitydetection/
 *  "Multiresolution community detection for megascale networks by information-based replica correlations," Physical Review E vol. 80 article 016109, by P. Ronhovde and Z. Nussinov
 *   ''solving a set of independent solutions ("replicas") for a range of resolutions using a local Potts model''
 *   
 *   Minor change  in MSL_List_Template.h to be able to compile with g++, line 615: 
 *   pn = new TListNode<T>(v);
 * 	//pn->n       = v;          // copy node data
 *
 * @param <V>
 * @param <E>
 */
public class PottsModel<V,E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath = "./execs/CM-Weighted-Potts-Ronhovde"; 

	{
		graphWriter = new GMLGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".gml";
		extentionWeighted = ".gml";
		startIndexId = 0; // integers starting from 0
//		parameters = new String[]{"-nsm:300", "-nim:100",  "-nrm:10", "-ntm:4", "-gstart:19", "-gstop:0.001", "-gsteps:20", "-v:2", "-rseed:8231593"};
		//	verbousMode = true;
		groupingReader = new ListGrouingReader<V>();
	}
	
	@Override
	public FileInputStream runCommunityMiningExecutable(Graph<V, E> graph,
			Transformer<String, V> vertexTransformer, String network,
			String path, boolean isWeighted) throws IOException,
			InterruptedException {

		//./rnMRAw -n:34 -nsm:300 -nim:100 -nrm:10 -ntm:4 -gstart:19 -gstop:0.001 -gsteps:20 -v:2 -rseed:8231593 
		//-inf:karate_weighted_nooffset.gml -infans:karate_answer_nooffset.txt -outf:MRAKarateTestW
		//./rnCDw -n:34 -nsm:300 -nim:1000 -ntm:10 -g:0.1 -v:2 -rseed:2436083 -inf:karate_weighted_nooffset.gml -infans:karate_answer_nooffset.txt -outf:APMCDTestW
		
		String args = " -nsm:300 -nim:100 -nrm:10 -ntm:4 -gstart:19 -gstop:0.001 -gsteps:20 -v:2 -rseed:8231593 ";
		runExecCommand(exePath+ "/rnMRA"+(isWeighted?"w ":" ")+ 
				"-n:"+ graph.getVertexCount()+ args +"-inf:" + path + "/" + network + extention + " -outf:" + path + "/" + network  );

		return  new FileInputStream(path + "/" + network +"Best.txt");	//MRAKarateTestWBest.txt
	}

	@Override
	public FileInputStream runCommunityMiningExecutable(
			Transformer<String, V> vertexTransformer, String network,
			String path, boolean isWeighted) throws IOException,
			InterruptedException {
		throw new IOException("Call the runCommunityMiningExecutable instead, it needs graph meta data!");
	}
	
//	{
//		groupingReader = new GroupingReader<V>() {
//	public Grouping<V> readPartitioning(InputStream is,	 Transformer<String,V> vertexTransformer) throws IOException {
//		Grouping<V> partitioning =  new Grouping<V>();
//		BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));
//
//		String tmp;
//		int clusterId = 0;
//		
//		while ((tmp = bufferedInputStream.readLine()) != null) if(tmp.length()>0){
//			partitioning.addGroup(new HashSet<V>());
//			
//			for (String vId :tmp.split("[,\\s]+"))
//				partitioning.getGroup(clusterId).add(vertexTransformer.transform(vId));
//			clusterId++;
//		}
//		
//		bufferedInputStream.close();
//		
//		return partitioning;
//	}};}

	@Override
	public String getName() {
		return "PottsModel";
	}

	@Override
	public String getShortName() {
		return "RN";
	}


}
