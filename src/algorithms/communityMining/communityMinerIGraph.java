package algorithms.communityMining;

import io.graph.gml.GMLGraphWriter;
import io.graph.pajek.PajekGraphWriter;
import io.group.ListGrouingReader;
import io.group.PairGrouingReader;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.collections15.Transformer;
import org.python.modules.synchronize;

import algorithms.communityMining.data.Grouping;
import edu.uci.ics.jung.graph.Graph;
//
//import org.python.*;
//import org.python.core.PyObject;
//import org.python.util.PythonInterpreter;

/*
 *
 *
iGraph in R: http://igraph.sourceforge.net/download.html
http://igraph.sourceforge.net/documentation.html

Graph.community_edge_betweenness           
Graph.community_leading_eigenvector        
Graph.community_spinglass
Graph.community_fastgreedy                 
Graph.community_leading_eigenvector_naive  
Graph.community_walktrap
Graph.community_infomap                    
Graph.community_multilevel                 
Graph.community_label_propagation          
Graph.community_optimal_modularity         

https://wiki.python.org/jython/NewUsersGuide
 *
 */
public  class communityMinerIGraph<V,E> extends CommunityMinerExecutableWrapper<V, E>{
	public enum Method { edge_betweenness, leading_eigenvector, spinglass, fastgreedy, leading_eigenvector_naive,
		walktrap, infomap, multilevel, label_propagation, optimal_modularity    }
	Method method  = Method.label_propagation;
	{
//		verbousMode = true;
	}
	public communityMinerIGraph() {
		super();
	}
	public communityMinerIGraph(Method method) {
		super();
		this.method = method;
	}

	{
		graphWriter = new GMLGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".gml";
		extentionWeighted = ".gml";
		startIndexId = 0; // integers starting from 0
		verbousMode = true;
		groupingReader = new ListGrouingReader< V>();
	}
	
	public synchronized FileInputStream runCommunityMiningExecutable(
			Transformer<String, V> vertexTransformer, String network,
			String path, boolean isWeighted) throws IOException,
			InterruptedException {
		String exePath = "./src/algorithms/communityMining/methods/iGraph/";
		String args ="";
		runExecCommand("python " + exePath+ "communityMinerInterface.py  "
				+" -i "+ path+ "/" + network + extention+
				" -o "+ path + "/" + network +".com  "+
				" -m "+ method+args);
		FileInputStream res =null;
		try{
		 res = new FileInputStream(path + "/" + network +".com");	
		}catch (Exception e){
			e.printStackTrace();
		}
		return  res;
		}
	
	
	@Override
	public String getName() {
		String name =  method.toString();
		return name.substring(0,1).toUpperCase()+name.substring(1);
//		return "iGraph_" +method;

	}

	@Override
	public String getShortName() {
		return "i"+method;
	}

}
