package io.graph.pajek;

import io.graph.GraphOutputStream;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;

public class PajekGraphWriter<V,E> extends GraphOutputStream<V, E>{
	@Override
	protected String verticesMetaData(Vector<V> vertexes){
		//*Vertices     26
		return ("*Vertices     " + (vertexes.size()) + "\n");
	}
	@Override
	protected String formatVetice(V v1,Transformer<V, Integer> vertex_Ids,Map<V, HashMap<Object,Object>> vertex_attributes){
		//     1 "pe0"   //TODO: write the attributes
		return (vertex_Ids.transform(v1) +"\t\"" + (vertex_Ids.transform(v1)) + "\"\n");
	}
	//*Edges
	@Override
	protected String edgeMetaData(){
		return ("*Edges     "  + "\n");
	}
	@Override
	protected String formatEdge(int v1, int v2, String weight) {
		//16     21       4.5
		return ((v1) + "\t" + (v2) +(weight==null?"": ("\t"+weight ))+ "\n");
	}
	


	
}
