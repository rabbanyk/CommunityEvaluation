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
	
	//Attributes with more than one value, will be written as attKey {val1,val2,...}
	@Override
	protected String formatVetice(V v1,Transformer<V, String> vertex_Ids,Map<V, HashMap<Object,Vector<Object>>> vertex_attributes){
		String vid = ""+(vertex_Ids!=null?vertex_Ids.transform(v1):v1.hashCode());
		String vlabel = "\"" +vid+ "\"" ;
		
		String res = "";
		String attKey, attValue;
		if(vertex_attributes!=null){
			Map<Object, Vector<Object>> attributes = vertex_attributes.get(v1);
			for(Object att : attributes.keySet()){
				attKey = att.toString();
				attValue = getAttValueString(attributes.get(att));
				
				if(attKey.equals("id"))// && attValue.trim().equals(vid)) 
					continue;
				if(attKey.equals("label")){
					vlabel = attValue;
					continue;
				}
				res += "\t"+ att +attValue ;
			}
		}
		return vid + "\t" + vlabel  + res + "\n";
	}
	//*Edges
	@Override
	protected String edgeMetaData(){
		return ("*Edges     "  + "\n");
	}
	@Override
	protected String formatEdge(String v1, String v2, String weight) {
		//16     21       4.5
		return ((v1) + "\t" + (v2) +(weight==null?"": ("\t"+weight ))+ "\n");
	}
	


	
}
