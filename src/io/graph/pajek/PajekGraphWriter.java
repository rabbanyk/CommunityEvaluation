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
	protected String formatVetice(V v1,Transformer<V, Integer> vertex_Ids,Map<V, HashMap<Object,Vector<Object>>> vertex_attributes){
		String res = "";
		int vid = (vertex_Ids!=null?vertex_Ids.transform(v1):v1.hashCode());
		res+= vid + "\t\"" + vid + "\"\n";
		if(vertex_attributes!=null){
			Map<Object, Vector<Object>> attributes = vertex_attributes.get(v1);
			for(Object att : attributes.keySet()){
				res += "\t"+ att;
				Vector<Object> values = attributes.get(att);
				if (values.size()>1) res +=" {";
				for (Object value :values ){
					String marker = (value instanceof String)?"\"":"";
					res += ((values.size()>1)?",":" ") +  marker + value  +marker;
				}
				if (values.size()>1) res +="} ";
			}
		}
		res +="\n";
 		return res; 
 		
//		return (vertex_Ids.transform(v1) +"\t\"" + (vertex_Ids.transform(v1)) + "\"\n");
	
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
