package io.graph.pajek;

import io.graph.GraphOutputStream;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;

public class PajekGraphWriter<V,E> extends GraphOutputStream<V, E>{
	
	protected String tokenizer ="\t";
	protected String commentIndicator = "*";
	protected boolean labeled = true;
	@Override
	protected String verticesMetaData(Vector<V> vertexes){
		//*Vertices     26
		return (commentIndicator+"Vertices     " + (vertexes.size()) + "\n");
	}
	public PajekGraphWriter() {
		super();
		this.tokenizer ="\t";
		this.labeled = true;
	}

	public PajekGraphWriter(String tokenizer, boolean labeled, String commentIndicator) {
		super();
		this.tokenizer = tokenizer;
		this.labeled = labeled;
		this.commentIndicator = commentIndicator;
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
				res += (labeled?(tokenizer+ att):"") +tokenizer +attValue ;
			}
		}
		return vid + (labeled?(tokenizer + vlabel):"")  + res + "\n";
	}
	//*Edges
	@Override
	protected String edgeMetaData(){
		return (commentIndicator+"Edges     "  + "\n");
	}
	@Override
	protected String formatEdge(String v1, String v2, String weight) {
		//16     21       4.5
		return ((v1) + tokenizer + (v2) +(weight==null?"": (tokenizer+weight ))+ "\n");
	}
	


	
}
