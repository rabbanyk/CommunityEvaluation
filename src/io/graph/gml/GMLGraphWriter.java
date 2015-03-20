package io.graph.gml;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import javax.swing.text.html.HTMLDocument.HTMLReader.IsindexAction;

import org.apache.commons.collections15.Transformer;
import org.python.modules.errno;

import edu.uci.ics.jung.graph.Graph;
import io.graph.GraphOutputStream;


/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 *
 * http://www.fim.uni-passau.de/en/fim/faculty/chairs/theoretische-informatik/projects.html
 * https://gephi.org/users/supported-graph-formats/gml-format/
 * 
graph
[
  node
  [
   id A
   label "Node A"
  ]
  node
  ...
   edge
  [
   source B
   target A
   label "Edge B to A"
  ]
  edge
  ...
]
 * 
 * 
 */
public class GMLGraphWriter<V,E> extends GraphOutputStream<V, E>{
	String tokenizer ="\t";
	/*
	graph
		[
		...
		]
	 */
	protected String graphMetaData(Graph<V, E> graph){
		return "graph\n[\n";
	}
	protected String graphEnd(){
		return "]\n";
	}
	
/*
 *  node
  [
   id A
   label "Node A"
  ]
 */
	//Attributes with more than one value, will be written as attKey {val1,val2,...}
	protected String formatVetice(V v1,Transformer<V, String> vertex_Ids,Map<V, HashMap<Object,Vector<Object>>> vertex_attributes){
		String vid = ""+(vertex_Ids!=null?vertex_Ids.transform(v1):v1.hashCode());

		String res = "";
		String attKey, attValue;
		if(vertex_attributes!=null){
			Map<Object, Vector<Object>> attributes = vertex_attributes.get(v1);
			for(Object att : attributes.keySet()){
				attKey = att.toString();
				attValue = getAttValueString(attributes.get(att));
				
				if(attKey.equals("id")) continue;// && attValue.trim().equals(vid)) continue;
				
				res += "\n" + tokenizer + att +tokenizer+attValue ;
			}
		}
 		return "node\n[\n"+tokenizer+"id"+tokenizer+ vid + res + "\n]\n";
	
	}
	
	
/*
 edge
  [
   source B
   target A
   label "Edge B to A"
  ]
  */
	@Override
	protected	String formatEdge(String v1, String v2 ,String weight){
		return ("edge\n[\n\tsource "+(v1) + "\n\ttarget " + (v2)+
				(weight==null?"": ("\n\tvalue " +weight)) 	+ "\n]\n");
	}

	
}
