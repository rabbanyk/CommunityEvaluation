package io.graph.gml;

import java.util.HashMap;
import java.util.Map;

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
	protected String formatVetice(V v1,Transformer<V, Integer> vertex_Ids,Map<V, HashMap<Object,Object>> vertex_attributes){
		String res = "node\n[\n\tid "+ (vertex_Ids!=null?vertex_Ids.transform(v1):v1.hashCode());
		if(vertex_attributes!=null){
			Map<Object, Object> attributes = vertex_attributes.get(v1);
			for(Object att : attributes.keySet()){
				Object value = attributes.get(att);
				String tmp = (value instanceof String)?"\"":"";
				res += "\n\t"+ att+" " +  tmp + value  +tmp;
			}
		}
		res +="\n]\n";
 		return res; 
 	}
	
/*
 edge
  [
   source B
   target A
   label "Edge B to A"
  ]
  */
	protected	String formatEdge(int v1, int v2 ,String weight){
		return ("edge\n[\n\tsource "+(v1) + "\n\ttarget " + (v2)+
				(weight==null?"": ("\n\tvalue " +weight)) 	+ "\n]\n");
	}

	
}
