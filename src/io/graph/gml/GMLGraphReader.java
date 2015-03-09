package io.graph.gml;

import io.graph.GraphInputStream;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import javax.print.attribute.standard.Compression;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;

import com.kenai.jaffl.provider.jffi.NumberUtil;

import algorithms.communityMining.topleaders.data.Partitioning;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * @author Reihaneh Rabbany (rabbanyk@ualberta.ca)
 *
 * limited reader, will not support all inputs
 */

public class GMLGraphReader <V,E> extends GraphInputStream<V, E> {
	boolean directed = false;
	
	private Object parseValue(String v){
		try{Integer d = Integer.parseInt(v);return d;}catch(Exception exception){}
		try{Double d = Double.parseDouble(v);return d;}catch(Exception exception){}
		if(v.startsWith("\"")) return v.substring(1, v.length()-1);
		return v;
	}

	

	@Override
	public Graph<V, E> readGraph(InputStream inputStream, Factory<V> nodeFactory,	Factory<E> edgeFactory) throws IOException {
		Graph<V, E> graph = new SparseMultigraph<V, E>();
		BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
		int edgeId = 1, nodeId=1;
		Object w;
		String tmp;
		while ((tmp = br.readLine()) != null) if(tmp.length()>0){		
			if (tmp.contains("node")){
				HashMap<Object, Vector<Object>> attributes = new HashMap<>();
				while ((tmp = br.readLine()) != null && !tmp.contains("]")) if(tmp.length()>0 &&  !tmp.contains("[")){	
					String[] vals = tmp.trim().split("[,{}\\s]+");
//					attributes.put(vals[0], parseValue(tmp.substring(tmp.indexOf(vals[0])+vals[0].length()).trim())); 
					if (attributes.get(vals[0])==null)
						attributes.put(vals[0], new Vector<Object>());
					for (int i=1; i< vals.length; i++)
						attributes.get(vals[0]).add(vals[i].trim());
				}
				V v;
				if((nodeFactory!= null))
					v = nodeFactory.create();
				else {
					Vector<Object> ids = (attributes.get("id"));	
					v = (V) (ids.size()>1?ids:(ids.size()>0?ids.get(0):nodeId++));
					}
				vertex_labels.put(v, v.toString());
				labels_vertices.put(v.toString(),v);
				nodeAttributes.put(v, attributes);
				graph.addVertex(v);
			}else if (tmp.contains("edge")){
				HashMap<String, Object> attributes = new HashMap<>();
				while ((tmp = br.readLine()) != null && !tmp.contains("]"))if(tmp.length()>0 &&  !tmp.contains("[")){	
//						System.err.println("parsing: "+tmp);
						String[] vals = tmp.trim().split("[,\\s]+"); 
						attributes.put(vals[0], parseValue(tmp.substring(tmp.indexOf(vals[0])+vals[0].length()).trim())); 
					}
				E e = (edgeFactory!= null)?edgeFactory.create():(E)(nodeAttributes.get("id")!=null?nodeAttributes.get("id"):edgeId++);	
				if( (w =  attributes.get("value"))!=null) {
					if (weights ==null) weights = new HashMap<E, Double>();
					weights.put(e, new Double(w.toString()));
				}

				graph.addEdge(e, labels_vertices.get(attributes.get("source").toString()),labels_vertices.get(attributes.get("target").toString()));
		
			}else if(tmp.contains("directed"))
				directed = (tmp.contains("1"))?true:false;
		}
		
		return graph;		
	}
	
	public Partitioning<V> groundTruth;

	
	
}
