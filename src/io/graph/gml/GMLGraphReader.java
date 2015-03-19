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
	
	protected void parseEdge() throws IOException{
		HashMap<String, Object> attributes = new HashMap<>();
		String tmp;
		while ((tmp = br.readLine()) != null && !tmp.contains("]"))if(tmp.length()>0 &&  !tmp.contains("[")){	
				String[] vals = tmp.trim().split("[,\\s]+"); 
				attributes.put(vals[0], parseValue(tmp.substring(tmp.indexOf(vals[0])+vals[0].length()).trim())); 
			}
		
		getAddEdge(labels_vertices.get(attributes.get("source").toString()), 
					labels_vertices.get(attributes.get("target").toString()), 
					attributes.get("value") , attributes.get("id"));

	}
	
	protected void parseNode() throws IOException{
		String tmp;
		HashMap<Object, Vector<Object>> attributes = new HashMap<>();
		while ((tmp = br.readLine()) != null && !tmp.contains("]")) if(tmp.length()>0 &&  !tmp.contains("[")){	
			String[] vals = tmp.trim().split("[,{}\\s]+");
			
			if (attributes.get(vals[0])==null)
				attributes.put(vals[0], new Vector<Object>());
			for (int i=1; i< vals.length; i++)
				attributes.get(vals[0]).add(parseValue(vals[i].trim()));
		}
		Vector<Object> ids = (attributes.get("id"));
		
		V v = getAddVertex( ((ids!=null && ids.size()>0)?ids.get(0):null).toString());
		nodeAttributes.put(v, attributes);
	}

	@Override
	protected void parse(String tmp) throws IOException {
		if (tmp.contains("node")){
			parseNode();
		}else if (tmp.contains("edge")){
			parseEdge();
		}else if(tmp.contains("directed"))
			directed = (tmp.contains("1"))?true:false;		
	}
	
	
}
