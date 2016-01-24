package io.graph.pajek;

import io.graph.GraphInputStream;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Vector;


import org.apache.commons.collections15.Factory;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;


public class PajekGraphReader<V,E> extends GraphInputStream<V, E> {
	protected boolean labeled ;
	protected String edgeListStartIndicator= "Edges";
	boolean edgeMode  = false;
	
	public PajekGraphReader(){
		commentIndicator = "*";
		tokenizationPattern = "[\\s]+";
		labeled = true;
	}
	
	public PajekGraphReader(boolean labeled, String tokenizationPattern, String commentIndicator) {
		super();
		this.labeled = labeled;
		this.tokenizationPattern = tokenizationPattern;
		this.commentIndicator = commentIndicator;
	}

	@Override
	public void parse(String tmp) {
		if (tmp.startsWith(commentIndicator)){
			if (tmp.contains(edgeListStartIndicator))
				edgeMode =true;
			return;
		}
		if (edgeMode)
			parseEdge(tmp);
		else  
			parseNode(tmp);
	}
	
	@Override
	protected void init() {
		super.init();
		edgeMode = false;
	}
	
	// Assuming  if labeled format is : id "label" att1 att1val att2 att2val ...
	// if not labeled : id att1val att2val ...
	protected void parseNode(String line){
		String[] vals = line.split(tokenizationPattern);
		if (vals.length<1) return;
		V v= getAddVertex(vals[0]);

		HashMap<Object, Vector<Object>> attributes = new HashMap<>();
		
		if (vals.length>1 && labeled){
			attributes.put("label", new Vector<Object>());
			attributes.get("label").add(parseValue(vals[1]));
		}
		
		Object attKey, attVal;
		String[] attVals; 
		for (int ia=labeled?2:1 ; ia +(labeled?1:0)<vals.length;ia+=(labeled?2:1)){
			attKey = labeled?parseValue(vals[ia]):ia;
			attVals = vals[ia+(labeled?1:0)].split("[,{}\\s]+");
//			if (attributes.get(attKey)==null)
			attributes.put(attKey, new Vector<Object>());
			for (String val: attVals){
				attVal = parseValue(val);
				if(attVal.toString().length()>0)
					attributes.get(attKey).add(attVal);
			}
		}
		
		nodeAttributes.put(v, attributes);
	}

	protected void parseEdge(String line){
		String[] vals = line.split(tokenizationPattern);
		if (vals.length<2) return;
		V v1= getAddVertex(vals[0]);
		V v2 = getAddVertex(vals[1]);
		
		getAddEdge(v1, v2, vals.length>=3? vals[2]: null);
		
	}
}
