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
	{
		commentIndicator = "*";
		tokenizationPattern = "[\\s]+";
	}
	protected String edgeListStartIndicator = "Edges";
	boolean edgeMode = false;
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
//		if (attributes.get(vals[0])==null)
//			attributes.put(vals[0], new Vector<Object>());
//		for (int i=1; i< vals.length; i++)
//			attributes.get(vals[0]).add(parseValue(vals[i].trim()));
//	}
//	Vector<Object> ids = (attributes.get("id"));
//	
//	V v = getAddVertex( ((ids!=null && ids.size()>0)?ids.get(0):null).toString());
	
	
	// Assuming format is : id "label" att1 att1val att2 att2val ...
	protected void parseNode(String line){
		String[] vals = line.split(tokenizationPattern);
		if (vals.length<1) return;
		V v= getAddVertex(vals[0]);

		HashMap<Object, Vector<Object>> attributes = new HashMap<>();
		
		if (vals.length>1){
			attributes.put("label", new Vector<Object>());
			attributes.get("label").add(parseValue(vals[1]));
		}
		
		Object attKey, attVal;
		String[] attVals; 
		for (int ia=2;ia+1<vals.length;ia+=2){
			attKey = parseValue(vals[ia]);
			attVals = vals[ia+1].split("[,{}\\s]+");
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
