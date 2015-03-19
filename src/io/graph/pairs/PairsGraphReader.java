package io.graph.pairs;

import io.graph.GraphInputStream;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.SparseMultigraph;


public  class PairsGraphReader<V,E> extends GraphInputStream<V, E>{
	{
		commentIndicator = "*";
		tokenizationPattern = "[\\s]+";
	}

	@Override
	protected void parse(String tmp) {
		String[] vals = tmp.split(tokenizationPattern);
		
		if (vals.length<2) return;
		
		V v1= getAddVertex(vals[0]);
		V v2 = getAddVertex(vals[1]);
	
		E e =getAddEdge(v1, v2, vals.length>=3? vals[2]: null);
	}

}
