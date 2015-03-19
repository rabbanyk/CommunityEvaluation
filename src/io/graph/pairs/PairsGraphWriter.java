package io.graph.pairs;
import io.graph.GraphOutputStream;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;


public  class PairsGraphWriter<V,E> extends GraphOutputStream<V, E>{

	@Override
	protected String formatEdge(String v1, String v2, String weight) {
		return ((v1) + "\t" + (v2) + (weight==null?"": "\t"+weight)+ "\n");
	}

}
