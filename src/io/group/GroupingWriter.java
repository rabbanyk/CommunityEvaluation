package io.group;

import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;

import edu.uci.ics.jung.graph.Graph;
import algorithms.communityMining.data.Grouping;

public abstract class GroupingWriter <V> {

	public void  writeGrouping(String path, Grouping<V> grouping
			, Map<V, String> vertex_Ids )  throws IOException{
			writeGrouping(path, grouping.getGroups(), vertex_Ids);
	}

//	public void writeGraph(String path, Vector<Set<V> >  grouping,
//			final Map<V, Integer> vertex_Ids) throws IOException {
//			writeGrouping(path, grouping , 	vertex_Ids==null?null: TransformerUtils.mapTransformer(vertex_Ids));
//	}
	public abstract void  writeGrouping(String path, Vector<Set<V> > grouping
			, Map<V, String> vertex_Ids ) throws IOException ;

}
