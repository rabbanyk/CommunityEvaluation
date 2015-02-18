package io.group;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.data.Grouping;

public class ListGrouingReader<V> extends GroupingReader<V>{

	String  ignoreline;
	public ListGrouingReader() {
		
	}
	
	public ListGrouingReader(String ignoreline) {
		super();
		this.ignoreline = ignoreline;
	}


	@Override
	public Grouping<V> readPartitioning(InputStream is,
			Transformer<String, V> vertexTransformer) throws IOException {
		/*	
		 * nodeId nodeId2 nodeId3 ...
		 * nodeId nodeId2 nodeId3 ...
		 */
		Grouping<V> partitioning =  new Grouping<V>();
		BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));

		String tmp;
		int clusterId = 0;
		
		while ((tmp = bufferedInputStream.readLine()) != null) 
			if(tmp.length()>0 && (ignoreline==null || !tmp.startsWith(ignoreline))){
				partitioning.addGroup(new HashSet<V>());
				if (tmp.startsWith("[")) tmp = tmp.substring(1, tmp.length()-1);
	//			System.err.println(tmp);
				for (String vId :tmp.split("[,\\s]+"))
					partitioning.getGroup(clusterId).add(vertexTransformer.transform(vId));
				clusterId++;
		}
		
		bufferedInputStream.close();
		
		return partitioning;
	}

}
