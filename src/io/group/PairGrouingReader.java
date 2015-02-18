package io.group;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.data.Grouping;

public class PairGrouingReader<V> extends GroupingReader<V>{
	public boolean ComIdfirst = false;
	public PairGrouingReader() {
		this(false);
	}
	public PairGrouingReader(boolean comIdfirst) {
		super();
		ComIdfirst = comIdfirst;
	}

	@Override
	public Grouping<V> readPartitioning(InputStream is,
			Transformer<String, V> vertexTransformer) throws IOException {
		/*	
		 * nodeId commId
		 */
		Grouping<V> partitioning =  new Grouping<V>();
		BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));

		String tmp;
		int clusterId = -1;
		
		while ((tmp = bufferedInputStream.readLine()) != null) if(tmp.length()>0)
		{
			String[] row =  tmp.split("\\s+");
			
			clusterId = Integer.parseInt(row[ComIdfirst?0:1]);
			tmp = row[ComIdfirst?1:0];
//				clusterId = Integer.parseInt( tmp.split("\\s+")[1]);
//				tmp = tmp.split("\\s+")[0];
			
			
			while(clusterId >= partitioning.getNumberOfGroups())
				partitioning.addGroup(new HashSet<V>());
			
			partitioning.getGroup(clusterId).add(vertexTransformer.transform(tmp));
		}
		
		bufferedInputStream.close();
		
		return partitioning;
	}

}
