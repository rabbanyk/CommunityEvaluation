package io.group;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.MapTransformer;

import algorithms.communityMining.data.Grouping;

public abstract class GroupingReader<V> {

	public abstract Grouping<V>  readPartitioning(InputStream is,
			Transformer<String, V> vertexTransformer) throws IOException ;

	public  Grouping<V>  readPartitioning(InputStream is,
			Map<String, V> vertexTransformer) throws IOException {
		return readPartitioning(is,   MapTransformer.getInstance(vertexTransformer));
	};
}
