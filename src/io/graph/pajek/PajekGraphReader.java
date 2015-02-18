package io.graph.pajek;

import io.graph.GraphInputStream;

import java.io.IOException;
import java.io.InputStream;

import org.apache.commons.collections15.Factory;

import edu.uci.ics.jung.graph.Graph;


public class PajekGraphReader<V,E> extends GraphInputStream<V, E> {

	@Override
	public Graph<V, E> readGraph(InputStream inputStream, Factory<V> nodeFactory,
			Factory<E> edgeFactory) throws IOException {
		return null;
	}


	

}
