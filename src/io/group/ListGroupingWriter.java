package io.group;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

public class ListGroupingWriter<V> extends GroupingWriter<V>{

	@Override
	public void writeGrouping(String filename, Vector<Set<V>> grouping, Map<V, Integer> vertex_Ids) throws IOException {
		FileOutputStream file = new FileOutputStream(filename);

		for (int i =0; i< grouping.size() ; i++){
			
			for (V v : grouping.get(i)) {				
				file.write(((vertex_Ids!=null?vertex_Ids.get(v):v)+" ").getBytes());
			}
			file.write(( "\n").getBytes());
		}
		
		file.close();
	}

}
