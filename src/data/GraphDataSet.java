package data;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;
import org.apache.commons.collections15.comparators.ComparableComparator;
import org.apache.commons.collections15.functors.MapTransformer;
import org.python.modules.synchronize;

import algorithms.communityMining.data.Grouping;
import edu.uci.ics.jung.graph.Graph;

public class GraphDataSet<V,E> {
		public String name;
		public Map<E, Double> weights;
		public Map<V, HashMap<Object, Object>> attributes;
		public Graph<V, E> graph;
		
		public GraphDataSet(String name) {
			this.name = name;
		}
		public Transformer<E, Double> getWeights(){
			return weights==null?null: TransformerUtils.mapTransformer(weights);
		}
		
		
		public Grouping<V> getGrouping(Object attribute, Transformer<Object, Integer> groupId){
			Grouping<V> res = new Grouping<>();
			for (V v : attributes.keySet()) {
				if(attributes.get(v).containsKey(attribute)){
//					System.err.println(attributes.get(v));
					int gId = groupId.transform(attributes.get(v).get(attribute));
					while (gId>=res.getNumberOfGroups()) res.addGroup(new HashSet<V>());
					res.getGroup(gId).add(v);
				}
			}
			return res;
		}
		
		public void addPartitioingAttribute (String key,Vector<Set<V>> grouping){
//			System.err.println("---------------adding:::  "+grouping.size());
			for (int groupid = 0 ; groupid< grouping.size(); groupid++) {
				for (V v : grouping.get(groupid)) {
//					System.err.println(v);
					if(attributes.get(v)==null)
						attributes.put(v,new HashMap<>());
					attributes.get(v).put(key, groupid);
				}
			}
		}
		public Map<V, String> getAttMap(Object attribute){
			Map<V, String> res = new HashMap<V, String> ();
			for (V v : attributes.keySet()) {
				if(attributes.get(v).containsKey(attribute)){
					res.put(v, attributes.get(v).get(attribute).toString());
				}
			}
			return res;
		}
		
		public synchronized void printStats(){
			System.out.print(">>>>>>>>>>>>>>>>>> Graph "+name+" : ");
			System.out.print(graph.getVertexCount()+" nodes and " + graph.getEdgeCount() + " Edges");
			if (weights!= null) System.out.print(" weighted ");
			if (attributes!= null) System.out.print(" attributed ");
			System.out.println();
			
//			if(graph.getVertexCount()<500){
//				System.out.println("Graph: " + graph);
//				if (weights!= null) System.out.println("-- weights: " + weights);
//				if (attributes!= null) System.out.println("-- attributes: " + attributes);
//			}
		}
}
