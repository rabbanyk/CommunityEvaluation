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
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class GraphDataSet<V,E> {
		public String name;
		public Map<E, Double> weights;
		public Map<V, HashMap<Object, Object>> attributes;
		public Graph<V, E> graph;
		public Map<String, V> labels_vertices;

		public GraphDataSet(String name) {
			this.name = name;
		}
		public Transformer<E, Double> getWeights(){
			return weights==null?null: TransformerUtils.mapTransformer(weights);
		}
		
		//Attribute values are the groupIds
		//Or that each unique value determines a new group
		public Grouping<V> getGrouping(Object attribute, final Object missingValue){
			 return getGrouping(attribute,new Transformer<Object, Integer>() {
					Vector<Object> clusters = new Vector<>();
					@Override
					public Integer transform(Object obj) {
						if (obj.equals(missingValue)||obj.toString().equals("null")) return null;
						if (!clusters.contains(obj)) clusters.add(obj);
						return clusters.indexOf(obj);
					}});
		}
		//Attribute values are the groupIds
		//Or that each unique value determines a new group
		public Grouping<V> getGrouping(Object attribute){
			 return getGrouping(attribute,new Transformer<Object, Integer>() {
					Vector<Object> clusters = new Vector<>();
					@Override
					public Integer transform(Object obj) {
						if (!clusters.contains(obj)) clusters.add(obj);
						return clusters.indexOf(obj);
					}});
		}
		// groupId transformers maps the attribute value to a groupId
		public Grouping<V> getGrouping(Object attribute, Transformer<Object, Integer> groupId){
			Grouping<V> res = new Grouping<>();
			for (V v : attributes.keySet()) {
				if(attributes.get(v).containsKey(attribute)){
//					System.err.println(attributes.get(v));
					Integer gId = groupId.transform(attributes.get(v).get(attribute));
					if (gId != null){ //missing
						while (gId>=res.getNumberOfGroups()) res.addGroup(new HashSet<V>());
						res.getGroup(gId).add(v);
					}//else System.err.println("miss");
				}
			}
			//Not good idea since if two nodes are put into two different cluster, we will lose this info
//			res.removeGroupsSmallerThan(2);
			return res;
		}
		
		
		public Map<V,String> getLabels(){
			Map<V, String> res = new HashMap<>();
			if(labels_vertices!=null)
			for(Map.Entry<String,V> entry : labels_vertices.entrySet()){
				res.put(entry.getValue(), entry.getKey());
			}
			else for(V v :graph.getVertices()){
				res.put(v, v.toString());
			}
			return res;
		}
		
		Boolean isConnected = null;
		public boolean isConnected (Graph<V,E> graph){
			if (isConnected!=null) return isConnected;
			WeakComponentClusterer<V, E> clusterer= new WeakComponentClusterer<V, E>();
			if(clusterer.transform(graph).size()>1) {
				isConnected =false;
			}
			isConnected =true;

			return isConnected;
		}
		public void addPartitioingAttribute (String key,Vector<Set<V>> grouping){
//			System.err.println("---------------adding:::  "+grouping.size());
			key  = key.toLowerCase();
			for (int groupid = 0 ; groupid< grouping.size(); groupid++) {
				for (V v : grouping.get(groupid)) {
//					System.err.println(v);
					if(attributes.get(v)==null)
						attributes.put(v,new HashMap<>());
					attributes.get(v).put(key, groupid+1);
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
