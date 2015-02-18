package io.graph;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.EdgeType;

public abstract class GraphOutputStream<V,E> {

	protected  String  graphMetaData(Graph<V,E> graph){return "";}
	protected  String  graphEnd(){return "";}
	protected  String  verticesMetaData(Vector<V> vertexes){return "";}
	protected  String  formatVetice(V v1,Transformer<V, Integer> vertex_Ids,Map<V, HashMap<Object,Object>> vertex_attributes){return "";}
	protected  String  edgeMetaData(){return "";}
	protected  String  formatEdge(int v1, int v2 ,String weight){return "";}

	public void writeGraph(String path, Graph<V,E> graph, Transformer<V, Integer> vertex_Ids ,
Transformer<E, ? extends Number> weights , String weightsOutputFormat, 
Map<V, HashMap<Object,Object>> vertex_attributes) throws IOException{
		final Vector<V> vertexes = new Vector<V>(graph.getVertices());
		FileOutputStream outputStream = new FileOutputStream(path);
		outputStream.write(graphMetaData(graph).getBytes());
		outputStream.write(verticesMetaData(vertexes).getBytes());
		for (int i = 0; i < vertexes.size(); i++) {
			V v1 = vertexes.get(i);
			outputStream.write( formatVetice(v1,vertex_Ids,vertex_attributes).getBytes());
		}
		
		outputStream.write(edgeMetaData().getBytes());

			for (int i = 0; i < vertexes.size(); i++) {
				Vector<Integer> neighs = new Vector<Integer>();
				Vector<E> neigEdges = new Vector<E>();

				V v1 = vertexes.get(i);
				int id1 = vertex_Ids!=null?vertex_Ids.transform(v1):v1.hashCode();
//				System.err.println(id1);
					
				for (E e : graph.getIncidentEdges(v1)){
					//Sorting neighbors 
					V v2 = graph.getOpposite(v1, e);
					int id2 = vertex_Ids!=null?vertex_Ids.transform(v2):v2.hashCode();
					
					
					if(graph.getEdgeType(e)==EdgeType.UNDIRECTED & id1>id2) continue; //eliminating repeated egdes 4 20 ... 20 4

					int index = 0;
					while (index < neighs.size() && neighs.elementAt(index) < id2)
						index++;
					
					neighs.insertElementAt(id2, index);
					if (weights!=null)
						neigEdges.insertElementAt(e, index);
				}
				
				for (int j = 0; j < neighs.size(); j++)
					outputStream.write(formatEdge(id1, neighs.get(j), 
							weights==null?null:String.format(weightsOutputFormat!=null?weightsOutputFormat:"%s",weights.transform(neigEdges.get(j)))).getBytes());
					//((id1) + "\t" + (neighs.get(j)) +(weights==null?"": "\t"+String.format(weightsOutputFormat,weights.transform(neigEdges.get(j)))) + "\n");
			}
			outputStream.write(graphEnd().getBytes());

			outputStream.flush();
			outputStream.close();
	}
	public void writeGraph(String path, Graph<V,E> graph, Transformer<V, Integer> vertex_Ids ,Transformer<E, ? extends Number> weights , String weightsOutputFormat) throws IOException{
		 writeGraph(path,  graph, vertex_Ids , weights, weightsOutputFormat, null);

	}
	public void writeGraph(String path, Graph<V,E> graph, Transformer<V, Integer> vertex_Ids ,Transformer<E, ? extends Number> weights ) throws IOException{
		 writeGraph(path,  graph, vertex_Ids , weights, "%s");
	}
	public void writeGraph(String path, Graph<V,E> graph ) throws IOException{
		 writeGraph(path,  graph, null );
	}
	public void writeGraph(String filepath, Graph<V, E> graph, final Map<V, Integer> vertex_Ids) throws IOException {
		writeGraph(filepath, graph , vertex_Ids==null?null: TransformerUtils.mapTransformer(vertex_Ids) , null);
	}
	public void writeGraph(String filepath, Graph<V, E> graph, final Map<V, Integer> vertex_Ids, final  Map<E, ? extends Number> weights) throws IOException {
		writeGraph(filepath, graph , vertex_Ids==null?null: TransformerUtils.mapTransformer(vertex_Ids),  weights==null?null:TransformerUtils.mapTransformer(weights));
	}
	public void writeGraph(String filepath, Graph<V, E> graph,
			final Map<V, Integer> vertex_Ids, final  Map<E, ? extends Number> weights,	 String weightsOutputFormat) throws IOException {
		writeGraph(filepath, graph , vertex_Ids==null?null: TransformerUtils.mapTransformer(vertex_Ids), weights==null?null:TransformerUtils.mapTransformer(weights), weightsOutputFormat );
	}
	public void writeGraph(String filepath, Graph<V, E> graph,
			final Map<V, Integer> vertex_Ids, final  Map<E, ? extends Number> weights,	 String weightsOutputFormat, Map<V, HashMap<Object,Object>> vertex_attributes) throws IOException {
		writeGraph(filepath, graph ,vertex_Ids==null?null: TransformerUtils.mapTransformer(vertex_Ids), 
				weights==null?null:TransformerUtils.mapTransformer(weights), weightsOutputFormat ,vertex_attributes);
	}
//	public  void writeGraph(String path, Graph<V,E> graph, Transformer<V, Integer> vertex_Ids ) throws IOException{
//		 writeGraph(path,  graph, vertex_Ids , null);
//	}
	
}
