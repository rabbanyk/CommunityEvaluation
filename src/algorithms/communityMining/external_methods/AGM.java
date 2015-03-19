package algorithms.communityMining.external_methods;

import io.graph.GraphOutputStream;
import io.graph.pairs.PairsGraphWriter;
import io.group.ListGrouingReader;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;

public class AGM <V, E> extends CommunityMinerExecutableWrapper<V, E> {
	String exePath ="./execs/CM-Overlapping-AGM-SNAP/agmfit/";
	DecimalFormat df = new DecimalFormat("#.###");
	{
		graphWriter = new AGAMPairsGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".edgelist";
		extentionWeighted = ".edgelist";
		startIndexId = 0;
		groupingReader = new ListGrouingReader<>();
//		verbousMode = true;
	}
	
	@Override
	public FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted)
			throws IOException, InterruptedException {
		
		int c = 0;
		double e = 0 ;// Math.min(1.0/ (n*n), 0.001);// 0.1;
		runExecCommand(exePath+"agmfitmain "
				+ " -i:" +  path + "/" + network + extention
				+ " -l:" +  path + "/" + network + extention+".labels"  
				+ " -o:" +  path + "/" + network +".res" 
				+ " -c:" + c 
				+ " -e:" + df.format(e));
//		agmfitmain -i:football.edgelist -l:football.labels -c:12 -e:0.1
		return new FileInputStream(path + "/" +  network + ".rescmtyvv.txt");
	}
	
	public  class AGAMPairsGraphWriter<V,E> extends PairsGraphWriter<V, E>{

		protected  String  formatVeticeLabel(V v1,Transformer<V, String> vertex_Ids,Map<V, HashMap<Object,Vector<Object>>> vertex_attributes){
			return vertex_Ids.transform(v1) +"\t" + (vertex_Ids.transform(v1)) + "\n";
		}

		@Override
		public void writeGraph(String path, Graph<V, E> graph, Transformer<V, String> vertex_Ids, Transformer<E, ? extends Number> weights,
				String weightsOutputFormat, Map<V, HashMap<Object, Vector<Object>>> vertex_attributes) throws IOException {
			super.writeGraph(path, graph, vertex_Ids, weights, weightsOutputFormat, vertex_attributes);
		
			
			final Vector<V> vertexes = new Vector<V>(graph.getVertices());
			FileOutputStream outputStream = new FileOutputStream(path+".labels");
			for (int i = 0; i < vertexes.size(); i++) {
				V v1 = vertexes.get(i);
				outputStream.write( formatVeticeLabel(v1,vertex_Ids,vertex_attributes).getBytes());
			}
			outputStream.flush();
			outputStream.close();
		
		}

		
		
	}
	
	@Override
	public String getName() {
		return "AGM";
	}
	
	@Override
	public String getShortName() {
		return "AGM";
	}
	
}
