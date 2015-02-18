package measure.cluster.agreement;

import io.group.GroupingWriter;
import io.group.ListGroupingWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.MapTransformer;

import algorithms.communityMining.data.Grouping;

public abstract class ClusteringAgreementExecWrapper<V> extends ClusteringAgreement<V>{
	int startIndexId =0;
	String exePath;
	String infoMarker ;

	@Override
	public synchronized double getAgreement(Vector<Set<V>> U, Vector<Set<V>> V) {
		String name = this.toString(), path = "./temp";
		new File(path).mkdir();
		
		Map<V, Integer> vertexLabels = new HashMap<V, Integer>();
		final Vector<V> vertexes = new Vector<V>();
		for(Set<V> c:U)
			vertexes.addAll(c);
		for(Set<V> c:V)
			vertexes.addAll(c);
		
		for(int i = 0; i< vertexes.size(); i++){
			vertexLabels.put(vertexes.elementAt(i), i+startIndexId);
		}
		
		GroupingWriter<V> groupingWriter = new ListGroupingWriter<>();
		try {
			groupingWriter.writeGrouping(path +"/"+ name +"_U", U, vertexLabels); 
			groupingWriter.writeGrouping(path +"/"+ name +"_V", V, vertexLabels); 
			
			String line = runExecCommandAndGrabFromSTO(exePath + " " + path +"/"+ name +"_U" + " " 
			+  path +"/"+ name +"_V", infoMarker);
			double res= Double.parseDouble(line.substring(infoMarker.length()).trim());
			runExecCommandAndGrabFromSTO("rm -rf " + path, null);
			return res;
		} catch (Exception e) {
			e.printStackTrace();
			return -10;
		}
	}
	
//			double res = runAgreementExecutable( path +"/"+ name +"_U",  path +"/"+ name +"_V");
//	protected double runAgreementExecutable(String U, String V) throws IOException, InterruptedException{
//		runExecCommand(exePath + " " + U + " " + V, V + "res.tmp");
//			double res = -10;
//			BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(V + "res.tmp")));
//			String line = null;
//			while ((line = input.readLine()) != null) {
//				if (line.startsWith(InfoMarker))
//					res = Double.parseDouble(line.substring(InfoMarker.length()).trim());
//			}
//			input.close();
//	}
	
	
	protected synchronized String runExecCommandAndGrabFromSTO(String command, String grabIndicator) throws IOException, InterruptedException{
		String result=null;
		System.out.println(command);
		Process child = Runtime.getRuntime().exec(command);
		Scanner sc = new Scanner(child.getInputStream());    		
		String output ;
		while ( sc.hasNext()) {
			output =sc.nextLine();
			if(output!=null && output.startsWith(grabIndicator))
				result = output;
		}
		child.waitFor();
		sc.close();
		sc = new Scanner(child.getErrorStream());    		
		while (sc.hasNext()){
			output =sc.nextLine();
		}
		sc.close();
		return result;
	}
//	
//	
//	protected synchronized void runExecCommand(String command, String outputPath) throws IOException, InterruptedException{
//		System.out.println(command);
//		Process child = Runtime.getRuntime().exec(command);
//		Scanner sc = new Scanner(child.getInputStream());    		
//		FileOutputStream outputStream = null;
//		String output ;
//		if(outputPath!=null)
//			outputStream = new FileOutputStream(outputPath);
//		while ( sc.hasNext()) {
//			output =sc.nextLine();
//			if(outputPath!=null)
//				outputStream.write((output+"\n").getBytes());
//		}
//		if(outputPath!=null){
//			outputStream.flush();
//			outputStream.close();
//		}
//		child.waitFor();
//		sc.close();
//		sc = new Scanner(child.getErrorStream());    		
//		while (sc.hasNext()){
//			output =sc.nextLine();
//		}
//		sc.close();
//
//	}
//	
//	

}
