package algorithms.communityMining.external_methods;

import io.graph.pajek.PajekGraphWriter;
import io.group.GroupingReader;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.data.Grouping;

/**
 * @author reihaneh
 *
 *  Wrapper for InfoMap/Infomod Community mining approach provided here:
 *	http://www.tp.umu.se/~rosvall/code.html
 *	Newer version here: http://www.mapequation.org/code.html
 *
 * 
 * Minor change in Infomod code to print the correctly find the name of input file if not in the same directory
 * string networkName(infile.begin(),infile.begin() + infile.find("."));
 * >>
 * string networkName(infile.begin(),infile.begin() + infile.find_last_of("."));
 *
 */
public class Infomap<V,E> extends CommunityMinerExecutableWrapper<V, E> {
	
	public enum Version{InfoMod, InfoMapOriginalUndir, InfoMapOriginalDir, InfoMap};
	Version version = Version.InfoMap;
	
	String[] execPaths ={"./execs/CM-Overlapping-InfoMap-Rosvall/infomod/infomod", "./execs/CM-Overlapping-InfoMap-Rosvall/infomap_undir/infomap",
		"./execs/CM-Overlapping-InfoMap-Rosvall/infomap_dir/infomap", "./execs/CM-Overlapping-InfoMap-Rosvall/Infomap-0.13.5/Infomap"};
	String execPath = execPaths[version.ordinal()];
	{
		graphWriter = new PajekGraphWriter< V, E>();
		weightFormat = "%s"; 
		extention = ".net";
		extentionWeighted = ".net";
		startIndexId = 1;
//		verbousMode = true;
	}
	public Infomap() {
	}
	
	public Infomap(Version version) {
		this.version = version;
//		weightFormat = (version == Version.InfoMod)? "%.0f":"%s";
		this.execPath = execPaths[version.ordinal()];
	}
	public FileInputStream runCommunityMiningExecutable(Transformer<String, V> vertexTransformer, String network, String path, boolean isWeighted) throws IOException,
	InterruptedException  {
		return findCommunities(vertexTransformer,network,path, 345234,  10); //default values used by authors 

	}	

	public synchronized FileInputStream findCommunities(Transformer<String, V> vertexTransformer, String network, String path, int random_seed, int numberOfAttempts) throws IOException, InterruptedException {
			
			if (version == Version.InfoMap)
				//TODO: the new InfoMap has many new features to set, including  --directed, and gives a hierarchy in map,  
				 runExecCommand(execPath +" " + path + "/" + network + ".net "+path + "/"  + " -N "+  numberOfAttempts + " --clu");
			else runExecCommand(execPath + " "+random_seed+" " + path + "/" + network + ".net "+numberOfAttempts); 

			//Removing files
//			String[] generatedFiles = {".clu", ".map" , ".net" , ".tree" ,"_map.net","_map.vec"};
//			for(String gf: generatedFiles)	runExecCommand("rm " + path + "/" + network + gf);

		return new FileInputStream(path + "/" + network + ".clu");
	}
	{
		groupingReader = new GroupingReader<V>() {
	/*
	 * Output format: start with number of nodes, then clusterId per line 
*Vertices 34
1
2
1
2
1
2
*/
public Grouping<V> readPartitioning(InputStream is,	 Transformer<String,V> vertexTransformer) throws IOException {
		Grouping<V> partitioning =  new Grouping<V>();
		BufferedReader bufferedInputStream = new BufferedReader(new InputStreamReader(is));

		String tmp;
		int nodeId = 1;
		int clusterId = -1;
		
		bufferedInputStream.readLine(); //"*Vertices 34"
		
		while ((tmp = bufferedInputStream.readLine()) != null) if(tmp.length()>0){
			clusterId = Integer.parseInt( tmp.trim());

			while(clusterId > partitioning.getNumberOfGroups())
				partitioning.addGroup(new HashSet<V>());
			
			partitioning.getGroup(clusterId-1).add(vertexTransformer.transform((nodeId++)+""));
		}
		
		bufferedInputStream.close();
		
		return partitioning;
	}};}
	@Override
	public String getName() {
		return version+"";
	}

	@Override
	public String getShortName() {
		return "IM"+version.ordinal();
	}
	
	
	
}
