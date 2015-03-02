package algorithms.dev_topleaders;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

public class Partitioning <V> {

	Vector<Set<V>> communities;
	Set<V> hubs;
	Set<V> outliers;
//	ArrayList<V> centers;
//	
//	
//
//	public ArrayList<V> getCenters() {
//		return centers;
//	}
//
//	public void setCenters(ArrayList<V> centers) {
//		this.centers = centers;
//	}

	public Partitioning(){
		communities = new Vector<Set<V>>();
		hubs = new HashSet<V>();
		outliers = new HashSet<V>();
	}
	
	public Partitioning(Vector<Set<V>> clusters) {
		this.communities = clusters;
		hubs = new HashSet<V>();
		outliers = new HashSet<V>();
	}

	public boolean isHub(V e){
		return hubs.contains(e);
	}
	
	public Set<V> getCluster(int ind){
		return communities.get(ind);
	}
	
	public boolean isOutlier(V e){
		return outliers.contains(e);
	}
	
	//TODO: not efficient for large clusters!!! change it if it is being called in a loop
	public Set<V> getCommunity(V e){
		for (Set<V> cluster : communities) {
			if(cluster.contains(e))
				return cluster;
		}
		if (hubs.contains(e)) return hubs;
		return outliers;
	}
	
	public Vector<Set<V>> getCommunities() {
		return communities;
	}

	public int getNumberOfClusters(){
		return communities.size();
	}
	
	public Set<V> getHubs() {
		return hubs;
	}
	public boolean removeHub(V e){
		 return hubs.remove(e);
	}
	
	public void setHubs(Set<V> hubs) {
		this.hubs = hubs;
	}

	public Set<V> getOutliers() {
		return outliers;
	}

	public boolean addCluster(Set<V> e) {
		return communities.add(e);
	}

	public boolean addOutlier(V e) {
		return outliers.add(e);
	}
	
	public boolean addHub(V e) {
		return hubs.add(e);
	}
	public void setOutliers(Set<V> outliers) {
		this.outliers = outliers;
	}

	@SuppressWarnings("unchecked")
	public boolean equals(Object o) {
		if (! (o instanceof Partitioning)) return false;
		return (new HashSet<Set<V>>(communities)).equals(new HashSet<Set<V>>(((Partitioning<V>)o).getCommunities()));
	}

	public String toString() {
		String res = "[";
		res +=" Communities("+communities.size()+"): " + communities.toString();
		res +="\n Hubs: "+ hubs.toString(); 
		res +="\n Outliers: "+ outliers.toString() + "]";
		return res;
	}

	public String getStatistics() {
		String res = " Number of Clusters: " + communities.size();
		
		double avgClusterSize = 0,	sdClusterSize = 0, maxClusterSize = 0, minClusterSize =-1;
		
		for (Set<V> cluster  : communities) {
			avgClusterSize += cluster.size();
			sdClusterSize += cluster.size()*cluster.size();
			if(cluster.size()>maxClusterSize) maxClusterSize = cluster.size();
			if(minClusterSize == -1 || minClusterSize > cluster.size()) minClusterSize = cluster.size();
		}
		
		avgClusterSize /= communities.size();
		sdClusterSize /= communities.size();
		
		res += " , size: " + avgClusterSize+" +/- " + Math.sqrt(sdClusterSize - avgClusterSize*avgClusterSize)+ " in [" + minClusterSize +","+maxClusterSize+"]" ;
		//res +=" Clusters: " + clusters.toString();
		
		res +=" , number of Hubs: "+ hubs.size(); 
		res +=" , Outliers: "+ outliers.size() ;
		return res;
	}



	public void saveToFile(String filename, Transformer<V,String> vertexTransformer) throws IOException{
		FileOutputStream file = new FileOutputStream(filename);

		for (int i =0; i< communities.size() ; i++){
			file.write(("Community " + i + "\n").getBytes());
			
			for (V v : communities.get(i)) {				
				file.write((vertexTransformer.transform(v)+"\n").getBytes());
			}
			file.write(( "\n").getBytes());
		}
		
		file.write(("Hubs\n").getBytes());
		for (V v : hubs) {
			file.write((vertexTransformer.transform(v)+"\n").getBytes());
		}
		file.write(( "\n").getBytes());
		
		file.write(("Outliers\n").getBytes());
		for (V v : outliers) {
			file.write((vertexTransformer.transform(v).toString()+"\n").getBytes());
		}
		file.write(( "\n").getBytes());
	}
	
	public void saveToFile(String filename) throws IOException{
		FileOutputStream file = new FileOutputStream(filename);

		for (int i =0; i< communities.size() ; i++){
			file.write(("Community " + i + "\n").getBytes());
			
			for (V v : communities.get(i)) {				
				file.write((v+"\n").getBytes());
			}
			file.write(( "\n").getBytes());
		}
		
		file.write(("Hubs\n").getBytes());
		for (V v : hubs) {
			file.write((v+"\n").getBytes());
		}
		file.write(( "\n").getBytes());
		
		file.write(("Outliers\n").getBytes());
		for (V v : outliers) {
			file.write((v.toString()+"\n").getBytes());
		}
		file.write(( "\n").getBytes());
	}

	public Partitioning<V> clone() {
		Partitioning<V> partitioning = new Partitioning<V>();
		
		partitioning.communities = (Vector<Set<V>>) communities.clone();
		partitioning.outliers = (HashSet<V>)((HashSet<V>)outliers).clone();
		partitioning.hubs = (HashSet<V>)((HashSet<V>)hubs).clone();
		
		
		return partitioning;
	}
	

}
