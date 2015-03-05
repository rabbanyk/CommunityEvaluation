package algorithms.communityMining.data;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.topleaders.data.Partitioning;

public  class Grouping <V> {
	Vector<Set<V>> groups;
	
	public Grouping(){
		groups = new Vector<Set<V>>();
	}
	@SuppressWarnings("unchecked")
	public Grouping(Vector<Set<V>> groups){
		this.groups = (Vector<Set<V>>) groups.clone();
	}
	public Set<V> getGroup(int ind){
		return groups.get(ind);
	}
	public Set<V> getLastGroup(){
		return groups.lastElement();
	}
	
	public Set<V> getGroup(V e){
		for (Set<V> cluster : groups) {
			if(cluster.contains(e))
				return cluster;
		}
		return null;
	}
	
	public void removeEmptyGroups(){
		this.removeGroupsSmallerThan(1);
//		Vector<Set<V>> tmp = new Vector<Set<V>>();
//		for (Set<V> set : groups) {
//			if(set.size()>0) tmp.add(set);
//		}
//		groups = tmp;
	}
	public void removeGroupsSmallerThan(int minSize){
		Vector<Set<V>> tmp = new Vector<Set<V>>();
		for (Set<V> set : groups) {
			if(set.size()>=minSize) tmp.add(set);
		}
		groups = tmp;
	}
	public Vector<Set<V>> getGroups() {
		return (Vector<Set<V>> )groups.clone();
	}
	public Set<V> getDataPoints() {
		Set<V> res = new HashSet<V>();
		for (Set<V> set : groups) {
//			System.err.println(set.size());
			res.addAll(set);
		}
		return res;
	}

	public int getNumberOfGroupsOfSizeAtLeast( int minSize){
		int res = 0;
		for (Set<V> set : groups) {
			if (set.size()>=minSize) res++;
		}
		return res;
	}

	public int getNumberOfGroups(){
		return groups.size();
	}

	public boolean addGroup(Set<V> e) {
		return groups.add(e);
	}

	public String toString() {
		String res = "[";
		res +=" Groups("+groups.size()+"): " + groups.toString();
		return res;
	}

	public String getStatistics() {
		String res = " Number of Groups: " + groups.size();
		
		double avgClusterSize = 0,	sdClusterSize = 0, maxClusterSize = 0, minClusterSize =-1;
		
		for (Set<V> group  : groups) {
			avgClusterSize += group.size();
			sdClusterSize += group.size()*group.size();
			if(group.size()>maxClusterSize) maxClusterSize = group.size();
			if(minClusterSize == -1 || minClusterSize > group.size()) minClusterSize = group.size();
		}
		avgClusterSize /= groups.size();
		sdClusterSize /= groups.size();
		
		res += " , size: " + avgClusterSize+" +/- " + Math.sqrt(sdClusterSize - avgClusterSize*avgClusterSize)+ " in [" + minClusterSize +","+maxClusterSize+"]" ;
		
		return res;
	}
	
	
	
	public void saveToFile(String filename) throws IOException{
		FileOutputStream file = new FileOutputStream(filename);

		for (int i =0; i< groups.size() ; i++){
			file.write(("Community " + i + "\n").getBytes());
			
			for (V v : groups.get(i)) {				
				file.write((v+"\n").getBytes());
			}
			file.write(( "\n").getBytes());
		}
		
		file.write(( "\n").getBytes());
		file.close();
	}

	@SuppressWarnings("unchecked")
	public boolean equals(Object o) {
		if (! (o instanceof Grouping)) return false;
		return (new HashSet<Set<V>>(groups)).equals(new HashSet<Set<V>>(((Grouping<V>)o).getGroups()));
	}

	
	public void saveToFile(String filename, Transformer<V,String> vertexTransformer) throws IOException{
		FileOutputStream file = new FileOutputStream(filename);

		for (int i =0; i< groups.size() ; i++){
			file.write(("Community " + i + "\n").getBytes());
			
			for (V v : groups.get(i)) {				
				file.write((vertexTransformer.transform(v)+"\n").getBytes());
			}
			file.write(( "\n").getBytes());
		}
		
		file.write(( "\n").getBytes());
		file.close();
	}
	
	@SuppressWarnings("unchecked")
	public Grouping<V> clone() {
		Grouping<V> partitioning = new Grouping<V>();
		partitioning.groups = (Vector<Set<V>>) groups.clone();
		return partitioning;
	}
	
}
