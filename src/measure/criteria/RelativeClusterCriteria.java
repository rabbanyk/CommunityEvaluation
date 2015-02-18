package measure.criteria;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;
import measure.base.Centroid;
import measure.base.Proximity;
import measure.base.RelativeCriteria;

public abstract class RelativeClusterCriteria<V> implements RelativeCriteria<V>{
	protected static enum Order { MAXIMIZER, MINIMIZER };
	protected Order order = Order.MAXIMIZER;
	
	public static enum Between{SINGLE, COMPLETE, AVERAGE, CENTROID,SUM} 
	public static enum Within{MAX, AVG, AVG2C, SUM, SUM2C}
	protected Between betweenSetMethod;
	protected Within withinSetMethod;
	
	
	/*
	 * Centrality of a graph or even a community can be calculated in different ways. 
	 */
//	protected Centroid<V> centroidMethod;
	/*
	 * the distance between two data points (i.e nodes) can be calculated in different ways.
	 */
	protected Proximity<V> proximityMeasure;
	
	protected Vector<Set<V>> clusters;
	
	//to cache centroids
	protected Vector<V> centroids;
	//to cache within group distances 
	protected Vector<Double> wd;
	//Number of data points, number of clusters
	protected int N,k;
	
	protected Set<V> allPoints;

	public RelativeClusterCriteria(){		
	}
	
//	public RelativeClusterCriteria(Distance<V> distanceMethod){
//		this.distanceMethod = distanceMethod;
//	}
//	
//	public RelativeClusterCriteria(Distance<V> distanceMethod, Centroid<V> centroidMethod ) {
//		this.distanceMethod = distanceMethod;
//		this.centroidMethod = centroidMethod;
//	}

	
	/**
	 * @param clusters
	 * @return The value of this Criteria over the given communities
	 */
	public double evaluate(Vector<Set<V>> clusters) {
		allPoints = new HashSet<V>();
		for (Set<V> v : clusters) {
			allPoints.addAll(v);
		}
		
		return evaluate(allPoints, clusters);
	}
	
	
	/**
	 * Same as {@link #evaluate(clusters)} but more efficient 
	 * @param allPoints
	 * @param clusters
	 * @return The value of this Criteria over the given communities
	 */
	public double evaluate(Set<V> allPoints, Vector<Set<V>> clusters) {
//		if (clusters.size() <= 1 ) throw new IllegalArgumentException("Number of clusters/communities should be at least 2");
		this.clusters = clusters;
		this.allPoints = allPoints;
		
		k = clusters.size();
		N = allPoints.size();
		
		centroids = null;
		wd = null;
		
		return evaluate();
	}
	
	protected abstract double evaluate();

	/**
	 * @param i
	 * @return cluster with index <tt>i<tt>
	 */
	protected Set<V> getCluster(int i){
		return clusters.get(i);
	}
	
	/**
	 * @param v1
	 * @param v2
	 * @return distance from <tt>v1<tt> to <tt>v2<tt> based on {@link #proximityMeasure}
	 */
	
	protected boolean closnessCompatible = false;
	
	public boolean isSimInUse (){
		return proximityMeasure.isSimilarity() && closnessCompatible;
	}
	protected double getDistance(V v1, V v2){
		return (isSimInUse()?proximityMeasure.getSimilarity(v1, v2):proximityMeasure.getDistance(v1, v2)).doubleValue();
	}

	/**
	 * @param i
	 * @return Centroid of cluster with index <tt>i</tt> 
	 * 
	 */
	protected V getCentroid(int i){
		if(centroids == null) computeAllCentroids();
		return centroids.get(i);
	}
	
	private V computeCentroid(int i){
		return computeCentroid(clusters.get(i));
	}
	
	/**
	 * @param cluster
	 * @return Centroid of the given <tt>cluster</tt>
	 * 
	 */
	protected V computeCentroid(Set<V> X){
		V center = null;
		double sum = 0, sumd; 
		for (V x : X) {
			sumd = 0;
			for (V y : X) if(x!=y){
				sumd += getDistance(y, x);//proximityMeasure.getProximity(y, x).doubleValue();
			}
			if(center == null || (isSimInUse()?(sum < sumd):(sum > sumd))){
				center = x;
				sum = sumd;
			}
		}

		if (center == null){
			System.err.println(" Medoid Null Exception: "+ X);
		}
		return center;
			
//		V centroid = centroidMethod.findCentroid(cluster);
//		if (centroid == null){
//			System.err.println(cluster);
//			throw new CentroidFailedExecption(centroidMethod);
//		}
//			return centroid;
	}

	/**
	 * @return centroid of all data points
	 */
	protected V getCentroid(){
		return computeCentroid(allPoints);
	}
	
	
	private void computeAllCentroids(){
		centroids = new Vector<V>();
		for (int i = 0; i < clusters.size(); i++) {
			centroids.add(computeCentroid(i));
		}
	}
	
	/**
	 * computes all the within scores at once for computational efficiency 
	 */
	private void computeAllWithinClassScores(){
		wd = new Vector<Double>();
		for (int i = 0; i < clusters.size(); i++) {
			wd.add(computeWithin (i));
		}
	}
	
	protected double getWithin(int i){
		if(wd == null) computeAllWithinClassScores();
		return wd.get(i);
	}
	
	/**
	 * @param v
	 * @param i
	 * @return distance of node <tt>v<tt> to cluster with index <tt>i<tt> according to {@link #betweenSetMethod}
	 */
	protected double getDistance(V v, int i){
		double d;
		Double res = null;
		
		switch(betweenSetMethod){
		case CENTROID:
			res = getDistance(v, getCentroid(i));
			break;
		case SINGLE:
			for (V v2: clusters.get(i)) {
					d = getDistance(v, v2);
					if(res == null || (isSimInUse()? (res < d):(res > d)) ) res = d;
				}
			break;
		case COMPLETE:
				for (V v2: clusters.get(i)) {
					d = getDistance(v, v2);
					if(res==null || (isSimInUse()? (res > d):(res < d)) ) res = d;
				}
			break;
		case AVERAGE:
			res = 0.;
				for (V v2: clusters.get(i)) {
					res += getDistance(v, v2);
				}
			res /= N(i) ;
			break;
		}	

		return res;
	}
	
	//get  
	/**
	 * @param i
	 * @param j
	 * @return distance between Cluster <tt>i</tt> and Cluster <tt>j</tt> according to {@link #betweenSetMethod}
	 */
	protected double getBetween(int i, int j){
		double d;
		Double res = null;
		
		switch(betweenSetMethod){
		case CENTROID:
			res = getDistance(getCentroid(i), getCentroid(j));
			break;
		case SINGLE:
//			d = Double.MAX_VALUE;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(j)) {
					d = getDistance(v1, v2);
					if(res == null ||  (isSimInUse()? (res < d):(res > d))) res = d;
				}
			break;
		case COMPLETE:
//			d = Double.MIN_VALUE;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(j)) {
					d = getDistance(v1, v2);
					if(res==null || (isSimInUse()? (res > d):(res < d))) res = d;
				}
			break;
		case SUM:
			res = 0.;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(j)) {
					res += getDistance(v1, v2);
				}
			break;
		case AVERAGE:
			res = 0.;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(j)) {
					res += getDistance(v1, v2);
				}
			res /= N(i) * N(j);
			break;
		}	

		return res;
	}
	
	/**
	 * @param i
	 * @return distance between elements within cluster <tt>i</tt> according to {@link #withinSetMethod}
	 */
	protected double computeWithin(int i){
		if(N(i) == 1) return getDistance(getCentroid(i),getCentroid(i));

		double d;
		Double res = null;		 
		
		switch (withinSetMethod) {
		case MAX:
			res = 0.;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(i)) 
					if(v1!=v2)
					{
						d = getDistance(v1, v2);
						if(res ==null || (isSimInUse()? (res < d):(res > d))) res = d; 
					}
			break;
		case SUM:
			res = 0.;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(i)) 
					if(v1!=v2)
					{
						d = getDistance(v1, v2);
						res += d; 
					}
			break;
		case AVG:
			res = 0.;
			for (V v1: clusters.get(i)) 
				for (V v2: clusters.get(i)) 
					if(v1!=v2)
					{
						d = getDistance(v1, v2);
						res += d; 
					}
			if(N(i)>1)
				res /= N(i)*(N(i)-1);	
			break;
		case AVG2C:
			res = 0.;
			for (V v1: clusters.get(i)) 
				if(v1!=  getCentroid(i))
				{
					res += getDistance(v1, getCentroid(i));
				}
			res /= N(i); 
			break;
		case SUM2C:
			res = 0.;
			for (V v1: clusters.get(i)) 
				if(v1!=  getCentroid(i))
				{
					res += getDistance(v1, getCentroid(i));
				}
			break;
		}
		return res;
	}
	
	
	/**
	 * @param i
	 * @return size of cluster <tt>i</tt>
	 */
	protected double N(int i){
		return clusters.get(i).size();
	}
	
//	/**
//	 * @return whether this criteria is defined for a global optimization or for computing difference/improvement between two clusterings  
//	 */
//	public boolean isDifferenceLike() {
//		return (type == Type.DIFFERENCE_LIKE);
//	}
//	
	/**
	 * @return whether higher values of this criteria is better 
	 */
	public boolean isMaximizer() {
		return (isSimInUse()? (order==Order.MINIMIZER):(order==Order.MAXIMIZER));// ): (order == Order.MAXIMIZER);
	}

	
	public RelativeClusterCriteria<V> setMetrics(Proximity<V> distanceMethod, Centroid< V> centroidMethod){
		setCentroidMethod(centroidMethod);
		setDistanceMethod(distanceMethod);
		return this;
	}
	/**
	 * @param centroidMethod
	 * sets the method that this criteria should use for calculating centroid of clusters
	 */
	public void setCentroidMethod(Centroid< V> centroidMethod){
	//	this.centroidMethod = centroidMethod;
		centroids = null;
		if(withinSetMethod == Within.SUM2C || withinSetMethod == Within.AVG2C) wd = null;
	}
	
	/**
	 * @param distanceMethod
	 * 	sets the method that this criteria should use for calculating distances between data points
	 */
	public void setDistanceMethod(Proximity<V> distanceMethod){
		this.proximityMeasure = distanceMethod;
		wd = null;
	}
//	
	
//	/**
//	 * @return the method that this criteria is using for calculating distances between data points
//	 */
//	public Distance<V> getDistanceMethod(){
//		return this.distanceMethod;
//	}
	
//	/**
//	 * @return the name of the method that this criteria is using for calculating distances between data points
//	 */
//	public String getDistanceStr(){
//		return this.distanceMethod.toString();
//	}
	
//	/**
//	 * @return the name of the method that this criteria is using for calculating calculating centroid of clusters
//	 */
//	public String getCentralityStr(){
//		return this.centroidMethod.toString();		
//	}
	
	public String toString(){
		return proximityMeasure.toString() ;//+ " " + centroidMethod;
	}
	
	
}
