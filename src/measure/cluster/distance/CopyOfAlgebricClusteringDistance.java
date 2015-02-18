package measure.cluster.distance;

import java.util.Collection;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.test;

import org.apache.commons.collections15.Transformer;
import org.la4j.Matrix;
import org.la4j.matrix.sparse.CRSMatrix;
import org.la4j.Vectors;
import org.la4j.vector.dense.BasicVector;
import org.python.modules.math;

import edu.uci.ics.jung.graph.Graph;


public class CopyOfAlgebricClusteringDistance<d> {
	
	class D{
		double raw;
		double normalized;
		double adjusted;
		public void reverse(){
			raw = 1/raw;
			normalized = 1- normalized;
			adjusted = 1- adjusted;
		}
	}
	
	public D getDistanceFromOverlap(Matrix U, Matrix V){
		return getDistanceFromOverlap(U, V, 
				new Transformer<Double, Double>() {
					@Override
					public Double transform(Double x) {
						return x*(x-1)*0.5;
					}
				});
	}	
	
	public D getDistanceFromOverlap(Matrix U, Matrix V, final Transformer<Double, ? extends Number> f){
		Matrix N = U.transpose().multiply(V);
		
//		b2 = np.ones((N.shape[1],1),dtype='float')
		BasicVector b2 = new BasicVector(N.columns());
		b2.setAll(1);
//	    b1 = np.ones((1, N.shape[0]),dtype='float')
		BasicVector b1 = new BasicVector(N.rows());
		b1.setAll(1);
		
//	    m1=b1.dot(N)
		org.la4j.Vector m1 = b1.multiply(N);
//	    m2=N.dot(b2)
		org.la4j.Vector m2 = N.multiply(b2);
//	    m = b1.dot(N).dot(b2)
//	    s1 = (b1.dot(f(m2))) 
//	    s2 = (f(m1).dot(b2)) 
//	    n = f(m)
//	    I = (b1.dot(f(N)).dot(b2))

		return new D();
	}

	
	public static<V,E> Matrix getAdjecency(Graph<V,E> G,final Transformer<E, ? extends Number> weights ){
		CRSMatrix A = new CRSMatrix(G.getVertexCount(),G.getVertexCount());
		
		for (E e : G.getEdges()) {
			double w = (Double) ( weights==null? 1 : weights.transform(e));

			Vector<V> nodes =  (Vector<V>) G.getIncidentVertices(e);
			//TODO: fix latter: assuming a node is always int and starts from 0!
			A.set((Integer)nodes.get(0), (Integer)nodes.get(1), w );
		} 
		return A;
	}
	public static<V,E> Matrix getIncidence(Graph<V,E> G,final Transformer<E, ? extends Number> weights  ){
		return getIncidence(G, weights, false);
	}

	public static<V,E> Matrix getIncidence(Graph<V,E> G,final Transformer<E, ? extends Number> weights , boolean directed ){
		CRSMatrix N = new CRSMatrix(G.getVertexCount(),G.getEdgeCount());
		
		Vector<E> edges = (Vector<E>) G.getEdges();
		for (int i = 0; i < edges.size(); i++) {
			E e = edges.get(i);
			double w =  ( weights==null? 1 : math.sqrt((Double)weights.transform(e)));
			Vector<V> nodes =  (Vector<V>) G.getIncidentVertices(e);
			//TODO: assuming an edge is always int 
			N.set((Integer)nodes.get(0), i, w  );
			N.set((Integer)nodes.get(1), i, w * (directed?-1:1) );
		} 
		return N;
	}
	
	public static<d> Matrix getClusteringMatrix(int numberOfdatapoints,Vector<Set<d>> clustering){
		CRSMatrix U = new CRSMatrix(numberOfdatapoints,clustering.size());
		for (int k = 0; k < clustering.size(); k++) {
			for (d i : clustering.get(k)) {
				U.set((Integer)i,k , 1);
			}
		}
		
		return U;
	}
	
	
	
	public static void main(String[] params){
		CopyOfAlgebricClusteringDistance<Double> algebricClusteringDistance = new CopyOfAlgebricClusteringDistance<>();
		
		Vector<Set<Integer>> groundTruth = test.createClusteringFromArray(new int[][]{ {0,1,2,3},{4,5,6}});
		Vector<Set<Integer>> p1 =  test.createClusteringFromArray(new int[][]{ {0,1,2},{3,4,5,6}});
		Vector<Set<Integer>> p2 =  test.createClusteringFromArray(new int[][]{ {0},{1,2,3,4,5,6}});

		Matrix U = getClusteringMatrix(7, groundTruth);
		Matrix V = getClusteringMatrix(7, p1);
		Matrix V2 = getClusteringMatrix(7, p2);

		System.out.println(U);
		algebricClusteringDistance.getDistanceFromOverlap(U, V);
		
		
	}
}
