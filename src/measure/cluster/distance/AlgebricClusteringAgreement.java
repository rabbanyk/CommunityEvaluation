package measure.cluster.distance;

import io.Logger;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.MapTransformer;

import algorithms.communityMining.data.Grouping;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.distance.AlgebricClusteringDistance.D;
import data.AbstractMatrix;
import data.AbstractVector;
import data.MatrixVectorZ;
import data.Matrixla4j;
import data.Pair;
import edu.uci.ics.jung.graph.Graph;

public class AlgebricClusteringAgreement<T> extends ClusteringAgreement<T>{
	public static enum AType { 	COMEMEBR_RI, ALT_NORM, ALT_TRACE, OVERLAP_RI, OVERLAP_VI, OMEGA, NMI};
	AType type = AType.COMEMEBR_RI;
	boolean adjusted = true, normalized = true, same_nodes=true;
	
	public static enum StructureType {INDEPENDENT, DEP_TRANS, DEP_SUM };
	Transformer<T,Integer> node_ids ;
	int numberofdatapoints;
	
	public static <M extends AbstractMatrix> M createMatrix(int vertexCount, int edgeCount){
		return (M) new  Matrixla4j(vertexCount,edgeCount);
//		return (M) new MatrixVectorZ(vertexCount, edgeCount);
	}
	public static<T,E> AbstractMatrix getIncidence(Graph<T,E> G,final Transformer<E, ? extends Number> weights , Transformer<T,Integer> node_ids ){
		return getIncidence(G, weights,node_ids, false);
	}
	public static <T,E>  AbstractMatrix getIncidence(Graph<T,E> G,final Transformer<E, ? extends Number> weights , Transformer<T,Integer> node_ids, boolean directed ){
		if (G==null) return null;
		AbstractMatrix N = createMatrix(G.getVertexCount(),G.getEdgeCount());
		Vector<E> edges = new Vector<E> (G.getEdges());
		for (int i = 0; i < edges.size(); i++) {
			E e = edges.get(i);
			double w =  ( weights==null? 1 : Math.sqrt((Double)weights.transform(e)));
			Vector<T> nodes = new Vector<T> ( G.getIncidentVertices(e));
			N.set(node_ids.transform(nodes.get(0)), i, w  );
			N.set(node_ids.transform(nodes.get(1)), i, w * (directed?-1:1) );
		} 
		return N;
	}
	
	public class StructureBasedClusteringAgreement<E>  extends ClusteringAgreement<T>{
		AbstractMatrix N;
		int n; //number of data points
//		SparseMatrix NT;
		StructureType structureType;
		
		public StructureBasedClusteringAgreement(Graph<T, E> G,  final Transformer<E, ? extends Number> weights ){
			this(StructureType.DEP_TRANS, G, weights);
		}
		public StructureBasedClusteringAgreement(StructureType structureType, Graph<T, E> G,  final Transformer<E, ? extends Number> weights ){
			this.structureType = structureType;
			this.n = G!=null?G.getVertexCount():0;
			if (G!=null && structureType!=StructureType.INDEPENDENT){
				node_ids = getNodeIds(G.getVertices());
				numberofdatapoints = G.getVertices().size();
				N = getIncidence(G, weights, node_ids);
//				System.err.println(">>>>>> N:    "+ N.rows()+"x"+N.columns()+" : " + N.getNonZeroCount());
//				NT = N.transpose();
			}
		}
		
		@Override
		public double getAgreement(Vector<Set<T>> U, Vector<Set<T>> V) {
				return getAgreement(getClusteringMatrix(n, U, node_ids), 
					getClusteringMatrix(n ,V,node_ids));
		}
		public double getAgreement(AbstractMatrix U, AbstractMatrix V) {
			switch (structureType) {
			case INDEPENDENT:
				return AlgebricClusteringAgreement.this.getAgreement(U,V);
			case DEP_TRANS:
				return AlgebricClusteringAgreement.this.getAgreement(N.transposeProduct(U), N.transposeProduct(V));
			case DEP_SUM:
				D d0 = AlgebricClusteringAgreement.this.getDistance(U,V);
				D d1 = AlgebricClusteringAgreement.this.getDistance(U,N);
				D d2 = AlgebricClusteringAgreement.this.getDistance(V,N);
	            double alpha =0.5;
	            return AlgebricClusteringAgreement.this.toAgreement(
	            		d0.multiply(alpha).add(d1.absSub(d2).multiply(1-alpha)));
			default:
				return AlgebricClusteringAgreement.this.getAgreement(N.transposeProduct(U), N.transposeProduct(V));	
			}
		}

		
		public double getAgreement(Vector<Set<T>> U) {
				return getAgreement(getClusteringMatrix(N.rows(),U,node_ids));
				
		}
		public double getAgreement(AbstractMatrix U) {
			return AlgebricClusteringAgreement.this.getAgreement(U,N);
			
	}
		public String toString(){
			return structureType+"_" +AlgebricClusteringAgreement.this.toString();
		}
		public String toLatexString(){
			String res ="";
			switch (structureType) {
			case DEP_SUM:
				res+="_{+}";
				break;
			case DEP_TRANS:
				res+="_{\\bot}";
				break;
			default:
				break;
			}
			res+=AlgebricClusteringAgreement.this.toLatexString();
			return res;
		}
		
		
	}
	
	
	public static void addAll(Pair<String[], double[]> tmp ,Vector<Double> measures,Vector<String> names,  String prefix, boolean postfix ){
		for (int i = 0; i < tmp.first.length; i++) {
			measures.add(1- tmp.second[i]);
			names.add( prefix+tmp.first[i]+(postfix?"'":""));
		}
	}
	public static<T,E> Pair<Vector<String>,Vector<Double>> getAllAgreements(
			Graph<T, E> G,  final Transformer<E, ? extends Number> weights,Vector<Set<T>> U, Vector<Set<T>> V){
		return getAllAgreements(G, weights , U, V, false, true);
	}
	public static<T,E> Pair<Vector<String>,Vector<Double>> getAllAgreements(
			Graph<T, E> G,  final Transformer<E, ? extends Number> weights,Vector<Set<T>> U, Vector<Set<T>> V,boolean doAlphas , boolean doTrans) {
//		Set<T> allDataPoints =getAllDatapoints(U, V);
//		int numberofdatapoints = allDataPointsPoints.size();
//		Transformer<T, Integer> node_ids =  getNodeIds(allDataPoints);
		Transformer<T, Integer> node_ids = getNodeIds(G!=null? G.getVertices():null);
		int numberofdatapoints = G!=null? G.getVertexCount(): 0;
		AbstractMatrix N = getIncidence(G, weights,node_ids);

		return getAllAgreements(N, getClusteringMatrix(numberofdatapoints,U,node_ids), 
							getClusteringMatrix(numberofdatapoints,V,node_ids), doAlphas, doTrans);
	}
	
	public static<T,E> Pair<Vector<String>,Vector<Double>> getAllAgreements(AbstractMatrix N ,AbstractMatrix U, AbstractMatrix V){
		return getAllAgreements(N, U, V, false , true);
	}
	public static<T,E> Pair<Vector<String>,Vector<Double>> getAllAgreements(
			AbstractMatrix N ,AbstractMatrix U, AbstractMatrix V,boolean doAlphas , boolean doTrans) {
//		int n = G!=null?G.getVertexCount():0;
		Vector<Double> measures = new Vector<Double>();
		Vector<String> names = new Vector<String>();
		for (boolean same_node : new boolean[]{true}){//,false}) {
	//		case INDEPENDENT:
			Pair<String[], double[]> dInd = AlgebricClusteringDistance.getAll(U, V, same_node);
			addAll(dInd , measures, names,  "C",same_node);
	//		case DEP_TRANS:
			if (doTrans){
				if (N!= null){
					long startTime = System.currentTimeMillis();
					AbstractMatrix NU = N.transposeProduct(U);
					System.err.println("\t\t      >>  -- N.transposeProduct(U) finished in: " + (System.currentTimeMillis() - startTime) +" milisecond");
					System.err.println("\t\t      >>  -- res in: " + NU.getNonZeroCount() +" non zeros in a " + NU.rows()+"x"+NU.columns()+" matrix "+ (NU.getNonZeroCount()*100.0) / (NU.rows()*NU.columns()) + "%");

					startTime = System.currentTimeMillis();
					AbstractMatrix NV = N.transposeProduct(V);
					System.err.println("\t\t      >>  -- N.transposeProduct(V) finished in: " + (System.currentTimeMillis() - startTime) +" milisecond");
					System.err.println("\t\t      >>  -- res in: " + NV.getNonZeroCount() +" non zeros in a " + NV.rows()+"x"+NV.columns()+" matrix "+ (NV.getNonZeroCount()*100.0)  / (NV.rows()*NV.columns()) + "%");

					startTime = System.currentTimeMillis();
					Pair<String[], double[]> dTran= 
							AlgebricClusteringDistance.getAll(NU, NV, same_node);
					System.err.println("\t\t      >>  -- TrasposeAgg finished in: " + (System.currentTimeMillis() - startTime) +" milisecond");
	
					addAll( dTran, measures, names,  "C_{\\bot}",same_node);
				}else {
					addAll(AlgebricClusteringDistance.getAll(null, null, same_node)
								, measures, names,  "C_{\\bot}",same_node);

				}
			}
	//		case DEP_SUM:
	//			D d0 = AlgebricClusteringDistance.getAll(U,V);
			Pair<String[], double[]> d1= AlgebricClusteringDistance.getAll(U,N,same_node);
			addAll( d1, measures, names,  "C_{1}",same_node);
			Pair<String[], double[]> d2= AlgebricClusteringDistance.getAll(V,N,same_node);
			addAll( d2, measures, names,  "C_{2}",same_node);
			
			double alpha = 0.5;
			for (int a =0; a<=100; a+=doAlphas?5:110,alpha = a/100.0){
				double[] dsum = new double[d1.second.length];
		        for (int i = 0; i < dInd.second.length; i++) {
		        	dsum[i] = alpha*dInd.second[i] + Math.abs(d1.second[i]-d2.second[i])*(1-alpha);
				}
				addAll(new Pair<String[], double[]>(dInd.first, dsum) , measures, names,  "C_{+\\alpha=" + alpha+"}" ,same_node);
			}
			}
        return new Pair<Vector<String>, Vector<Double>>(names, measures);
	}
	
	public AlgebricClusteringAgreement(AType type){
		this(type, true, true, true);
	}
	public AlgebricClusteringAgreement(AType type, boolean adjusted){
		this(type, adjusted, true, true);
	}
	public AlgebricClusteringAgreement(AType type, boolean adjusted,boolean normalized ){
		this(type, adjusted, normalized, true);
	}
	public AlgebricClusteringAgreement(AType type, boolean adjusted,boolean normalized ,boolean same_nodes){
		this.type = type;
		this.adjusted = adjusted;
		this.normalized = normalized;
		this.same_nodes = same_nodes;	}
	
	private static <T> Transformer<T,Integer> getNodeIds(Collection<T> allDataPoints){
		Transformer<T, Integer> node_ids =null;
//		int numberofdatapoints = allDataPoints.size();
		if (allDataPoints==null) return null;
		Iterator iterator = allDataPoints.iterator(); 
		T t = null;
		if (iterator.hasNext())
			t = (T) iterator.next();
//		System.err.println(allDataPoints.getClass());
		if (t instanceof Integer){
			int min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
			for (T num : allDataPoints) {
				int integer = (Integer)num;
				if (integer<min) min =integer;
				if (integer>max) max =integer;
			}
			final int minVal = min;
			if (max-min +1 ==allDataPoints.size()) // the data itself is inform of ids, so we use that 
				node_ids = new Transformer<T,Integer>(){
					@Override
					public Integer transform(T input) {
						return (Integer)input-minVal;
					}
					};
		}
		// otherwise assign ids to nodes
		if (node_ids ==null){
			HashMap<T, Integer> ids = new HashMap<>(allDataPoints.size());
			int id =0;
			for (T datapoint : allDataPoints) {
				ids.put(datapoint,id++);
			}
			node_ids =  MapTransformer.getInstance(ids);
		}
		return node_ids; 
	}
	public static <T> Set<T> getAllDatapoints(Vector<Set<T>> U, Vector<Set<T>> V){
		HashSet<T> allDataPoints = new HashSet<>();
		if(U!=null)for (Set<T> x : U) allDataPoints.addAll(x);
		if(V!=null)for (Set<T> x : V) allDataPoints.addAll(x);
		return allDataPoints;
	}
	@Override
	public double getAgreement(Vector<Set<T>> U, Vector<Set<T>> V) {
		if (node_ids ==null || numberofdatapoints ==0) {
			Set<T> allDataPoints =getAllDatapoints(U, V);
			numberofdatapoints = allDataPoints.size();
			node_ids =  getNodeIds(allDataPoints);
		}

		return getAgreement(getClusteringMatrix(numberofdatapoints,U,node_ids), 
							getClusteringMatrix(numberofdatapoints,V,node_ids));
	}
	
//	
//	public SparseRowMatrix getClusteringMatrix(Vector<Set<d>> clustering){
//		return getClusteringMatrix(clustering);
//	}
	public static <T> AbstractMatrix getClusteringMatrix(int numberOfdatapoints,Vector<Set<T>> clustering, Transformer<T, Integer> node_ids){
		if(clustering==null) return null;
		AbstractMatrix U = createMatrix(numberOfdatapoints,clustering.size());
		for (int k = 0; k < clustering.size(); k++) {
			for (T i : clustering.get(k)) {
				U.set(node_ids.transform(i) , k , 1);
			}
			}
		return U;
	}
	
	
	protected double toAgreement(D d){
		switch (type) {
		case OMEGA: case NMI:	
			return (adjusted)?d.adjusted: (normalized?d.normalized:1-d.raw) ;		
		case ALT_NORM: 
			adjusted =false;
		default:
			d.reverse();
			return (adjusted)?d.adjusted:d.normalized ;		
		}
	}
	public double getAgreement(AbstractMatrix U, AbstractMatrix V) {
		return toAgreement( getDistance(U, V));
	}
	
	public String toString(){
		String res = "a_";
		switch (type) {
		case NMI:
			res += (adjusted)?"NMI_sqrt": (normalized?"NMI_sum":"VI") ;	
			break;
		case ALT_TRACE:
			res += (adjusted)?"Tr_sqrt": (normalized?"Tr":"Tr_un_normalized") ;	
			break;
		case ALT_NORM:
			res += "Norm";
			break;
		case OMEGA:
			res = ((adjusted)?"A":"")+"OMEGA";
			break;
		default:
			res += ((adjusted)?"A":"" )+type ;		
		}
		res+= (same_nodes&&type!=AType.NMI)?"'":"";
		return res;
	}
	public String toLatexString(){
		String res = "";
		switch (type) {
		case NMI:
			res += (adjusted)?"NMI_\\sqrt{.}^a": (normalized?"NMI_{+}^a":"VI^a") ;	
			break;
		case ALT_TRACE:
			res += (adjusted)?"I_\\sqrt{tr}": (normalized?"I_{tr}":"Tr_un_normalized") ;	
			break;
		case ALT_NORM:
			res += "I_{norm}";
			break;
		default:
			res += ((adjusted)?"A":"" ) ;
			switch (type) {
			case COMEMEBR_RI:
				res+="RI_\\delta";
				break;
			case OVERLAP_RI:
				res+="RI_{o}";
				break;
			case OVERLAP_VI:
				res+="VI_{o}";
				break;
			case OMEGA:
				res+="\\omega";
				break;
			default:
				break;
			}
		}
		res+= (same_nodes&&type!=AType.NMI)?"'":"";
		return res;
	}
	public D getDistance(AbstractMatrix U, AbstractMatrix V) {
		switch (type) {
		case OVERLAP_RI:
			if (same_nodes)
				return AlgebricClusteringDistance.getDistanceFromOverlap(U, V, false, 
						new Transformer<Double, Double>() {
							public Double transform(Double x) {	return x*x;	}  });
			else 
				return AlgebricClusteringDistance.getDistanceFromOverlap(U, V, true);
		case OVERLAP_VI:
			return AlgebricClusteringDistance.getDistanceFromOverlap(U, V, !same_nodes,	
					new Transformer<Double, Double>() {
						@Override
						public Double transform(Double x) {
							return x==0 ? 0 : x*Math.log(x);
						}
			} );
		case COMEMEBR_RI:
			return AlgebricClusteringDistance.getDistanceFromCoMemberships(U, V, same_nodes);
		case OMEGA:
			return AlgebricClusteringDistance.omega(U, V, same_nodes);
		case NMI:
			return AlgebricClusteringDistance.nmi(U, V);
		case ALT_NORM:
			return AlgebricClusteringDistance.getDistanceFromMatrixNorms(U, V, same_nodes);
		case ALT_TRACE:
			return AlgebricClusteringDistance.getDistanceTarceForms(U, V, same_nodes);
		default:
			return AlgebricClusteringDistance.getDistanceFromCoMemberships(U, V, same_nodes);
		}
	}
	
		
	
}