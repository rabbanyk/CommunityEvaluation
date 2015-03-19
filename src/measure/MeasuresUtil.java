package measure;

import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import data.GraphDataSet;
import edu.uci.ics.jung.graph.Graph;
import measure.base.Centroid;
import measure.base.Medoid;
import measure.base.Proximity;
import measure.base.RelativeCriteria;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.agreement.NMI_Overlapping_LFR;
import measure.cluster.agreement.NMI_Overlapping_MacDaid;
import measure.cluster.agreement.partitioning.PartiotioningAgreement;
import measure.cluster.agreement.partitioning.classics.AMI;
import measure.cluster.agreement.partitioning.classics.ARI;
import measure.cluster.agreement.partitioning.classics.FMeasure;
import measure.cluster.agreement.partitioning.classics.Jaccard;
import measure.cluster.agreement.partitioning.classics.NMI;
import measure.cluster.agreement.partitioning.classics.VI;
import measure.cluster.agreement.partitioning.generalization.AGAM;
import measure.cluster.agreement.partitioning.generalization.GAM;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM;
import measure.cluster.agreement.partitioning.generalization.GAM.Type;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.AdjustionMethod;
import measure.cluster.agreement.partitioning.generalization.GraphAGAM.ExternalOverlap;
import measure.cluster.distance.AlgebricClusteringAgreement;
import measure.cluster.distance.AlgebricClusteringAgreement.AType;
import measure.cluster.distance.AlgebricClusteringAgreement.StructureType;
import measure.criteria.CIndex;
import measure.criteria.DaviesBouldin;
import measure.criteria.Dunn;
import measure.criteria.PBM;
import measure.criteria.PointBiserial;
import measure.criteria.Silhouette;
import measure.criteria.VarianceRatio;
import measure.criteria.ZIndex;
import measure.criteria.RelativeClusterCriteria.Between;
import measure.criteria.RelativeClusterCriteria.Within;
import measure.graph.RelativeCommunityCriteria;
import measure.graph.criteria.Modularity;
import measure.graph.distance.Adjacency;
import measure.graph.distance.AdjacencyPearsonCorrelation;
import measure.graph.distance.AdjacencyPearsonOverlap;
import measure.graph.distance.AdjacencyRelation;
import measure.graph.distance.ICloseness;
import measure.graph.distance.IClosenessVariation;
import measure.graph.distance.ModSimilarity;
import measure.graph.distance.NeighborOverlap;
import measure.graph.distance.NeighborOverlapVariation;
import measure.graph.distance.NumberOfPathes;
import measure.graph.distance.ShortestPath;
import measure.graph.distance.TopologicalOverlap;
import measure.graph.distance.ModSimilarity.Norm;

public class MeasuresUtil {
	public static<V,E> Vector< RelativeCriteria<V>> getRelativeAlternatives(Graph<V, E> graph, Transformer<E, Double> weights){
		Vector< RelativeCriteria<V>> criterias = new Vector<>();
		criterias.add(new Modularity<V, E>(graph,weights));
		Vector<Proximity<V>> distances = getGraphDistanceAlternatives(graph, weights);
		for (Proximity<V> distanceMethod: distances ) {
			Centroid<V> centroidMethod = new Medoid<V>(distanceMethod);
			if(distanceMethod.isSimilarity()) {
				criterias.add(new ZIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
//				criterias.add(new CIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
			}else{
				criterias.add(new ZIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
//				criterias.add(new CIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
			}
//			criterias.add(new VarianceRatio<V>(true).setMetrics(distanceMethod, centroidMethod));
//			criterias.add(new Silhouette<V>(Between.CENTROID,false,true).setMetrics(distanceMethod, centroidMethod)); //SSWC
		}
		return criterias;
	}
	public static<V,E> Vector< RelativeCriteria<V>> getAllRelativeAlternatives(Graph<V, E> graph, Transformer<E, Double> weights){
		Vector< RelativeCriteria<V>> criterias = new Vector<>();
		criterias.add(new Modularity<V, E>(graph,weights));
		Vector<Proximity<V>> distances = getGraphDistanceAlternatives(graph, weights);
		for (Proximity<V> distanceMethod: distances ) {
			Centroid<V> centroidMethod = new Medoid<V>(distanceMethod);

			criterias.add(new VarianceRatio<V>(false).setMetrics(distanceMethod, centroidMethod));
			if(distanceMethod.isSimilarity())
				criterias.add(new VarianceRatio<V>(true).setMetrics(distanceMethod, centroidMethod));
			
			if(distanceMethod.isSimilarity())
				criterias.add(new DaviesBouldin<V>(true).setMetrics(distanceMethod, centroidMethod));
			criterias.add(new DaviesBouldin<V>(false).setMetrics(distanceMethod, centroidMethod));
//				// Different variation of Dunn
			if(!distanceMethod.isSimilarity())
			for (Between between : Between.values()) {
				for (Within within : Within.values()) {
					if(between!=Between.SUM && within!=Within.SUM &&  within!=Within.SUM2C ){ // ignoring unnecessary ones
							criterias.add(new Dunn<V>(between, within,false).setMetrics(distanceMethod, centroidMethod));
						//Dunn is not sim friendly!!!
					}
				}
			}
//			
			if(distanceMethod.isSimilarity()) 		
				criterias.add(new PBM<V>(true).setMetrics(distanceMethod, centroidMethod));
			criterias.add(new PBM<V>(false).setMetrics(distanceMethod, centroidMethod));
////////	//Different Variation of Silhouette
			criterias.add(new Silhouette<V>(Between.AVERAGE,false,false).setMetrics(distanceMethod, centroidMethod)); //SWC
			criterias.add(new Silhouette<V>(Between.CENTROID,false,false).setMetrics(distanceMethod, centroidMethod)); //SSWC
			if(distanceMethod.isSimilarity()){
				criterias.add(new Silhouette<V>(Between.AVERAGE,false,true).setMetrics(distanceMethod, centroidMethod)); //SWC
				criterias.add(new Silhouette<V>(Between.CENTROID,false,true).setMetrics(distanceMethod, centroidMethod)); //SSWC
			}
			criterias.add(new Silhouette<V>(Between.AVERAGE,true,false).setMetrics(distanceMethod, centroidMethod));//ASWC
			criterias.add(new Silhouette<V>(Between.CENTROID,true,false).setMetrics(distanceMethod, centroidMethod));//ASSWC
			if(distanceMethod.isSimilarity()){
				criterias.add(new Silhouette<V>(Between.AVERAGE,true,true).setMetrics(distanceMethod, centroidMethod));//ASWC
				criterias.add(new Silhouette<V>(Between.CENTROID,true,true).setMetrics(distanceMethod, centroidMethod));//ASSWC
			}
//			
			if(distanceMethod.isSimilarity()) 					
				criterias.add(new PointBiserial<V>(true).setMetrics(distanceMethod, centroidMethod));
			else 
				criterias.add(new PointBiserial<V>(false).setMetrics(distanceMethod, centroidMethod));
			
			if(distanceMethod.isSimilarity())
				criterias.add(new CIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
			else 
				criterias.add(new CIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
//			
			if(distanceMethod.isSimilarity()) 
				criterias.add(new ZIndex<V>(true).setMetrics(distanceMethod, centroidMethod));
			else
				criterias.add(new ZIndex<V>(false).setMetrics(distanceMethod, centroidMethod));
			
		}
		return criterias;
	}
	
	/**
	 * @param graph
	 * @param weights
	 * @return possible alternative distance measures that could compute distance between two nodes in the given graph
	 */
	public static <V,E> Vector<Proximity<V>>  getGraphDistanceAlternatives(Graph<V, E> graph,  Transformer<E,  Double> weights){
		Vector<Proximity<V>> distances = new Vector<Proximity<V>>();
//		distances.add(new Adjacency<V, E>(graph ,weights,true));
//		distances.add(new AdjacencyRelation<V, E>(graph ,weights,true));
		distances.add(new NeighborOverlap<V, E>(graph ,weights,true));
		distances.add(new TopologicalOverlap<V, E>(graph ,weights,true));
//		distances.add(new AdjacencyPearsonCorrelation<V, E>(graph ,weights,true,true)); //Normalized Augmented
//		distances.add(new AdjacencyPearsonOverlap<V, E>(graph ,weights)); //Normalized Augmented
//		distances.add(new ShortestPath<V, E>(graph ,weights));
//		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.UNI,3,true));
//		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.LIN,3,true));
//		distances.add(new NumberOfPathes<V, E>(graph ,weights,NumberOfPathes.Coeff.EXP,3,true));
//		distances.add(new ModSimilarity<V, E>(graph ,weights,true,Norm.DIVIDE));		
//		distances.add(new ModSimilarity<V, E>(graph ,weights,true,Norm.MINUS));	
//		distances.add(new ICloseness<V, E>(graph ,weights)); 
//		distances.add(new ICloseness<V, E>(3, graph ,weights)); 
//		distances.add(new ICloseness<V, E>(1, graph ,weights)); 
//		distances.add(new IClosenessVariation<V, E>(3, graph ,weights)); 
//		System.err.println(" resulted in "+distances.size()+" distance variations ");
		return distances;
	}
	
	public static<V> Vector< ClusteringAgreement<V>> getAgreementAlternatives(){
		return getAgreementAlternatives(null,null);
	}

	public static<V,E> Vector< ClusteringAgreement<V>> getAgreementAlternatives(Graph<V, E> graph, Transformer<E, Double> weights){
		return getAgreementAlternatives(graph, weights, false);
	}
		

	public static<V,E> Vector< ClusteringAgreement<V>> getAgreementAlternatives(Graph<V, E> graph, Transformer<E, Double> weights, boolean overlapping){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		if (overlapping){
			res.addAll(MeasuresUtil.<V>getOverlappingAlternatives());
		}else {
			res.addAll(MeasuresUtil.<V>getClassicAgreements());
			if(graph!=null){
				res.add(new GraphAGAM<V,E>(graph, weights,Type.X2,ExternalOverlap.NodesWeightedByDegree ,AdjustionMethod.ARI_ADJUSTED));
				res.add(new GraphAGAM<V,E>(graph, weights,Type.X2,ExternalOverlap.Edges ,AdjustionMethod.ARI_ADJUSTED));
			}
//			res.addAll(MeasuresUtil.<V,E>getAlgebricAgreementAlternatives(graph,weights));
		}
		return res;
	}
	
	public enum Implementation{ORIGINAL, GAM, ALGEBRIC_OVERLAP, ALGEBRIC_DELTA};
	public static<V,E> Vector< ClusteringAgreement<V>> getDisjointAgreements( GraphDataSet<V, E> dataset){
		boolean doExacts, doRI, doNorm, doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega,
		doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased;
		
		doRI=  doVI= doNMIsum=doJacc= doFM=
		doDegreeWeightedStructureBased= doEdgeStructureBased=true;
		
		doExacts=  doNorm= doTrace=   doNMIsqrt=   doAMI=   doNMIvars= doOmega=
		  doTransStructureBased= doSumStructureBased=false;
		return getAgreements(dataset, false, doExacts, doRI, doNorm, doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega,
				doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased, Implementation.GAM);
		
	}
	public static<V,E> Vector< ClusteringAgreement<V>> getAgreements( GraphDataSet<V, E> dataset, boolean overlapping, Implementation imp ){
		boolean doExacts, doRI, doNorm, doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega,
		doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased;
		doExacts= doRI= doNorm= doTrace= doVI= doNMIsum= doNMIsqrt= doAMI= doJacc= doFM= doNMIvars= doOmega=
		doDegreeWeightedStructureBased= doEdgeStructureBased= doTransStructureBased= doSumStructureBased =true;
		return getAgreements(dataset, overlapping, doExacts, doRI, doNorm, doTrace, doVI, doNMIsum, doNMIsqrt, doAMI, doJacc, doFM, doNMIvars, doOmega,
				doDegreeWeightedStructureBased, doEdgeStructureBased, doTransStructureBased, doSumStructureBased, imp);
		
	}
	public static<V,E> Vector< ClusteringAgreement<V>> getAgreements(GraphDataSet<V, E> dataset, 
			boolean overlapping, boolean doExacts , boolean doRI, boolean doNorm , boolean doTrace ,
			boolean doVI , boolean doNMIsum , boolean doNMIsqrt ,boolean  doAMI,
			boolean doJacc ,boolean  doFM, 
			boolean doNMIvars , boolean doOmega ,
			boolean doDegreeWeightedStructureBased  , 
			boolean doEdgeStructureBased  , boolean doTransStructureBased  , boolean doSumStructureBased  , Implementation imp ){
		if(imp==null) {
			if(overlapping)
				imp =Implementation.ALGEBRIC_DELTA;
			else 
				imp =Implementation.GAM;

		}
//		System.err.println(imp);
		Vector<ClusteringAgreement<V>> measures = new Vector<>(); 
		Vector< AlgebricClusteringAgreement<V>> algebricMeasures = new Vector< AlgebricClusteringAgreement<V>>();
		// RI & ARI
		switch (imp) {
		case ORIGINAL:
			measures.add(new ARI<V>());
			break;
		case GAM:
			measures.add(new GraphAGAM<V,E>(null,null,doExacts?Type.RI :Type.X2,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED));
			if(doRI) measures.add(new GraphAGAM<V,E>(null,null,doExacts?Type.RI :Type.X2,ExternalOverlap.Nodes ,AdjustionMethod.NORMALIZE));
			break;
		case ALGEBRIC_OVERLAP:
			measures.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, true, true , !doExacts));
			if(doRI) measures.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, false, true , !doExacts));
			break;
		case ALGEBRIC_DELTA:
			algebricMeasures.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , !doExacts));
			if(doRI) algebricMeasures.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , !doExacts));
			break;
		}
		if(doNorm && imp==Implementation.ALGEBRIC_DELTA)algebricMeasures.add(new AlgebricClusteringAgreement<V>(AType.ALT_NORM, true, true , !doExacts));
		if(doTrace && imp==Implementation.ALGEBRIC_DELTA)algebricMeasures.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE, true, true , !doExacts));
		
		if(doJacc) measures.add(new Jaccard<V>());
		if(doFM) measures.add(new FMeasure<V>());
	
		
		// VI & NMI_sum
		switch (imp) {
		case ORIGINAL:
			if(doVI) measures.add(new VI<V>());
			if(doNMIsum) measures.add(new NMI<V>());
			break;
		case GAM:
			if(doVI) measures.add(new GraphAGAM<V,E>(null,null,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.NORMALIZE));
			if(doNMIsum) measures.add(new GraphAGAM<V,E>(null,null,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED));
			break;
		case ALGEBRIC_OVERLAP: 
			if(doVI) measures.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, false, true , !doExacts));
			if(doNMIsum) measures.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, true, true , !doExacts));
			break;
		case ALGEBRIC_DELTA: //THESE are not really based on delta, but are alternative algebraic implementations
			if(doVI) measures.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, false , !doExacts));
			if(doNMIsum) measures.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, true , !doExacts));
			if(doNMIsqrt) measures.add(new AlgebricClusteringAgreement<V>(AType.NMI, true, true , !doExacts));
			break;
		}
		if(!(imp== Implementation.ALGEBRIC_DELTA) && doNMIsqrt) measures.add(new NMI<V>(NMI.Mode.SQRT));
		if(doAMI) measures.add(new AMI<V>());
			
		if(doNMIvars){
			measures.add(new NMI_Overlapping_LFR<V>());
			measures.add(new NMI_Overlapping_MacDaid<V>());
		}
		if(doOmega) {
			measures.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, false ));
			measures.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, false ));
		}
		
		measures.addAll(0,algebricMeasures);
		if(dataset==null) dataset = new GraphDataSet<V, E>("Null");
//			return measures;
		
		if(doDegreeWeightedStructureBased)
			measures.add(new GraphAGAM<V,E>(dataset.graph, dataset.getWeightsTransformer(),Type.X2,ExternalOverlap.NodesWeightedByDegree , AdjustionMethod.ARI_ADJUSTED));
		if(doEdgeStructureBased)   
			measures.add(new GraphAGAM<V,E>(dataset.graph, dataset.getWeightsTransformer(),Type.X2,ExternalOverlap.Edges , AdjustionMethod.ARI_ADJUSTED));
	
		if(doTransStructureBased)
			for (AlgebricClusteringAgreement<V> algebricClusteringAgreement :algebricMeasures) 
				measures.add(algebricClusteringAgreement.new StructureBasedClusteringAgreement<E>(StructureType.DEP_TRANS,dataset.graph, dataset.getWeightsTransformer()));
		if(doSumStructureBased) 
			for (AlgebricClusteringAgreement<V> algebricClusteringAgreement :algebricMeasures) 
				measures.add(algebricClusteringAgreement.new StructureBasedClusteringAgreement<E>(StructureType.DEP_SUM,dataset.graph, dataset.getWeightsTransformer()));

		return measures;
	}
	
	
	private static<V> Vector< ClusteringAgreement<V>> getClassicAgreements(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new Jaccard<V>());
		res.add(new FMeasure<V>());
//		res.add(new ARI<V>());
		res.add(new VI<V>());
		res.add(new NMI<V>());
		res.add(new NMI<V>(NMI.Mode.SQRT));
		res.add(new AMI<V>());
		return res;
	}
	
	private static<V,E> Vector< ClusteringAgreement<V>> getAlgebricAgreementAlternatives(Graph<V, E> graph, Transformer<E, Double> weights){
		System.err.println("getAgreementAlternatives");
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		for (boolean same_node : new boolean[]{true,false}) {
			Vector<AlgebricClusteringAgreement<V>> algAgrMeas=  MeasuresUtil.<V>getAlgebricAgreementAlternatives(same_node);
			for (AlgebricClusteringAgreement<V> algebricClusteringAgreement :algAgrMeas) 
				res.add(algebricClusteringAgreement.new StructureBasedClusteringAgreement<E>(StructureType.INDEPENDENT,graph, weights));
			for (AlgebricClusteringAgreement<V> algebricClusteringAgreement :algAgrMeas) 
				res.add(algebricClusteringAgreement.new StructureBasedClusteringAgreement<E>(StructureType.DEP_TRANS,graph, weights));
			for (AlgebricClusteringAgreement<V> algebricClusteringAgreement :algAgrMeas) 
				res.add(algebricClusteringAgreement.new StructureBasedClusteringAgreement<E>(StructureType.DEP_SUM,graph, weights));
		}
		return res;
	}
	
	private static<V> Vector< ClusteringAgreement<V>> getOverlappingAlternatives(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new NMI_Overlapping_LFR<V>());
		res.add(new NMI_Overlapping_MacDaid<V>());
		
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, false ));
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, false ));

//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, true ));
//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_NORM));
//		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,false));
//		
		return res;
	}
	
	

	private static<V> Vector< ClusteringAgreement<V>> getPartitioningAlternatives(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new Jaccard<V>());
		res.add(new FMeasure<V>());
		
//		res.add(new GAM<V>(Type.RI));
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, false, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, false ));
//		
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, false, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, true ));

		
//		res.add(new ARI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, true, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , false));
//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, false ));
//
//		res.add(new AGAM<V>(Type.RI)); 	
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, true, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, true));

		
//		res.add(new VI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, false, true , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, false ));

		
//		res.add(new NMI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, true, false , true));
//		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, true ));

//		res.add(new NMI<V>(NMI.Mode.SQRT));
		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, true, false));

		
//		res.add(new AMI<V>());
		
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_NORM));
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,false));
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,true));
		
		return res;
	}
	
	public static<V> Vector< ClusteringAgreement<V>> getAllAgreementImplementations(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new Jaccard<V>());
		res.add(new FMeasure<V>());
		//ADD THIS \sum_{ij} [eij - dj.dj/2E]

		res.add(new GAM<V>(Type.RI));
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, false, true , false));
		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , false));
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, false ));
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, false, true , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, false, true, true ));

		res.add(new ARI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, true, true , false));
		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , false));
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, false ));
		
		res.add(new AGAM<V>(Type.RI)); 	
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_RI, true, true , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.OMEGA, true, true, true));

		res.add(new VI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, false, true , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, false ));

		res.add(new NMI<V>());
		res.add(new AlgebricClusteringAgreement<V>(AType.OVERLAP_VI, true, false , true));
		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, false, true ));

		res.add(new NMI<V>(NMI.Mode.SQRT));
		res.add(new AlgebricClusteringAgreement<V>(AType.NMI, true, false));
//		res.add(new AMI<V>());
		
		
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_NORM));
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,false));
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE,true));
		
		return res;
	}
	

	private static<V> Vector< AlgebricClusteringAgreement<V>> getAlgebricAgreementAlternatives(boolean same_node){
//		System.err.println("getAgreementAlternatives");
		Vector< AlgebricClusteringAgreement<V>> res = new Vector< AlgebricClusteringAgreement<V>>();
		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , same_node));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, false, true , true));

		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , same_node));
//		res.add(new AlgebricClusteringAgreement<V>(AType.COMEMEBR_RI, true, true , true));

		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_NORM, false, true , same_node));
//		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE, false, true , same_node));
		res.add(new AlgebricClusteringAgreement<V>(AType.ALT_TRACE, true, true , same_node));
	
		return res;
	}
	
	
	private static<V> Vector< ClusteringAgreement<V>> getAGAMAlternatives(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new GAM<V>(Type.RI));
		res.add(new GAM<V>(Type.VI));
		res.add(new AGAM<V>(Type.RI)); 	
		
		return res;
	}
	
	
	private static<V,E> Vector< ClusteringAgreement<V>> getAGAMAlternatives(Graph<V, E> graph, Transformer<E, Double> weights){
		System.err.println("getAgreementAlternatives");
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		//		res.add(new GAM<V>(Type.RI));
		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.Nodes ,AdjustionMethod.NORMALIZE));
		res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED));
//		res.add(new AGAM<V>(Type.RI)); 	
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.RI,ExternalOverlap.Nodes ,AdjustionMethod.NONE));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.NONE));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.NORMALIZE));
//		res.add(new GraphAGAM<V,E>(graph, weights,Type.VI,ExternalOverlap.Nodes ,AdjustionMethod.ARI_ADJUSTED));
		for (ExternalOverlap eta: ExternalOverlap.values() ) if (eta !=ExternalOverlap.Nodes && 
				eta != ExternalOverlap.AdjustedNodes &&
				eta != ExternalOverlap.CommonClusteringCoefficientRatio &&
				eta != ExternalOverlap.CommonTrianglesRatio 
				)
		{
//			if 	(	eta == 	ExternalOverlap.NodesWeightedByClusteringCoefficient||
//					eta == 	ExternalOverlap.CommonClusteringCoefficient||
//							eta == 	ExternalOverlap.EdgesRatio||
//									eta == 	 ExternalOverlap.CommonTrianglesRatio||
//											eta == 	 ExternalOverlap.CommonClusteringCoefficientRatio 
//											||eta == 	  ExternalOverlap.AdjustedNodes
//											||eta == 	  ExternalOverlap.AdjustedEdgesGlobal
//											||eta == 	  ExternalOverlap.AdjustedEdgesInter
//											||eta == 	  ExternalOverlap.AdjustedEdgesUnion
//											){
//				for (Type phi: Type.values() ) if(phi == Type.X2){//NVI//if(eta != ExternalOverlap.AdjustedNodes){ 
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
//				}
////				for (Type phi: Type.values() ) if(phi == Type.NVI){//NVI//if(eta != ExternalOverlap.AdjustedNodes){ 
////					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
////					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
////				}
//			}else {
//				for (Type phi: Type.values() ) if(phi == Type.NVI 	  ){//NVI//if(eta != ExternalOverlap.AdjustedNodes){ 
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
////					if(eta == ExternalOverlap.AdjustedNodes)
////						res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NONE));
//				}
				for (Type phi: Type.values() ) if(phi == Type.X2){//NVI//if(eta != ExternalOverlap.AdjustedNodes){ 
					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
//					if(eta == ExternalOverlap.AdjustedNodes)
//						res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NONE));
				}
//				for (Type phi: Type.values() ) if(phi == Type.VI){//NVI//if(eta != ExternalOverlap.AdjustedNodes){ 
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
////					if(eta == ExternalOverlap.AdjustedNodes)
////						res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NONE));
//				}
//				for (Type phi: Type.values() ) if(phi == Type.RI){//X2){//if(eta != ExternalOverlap.AdjustedNodes){ 
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NORMALIZE));
//					res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.ARI_ADJUSTED));
////					if(eta == ExternalOverlap.AdjustedNodes)
////						res.add(new GraphAGAM<V,E>(graph, weights,phi,eta,AdjustionMethod.NONE));
//				}
//			}
		}
		
		return res;
	}
}




