package measure;

import java.util.Vector;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
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

public class MeasuresUtil {
	
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
			res.add(new GraphAGAM<V,E>(graph, weights,Type.X2,ExternalOverlap.NodesWeightedByDegree ,AdjustionMethod.ARI_ADJUSTED));
			res.add(new GraphAGAM<V,E>(graph, weights,Type.X2,ExternalOverlap.Edges ,AdjustionMethod.ARI_ADJUSTED));

//			res.addAll(MeasuresUtil.<V,E>getAlgebricAgreementAlternatives(graph,weights));
		}
		return res;
	}
	
	private static<V> Vector< ClusteringAgreement<V>> getClassicAgreements(){
		Vector< ClusteringAgreement<V>> res = new Vector< ClusteringAgreement<V>>();
		res.add(new Jaccard<V>());
		res.add(new FMeasure<V>());
//		res.add(new ARI<V>());
		res.add(new VI<V>());
		res.add(new NMI<V>());
		res.add(new NMI<V>(NMI.Mode.SQRT));
//		res.add(new AMI<V>());
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




