package algorithms;

import java.util.Vector;
import java.util.stream.IntStream;

import algorithms.communityMining.CommunityMiner;
import algorithms.communityMining.external_methods.*;
import algorithms.communityMining.external_methods.Infomap.Version;
import algorithms.communityMining.topleaders.TopLeaders;

public class AlgorithmUtils {

	public  static  enum Method {
		 BP, Donetti, Louvain, PottsModel, FastModularity, WalkTrap, Infomap
		 , OSLOM, MOSES, COPRA, AGM, LINKCLUST, BIGCLAM
	}//, Infomod //Too costly!!
	public  static < V, E> CommunityMiner<V,E> getCommunititMiner (Method method){
		switch (method) {
//		case BP: return new BP<V,E>();
		case Donetti: return new Donetti<V,E>();
		case FastModularity: return new FastModularity<V,E>();
		case Louvain: return new Louvain<V,E>();
		case PottsModel: return new PottsModel<V,E>();
		case WalkTrap: return new WalkTrap<V,E>();
		case Infomap: return new Infomap<V,E>();
//		case Infomod: return new Infomap<V,E>(Version.InfoMod);

		case OSLOM: return new OSLOM<V,E>();
		case MOSES: return new MOSES<V,E>();
		case COPRA: return new COPRA<V,E>();
		case AGM: return new AGM<V,E>();
		case BIGCLAM: return new BIGCLAM<V,E>();

		default:
			try {
				System.err.println(method.toString());
				return (CommunityMiner< V, E> ) Class.forName("algorithms.communityMining.methods."+method.toString()).newInstance();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return null;
	}
	//variation of algorithm for different parametes
	public   static < V, E> Vector<CommunityMiner< V, E>> getVariant (Method method){
		Vector<CommunityMiner< V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		switch (method) {
//		case Donetti: return new Donetti<V,E>();
//		case FastModularity: return new FastModularity<V,E>();
//		case Louvain: return new Louvain<V,E>();
//		case PottsModel: return new PottsModel<V,E>();
//		case WalkTrap: return new WalkTrap<V,E>();
//		case Infomap: return new Infomap<V,E>();			
		}
		return communityMiners;
	}
	public   static < V, E> Vector<CommunityMiner< V, E>> getAllCommunititMiners (){
		Vector<CommunityMiner< V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		for (Method method: Method.values()){
			communityMiners.add(AlgorithmUtils.<V,E>getCommunititMiner(method));
		}
		return communityMiners;
	}
	public   static < V, E> Vector<CommunityMiner< V, E>> getModularityBasedCommunititMiners (){
		Vector<CommunityMiner< V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		communityMiners.add(AlgorithmUtils.<V,E>getCommunititMiner(Method.Louvain));
//		communityMiners.add(getCommunititMiner(Method.PottsModel));
		communityMiners.add(AlgorithmUtils.<V,E>getCommunititMiner(Method.FastModularity));
		communityMiners.add(AlgorithmUtils.<V,E>getCommunititMiner(Method.BP));

		return communityMiners;
	}

	public   static < V, E>Vector<CommunityMiner< V, E>> getIgraphCommunititMiners (){
		Vector<CommunityMiner<V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.fastgreedy));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.multilevel));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.infomap));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.walktrap));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.spinglass));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.label_propagation));		
		communityMiners.add(new communityMinerIGraph<V,E>(communityMinerIGraph.Method.leading_eigenvector));		
		return communityMiners;
	}
	public  static < V, E> Vector<CommunityMiner< V, E>> getExecCommunititMiners (){
		Vector<CommunityMiner<V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		communityMiners.add(new WalkTrap<V,E>());//
		communityMiners.add(new Infomap<V,E>());//
		communityMiners.add(new FastModularity<V,E>());//w 
		communityMiners.add(new Louvain<V,E>());//w
//		communityMiners.add(new PottsModel<V,E>());//A bit slow
//		communityMiners.add(new Donetti<V,E>()); needs a parameter, otherwise performs poorly

		return communityMiners;
	}
	public  static < V, E> Vector<CommunityMiner< V, E>> getExecOverlappingCommunititMiners (){
		Vector<CommunityMiner<V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();
		communityMiners.add(new OSLOM<V,E>());// 
		communityMiners.add(new MOSES<V,E>());//
		communityMiners.add(new COPRA<V,E>());//
//		communityMiners.add(new AGM<V,E>());// Veeery slow!!!
		communityMiners.add(new BIGCLAM<V,E>());//

		return communityMiners;
	}
	
	public  static < V, E> Vector<CommunityMiner< V, E>> getSelectedCommunititMiners (){
		return getSelectedCommunititMiners(false);
	}
	public  static < V, E> Vector<CommunityMiner< V, E>> getSelectedCommunititMiners (boolean overlapping){
		return getSelectedCommunititMiners(overlapping,-1);
	}
		
	public  static < V, E> Vector<CommunityMiner< V, E>> getSelectedCommunititMiners (boolean overlapping, int maxK){
		return getSelectedCommunititMiners(overlapping, IntStream.range(2, maxK).toArray() );
	}
	
	public  static < V, E> Vector<CommunityMiner< V, E>> getSelectedCommunititMiners (boolean overlapping, int[] Tlks){
		Vector<CommunityMiner<V, E>> communityMiners = new Vector<CommunityMiner<V,E>>();

		if (!overlapping){
			for (int k :Tlks) {
				communityMiners.add(new TopLeaders<V,E>(k));
			}
			communityMiners.addAll(AlgorithmUtils.<V,E>getExecCommunititMiners());
//		communityMiners.add(minerUtils.getCommunititMiner(Method.FastModularity));
//		communityMiners.addAll(minerUtils.getModularityBasedCommunititMiners());

//		communityMiners.addAll(AlgorithmUtils.<V,E>getIgraphCommunititMiners());
		}else{
			communityMiners.addAll(AlgorithmUtils.<V,E>getExecOverlappingCommunititMiners());
		}
		
		return communityMiners;
	}
	
	
	
	
//	public static void main(String[] args){
//		AlgorithmUtils<Integer, Integer> minerUtils = new AlgorithmUtils<Integer, Integer>();
////		CommunityMiner<Integer, Integer> miner = minerUtils.getCommunititMiner(Method.Donetti);
//		for ( CommunityMiner< Integer, Integer> miner : minerUtils.getAllCommunititMiners())
//			System.err.println(miner.getName());
//	}
	
	
}
