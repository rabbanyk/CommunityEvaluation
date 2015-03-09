package measure.cluster.agreement.partitioning.generalization;

import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import measure.cluster.agreement.partitioning.PartiotioningAgreement;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;

public class GraphAGAM<V,E>  extends PartiotioningAgreement<V> {
	public enum ExternalOverlap{Nodes,
		NodesWeightedByDegree,NodesWeightedByTRIANGLES,NodesWeightedByClusteringCoefficient,
		Edges, CommonTriangles, CommonClusteringCoefficient, 
		EdgesRatio, CommonTrianglesRatio, CommonClusteringCoefficientRatio, 
		AdjustedEdgesGlobal, AdjustedEdgesInter, AdjustedEdgesUnion, 
		AdjustedNodes};		
	ExternalOverlap externalOverlapMode = ExternalOverlap.Nodes;
	GAM.Type externalType;
	public enum AdjustionMethod{ NONE, NORMALIZE, ARI_ADJUSTED };
	AdjustionMethod adjusted = AdjustionMethod.ARI_ADJUSTED;
//	Boolean adjusted = true;
	//	AGAM<V> agam ;
	GAM<V> gam ;

	Transformer<Pair<Set<V>>, Double> clusteringOverlapWeighter;
	

	double maxWeight = 1;
	protected double getWeight(E e, final Transformer<E, ? extends Number> weights ){
		if (e==null) return 0;
		return  (weights==null)? 1: (weights.transform(e).doubleValue()/maxWeight);
	}
	protected double getWeight(V v1, V v2, final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights ){
		if(curGraph.findEdge(v1, v2)==null) return 0;
		return  (weights==null)? 1: (weights.transform(curGraph.findEdge(v1, v2)).doubleValue()/maxWeight);
	}
	
	public GraphAGAM(final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights){
		this(curGraph,weights,GAM.Type.RI , ExternalOverlap.Nodes);
	}
	public GraphAGAM(final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights , GAM.Type type , ExternalOverlap mode,Transformer<Double, Double> phi, AdjustionMethod adjusted  ){
		maxWeight = 1;
		if (weights!=null )
			for (E e : curGraph.getEdges())
				maxWeight = Math.max(maxWeight,weights.transform(e).doubleValue());
				externalOverlapMode = mode;
		externalType = type;
		this.adjusted = adjusted;
	
		clusteringOverlapWeighter = new Transformer<Pair<Set<V>>, Double>() {
			public Double transform(Pair<Set<V>> uv) {
				double res = 0, tmp,tmp2;
				V v1,v2;
//				System.err.println(uv.getFirst()+" , "+uv.getSecond());
				Set<V> cap = new HashSet<V>(uv.getFirst());
				cap.retainAll(uv.getSecond());
//				return (double) cap.size();
				Set<V> cup = new HashSet<V>(uv.getFirst());
				cup.addAll(uv.getSecond());
		
				switch (externalOverlapMode) {
				case  Nodes:
					res+=cap.size();
					break;
				case AdjustedNodes:
					double e = (uv.getFirst().size()*uv.getSecond().size())*1./curGraph.getVertexCount();
					res+= 1+ (cap.size() - e)/(curGraph.getVertexCount() -e);
					break;
				case Edges:
					for (V v: cap) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cap.contains(v1))
								res+= getWeight(e1,weights);
						}
					res+=cap.size();
					break;
				case NodesWeightedByDegree:
					for (V v: cap) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							res+=getWeight(e1,weights);
						}
					res+=cap.size();
					break;
				case NodesWeightedByClusteringCoefficient:
					for (V v: cap) {
						tmp= 0;
						tmp2 = 0;
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							for (E e2: curGraph.getIncidentEdges(v)) if(e1!=e2){
								v2 = curGraph.getOpposite(v, e2);
								tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
										getWeight(e1,weights)*
										getWeight(e2,weights);
								
								tmp2+=getWeight(e1,weights)*
										getWeight(e2,weights);
							}
						}
						tmp+=cap.size();
						tmp2+=cap.size();
						res+=(tmp2==0?0: (tmp/tmp2));
					}
//					res+=cap.size();
					break;
				case CommonClusteringCoefficient:
					tmp= 0;
					tmp2 = 0;
					for (V v: cap) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cap.contains(v1 ))
								for (E e2: curGraph.getIncidentEdges(v))  if(e1!=e2){
									v2 = curGraph.getOpposite(v,e2);
									if(cap.contains(v2)){
										tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
												getWeight(e1,weights)*
												getWeight(e2,weights);
										
										tmp2+=getWeight(e1,weights)*
												getWeight(e2,weights);
									}
								}
						}
					tmp+=cap.size();
					tmp2+=cap.size();
					res+=(tmp2==0?0: (tmp/tmp2));
//					res= (res+1)*cap.size();
//					res+=1;
					break;
				case CommonClusteringCoefficientRatio:
					tmp= 0;
					tmp2 = 0;
					for (V v: cup) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cup.contains(v1 ))
								for (E e2: curGraph.getIncidentEdges(v))  if(e1!=e2){
									v2 = curGraph.getOpposite(v,e2);
									if(cup.contains(v2)){
										if(cap.contains(v) && cap.contains(v1) && cap.contains(v2))
											tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
													getWeight(e1,weights)*
													getWeight(e2,weights);
											
										tmp2+=getWeight(e1,weights)*
												getWeight(e2,weights);
									}
								}}
					tmp+=cap.size();
					tmp2+=cap.size();
					res+=(tmp2==0?0: (tmp/tmp2));//*cap.size();
					break;
				case NodesWeightedByTRIANGLES:
					tmp=0;
					for (V v: cap) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							for (E e2: curGraph.getIncidentEdges(v))  if(e1!=e2){
								v2 = curGraph.getOpposite(v,e2);
									tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
											getWeight(e1,weights)*
											getWeight(e2,weights);
							}
						}
					res+=tmp;
					res+=cap.size();
					break;
				case CommonTriangles:
					tmp=0;
					for (V v: cap) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cap.contains(v1 ))
								for (E e2: curGraph.getIncidentEdges(v))  if(e1!=e2){
									v2 = curGraph.getOpposite(v,e2);
									if(cap.contains(v2))
										tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
												getWeight(e1,weights)*
												getWeight(e2,weights);
									
								}
							}
					res+=tmp;
					res+=cap.size();

					break;
				case CommonTrianglesRatio:
					tmp=0;
					tmp2 = 0;
					for (V v: cup) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cup.contains(v1 ))
								for (E e2: curGraph.getIncidentEdges(v))  if(e1!=e2){
									v2 = curGraph.getOpposite(v,e2);
									if(cup.contains(v2)){
										if(cap.contains(v) && cap.contains(v1) && cap.contains(v2))
											tmp+=getWeight(curGraph.findEdge(v1, v2),weights) * 
													getWeight(e1,weights)*
													getWeight(e2,weights);
										tmp2+=getWeight(curGraph.findEdge(v1, v2),weights) * 
												getWeight(e1,weights)*
												getWeight(e2,weights);
									}
								}
						}
					tmp+=cap.size();
					tmp2+=cup.size();
					res+=(tmp2==0?0: (tmp/tmp2));//*cap.size();
					break;
				
				case EdgesRatio:
					tmp=0;
					tmp2=0;
					for (V v: cup) 
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							if(cup.contains(v1)){
								if(cap.contains(v1) && cap.contains(v))
									tmp+=getWeight(e1,weights);
								tmp2+= getWeight(e1,weights);
							}
						}
					tmp+=cap.size();
					tmp2+=cup.size();
					res+=(tmp2==0?0: (tmp/tmp2));//cap.size();
					break;
				case AdjustedEdgesGlobal:
					for (V v: cap){ 
						tmp=0;
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							tmp+=getWeight(e1,weights);

							if(cap.contains(v1)){
								double tmp3=0, tmp4=0;
								for(E e2: curGraph.getIncidentEdges(v)){ 
									tmp3+= getWeight(e2,weights);
								}
								for(E e2: curGraph.getIncidentEdges(v1)){ 
									tmp4+= getWeight(e2,weights);
								}
								tmp2= Math.sqrt((tmp3+1)*(tmp4+1));
								res+= tmp2==0?0: (getWeight(e1,weights)/tmp2);
							}
						}
						
						res+= 1./ tmp;
					}
//					res*=cap.size();
					break;
				case AdjustedEdgesInter:
					for (V v: cap){ 
						tmp=0;
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							tmp+=getWeight(e1,weights);

							if(cap.contains(v1)){
								double tmp3=0, tmp4=0;
								for(E e2: curGraph.getIncidentEdges(v)){ 
									v2 = curGraph.getOpposite(v,e2);
									if(cap.contains(v2))
									tmp3+= getWeight(e2,weights);
								}
								for(E e2: curGraph.getIncidentEdges(v1)){ 
									v2 = curGraph.getOpposite(v1,e2);
									if(cap.contains(v2))
									tmp4+= getWeight(e2,weights);
								}
								tmp2= Math.sqrt((tmp3+1)*(tmp4+1));
								res+= tmp2==0?0: (getWeight(e1,weights)/tmp2);
							}
						}
						res+= 1./ tmp;
					}

					break;
				case AdjustedEdgesUnion:
					for (V v: cap){ 
						tmp=0;
						for (E e1: curGraph.getIncidentEdges(v)){ 
							v1 = curGraph.getOpposite(v,e1);
							tmp+=getWeight(e1,weights);

							if(cap.contains(v1)){
								double tmp3=0, tmp4=0;
								for(E e2: curGraph.getIncidentEdges(v)){ 
									v2 = curGraph.getOpposite(v,e2);
									if(cup.contains(v2))
									tmp3+= getWeight(e2,weights);
								}
								for(E e2: curGraph.getIncidentEdges(v1)){ 
									v2 = curGraph.getOpposite(v1,e2);
									if(cup.contains(v2))
									tmp4+= getWeight(e2,weights);
								}
								tmp2= Math.sqrt((tmp3+1)*(tmp4+1));
								res+= tmp2==0?0: (getWeight(e1,weights)/tmp2);
							}
						}
						res+= 1./ tmp;
					}
					break;
				default:
					break;
				}
//				if (res>1 || res<0){
//					if 	(externalOverlapMode == ExternalOverlap.NodesWeightedByClusteringCoefficient ||
//							externalOverlapMode == 	ExternalOverlap.CommonClusteringCoefficient||
//									externalOverlapMode == 	ExternalOverlap.EdgesRatio||
//											externalOverlapMode == 	 ExternalOverlap.CommonTrianglesRatio||
//													externalOverlapMode == 	 ExternalOverlap.CommonClusteringCoefficientRatio
//													//||externalOverlapMode == 	  ExternalOverlap.AdjustedNodes
//													)
//					System.err.println(externalOverlapMode+":::: " + res);
//				}
				
				if 	(	externalOverlapMode == 		ExternalOverlap.NodesWeightedByClusteringCoefficient||
						externalOverlapMode == 	ExternalOverlap.CommonClusteringCoefficient||
						externalOverlapMode == 	ExternalOverlap.EdgesRatio||
								externalOverlapMode == 	 ExternalOverlap.CommonTrianglesRatio||
										externalOverlapMode == 	 ExternalOverlap.CommonClusteringCoefficientRatio 
										||externalOverlapMode == 	  ExternalOverlap.AdjustedNodes
										||externalOverlapMode == 	  ExternalOverlap.AdjustedEdgesGlobal
										||externalOverlapMode == 	  ExternalOverlap.AdjustedEdgesInter
										||externalOverlapMode == 	  ExternalOverlap.AdjustedEdgesUnion
										){
					 if (res<0)
						System.err.println(externalOverlapMode+"::::::::::NEGATIVE " + res);
				}else if (res<1&&res!=0)
					System.err.println(externalOverlapMode+"::::::::::SMALLER THAN ONE!! " + res);


				
				return res;//==0?0:res+1;//(res==0 || externalOverlapMode==ExternalOverlap.Nodes)?res:res+1;
			}
		};
		
	//	GAM<Integer> gam = new GAM<Integer>();
	
//		Transformer<Double, Double> phiFunction = null;
		switch (type) {
		case NVI:
			phi = GAM.getNegativeXlogx();
			break;
		case VI:
			phi = GAM.getVIPhi();
			break;
		case RI:
			phi = GAM.getRIPhi();
			break;
		case X2:
			phi = GAM.getX2();
			break;
//		case ELSE:
//			phi = phi;
//			break;
//		default:
//			break;
		}
		
//		if(type==Type.ELSE) 	{
		switch (adjusted) {
		case NONE:
			gam = new AGAMA<V>(phi,clusteringOverlapWeighter);
			break;
		case ARI_ADJUSTED:
			gam = new AGAM<V>(phi,clusteringOverlapWeighter);
			break;
		case NORMALIZE:
			gam = new GAM<V>(phi,clusteringOverlapWeighter);
			break;
		default:
			break;
		}
//			gam = adjusted? new AGAM<V>(phi,clusteringOverlapWeighter):new AGAMA<V>(phi,clusteringOverlapWeighter);
//		}
//		else{
//			switch (adjusted) {
//			case NONE:
//				gam = new AGAMA<V>(  ((type == Type.VI) ? GAM.getVIPhi() : ((type == Type.RI)? GAM.getRIPhi(): GAM.getRIPhi())),clusteringOverlapWeighter);
//				break;
//			case ARI_ADJUSTED:
//				gam = (type == Type.VI) ? new AGAM<V>(GAM.getVIPhi(),clusteringOverlapWeighter):new AGAM<V>(GAM.getRIPhi(),clusteringOverlapWeighter);
//				break;
//			case NORMALIZE:
//				gam = (type == Type.VI) ? new GAM<V>(GAM.getVIPhi(),clusteringOverlapWeighter):new GAM<V>(GAM.getRIPhi(),clusteringOverlapWeighter);
//				break;
//			default:
//				break;
//			}	
//		}
	}
	public GraphAGAM(final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights , GAM.Type type , ExternalOverlap mode,Transformer<Double, Double> phi ){
		this( curGraph, weights ,  type ,  mode, phi ,AdjustionMethod.ARI_ADJUSTED);
	}
	public GraphAGAM(final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights , GAM.Type type , ExternalOverlap mode ){
		this( curGraph, weights ,  type ,  mode, null ,AdjustionMethod.ARI_ADJUSTED);
	}
	public GraphAGAM(final Graph<V,E> curGraph, final Transformer<E, ? extends Number> weights , GAM.Type type , ExternalOverlap mode , AdjustionMethod adjusted){
		this( curGraph, weights ,  type ,  mode, null,adjusted );
	}
	
	public double getAgreement(Vector<Set<V>> partitioning, Vector<Set<V>> groundTruth) {
		return gam.getAgreement(partitioning, groundTruth);
	}
	
	
	public String toLatexString(){
		String tail = (externalOverlapMode)+"" ;//+ "-"+externalType;
		String head ="";
		
		switch (externalType) {
		case RI:
			head = "x(x-1)";
			break;
		case VI:
			head = "xlogx";
			break;
		case X2:
			head = "x^2";
			break;
		case NVI:
			head = "(x+1)log(x+1)";
			break;
		default:
			externalType.toString();
			break;
		}
//		head = (externalType == Type.VI) ? "xlogx":((externalType == Type.RI) ? "x^2" : externalType.toString());
		
		switch (externalOverlapMode) {
		case Edges:
			tail = "\\xi";
			break;
		case Nodes:
			tail = "|\\cap|";
			break;
		case NodesWeightedByDegree:
			tail = "\\Sigma d";
			break;
		case NodesWeightedByClusteringCoefficient:
			tail = "\\Sigma c";
			break;
		case NodesWeightedByTRIANGLES:
			tail = "\\Sigma t";
			break;
	
		case CommonClusteringCoefficient:
			tail = "\\Sigma cc";
			break;
		case CommonTriangles:
			tail = "\\Sigma ct";
			break;
		case EdgesRatio:
			tail = "\\xi\\%";
			break;
		case CommonTrianglesRatio:
			tail = "\\Sigma cc\\%";
			break;
		case CommonClusteringCoefficientRatio:
			tail = "\\Sigma ct\\%";
			break;
		case AdjustedEdgesGlobal:
			tail = "\\xi^\\forall";
			break;
		case AdjustedEdgesInter:
			tail = "\\xi^\\cap";
			break;
		case AdjustedEdgesUnion:
			tail = "\\xi^\\cup";
			break;
		case AdjustedNodes:
			tail = "\\hat{n}";
			break;
		}
		
		String modifier ="";
		switch (adjusted) {
		case ARI_ADJUSTED:
			modifier = "A";
			break;
		case NORMALIZE:
					
			break;
		case NONE:
			modifier = "U";

			break;
		default:
			break;
		}
		
		String euivalentTo ="";
		if (externalOverlapMode == ExternalOverlap.Nodes) {
			if (externalType==GAM.Type.X2) euivalentTo= adjusted==AdjustionMethod.ARI_ADJUSTED?"ARI'":(adjusted==AdjustionMethod.NORMALIZE?"RI'":"");
			if (externalType==GAM.Type.VI) euivalentTo= adjusted==AdjustionMethod.ARI_ADJUSTED?"NMI_{+}":(adjusted==AdjustionMethod.NORMALIZE?"VI":"");
			if (externalType==GAM.Type.RI) euivalentTo= adjusted==AdjustionMethod.ARI_ADJUSTED?"ARI":(adjusted==AdjustionMethod.NORMALIZE?"RI":"");
		}
		String res = euivalentTo+" : "+(modifier)+"R_{"+head+ "}^{" +tail+"}";
		
		return res+" ";
	}
	
	public String toString(){
		String tail = (externalOverlapMode)+"" ;//+ "-"+externalType;
		String head =" ";
		
		head = (externalType == GAM.Type.VI) ? "xlogx":((externalType == GAM.Type.RI) ? "x\u00B2" : externalType.toString());
		
		switch (externalOverlapMode) {
		case Edges:
			tail = "\u2211 A\u1D62\u2C7C";
			break;
		case Nodes:
			tail = "\u22C2";//"\u207F";
			break;
		case NodesWeightedByDegree:
			tail = "\u2211 d\u1D62";
			break;
		case NodesWeightedByClusteringCoefficient:
			tail = "\u2211 cc\u1D62";
			break;
		case NodesWeightedByTRIANGLES:
			tail = "\u2211 \u22D5t\u1D62";
			break;
		default:
			tail = externalOverlapMode.toString();
		}
		String modifier ="";
		switch (adjusted) {
		case ARI_ADJUSTED:
			modifier = "A";
			break;
		case NORMALIZE:
					
			break;
		case NONE:
			modifier = "U";

			break;
		default:
			break;
		}
		
		String euivalentTo ="";
		if (externalOverlapMode == ExternalOverlap.Nodes) {
			if (externalType==GAM.Type.RI) euivalentTo= adjusted==AdjustionMethod.ARI_ADJUSTED?"=ARI":(adjusted==AdjustionMethod.NORMALIZE?"=RI":"");
			//else if(externalType == Type.VI) res+= "\u2243VI";
		}
		
		String res = modifier+"{"+head+ "," +tail+"}"+euivalentTo+"";
		
		
		return res+" ";
	}
	
//	public String toSNAMString(){
//		String tail = (externalOverlapMode)+"" ;//+ "-"+externalType;
//		String head =" ";
//		switch (externalOverlapMode) {
//		case Edges:
////			tail = "\u1D5E";
//			head = "\u03BE";//&#958;	&#x03BE 	&xi;
//			tail = "";
//			break;
//		case Nodes:
//			tail = "ARI";//"\u207F";
//			break;
//		case NodesWeightedByDegree:
//			head = "\u03B7";
//			//tail = "\u1D51, w\u1D62=d\u1D62";
//			tail = "(w\u1D62=d\u1D62)";
//			break;
//		case NodesWeightedByClusteringCoefficient:
//			head = "\u03B7";
//			tail = "(w\u1D62=c\u1D62)";
////			tail = "_w\u01D6=c\u01D6";
//			break;
//		case NodesWeightedByTRIANGLES:
//			tail = "(w\u1D62=t\u1D62)";
//			head = "\u03B7";
////			tail = "_w\u01D6=\u01D6";
//			break;
//		}
//		
//	
//		return head+ "" +tail+" ";
//	}
}
