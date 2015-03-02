package measure.cluster.distance;

import java.util.Set;
import java.util.Vector;

import measure.MeasuresUtil;
import measure.cluster.agreement.ClusteringAgreement;
import measure.cluster.agreement.test;

import org.apache.commons.collections15.Transformer;
import org.python.modules.math;

import util.DatasetUtils;
import data.AbstractMatrix;
import data.AbstractVector;
import data.LambertW;
import data.Pair;


public class AlgebricClusteringDistance {
	
	
	static class D{
		public double raw=0;
		public double normalized=0;
		public double adjusted=0;
		public void reverse(){
			raw = raw!=0? (1/raw): (1/0.000001);
			normalized = 1- normalized;
			adjusted = 1- adjusted;
		}
		public D(){}
		public D(double first, double second, double third){
			raw = first;
			normalized = second;
			adjusted = third;
		}
		public D add(D d){
			raw +=d.raw;
			normalized +=d.normalized;
			adjusted += d.adjusted;
			return this;
		}
		public D sub(D d){
			raw -= -d.raw;
			normalized -=d.normalized;
			adjusted -= d.adjusted;
			return this;
		}
		public D abs(){
			raw = Math.abs(raw);
			normalized =Math.abs(normalized);
			adjusted = Math.abs(adjusted);
			return this;
		}
		public D absSub(D d){
			return sub(d).abs();
		}
		public D multiply(double k){
			raw *= k;
			normalized *= k;
			adjusted *= k;
			return this;
		}
		public String toString(){
//			return raw+","+normalized +","+adjusted;
			return normalized +","+adjusted;

		}
	}
	static double log(double x){
		return x==0?0:Math.log(x+1); 
	}
	static double invxlogx(double y){
		//https://books.google.ca/books?id=1aAOdzK3FegC&pg=PA307&lpg=PA307&dq=inverse+of+x+log+x&source=bl&ots=3jOiG7Bsng&sig=7y3PNLHe8PpzIIdtCMgCxOViwec&hl=en&sa=X&ei=diDIVPrUBIOxyQSorIKwAw&ved=0CB0Q6AEwAA#v=onepage&q=inverse%20of%20x%20log%20x&f=false
		return y/LambertW.branch0(y) -1;
	}
	static Transformer<Pair<Double,Double>, Double> x2 = new Transformer<Pair<Double,Double>, Double>() {
		@Override
		public Double transform(Pair<Double, Double> input) {
			return input.first*input.second;
		}
	};
	static Transformer<Pair<Double,Double>, Double> xlogx = new Transformer<Pair<Double,Double>, Double>() {
		@Override
		public Double transform(Pair<Double, Double> input) {
			return input.first*Math.log(input.second);
		}
	};
	
	
	public static Pair<String[],double[]> getAll(AbstractMatrix U, AbstractMatrix V, boolean same_nodes ){
		double[] measures = new double[5];//+5];
		String[] names = new String[]{"RI", "ARI", "I_{Norm}", "I_{Tr}", "I_{\\sqrt{Tr}}"};
//					,"VI", "AVI", "VI_{Norm}", "VI_{Tr}", "VI_{\\sqrt{Tr}}"};
		if(U!=null && V!=null){
			if (U.rows() !=V.rows() ) System.err.println("ERROR==========  U.rows() !=V.rows() ");
			double[] ress = U.pByTSumPowSumDiffPowSum(V, same_nodes);
			//0:max1, 1:sum1, 2:powsu1, 3:max2, 4:sum2, 5:powsum2, 6:powsumDiff12, 
			//7:Tr1, 8:Tr2, 9:timessum12,
			//10:logsu1, 11: logsu2, 12: logsumdiff12, 13: timeslog12
//			for (int i = 0; i < ress.length; i++) {
//				System.err.println(i+": "+ress[i]+", ");
//			}
			double pairCount =  V.rows()*V.rows() - (same_nodes?0 :V.rows());
		   
		    if (ress[6]!=0){
		    	 measures[0]= ress[6]/( pairCount*Math.pow(Math.max(ress[0],ress[3]),2));
		    	 measures[1]= ress[6] /((ress[2]+ress[5]) - 2*(ress[1]*ress[4])/(pairCount));
		    	 measures[2]= Math.sqrt(ress[6]) /((Math.sqrt(ress[2])+Math.sqrt(ress[5])));
		    }else {}
		    
//		    if (ress[12]!=0){
//		    	 measures[5]= ress[12]/( pairCount*Math.pow(Math.max(ress[0],ress[3]),2));
//		    	 measures[6]= ress[12] /((ress[10]+ress[11]) - (ress[1]*log(ress[4])+ress[4]*log(ress[1]))/(pairCount));
//		    	 measures[7]= invxlogx(ress[12]) /((invxlogx(ress[10])+invxlogx(ress[11])));
//		    }
		    
			if(ress[9]!=0){
				if(same_nodes){ //only defined in this case
					measures[3] = ress[9]/ ( ress[7] *ress[8]);
				}
				measures[4] = ress[9]/ Math.sqrt( ress[2] *ress[5] );
			}
			measures[3] = 1 - measures[3];
			measures[4] = 1 - measures[4];
//			
//			if(ress[13]!=0){
//				if(same_nodes){ //only defined in this case
//					measures[8] = ress[13]/ (ress[7] *log(ress[8]) + ress[8] *log(ress[7]));
//				}
//				double tmp = invxlogx(ress[10] *log(ress[11]) + ress[11] *log(ress[10]) );
//				measures[9] = tmp!=0?ress[13]/tmp:0 ;
//			}
//			measures[8] = 1 - measures[8];
//			measures[9] = 1 - measures[9];
			
			//TODO: Investigate relation between Tr & RI:: since   res[6] = res[2]+res[5] -2*res[9] ;
			//TODO: add NMI OMEGA?
	//		return new Pair(new String[]{"COMEMEBR_RI", "Adj_COMEMEBR_RI", "Norm", "Tr", "Tr_sqrt"},measures);
		}
		return new Pair<String[], double[]>( names,	measures);

	}
	
	public static D getDistanceFromOverlap(AbstractMatrix U, AbstractMatrix V){
		return getDistanceFromOverlap(U, V, false);
	}

	public static D getDistanceFromOverlap(AbstractMatrix U, AbstractMatrix V, boolean exact){
		return getDistanceFromOverlap(U, V, exact,
				new Transformer<Double, Double>() {
					@Override
					public Double transform(Double x) {
						return x*(x-1)*0.5;
					}
				});
	}	
	
	public static D getDistanceFromOverlap(AbstractMatrix U, AbstractMatrix V, boolean exact, final Transformer<Double, Double> f){

		AbstractMatrix N = U.transposeProduct(V);
		AbstractVector m1 = N.getElemetMarginalSum(1);// b1.innerProduct(N);
		AbstractVector m2 = N.getElemetMarginalSum(0);//N.innerProduct(b2);

		double m = N.sum();
		double n = f.transform(m); 

		AbstractMatrix m2m1 = m2.outerProduct(m1);
		m2.inPlaceApplyFunction(f);
		double s1 = m2.sum();
		m1.inPlaceApplyFunction(f);
		double s2 = m1.sum();
		N.applyFunctionInPlace(f);
		double I = N.sum();
		
		D res = new D();
	    res.raw = s1+s2- 2* I ;
	    res.normalized = res.raw/n;
	    
	    double exp =0;
	    if (exact){
	    	exp = m2.outerProduct(m1).sum()/n;
	    }else {
	    	m2m1.divideInPlace(m);
	    	m2m1.applyFunctionInPlace(f);
	    	exp = m2m1.sum();
	    }
	    
	    res.adjusted = res.raw;
	    if (res.raw!=0)
	    	res.adjusted /= (s1+s2 - 2*exp);
	    		
		return res;
	}
	
	public static D getDistanceFromCoMemberships(AbstractMatrix U, AbstractMatrix V ){
		return getDistanceFromCoMemberships(U, V, true);
	}
	//Optimized
	public static D getDistanceFromCoMembershipsOptimized(AbstractMatrix U, AbstractMatrix V, boolean same_nodes ){
		double[] ress = U.pByTSumPowSumDiffPowSum(V, same_nodes);
		//max, sum, powsu, max, sum, powsum, powsumDiff
//		System.err.println(U.rows() +"x"+U.columns() + " , "+V.rows() +"x"+V.columns());
		double pairCount =  U.rows()*U.rows() - (same_nodes?0 :U.rows());
        double ssU = ress[1];
        double ssV = ress[4];
        double m = Math.max(ress[0],ress[3]);

        double sU = ress[2];
        double sV = ress[5];
        double sUV = ress[6];

//        System.err.println(ssU+" "+ssV+" "+m+" "+ sU+" "+ sV+" "+ sUV);
        D res = new D();
	    res.raw =res.normalized =res.adjusted = sUV;
	    if (res.raw!=0){
		    res.normalized /=  pairCount*m*m;
		    res.adjusted /=  (sU+sV) - 2*(ssU*ssV)/(pairCount);
	    }
//	    System.err.println(ress);
//	    System.err.println(res);
	    return res;
	}
	public static D getDistanceFromCoMemberships(AbstractMatrix U, AbstractMatrix V, boolean same_nodes ){
		AbstractMatrix CU = U.productByItsTranspose();
		AbstractMatrix CV = V.productByItsTranspose();

	    if (!same_nodes){
	        CU.zeroDiagonalsInPlace();
	        CV.zeroDiagonalsInPlace();
	    }
	    
        double pairCount =  CU.rows()*CU.rows() - (same_nodes?0 :CU.rows());
        double ssU = CU.sum();
        double ssV = CV.sum();
        double m = Math.max(CU.max(),CV.max());

        double sU = CU.powSum(2);
        double sV = CV.powSum(2);
        double sUV = (CU.subtract(CV)).powSum(2);
//        System.err.println(ssU+" "+ssV+" "+m+" "+ sU+" "+ sV+" "+ sUV);

        D res = new D();
	    res.raw = sUV;
	    if (res.raw!=0){
		    res.normalized =  res.raw /(pairCount*m*m);
		    res.adjusted = res.raw / ((sU+sV) - 2*(ssU*ssV)/(pairCount));
	    }
		return res;
	}
	
	public static D omega(AbstractMatrix U, AbstractMatrix V, boolean same_nodes ){
		AbstractMatrix CU = U.productByItsTranspose();
		AbstractMatrix CV = V.productByItsTranspose();
	    if (!same_nodes){
	        CU.zeroDiagonalsInPlace();
	        CV.zeroDiagonalsInPlace();
	    }
        double pairCount =  CU.rows()*CU.rows() - (same_nodes?0 :CU.rows());
        
        AbstractVector tU = CU.getValueFrequencies(same_nodes);
        AbstractVector tV = CV.getValueFrequencies(same_nodes);

        D res = new D();
	    res.raw = pairCount-CU.subtract(CV).getNonZeroCount();
	    res.normalized =  res.raw /(pairCount);
		
	    double ew = tU.hadamardProduct(tV).sum();
//	    AbstractVector tmp = tU.hadamardProduct(tV);
	    ew = ew / (pairCount*pairCount);
	    res.adjusted = (res.normalized -ew);
	    if (res.adjusted !=0)
	    		res.adjusted /=  (1-ew);
	    
	    return res;
	}
	
	public static D getDistanceFromMatrixNorms(AbstractMatrix U,AbstractMatrix V, boolean same_nodes){
		return getDistanceFromMatrixNorms(U, V, 
			new Transformer<AbstractMatrix, Double>() {
				@Override
				public Double transform(AbstractMatrix M) {
					return Math.sqrt(M.powSum(2));
				}
			},same_nodes);
	}
	public static D getDistanceFromMatrixNorms(AbstractMatrix U,AbstractMatrix V,Transformer< AbstractMatrix, Double > norm, boolean same_nodes){
		AbstractMatrix CU = U.productByItsTranspose();
		AbstractMatrix CV = V.productByItsTranspose();
		   if (!same_nodes){
		        CU.zeroDiagonalsInPlace();
		        CV.zeroDiagonalsInPlace();
		    }
		    
		double NF = norm.transform(CU)+norm.transform(CV);
//		CU.subtract(CV);
	    D res = new D();
	    res.raw = norm.transform(CU.subtract(CV));
	    res.normalized = res.raw;
	    if (res.normalized!=0)
	    	res.normalized/=NF;
	    return res;
	}
	
	
	public static D getDistanceTarceFormsOptimized(AbstractMatrix U,AbstractMatrix V, boolean same_nodes){
		double[] ress = U.pByTSumPowSumDiffPowSum(V, same_nodes);
		//max, sum, powsu, max, sum, powsum, powsumDiff, tr, tr
//		System.err.println(U.rows() +"x"+U.columns() + " , "+V.rows() +"x"+V.columns());

		double CUTrace = ress[7]; // = CU.getTrace()
		double CVTrace = ress[8];
		double CU2Trace = ress[2];
		double CV2Trace = ress[5];
		double QTrace = ress[9];
		
		D res = new D();
		res.raw =res.normalized =res.adjusted = QTrace;
		if(res.raw!=0){
			if (same_nodes)//Not defined otherwise
				res.normalized /= (CUTrace * CVTrace);
			else res.normalized=0;
			res.adjusted /=  Math.sqrt( CU2Trace * CV2Trace );
		}
	    res.reverse();
		return res; 
		
	}
	
	
	public static D getDistanceTarceForms(AbstractMatrix U,AbstractMatrix V, boolean same_node){
//		AMatrix CU = U.productByItsTranspose();
//		AMatrix CV = V.productByItsTranspose();
//		AMatrix CU2 = CU.productByItsTranspose();
//		AMatrix CV2 = CV.productByItsTranspose();
//		    
//		AMatrix Q = CU.product(CV);
//		D res = new D();
//		res.raw = Q.getTrace();
//		res.normalized = res.raw / (CU.getTrace() * CV.getTrace());
//		res.adjusted = res.raw / Math.sqrt( CU2.getTrace() * CV2.getTrace() );
//	    res.reverse();
//		return res; 
		// tr(U.T x V) = sum (U o V) ==> so
		
		
		AbstractMatrix CU = U.productByItsTranspose();
		AbstractMatrix CV = V.productByItsTranspose();
		
		double CUTrace = U.elementwiseProduct(U).sum(); // = CU.getTrace()
		double CVTrace = V.elementwiseProduct(V).sum();
		double CU2Trace = CU.powSum(2);
		double CV2Trace = CV.powSum(2);;
		
		double QTrace = CU.elementwiseProduct(CV).sum();
		
		D res = new D();
		res.raw =res.normalized =res.adjusted = QTrace;
		if(res.raw!=0){
			res.normalized /= (CUTrace * CVTrace);
			res.adjusted /=  Math.sqrt( CU2Trace * CV2Trace );
		}
	    res.reverse();
		return res; 
		
	}
	

	public static D nmi(AbstractMatrix U,AbstractMatrix V){
		Transformer<Double, Double> vectorizef = new Transformer<Double, Double>() {
			@Override
			public Double transform(Double x) {
				return x<=0? 0:x*Math.log(x);
			}
		}; 
	//    # N = np.dot(U.T,V) # the overlaps, a.k.a contingency table
		AbstractMatrix N = U.transposeProduct(V);
//		mikera.vectorz.Vector tb2 =  mikera.vectorz.Vector.createLength(N.columnCount());
//		tb2.fill(1);
//		ColumnMatrix b2 = new  ColumnMatrix(tb2);
//
//		mikera.vectorz.Vector tb1 =  mikera.vectorz.Vector.createLength(N.rowCount());
//		tb1.fill(1);
//		RowMatrix b1 = new  RowMatrix(tb1);

		double n = N.sum();
		N.divideInPlace(n);
		AbstractVector nv = N.getElemetMarginalSum(1);
		AbstractVector nu = N.getElemetMarginalSum(0);
		N.applyFunctionInPlace(vectorizef);
		nu.inPlaceApplyFunction(vectorizef);
		nv.inPlaceApplyFunction(vectorizef);
		
	    double Hu = -nu.sum();
	    double Hv = -nv.sum();
	    double Huv = -N.sum();
	    double Iuv = Hu+Hv-Huv;
	    
	    double vi = (2*Huv-(Hu+Hv) )/Math.log(n);
	    double nmi_sum = Iuv==0?0: 2*Iuv/(Hu+Hv); // VM in sklearn
	    double nmi_sqrt = Iuv==0?0:Iuv/math.sqrt(Hu*Hv); // NMI in sklearn
	    
	    return new D(vi,nmi_sum, nmi_sqrt);
	}
//	public static D Q_modularity_with_incidence(AMatrix U,AMatrix N){
//		return Q_modularity_incidence(U, N, false);
//	}
//	public static D Q_modularity_incidence(AMatrix U,AMatrix N, boolean same_nodes ){
//		AMatrix CU = U.productByItsTranspose();
//		AMatrix A = N.productByItsTranspose();
//	    if (!same_nodes){
//	    	CU.inPlaceZero_diagonal();
//	    	A.inPlaceZero_diagonal();
//	    }
//	    double ssA = A.getElementsSum();
//	    
//		AVector D = A.getElemetMarginalSum(1);
//		AMatrix E = D.outerProduct(D);
//		E.inPlaceDivide(ssA);
//
//		A.inPlaceSub(E);
//
//		AMatrix Q = CU.product(A);
//        D res = new D();
//	    res.raw = Q.getElementsSum();
//	    res.normalized =  res.raw /ssA;
//		return res;
//	}
//	
	
	public void test (){
		Vector<Set<Integer>> groundTruth = DatasetUtils.createClusteringFromArray(new int[][]{ {0,1,2,3},{4,5,6}});
		Vector<Set<Integer>> p1 =  DatasetUtils.createClusteringFromArray(new int[][]{ {0,1,2},{3,4,5,6}});
		Vector<Set<Integer>> p2 =  DatasetUtils.createClusteringFromArray(new int[][]{ {0},{1,2,3,4,5,6}});

		for (ClusteringAgreement<Integer> meas : MeasuresUtil.<Integer>getAgreementAlternatives()) {
//			System.err.println(meas + "(t,t) =" +meas.getAgreement(t, t));
//			System.err.println(meas + "(V,V) =" +meas.getAgreement(groundTruth, groundTruth));
			System.err.println(meas + "(U1,V) =" +meas.getAgreement(p1, groundTruth));
			System.err.println(meas + "(U2,V) =" + (meas.getAgreement(p2, groundTruth)));
			System.err.println();
		}
		
		Transformer<Integer, Integer> ids = new Transformer<Integer, Integer>() {
			@Override
			public Integer transform(Integer input) {
				return (Integer)input;
			}
		};
		AbstractMatrix U = AlgebricClusteringAgreement.getClusteringMatrix( 7, groundTruth, ids);
		AbstractMatrix V = AlgebricClusteringAgreement.getClusteringMatrix( 7,p1,ids);
		AbstractMatrix V2 = AlgebricClusteringAgreement.getClusteringMatrix(7, p2,ids);

		System.out.println(U);
//		System.out.println(U.columnCount());
		System.out.println(V);
//		System.out.println(V.columnCount());
		D res = getDistanceFromOverlap(U, V,true);
		res.reverse();
		System.err.println(res);
		res = getDistanceFromOverlap(U, V2,true);
		res.reverse();
		System.err.println(res);
		System.err.println(nmi(U, V));
		System.err.println(nmi(U, V2));
		res = getDistanceFromCoMemberships(U, V,false);
		res.reverse();
		System.err.println(res);
		res = getDistanceFromCoMemberships(U, V2,false);
		res.reverse();
		System.err.println(res);
		res = omega(U, V,false);
//		res.reverse();
		System.err.println(res);
//		res = omega(U, V2,false);
//		res.reverse();
//		System.err.println(res);
	}
	public static void main(String[] params){
		AlgebricClusteringDistance algebricClusteringDistance = new AlgebricClusteringDistance();
		
		algebricClusteringDistance.test();
	}
}
