package data;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections15.Transformer;
import org.la4j.iterator.MatrixIterator;
import org.la4j.Matrix;
import org.la4j.matrix.ColumnMajorSparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;
import org.la4j.matrix.sparse.CRSMatrix;

public class Matrixla4j extends CRSMatrix implements AbstractMatrix{

	public Matrixla4j(int vertexCount, int edgeCount) {
		super(vertexCount, edgeCount);
	}

	@Override
    public Matrix blankOfShape(int rows, int columns) {
        return new Matrixla4j(rows, columns);
    }
	static double log(double x){
		return x==0?0:Math.log(x+1); 
	}

	public <M extends AbstractMatrix> double[] pByTSumPowSum() {
		double[] res = {0,0};
	        List<Integer> nzRows = new ArrayList<Integer>();
	        Iterator<Integer> it = this.iteratorOfNonZeroRows();

	        while (it.hasNext()) {
	            nzRows.add(it.next());
	        }

	        double tmp;
	        for (int i: nzRows) {
	            for (int j: nzRows) {
	                tmp = nonZeroIteratorOfRow(i)
                            .innerProduct(nonZeroIteratorOfRow(j));
	                res[0] += tmp;
	                res[1] += tmp*tmp;
//	            	result.set(i, j, a.nonZeroIteratorOfRow(i)
//	                                  .innerProduct(a.nonZeroIteratorOfRow(j)));
	            }
	        }
		return res;
	}
	
	@Override
	public <M extends AbstractMatrix> double[] pByTSumPowSumDiffPowSum(M b, boolean diag) {
		double[] res = {0,0,0,0,0,0,0,
						0,0,0 };//,0,0,0,0};
		//0:max1, 1:sum1, 2:powsu1, 3:max2, 4:sum2, 5:powsum2, 6:powsumDiff12, 
		//7:Tr1, 8:Tr2, 9:timessum12,
		//10:logsu1, 11: logsu2, 12: logsumdiff12, 13: timeslog12

		Matrixla4j B = (Matrixla4j)b;
		
		long startTime = System.currentTimeMillis();

        Set<Integer> nzRows = new HashSet<Integer>();
        Iterator<Integer> it = this.iteratorOfNonZeroRows();
        while (it.hasNext()) {
            nzRows.add(it.next());
        }
        
        double tmp;
        for (int i: nzRows) {
            for (int j: nzRows) if(i!=j || diag)
            {
                tmp = nonZeroIteratorOfRow(i)
                        .innerProduct(nonZeroIteratorOfRow(j));
                if (tmp>res[0]) res[0]=tmp;
                res[1] += tmp;
                res[2] += tmp*tmp;
//                res[10] += tmp*log(tmp);
                if (i==j) res[7]+=tmp;
            }
        }
        System.err.println("\t\t      >>  -- su2(U): " + (System.currentTimeMillis() - startTime) +" milisecond");
		startTime = System.currentTimeMillis();
        Set<Integer> nzRowsB = new HashSet<Integer>();
        Iterator<Integer> itB = (B).iteratorOfNonZeroRows();
        while (itB.hasNext()) {
        	nzRowsB.add(itB.next());
        }
        for (int i: nzRowsB) {
            for (int j: nzRowsB) if(i!=j || diag)
            	{
                tmp = B.nonZeroIteratorOfRow(i)
                        .innerProduct(B.nonZeroIteratorOfRow(j));
                if (tmp>res[3]) res[3]=tmp;
                res[4] += tmp;
                res[5] += tmp*tmp;
//                res[11] += tmp*log(tmp);
                if (i==j) res[8]+=tmp;
            }
        }
        System.err.println("\t\t      >>  -- su2(V): " + (System.currentTimeMillis() - startTime) +" milisecond");
		startTime = System.currentTimeMillis();
        
		nzRows.retainAll(nzRowsB);
        double tmp2= 0;
        for (int i: nzRows) {
            for (int j: nzRows) if(i!=j || diag)
            	{
                tmp = nonZeroIteratorOfRow(i)
                        .innerProduct(nonZeroIteratorOfRow(j));
                tmp2 = B.nonZeroIteratorOfRow(i)
                        .innerProduct(B.nonZeroIteratorOfRow(j));
                
                res[9] += tmp*tmp2;
//                res[13] += tmp*log(tmp2)+tmp2*log(tmp);
        	}
        }
        System.err.println("\t\t      >>  -- su2(U.V): " + (System.currentTimeMillis() - startTime) +" milisecond");

        res[6] = res[2]+res[5] - 2*res[9] ;
//        res[12] = res[10]+res[11] - res[13] ;
        
		return res;
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public <M extends AbstractMatrix> M productByItsTranspose() {
		M res;

//		System.err.println(">>>>>>M    "+ rows()+"x"+columns()+" : " + getNonZeroCount());
//		long tmp = System.nanoTime();
//		ColumnMajorSparseMatrix Tr = (ColumnMajorSparseMatrix) this.transpose();
//		res =  (M) this.multiply(Tr);
//		long t1 = (System.nanoTime() - tmp);
//		tmp= System.nanoTime();
		res =  (M) this.multiplyByItsTranspose();
//		long t2 = (System.nanoTime() - tmp);
//		System.err.println("----- MT toook: " +  "  " + t2  );
		
		return res;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <M extends AbstractMatrix> M elementwiseProduct(M b) {
		return (M) this.hadamardProduct((Matrix)b);
	}

	@SuppressWarnings("unchecked")
	@Override
	public <M extends AbstractMatrix> M transposeProduct(M b) {
		ColumnMajorSparseMatrix columnMajorSparseMatrix= (ColumnMajorSparseMatrix) this.transpose();
		columnMajorSparseMatrix.multiply((Matrix)b);
		return (M) this.transpose().multiply((Matrix)b);
	}
	
	@Override
	public int getNonZeroCount() {
		return this.cardinality();
	}

	
	@Override
	public void divideInPlace(double m) {
		MatrixIterator it = nonZeroIterator();
		while (it.hasNext()) {
		  double x = it.next();
		  it.set(x/m);
		}	
	}

	@Override
	public void applyFunctionInPlace(final Transformer<Double, Double> f) {
		MatrixIterator it = nonZeroIterator();
		while (it.hasNext()) {
		  double x = it.next();
		  it.set(f.transform(x));
		}
//		MatrixProcedure procedure = new MatrixProcedure() {
//			@Override
//			public void apply(int i, int j, double value) {
//				set(i,j,f.transform(value));
//			}
//		};
//		this.eachNonZero(procedure );
	}

	@Override
	public void zeroDiagonalsInPlace() {
//      Is better if number of non zero elements < n log n	
//		MatrixIterator it = nonZeroIterator();
//		while (it.hasNext()) {
//			double x= it.next();
//			int i = it.rowIndex(), j = it.columnIndex();
//			if (i==j) {
//				it.set(0);
//			}
//		}
	
		for(int i =0; i<this.rows();i++)
			this.set(i, i, 0);

	}
	
	@Override
	public double powSum(int p) {
		double res =0;
		MatrixIterator it = nonZeroIterator();
		while (it.hasNext()) {
		  double x = it.next();
		  res+=Math.pow(x, p);
		}
		return res;
	}


	@SuppressWarnings("unchecked")
	@Override
	public <V extends AbstractVector> V getValueFrequencies(boolean countDiagonals) {
		double[] freq = new double[(int)max()+1];
		MatrixIterator it = nonZeroIterator();
		while (it.hasNext()) {
		  double x = it.next();
		  freq[(int)x]++;
		}

		Vectorla4j res = new Vectorla4j(freq);
		double pairCount =  rows()*columns() - (countDiagonals?0 :rows());
		res.set(0,  pairCount-res.sum());
		return (V) res;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <V extends AbstractVector> V getElemetMarginalSum(int dimension) {
		V res =null ;
		if(dimension ==1){
			res = (V) new Vectorla4j(this.rows());
			for (int i = 0; i < res.length(); i++) {
				res.set(i,getRow(i).sum());
			}
		}else if (dimension ==0){
			Matrix tmp = this.toColumnMajorSparseMatrix();
			res = (V) new Vectorla4j(tmp.columns());
			for (int i = 0; i < res.length(); i++) {
				res.set(i,tmp.getColumn(i).sum());
			}
		}
		return  res;
	}

	@Override
	public <M extends AbstractMatrix> M subtract(M b) {
//		System.err.println("  A = \n"+ this);
//		System.err.println("  B = \n"+ b);
		Matrix result = this.subtract((Matrix)b);
//		System.err.println("  A-B = \n"+ result);
		return (M)result;
//		return (M) this.add(((Matrix)b).multiply(-1));
//		System.err.println("  A = \n"+ this);
//		System.err.println("  B = \n"+ b);
//		
//		if (b instanceof Matrixla4j){
//			Matrixla4j B = (Matrixla4j)b;
//			Matrix result = this.blank();
//	        MatrixIterator these = this.nonZeroRowMajorIterator();
//	        MatrixIterator those = B.nonZeroRowMajorIterator();
//	        MatrixIterator both = these.orElseSubtract(those);
//	
//	        while (both.hasNext()) {
//	            double x = both.next();
//	            int i = both.rowIndex();
//	            int j = both.columnIndex();
//	            result.set(i, j, x);
//	        }
//	        System.err.println("--------> "+this.get(1, 4)+" - "+B.get(1, 4)+" = "+result.get(1, 4));
//	      
//			System.err.println("  A-B = \n"+ result);
//
//	        return (M) result;
//		}
//        return null;//(M) this.subtract((Matrix)b);
	}

	
}


//
//	@Override
//	public AVector getValueFrequencies() {
//		double[] freq = new double[(int)max()+1];
//		MatrixIterator it = nonZeroIterator();
//		while (it.hasNext()) {
//		  double x = it.next();
//		  freq[(int)x]++;
//		}
//		return new Vectorla4j(freq);
////		class FreqProce implements MatrixProcedure{
////			public int[] freq = new int[(int)max()+1];
////			@Override
////			public void apply(int i, int j, double value) {
////				freq[(int)value]++;
////			}
////		};
////		FreqProce procedure = new FreqProce() ;
////		this.eachNonZero(procedure );		
////		return procedure.freq;
//	}
//@Override
//public AVector getElemetMarginalSum(int dimension) {
//	double[] freq = new double[(int)max()+1];
//	MatrixIterator it = nonZeroIterator();
//	while (it.hasNext()) {
//	  double x = it.next();
//	  freq[(int)x]++;
//	}
//	return new Vectorla4j(freq);
//}

//@Override
//public Matrixla4j elementwiseProduct(Matrixla4j b) {
//	return (Matrixla4j) this.hadamardProduct(b);
//}

//@Override
//public Matrixla4j productByItsTranspose() {
//	return (Matrixla4j) 
//}

//@Override
//public Matrixla4j getTranspose() {
//	return (Matrixla4j) this.transpose();
//}
//@Override
//public Matrixla4j transposeProduct(Matrixla4j b) {
//	return (Matrixla4j) this.getTranspose().multiply(b);
//}

//@Override
//public int rows() {
//	return this.rows();
//}

//@Override
//public int columns() {
//	return this.columns();
//}

//
//@Override
//public void set(int i, int j, double value) {
//	this.set(i, j, value);
//}


//
//@Override
//public double trace() {
//	return this.trace();
//}
//
//@Override
//public double sum() {
//	return this.sum();
//}
//
//@Override
//public double max() {
//	return this.max();
//}

