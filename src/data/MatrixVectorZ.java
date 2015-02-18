package data;

import java.util.Iterator;
import java.util.List;

import org.apache.commons.collections15.Transformer;
import org.la4j.Matrix;
import org.la4j.matrix.ColumnMajorSparseMatrix;

import mikera.arrayz.INDArray;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.matrixx.impl.ColumnMatrix;
import mikera.matrixx.impl.IdentityMatrix;
import mikera.matrixx.impl.SparseColumnMatrix;
import mikera.matrixx.impl.SparseRowMatrix;
import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.Vectorz;
import mikera.vectorz.ops.AFunctionOp;


public class MatrixVectorZ extends SparseRowMatrix implements AbstractMatrix{


	public MatrixVectorZ(int vertexCount, int edgeCount) {
		super(vertexCount, edgeCount);
	}
	
	public MatrixVectorZ(AVector[] data, int rowCount, int columnCount) {
		super(data, rowCount, columnCount);
	}

	public MatrixVectorZ(List<AVector> data, int rowCount, int columnCount) {
		super(data, rowCount, columnCount);
	}
	public static MatrixVectorZ create(int rows, int cols) {
		return new MatrixVectorZ(rows, cols);
	}
	public MatrixVectorZ wrap(SparseRowMatrix a){
		return new MatrixVectorZ(a.getRows(),a.rowCount(),a.columnCount());
	}
	
	public static MatrixVectorZ create(AMatrix source) {
		if (source instanceof SparseColumnMatrix) 
			source = ((SparseColumnMatrix)source).toSparseRowMatrix();
		int rc = source.rowCount();
		int cc = source.columnCount();
		AVector[] data = new AVector[rc];
		List<AVector> rows=source.getRows();
		for (int i = 0; i < rc; i++) {
			AVector row = rows.get(i);
			if (!row.isZero()) {
			    data[i] = Vectorz.createSparse(row);
			}
		}
		return new MatrixVectorZ(data,rc,cc);
	}
	@Override
	public MatrixVectorZ exactClone() {
		MatrixVectorZ result = new MatrixVectorZ(rows, cols);
        for (int i = 0; i < rows; ++i) {
			AVector row = unsafeGetVector(i);
			if (row != null)
                result.replaceRow(i, row.exactClone());
		}
		return result;
	}
	
	public MatrixVectorZ innerProduct(SparseColumnMatrix a) {
		MatrixVectorZ r = create(rows, a.columnCount());
		for (int i = 0; i < rows; ++i) {
			AVector row = unsafeGetVector(i);
            if (row != null) {
			    r.replaceRow(i,row.innerProduct(a));
            }
		}
		return r;
	}

	
	public MatrixVectorZ getCloneTransposeView() {
		MatrixVectorZ clone =  this.exactClone();
		return MatrixVectorZ.create(clone.getTransposeView().toSparseRowMatrix());
	}
	@Override
	public <M extends AbstractMatrix> M elementwiseProduct(M b) {
		return (M) this.wrap(super.wrap((this.multiplyCopy((INDArray)b).toVector())));
	}

	@Override
	public <M extends AbstractMatrix> M productByItsTranspose() {
		
		MatrixVectorZ clone =  this.exactClone();
		SparseColumnMatrix Tr = clone.getTransposeView();
				
		return (M) this.innerProduct(Tr);
//		return (M) this.innerProduct(this.getCloneTransposeView());
	}

	@Override
	public <M extends AbstractMatrix> M transposeProduct(M b) {
		return (M) wrap((SparseRowMatrix) this.getCloneTransposeView().innerProduct((AMatrix)b));

	}

	@Override
	public int rows() {
		return rowCount();
	}

	@Override
	public int columns() {
		return columnCount();
	}

	@Override
	public int getNonZeroCount() {
		return (int)nonZeroCount();
	}

	@Override
	public <V extends AbstractVector> V getValueFrequencies(boolean countDiagonals) {
		VectorVectorZ tU = new VectorVectorZ(Vectorz.newVector((int)elementMax() + 1));

		for (Iterator iterator = elementIterator(); iterator.hasNext();) {
			Double val = (Double) iterator.next();
			tU.addAt(val.intValue(), 1);
		}
	    tU.subAt(0, countDiagonals?0:rowCount());

		return (V) tU;
	}

	@Override
	public <V extends AbstractVector> V getElemetMarginalSum(int dimension) {
		V res =null ;
		if(dimension ==1){
			res = (V) VectorVectorZ.wrap(new double[this.rows()]);
			for (int i = 0; i < res.length(); i++) {
				res.set(i,getRow(i).elementSum());
			}
		}else if (dimension ==0){
			AMatrix tmp = this.toSparseColumnMatrix();
			res = (V) VectorVectorZ.wrap(new double[tmp.columnCount()]);
			for (int i = 0; i < res.length(); i++) {
				res.set(i,tmp.getColumn(i).elementSum());
			}
		}
		return  res;
		
//		
//		mikera.vectorz.Vector tb =  mikera.vectorz.Vector.createLength(dimension==0?columnCount():rowCount());
//		tb.fill(1);
//		ColumnMatrix b = new  ColumnMatrix(tb);
//		AMatrix res = dimension==0? innerProduct(b):b.innerProduct(this);
//		return (V) VectorVectorZ.wrap(res.toDoubleArray());
	}

	@Override
	public double sum() {
		return elementSum();
	}

	@Override
	public double powSum(int p) {
		return elementPowSum(p);
	}

	@Override
	public double max() {
		return elementMax();
	}

	@Override
	public <M extends AbstractMatrix> M subtract(M b) {
		MatrixVectorZ res = this.exactClone();
		res.sub((AMatrix) b);
		return (M) res;		
	}

	@Override
	public void divideInPlace(double m) {
		divide(m);
	}

	@Override
	public void applyFunctionInPlace(final Transformer<Double, Double> f) {
		Op vectorizef = new AFunctionOp() {
			@Override
			public double apply(double x) {
				return f.transform(x);
			}
		};	
		this.applyOp(vectorizef);		
	}

	@Override
	public void zeroDiagonalsInPlace() {
		IdentityMatrix I = IdentityMatrix.create(rowCount());
    	sub(multiplyCopy(I));		
	}

	@Override
	public <M extends AbstractMatrix> double[] pByTSumPowSumDiffPowSum(M b, boolean diag) {
		// TODO Auto-generated method stub
		return null;
	}

}
