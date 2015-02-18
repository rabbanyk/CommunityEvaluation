package data;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;

import org.apache.commons.collections15.Transformer;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;


public class SparseMatrixMTJ<T>{// implements AMatrix<AMatrix>{
	Matrix matrix ;
	
	public SparseMatrixMTJ(Matrix M) {
		this.matrix = M;
	}

//	@Override
//	public AMatrix product(AMatrix b) {
//		if(b instanceof SparseMatrixMTJ){
//			@SuppressWarnings("unchecked")
//			SparseMatrixMTJ<T> B = (SparseMatrixMTJ<T>)b;
//			return new SparseMatrixMTJ<Double>(matrix.mult(matrix, B.matrix));
//		}
//		return null;
//	}

//	@Override
//	public AMatrix transposeProduct(AMatrix b) {
//		if(b instanceof SparseMatrixMTJ){
//			@SuppressWarnings("unchecked")
//			SparseMatrixMTJ<T> B = (SparseMatrixMTJ<T>)b;
//			return new SparseMatrixMTJ<Double>(matrix.transAmult(matrix, B.matrix));
//		}
//		return null;
//	}

	
}
