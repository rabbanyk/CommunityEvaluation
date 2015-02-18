package data;

import org.apache.commons.collections15.Transformer;

public interface AbstractMatrix {

	// public abstract M product(M b);
	public abstract <M extends AbstractMatrix> M elementwiseProduct(M b);

	public abstract <M extends AbstractMatrix> M productByItsTranspose();
	
	// So to not store the matrix, just sum over
//	public abstract <M extends AbstractMatrix> double[] pByTSumPowSumDiffPowSum(M b, boolean diag, Transformer<Pair<Double,Double>, Double> f);

	public abstract <M extends AbstractMatrix> double[] pByTSumPowSumDiffPowSum(M b, boolean diag);

	public abstract <M extends AbstractMatrix> M transposeProduct(M b);
	
//	public abstract <M extends AMatrix> M getTranspose();
	// public abstract AMatrix productTranspose(AMatrix b);

	public abstract int rows();

	public abstract int columns();

	public abstract int getNonZeroCount();

	public abstract void set(int i, int j, double value);

	public abstract double trace();

	public abstract <V extends AbstractVector> V getValueFrequencies(boolean countDiagonal);

	public abstract <V extends AbstractVector> V getElemetMarginalSum(int dimension);

	public abstract double sum();

	public abstract double powSum(int p);

	public abstract double max();


	public abstract <M extends AbstractMatrix> M subtract(M b);

	public abstract void divideInPlace(double m);

	public abstract void applyFunctionInPlace(Transformer<Double, Double> f);

	public abstract void zeroDiagonalsInPlace();


}
