package data;

import org.apache.commons.collections15.Transformer;

public interface AbstractVector {

	/**
	 * returns the inner product with given factor, if length of the vectors are not the same, the minimum is considered instead of error
	 * 
	 * @param tV
	 * @return
	 */
	public abstract  <V extends AbstractVector> V hadamardProduct(V b);

	public abstract  <V extends AbstractVector,M extends AbstractMatrix> M outerProduct(V b);

	public abstract  double sum();

	public abstract  void inPlaceApplyFunction(Transformer<Double, Double> f);

	public int length();
	public double get(int i) ;
	public void set(int i, double value);
}

//	public abstract  <M extends AMatrix> M productByItsTranspose();