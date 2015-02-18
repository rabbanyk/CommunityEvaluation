package data;

import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.ops.AFunctionOp;
import mikera.vectorz.Vector;

import org.apache.commons.collections15.Transformer;

public class VectorVectorZ extends Vector implements AbstractVector{

//	public VectorVectorZ(int length) {
//        this(new double[length]);
//	}
	
	public VectorVectorZ(AVector source) {
		super(source);
	}
	/**
	 * Wraps a double array into a Vector, does *no defensive copy* so use with caution
	 * @param source
	 * @return
	 */
	public static VectorVectorZ wrap(AVector source) {
		return new VectorVectorZ(source);
	}
	
	public static VectorVectorZ wrap(double[] source) {
		return new VectorVectorZ(Vector.wrap(source));
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public <V extends AbstractVector> V hadamardProduct(V b) {
		return (V) VectorVectorZ.wrap(this.multiplyCopy((AVector) b));
	}


	@Override
	public <V extends AbstractVector, M extends AbstractMatrix> M outerProduct(V b) {
		@SuppressWarnings("unchecked")
		M result = (M) new MatrixVectorZ(this.length(), b.length());

        for (int i = 0; i < this.length(); i++) {
            for (int j = 0; j < b.length(); j++) {
                result.set(i, j, this.get(i) * b.get(j));
            }
        }
        return result;
//		return (M) this.outerProduct((Vector)b);
	}

	@Override
	public void inPlaceApplyFunction(final Transformer<Double, Double> f) {
		Op vectorizef = new AFunctionOp() {
			@Override
			public double apply(double x) {
				return f.transform(x);
			}
		};	
		this.applyOp(vectorizef);
	}
	@Override
	public double sum() {
		return elementSum();
	}


}
