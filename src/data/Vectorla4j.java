package data;

import org.apache.commons.collections15.Transformer;
import org.la4j.iterator.MatrixIterator;
import org.la4j.iterator.VectorIterator;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.vector.dense.BasicVector;

public class Vectorla4j extends BasicVector implements AbstractVector{

	public Vectorla4j(double[] freq) {
		super(freq);
	}
	public Vectorla4j(int length) {
        this(new double[length]);
	}
	
	@Override
    public Vector blankOfLength(int length) {
        return new Vectorla4j(length);
    }
	
	@SuppressWarnings("unchecked")
	@Override
	public <V extends AbstractVector> V hadamardProduct(V b) {
		if (this.length<b.length()){
			return (V) this.hadamardProduct(((Vector) b).sliceLeft(this.length()));
		}else if (b.length()<this.length()){
			return (V) this.sliceLeft(b.length()).hadamardProduct((Vector) b);
		}
		return (V) this.hadamardProduct((Vector) b);
	}


	@Override
	public <V extends AbstractVector, M extends AbstractMatrix> M outerProduct(V b) {
		@SuppressWarnings("unchecked")
		M result = (M) new Matrixla4j(this.length(), b.length());

        for (int i = 0; i < this.length(); i++) {
            for (int j = 0; j < b.length(); j++) {
                result.set(i, j, this.get(i) * b.get(j));
            }
        }
        return result;
//		return (M) this.outerProduct((Vector)b);
	}

	@Override
	public void inPlaceApplyFunction(Transformer<Double, Double> f) {
		VectorIterator it = this.iterator();
		while (it.hasNext()) {
		  double x = it.next();
		  it.set(f.transform(x));
		}		
	}


}
