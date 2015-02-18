package measure.base;

import java.util.Set;
import java.util.Vector;

public interface RelativeCriteria<V> {
	public double evaluate(Vector<Set<V>> clusters);
	public String getName();
}
