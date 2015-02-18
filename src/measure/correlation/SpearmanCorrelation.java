package measure.correlation;

import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

//http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
public class SpearmanCorrelation extends Correlation<Vector<? extends Number>>{

	@Override
	public Number compute(Vector<? extends Number> X, Vector<? extends Number> Y) {
		if (X.size() != Y.size()) throw new IllegalArgumentException("Size of given vectors should be the same");
	
		//Ranking
		Vector<Double> rankedX = rank(X, sort(X));
		Vector<Double> rankedY = rank(Y, sort(Y));
		return (new PearsonCorrelation()).compute(rankedX, rankedY);
	}
	
	private Vector<Double> rank(final Vector<? extends Number> X, Vector<Integer> sortedIndexes){
		Vector<Double> ranked = new Vector<Double>();
		for (int i = 0; i < sortedIndexes.size(); i++) ranked.add(0.);
		for (int i=0; i < sortedIndexes.size(); ) {
			int j = i+1;
			while(j < sortedIndexes.size() && X.get(sortedIndexes.get(i)).equals(X.get(sortedIndexes.get(j))))	j++; 
			double n = j-i;
			double avgR = (n-1)/2 + i+1;
			
			for (; i < j; i++) ranked.set(sortedIndexes.get(i), avgR);
		}
		return ranked;
	}
	
	private Vector<Integer> sort(final Vector<? extends Number> X){
		Vector<Integer> res = new Vector<Integer>();
		for (int i = 0; i < X.size() ; i++) {
			res.add(i);
		}
		
		Collections.sort(res , new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				return Double.compare(X.get(o1).doubleValue(),X.get(o2).doubleValue());
			}
		});
		return res;
	}
	
	public String toString(){
		return "Spearman's Correlation";
	}
	
	public static void main(String args[]){
		SpearmanCorrelation correlation = new SpearmanCorrelation();
		Vector<Double> X = new Vector<Double>();
		X.add(.3);
		X.add(-.1);
		X.add(-.2);
		X.add(.4);
		X.add(.5);
		Vector<Double> Y = new Vector<Double>();
		Y.add(-.1);
		Y.add(-.1);
		Y.add(.3);
		Y.add(.4);
		Y.add(.5);
		System.err.println(correlation.compute(X, Y));
		
	}

}
