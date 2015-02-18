package measure.correlation;

import java.util.Vector;


/**
 * @author rabbanyk
 * @see {@link <a href="http://en.wikipedia.org/wiki/Correlation_and_dependence"> Correlation_and_dependence from wiki</a> }
 */
public class PearsonCorrelation extends Correlation<Vector<? extends Number>>{

	@Override
	public Number compute(Vector<? extends Number> X, Vector<? extends Number> Y) {
		if (X.size() != Y.size()) throw new IllegalArgumentException("Size of given vectors should be the same");
		
		
		double nom = 0, domX =0, domY =0 ;
		double avgX = getAvg(X) , avgY = getAvg(Y);
		
		for (int i = 0; i < X.size(); i++) {
			nom += (X.get(i).doubleValue()-avgX)*(Y.get(i).doubleValue()-avgY);
			
			domX += (X.get(i).doubleValue()-avgX)*(X.get(i).doubleValue()-avgX);
			domY += (Y.get(i).doubleValue()-avgY)*(Y.get(i).doubleValue()-avgY);
		}
		
	//	System.err.println(nom +" "+  domX * domY);
		return (nom==0 &&  domX * domY ==0)? 0 : nom / Math.sqrt(domX * domY);
	}

	
	private double getAvg(Vector<? extends Number> X){
		double res = 0;
		
		for (Number number : X) {
			res += number.doubleValue();
		}
		
		return res/X.size();
	}
	
	public String toString(){
		return "Pearson Correlation";
	}
}
