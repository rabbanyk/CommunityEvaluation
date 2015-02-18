package measure.graph.distance;

import java.util.Set;
import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.graph.Graph;
import measure.graph.GraphProximity;

public class NumberOfPathes<V, E> extends GraphProximity<V,E> {
	public enum Coeff{UNI, EXP, LIN};
	Coeff coeff = Coeff.UNI;
	double MaxLenght = 3;
	
	{
		type = Type.SIMILARITY;
	}
	public NumberOfPathes(Graph<V, E> g) {
		this(g, Coeff.UNI);
	}
	public NumberOfPathes(Graph<V, E> g, Coeff coeff) {
		this(g, coeff, 3);
	}
	public NumberOfPathes(Graph<V, E> g, Coeff coeff,double MaxLenght) {
		this(g, null, coeff, MaxLenght);
	}
	public NumberOfPathes(Graph<V, E> g, Transformer<E, ? extends Number> weights) {
		this(g, weights, Coeff.UNI);
	}
	public NumberOfPathes(Graph<V, E> g, Transformer<E, ? extends Number> weights, Coeff coeff) {
		this(g, weights,coeff, 3);
	}
	public NumberOfPathes(Graph<V, E> g, Transformer<E, ? extends Number> weights, Coeff coeff,double MaxLenght) {
		this(g, weights,coeff, MaxLenght, false);
	}
	public NumberOfPathes(Graph<V, E> graph,Transformer<E, ? extends Number> weights, Coeff coeff , boolean selfWeight) {
		this(graph, weights,coeff, 3, selfWeight);
	}

	public NumberOfPathes(Graph<V, E> graph,Transformer<E, ? extends Number> weights, Coeff coeff,double MaxLenght , boolean selfWeight) {
		this(graph, weights,coeff, MaxLenght, selfWeight, false);
	}

	public NumberOfPathes(Graph<V, E> graph,
			Transformer<E, ? extends Number> weights, Coeff coeff,double MaxLenght, boolean selfWeight,	boolean normalized) {
		super(graph, weights, selfWeight, normalized);
		this.coeff = coeff;
		this.MaxLenght=MaxLenght;
	}

	public Number computeMeasure(V source, V target) {
		double sp1=0,sp2=0,sp3=0;
		//Paths length 1
		sp1 =  W(source, target);

		//paths length 2
		if(MaxLenght>=2){
			Set<V> cNeigh = getUnionNeighbours(source, target);
			
			double sk, tk; 
			for (V k : cNeigh) {
				sk = W(source, k); tk = W(target, k);
				sp2 += sk*tk;
			}
		
			//paths length 3
			if(MaxLenght>=3){
				double sp,tq, pq; 
				for (V p : cNeigh) for (V q : cNeigh) {//if(graph.getNeighbors(p).contains(q)){
					sp = W(source, p);
					tq = W(target, q);
					pq = W(p, q);
					sp3 += sp*tq*pq;
				}
			}
		}
		
		double sim=0;	
		if( coeff == Coeff.UNI)
			sim = sp1+sp2+sp3;
		else if( coeff == Coeff.LIN)
			sim = sp1+sp2/2+sp3/3;
		else if( coeff == Coeff.EXP)
			sim = sp1+Math.pow(sp2, 1./2.) +Math.pow(sp3, 1./3.);
			
		return sim;
	}
	public Number getDistance(V source, V target){
		return MaxLenght*Math.max(Math.pow(MaxWeight,MaxLenght),MaxWeight)*Math.pow(graph.getVertexCount(),MaxLenght) - getProximity(source, target).doubleValue();///3*Math.max(Math.pow(MaxWeight,3),MaxWeight)*graph.getVertexCount() ;
	}
	public String toString(){
		return "NumberOfPathes"+type+coeff+MaxLenght+ ((selfWeight)?"Aug":"");
	}
	@Override
	public String getName() {
		if (selfWeight)
			return "$\\hat{NP"+(coeff==Coeff.LIN?"L":(coeff==Coeff.EXP?"E":""))+MaxLenght+"}$";
		else return "NP"+(coeff==Coeff.LIN?"L":(coeff==Coeff.EXP?"E":""))+MaxLenght;
	}

}
