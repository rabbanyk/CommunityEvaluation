package measure.criteria;

public class Dunn<V> extends RelativeClusterCriteria<V> {
	public Dunn() {
		this(Between.SINGLE,Within.MAX,false);
	}
	
	public Dunn(Between setDIST, Within diameter,boolean closnessCompatible) {
		super();
		this.betweenSetMethod = setDIST;
		this.withinSetMethod = diameter;
		this.closnessCompatible = closnessCompatible;
	}
	
	public boolean isMaximizer() {
		return (order==order.MAXIMIZER);
	}
	// Want to maximize
	protected double evaluate(){
		Double mindelta = null; double d ;
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) 
				if (i!=j){
					d = getBetween(i, j);
					if(mindelta==null || (isSimInUse()? (d>mindelta): (d < mindelta))) mindelta = d;
				}	
		}
		
		Double maxDelta = null;
		for (int i = 0; i < k; i++) {
			d = getWithin(i);
			if(maxDelta==null ||(isSimInUse()? (d<maxDelta):( d > maxDelta))) maxDelta = d;
		}

		if(isSimInUse()) {
			double tmp= mindelta;
			mindelta= maxDelta;
			maxDelta = mindelta;
		}

//		System.err.println(isSimInUse()+" "+mindelta + " "  + maxDelta);
		
		return (maxDelta!=0? ( mindelta / maxDelta) : 0) ;
	}
	 
	public String toString(){
		return "Dunn " + betweenSetMethod + " " +withinSetMethod+" "  +closnessCompatible+" " + super.toString();
	}
	public String getName() {
		return "Dunn" + betweenSetMethod.ordinal() + withinSetMethod.ordinal()+(closnessCompatible?"'":"")+ " "+ proximityMeasure.getName();
	}

}
