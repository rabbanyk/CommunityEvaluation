package measure.graph.centrality;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;

import org.apache.commons.collections15.Transformer;
import edu.uci.ics.jung.algorithms.util.MapBinaryHeap;

//TODO: Justins implementation, haven't checked it
public class dev_EstimatedBetweenness<V,E>  extends GraphVertexScorer<V,E>{

    private class BetweennessData {
        double distance;
        double numSPs;
        List<E> incomingEdges;
        double dependency;

        BetweennessData() {
            distance = -1;
            numSPs = 0;
            incomingEdges = new ArrayList<E>();
            dependency = 0;
        }
       
        public String toString() {
        	return "[d:" + distance + ", sp:" + numSPs + 
        		", p:" + incomingEdges + ", d:" + dependency + "]\n";
        }
    }
    
    private class BetweennessComparator implements Comparator<V>
    {
		public int compare(V v1, V v2) 
		{
			return vertex_data.get(v1).distance > vertex_data.get(v2).distance ? 1 : -1;
		}
    }
	

	protected void clear(){
		super.clear();
	//	vertex_data = null;
	}
    //public Map<E, Double> edge_scores;
	public Map<V, BetweennessData> vertex_data;

	protected void initializeAllScores() {
		
		//to normalized the betweenness centrality 
		//int n = graph.getVertexCount();
		//double divisor = (n-1)*(n-2)/2;
		vertexScores = new HashMap<V, Double>();
		
		if (weights == null) {
			weights = new Transformer<E, Double>() {
				public Double transform(E arg0) {
					if(arg0 == null) return new Double(0);
					return new Double(1);
				}
			};
		}
		
		computeBetweenness(new MapBinaryHeap<V>(new BetweennessComparator()));
		
		/*for(V v : this.graph.getVertices()){
			double score = bc.getVertexScore(v);
			this.vertexScores.put(v, score / divisor);	
		}*/
				
	}

	public void computeBetweenness(Queue<V> queue){
		for (V v : graph.getVertices()){
			// initialize the betweenness data for this new vertex
			for (V s : graph.getVertices()) 
				this.vertex_data.put(s, new BetweennessData());


            vertex_data.get(v).numSPs = 1;
            vertex_data.get(v).distance = 0;

            Stack<V> stack = new Stack<V>();
            queue.offer(v);

            while (!queue.isEmpty()) 
            {
            	V w = queue.poll();
                stack.push(w);
            	BetweennessData w_data = vertex_data.get(w);
                
                for (E e : graph.getOutEdges(w))
                {
                	V x = graph.getOpposite(w, e);
                	if (x.equals(w))
                		continue;
                	double wx_weight = weights.transform(e).doubleValue();
                	                	
                	BetweennessData x_data = vertex_data.get(x);
                	double x_potential_dist = w_data.distance + wx_weight;
                	
                    if (x_data.distance < 0) 
                    {
                    	x_data.distance = x_potential_dist;
                      	queue.offer(x);
                    }
                    
                    // note:
                    // (1) this can only happen with weighted edges
                    // (2) x's SP count and incoming edges are updated below 
                    if (x_data.distance > x_potential_dist)
                    {
                    	x_data.distance = x_potential_dist;
                    	// invalidate previously identified incoming edges
                    	// (we have a new shortest path distance to x)
                    	x_data.incomingEdges.clear(); 
                        // update x's position in queue
                    	((MapBinaryHeap<V>)queue).update(x);
                    }

                }
                for (E e: graph.getOutEdges(w))
                {
                	V x = graph.getOpposite(w, e);
                	if (x.equals(w))
                		continue;
                	double e_weight = weights.transform(e).doubleValue();
                	BetweennessData x_data = vertex_data.get(x);
                	double x_potential_dist = w_data.distance + e_weight;
                    if (x_data.distance == x_potential_dist) 
                    {
                        x_data.numSPs += w_data.numSPs;
                        x_data.incomingEdges.add(e);
                    }
                }
            }
    		while (!stack.isEmpty()) 
    		{
    		    V x = stack.pop();

    		    for (E e : vertex_data.get(x).incomingEdges)
    		    {
    		    	V w = graph.getOpposite(x, e);
    		        double partialDependency = 
    		        	vertex_data.get(w).numSPs / vertex_data.get(x).numSPs *
    		        	(1.0 + vertex_data.get(x).dependency);
    		        vertex_data.get(w).dependency +=  partialDependency;
    		    }
    		    if (!x.equals(v)) 
    		    {
    		    	double x_score = vertexScores.get(x).doubleValue();
    		    	x_score += vertex_data.get(x).dependency;
    		    	vertexScores.put(x, x_score);
    		    }
    		}
        }
        vertex_data.clear();
	}
	
	public Double getVertexScore(V v) {
		if(vertexScores == null) initializeAllScores();
		return vertexScores.get(v);
	}

	public String toString(){
		return "EstimatedBetweenness";
	}


}
