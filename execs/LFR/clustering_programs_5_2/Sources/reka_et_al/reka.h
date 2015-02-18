


#include "static_network.h"

class L_module {

	
public:
	
	
	L_module(double a) {  nc=1; kin=0; ktot=a; };
	~L_module(){};
	
	int nc;
	double kin;
	double ktot;
	
};



typedef map<int, double > mapip;
typedef map<int, L_module> map_int_om;



void prints(map_int_om & M) {

	
	for(map_int_om :: iterator itm = M.begin(); itm!=M.end(); itm++)
		cout<<"module: "<<itm->first<<"\t\t\t\tnc= "<<itm->second.nc<<"\t ktot= "<<itm->second.ktot<<"\t kin= "<<itm->second.kin<<endl;


}




class siglouvain : public static_network {
	
public:
	
		
	siglouvain(): static_network(){};
	~siglouvain(){};
	
	
	double collect_raw_groups_once(deque<deque<int> > & );
	double collect_raw_groups(int_matrix & best_partition, int reals);

	
private:
	
	int single_pass_weighted();	
	void set_partition_collected(deque<deque<int> > & ten2);
	
	
	
	// data
	DI vertex_order;
	DI vertex_label;
	int nodes_changed;
	
	

};






int siglouvain::single_pass_weighted() {
	
	
	nodes_changed=0;
	shuffle_s(vertex_order);		
	
	//cout<<"vertex_order:"<<endl;
	//prints(vertex_order);
	
	RANGE_loop(i, vertex_order) {
		
		
		//*********************************************************************************************
		int node=vertex_order[i];
		mapid Mlabel_kin;		// M is a map module_label -> kin (internal weight)
		
		for(int j=0; j<vertices[node]->links->size(); j++)
			int_histogram(vertex_label[vertices[node]->links->l[j]], Mlabel_kin, vertices[node]->links->w[j]);
		
		//cout<<"Mlabel_kin"<<endl;
		//prints(Mlabel_kin);
		
		DI most_frequent;
		double w_max=0;
		
		for(mapip:: iterator itM= Mlabel_kin.begin(); itM!=Mlabel_kin.end(); itM++) {
			
			if (itM->second> w_max *(1.000001) ) {
				most_frequent.clear();
				most_frequent.push_back(itM->first);
				w_max=itM->second;
			}
			else if(fabs(itM->second-w_max)<1e-7*w_max) {
				most_frequent.push_back(itM->first);
			}
		}
		
		//cout<<"most_frequent "<<endl;
		//prints(most_frequent);
		
		if(most_frequent.size()>0) {
			
			int new_label=most_frequent[irand(most_frequent.size()-1)];
			if(new_label!=vertex_label[node]) {
				++nodes_changed;
				vertex_label[node]=new_label;
			}
		}
		//*********************************************************************************************

	}
	
	
	set<int> s;
	deque_to_set(vertex_label,s);
	return s.size();
	
}








double siglouvain::collect_raw_groups_once(deque<deque<int> > & P) {

	
	vertex_label.clear();
	vertex_order.clear();
	for(int i=0; i<dim; i++) {
		vertex_label.push_back(i);
		vertex_order.push_back(i);
	}
	
	int stopper=0;	
	int iteration=0;
	int previous_nodes_changed=0;
	
	
	while (true) {
		
		int number_of_modules=single_pass_weighted();
				
		if(paras.print_flag_subgraph && iteration%1==0)
			cout<<"iteration: "<<iteration<<" number of modules: "<<number_of_modules<<"\t\t\t changing nodes...  "<<nodes_changed<<" vs "<<previous_nodes_changed<<endl;
		
		++iteration;
		
		/* this conditions means that at least max_iteration_convergence iterations are done and, after that, the number of nodes changed has to decrease (up to a facto 1e-3) */
		if(double(nodes_changed - previous_nodes_changed) > 1e-3 * previous_nodes_changed && iteration>paras.max_iteration_convergence)
			stopper++;
		
		if(stopper==paras.max_iteration_convergence || nodes_changed==0 || iteration==dim)
			break;
		
		previous_nodes_changed=nodes_changed;		
	}
	
	
	set_partition_collected(P);
	
	
	vertex_label.clear();
	vertex_order.clear();
	nodes_changed=0;
	
	return newman_modularity(P);
	
}







double siglouvain::collect_raw_groups(int_matrix & best_partition, int reals=paras.Or) {
	
	
	double qmax=-1e100;
	best_partition.clear();
	
	for (int r=0; r<reals; r++) {
		
		int_matrix P;
		double q=collect_raw_groups_once(P);
		cout<<"qs: "<<q<<" run "<<r<<" done!"<<endl<<endl<<endl;
;
		if (q>qmax) {
			qmax=q;
			best_partition=P;
		}
	}
	
	return qmax;
}





void siglouvain::set_partition_collected(deque<deque<int> > & ten2) {
	
	
	ten2.clear();
	
	deque<deque<int> > M;
	
	
	// take partition from vertex_label  //******************************
	map<int, int> mems;
	for(int i=0; i<dim; i++) {
		
		
		pair<map<int, int>::iterator, bool>  itm_bool= mems.insert(make_pair(vertex_label[i], mems.size()));
		
		if(itm_bool.second==true) {
			deque<int> first;
			M.push_back(first);
		}
		
		
		M[itm_bool.first->second].push_back(i);
	}
	
	
	
	// check if subgraphs are connected  //******************************
	RANGE_loop(i, M) {
		
		deque<DI> link_per_node;
		deque<DD> weights_per_node;
		set_subgraph(M[i], link_per_node, weights_per_node);
		static_network giovanni;
		giovanni.set_graph(link_per_node, weights_per_node, M[i]);
		deque<deque<int> > gM;
		giovanni.set_connected_components(gM);
		
		if(gM.size()==1)
			ten2.push_back(M[i]);
		else {
			for(int j=0; j<int(gM.size()); j++) {
				giovanni.deque_id(gM[j]);
				ten2.push_back(gM[j]);
			}
		}
		
		
	}
	
	
}
