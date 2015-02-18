


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
	double compute_rand_mod(int reals);
	void minimality(int_matrix & P);
	void check_merging(int_matrix & P);
	void check_merging_pairs(int_matrix & P);
	double collect_raw_groups(int_matrix & best_partition, int reals);

	
private:
	
	void module_initializing();
	void weighted_favorite_of(const int & node);	
	void single_pass_weighted();	
	void set_partition_collected(deque<deque<int> > & ten2);
	bool together_(DI & a, int_matrix & minimal_groups);
	
	
		
	double qf(const double kin_g, const double ktot_g) {
		return (kin_g - ktot_g*ktot_g/(two_M));
	};
		
		
		
		
	//int check_all();
	map<int, L_module> label_module;
	deque<int> vertex_label;
	deque<int> vertex_order;
	deque<bool> vertex_to_check;
	deque<bool> vertex_to_check_next;
	int nodes_changed;
	
	
	double qtot;
	double two_M;
	

};




void siglouvain::module_initializing() {
	
	qtot=0;
	two_M= 2 *tstrength;
	
	for(int i=0; i<dim; i++) {
		
		vertex_label.push_back(i);
		vertex_order.push_back(i);
		vertex_to_check.push_back(true);
		vertex_to_check_next.push_back(false);
		L_module newm(vertices[i]->strength);
		label_module.insert(make_pair(i, newm));
		qtot+=qf(0, vertices[i]->strength);
		
	}
}


void siglouvain::weighted_favorite_of(const int & node) {
	
	
	//cout<<qtot/two_M<<" <<<---- "<<endl;
	
	mapip Mlabel_kin;		// M is a map module_label -> kin (internal weight)
	
	for(int j=0; j<vertices[node]->links->size(); j++)
		int_histogram(vertex_label[vertices[node]->links->l[j]], Mlabel_kin, vertices[node]->links->w[j]);
	
	double dq1=-qf(0, vertices[node]->strength);	// at the beginning the node is alone
	double best_dq2=-1e200;
	int module_to_move=-1;	
	double kop=0;		// internal degree respect with the old module
	double kp=0;		// internal degree respect with the new module
	
	//cout<<"node: "<<node<<endl;
	//cout<<"dq1 "<<dq1<<endl;
	//prints(Mlabel_kin);
	
	for(mapip:: iterator itM= Mlabel_kin.begin(); itM!=Mlabel_kin.end(); itM++) {
				
		map_int_om :: iterator itOM = label_module.find(itM->first);
		
		if(itM->first==vertex_label[node]) {
						
			// this can be optimized a little bit, I suppose
			
			double kin_old=itOM->second.kin;
			double ktot_old=itOM->second.ktot;
			double qpartial1 = qf(kin_old, ktot_old);
			
			kin_old-=2*itM->second;
			ktot_old-=vertices[node]->strength;
			double qpartial1_new = qf(kin_old, ktot_old);
			dq1=qpartial1_new - qpartial1;
			kop=itM->second;
			
			
		} else {
			
			double kin_old=itOM->second.kin;
			double ktot_old=itOM->second.ktot;
			double qpartial1 = qf(kin_old, ktot_old);
			
			kin_old+=2*itM->second;
			ktot_old+=vertices[node]->strength;
			double qpartial1_new = qf(kin_old, ktot_old);
			double dq2=qpartial1_new - qpartial1;
			if (dq2 + 1e-10 * (ran4()-0.5) > best_dq2) {
				best_dq2=dq2;
				module_to_move=itM->first;
				kp=itM->second;
			}
			
		}
		
	}
	
	
	
	// updates modules
	if(dq1+best_dq2>0) {
		
		qtot+=dq1+best_dq2;
		nodes_changed++;
		
		for(int j=0; j<vertices[node]->links->size(); j++)
			vertex_to_check_next[vertices[node]->links->l[j]]=true;
		
		map_int_om :: iterator itm= label_module.find(vertex_label[node]);
		
		--(itm->second.nc);

		
		if(itm->second.nc==0)
			label_module.erase(itm);
		else {
			itm->second.kin-= 2 * kop;
			itm->second.ktot-= vertices[node]->strength;
		}
		

		itm= label_module.find(module_to_move);
		++(itm->second.nc);
		itm->second.kin+= 2 * kp;
		itm->second.ktot+= vertices[node]->strength;
		

		vertex_label[node] = module_to_move;
	}
	
	
}



void siglouvain::single_pass_weighted() {
	
	
	for(deque<int> :: iterator itd=vertex_order.begin(); itd!=vertex_order.end(); itd++) {		
		if(vertex_to_check[*itd]==true) {
			weighted_favorite_of(*itd);
			//check_all();
		}
		
		
	}
	
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



double siglouvain::collect_raw_groups(int_matrix & best_partition, int reals=paras.Or) {
	
	
	double qmax=-1e100;
	best_partition.clear();
	
	for (int r=0; r<reals; r++) {
		
		int_matrix P;
		double q=collect_raw_groups_once(P);
		if (q>qmax) {
			qmax=q;
			best_partition=P;
		}
	}
	
	return qmax;
}




double siglouvain::collect_raw_groups_once(deque<deque<int> > & P) {

	module_initializing();
	
	int stopper=0;	
	int previous_nodes_changed=dim;
	int iteration=0;
	
	while (true) {
		
		
		nodes_changed=0;
		
		for(int i=0; i<int(vertex_to_check.size()); i++)
			vertex_to_check_next[i]=false;
		shuffle_s(vertex_order);
		
		single_pass_weighted();
				
		if(paras.print_flag_subgraph && iteration%20==0)
			cout<<"iteration: "<<iteration<<" number of modules: "<<label_module.size()<<" q: "<<qtot/two_M<<endl;
		
		++iteration;
		
		/* this conditions means that at least max_iteration_convergence iterations are done and, after that, the number of nodes changed has to decrease (up to a facto 1e-3) */
		if(double(nodes_changed - previous_nodes_changed) > 1e-3 * previous_nodes_changed && iteration>paras.max_iteration_convergence)
			stopper++;
		
		if(stopper==paras.max_iteration_convergence || nodes_changed==0 || iteration==dim)
			break;
		
		
		vertex_to_check=vertex_to_check_next;
		previous_nodes_changed=nodes_changed;		
	}
	
	set_partition_collected(P);
	
	if(paras.print_flag_subgraph)
		cout<<"collection done "<<endl<<endl<<endl;
	
	label_module.clear();
	vertex_label.clear();
	vertex_order.clear();
	vertex_to_check.clear();
	vertex_to_check_next.clear();
	nodes_changed=0;
	
	return qtot/(2*edges());
	
}


double siglouvain::compute_rand_mod(int reals=10) {
	
	paras.print_flag_subgraph=false;
	
	double max_mod=-1e100;
	
	/*for (int r=0; r<reals; r++) {
		
		int_matrix ten;
		DI labs;
		set_rewired_matrix(ten, labs);
		siglouvain net2;
		net2.set_graph(ten, labs);
		int_matrix P;
		double qmod=net2.collect_raw_groups_once(P);
		max_mod=max(qmod, max_mod);
		
	}*/
	
	double sqrt_ave=0;
	double ave=0;
	
	
	RANGE_loop(i, vertices) {		
		sqrt_ave+=sqrt(vertices[i]->strength);
		ave+=vertices[i]->strength;
	}
		
	sqrt_ave/=dim;
	ave/=dim;
	
	double	qr=0.765 * sqrt_ave/ave;
	cout<<" "<<max_mod<<" "<<qr<<endl;
	return qr;
	
}



bool siglouvain::together_(DI & a, int_matrix & minimal_groups) {
	
	deque<DI> link_per_node;
	deque<DD> weights_per_node;
	set_subgraph(a, link_per_node, weights_per_node);
	siglouvain giovanni;
	giovanni.set_graph(link_per_node, weights_per_node, a);
	
	minimal_groups.clear();
	double qreal_mod= giovanni.collect_raw_groups(minimal_groups);
	double qrand=giovanni.compute_rand_mod();
	
	cout<<"qreal_mod: "<<qreal_mod<<" qrand: "<<qrand<<" density: "<<kin(a)/(a.size()*(a.size()-1))<<" dim "<<giovanni.size()<<endl;
	if(qreal_mod>qrand) {
		cout<<"split: "<<minimal_groups.size()<<endl;
		return false;
	}
	
	
	
	return true;
	


}

void erase_empty(int_matrix & P) {

	int_matrix B;
	RANGE_loop(i, P) {
		if(P[i].size()!=0)
			B.push_back(P[i]);
	}
	P=B;
	

}


void siglouvain::minimality(int_matrix & P) {
	
	/*
		this function is to check if each module is minimal
	 */
		
	int_matrix to_check_again;
	
	RANGE_loop(i, P) {
		
		
		int_matrix minimal_groups;
		if(together_(P[i], minimal_groups)==false) {
			P[i].clear();
			RANGE_loop(j, minimal_groups) {
				P.push_back(minimal_groups[j]);
			}
		
		}
	}
	
	erase_empty(P);
}



void siglouvain::check_merging(int_matrix & P) {
	
	map<int, map<int, double > > neigh_weight_f;
	set_upper_network(neigh_weight_f, P);
	siglouvain giovanni;
	giovanni.set_graph(neigh_weight_f);
	
	int_matrix big_groups;
	double qreal_mod= giovanni.collect_raw_groups(big_groups);
	printm(big_groups);
	
	cout<<"qreal_mod fusion: "<<qreal_mod<<" "<<giovanni.size()<<endl;
	giovanni.draw("upnet.dat");
	
	RANGE_loop(i, big_groups) {
		
		set<int> all;
		RANGE_loop(j, big_groups[i]) {
			deque_to_set_app(P[big_groups[i][j]], all);
		}
		DI alld;
		set_to_deque(all, alld);
		int_matrix minimal_groups;
		if(together_(alld, minimal_groups)) {
			
			cout<<"merging!"<<endl;
			RANGE_loop(j, big_groups[i]) {
				P[big_groups[i][j]].clear();
			}
			P.push_back(alld);
		}
	}
	
	
	erase_empty(P);
	
}


void siglouvain::check_merging_pairs(int_matrix & P) {
	
	map<int, map<int, double > > neigh_weight_f;
	set_upper_network(neigh_weight_f, P);
	siglouvain giovanni;
	giovanni.set_graph(neigh_weight_f);
	
	set<pair<int, int> > couples_to_merge;
	RANGE_loop(i, giovanni.vertices) {
		
		double maxw=0;
		int fav=0;
		for(int h=0; h<giovanni.vertices[i]->links->size(); h++) if(giovanni.vertices[i]->links->w[h]> maxw) {
			maxw=giovanni.vertices[i]->links->w[h];
			fav=giovanni.vertices[i]->links->l[h];
		}
		
		couples_to_merge.insert(make_pair(min(int(i), fav), max(int(i), fav)));
		cout<<" -> "<<min(int(i), fav) <<"  ...   "<<max(int(i), fav) <<endl;
	
	}
	
	
	set<int> already_merged;
	for(set<pair<int, int> >:: iterator its = couples_to_merge.begin(); its!=couples_to_merge.end(); its++) if(already_merged.find(its->first)==already_merged.end() 
																											and already_merged.find(its->second)==already_merged.end() ) {
		
		
		/*cout<<"--------------------"<<endl;
		print_id(P[its->first], cout);
		print_id(P[its->second], cout);*/
		set<int> all;
		deque_to_set_app(P[its->first], all);
		deque_to_set_app(P[its->second], all);
		DI alld;
		set_to_deque(all, alld);
		
		
		int_matrix minimal_groups;
		if(together_(alld, minimal_groups)) {
			already_merged.insert(its->first);
			already_merged.insert(its->second);
			cout<<"merging!"<<endl;
			P[its->first].clear();
			P[its->second].clear();
			P.push_back(alld);
		}
	}
	
	cout<<"already_merged: "<<already_merged.size()<<endl;
	erase_empty(P);
	
}



/*
int siglouvain::check_all() {
	
	
	bool print_stuff_oslom_local= false;
	
	deque<deque<int> > M;
	set_partition_collected(M);
	
	
	
	double one=0;
	double two=0;
	
	
	for(map_int_om :: iterator itm = label_module.begin(); itm!=label_module.end(); itm++) {
		one+=sin(itm->second.nc) + log(itm->second.ktot) + cos(itm->second.kin);
		if(print_stuff_oslom_local)
			cout<<"nc: "<<itm->second.nc<<" kin: "<<itm->second.kin<<" ktot: "<<(itm->second.ktot)<<" -----------<"<<endl;
	}
	
	
	for(UI i=0; i<M.size(); i++) {
		
		double _ktot= ktot(M[i]);
		double _kin= kin(M[i]);
		if(print_stuff_oslom_local)
			cout<<"nc: "<<M[i].size()<<" kin: "<<_kin<<" ktot: "<<_ktot<<" "<<i<<endl;
		two+=sin(M[i].size()) + log( _ktot) + cos(_kin);
		
		
	}
	
	
	cherr(one-two, 1e-6);
	if(print_stuff_oslom_local)
		cout<<"one, two: "<<one<<" "<<two<<endl;
	
	cout<<"check:: "<<qtot/two_M<<" "<<" "<<qtot/two_M-newman_modularity(M)<<endl;
	return 0;
	

}

*/
