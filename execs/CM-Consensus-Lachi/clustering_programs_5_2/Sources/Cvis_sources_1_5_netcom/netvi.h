




class frame  {
	
	
public:
	
	frame(){ number_of_nodes=-1; radius=-1; center_x=-1; center_y=-1; h_color_low=-1; h_color_hi=-1; };
	~frame(){};
	
	
	double center_x;
	double center_y;
	double radius;
	int number_of_nodes;
	double h_color_low;
	double h_color_hi;
	DI contains_these;
	DI original_nodes;
	
	string color_;
	
};





class netvi : public static_network {
	
public:
	
	netvi():static_network(){};
	~netvi(){};
	int print_gdf(deque<frame> & frames, ostream & pout, ostream & pout2, double radius_factor, netvi & original_network, int_matrix & level__number_of_mems);
	
	void set_frame_network(netvi & frame_network, int_matrix & next_partition);
	void set_frame_network(netvi & frame_network, DI & high_contains_these_i, int_matrix & next_partition);
	
	void edit_positioning(deque<frame*> & Fr);
	void edit_positioning(DI & subset, deque<frame*> & Fr);
	
	double square_distance_node(int i, deque<frame*> & Fr);
	
	bool try_to_swap(int node_one, int node_two, deque<frame*> & Fr);
	double global_square_distance(deque<frame*> & Fr);
	
	void partition_in_words(int_matrix & p, ostream & pout);
	map<int, string> id_word;
	
	
	int check_node_assignments(int_matrix & tool_matrix, set<int> & homeless_nodes);
	int check_node_assignments_power(int_matrix & tool_matrix, set<int> & homeless_nodes);
	
	void print_most_important_among(DI & a, string & newname);
	void deque_id(DI & a) { RANGE_loop(i, a) { a[i]=id_of(a[i]); } };
	
};






#include "some_functions.h"




void netvi::set_frame_network(netvi & frame_network, DI & high_contains_these_i, int_matrix & next_partition) {
	
	module_collection mall(dim);
	
	
	
	RANGE_loop(j, high_contains_these_i) {
		//cout<<"module"<<endl;
		//prints(next_partition[high_contains_these_i[j]]);
		mall.insert(next_partition[high_contains_these_i[j]], 0.01);
		
	}
	
	
	/* now we construct the community network */
	map<int, map<int, pair<int, double> > > neigh_weight_s;
	set_upper_network(neigh_weight_s, mall);
	
	frame_network.set_graph(neigh_weight_s);
	//cout<<"frame network: "<<frame_network.size()<<endl;
	
}



void netvi::set_frame_network(netvi & frame_network, int_matrix & next_partition) {
	
	module_collection mall(dim);
	

	
	RANGE_loop(j, next_partition) {
		mall.insert(next_partition[j], 0.01);

	}
	
	
	/* now we construct the community network */
	map<int, map<int, pair<int, double> > > neigh_weight_s;
	set_upper_network(neigh_weight_s, mall);
	
	frame_network.set_graph(neigh_weight_s);
	
}



int choose_saturation_and_brightness(DI & level__number_of_mems_nodei, double & saturat, double & brightn) {
	
	/*
	 I want to decide the saturat and brightn according to the number of times the node is homeless and overlapping
	 */
	
	int number_of_homeless_times=0;
	int number_of_overlap_times=0;
	RANGE_loop(i, level__number_of_mems_nodei) {
		
		if(level__number_of_mems_nodei[i]==0)
			++number_of_homeless_times;
		if(level__number_of_mems_nodei[i]>1)
			++number_of_overlap_times;
	
	}
	
	double normal_value=0.5;
	double unique_value=0;
	
	if(number_of_overlap_times>0) {		
		unique_value=normal_value + (1.-normal_value) * double(number_of_overlap_times-0.25)/level__number_of_mems_nodei.size();
		saturat=unique_value;
		brightn=unique_value;
		return 0;
	}
	
		
	unique_value=normal_value - (1.-normal_value) * double(number_of_homeless_times)/level__number_of_mems_nodei.size();
	saturat=unique_value;
	brightn=unique_value;
	return 0;

}




void netvi::edit_positioning(DI & subset, deque<frame*> & Fr) {

	
	int stopper=cast_int(pow(dim, paras.edit_positioning_exponent));
	
	// one quick improvement would be to make a more efficient sampling
	
	
	for(int i=0; i<stopper; i++) {
		
		int node_one=subset[irand(subset.size()-1)];
		int node_two=subset[irand(subset.size()-1)];
		
		if(node_one!=node_two)
			try_to_swap(node_one, node_two, Fr);
	}

	
}

void netvi::edit_positioning(deque<frame*> & Fr) {
	
	
	// this function can be used to swap the positions of modules with the same number of nodes
	
	
	int stopper=cast_int(pow(dim, paras.edit_positioning_exponent));
	
	// one quick improvement would be to make a more efficient sampling
	//cout<<"stopper: >"<<stopper<<endl;
	
	for(int i=0; i<stopper; i++) {
		
		int node_one=irand(dim-1);
		int node_two=irand(dim-1);
		
		if(node_one!=node_two)
			try_to_swap(node_one, node_two, Fr);
	}
	
}
 


double netvi::global_square_distance(deque<frame*> & Fr) {
	
	
	double t=0;
	for(int i=0; i<dim; i++)
		t+=square_distance_node(i, Fr);
	
	return t;
}




bool netvi::try_to_swap(int node_one, int node_two, deque<frame*> & Fr) {
	
	double distance_one=square_distance_node(node_one, Fr);
	double distance_two=square_distance_node(node_two, Fr);
	
	swap(Fr[node_one]->center_x, Fr[node_two]->center_x);
	swap(Fr[node_one]->center_y, Fr[node_two]->center_y);
	
	double distance_new_one=square_distance_node(node_one, Fr);
	double distance_new_two=square_distance_node(node_two, Fr);
	
	
	//cout<<"diff: "<<distance_new_one+distance_new_two -(distance_one+distance_two)<<endl;
	
	if(distance_new_one+distance_new_two<distance_one+distance_two)
		return true;
	
	swap(Fr[node_one]->center_x, Fr[node_two]->center_x);
	swap(Fr[node_one]->center_y, Fr[node_two]->center_y);
	
	return false;
	
	
}

void netvi::partition_in_words(int_matrix & p, ostream & pout) {
	
	RANGE_loop(i, p) {
		RANGE_loop(j, p[i]) {
			pout<<id_word[id_of(p[i][j])]<<" ";
		}
		pout<<endl;
	}
	
	
}






void netvi::print_most_important_among(DI & a, string & newname) {
	
	
	deque<pair<double, int> > kp_node;
	RANGE_loop(i, a) {
		kp_node.push_back(make_pair(-vertices[a[i]]->kplus_wsum(a), a[i]));		
	}
	
	
	sort(kp_node.begin(), kp_node.end());
	
	newname.clear();
	
	for(int i=0; i<min(int(kp_node.size()),3); i++) {
		
		int nodei=kp_node[i].second;
		
		newname+=id_word[id_of(nodei)]+" ";
		if(id_word[id_of(nodei)].size()==0) {
			
			char nn[100];
			sprintf(nn, "%d ", id_of(nodei));
			newname+=nn;
		
		}
	}
	
	if(kp_node.size()>=3)
		newname+="...";
	
	/*RANGE_loop(i, kp_node) {
		
		int nodei=kp_node[i].second;
		po<<id_word[id_of(nodei)]<<"-"<<id_of(nodei)<<" "<<-kp_node[i].first<<endl;
		
	}*/
	
	
}


