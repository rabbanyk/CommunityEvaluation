


double count_dash(string s) {
	
	double a=0;
	RANGE_loop(i, s) if(s[i]=='-')
		a++;
	
	return a/2;
}


void convert_set_to_string(set<int> & A, string & c, int i) {
	
	deque<string> AA;
	
	char b[100];
	
	for(set<int>::iterator its=A.begin(); its!=A.end(); its++) {
		sprintf(b, "-%d-", *its);
		string p(b);
		AA.push_back(p);
	}
	
	
	c.clear();
	RANGE_loop(i, AA)
	c+=AA[i];
	
	
}




class not_overlapping_tree {

	
	//		to visualize clusters we need a dendrogram, meaning nodes cannot overlap
	//		information about overlap will be used in considering internal links between clusters within a certain level. 
	
	
public:

	not_overlapping_tree(){};
	~not_overlapping_tree(){};
	

	int read(netvi & original_network, string filename);
	
	
	//  *********************************   DATA   **********************************
	
	deque<int_matrix> short_pts;
	deque<int_matrix> orig_pts_node_no_overlap;		
	// module x= short_pts[1][i] corresponds to module orig_pts_node_no_overlap[1][i]	
	
	/* 
	 short_pts[0] contains the modules at level 0.
	 this means that short_pts[0][0] is the module 0 at level 0 (the lowest)
	 node_x = short_pts[1][0][0] is a node at level 1. this node (node_x) corresponds to the module short_pts[0][node_x]
	 */
	
	int_matrix level__number_of_mems;
	/*
	 level__number_of_mems[x][i]=3 means: node x has three memberships at level i
	 level 0 doesn't mean the nodes! but the first modules
	 for the memberships I am using tool_matrix. so it means that the overlapping super-nodes are not accounted
	 */
	
	/*deque< deque<pair<int, int> >  > mates;
	 call (for some i)
	 x = mates[0][i].first
	 y = mates[0][i].second
	 this means that links between "modules" short_pts[0][x] and short_pts[0][y] are internal
	 */

	

	
private:
	
	deque<int_matrix> orig_pts_node;
	deque<string> pt_names; // pts names
	int dim;


	void get_pts(string filename);
	void set_all_pts(int_matrix & highest_pt, int level, netvi & original_network, bool last_level);
	
	
	void set_short_pts(int_matrix & tool_matrix, int_matrix & previous_pno_overlap, int_matrix & pno_overlap);	
	int separate_overlaps(int_matrix & ten_temp);
	int separate_overlaps_original(int_matrix & ten_temp);
	void get_original_pt(int_matrix & original_part, int_matrix & ten_temp, int_matrix & previous_original_part);
	
	
	void update_level__number_of_mems(int_matrix & tool_matrix);
	void select_homeless_this_level(set<int> & homeless_nodes);

	
};




void not_overlapping_tree::get_pts(string filename) {
	
	
	/*
	 file formats: tp can have whatever labels
	 short_tp_ must have consecutive labels starting from 0
	 homeless nodes present!
	 */
	
	char b[1000];
	char f[1000];
	cast_string_to_char(filename, f);	

	int h=0;
	while(true) {
		
		if(h==0)
			sprintf(b, "%s/tp", f);
			else
				sprintf(b, "%s/short_tp%d", f, h);
				ifstream gin(b);
				if(gin.is_open()) {
					string s(b);
					pt_names.push_back(s);
				}
				else
					break;
		
		++h;
		
	}
	
	cout<<"pts found: "<<pt_names.size()<<endl;
	prints(pt_names);
		
}



bool all_singleton_modules(int_matrix & ten_temp){

	RANGE_loop(i, ten_temp) {	
		if(ten_temp[i].size()>1)
			return false;
	}
	
	return true;
}


int not_overlapping_tree::read(netvi & original_network, string filename) {
	
	dim=original_network.size();
	
	level__number_of_mems.clear();
	for(int i=0; i<dim; i++) {
		DI f;
		level__number_of_mems.push_back(f);
	}
	
	get_pts(filename);

	if(pt_names.size()==0) {
		cout<<"no levels found.. exiting"<<endl;
		return -1;
	}

	for(UI ii=0; ii<pt_names.size(); ii++) {
		
		int_matrix ten_temp;
		cout<<"getting partition from file: "<<pt_names[ii]<<" level: "<<ii<<endl;
		get_partition_from_file(pt_names[ii], ten_temp);
		
		if(ii==0) {
			original_network.translate(ten_temp);
		}
				
		if(ii!=0 && ii==pt_names.size()-1) {		// check if the last partition still contains something
			if(all_singleton_modules(ten_temp))
				return 0;			
		}
		
		
		if(ii==pt_names.size()-1)
			set_all_pts(ten_temp, ii, original_network, true);
		else
			set_all_pts(ten_temp, ii, original_network, false);
	}
	/*cout<<"level__number_of_mems"<<endl;
	printm(level__number_of_mems, cout);*/
	
	return 0;
}


void not_overlapping_tree::select_homeless_this_level(set<int> & homeless_nodes) {

	RANGE_loop(i, level__number_of_mems) if(level__number_of_mems[i][level__number_of_mems[i].size()-1]==0)
		homeless_nodes.insert(i);

}

void not_overlapping_tree::set_all_pts(int_matrix & ten_temp, int level, netvi & original_network, bool last_level) {
	
	// this function is to reset pts in order to have only the modules with the same memberships together
	
	int_matrix tool_matrix;;
	
	if(level==0) {
		orig_pts_node.push_back(ten_temp);
		tool_matrix= ten_temp;
	}
	else {
		
		int_matrix original_part;
		int_matrix & previous_original_part=orig_pts_node[orig_pts_node.size()-1];
		
		get_original_pt(original_part, ten_temp, previous_original_part);
		
		orig_pts_node.push_back(original_part);
		
		//cout<<"-> level: "<<level<<endl;
		//printm(ten_temp);
		
		separate_overlaps(ten_temp);	// this is here because there is a difference between overlapping nodes and modules
		
		get_original_pt(tool_matrix, ten_temp, previous_original_part);		// so in tool_matrix overlapping modules are now separated
		
	}
	

	// here you can compute the homeless and overlapping nodes at this level
	update_level__number_of_mems(tool_matrix);
	//cout<<"-> tool: "<<level<<endl;
	//printm(tool_matrix);
	separate_overlaps(tool_matrix);
	//cout<<"-> tool after: "<<level<<endl;
	//printm(tool_matrix);

	
	
	set<int> homeless_nodes;
	select_homeless_this_level(homeless_nodes);
	original_network.check_node_assignments(tool_matrix, homeless_nodes);
	
	if(last_level)
		original_network.check_node_assignments_power(tool_matrix, homeless_nodes);

	orig_pts_node_no_overlap.push_back(tool_matrix);
	//cout<<"orig_pts_node_no_overlap after ASSIGNING"<<endl;
	//printm(orig_pts_node_no_overlap[orig_pts_node_no_overlap.size()-1]);
	
	if(level==0) {
		short_pts.push_back(tool_matrix);
	}
	else {
		int_matrix & previous_pno_overlap= orig_pts_node_no_overlap[orig_pts_node_no_overlap.size()-2];
		int_matrix & pno_overlap= orig_pts_node_no_overlap[orig_pts_node_no_overlap.size()-1];
		set_short_pts(tool_matrix, previous_pno_overlap, pno_overlap);
		short_pts.push_back(tool_matrix);
	}
	
	
	/*cout<<"short_pt last"<<endl;
	printm(tool_matrix);
	
	cout<<"orig_pts_node_no_overlap last "<<endl;
	printm(orig_pts_node_no_overlap[orig_pts_node_no_overlap.size()-1]);*/


	
	//cout<<"************** END *********************************"<<endl;
	
}

void not_overlapping_tree::set_short_pts(int_matrix & tool_matrix, int_matrix & previous_pno_overlap, int_matrix & pno_overlap) {
	
	/*cout<<"PNO overlap"<<endl;
	printm(pno_overlap);
	
	cout<<"PNO overlap pre"<<endl;
	printm(previous_pno_overlap);*/
	
	module_collection mall(dim);
	for(UI i=0; i<previous_pno_overlap.size(); i++)
		mall.insert(previous_pno_overlap[i], 0.001);
	
	tool_matrix.clear();
	RANGE_loop(i, pno_overlap) {
		
		set<int> c;
		RANGE_loop(j, pno_overlap[i]) {
			int mem_node=*(mall.memberships[pno_overlap[i][j]].begin() );
			c.insert(mem_node);
		}
		
		DI f;
		set_to_deque(c, f);
		tool_matrix.push_back(f);
	}
	
	
	
	
	separate_overlaps(tool_matrix);
	//cout<<"short_pt separated"<<endl;
	//printm(tool_matrix);

	
	// this happens if the pno_overlap does not perfectly cover previous_pno_overlap
	// in this case I sepated the overlap. now I reset pno_overlap
	get_original_pt(pno_overlap, tool_matrix, previous_pno_overlap);
	
	
	
	
	
	
	
}



int not_overlapping_tree::separate_overlaps_original(int_matrix & ten_temp) {
	
	
	
	//cout<<"dim: "<<dim<<endl;
	
	module_collection mall(dim);
	for(UI i=0; i<ten_temp.size(); i++)
		mall.insert(ten_temp[i], 0.001);
	
	
	map<string, int> mems_module;
	map<int, string> mems;
	
	RANGE_loop(i, mall.memberships) if(mall.memberships[i].size()>0) {
		
		string s;
		convert_set_to_string(mall.memberships[i], s, i);
		mems[i]=s;
		mems_module.insert(make_pair(s, mems_module.size()));
	}
	
	
	ten_temp.clear();
	
	DI first_;
	RANGE_loop(i, mems_module) {
		ten_temp.push_back(first_);
	}
	
	
	for(map<int, string>:: iterator itm=mems.begin(); itm!=mems.end(); itm++) {
		ten_temp[mems_module[itm->second]].push_back(itm->first);
	}
	
	return 0;
	
	
}

int not_overlapping_tree::separate_overlaps(int_matrix & ten_temp) {
	
	
	
	//return separate_overlaps_original(ten_temp);
	//cout<<"dim: "<<dim<<endl;
	
	module_collection mall(dim);
	for(UI i=0; i<ten_temp.size(); i++)
		mall.insert(ten_temp[i], 0.001);
	
	
		
	DI node_community;
	for(int i=0; i<dim; i++)
		node_community.push_back(-2);
	
	
	RANGE_loop(i, mall.memberships) {
		if(mall.memberships[i].size()==1) node_community[i]=*(mall.memberships[i].begin());
		else if(mall.memberships[i].size()>1) node_community[i]=-1;
	}
	
	int maxmem=0;
	RANGE_loop(i, node_community) maxmem=max(maxmem, node_community[i]);
	++maxmem;
	
	bool overlap_exists=false;
	RANGE_loop(i, node_community) if(node_community[i]==-1) { node_community[i]=maxmem; overlap_exists=true; }
	
	
	
	set_partition_from_memberships(node_community, ten_temp);
	
	if(overlap_exists) {
		
		DI ff=ten_temp[ten_temp.size()-1];
		ten_temp.pop_back();
		
		
		RANGE_loop(i, ff) {
			
			DI pu;
			pu.push_back(ff[i]);
			ten_temp.push_front(pu);
			//ten_temp.push_back(pu); 
		
		}
		
		
	}
	
	
	return 0;
	
}

void not_overlapping_tree::get_original_pt(int_matrix & original_part, int_matrix & ten_temp, int_matrix & previous_original_part) {
	
	original_part.clear();
	
	RANGE_loop(i, ten_temp) {
		
		set<int> a;
		RANGE_loop(j, ten_temp[i]) {
			
			deque_to_set_app(previous_original_part[ten_temp[i][j]], a);
		}
		
		DI b;
		set_to_deque(a, b);
		original_part.push_back(b);
	}
	
}



void not_overlapping_tree::update_level__number_of_mems(int_matrix & tool_matrix) {
	
	module_collection mall(dim);
	for(UI i=0; i<tool_matrix.size(); i++)
		mall.insert(tool_matrix[i], 0.001);
	
	mall.put_gaps();
	
	RANGE_loop(i, mall.memberships) {
		level__number_of_mems[i].push_back(mall.memberships[i].size());
	}

}


