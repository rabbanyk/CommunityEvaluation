


void all_levels::print_all_frames_inclusions() {
	
	
	RANGE_loop(i, all_level_frames) {
		
		
		cout<<"LEVEL::: "<<i<<" frames: "<<all_level_frames[i].size()<<" ************************"<<endl;
		
		
		/*if(int(all_level_frames[i].size())!=original_network.size())*/ RANGE_loop(j, all_level_frames[i]) {
			
			cout<<"#"<<j<<"   number_of_nodes: "<<all_level_frames[i][j].number_of_nodes<<" radius: "<<all_level_frames[i][j].radius<<" contains: "<<endl;
			prints(all_level_frames[i][j].contains_these);
			cout<<"original nodes"<<endl;
			original_network.print_id(all_level_frames[i][j].original_nodes, cout);
			
		}
		
		
	}
}



int all_levels::print_all_frames_radii() {
	
	
	
	RANGE_loop(i, all_level_frames) {
		
		char file1[1000];
		sprintf(file1, "%s/frames_%d", folder_output, i);		
		ofstream pout(file1);
		
		
		for(UI j=0; j<all_level_frames[i].size(); j++) {
			
			//cout<<"frame: "<<i<<" "<<all_level_frames[i][j].center_x<<" "<<all_level_frames[i][j].center_y<<endl;
			
			pout<<all_level_frames[i][j].center_x<<" "<<all_level_frames[i][j].center_y<<endl;
			print_circle(pout, all_level_frames[i][j].center_x, all_level_frames[i][j].center_y, all_level_frames[i][j].radius);
			
		}
	}
	
	
	return 0;
	
}



int all_levels::print_level_gdf(int level) {
	
	
	char file1[1000];
	sprintf(file1, "%s/cvis_%d.gdf", folder_output, int(all_level_frames.size()) -1 -level);
	ofstream pout(file1);
	
	
	cout<<"printing level:: "<<level<<" in FILE: "<<file1<<endl;
	
	char file2[1000];
	sprintf(file2, "%s/modules_%d.dat", folder_output, max(1, level));		/*modules_0.dat does not exist*/
	ofstream pout2(file2);
	
	
	if(level==0) 
		return original_network.print_gdf(all_level_frames[level], pout, pout2, 2., original_network, node_mems_across_levels);
	
	
	int_matrix modules_level;
	
	RANGE_loop(i, all_level_frames[level]) {
		modules_level.push_back(all_level_frames[level][i].original_nodes);
	}

	
	netvi frame_network;
	original_network.set_frame_network(frame_network, modules_level);
	
	
	
	
	return frame_network.print_gdf(all_level_frames[level], pout, pout2, 2., original_network, node_mems_across_levels);
	
}






int all_levels::print_levels_gdf() {
	
	
	RANGE_loop(i, all_level_frames) {
		print_level_gdf(i);
	}
	
	return 0;
	
}


void get_link_weight(map<int, deque<pair<int, double> > > & link_weight, int level) {

	link_weight.clear();
	
	char file1[1000];
	sprintf(file1, "%s/cvis_%d.gdf", folder_output, level);
	ifstream gin(file1);
	
	string s;
	bool check_=false;
	while(getline(gin, s)) {
		
		if(s.size()>5 and s[0]=='e' and s[1]=='d' and s[2]=='g' and s[3]=='e' ) {
			
			check_=true;
			
		}
		
		if(check_ and s[0]!='e') {
			
			DD a;
			cast_string_to_doubles(s, a);
			
			
			if(link_weight.find(a[0])==link_weight.end()) {
				deque<pair<int, double> > newd;
				link_weight.insert(make_pair(a[0], newd));
			}
			
			link_weight[a[0]].push_back(make_pair(a[1], a[2]));
		}
	}
	

}

void print_pairs(int a, deque<pair<int, double> > b, ostream & node_prou, set<int> & set_contains, int & edges_) {


	RANGE_loop(i, b) if(set_contains.find(b[i].first)!= set_contains.end() ) {
		node_prou<<a<<","<<b[i].first<<","<<b[i].second<<endl;
		++edges_;
	}



}



void all_levels::reprint_with_single_nodes(DI & nodes, string filename, map<int, deque<pair<int, double> > > & link_weight ) {
	
	// this function is to rewrite everything at the last level
	
	ofstream node_prou(filename.c_str());
	
	node_prou<<"nodedef>name VARCHAR,label VARCHAR, visible BOOLEAN,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR"<<endl;
	
	int level=1;
	
	RANGE_loop(ii, nodes) {
		
		
		int nodei_=nodes[ii];
		frame & sub_frame= all_level_frames[level-1][nodei_];
		
		double unique_value=0.8;
		string color= HSV_to_RGB((sub_frame.h_color_low + sub_frame.h_color_hi)/2, unique_value, unique_value);
		
		if(sub_frame.number_of_nodes==1) {
			
			int nodei=sub_frame.original_nodes[0];
			
			
			double saturat, brightn;
			choose_saturation_and_brightness(node_mems_across_levels[nodei], saturat, brightn);
			
			if(node_mems_across_levels[nodei][node_mems_across_levels[nodei].size()-1]==0)			// if the node is homeless even here it is white
				brightn=1;
			color= HSV_to_RGB((sub_frame.h_color_low + sub_frame.h_color_hi)/2, saturat, brightn);
		}
		
		string new_n;
		original_network.print_most_important_among(sub_frame.original_nodes, new_n);
		double size_node= sub_frame.radius*2;
		node_prou<<original_network.id_of(nodei_)<<",'"<<new_n<<"',"<<"true"<<","<<"true"<<","<<size_node<<","<<size_node<<","<<sub_frame.center_x<<","<<sub_frame.center_y<<","<<color<<endl;
	
	}
	
	
	
	set<int> set_contains;
	DI  new_nonn=nodes;
	if(level==1) {
		original_network.deque_id(new_nonn);
	}
	
	deque_to_set(new_nonn, set_contains);
	node_prou<<"edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE"<<endl;
	int edges_=0;
	RANGE_loop(k, new_nonn) if(link_weight.find(new_nonn[k])!=link_weight.end() ) {
		print_pairs(new_nonn[k], link_weight[new_nonn[k]], node_prou, set_contains, edges_);
	}
	
	
	
}


void all_levels::print_levels_multiple_files() {

	
	int level_renamed=0;
	
	
	
	
	for(int level=int(all_level_frames.size())-1; level>=0; level--) {
		
		//cout<<" > "<<level<<" frames: "<<all_level_frames[level].size()<<endl;
	
	}
	
	map<int, deque<pair<int, double> > > link_weight_last_level;
	get_link_weight(link_weight_last_level, all_level_frames.size() -1);
	
	
	for(int level=int(all_level_frames.size())-1; level>0; level--) {
		
		cout<<"label renamed: "<<level_renamed<<" level: "<<level<<" link_weight: "<<all_level_frames.size() -level<<endl;
		
		map<int, deque<pair<int, double> > > link_weight;
		get_link_weight(link_weight, all_level_frames.size() -level);
		
		deque<frame> & frames=all_level_frames[level];
		
		for(UI ii_=0; ii_<frames.size(); ii_++) {
			
			
			
			char node_prou_c[10000];
			sprintf(node_prou_c, "%s/level_%d_cluster_%d.dat", folder_output, level_renamed, ii_);
			//cout<<node_prou_c<<" <---"<<endl;
			string filename_(node_prou_c);
			ofstream node_prou(node_prou_c);
			
			//node_prou<<"number of nodes: "<<frames[ii_].number_of_nodes<<" size: "<<size_node<<" pos: "<<frames[ii_].center_x<<","<<frames[ii_].center_y<<" color: "<<frames[ii_].color_<<endl;
			//node_prou<<"contains clusters:"<<endl;
			//prints(frames[ii_].contains_these, node_prou);
			node_prou<<"nodedef>name VARCHAR,label VARCHAR, visible BOOLEAN,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR"<<endl;

			RANGE_loop(k, frames[ii_].contains_these) {
				
				frame & sub_frame= all_level_frames[level-1][frames[ii_].contains_these[k]];
				
				double unique_value=0.8;
				string color= HSV_to_RGB((sub_frame.h_color_low + sub_frame.h_color_hi)/2, unique_value, unique_value);
				
				if(sub_frame.number_of_nodes==1) {
					
					int nodei=sub_frame.original_nodes[0];
					
					//cout<<"node: "<<nodei<<endl;
					
					double saturat, brightn;
					choose_saturation_and_brightness(node_mems_across_levels[nodei], saturat, brightn);
					
					if(node_mems_across_levels[nodei][node_mems_across_levels[nodei].size()-1]==0)			// if the node is homeless even here it is white
						brightn=1;
					color= HSV_to_RGB((sub_frame.h_color_low + sub_frame.h_color_hi)/2, saturat, brightn);
				}
				
				string new_n;
				original_network.print_most_important_among(sub_frame.original_nodes, new_n);
				double size_node= sub_frame.radius*2;
				if(level==1)
					node_prou<<original_network.id_of(frames[ii_].contains_these[k])<<",'"<<new_n<<"',"<<"true"<<","<<"true"<<","<<size_node<<","<<size_node<<","<<sub_frame.center_x<<","<<sub_frame.center_y<<","<<color<<endl;
				else
					node_prou<<frames[ii_].contains_these[k]<<",'"<<new_n<<"',"<<"true"<<","<<"true"<<","<<size_node<<","<<size_node<<","<<sub_frame.center_x<<","<<sub_frame.center_y<<","<<color<<endl;
			}
			
			
			//node_prou<<"original nodes: "<<endl;
			//original_network.print_id(frames[ii].original_nodes, node_prou);
			//node_prou<<"ranking:"<<endl;
			set<int> set_contains;
			DI  new_nonn=frames[ii_].contains_these;
			if(level==1) {
				original_network.deque_id(new_nonn);
			}
			
			
			int edges_=0;
			deque_to_set(new_nonn, set_contains);
			node_prou<<"edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE"<<endl;
			RANGE_loop(k, new_nonn) if(link_weight.find(new_nonn[k])!=link_weight.end() ) {
				print_pairs(new_nonn[k], link_weight[new_nonn[k]], node_prou, set_contains, edges_);
			}
			
			//cout<<"edges "<<edges_<<" filename "<<filename_<<endl;
			
			if(edges_==0) {
				
				reprint_with_single_nodes(frames[ii_].original_nodes, filename_, link_weight_last_level);
			
			}
			
			
		}
		
		
		++level_renamed;
		
	}
	

}





