

double all_levels::minimize_distance_frame(deque<frame*> & Fr, DI & ranking, netvi & frame_network) {
	
	// ranking must contain the names of the frame in the right order
	// in Fr there are some frames. Each frame corresponds to a group of nodes
	// which are written in modules_correspondent_to_these_frames
	
	
	frame_network.edit_positioning(Fr);
	
	DD distances_from_the_center;
	RANGE_loop(i, Fr) {
		distances_from_the_center.push_back(Fr[i]->center_x*Fr[i]->center_x + Fr[i]->center_y*Fr[i]->center_y);
	}	
	
	set_ranking_size_based(ranking, distances_from_the_center, 1);	
	return frame_network.global_square_distance(Fr);
}

void all_levels::frame_location(deque<frame*> & Fr, DI & ranking) {
	
	DD nodes_x;
	DD nodes_y;
	
	DD radius_sizes;
	RANGE_loop(j, Fr) {
		radius_sizes.push_back(0.);
	}
	RANGE_loop(i, Fr) {
		radius_sizes[i]=Fr[ranking[i]]->radius;
	}
	
	// radius_sizes contains the radii of the modules sorted from the bulk to the skirts
	locate_nodes(nodes_x, nodes_y, radius_sizes);
	
	
	for(UI i=0; i<nodes_x.size(); i++) {
		
		Fr[ranking[i]]->center_x=nodes_x[i];
		Fr[ranking[i]]->center_y=nodes_y[i];
		
	}
	
	
	set_zero_baricenter(Fr);
	set_colors(Fr, paras.color_reduction_factor);	// should be set in according to the ranking  ********* IMPORTANT ***********
	
	//cout<<"set single_frame_positions done"<<endl;
}



int all_levels::swap_homeless_frames(DI & homeless_frames, deque<frame*> & Fr) {
	
	
	if(homeless_frames.size()==0)
		return 0;
	
	netvi frame_network;
	int_matrix modules_correspondent_to_these_frames;
	RANGE_loop(i, Fr) {
		modules_correspondent_to_these_frames.push_back(Fr[i]->original_nodes);
	}
	
	
	original_network.set_frame_network(frame_network, modules_correspondent_to_these_frames);	
	frame_network.edit_positioning(homeless_frames, Fr);

	return 0;

}



int all_levels::set_single_frame_positions_no_homeless(deque<frame*> & Fr, int level) {
	
	// idea 1: removing homeless nodes in this frame and add them later (and far)
	// idea 2: multiply links between friend frames
	
	/*
	 idea 1:
	 step a) remove homeless nodes from Fr -> Fr_nohomeless and compute the ranking for them
	 step b) do the usual things to locate Fr_nohomeless
	 step c) resume Fr. add the homeless at the end of the ranking and build up a special function to swap ONLY the new nodes	 
	 */
	
	
	DI ranking;
	
	if(level==0)
		return set_single_frame_positions(Fr, ranking);
		

	// create a new frame without homeless nodes

	deque<frame*> Fr_nohomeless;
	
	DI homeless_frames;
	map<int, int> new_names;
	
	RANGE_loop(i, Fr) {
	
		if(Fr[i]->number_of_nodes>1 || node_mems_across_levels[Fr[i]->original_nodes[0]][level-1]>0) {
			
			Fr_nohomeless.push_back(Fr[i]);
			new_names[Fr_nohomeless.size()-1]=i;
			
			//cout<<i<<" level: "<<level<<" module:"<<endl;
			//original_network.print_id(Fr[i]->original_nodes, cout);
		
		}
		else {
			homeless_frames.push_back(i);
		}
	}
	
	
	
	// set ranking for Fr_nohomeless
	set_single_frame_positions(Fr_nohomeless, ranking);
	
	
	/*cout<<"no homeless->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>: "<<Fr_nohomeless.size()<<endl;
	cout<<"homeless frames: "<<homeless_frames.size()<<endl;
	prints(homeless_frames);
	cout<<"no homeless"<<endl;
	prints(new_names);*/
	
	RANGE_loop(i, ranking) {
		ranking[i]=new_names[ranking[i]];
	}
	
	//cout<<"ranking before"<<endl;
	//prints(ranking);
	
	RANGE_loop(i, homeless_frames) {
		ranking.push_back(homeless_frames[i]);
	}
	
	//cout<<"ranking after"<<endl;
	//prints(ranking);
	
	deque<frame> fakes;
	int ho=homeless_frames.size();
	for(int i=0; i<ho*paras.percentage_fake_homeless; i++) {
	
		frame fake;
		fake.number_of_nodes=0;
		fake.radius=paras.node_radius;
		fakes.push_back(fake);
		
		Fr.push_back(&(fakes[fakes.size()-1]));
		
		int new_frame_name=ranking.size();
		
		ranking.push_back(new_frame_name);
		homeless_frames.push_back(new_frame_name);
	}
	
	//cout<<"frame location"<<endl;
	//prints(ranking);
	frame_location(Fr, ranking);
	
	//cout<<"swapping homeless nodes: "<<ho<<endl;
	//prints(homeless_frames);
	
	
	
	swap_homeless_frames(homeless_frames, Fr);
	
	

	return 0;
	
}

int all_levels::set_single_frame_positions(deque<frame*> & Fr, DI & ranking) {
	
	
	// this functions is to set the positions of a given frame, after you set number_of_nodes and radii
	
	
	if(Fr.size()==0)
		return -1;
	
	DD module_radius;
	RANGE_loop(j, Fr) {
		module_radius.push_back(Fr[j]->radius);
	}
	
	ranking.clear();
	set_ranking_size_based(ranking, module_radius);

	
	double total_square_length= 1e+300;

	netvi frame_network;	// one idea is that frame_network should account for the friend frames. meaning if they are actually overlapping
	int_matrix modules_correspondent_to_these_frames;
	RANGE_loop(i, Fr) {
		modules_correspondent_to_these_frames.push_back(Fr[i]->original_nodes);
	}
	
	
	original_network.set_frame_network(frame_network, modules_correspondent_to_these_frames);
	
	while (true) {
		
		double total_square_length_temp=total_square_length;
		
		frame_location(Fr, ranking);	
		
		if(Fr.size()>2)
			total_square_length_temp=minimize_distance_frame(Fr, ranking, frame_network);
		
		//cout<<"total_square_length_temp: "<<total_square_length_temp<<" "<<Fr.size()<<endl;
		
		if(total_square_length_temp>=total_square_length*0.9) {
			
			if(total_square_length!=total_square_length_temp)
				frame_location(Fr, ranking);
			break;			
		}
		total_square_length=total_square_length_temp;
		
	}
	
	return 0;
	
}

int all_levels::set_positions(UI current_level) {
	
	
	
	cout<<"############### position level: "<<current_level<<endl;
	
	// this function loops on one higher level and locates nodes with positions relative to the previous centers
	
	UI higher_level= current_level+1;
	
	if(higher_level==all_level_frames.size())		// highest level. no higher level to loop on
		return highest_level_set_positions();
	
	RANGE_loop(i, all_level_frames[current_level]) {
		all_level_frames[current_level][i].center_x=-1;
		all_level_frames[current_level][i].center_y=-1;
	}
	
	RANGE_loop(i, all_level_frames[higher_level]) {		// loop on higher level
				
		DI & frames_to_put_together=all_level_frames[higher_level][i].contains_these;
		
		deque<frame*> Fr;
		RANGE_loop(j, frames_to_put_together) {	// contains_these is a module of current level
			Fr.push_back(& all_level_frames[current_level][frames_to_put_together[j]]);
		}
		
		// Fr groups all the current_level frames which are together at higher level
		set_single_frame_positions_no_homeless(Fr, current_level);
	}
	
	return 0;
}

int all_levels::highest_level_set_positions() {
	
	deque<frame*> Fr;
	
	UI current_level=all_level_frames.size()-1;
	
	RANGE_loop(i, all_level_frames[current_level])	{
		Fr.push_back(& all_level_frames[current_level][i]);
	}
	
	
	
	set_single_frame_positions_no_homeless(Fr, current_level);	
	
	/*DD module_radius;
	RANGE_loop(j, Fr) {
		module_radius.push_back(Fr[j]->radius);
	}
	
	DI ranking;
	set_ranking_size_based(ranking, module_radius);
	frame_location(Fr, ranking);*/
		
	return 0;
	
}

void all_levels::set_radii(frame* A, int level) {
	
	// this function works supposing contains_these have the right radii
	// level is the level A belongs to
	
	double max_distance=0;
	
	//cout<<"reset_radii level: "<<level<<" nodes: "<<A->number_of_nodes<<endl;
	
	
	RANGE_loop(i, A->contains_these) {
		
		double & x=all_level_frames[level-1][A->contains_these[i]].center_x;
		double & y=all_level_frames[level-1][A->contains_these[i]].center_y;
		
		double dist=sqrt(x*x+y*y)+all_level_frames[level-1][A->contains_these[i]].radius;
		max_distance=max(max_distance,dist);
		
		//cout<<"max_distance: "<<max_distance<<endl;
		
	}
	
	A->radius=max_distance;
}


int all_levels::set_radii_and_positions() {
	
	
	// this function is to compute the right radii of all the frames
	// it requires that number_of_nodes has been set and also the radii of the lowest level
	
	if(all_level_frames.size()<2)
		return -1;
		
	set_positions(0); // set positions of single nodes. 
	
	for(UI i=1; i<all_level_frames.size(); i++) {
		
		cout<<"************************* Setting radii and positions: LEVEL "<<i<<" modules: "<<all_level_frames[i].size()<<endl;
		
		RANGE_loop(j, all_level_frames[i])
			set_radii(&all_level_frames[i][j], i);
		
		set_positions(i);
	}
	
	return 0;
}



int all_levels::reset_frames() {
	
	for(int level=all_level_frames.size()-1; level>0; level--) {
		
		cout<<"************************* Resetting positions LEVEL "<<level<<endl;
		
		cascade(level);
		
	}
	
	
	return 0;
	
}


void all_levels::cascade(int level) {
	
	
	RANGE_loop(i, all_level_frames[level]) {
		
		frame *fr_dad = & all_level_frames[level][i];
		DI & sub_frames = fr_dad->contains_these;
		
		RANGE_loop(j, sub_frames) {
			
			frame *fr = & all_level_frames[level-1][sub_frames[j]];
			
			fr->center_x+=fr_dad->center_x;
			fr->center_y+=fr_dad->center_y;
			
			double prev_col_low= fr_dad->h_color_low;
			double prev_col_delta= fr_dad->h_color_hi - prev_col_low;
			
			fr->h_color_low = fr->h_color_low * prev_col_delta + prev_col_low;
			fr->h_color_hi = fr->h_color_hi * prev_col_delta + prev_col_low;
			
			if(level==1) {
				fr->h_color_low = fr_dad->h_color_low;
				fr->h_color_hi = fr_dad->h_color_hi;
			}
			
		}
		
	}
	
}

