

void all_levels::from_partitions_to_frames_no_overlap(int_matrix & short_pt, int_matrix & orig_pt_node) {
	
	// this function sets contains_these, number_of_nodes and radius if number_of_nodes=1
	// all_level_frames[0] is the node level
	
	deque<frame> new_frames;

	if(all_level_frames.size()==0) {	// node level
		short_pt.clear();
		for(int i=0; i<original_network.size(); i++) {		
			DI first;
			first.push_back(i);
			short_pt.push_back(first);	
		}
		orig_pt_node=short_pt;		
	}
	// ******************************************
	
	// setting frames ***************************
	RANGE_loop(i, short_pt) {
		frame F;
		F.number_of_nodes=orig_pt_node[i].size();
		if(F.number_of_nodes==1)
			F.radius=paras.node_radius;
		new_frames.push_back(F);	
	}
		
	all_level_frames.push_back(new_frames);
	
	RANGE_loop(i, short_pt) {
		all_level_frames[all_level_frames.size()-1][i].contains_these=short_pt[i];
		all_level_frames[all_level_frames.size()-1][i].original_nodes=orig_pt_node[i];
	}
	// setting frames ***************************
	
}


void all_levels::original_nodes_in_these_frames(DI & these_frames, DI & original_n, deque<frame> & last_frames) {
	
	set<int> a;
	RANGE_loop(i, these_frames) {
		deque_to_set_app(last_frames[these_frames[i]].original_nodes, a);
	}
	
	set_to_deque(a, original_n);
}


int all_levels::set_frames(deque<int_matrix> & short_parts, deque<int_matrix> & orig_pts_node_no_overlap ) {
	
	all_level_frames.clear();
	
	int_matrix void_matrix;
	int_matrix void_matrix2;
	from_partitions_to_frames_no_overlap(void_matrix, void_matrix2);
	
	for(UI ii=0; ii<short_parts.size(); ii++)
		from_partitions_to_frames_no_overlap(short_parts[ii], orig_pts_node_no_overlap[ii]);
	
	return 0;

}




