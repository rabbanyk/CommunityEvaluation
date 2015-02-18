




/*
 things I might want to improve:
 1. the swap could be improved quite a lot, especially swap_homeless_frames (add a self-avoiding constraint?)
 2. some_functions.h -> locate_nodes is slow. it could be made more efficient
 */


bool separate_strings_cvis(string &b, deque<string> & v) {		
	
	
	
	v.clear();
	string s1;
	
	for(int i=0; i<int(b.size()); i++) {
		
		
		if(b[i]==' ' || b[i]=='\t' || b[i]=='\n' || b[i]==',' || b[i]=='\"' || b[i]=='\'') {
			
			if(s1.size()>0)
				v.push_back(s1);
			
			s1.clear();
			
			
		} else
			s1.push_back(b[i]);
		
		
		if(i==int(b.size())-1) {
			
			
			if(s1.size()>0)
				v.push_back(s1);			
			s1.clear();
			
			
		}
		
		
		
	}
	
	
	
	return true;
	
	
	
	
}




int main(int argc, char * argv[]) {		
	
	
	
	if(argc<2) {
		program_statement(argv[0]);
		return -1;
	}
	
	if(paras._set_(argc, argv)==false)
		return -1;
	
	paras.print();
	
	string netfile=paras.file1;
	
	
	/* check if file_name exists */
	if(check_if_file_exists(netfile)==false)
		return -1;
	
	all_levels AL;
	AL.original_network.set_graph(netfile);
	
	if(AL.original_network.size()==0 || AL.original_network.stubs()==0) {
		cerr<<"network empty"<<endl;
		return -1;
	}
	
	
	if(paras.labels_set) {
		
		char labels_file_c[1000];
		cast_string_to_char(paras.label_file, labels_file_c);
		ifstream gin(labels_file_c);
		
		if(gin.is_open()==false) {
			cout<<"label file not found"<<endl;		
		}
		
		string ww;
		while(getline(gin, ww)) {
			
			
			deque<string> ws;
			separate_strings_cvis(ww, ws);
			string news;
			
			if(ws.size()<2) 
				cout<<"label file has something wrong. (not specified?)"<<endl;
			else {				
				for(int i=0; i<int(ws.size())-1; i++)
					news+=ws[i];
				int idd=cast_int(cast_string_to_double(ws[ws.size()-1]));
				AL.original_network.id_word[idd]=news;
				
			}
		}
		
	}
	
	
	cout<<"network:: "<<AL.original_network.size()<<" nodes and "<<AL.original_network.stubs()<<" stubs;\t average degree = "<<AL.original_network.stubs()/AL.original_network.size()<<endl;
	
	if(paras.tp_flag && paras.partitions_file.size()>0 && paras.partitions_file[paras.partitions_file.size()-1]=='/') {
		
		string news;
		for(int i=0; i<int(paras.partitions_file.size())-1; i++)
			news.push_back(paras.partitions_file[i]);
		
		paras.partitions_file=news;
	
	}
		
		
		
	char directory_char[1000];
	cast_string_to_char(netfile + "_cvis", directory_char);
	
	folder_output=directory_char;
	
	not_overlapping_tree NTree;
	if(paras.tp_flag) {
		
		
		if(NTree.read(AL.original_network, paras.partitions_file)==-1)
			return -1;
		
		{
			char char_to_use[1000];
			sprintf(char_to_use, "mkdir %s", directory_char);
			int sy=system(char_to_use);
			sprintf(char_to_use, "rm -r %s/*", directory_char);
			sy=system(char_to_use);		
		}
		
		
	}
	
	
	
	if(paras.tree_flag) {
		
		
		
		/* check if file_name exists */
		if(check_if_file_exists(paras.partitions_file)==false)
			return -1;
		
		
		
		convert_tree_format Martin_tree;
		Martin_tree.set_orig_and_short_pts(paras.partitions_file);
		
		string folder_output_string(folder_output);
		
		Martin_tree.print(folder_output_string);
		if(NTree.read(AL.original_network, folder_output_string)==-1)
			return -1;
	
	
	}
	
	
	
	
	AL.node_mems_across_levels=NTree.level__number_of_mems;
	
	
	
		
	/*RANGE_loop(i, NTree.short_pts) printm(NTree.short_pts[i]);
	RANGE_loop(i, NTree.orig_pts_node_no_overlap) printm(NTree.orig_pts_node_no_overlap[i]);
	printm(NTree.level__number_of_mems);*/
	
	AL.set_frames(NTree.short_pts, NTree.orig_pts_node_no_overlap);
	
	//AL.print_all_frames_inclusions();
	
	
	AL.set_radii_and_positions();
	
	AL.reset_frames();
	cout<<"recording... "<<endl;
	
	//AL.print_all_frames_radii();
	AL.print_levels_gdf();
	//AL.print_levels_multiple_files();
	
	return 0;
		
	
	

}



