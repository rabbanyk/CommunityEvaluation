
#include "standard_package/standard_include.cpp"


class convert_tree_format  {
	
	
public:
	
	convert_tree_format(){};
	~convert_tree_format(){};
	
	deque<int_matrix> original_pts;
	deque<int_matrix> short_pts;
	
	int set_orig_and_short_pts(string filename);
	void print(string file);
	
	
	
private:
	
	bool separate_strings_special(string &b, deque<string> & v);	
	UI set_node_mems(string file);	
	int set_partition_from_node_mems(int_matrix & partition_this_level, UI level);
	void translate(deque<int_matrix> & A, int_matrix & b);
	void set_short_pts(int_matrix & tool_matrix, int_matrix & 	lower_matrix, int_matrix & upper_matrix);
	
	void compute_simple_mems(int_matrix & memberships_with_previous, DI & memberships);
	
	map<int, DI> node_mems;
	
};


bool convert_tree_format::separate_strings_special(string &b, deque<string> & v) {		
	
	
	
	v.clear();
	string s1;
	
	for(int i=0; i<int(b.size()); i++) {
		
		
		if(b[i]==' ' || b[i]=='\t' || b[i]=='\n' || b[i]==',' || b[i]=='"') {
			
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

void convert_tree_format::translate(deque<int_matrix> & A, int_matrix & b) {
	
	DI node_id;
	for(map<int, DI> :: iterator itm=node_mems.begin(); itm!=node_mems.end(); itm++)
		node_id.push_back(itm->first);
	
	RANGE_loop(ii, A) RANGE_loop(i, A[ii]) RANGE_loop(j, A[ii][i]) {
		A[ii][i][j]=node_id[A[ii][i][j]];
	}
	
	RANGE_loop(i, b) RANGE_loop(j, b[i]) {
		b[i][j]=node_id[b[i][j]];
	}
	
	
	
}




void convert_tree_format::compute_simple_mems(int_matrix & memberships_with_previous, DI & memberships) {
	
	
	memberships.clear();
	map<DI, int> prev_simple;
	RANGE_loop(i, memberships_with_previous) {
		
		prev_simple.insert(make_pair(memberships_with_previous[i], prev_simple.size()));	
		memberships.push_back(prev_simple[memberships_with_previous[i]]);
		
	}
	
}




int convert_tree_format::set_partition_from_node_mems(int_matrix & partition_this_level, UI level) {
	
	
	partition_this_level.clear();
	
	
	
	
	//cout<<" ...LEVEL::: "<<level<<endl;
	
	int_matrix memberships_with_previous;
	
	for(map<int, DI> :: iterator itm=node_mems.begin(); itm!=node_mems.end(); itm++) {
		
		DI & mems=itm->second;
		if(level>=mems.size())
			mems.push_back(1);
		
		/*cout<<"node: "<<itm->first<<endl;
		prints(mems);*/
		
		
		DI allmems;
		for(UI i=0; i<=level; i++)
			allmems.push_back(mems[i]);
		
		memberships_with_previous.push_back(allmems);
		
	}
	
	
	DI memberships;
	compute_simple_mems(memberships_with_previous, memberships);
	
	/*cout<<"memberships_with_previous"<<endl;
	 printm(memberships_with_previous);
	 
	 cout<<"memberships"<<endl;
	 prints(memberships);*/
	
	set_partition_from_memberships(memberships, partition_this_level);
	
	
	return 0;
	
}

UI convert_tree_format::set_node_mems(string file) {
	
	
	
	//node_mems[i] are simply the number correspondent to node i
	
	string s;
	char b[1000];
	cast_string_to_char(file, b);
	ifstream gin(b);
	
	
	UI max_depth=0;
	while(getline(gin, s)) if(s.size()>0 && s[0]!='#') {
		
		//cout<<"s-> "<<s<<endl;
		
		deque<string> seps;
		separate_strings_special(s, seps);		
		//prints(seps);
		
		int node=cast_int(cast_string_to_double(seps[seps.size()-1]));
		//cout<<"node: "<<node<<endl;
		
		DI mems;
		cast_string_to_doubles(seps[0], mems);
		
		if(mems.size()>0)
			mems.pop_back();
		if(mems.size()>max_depth)
			max_depth=mems.size();
		
		//prints(mems);
		node_mems[node]=mems;
		
		
	}
	
	//cout<<"max_depth: "<<max_depth<<endl;
	
	return max_depth;
	
}



void convert_tree_format::set_short_pts(int_matrix & tool_matrix, int_matrix & 	lower_matrix, int_matrix & upper_matrix) {
	
	
	DI memberships_lower_level;
	for(UI i=0; i<node_mems.size(); i++)
		memberships_lower_level.push_back(0);
	
	RANGE_loop(i, lower_matrix) RANGE_loop(j, lower_matrix[i]) {
		memberships_lower_level[lower_matrix[i][j]]=i;
	}
	
	
	tool_matrix.clear();
	RANGE_loop(i, upper_matrix) {
		
		set<int> c;
		RANGE_loop(j, upper_matrix[i]) {
			int mem_node=memberships_lower_level[upper_matrix[i][j]];
			c.insert(mem_node);
		}
		
		DI f;
		set_to_deque(c, f);
		tool_matrix.push_back(f);
	}
	
}






int convert_tree_format::set_orig_and_short_pts(string filename) {
	
	
	
	UI max_depth= set_node_mems(filename);
	
	
	//cout<<"max_depth:************************************ "<<max_depth<<endl;
	
	if(max_depth==0)
		max_depth=1;
	
	for(UI i=0; i<max_depth; i++) {
		
		int_matrix partition_this_level;
		set_partition_from_node_mems(partition_this_level, i);
		original_pts.push_back(partition_this_level);		
	}
	
	
	// if level+1==original_pts.size() -> level is the last level
	
	RANGE_loop(level, original_pts)  {
		
		
		if(level+1!=original_pts.size()) {
			int_matrix & lower_matrix= original_pts[level+1];
			int_matrix & upper_matrix= original_pts[level];	
			int_matrix tool_matrix;
			set_short_pts(tool_matrix, lower_matrix, upper_matrix);
			short_pts.push_back(tool_matrix);
			
		}
		else
			short_pts.push_back(original_pts[level]);
		
	}
	
	
	translate(original_pts, short_pts[short_pts.size()-1]);
	
	
	return 0;
}


void print_shortly(int_matrix & A, ostream & po) {
	
	
	RANGE_loop(i, A) {
		
		RANGE_loop(j, A[i]) {
			po<<A[i][j]<<" ";
		}
		
		po<<endl;
	}
	
	
	
}

void convert_tree_format::print(string file) {
	
	string s;
	char directory_char[1000];
	cast_string_to_char(file, directory_char);
	
	char char_to_use[1000];
	sprintf(char_to_use, "mkdir %s", directory_char);
	int sy=system(char_to_use);
	sprintf(char_to_use, "rm -r %s/*", directory_char);
	sy=system(char_to_use);	
	
	RANGE_loop(i, original_pts) {
		sprintf(char_to_use, "%s/original_pt_%d.dat", directory_char, int(original_pts.size()) -1 -int(i) );
		ofstream go(char_to_use);
		print_shortly(original_pts[i], go);
	}
	
	RANGE_loop(i, short_pts) {
		
		if(i+1==short_pts.size()) {
			sprintf(char_to_use, "%s/tp", directory_char);			
		} else {
			sprintf(char_to_use, "%s/short_tp%d", directory_char, int(original_pts.size()) -1 -int(i) );
		}
		ofstream go(char_to_use);
		print_shortly(short_pts[i], go);
	}
	
}







