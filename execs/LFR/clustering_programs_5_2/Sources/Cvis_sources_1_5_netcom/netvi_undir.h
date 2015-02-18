



#if !defined(frames_included)
#define frames_included


#include "netvi.h"
#include "all_levels.h"




int netvi::print_gdf(deque<frame> & frames, ostream & pout, ostream & pout2, double radius_factor, netvi & original_network, int_matrix & level__number_of_mems) {
	
	
	DD nodes_x; 
	DD nodes_y;
	
	for(UI i=0; i<frames.size(); i++) {
		nodes_x.push_back(frames[i].center_x);
		nodes_y.push_back(frames[i].center_y);
	}
	
	pout<<"nodedef>name VARCHAR,label VARCHAR, visible BOOLEAN,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR"<<endl;
	
	
	//cout<<"dims: "<<dim<<" "<<original_size<<endl;
	/*cout<<"-> level__number_of_mems[i].size()"<<level__number_of_mems[0].size()<<endl;
	 printm(level__number_of_mems);*/
	
	if(dim==original_network.size()) {
		
		
		char file_pat[1000];
		sprintf(file_pat, "%s/node_memberships.dat", folder_output);
		ofstream pat(file_pat);
		
		//********************************************* network level **************************************************************
		
		for(UI i=0; i<frames.size(); i++) {	
			
			
			pat<<"node: "<<id_of(i)<<" memberships:\t";
			prints(level__number_of_mems[i], pat);
			
			double saturat, brightn;
			choose_saturation_and_brightness(level__number_of_mems[i], saturat, brightn);
			if(level__number_of_mems[i][level__number_of_mems[i].size()-1]==0)			// if the node is homeless even here it is white
				brightn=1;
			string color= HSV_to_RGB((frames[i].h_color_low + frames[i].h_color_hi)/2, saturat, brightn);
			
			
			double width= radius_factor* paras.node_radius;
			
			string new_n;
			original_network.print_most_important_among(frames[i].original_nodes, new_n);
			
			//cout<<"label: "<<id_word[id_of(i)]<<" memeberships: "<<mall.memberships[i].size()<<endl;
			pout<<id_of(i)<<",'"<<new_n<<"',"<<"true"<<","<<"true"<<","<<width<<","<<width<<","<<frames[i].center_x<<","<<frames[i].center_y<<","<<color<<endl;
		}
		
		
	}
	
	//********************************************* higher levels **************************************************************
	
	
	else for(UI ii=0; ii<frames.size(); ii++) {
		
		double unique_value=0.8;
		string color= HSV_to_RGB((frames[ii].h_color_low + frames[ii].h_color_hi)/2, unique_value, unique_value);
		
		if(frames[ii].number_of_nodes==1) {
			
			int nodei=frames[ii].original_nodes[0];
			
			//cout<<"node: "<<nodei<<endl;
			
			double saturat, brightn;
			choose_saturation_and_brightness(level__number_of_mems[nodei], saturat, brightn);
			
			if(level__number_of_mems[nodei][level__number_of_mems[nodei].size()-1]==0)			// if the node is homeless even here it is white
				brightn=1;
			
			color= HSV_to_RGB((frames[ii].h_color_low + frames[ii].h_color_hi)/2, saturat, brightn);
		}
		
		frames[ii].color_=color;
		//double size_node=paras.node_radius*(1+log(frames[ii].radius*radius_factor));
		double size_node= frames[ii].radius*radius_factor ;
		
		string new_n;
		original_network.print_most_important_among(frames[ii].original_nodes, new_n);
		
		pout<<ii<<",'"<<new_n<<"',"<<"true"<<","<<"true"<<","<<size_node<<","<<size_node<<","<<frames[ii].center_x<<","<<frames[ii].center_y<<","<<color<<endl;

		//pout2<<"#module: "<<ii<<" nodes: "<<frames[ii].number_of_nodes<<endl;
		//prints(frames[ii].contains_these, pout2);
		
	}
	
	cout<<"recording links... "<<endl;
	pout<<"edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE"<<endl;
	
	if(dim==original_network.size()) {
		if(paras.weighted)	{
			for (int i=0; i<int(vertices.size()); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(i < vertices[i]->links->l[j])
					pout<<id_of(i)<<","<<id_of(vertices[i]->links->l[j])<<","<<vertices[i]->links->w[j].second<<endl;
		} else {
			for (int i=0; i<int(vertices.size()); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(i < vertices[i]->links->l[j])
					pout<<id_of(i)<<","<<id_of(vertices[i]->links->l[j])<<","<<vertices[i]->links->w[j].first<<endl;
		}
	}
	else {
		
		
		if(paras.weighted)	{
			for (int i=0; i<int(vertices.size()); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(i < vertices[i]->links->l[j])
					pout<<i<<","<<vertices[i]->links->l[j]<<","<<vertices[i]->links->w[j].second<<endl;
		} else {
			for (int i=0; i<int(vertices.size()); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(i < vertices[i]->links->l[j])
					pout<<i<<","<<vertices[i]->links->l[j]<<","<<vertices[i]->links->w[j].first<<endl;
		}
	}
	return 0;
	
}




double netvi::square_distance_node(int i, deque<frame*> & Fr) {
	
	double tlength=0;
	
	for(int j=0; j<vertices[i]->links->size(); j++) {
		
		
		double dx=(Fr[i]->center_x - Fr[vertices[i]->links->l[j]]->center_x);
		double dy=(Fr[i]->center_y - Fr[vertices[i]->links->l[j]]->center_y);
		if(paras.weighted) {
			dx*=vertices[i]->links->w[j].second;
			dy*=vertices[i]->links->w[j].second;
		} else {
			dx*=vertices[i]->links->w[j].first;
			dy*=vertices[i]->links->w[j].first;
		}
		
		
		tlength+= dx*dx+dy*dy;
	}
	
	return tlength;
}




int netvi::check_node_assignments(int_matrix & tool_matrix, set<int> & homeless_nodes) {
	
	// this function checks nodes which have all their links towards a certain group
	// if this group is not their own group they are assigned to that
	
	// in tool_matrix there is no overlap
	
	
	DI node_community;
	for(int i=0; i<dim; i++)
		node_community.push_back(0);
	
	RANGE_loop(i, tool_matrix) RANGE_loop(j, tool_matrix[i]) {
		node_community[tool_matrix[i][j]]=i;
	}
	
	
	
	RANGE_loop(i, tool_matrix) RANGE_loop(j, tool_matrix[i]) {
		
		set<int> neigh_mems;
		
		bool homeless_to_homeless=false;
		if(vertices[tool_matrix[i][j]]->links->size()==1 && homeless_nodes.find(vertices[tool_matrix[i][j]]->links->l[0])!=homeless_nodes.end() )
			homeless_to_homeless=true;
		
		for(int h=0; h<vertices[tool_matrix[i][j]]->links->size(); h++)
			neigh_mems.insert(node_community[vertices[tool_matrix[i][j]]->links->l[h]]);
		
		
		if(homeless_to_homeless==false) if(neigh_mems.size()==1 && *(neigh_mems.begin())!= node_community[tool_matrix[i][j]])
			node_community[tool_matrix[i][j]]=*(neigh_mems.begin());		
	}
	
	
	set_partition_from_memberships(node_community, tool_matrix);	
	
	//cout<<"check_node_assignments tool matrix after"<<endl;
	//printm(tool_matrix);
	
	
	return 0;
	
}



int netvi::check_node_assignments_power(int_matrix & tool_matrix, set<int> & homeless_nodes) {
	
	// this function checks nodes which have all their links towards a certain group
	// if this group is not their own group they are assigned to that
	
	// in tool_matrix there is no overlap
	
	//cout<<"check_node_assignments tool matrix POWER"<<endl;
	//printm(tool_matrix);
	
	
	DI node_community;
	for(int i=0; i<dim; i++)
		node_community.push_back(0);
	
	RANGE_loop(i, tool_matrix) RANGE_loop(j, tool_matrix[i]) {
		node_community[tool_matrix[i][j]]=i;
	}
	
	
	int change=1;
	
	set<int> check_this=homeless_nodes;
	
	while(change>0) {
		
		
		change=0;
		
		//cout<<"homeless "<<endl;
		//prints(homeless_nodes);
		
		for(set<int>:: iterator it=check_this.begin(); it!=check_this.end(); it++) {
			
			
			
			int nodei=*it;
			//cout<<"-------------- nodei: "<<nodei<<endl;
			set<int> neigh_mems;
			
			for(int h=0; h<vertices[nodei]->links->size(); h++) if(homeless_nodes.find(vertices[nodei]->links->l[h])==homeless_nodes.end() )
				neigh_mems.insert(node_community[vertices[nodei]->links->l[h]]);
			
			
			/*cout<<"neigh_mems"<<endl;
			 prints(neigh_mems);*/
			
			if(neigh_mems.size()==1 && *(neigh_mems.begin())!= node_community[nodei]) {
				node_community[nodei]= *(neigh_mems.begin());
				++change;
				homeless_nodes.erase(nodei);				
				
			}			
		}
		
		//cout<<"last level: "<<change<<endl;
		check_this=homeless_nodes;
		
	}
		
	set_partition_from_memberships(node_community, tool_matrix);	
	
	/*cout<<"check_node_assignments tool matrix after power"<<endl;
	 printm(tool_matrix);*/
	
	return 0;
	
}


#endif



