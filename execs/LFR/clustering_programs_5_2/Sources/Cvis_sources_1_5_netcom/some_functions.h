



void print_circle(ostream & pout, double center_x, double center_y, double radius) {
	
	
	for(double theta=0; theta<=twopi; theta+=0.5/radius)
		pout<<center_x+radius*cos(theta)<<" "<<center_y+radius*sin(theta)<<endl;
	
}



void set_zero_baricenter(deque<frame*> & Fr) {
	
	double barix=0;
	double bariy=0;
	
	double totradii=0;
	
	for(UI i=0; i<Fr.size(); i++) {
		barix+= Fr[i]->center_x * Fr[i]->radius;
		bariy+= Fr[i]->center_y * Fr[i]->radius;
		totradii+=Fr[i]->radius;
	}
	
	
	//cout<<"barix: "<<barix<<endl;
	//cout<<"bariy: "<<barix<<endl;
	
	barix/=totradii;
	bariy/=totradii;

	for(UI i=0; i<Fr.size(); i++) {
		Fr[i]->center_x-= barix;
		Fr[i]->center_y-= bariy;
	}
	
	

}



void set_ranking_size_based(DI & ranking, DD & module_number_of_nodes, int _sign_=-1) {
	
	ranking.clear();
	multimap<double, int> A;
	
	for(UI ii=0; ii<module_number_of_nodes.size(); ii++)
		A.insert(make_pair( _sign_ * module_number_of_nodes[ii], ii));
	
	for(multimap<double, int>:: iterator itm= A.begin(); itm!= A.end(); itm++)
		ranking.push_back(itm->second);
	
	/*cout<<"set_ranking_size_based "<<endl;
	prints(module_number_of_nodes);
	prints(ranking);*/

}




void set_factors(double & node_separation, double & first_circle_radius, double & spiral_separation, double node_separation_factor, double nodes_first_circle, double spiral_separation_factor) {

	
	node_separation= node_separation_factor * paras.node_radius;
	if(spiral_separation_factor==0)
		node_separation*=1.1*paras.gap_between_frames;
	
	first_circle_radius=  nodes_first_circle * node_separation /twopi;		
	spiral_separation= spiral_separation_factor * node_separation/twopi;
	
		
}





bool compatible(double theta, double r, DD & radius_sizes, DD & theta_nodes, DD & radius_nodes) {


	/* if this function will be slow, it can be optimized, for instance returning the first false or other stuff */
	
	
	int H=theta_nodes.size()-1;
	
	double x= r * sin(theta);
	double y= r * cos(theta);
	
	
	while(H>=0) {		// here I should add that 
		
		
		
		double x1= radius_nodes[H]* sin(theta_nodes[H]);
		double y1= radius_nodes[H]* cos(theta_nodes[H]);
		//cout<<H<<" x1:  "<<x1<<" y1:  "<<y1<<endl;
		//cout<<"distance: "<<sqrt((x1-x)*(x1-x) +(y1-y)*(y1-y))<<endl;
		
		if(sqrt((x1-x)*(x1-x) +(y1-y)*(y1-y))< paras.gap_between_frames * (radius_sizes[theta_nodes.size()] + radius_sizes[H]))
			return false;
		
		--H;
	}
	
	return true;
}




double locate_nodes(DD & nodes_x, DD & nodes_y, DD & radius_sizes,
					double node_separation_factor=paras.gap_between_frames, double spiral_separation_factor=paras.Spi_Sep, double nodes_first_circle=paras.Nodes_fir_C ) {
	
	
	//int clockwise= 2*irand(1)-1;
	
	
	//cout<<"locate_nodes"<<endl;
	//prints(radius_sizes);
	
	
	int dim=radius_sizes.size();
	nodes_x.clear();
	nodes_y.clear();
	
	if(dim==0)
		return 0;
	
	if(dim==1) {
		
		nodes_x.push_back(0.);
		nodes_y.push_back(0.);
		return paras.node_radius;
		
	}
	
	// this function finds positions of nodes in group
	// nodes_x[i] gives the x coordinate of node i
	
	
	// ******************** SPIRAL PARAMETERS ************************** //
	
	double theta=0;
	double r;
	
	double node_separation, first_circle_radius, spiral_separation;
	set_factors(node_separation, first_circle_radius, spiral_separation, node_separation_factor, nodes_first_circle, spiral_separation_factor);

	
	
	// ******************** SPIRAL PARAMETERS ************************** //
	
	
	// ****************** add your sorting here ************************ //
	// ****************** add your sorting here ************************ //
	
	
	DD theta_nodes;
	DD radius_nodes;
	
	
	
	
	while(true) {
		
		
		r=spiral_separation * theta + first_circle_radius;
		
		//cout<<radius_sizes[radius_nodes.size()]<<" <--- radius"<<endl;
		
		if(compatible(theta, r, radius_sizes, theta_nodes, radius_nodes) || spiral_separation_factor<1e-7) {
			theta_nodes.push_back(theta);
			radius_nodes.push_back(r);
		}
		
		
		
		/*for(UI j=0; j<all_level_frames[i].size(); j++) {
			
			//cout<<"frame: "<<i<<" "<<all_level_frames[i][j].center_x<<" "<<all_level_frames[i][j].center_y<<endl;
			
			pout<<all_level_frames[i][j].center_x<<" "<<all_level_frames[i][j].center_y<<endl;
			print_circle(pout, all_level_frames[i][j].center_x, all_level_frames[i][j].center_y, all_level_frames[i][j].radius);
			
		}*/
		//cout<<"theta... "<<theta<<" ...r "<<r<<" spiral_separation: "<<spiral_separation<<endl;
		
		theta+=node_separation/r;
		
		if(int(theta_nodes.size())==dim)
			break;		
	}
	
	
	
	for(int i=0; i<dim; i++) {
		nodes_x.push_back(radius_nodes[i]*sin(theta_nodes[i]));
		nodes_y.push_back(radius_nodes[i]*cos(theta_nodes[i]));
	}
	
	//cout<<"loops: "<<theta/twopi<<" radius: "<<r<<endl;
	
	return r+paras.node_radius;
	
}


string RETURN_RGB(double v, double  n, double m) {

	
	char b[1000];
	//sprintf(b, "r=\"%d\" g=\"%d\" b=\"%d\"", cast_int(v*255), cast_int(n*255), cast_int(m*255));
	sprintf(b, "'%d,%d,%d'", cast_int(v*255), cast_int(n*255), cast_int(m*255));
	string s(b);
	
	//cout<<"s: "<<s<<endl;
	
	return s;
}



string HSV_to_RGB(double h_zero_one, double s, double v) {
	
	
	
	double h= h_zero_one * 6;
	

	
	// h is given on [0, 6]. s and v are given on [0, 1].
	// RGB are each returned on [0, 1].
	
	
	double m, n, f;
	int i;
	i = floor(h);
	f = h - i;
	if ( !(i&1) ) f = 1 - f; // if i is even
	m = v * (1 - s);
	n = v * (1 - s * f);
	switch (i) {
		case 0: return RETURN_RGB(v, n, m);
		case 1: return RETURN_RGB(n, v, m);
		case 2: return RETURN_RGB(m, v, n);
		case 3: return RETURN_RGB(m, n, v);
		case 4: return RETURN_RGB(n, m, v);
		case 5: return RETURN_RGB(v, m, n);
		default : return RETURN_RGB(v, v, v);
	}
	
	
}


void set_colors(deque<frame*> & Fr, double c_fact) {
	
	
	double total_number_of_nodes=0;
	for(UI i=0; i<Fr.size(); i++) {
		
		if(Fr[i]->number_of_nodes>1)
			total_number_of_nodes+=Fr[i]->radius;
		else
			++total_number_of_nodes;
	}
			
		

	double right_pointer=0;

	for(UI i=0; i<Fr.size(); i++) {
		
		
		double step=0;
		if(Fr[i]->number_of_nodes>1)
			step=Fr[i]->radius/total_number_of_nodes;
		else {
			step+=1./total_number_of_nodes;		// if there is only one node the color step is smaller
		}

		
		Fr[i]->h_color_low= right_pointer + step/2  - c_fact * step/2; 
		Fr[i]->h_color_hi= right_pointer + step/2  + c_fact * step/2;
		right_pointer+=step;
		//cout<<i<<" -> colors: "<<Fr[i]->h_color_low<<" "<<Fr[i]->h_color_hi<<" "<<right_pointer<<endl;
	
	}
	
	
	
	
	
	

}




