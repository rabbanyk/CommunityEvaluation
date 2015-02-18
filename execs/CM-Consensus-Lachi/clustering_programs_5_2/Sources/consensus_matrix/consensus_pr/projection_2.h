





class link_collection {
	
	
public:
	
	
	link_collection(int dim_) { init_(dim_); partitions=0; };
	~link_collection(){};
	
	void add_partition(int_matrix & a);
	void init_(int  a);
	
	int projection(DI & a);
	void add_link(int a, int b, double w);
	
	void print(ostream & po, DI & id_of, double );
	void print_disconnected_nodes(ostream & po, DI & id_of, double );
	double average_func_partitions(double a);
	//double std_func_partitions(DD & a);
	
	//bool check_errors(double rew_probability, double & max_err__);
	double nmi();

	
	deque<map<int, double> > link_per_node;
	
	
	int partitions;
	deque<int_matrix> Ps;
	
	
};



void link_collection::init_(int  a) {
	
	
	link_per_node.clear();
	map<int, double> b;
	for(int i=0; i<a; i++)
		link_per_node.push_back(b);
	
}


int coverage(int_matrix & A) {

	set<int> s;

	
	RANGE_loop(i, A) {
		deque_to_set_app(A[i], s);
	}

	return int(s.size());

}


double link_collection::nmi() {

	if(Ps.size()<2)
		return 1;
	
	
	
	DD h;
	while(true) {
		
		int one=irand(Ps.size()-1);
		int two=irand(Ps.size()-1);
		
		if(one!=two) {
			
			double nmi=mutual3(Ps[one], Ps[two]);
			if(nmi!=nmi)
				nmi=0;
			h.push_back(nmi);
		}
		if(h.size()>1000)
			break;
	}
	
	DD coverages;
	RANGE_loop(i, Ps) {		
		coverages.push_back(double(coverage(Ps[i])) / link_per_node.size());
	}
	
	//double covo=average_func(coverages);
	
	//av_nmi=average_func(h);
	//cout<<"average_ nmi: "<<av_nmi<<" "<<covo<<endl;
	//av_nmi=av_nmi*covo;
	
	
	//prints(h);
	
	return average_func(h);
	

}

void link_collection::add_partition(int_matrix & A) {
	
	
	if(A.size()>1) {
		
		RANGE_loop(i, A) {
			projection(A[i]);
		}
	
		++partitions;
		Ps.push_back(A);
	}
	

}




int link_collection::projection(DI & A) {
	
	
	if(A.size()>0.9*link_per_node.size())
		return -1;
	
	RANGE_loop(i, A) for(UI j=i+1; j<A.size(); j++) {
		//add_link(A[i], A[j], 1./(A.size()-1));
		//add_link(A[j], A[i], 1./(A.size()-1));
		//cout<<"--> "<<A[i]<<" "<<A[j]<<endl;
		add_link(A[i], A[j], 1.);
		add_link(A[j], A[i], 1.);
	}
	
	return 0;
	
}



void link_collection::add_link(int a, int b, double w) {
		
	// here we go
	if(link_per_node[a].find(b)==link_per_node[a].end()) {
		double f=0;
		link_per_node[a].insert(make_pair(b, f));
	}
	link_per_node[a][b]+=w;
}


void link_collection::print(ostream & po, DI & id_of, double thre) {
	
	RANGE_loop(i, link_per_node) {
		
		for(map<int, double>::iterator itm= link_per_node[i].begin(); itm!=link_per_node[i].end(); itm++) {
			if(itm->first>int(i)) {
				
				double av=average_func_partitions(itm->second);
				if(av>thre)
					po<<id_of[i]<<" "<<id_of[itm->first]<<" "<<av<<endl;				
			}
		}		
	}
}




void link_collection::print_disconnected_nodes(ostream & po, DI & id_of, double thre) {
	
	// this function is to print nodes which would have degree =0 because of the threshold
	// to fix that, the idea is simply to look at them and print the best weighted links they have 
	
	DI disconnected_nodes;
	
	RANGE_loop(i, link_per_node) {
		bool dis=true;
		for(map<int, double>::iterator itm= link_per_node[i].begin(); itm!=link_per_node[i].end(); itm++) {
			double av=average_func_partitions(itm->second);
			if(av>thre)
				dis=false;				
		}		
		if(dis==true)
			disconnected_nodes.push_back(i);
	}
	
	RANGE_loop(i, disconnected_nodes) {
		
		double maxw=0;
		deque<pair<int, double> > max_weighted_links;
		for(map<int, double>::iterator itm= link_per_node[disconnected_nodes[i]].begin(); itm!=link_per_node[disconnected_nodes[i]].end(); itm++) {
			
			double av=average_func_partitions(itm->second);
			
			if(av>maxw+1e-6) {
				max_weighted_links.clear();
				maxw=av;
			}
			if(av>maxw-1e-6)
				max_weighted_links.push_back(make_pair(itm->first, av));
		}
		
		RANGE_loop(j, max_weighted_links) {
			po<<id_of[disconnected_nodes[i]]<<" "<<id_of[max_weighted_links[j].first]<<" "<<max_weighted_links[j].second<<endl;			
		}
		
	
	}
	
	
	
}




double link_collection::average_func_partitions(double  a) {

	
	double s=a;
	s/=partitions;
	return s;

}


/*
double link_collection::std_func_partitions(DD & a) {
	
	
	double av=0;
	double std=0;
	RANGE_loop(i, a) {
		av+=a[i];
		std+=a[i]*a[i];
	}
	
	av/=partitions;
	std/=partitions;
	
	std= sqrt((std - av*av)/partitions);
	return std;
	
}*/






/*
 bool link_collection::check_errors(double rew_probability, double & max_err__) {
 
 // so, I want to say that the errors on the averages should be smaller than the rew_probability_precision * rew_probability.
 
 
 DD errors;
 RANGE_loop(i, link_per_node) {
 
 for(map<int, DD>::iterator itm= link_per_node[i].begin(); itm!=link_per_node[i].end(); itm++) {
 errors.push_back(std_func_partitions(itm->second));
 }
 
 }
 
 
 //cout<<"errors"<<endl;
 //prints(errors);
 
 double max_err=0;
 RANGE_loop(i, errors) max_err=max(max_err, errors[i]); 
 
 cout<<"max_err: "<<max_err<<" should be less than "<<rew_probability_precision * rew_probability<<endl;
 max_err__=max_err;
 
 if(max_err<rew_probability_precision * rew_probability)
 return true;
 
 return false;
 
 
 }
 */

