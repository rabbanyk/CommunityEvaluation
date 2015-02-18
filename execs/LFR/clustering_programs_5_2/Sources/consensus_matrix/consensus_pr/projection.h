





class link_collection {
	
	
public:
	
	
	link_collection(int dim_) { init_(dim_); partitions=0; };
	~link_collection(){};
	
	void add_partition(int_matrix & a);
	void init_(int  a);
	
		
	void print(ostream & po, DI & id_of, double );
	void print_disconnected_nodes(ostream & po, DI & id_of, double );
	double average_func_partitions(double a);
	
	double nmi();

	
	int_matrix memberships_per_node;
	int_matrix partition_per_node;
	int partitions;
	int_matrix modules;
	deque<int_matrix> Ps;
	
	
};



void link_collection::init_(int  a) {
	
	Ps.clear();
	modules.clear();
	memberships_per_node.clear();
	partition_per_node.clear();
	DI b;
	for(int i=0; i<a; i++) {
		memberships_per_node.push_back(b);
		partition_per_node.push_back(b);
	}
	
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
		if(h.size()>100)
			break;
	}
	
	DD coverages;
	RANGE_loop(i, Ps) {		
		coverages.push_back(double(coverage(Ps[i])) / memberships_per_node.size());
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
		
		++partitions;
		Ps.push_back(A);
		RANGE_loop(i, A) {
			
			modules.push_back(A[i]);
			RANGE_loop(j, A[i]) {
				
				memberships_per_node[A[i][j]].push_back(modules.size()-1);
				partition_per_node[A[i][j]].push_back(partitions);
				
			}
		
		}
	}
	

}


double set_inter(DI & a, DI & b) {

	
	double inter=0;
	
	
	
	
	RANGE_loop(i, a) {
		if(binary_search(b.begin(), b.end(), a[i]))
			++inter;
	}

	return inter;

}



void link_collection::print(ostream & po, DI & id_of, double thre) {
	
	
	
	RANGE_loop(i, partition_per_node) {		
		set<int> pset;
		deque_to_set(partition_per_node[i], pset);
		set_to_deque(pset, partition_per_node[i]);
	
	}
	
	
	
	set<pair<int, int> > additional_links;
	
	RANGE_loop(i, memberships_per_node) {
		
		mapii node_occurrences;
		RANGE_loop(j, memberships_per_node[i]) RANGE_loop(k, modules[memberships_per_node[i][j]]) if(i!=modules[memberships_per_node[i][j]][k]) {
			int_histogram(modules[memberships_per_node[i][j]][k], node_occurrences);
		}
	
		
		
		
		//######################################################

		int links=0;
		IT_loop(mapii, itm, node_occurrences) {
			
			
			//cout<<"---------------------------------------"<<endl;			
			double partitions_here=set_inter(partition_per_node[i], partition_per_node[itm->first]);
			//cout<<"-> "<<id_of[i]<<" "<<id_of[itm->first]<<" "<<partitions_here<<endl;			
			//cout<<"---------------------------------------"<<endl;
			
			if(itm->second>thre*partitions_here) {
				
				if (id_of[i]<id_of[itm->first])
					po<<id_of[i]<<" "<<id_of[itm->first]<<" "<<double(itm->second)/partitions_here<<endl;
				++links;
			}
		}
		
		
		//######################################################
		//cout<<"links: "<<links<<endl;
		if(links==0) {
			
			double maxw=0;
			deque<pair<int, double> > max_weighted_links;
			IT_loop(mapii, itm, node_occurrences) {
				
				double partitions_here=set_inter(partition_per_node[i], partition_per_node[itm->first]);
				double av=double(itm->second)/partitions_here;
				
				if(av>maxw+1e-6) {
					max_weighted_links.clear();
					maxw=av;
				}
				if(av>maxw-1e-6)
					max_weighted_links.push_back(make_pair(itm->first, av));
			}
			
			
			RANGE_loop(j, max_weighted_links) {
				
				pair<int, int> p;
				p.first=min(id_of[i], id_of[max_weighted_links[j].first]);
				p.second=max(id_of[i], id_of[max_weighted_links[j].first]);
				if(additional_links.find(p)==additional_links.end()) {
					po<<p.first<<" "<<p.second<<" "<<max_weighted_links[j].second<<endl;
					additional_links.insert(p);
				}
			}
		}
		
		
	}
	
	
	
}








double link_collection::average_func_partitions(double  a) {

	
	double s=a;
	s/=partitions;
	return s;

}

