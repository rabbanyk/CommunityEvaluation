





class link_collection {
	
	
public:
	
	
	link_collection(int dim_) { init_(dim_); partitions=0; };
	~link_collection(){};
	
	void add_partition(int_matrix & a);
	void init_(int  a);
	
	int projection(DI & a);
	void add_link(int a, int b, double w);
	
	void print(ostream & po, DI & id_of);
	
	double average_func_partitions(DD & a);
	double std_func_partitions(DD & a);
	
	//bool check_errors(double rew_probability, double & max_err__);
	double nmi();

	
	deque<map<int, DD> > link_per_node;
	// this could be actually better designed because I'm just interested in averages and standard deviations
	
	
	int partitions;
	deque<int_matrix> Ps;
	
	
};



void link_collection::init_(int  a) {
	
	
	link_per_node.clear();
	map<int, DD> b;
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
	
	RANGE_loop(i, A) {
		projection(A[i]);
	}
	
	++partitions;
	Ps.push_back(A);
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
		DD f;
		link_per_node[a].insert(make_pair(b, f));
	}
	
	link_per_node[a][b].push_back(w);
	
	
}


void link_collection::print(ostream & po, DI & id_of) {
	
	RANGE_loop(i, link_per_node) {
		
		for(map<int, DD>::iterator itm= link_per_node[i].begin(); itm!=link_per_node[i].end(); itm++) {
			if(itm->first>i) {
				po<<id_of[i]<<" "<<id_of[itm->first]<<" "<<average_func_partitions(itm->second)<<endl;				
			}
		}		
	}
}




double link_collection::average_func_partitions(DD & a) {

	
	double s=0;
	RANGE_loop(i, a) {
		s+=a[i];
	}

	s/=partitions;
	return s;

}



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
	
}






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

