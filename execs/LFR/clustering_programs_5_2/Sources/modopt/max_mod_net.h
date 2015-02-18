
class group {
	
	
	public:
	
		group(double a, double b, set<int> & d){
			
			kin_g=a;
			ktot_g=b;
			s=d;
		
		};
		
		~group(){};
		
		double kin_g;
		double ktot_g;
		set <int> s;




};




class max_mod_net : public static_network {
	
	
	public:
				
		max_mod_net(string a) : static_network(a) {realizations=1; initial_temp=0.01; tstep=0.99;};
		~max_mod_net() {};
				
				
		double qf(const double kin_g, const double ktot_g) {
			return ((kin_g)/(2.*tstrength)- LAMBDA * pow((ktot_g)/(2.*tstrength),2));
		};
		
		void set_real(int a) {realizations=a;};
		void set_temp(double t) {initial_temp=t;};
		void set_temp_step(double t) {tstep=t;};


		
		
		deque<deque<int> > best_partition;
		double iter_simannealing(bool print_values, deque<deque<int> > & ten, string partfile);

		double LAMBDA;
		
	private:
	
		
				
		double initial_temp;
		int realizations;
		double tstep;

		
		double qtot;
		
		double check_members(deque<group> &, deque<int> &);
		double simannealing(bool , deque<deque<int> > & );
		double move(int , deque<group> & , deque<int> & , int);
		int divide(deque<group> & , deque<int> & );
		int smart_merge(deque<group> & , deque<int> & );
		
		double temp;
		deque<deque<int> > partition_work;

		
		
};

double max_mod_net::check_members(deque<group> & members, deque<int> & membership) {
	
	double qcheck=0;
	
	for (int i=0; i<members.size(); i++) {
		
		double Kin_g=kin(members[i].s);
		double Ktot_g=ktot(members[i].s);
			
		if (fabs(Kin_g - members[i].kin_g)>1e-6 || fabs(Ktot_g - members[i].ktot_g)>1e-6 ) {
			
			int err;
			cerr<<"error in check kin or ktot"<<endl;
			cin>>err;
		
		
		}
		
		for(set<int>::iterator its= members[i].s.begin(); its!=members[i].s.end(); its++)
			if(membership[*its]!=i) {
			
				int err;
				cerr<<"error in membership"<<endl;
				cin>>err;
		
		
			}
		
		
		qcheck+=qf(Kin_g, Ktot_g);	
	
	}
	
	
	if(fabs(qtot - qcheck)>1e-6 ) {
			
		int err;
		cerr<<"error in check qtot"<<endl;
		cin>>err;
		
		
	}
	
	cerr<<"passed"<<endl;
	
	return 0;


}


double max_mod_net::iter_simannealing(bool print_values, deque<deque<int> > & ten, string partfile) {
	
	// print_values is to print the values of the nodes and can be used to set_partition as starting point
	
	
	best_partition.clear();
	double qmax=-1e200;
	
	for (int i=0; i<realizations; i++) {
		
		double qs=simannealing(print_values, ten);
		if (qs>qmax or i==0) {
			best_partition=partition_work;
			qmax=qs;
        }
		cout<<"qs "<<qs<<endl;
	}
	
	
	cout<<"qmax "<<qmax<<endl;
	
	char b[10000];
	cast_string_to_char(partfile, b);
	
	ofstream partition_out(b);
	print_id(best_partition, partition_out);
	
	return qmax;

}





double max_mod_net::simannealing(bool print_values, deque<deque<int> > & ten) {

	
	
	// initialization ********************************************************************************************************
	
	qtot=0;
	
	deque<int> membership(dim);
	deque<group> members;
	
	if (!print_values) {
		for (int i=0; i<dim; i++) {
		
			set<int> first__;
			first__.insert(i);
			group g(0, vertices[i]->strength, first__);
			members.push_back(g);
			membership[i]=i;
			qtot+=qf(0, vertices[i]->strength);
			
		}
	}
	
	else {
		
		

		for (int i=0; i<ten.size(); i++) {
			
			set<int> first__;
			
			for(int j=0; j<ten[i].size(); j++) {
				first__.insert(ten[i][j]);
				membership[ten[i][j]]=i;
			}
			
			double Kin_g=kin(first__);
			double Ktot_g=ktot(first__);
			
			group g(Kin_g, Ktot_g, first__);
			members.push_back(g);
			qtot+=qf(Kin_g, Ktot_g);
			
			
		}
		
		
		//cout<<"qtot from the starting partition "<<qtot<<endl;
		
		for (int i=ten.size(); i<dim; i++) {		// empty modules
			
			set<int> first__;
			group g(0., 0., first__);
			members.push_back(g);
			
			
		}

		
		
	
	
	
	
	
	}
		
		
	

	// initialization ********************************************************************************************************
	
	
	int max_loops=50;
	
	
	
	
	temp=initial_temp;
	double qprevious=0;
	int stopper=0;
	double qmax=-1e200;

	while (true) {
		
		
		//check_members(members, membership);

		
		temp=temp*tstep;
		qprevious=qtot;
		
		for (int h=0; h<dim; h++) {
			
			int node=irand(dim-1);
			int new_com=0;
			
			new_com=membership[irand(membership.size()-1)];
			move(node, members, membership, new_com);

		}
		
		
		smart_merge(members, membership);
		divide(members, membership);
		
		cout<<"qtot: "<<qtot<<"\ttemp: "<<temp<<endl;		
		
		if (qtot>qmax) {
				
			qmax=qtot;
			partition_work.clear();
			
			for (int i=0; i<members.size(); i++) if(members[i].s.size()>0) {
					
				deque <int> v_;
				for (set<int>::iterator it=members[i].s.begin(); it!=members[i].s.end(); it++)
					v_.push_back(*it);
				partition_work.push_back(v_);
				
			}
        
		}
				
		if (fabs(qtot-qprevious)<1e-10)
			stopper++;
		else
			stopper=0;
		
		if (stopper==max_loops)
			break;


	}
	
	
	return qmax;
	
	

}







double max_mod_net::move(int node, deque<group> & members, deque<int> &membership, int new_com) {

	
	
	
	
	int old_com=membership[node];
	
	
		
	if (membership[node]==new_com)
		return 0;
	
	
	double dq=0;
	
	
	double kin_node_new=vertices[node]->kplus(members[new_com].s);
	
	
	
	{

		double kin_old=members[new_com].kin_g;
		double ktot_old=members[new_com].ktot_g;

		double qpartial1=qf(kin_old, ktot_old);

		kin_old+=2*kin_node_new;
		ktot_old+=vertices[node]->strength;

		double qpartial1_new=qf(kin_old, ktot_old);
		
		dq+=qpartial1_new - qpartial1;

	
	}
	
	
	
	
	double kin_node_old=vertices[node]->kplus(members[old_com].s);
	
		
	

	
	{


		double kin_old=members[old_com].kin_g;
		double ktot_old=members[old_com].ktot_g;

		double qpartial1=qf(kin_old, ktot_old);

		kin_old-=2*kin_node_old;
		ktot_old-=vertices[node]->strength;

		double qpartial1_new=qf(kin_old, ktot_old);
		
		dq+=qpartial1_new - qpartial1;

	
	}
	
	

	double bol=exp(dq/temp);
	
	if (ran4()<bol) {
		
		
		//cout<<dq<<" "<<qtot<<endl;
		//cout<<"move"<<endl;
		
		(members[old_com].s).erase(node);
		members[old_com].kin_g-=2*kin_node_old;
		members[old_com].ktot_g-=vertices[node]->strength;
		
		
		
		
		(members[new_com].s).insert(node);
		members[new_com].kin_g+=2*kin_node_new;
		members[new_com].ktot_g+=vertices[node]->strength;
		
		membership[node]=new_com;
		qtot+=dq;
		

	}


	return 0;


}







int max_mod_net::smart_merge(deque<group> & members, deque<int> &membership) {
	
	
	
	
	
	// to update: members, membership, qtot
	
	
	int random_g1=membership[irand(dim-1)];
	
	set<int> & merge_s1=members[random_g1].s;

	DI good_mems;
	for(set<int>::iterator its=merge_s1.begin(); its!=merge_s1.end(); its++) {
		for(int j=0; j<vertices[*its]->links->size(); j++) if(membership[vertices[*its]->links->l[j]]!=random_g1)
			good_mems.push_back(membership[vertices[*its]->links->l[j]]);
	}
	
	
	//cout<<"random_g1 "<<random_g1<<endl;
	//prints(good_mems);
	
	if(good_mems.empty())
		return 0;
	
	int random_g2=good_mems[irand(good_mems.size()-1)];	
	if (random_g1==random_g2 || members[random_g1].s.empty() || members[random_g2].s.empty())
		return 0;
	
	
	
	double qpartial1=qf(members[random_g1].kin_g, members[random_g1].ktot_g);
	double qpartial2=qf(members[random_g2].kin_g, members[random_g2].ktot_g);
	
	
	set<int> merge_s=members[random_g1].s;
	
	for (set<int>::iterator its= members[random_g2].s.begin(); its!=members[random_g2].s.end(); its++)
		merge_s.insert(*its);
	
	
	double Kin_merge= kin(merge_s);
	double Ktot_merge= ktot(merge_s);
	
	
	double qpartial3=qf(Kin_merge, Ktot_merge);
	
	double dq= - qpartial1 - qpartial2 + qpartial3;
	
	
	double bol=exp(dq/temp);
	
	if (ran4()<bol) {
		
		
		for (set<int>::iterator its= members[random_g1].s.begin(); its!=members[random_g1].s.end(); its++)
			membership[*its]=random_g2;
		
		
		members[random_g1].s.clear();
		members[random_g1].kin_g=0;
		members[random_g1].ktot_g=0;
		
		
		members[random_g2].s=merge_s;
		members[random_g2].kin_g=Kin_merge;
		members[random_g2].ktot_g=Ktot_merge;
		
		
		qtot+=dq;
		
		
	}
	
	return 0;
	
	
}


int max_mod_net::divide(deque<group> & members, deque<int> & membership) {
	
	

	//------------------------------------------------------------------------------------------ old 
	
	
	// to update: members, membership, qtot
	
	
	int random_g=membership[irand(dim-1)];
	
	if (members[random_g].s.size()<2)
		return 0;
	
	
	double kin_old=members[random_g].kin_g;
	double ktot_old=members[random_g].ktot_g;

	double qtot_p=qtot;
	double qpartial1=qf(kin_old, ktot_old);

	
	set <int> s_old_=members[random_g].s;
	
	
	
	
	
	//------------------------------------------------------------------------------------------ old 
	
	
	
	int empty_g;
	for (int mi=0; mi<members.size(); mi++) 
		if(members[mi].s.size()==0) {
			empty_g=mi;
			break;
		}
	
	
	
	for (set<int>::iterator it=s_old_.begin(); it!=s_old_.end(); it++) {
		
		int node=*it;
		
		
		if (ran4()<0.5) {
		
			members[empty_g].s.insert(node);
			members[random_g].s.erase(node);
			
			membership[node]=empty_g;
			
		
		}
		
	
	}
	
	

	
	
	members[random_g].kin_g=kin(members[random_g].s);
	members[random_g].ktot_g=ktot(members[random_g].s);
	
	members[empty_g].kin_g=kin(members[empty_g].s);
	members[empty_g].ktot_g=ktot(members[empty_g].s);
	
	
	double qpnew=qf(members[random_g].kin_g, members[random_g].ktot_g) + qf(members[empty_g].kin_g, members[empty_g].ktot_g);
	
	qtot+=qpnew-qpartial1;
	
	for(set<int>::iterator it=s_old_.begin(); it!=s_old_.end(); it++)
		if (membership[*it]==random_g)
			move(*it, members, membership, empty_g);
		else
			move(*it, members, membership, random_g);
	
	
	
	double dq=qtot - qtot_p;
	
	double bol=exp(dq/temp);
	
	if (ran4()>bol) {
		
		
		
		for (set<int>::iterator it=members[empty_g].s.begin(); it!=members[empty_g].s.end(); it++)
			membership[*it]=random_g;
			
		
		members[random_g].kin_g=kin_old;
		members[random_g].ktot_g=ktot_old;
		
		members[empty_g].kin_g=0;
		members[empty_g].ktot_g=0;
		
		members[empty_g].s.clear();
		members[random_g].s=s_old_;
		
		qtot=qtot_p;
		
	}
	
	
	return 0;
		
		
}





