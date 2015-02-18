
#include "./standard_package/standard_include.cpp"




class module_collection {		

	/* all the labels refers to the index in int_matrix modules */

	
	public:
	
		module_collection(int d);
		~module_collection(){};

		int size() { return module_bs.size(); };		
		
		bool insert(deque<int> & c, double bs, int & new_name);
		bool insert(deque<int> & c, double bs);
		bool erase(int);			
		void print(ostream & outt, deque<int> & netlabels, bool);		
		
		void fill_gaps();
		void put_gaps();
		void homeless(deque<int> & h);
		int coverage();
		int effective_groups();

		void set_partition(deque<deque<int> > & A);
		void set_partition(deque<deque<int> > & A, deque<double> & b);
		
	
		/*************************** DATA ***************************/
		
		
		deque<set<int> > memberships;					
		int_matrix modules;		
		map<int, double> module_bs;						/* it maps the module id into the b-score */
		
		/***********************************************************/

	private:
		
		void _set_(int dim);
		bool check_already(const deque<int> & c);

};




module_collection::module_collection(int dim) {

	_set_(dim);
}


void module_collection::_set_(int dim) {

	
	set<int> first;
	for(int i=0; i<dim; i++)
		memberships.push_back(first);
	
		


}



bool module_collection::insert(deque<int> & c, double bs) {


	int new_name;
	return insert(c, bs, new_name);
	
}




bool module_collection::insert(deque<int> & c, double bs, int & new_name) {

	
	if(bs==0)
		bs=ran4() * 1e-100;
	
	
	sort(c.begin(), c.end());
	new_name=-1;
	
	if(check_already(c)==true) {	
		
		new_name=modules.size();
		for(int i=0; i<int(c.size()); i++)
			memberships[c[i]].insert(new_name);	
		
		
		modules.push_back(c);
		module_bs[new_name]=bs;		
		return true;
		
	}
	
	return false;
	
}





bool module_collection::erase(int a) {
	
	// it erases module a 
	
	
	
	if(module_bs.find(a)==module_bs.end())		// it only erases not empty modules
		return false;
	
	deque<int> & nodes_a = modules[a];
	
	for(int i=0; i<int(nodes_a.size()); i++)
		memberships[nodes_a[i]].erase(a); 
	
			
	modules[a].clear();
	module_bs.erase(a);
	
	
	
	return true;

}











void module_collection::print(ostream & outt, deque<int> & netlabels, bool not_homeless) {
	
	
		
	
	int nmod=0;
	for(map<int, double >::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(not_homeless==false || modules[itm->first].size() > 1) {
		
		
		nmod++;
		
		
		deque<int> & module_nodes= modules[itm->first];
		outt<<"#module "<<itm->first<<" size: "<<modules[itm->first].size()<<" bs: "<<module_bs[itm->first]<<endl;
		
		deque<int> labseq;
		for(int i=0; i<int(module_nodes.size()); i++) {
			labseq.push_back(netlabels[module_nodes[i]]);
		}
		
		sort(labseq.begin(), labseq.end());
		
		for(int i=0; i<int(labseq.size()); i++) {
			outt<<labseq[i]<<" ";
		}
		outt<<endl;
		
		

	}
	
	

}





void module_collection::fill_gaps() {

	for(int i=0; i<int(memberships.size()); i++)
		if(memberships[i].size()==0) {
			
			deque<int> new_d;
			new_d.push_back(i);
			insert(new_d, 1.);
			
		
		}

}



void module_collection::put_gaps() {
		
	deque<int> to_erase;
	
	
	for(int i=0; i<int(modules.size()); i++) {
		
		if(modules[i].size()==1)
			to_erase.push_back(i);
	}
	
	
	for(int i=0; i<int(to_erase.size()); i++)
		erase(to_erase[i]);
	
}



//*/



void module_collection::homeless(deque<int> & h) {
	
	h.clear();
	
	for(int i=0; i<int(memberships.size()); i++)
		if(memberships[i].size()<1)
			h.push_back(i);
	
	for(int i=0; i<int(modules.size()); i++) {
		
		if(modules[i].size()==1)
			h.push_back(modules[i][0]);
	
	}
	
	
	sort(h.begin(), h.end());
		

}


int module_collection::coverage() {
	
	
	// this function returns the number of nodes which are covered by at least one module
	
	int cov=0;
	for(int i=0; i<int(memberships.size()); i++)
		if(memberships[i].size()>0)
			cov++;
	
	
	
	return cov;


}



int module_collection::effective_groups() {
	
	
	
	int nmod=0;
	for(map<int, double >::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(modules[itm->first].size() > 1)		
		nmod++;
	
	return nmod;
	
	
	
}






void module_collection::set_partition(deque<deque<int> > & A) {

	A.clear();

	for(map<int, double >::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(modules[itm->first].size()>1)
		A.push_back(modules[itm->first]);
}


void module_collection::set_partition(deque<deque<int> > & A, deque<double> & b) {


	A.clear();
	b.clear();
	
	for(map<int, double>::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(modules[itm->first].size()>1){
		A.push_back(modules[itm->first]);
		b.push_back(module_bs[itm->first]);
	}





}






bool module_collection::check_already(const deque<int> & c) {

	// returns false if the module is already present
	
	
	
	map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
	for(int i=0; i<int(c.size()); i++) {
		
		for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
			int_histogram(*itj, com_ol);
		
	
	}
	
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++) {
		
		if(itm->second==int(c.size()) && itm->second==int(modules[itm->first].size()))
			return false;
			
	}
	
		
	return true;



}




































