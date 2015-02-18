


#include "siglouvain.h"

#include <sstream>

template <class T>
inline string to_string (const T & t) {
	stringstream ss;
	ss<<t;
	return ss.str();
}

double qtot_higher_level(siglouvain & luca, int_matrix & P, int_matrix & short_tp, int_matrix & newP) {
	
	
	short_tp.clear();
	
	
	map<int, map<int, double > > UPnet;
	luca.set_upper_network(UPnet, P);
	siglouvain luca2;
	luca2.set_graph(UPnet);
	luca2.collect_raw_groups(short_tp);
	
	newP.clear();
	RANGE_loop(i, short_tp) {
		set<int> s;
		RANGE_loop(j, short_tp[i]) {
			deque_to_set_app(P[short_tp[i][j]], s);
		}
		DI pp;
		set_to_deque(s, pp);
		newP.push_back(pp);
	}
	
	
	return luca.newman_modularity(newP);

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
	
	
	{	/* check if file_name exists */
		char b[netfile.size()+1];
		cast_string_to_char(netfile, b);
		ifstream inb(b);
		if(inb.is_open()==false) {
			
			cout<<"File "<<netfile<<" not found"<<endl;
			return -1;
		
		}	
	}	/* check if file_name exists */

	
	
	siglouvain luca;
	luca.set_graph(netfile);
	
	
	if(luca.size()==0 || luca.edges()==0) {
		cerr<<"network empty"<<endl;
		return -1;
	}
	
	
	cout<<"network of "<<luca.size()<<" nodes, average degree: "<<2*luca.edges()/luca.size()<<endl;
	int sy=system(("rm "+paras.folder+"tp").c_str());
	sy=system(("rm "+paras.folder+"short_tp*").c_str());
	
	int_matrix P;
	double qtot=luca.collect_raw_groups(P);
	cout<<"qtot: "<<qtot<<endl;
	ofstream pout((paras.folder+"tp").c_str());
	luca.print_id(P, pout);
	
	int level=0;
	
	while (true) {
		
		cout<<"level "<<level<<" qtot: "<<qtot<<" modules: "<<P.size()<<endl;

		int_matrix short_tp;
		int_matrix newP;
		double q2=qtot_higher_level(luca, P, short_tp, newP);
		if(q2<qtot)
			break;
		else {
			
			++level;
			qtot=q2;
			P=newP;
			ofstream poo((paras.folder+"short_tp"+to_string(level)).c_str());
			RANGE_loop(i, short_tp) {
				prints(short_tp[i], poo);
			}
			
			ofstream pout2((paras.folder+"tp"+to_string(level)).c_str());
			luca.print_id(P, pout2);

			
			
		}
	}
	
		
	return 0;
    
    
   	
	
	
	

}


