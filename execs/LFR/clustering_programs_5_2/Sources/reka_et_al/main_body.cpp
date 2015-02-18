






#include "reka.h"


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
			return false;
		
		}	
	}	/* check if file_name exists */

	
	siglouvain luca;
	luca.set_graph(netfile);
	
	
	if(luca.size()==0 || luca.edges()==0) {
		cerr<<"network empty"<<endl;
		return -1;
	}
	
		
	cout<<"network of "<<luca.size()<<" nodes, average degree: "<<2*luca.edges()/luca.size()<<endl;
	
	
	int_matrix P;
	double qtot=luca.collect_raw_groups(P);
	cout<<"qtot: "<<qtot<<endl;
	ofstream pout(string(netfile+"part").c_str());
	luca.print_id(P, pout);
	
	return 0;
    
    
   	
	
	return 0;
	
	

}


