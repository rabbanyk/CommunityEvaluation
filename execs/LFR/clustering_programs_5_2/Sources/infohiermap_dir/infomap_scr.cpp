
#include "standard_package/standard_include.cpp"







int main(int argc, char * argv[]) {
	
	
	if(argc<5) {
		
		cerr<<argv[0]<<" "<<"[filename] [seed] [number_of_runs] [path to the bin folder]"<<" -- output is infomap_net.net infomap_net.tree"<<endl;
		return 0;
	}
	
	
	pajek_format(string(argv[1]), true);
	int syy=system("mv net.paj infomap_net.net");
	
	string bbs= string(argv[4])+ "/infohiermap_dir "+ string(argv[2])+" infomap_net.net "+ string(argv[3]);
	
	cout<<"running: "<<bbs<<endl;
	int sy=system(bbs.c_str());
	
	
	return 0;




}





