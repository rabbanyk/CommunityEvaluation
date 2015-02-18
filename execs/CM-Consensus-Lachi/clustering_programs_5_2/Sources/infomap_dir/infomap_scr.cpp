
#include "standard_package/standard_include.cpp"







int main(int argc, char * argv[]) {
	
	
	if(argc<5) {
		
		cerr<<argv[0]<<" "<<"[filename] [seed] [number_of_runs] [path to the bin folder]"<<" -- output is infomap_net.net infomap.part"<<endl;
		return 0;
	}
	
	
	pajek_format(string(argv[1]), true);
	int syy=system("mv net.paj infomap_net.net");
	
	string bbs= string(argv[4])+ "/infomap_dir "+ string(argv[2])+" infomap_net.net "+ string(argv[3]);
	
	cout<<"running: "<<bbs<<endl;
	int sy=system(bbs.c_str());
	
	
	
	
	
	
	deque<int> oldlabs;
	ifstream lin("labels.dat");
	int a;
	while(lin>>a)
		oldlabs.push_back(a);
	

		
	
	deque<deque<int> > one;
	get_partition_from_file_list_pajek("infomap_net.clu", one, oldlabs);
	ofstream oneout("infomap.part");
	for(int i=0; i<one.size(); i++)
		prints(one[i], oneout);

	
	
	
	return 0;




}





