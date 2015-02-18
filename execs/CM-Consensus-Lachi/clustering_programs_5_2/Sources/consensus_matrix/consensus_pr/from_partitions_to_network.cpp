
#include "standard_package/standard_include.cpp"
#include "projection.h"

int main(int argc, char * argv[]) {

	if(argc<3) {
		
		cerr<<argv[0]<<" [file with list of partitions to aggregate] [threshold] (strictly bigger links will be kept)"<<endl;
		return -1;
	}
	
	
	string thre_s(argv[2]);
	double thre=cast_string_to_double(thre_s);
	
	ifstream gin(argv[1]);
	deque<string> files;
	string s;
	while(gin>>s)
		files.push_back(s);
	
	if(files.size()==0) {
		cout<<"-> no partitions found"<<endl;
		return -1;
	}
		
	cout<<"files"<<endl;
	prints(files);
	
	deque<int_matrix> pts;
	map<int, int> all_elements;
	
	RANGE_loop(i, files) {
		int_matrix ten;
		get_partition_from_file(files[i], ten, 1);
		DI sizes;
		RANGE_loop(k, ten) {
			RANGE_loop(j, ten[k]) {
				if(all_elements.find(ten[k][j])==all_elements.end())
					all_elements.insert(make_pair(ten[k][j], all_elements.size()));
				
				ten[k][j]=all_elements[ten[k][j]];
			}
			sizes.push_back(ten[k].size());
		}
		
		pts.push_back(ten);
		cout<<"pt: "<<i<<" average_size: "<<average_func(sizes)<<endl;
	}
	
	
	DI ids(all_elements.size());
	IT_loop(mapii, itm, all_elements) {
		ids[itm->second]=itm->first;
	}
	

	link_collection built_net(all_elements.size());
	RANGE_loop(i, pts) 	built_net.add_partition(pts[i]);
	double nmi_ave=built_net.nmi();
	//cout<<"av nmi: "<<nmi_ave<<endl;
	
	string netfile(argv[1]);
	netfile=netfile+"_net.dat";
	ofstream neo(netfile.c_str());
	built_net.print(neo, ids, thre);
	//built_net.print_disconnected_nodes(neo, ids, thre);
	
	return 0;
	
}

