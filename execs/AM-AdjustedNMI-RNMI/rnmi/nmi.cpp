/*
 *   nmi version 2, release date 09/05/2014
 *   Copyright 2014 Pan Zhang
 *   nmi is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or 
 *   (at your option) any later version.

 *   nmi is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
*/
//{{{header files
#include "nmi.h"
//}}}

//{{{void show_help_short()
void show_help_short(char **argv)
{
	cout<<"Calculate normalized mutual information between two configurations"<<endl;
	cout<<"Usage: "<<argv[0]<<" file1 file2"<<endl;
	exit(0);
}
//}}}

//{{{ int main(int argc, char** argv)
int main(int argc, char** argv)
{
	//read the first argument to set the function.
	if(argc != 3) show_help_short(argv);
	long t_begin=get_cpu_time();
	
	ifstream fa(argv[1]);
	ifstream fb(argv[2]);
	assert(fa.good()&&fb.good()&&"I can not open the file that contains the partition.");

	vector <string> pas;//partition a in string
	vector <string> pbs;//partition b in string
	string tmpstr;
	while(fa >> tmpstr) pas.push_back(tmpstr);
	while(fb >> tmpstr) pbs.push_back(tmpstr);
	assert(pas.size() == pbs.size() && "Two partitions have different number of nodes !");

	vector <int> pa;//partition a in number
	vector <int> pb;//partition b in number
	int qa=ps2p(pas,pa);
	int qb=ps2p(pbs,pb);

	//for(int i=0;i<pas.size();i++){
	//	cout<<pas[i]<<" "<<pa[i]<<endl;
	//}
	//cout<<q<<" groups"<<endl;
	double nmi=0;
	if(qa>qb) nmi=compute_nmi(pa,pb);
	else nmi=compute_nmi(pb,pa);
	cout<<nmi<<endl;
//	cout<<"NMI="<<nmi<<endl;
//	cout<<"time used: "<<(get_cpu_time()-t_begin)/1000.0<<" seconds."<<endl;
}
//}}}

