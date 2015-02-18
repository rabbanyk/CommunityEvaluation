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
	if(argc < 3) show_help_short(argv);
	long t_begin=get_cpu_time();
	int randseed=int(time(NULL));
	ZRANDOMv3 rg(randseed);
	int nsample=10;
	if(argc>=4) nsample=atoi(argv[3]);
//	cout<<"nsample="<<nsample<<endl;
	
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

	double the_nmi=0;
	bool agtb=true; //qa is greater than qb
	if(qa<qb) agtb=false;
	if(agtb) the_nmi=compute_nmi(pa,pb);
	else the_nmi=compute_nmi(pb,pa);
//	cout<<"NMI="<<the_nmi<<endl;
	double tot_nmi=0;
	for(int sample=0;sample<nsample;sample++){
		double nmi=0;
		if(agtb){
			shuffle_seq(pb,rg);
			nmi=compute_nmi(pa,pb);
		}else{
			shuffle_seq(pa,rg);
			nmi=compute_nmi(pb,pa);
		}
//		cout<<nmi<<" "<<flush;
		tot_nmi += nmi;
	}
//	cout<<endl;
	tot_nmi /= nsample;
	cout<<the_nmi-tot_nmi<<endl;
//	cout<<"RNMI="<<the_nmi-tot_nmi<<endl;
//	cout<<"time used: "<<(get_cpu_time()-t_begin)/1000.0<<" seconds."<<endl;
}
//}}}

