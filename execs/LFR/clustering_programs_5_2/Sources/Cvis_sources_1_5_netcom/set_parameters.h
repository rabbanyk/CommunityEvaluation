


// to insert a new parameter there are four steps:
// 1- define it in the class
// 2- initialize it
// 3- set the flag
// 4- set the parameter






void general_program_statement(char * b) {
	
	
	cout<<"\nUSAGE 1: "<<b<<" -f network.dat -folder network.dat_oslo_files"<<endl<<endl;
	cout<<"USAGE 2: "<<b<<" -f network.dat -tree network.dat_tree"<<endl<<endl;
	
	cout<<"look at ReadMe.pdf for more information"<<endl;
	
	cout<<"\nOPTIONS"<<endl<<endl;
	cout<<"   [-labels filename]: filename should contain the labels: label - id"<<endl;
	cout<<"   [-opt U]: set the optimization parameter equal to U, default is U=2. bigger values lead to more accurate results and take more time"<<endl;
	cout<<"   [-seed U]: seed for the random number generator. default is read from file."<<endl;
	cout<<"   [-singlets U]: this is a parameter to locate singletons. default is 0.5 (must be U>=0). bigger values create a layout where singletons occupy more space (more CPU consuming)"<<endl;
	cout<<"   [-gap U]: this is to decide how far the modules should be. default is 1.4. reasonable values should be between 1 and 2"<<endl;
	
	cout<<"\nOUTPUT"<<endl<<endl;
	cout<<"all output files will be written in folder network.dat_cvis"<<endl;
	cout<<endl<<endl;
	
}




void error_statement(char * b) {

	cerr<<"\n\n************************************************************"<<endl;
	cerr<<"ERROR while reading parameters from command line... Please read program instructions or type: \n"<<b<<endl;
	cerr<<"************************************************************"<<endl;


}






class Parameters {
	
	
public:
	
	
	
	
	Parameters();
	~Parameters(){};
	void print();
	bool _set_(int argc, char * argv[]);
	bool set_flag_and_number(double & number_to_set, int & argct, int argc, char * argv[], double min_v, double max_v, string warning);
	bool set_flag_and_number(int & number_to_set, int & argct, int argc, char * argv[], int min_v, int max_v, string warning);
	

	//*******************************************************
	string file1;
	string label_file;
	string partitions_file;
	
	
	int seed_random;	
	
	double gap_between_frames;  	// this parameters is the factor for the minimum distance between the centers of f1 and f2: min_distance= gap_between_frames * (radius1 + radius2)

	
	double color_reduction_factor;
	double node_radius;
	
	double Spi_Sep;
	double Nodes_fir_C;
	
	
	bool weighted;
	bool labels_set;
	
	double edit_positioning_exponent;
	
	double percentage_fake_homeless;
	
	bool tp_flag;
	bool tree_flag;

	
	
private:
	
	map<string, int> command_flags;
		
	
};



bool Parameters::set_flag_and_number(int & number_to_set, int & argct, int argc, char * argv[], int min_v, int max_v, string warning) {
	
	argct++;
	if(argct==argc) {
		
		cout<<"you didn't set any number for the "<<warning<<endl;
		error_statement(argv[0]);
		return false;
	}
	
	string tt=argv[argct];
	double ttt;
	if(cast_string_to_double(tt, ttt)==false) {
		
		cout<<"you didn't set any number for the "<<warning<<endl;	
		error_statement(argv[0]);
		return false;
	}
	
	number_to_set=cast_int(ttt);
	
	if(number_to_set<min_v || number_to_set>max_v) {	
		cout<<"the "<<warning<<" must be between "<<min_v<<" and "<<max_v<<endl;
		error_statement(argv[0]);
		return false;
	}
	
	return true;
	
	
	
}


bool Parameters::set_flag_and_number(double & number_to_set, int & argct, int argc, char * argv[], double min_v, double max_v, string warning) {
	
	argct++;
	if(argct==argc) {
		
		cout<<"you didn't set any number for the "<<warning<<endl;
		error_statement(argv[0]);
		return false;
	}
	
	string tt=argv[argct];
	double ttt;
	if(cast_string_to_double(tt, ttt)==false) {
		
		cout<<"you didn't set any number for the "<<warning<<endl;	
		error_statement(argv[0]);
		return false;
	}
	
	number_to_set=ttt;
	
	if(number_to_set<min_v || number_to_set>max_v) {	
		cout<<"the "<<warning<<" must be between "<<min_v<<" and "<<max_v<<endl;
		error_statement(argv[0]);
		return false;
	}
	
	return true;
}




void Parameters::print() {
	
	
	cout<<"**************************************"<<endl;
	cout<<"Network file:\t\t\t"<<file1<<endl;
	
	
	
	
	if(labels_set)
		cout<<"file with labels: "<<label_file<<endl;

	cout<<"gap: "<<gap_between_frames<<endl;
	cout<<"optimization parameter: "<<edit_positioning_exponent<<endl;
	cout<<"singlets parameter: "<<percentage_fake_homeless<<endl;
	
	if(tp_flag)
		cout<<"tp folder: "<<partitions_file<<endl;
	
	if(tree_flag)
		cout<<"tree file: "<<partitions_file<<endl;
		
	if(seed_random!=-1)
		cout<<"Random number generator seed:\t\t\t"<<seed_random<<endl;
	
		
	cout<<"**************************************"<<endl<<endl;
	
	
	
	
	
}

Parameters::Parameters() {
	
	//**************************************************************************
	
	weighted=true;
	seed_random=-1;	
	
	gap_between_frames=1.4;	// this parameters is the factor for the minimum distance between the centers of f1 and f2: min_distance= gap_between_frames * (radius1 + radius2)	
	
	percentage_fake_homeless=0.5;
	
	
	color_reduction_factor= 0.99;
	node_radius=5;
	
	
	// spiral default parameters
	Spi_Sep=1;
	Nodes_fir_C=3;
	// spiral default parameters
	
	
	edit_positioning_exponent=0.5;
	
	labels_set=false;
	
	tp_flag=false;
	tree_flag=false;
	
	
	command_flags.insert(make_pair("-f", 3));
	command_flags.insert(make_pair("-labels", 4));	
	command_flags.insert(make_pair("-gap", 5));
	command_flags.insert(make_pair("-opt", 8));
	command_flags.insert(make_pair("-seed", 9));
	command_flags.insert(make_pair("-singlets", 10));
	command_flags.insert(make_pair("-folder", 11));
	command_flags.insert(make_pair("-tree", 12));

	
}




bool Parameters::_set_(int argc, char * argv[]) {
	
	int argct = 0;
	string temp;
	
	
	if (argc <= 1) {			/* if no arguments, return error_statement about program usage.*/
		error_statement(argv[0]);
		return false;
	}
	
	
	
	bool f_set=false;

	while (++argct < argc) {			// input file name
	
		cout<<"setting "<<argv[argct]<<endl;
		temp = argv[argct];
		map<string, int>::iterator itf=command_flags.find(temp);
		
		if(itf==command_flags.end()) {
			error_statement(argv[0]);
			return false;
		}
		
		int vp=itf->second;
		
		switch(vp) {
				
			case 3:
				argct++;
				if(argct==argc) {
					error_statement(argv[0]);
					return false;
				}
				file1=argv[argct];
				f_set=true;
				break;
			case 4:
				argct++;
				if(argct==argc) {
					error_statement(argv[0]);
					return false;
				}
				label_file=argv[argct];
				labels_set=true;
				break;
			case 5:
				if(set_flag_and_number(gap_between_frames, argct, argc, argv, 0.0001, 1000., "initial space between modules")==false)
					return false;
				break;
			case 8:
				if(set_flag_and_number(edit_positioning_exponent, argct, argc, argv, 0., 5., "optimization parameter")==false)
					return false;
				break;
			case 9:
				if(set_flag_and_number(seed_random, argct, argc, argv, 1, R2_IM2, "seed of the random number generator")==false)
					return false;
				break;
			case 10:
				if(set_flag_and_number(percentage_fake_homeless, argct, argc, argv, 0, 1000., " singletons parameter")==false)
					return false;
				break;
			case 11:
				argct++;
				if(argct==argc) {
					error_statement(argv[0]);
					return false;
				}
				partitions_file=argv[argct];
				tp_flag=true;
				break;
			case 12:
				argct++;
				if(argct==argc) {
					error_statement(argv[0]);
					return false;
				}
				partitions_file=argv[argct];
				tree_flag=true;
				break;
			default:
				error_statement(argv[0]);
				return false;		
		}

	
	
	
	
	}
	
	/*******************************************************************/
	
	
	if(f_set==false) {
		
		cerr<<"\n\n************************************************************"<<endl;
		cout<<"ERROR: you didn't set the file with the network.  Please read program instructions or type: \n"<<argv[0]<<endl;
		cerr<<"************************************************************"<<endl;
		
		return false;
		
	}
	
	
	if(tree_flag==false && tp_flag==false) {
		
		cerr<<"\n\n************************************************************"<<endl;
		cout<<"ERROR: you didn't set the file(s) with the partitions (option -folder OR -tree).  Please read program instructions or type: \n"<<argv[0]<<endl;
		cerr<<"************************************************************"<<endl;
		
		return false;
		
	}
	
	if(tree_flag==true && tp_flag==true) {
		
		cerr<<"\n\n************************************************************"<<endl;
		cout<<"ERROR: you cannot use both option -folder AND -tree.  Please read program instructions or type: \n"<<argv[0]<<endl;
		cerr<<"************************************************************"<<endl;
		
		return false;
		
	}
	
	
	
	
	if(seed_random==-1)
		srand_file();
	else
		srand5(seed_random);
	
	
	
	

	
	return true;
}





