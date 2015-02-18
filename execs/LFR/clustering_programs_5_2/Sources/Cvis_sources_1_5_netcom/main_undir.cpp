




#include "standard_package/standard_include.cpp"
#define twopi 6.28318


#include "set_parameters.h"

Parameters paras;
char * folder_output;


#include "module_collection.h"
#include "undirected_network.h"
#include "netvi_undir.h"
#include "get_no_overlapping_parts.h"
#include "convert_tree_format.h"



void program_statement(char * b) {
	
	
	cout<<"\n\n\n***************************************************************************************************************************************************"<<endl;
	cout<<"This program implements the hier-circle-visualization method for undirected networks"<<endl;	
	general_program_statement(b);

	
}


#include "main_body.cpp"







