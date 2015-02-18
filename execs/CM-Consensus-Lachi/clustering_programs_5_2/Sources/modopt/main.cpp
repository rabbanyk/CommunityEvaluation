
/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *	This program is free software; you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 22/10/08 (email: arg.lanci@gmail.com)     *
 *	Modified on 22/10/08                                                          *
 *  Location: ISI foundation, Turin, Italy                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */




#include "standard_package/standard_include.cpp"
#include "static_network.h"
#include "max_mod_net.h"





int main(int argc, char * argv[]) {		
	
		
	
	if (argc<2) {
		
		cerr<<argv[0]<<" [netfile] [random_seed] [lambda] [trials] [temp_step] [initial_temp]"<<endl;
		cerr<<"default random_seed = 1"<<endl;
		cerr<<"default lambda= 1"<<endl;
		cerr<<"default trials= 5"<<endl;
		cerr<<"default temp_step= 0.999"<<endl;
		cerr<<"default initial_temp= 1e-6"<<endl;
		return -1;
	}
	
	
	
	
	string netfile(argv[1]);
	srand5(1);
	int random_seed=1;
	if (argc>=3) {
		
		string rs(argv[2]);
		random_seed=cast_int(cast_string_to_double(rs));
		srand5(random_seed);
	}
	
	double lambda=1;
	if (argc>=4) {
		string rs(argv[3]);
		lambda=cast_string_to_double(rs);
	}
	
	
	int trials=5;
	if (argc>=5) {
		string rs(argv[4]);
		trials=cast_int(cast_string_to_double(rs));
	}
	
	double temp_step=0.999;
	if (argc>=6) {
		string rs(argv[5]);
		temp_step=cast_string_to_double(rs);
	}
	
	double initial_temp=1e-6;
	if (argc>=7) {
		string rs(argv[6]);
		initial_temp=cast_string_to_double(rs);
	}
	
	max_mod_net luca(netfile);
	
	cout<<"random_seed = "<<random_seed<<endl;
	cout<<"lambda= "<<lambda<<endl;
	cout<<"trials= "<<trials<<endl;
	cout<<"temp_step= "<<temp_step<<endl;
	cout<<"initial_temp= "<<initial_temp<<endl;
	
	
	cout<<"net:: "<<luca.size()<<" nodes and "<<luca.edges()<<" edges ;\t average degree = "<<2*luca.edges()/luca.size()<<endl;
	
	deque<deque<int> > ten;
	
		
		
	luca.set_real(trials);
	luca.set_temp(initial_temp);
	luca.set_temp_step(temp_step);
	luca.LAMBDA=lambda;
	
	
	string part= netfile+"part";
	cout<<"maximum modularity "<<luca.iter_simannealing(false, ten, part)<<endl;
	
	
	
	
	
	return 0;
	
	
	
	
}



