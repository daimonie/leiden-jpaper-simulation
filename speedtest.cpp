#include <iostream>
#include <string>
#include <ctime>
#include "omp.h"

#include "data.h"
#include "simulation.h"   
#include "symmetryctwo.h"
#include "symmetryctwoh.h"
#include "symmetryctwov.h"
#include "symmetrydtwod.h"
#include "symmetrydtwoh.h"
#include "symmetrystwo.h"
#include "symmetrysfour.h"
//includes from previous implementation
#include <iostream>
#include <fstream> 
#include <string> 
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <chrono> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <iomanip> 
#include <ctime>
#include <string>
#include "omp.h"
#include <vector> 
#include <memory> 
//random generation libraries
#include "dSFMT-src-2.2.3/dSFMT.h"
#include <random> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
using namespace std;

int main(int argc, char* argv[])
{
	/***
	 * 	Note: Always preface  with $$ for random information.
	 * 	That way, the python functions will *ignore* the lines.
	 ***/
	printf("$$ small speedtest on d2d 4/6/8. \n");
	auto time_start = std::chrono::high_resolution_clock::now();
	
	vector<vector<data>> results;
	 
	results.resize(3);
	 
	symmetry* gauge = new symmetry_d2d;
	
	simulation tiny = new simulation(4);
	simulation small = new simulation(6);
	simulation medium = new simulation(8);
	
	tiny.dice_mode = 2;
	tiny.generate_rotation_matrices ();  
	tiny.tau = 100; 
	tiny.build_gauge_bath(gauge);
	tiny.j_one = -0.1;
	tiny.j_two = -0.1;
	tiny.j_three = -1.0;
	tiny.random_initialization ();  
	tiny.mpc_initialisation ();  
	tiny.beta = 0.5;
	tiny.thermalization();
	tiny.e_total /= 2;
	tiny.e_ground = tiny.length_three*3*(tiny.j_one + tiny.j_two + tiny.j_three); 
	tiny.accuracy = 0.05;
	
	
	small.dice_mode = 2;
	small.generate_rotation_matrices ();  
	small.tau = 100; 
	small.build_gauge_bath(gauge);
	small.j_one = -0.1;
	small.j_two = -0.1;
	small.j_three = -1.0;
	small.random_initialization ();  
	small.mpc_initialisation ();  
	small.beta = 0.5;
	small.thermalization();
	small.e_total /= 2;
	small.e_ground = small.length_three*3*(small.j_one + small.j_two + small.j_three); 
	small.accuracy = 0.05;
	
	
	
	medium.dice_mode = 2;
	medium.generate_rotation_matrices ();  
	medium.tau = 100; 
	medium.build_gauge_bath(gauge);
	medium.j_one = -0.1;
	medium.j_two = -0.1;
	medium.j_three = -1.0;
	medium.random_initialization ();  
	medium.mpc_initialisation ();  
	medium.beta = 0.5;
	medium.thermalization();
	medium.e_total /= 2;
	medium.e_ground = medium.length_three*3*(medium.j_one + medium.j_two + medium.j_three); 
	medium.accuracy = 0.05;
	
	
	auto time_end = std::chrono::high_resolution_clock::now();        
	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
	printf("$$ Initialisation took %d microseconds. \n", microseconds);
	
	for( int s = 10; s < 2000; s += 10)
	{  
		auto tiny_start = std::chrono::high_resolution_clock::now();   
		tiny.samples = s;
		results[0].push_back(sweeps[i].calculate ());  
		auto tiny_end = std::chrono::high_resolution_clock::now();   
		auto tiny_microseconds = std::chrono::duration_cast<std::chrono::microseconds>( tiny_start - tiny_end).count();
		
		
		auto small_start = std::chrono::high_resolution_clock::now();   
		small.samples = s;
		results[0].push_back(sweeps[i].calculate ());  
		auto small_end = std::chrono::high_resolution_clock::now();   
		auto small_microseconds = std::chrono::duration_cast<std::chrono::microseconds>( small_start - small_end).count();
		
		
		auto medium_start = std::chrono::high_resolution_clock::now();   
		medium.samples = s;
		results[0].push_back(sweeps[i].calculate ());  
		auto medium_end = std::chrono::high_resolution_clock::now();   
		auto medium_microseconds = std::chrono::duration_cast<std::chrono::microseconds>( medium_start - medium_end).count();
		
		printf("%d\t%d\t%d\t%d\n", s, tiny_microseconds, small_microseconds, medium_microseconds);		
	} 
        auto time_end = std::chrono::high_resolution_clock::now();        
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
        printf("$$ Elapsed time %ld microseconds. \n", microseconds); 
	//gauge was new'd, so it should be deleted
	delete gauge; 
	return 0;
}
 
