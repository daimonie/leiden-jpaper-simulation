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
#include <iomanip> 
#include <ctime>
#include <string>
#include "omp.h"
#include <vector> 
#include <memory>  
#include <random> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
using namespace std;

int main(int argc, char* argv[])
{
	printf("$$ Welcome to speedtest. \n");
	if(argc > 0)
	{
		for(int j = 0; j < argc; j++)
		{
			string argument = string(argv[j]);
			printf("$$ Argument [%d]: %s \n", j, argument.c_str());
		}
		
	}
	auto time_start = std::chrono::high_resolution_clock::now();
	/***
	 * 	Note: Always preface  with $$ for random information.
	 * 	That way, the python functions will *ignore* the lines.
	 ***/
	printf("$$ Testing speed versus samples on d2d 4/6/8. \n");
	
	vector<vector<data>> results;
	 
	results.resize(3);
	 
	symmetry* gauge = new symmetry_d2d;
	 
	vector<vector <simulation>> sweeps;
	simulation tiny(4); 
	
	sweeps.resize(3);
	
	
	for(int i = 0; i < 3; i++)
	{
		int lattice_size = 2 * (i+1);
		for(int j = 0; j < 90; j++)
		{
			sweeps[i].emplace_back( simulation( lattice_size ) );
		}
	}
	
	for(unsigned int i = 0; i < sweeps.size(); i++)
	{
		for(unsigned int j = 0; j < sweeps.size(); j++)
		{
			sweeps[i][j].build_gauge_bath(gauge);
			sweeps[i][j].dice_mode = 2;
			sweeps[i][j].generate_rotation_matrices ();  
			sweeps[i][j].tau = 100; 
			sweeps[i][j].j_one = -0.1;
			sweeps[i][j].j_two = -0.1;
			sweeps[i][j].j_three = -1.0;
			sweeps[i][j].random_initialization ();  
			sweeps[i][j].mpc_initialisation ();   
			for (int ii = 0; ii < sweeps[i][j].length_three; ii++)
			{
				sweeps[i][j].e_total += sweeps[i][j].site_energy(ii);
			}
			sweeps[i][j].beta = 0.45;
			sweeps[i][j].e_total /= 2;
			sweeps[i][j].e_ground = sweeps[i][j].length_three*3*(sweeps[i][j].j_one + sweeps[i][j].j_two + sweeps[i][j].j_three); 
			
			
			sweeps[i][j].accuracy = 0.05;  
			sweeps[i][j].thermalization();
		}
	}
	
	auto time_end = std::chrono::high_resolution_clock::now();        
	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
	printf("$$ Initialisation took %ld microseconds. \n", microseconds);
	 
	for(unsigned int j = 0; j < sweeps.size(); j++)
	{
		vector<long int> timer_array;
		timer_array.resize(sweeps.size());
		for(unsigned int i = 0; i < sweeps.size(); i++)
		{
			time_start = std::chrono::high_resolution_clock::now();   
			
			sweeps[i][j].accuracy = 0.9 - 0.1 * j;  
			sweeps[i][j].thermalization();
			
			time_end = std::chrono::high_resolution_clock::now();   
			microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count(); 
			timer_array[i] = microseconds;
		}
		fprintf(stderr,"%d\t%ld\t%ld\t%ld\n", s, timer_array[0], timer_array[1], timer_array[2] );		
		printf("%.3f\t%ld\t%ld\t%ld\n", sweeps[0][j].accuracy, timer_array[0], timer_array[1], timer_array[2] );
	}
	//gauge was new'd, so it should be deleted
	delete gauge; 
	return 0;
}
 
