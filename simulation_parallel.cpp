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
	
        auto time_start = std::chrono::high_resolution_clock::now();
	//First, check if the arguments are proper.
	//Return 0 terminates program.  
	if(argc != 3)
	{
		printf("$$ Three arguments are required, yet %d were given. \n", argc);
		return 0;
	} 
	
	string arg_symmetry	= string(argv[1]);
	string arg_size		= string(argv[2]);
	
	int samples = 10;
	if(arg_size == "small")
	{ 
		samples = 5000;
		printf("$$ Will simulate small (8) lattice for point group %s. \n", arg_symmetry.c_str()); 
	}
	else if(arg_size == "large")
	{
		samples = 2000;
		printf("$$ Will simulate large (12) lattice for point group %s. \n", arg_symmetry.c_str());
	}
	else if(arg_size == "test")
	{
		samples = 10;
		printf("$$ Will simulate test (4) lattice for point group %s. \n", arg_symmetry.c_str());
	}
	else
	{
		printf("$$ Unknown lattice size. \n");
		return 0;
	}
	
	if(arg_symmetry != "c2" && arg_symmetry != "c2h" && arg_symmetry != "c2v" && arg_symmetry != "d2d"
		&& arg_symmetry != "d2h" && arg_symmetry != "s2" && arg_symmetry != "s4")
	{
		printf("$$ Unknown point group. \n");
		return 0;
	} 
	//if the code is not terminated before this point, the arguments are well formed and we can just start calculating.
	
	//This should be fine, given that we copy from this object. Right?
	symmetry* gauge;
	
	if(arg_symmetry == "c2")
	{
		gauge = new symmetry_c2;
	}
	else if(arg_symmetry == "c2v")
	{
		gauge = new symmetry_c2v;
	}
	else if(arg_symmetry == "c2h")
	{
		gauge = new symmetry_c2h;
	}
	else if(arg_symmetry == "d2d")
	{
		gauge = new symmetry_d2d;
	}
	else if(arg_symmetry == "d2h")
	{
		gauge = new symmetry_d2h;
	}
	else if(arg_symmetry == "s2")
	{
		gauge = new symmetry_s2;
	}
	else if(arg_symmetry == "s4")
	{
		gauge = new symmetry_s4;
	}
	else
	{
		printf("$$ Something went wrong. Terminating. \n");
		return 0;
	}
	printf("$$ Setting gauge group to %s. \n", gauge->label().c_str());
	
	/***
	 * We will simulate 0.1 < J < 1 at dJ = 0.01, and 1 < J < 2 at Dj=0.05.
	 * We will thus require 99+20 values.	  
	 ***/
	int imax = 119;

	vector<simulation> sweeps;
	vector<vector<data>> results;
	 
	results.resize(imax);
	 
	int lattice_size = 8;
	if(arg_size == "large")
	{ 
		lattice_size = 12;
	}   
	else if(arg_size == "test")
	{ 
		lattice_size = 4;
		imax = 64;
	}   
	
	for(int i = 0; i < imax; i++)
	{
		sweeps.emplace_back( simulation( lattice_size ) );
	}
	//Let's do this before OMP.
	for(unsigned int i = 0; i < sweeps.size(); i++)
	{
                sweeps[i].dice_mode = 2;
                sweeps[i].generate_rotation_matrices (); 
		
		sweeps[i].build_gauge_bath(gauge);
	}
	
        omp_set_num_threads(4);
        #pragma omp parallel for
        for(unsigned int i = 0; i < sweeps.size(); i++) 
        { 		
		if( i > 99)
		{
			sweeps[i].j_one = 1.0 + 0.05 * (i-98); 
		}
		else
		{ 
			sweeps[i].j_one = 0.01 * (i+1); 		
		}
		
                sweeps[i].j_two = sweeps[i].j_one;
                sweeps[i].j_three = 1.0;
		 
		sweeps[i].sample_amount = samples; 
			
                //start at beta=0 (T=inf), then go to Bmax, then go down again. Hence, random.
		sweeps[i].random_initialization ();
                sweeps[i].mpc_initialisation ();
		
                for (int ii = 0; ii < sweeps[i].length_three; ii++)
                {
                        sweeps[i].e_total += sweeps[i].site_energy(ii);
                }
                sweeps[i].e_total /= 2;
                sweeps[i].e_ground = sweeps[i].length_three*3*(sweeps[i].j_one + sweeps[i].j_two + sweeps[i].j_three); 
		
                sweeps[i].accuracy = 0.05;
		if(arg_size == "test")
		{
			sweeps[i].accuracy = 0.5;
		}
                
		double beta_max = (-8.5)/(1.99)*sweeps[i].j_one + 10.0; 
		if(beta_max > 10.0 || beta_max < 0)
		{
			printf("$$ Warning: beta_max=%.3f out of bounds.\n", beta_max);
		}
		//recall beta = inv T
		double beta = 0.0;  
		// cool it down (Tinf -> T finite)
		while( beta <= beta_max )
		{
			if(arg_size == "test")
			{
				beta = beta + 0.5;
			}
			else
			{ 
				beta = beta + 0.05;
			}
			sweeps[i].beta = beta;  
			sweeps[i].thermalization (); 
		
			results[i].push_back(sweeps[i].calculate ());  
		}
		//Heat it up (T finite -> T inf)
		while( beta > 0 )
		{
			if(arg_size == "test")
			{
				beta = beta - 0.5;
			}
			else
			{ 
				beta = beta - 0.05;
			}
			
			sweeps[i].beta = beta;  
			sweeps[i].thermalization (); 
		
			results[i].push_back(sweeps[i].calculate ());  
		}            
        }
        
        for(unsigned int i = 0; i < sweeps.size(); i++)
        {
                printf("$$ %.3f\t%s\t%d\n", sweeps[i].j_one, sweeps[i].u_label.c_str(), sweeps[i].u_order);
                for(unsigned int jj = 0; jj < results[i].size(); jj++)
                {
                        auto result = results[i][jj];
                        result.report ();
                }
        }
	//report time, end program 
        auto time_end = std::chrono::high_resolution_clock::now();        
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
        printf("$$ Elapsed time %ld microseconds. \n", microseconds); 
	//gauge was new'd, so it should be deleted
	delete gauge; 
	return 0;
}
 