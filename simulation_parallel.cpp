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
//random generation libraries
#include "dSFMT-src-2.2.3/dSFMT.h"
#include <random> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
using namespace std;

int main()
{ 
        auto time_start = std::chrono::high_resolution_clock::now();
        vector <simulation> sweeps;
        vector <vector<data>> results;
        
        int imax = 7;
	unsigned int j = 0;
        results.resize(imax);
	
        for(int i = 0; i < imax; i++)
        {
                sweeps.emplace_back( simulation(4) );
                
        }
        
        omp_set_num_threads(4);
        #pragma omp parallel for private(j)
        for(int i = 0; i < imax; i++)
        {
		
		
                sweeps[i].dice_mode = 2;
                sweeps[i].generate_rotation_matrices (); 
		
		if( i == 0 )
		{
			symmetry_c2 gauge; 
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 1 )
		{
			symmetry_c2h gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 2 )
		{
			symmetry_c2v gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 3 )
		{
			symmetry_d2d gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 4 )
		{
			symmetry_d2h gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 5 )
		{
			symmetry_s2 gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
		else if( i == 6 )
		{
			symmetry_s4 gauge;
			sweeps[i].build_gauge_bath (gauge);
		}
                
                sweeps[i].j_one = 0.5;
                
                sweeps[i].j_two = sweeps[i].j_one;
                sweeps[i].j_three = 1.0;
                sweeps[i].sample_amount = 10;
                sweeps[i].random_initialization ();
                sweeps[i].mpc_initialisation ();
                for (int ii = 0; ii < sweeps[i].length_three; ii++)
                {
                        sweeps[i].e_total += sweeps[i].site_energy(ii);
                }
                sweeps[i].e_total /= 2;
                sweeps[i].e_ground = sweeps[i].length_three*3*(sweeps[i].j_one + sweeps[i].j_two + sweeps[i].j_three);
                sweeps[i].accuracy = 0.5;
                
                for(j = 0; j < 20; j++)
                {
                        sweeps[i].beta = 5.0/20 * j; 
                        sweeps[i].thermalization (); 
                        
                        results[i].push_back(sweeps[i].calculate ());
                }
        }
        
        for(int i = 0; i < imax; i++)
        {
                printf("Report J1/J3 = %.3f .\n", sweeps[i].j_one);
                for(j = 0; j < results[i].size(); j++)
                {
                        auto result = results[i][j];
                        result.report ();
                }
        }
        
        
        
        auto time_end = std::chrono::high_resolution_clock::now();
        
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
        printf("Time taken is %ld microseconds. \n", microseconds); 
        
        return 0;
}
 