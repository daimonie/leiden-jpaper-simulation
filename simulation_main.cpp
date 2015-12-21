#include <iostream>
#include <string>
#include <ctime>
#include "omp.h"

#include "data.h"
#include "simulation.h"  
#include "symmetry.h"
#include "symmetrydtwod.h"
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

//random generation libraries
#include "dSFMT-src-2.2.3/dSFMT.h"
#include <random> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
using namespace std;
int main(int argc, char **argv)
{
        printf("Single object, this time. \n");
        
        simulation ares(4); //I like greek names
        ares.dice_mode = 2;
        
        auto time_start = std::chrono::high_resolution_clock::now();
         
        if (ares.dice_mode == 0)
        {
                dsfmt_init_gen_rand(&ares.dsfmt, time(0)); //seed dsfmt 
        }
//         omp_set_num_threads(4);
        ares.generate_rotation_matrices(); 
        
	symmetry_d2d gauge; 
        ares.build_gauge_bath(gauge);
        
        ares.j_one = -atof(argv[11]);
        ares.j_two = -atof(argv[12]);
        ares.j_three = -atof(argv[13]);
                 

        ares.sample_amount = atof(argv[14]);

        char *C_or_H = argv[1];
        double beta_step_big, beta_step_small;
        double beta_lower = atof(argv[2]);
        double beta_upper = atof(argv[3]);
        
        if (*C_or_H == 'C') // if cooling
        {                       
                ares.beta = beta_lower; // beta from small to big when cooling
                beta_step_small = atof(argv[4]);
                beta_step_big = atof(argv[5]);

                ares.random_initialization();
                ares.mpc_initialisation();
                
                for (int i=0; i < ares.length_three; i++)
                {
                        ares.e_total += ares.site_energy(i);
                }
             
        }
        else// if (*C_or_H == 'H') // if heating
        {                       
                ares.beta = beta_upper;
                beta_step_small = -atof(argv[4]); // the step is negative since beta is decreasing
                beta_step_big = -atof(argv[5]);

                ares.uniform_initialization();
                ares.mpc_initialisation();

                for (int i=0; i < ares.length_three; i++)
                {
                        ares.e_total += ares.site_energy(i);
                }
                                                  
        }
        ares.e_total /= 2;
        ares.e_ground = ares.length_three*3*(ares.j_one + ares.j_two + ares.j_three);
        

                
        /**** initialize beta and the determine the fine region****/            
        double beta_1 = atof(argv[6]);
        double beta_2 = atof(argv[7]);
        
        ares.accuracy = atof(argv[8]);
        
        char *choice = argv[9]; 
         
        printf("Maximum cores %d \n", omp_get_max_threads());
        if(*choice == 'E')
        {
                string cooling_heating = "heating";
                if(*C_or_H == 'C')
                {
                        cooling_heating = "cooling";
                }
                printf("Calculating heap capacity, total energy for temperatures [%s from] %.3f to %.3f, J=diag(%.3f, %.3f, %.3f), accuracy %.3f, samples %d, size %d. (argc %d) \n",
                        cooling_heating.c_str(), beta_lower, beta_upper, ares.j_one, ares.j_two, ares.j_three, ares.accuracy, ares.sample_amount, ares.length_one, argc);

                while( ares.beta <= beta_upper && ares.beta >=beta_lower )
                {
                        ares.thermalization ();
                        auto results = ares.calculate(); 
                        results.report();
                        
                        if ( (ares.beta >= beta_1) && (ares.beta <= beta_2) )
                        {
                                ares.beta += beta_step_small;
                        }
                        else
                        {
                                ares.beta += beta_step_big;
                        } 
                        
                }
        }
        
        auto time_end = std::chrono::high_resolution_clock::now();
        
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
        printf("Time taken is %ld microseconds. \n", microseconds); 
        return 0;
}