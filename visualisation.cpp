#include <iostream>
#include <string>
#include <ctime>
#include "omp.h"
#include <stdio.h>
#include "data.h"
#include "simulation.h"    
#include "symmetrydtwod.h" 
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

    if(argc != 1)
    {
        printf( "No arguments, please. \n"); 
        string arg_symmetry = string(argv[1]);
        string arg_size     = string(argv[2]);   
    }
    
    symmetry * gauge = new symmetry_d2d;
    
    int keep_at_temperature = 5;
    
    
    double beta_max = 3.0 ;
    double beta_number = beta_max / 0.25;
    
    
    int tmax = keep_at_temperature * (int) beta_number;
    
    if (tmax > 100)
    {
        fprintf(stderr, "Python will crash. Don't do this. \n");
        return 0;
    }
    
    printf("import numpy as np \n\n");
    printf("field_r = np.zeros((4, %d, 6, 6, 6, 3, 3)) \n", tmax);
    printf("beta = np.zeros((4, %d)) \n", tmax);
    printf("specific = np.zeros((4, %d)) \n", tmax);
    printf("jone = np.zeros((4)) \n");
     
    omp_set_num_threads(omp_get_max_threads());
    #pragma omp parallel for
    for(int i = 0; i < 4; i++) 
    { 	 
        simulation sweep( 6 ); 

        sweep.build_gauge_bath( gauge ); 
        
        sweep.tau = 100;
        sweep.dice_mode = 4;
        sweep.generate_rotation_matrices ();
        sweep.accuracy = 0.05;
 
        sweep.j_one = -2.0 / 4.0 * i; 
        sweep.j_two = sweep.j_one;
        sweep.j_three = -1.0;
        sweep.sample_amount = 100; 

        //start at beta=0 (T=inf), then go to Bmax, then go down again. Hence, random.
        sweep.random_initialization ();
        sweep.mpc_initialisation ();

        for (int ii = 0; ii < sweep.length_three; ii++)
        {
            sweep.e_total += sweep.site_energy(ii);
        }
        sweep.e_total /= 2;
        sweep.e_ground = sweep.length_three*3*(sweep.j_one + sweep.j_two + sweep.j_three); 


        // cool it down (Tinf -> T finite)
        sweep.beta = 0.0;
        int tnumber = 0;
        
        printf("jone[%d] = %.3f \n", i, sweep.j_one);
        while( sweep.beta <= beta_max && tnumber < tmax)
        { 	
            sweep.beta += beta_max / beta_number;
            for(int keep = 0; keep < keep_at_temperature; keep++)
            {
                sweep.thermalization ();   
                printf("beta[%d][%d] = %.3f \n", i, tnumber, sweep.beta);
                
                auto result = sweep.calculate();
                printf("specific[%d][%d] = %.3f \n", i, tnumber, result.heat_capacity);
                
                int x = 0;
                int y = 0; 
                int z = 0;
                for(int site = 0; site < sweep.length_three; site++)
                { 
                    if (x >= sweep.length_one)
                    {
                        x = 0;
                        y += 1;
                    }
                    if (y >= sweep.length_one)
                    { 
                        y = 0;
                        z += 1;
                    }
//                     fprintf(stderr, "%d %d %d\n", x,y, z); seems to work. Nobody cares about the order if it consistent, right?
//                         comes down to a global rotation anyway.

                    int matrix_index = 0;
                    for( int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            double matrix_value = sweep.field_r[site][matrix_index];
                            printf("field_r[%d][%d][%d][%d][%d][%d][%d] = %.3f \n", i, tnumber, x, y, z, m, n, matrix_value);
                            matrix_index += 1;
                        }
                    }
                    
                    x += 1;
                }
                tnumber += 1;
            }
        }    
    }  
    delete gauge; 
    auto time_end = std::chrono::high_resolution_clock::now();        
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
    fprintf(stderr, "Elapsed time %ld microseconds. \n", microseconds); 
    return 0;
}

