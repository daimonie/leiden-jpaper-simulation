#include <iostream>
#include <string>
#include <ctime>
#include "omp.h"

#include "data.h"
#include "simulation.h"
int main(int argc, char **argv)
{
        printf("Single object, this time. \n");
        
        simulation ares(4); //I like greek names
        
        auto time_start = std::chrono::high_resolution_clock::now();
         
        if (dice_mode == 0)
        {
                dsfmt_init_gen_rand(&dsfmt, time(0)); //seed dsfmt 
        }
        omp_set_num_threads(4); //I still want to use my computer :)
        generate_rotation_matrices(); 
        
        ares.build_gauge_bath();
        
        ares.j_one = -atof(argv[11]);
        ares.j_two = -atof(argv[12]);
        ares.j_three = -atof(argv[13]);
                 

        ares.sample_amount = atof(argv[14]);

        char *C_or_H = argv[1];
        
        beta_lower = atof(argv[2]);
        beta_upper = atof(argv[3]);
        
        if (*C_or_H == 'C') // if cooling
        {                       
                beta = beta_lower; // beta from small to big when cooling
                beta_step_small = atof(argv[4]);
                beta_step_big = atof(argv[5]);

                ares.random_initialization();
                ares.mpc_initialisation();
                
                for (int i=0; i < L3; i++)
                {
                        ares.e_total += site_energy(i);
                }

                ares.e_total /= 2;
                ares.e_g = L3*3*(J1+J2+J3);                          
        }
        else if (*C_or_H == 'H') // if heating
        {                       
                beta = beta_upper;
                beta_step_small = -atof(argv[4]); // the step is negative since beta is decreasing
                beta_step_big = -atof(argv[5]);

                uniform_initialization();
                mpc_initialisation();

                for (int i=0; i < L3; i++)
                {
                        ares.e_total += site_energy(i);
                }

                ares.e_total /= 2;
                ares.e_g = L3*3*(J1+J2+J3);                                                  
        }
        

                
        /**** initialize beta and the determine the fine region****/            
        beta_1 = atof(argv[6]);
        beta_2 = atof(argv[7]);
        
        ares.accuracy = atof(argv[8]);
        
        char *choice = argv[9]; 
        
        //I put in  argc here because of the annoying warning from the compiler.
        printf("Calculate from %2.3f to %2.3f, using {%2.3f, %2.3f, %2.3f}, accuracy %2.3f and %d samples. Argc=%d . \n", beta_lower, beta_upper, ares.j_one, ares.j_two, ares.j_three,
               ares.accuracy, ares.sample_amount, argc              
        );
        printf("Maximum cores %d \n", omp_get_max_threads());
        if(*choice == 'E')
        {
               while( beta <= beta_upper && beta >=beta_lower )
               {
                        auto results = ares.estimate_beta_c(); 
                        results.report();
                        
                        if ( (beta >= beta_1) && (beta <= beta_2) )
                        {
                                beta += beta_step_small;
                        }
                        else
                        {
                                beta += beta_step_big;
                        } 
                      
               }
        }
        
        auto time_end = std::chrono::high_resolution_clock::now();
        
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start).count();
        printf("Time taken is %ld microseconds. \n", microseconds);

        return 0;
}