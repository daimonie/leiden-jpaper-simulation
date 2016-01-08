#ifndef DATA_H
#define DATA_H
//just a data holder

#include <string>
#include <stdio.h>
#include <vector>
using namespace std;
class data
{
         public:
                double beta            = 0;
                double total_energy    = 0;
                double heat_capacity   = 0;
                double energy          = 0;
                double chi_energy      = 0;
                vector<double> chi_order;
                vector<double> order; 
                double j_one           = 0;
                double j_two           = 0;
                double j_three         = 0;  
                double accuracy        = 0;
                int sample_amount      = 0;
                string point_group     = "default";	
                void report ();
		void shout (FILE *);
};
#endif