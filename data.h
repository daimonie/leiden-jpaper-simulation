#ifndef DATA_H
#define DATA_H
//just a data holder

#include <string>
class data
{
         public:
                double beta            = 0;
                double total_energy    = 0;
                double heat_capacity   = 0;
                double energy          = 0;
                double chi_energy      = 0;
                double chi_order_one   = 0;
                double order_one       = 0;
                double chi_order_two   = 0;
                double order_two       = 0;
                double j_one           = 0;
                double j_two           = 0;
                double j_three         = 0;  
                double accuracy        = 0;
                int sample_amount      = 0;
                std::string point_group     = 0;	
                void report ();
};
#endif