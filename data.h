#ifndef DATA_H
#define DATA_H
//just a data holder
class data
{
         public:
                double beta            = 0;
                double total_energy    = 0;
                double heat_capacity   = 0;
                double energy          = 0;
                double chi_energy      = 0;
                double chi_order       = 0;
                double order           = 0;
                double j_one           = 0;
                double j_two           = 0;
                double j_three         = 0;  
                double accuracy        = 0;
                int sample_amount      = 0;
                
                void report ();
};
#endif