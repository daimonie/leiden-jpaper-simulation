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
                
                void report ();
};

void data::report ()
{
        printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n", beta, total_energy, j_one, j_two, j_three);
}
#endif