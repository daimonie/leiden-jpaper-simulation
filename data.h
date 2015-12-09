#ifndef DATA_H
#def DATA_H
//just a data holder
class data
{
         public:
                beta            = 0;
                total_energy    = 0;
                heat_capacity   = 0;
                energy          = 0;
                chi_energy      = 0;
                chi_order       = 0;
                order           = 0;
                j_one           = 0;
                j_two           = 0;
                j_three         = 0;  
                
                void report ();
} 

void data::report ()
{
        printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n", beta, total_energy, j_one, j_two, j_three);
}
#endif