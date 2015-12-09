#ifndef SIMULATION_H
#def SIMULATION_H

class simulation
{
        private:  
                //
                const int U_order = 8;
                const int tau = 100;
                const int rmc_number = 10; 
                const int rmc_number_total = rmc_number*rmc_number*rmc_number;
                const int dice_mode = 2;  
                
                dsfmt_t dsfmt; 
                std::mt19937_64 std_engine(0);
                boost::random::mt19937 boost_rng_mt; 
                boost::random::lagged_fibonacci44497 boost_rng_fib; 
                boost::random::uniform_01<> boost_mt;  
                std::uniform_real_distribution<double> std_random_mt (0.0, 1.0);
                
                int L;
                int L2;
                int L3;
                
                double E_total;
                double E_change;
                double E_g;
                double J1 = 0.1;
                double J2 = 0.1;
                double J3 = 1;
                double beta; 
                double accurate;
                int sample_amount;
                double R[L3][9] = {{0}};
                double Ux[L3][9] = {{0}};
                double Uy[L3][9] = {{0}};
                double Uz[L3][9] = {{0}};
                double U_bath[U_order][9]={{0}}; 
                double s[L3] = {0};
                double rmc_matrices[rmc_number_total][9] = {{0}};
                double mpc_urx[L3][9] = {{0}};
                double mpc_ury[L3][9] = {{0}};
                double mpc_urz[L3][9] = {{0}}; 
                
                void build_gauge_bath();
                void uniform_initialization();
                void random_initialization();
                void build_rotation_matrix(int i, double, double);
                double site_energy(int i);   
                void flip_R(int, double, double);
                void flip_Ux(int, double, double);
                void flip_Uy(int, double, double);
                void flip_Uz(int, double, double);
                void thermalization();
                void estimate_beta_c(); 
                double orderparameter_n(); 
                void print_matrix(string, double[]);
                void mpc_initialisation(); 
                void generate_rotation_matrices();
                void flipper(double, double, double, double, double);
                bool josko_diagnostics ();
};
#endif