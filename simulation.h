#ifndef SIMULATION_H
#def SIMULATION_H

class simulation
{
        private:  
                /***
                 * u_order defines the number of (manually set) arrays in the baths. 
                 * It is called u_order instead of u_samples, because it is not a sampling; this has to do with
                 *      the symmetry group we want to consider.
                 ***/ 
                const int u_order = 8;
                /***
                 *  There are L3 grid points, and 4 perturbations.
                 *  When sampling, each fluctuation is [the number of grid points, times the number
                 *      of perturbations possible, times tau] randomed perturbations on random sites. 
                 *  I.e. each grid point can possibly experience the 4 perturbations a hundred times each. 
                 ***/
                const int tau = 100;
                /***
                 *  A random quaternion is constructed from 3 numbers between 0 and 1.
                 *  We pre-generate a list of rotation matrices, the cache.
                 *  Each of the three numbers is split into rmc_number parts from 0 to 1.
                 *  This generates a list; we random one element of this list. The list has,
                 *      in total, rmc_number_total elements.
                 ***/
                const int rmc_number = 10; 
                const int rmc_number_total = rmc_number*rmc_number*rmc_number; 
                
                
                case 0:
                        return dsfmt_genrand_close_open(&dsfmt);
                break;
                case 1:
                        return std_random_mt(std_engine);
                break;
                case 2:
                        return boost_mt(boost_rng_fib);
                break;
                case 3:
                        return boost_mt(boost_rng_mt);
                break;
                default:
                        printf("Dice mode is not set... \n");
                        return 0.5;
                break;
                /***
                 * What random engine do we want to use?
                 * - 0, dsfmt
                 * - 1, std mt
                 * - 2, lagged lagged_fibonacci44497
                 * - 3, boost mt 
                 ***/
                const int dice_mode = 2;  
                
                /***
                 * The (seeded) generators for each random distribution.
                 * Also, the boost-uniform distribution.
                 ***/
                dsfmt_t dsfmt; 
                std::mt19937_64 std_engine(0);
                boost::random::mt19937 boost_rng_mt; 
                boost::random::lagged_fibonacci44497 boost_rng_fib; 
                boost::random::uniform_01<> boost_mt;  
                std::uniform_real_distribution<double> std_random_mt (0.0, 1.0);
                
                /***
                 * There are length points in each direction,
                 *      so length2 points in a plane,
                 *      so length3 points in a cube.
                 ***/
                int length_one = 1;
                int length_two = 1;
                int length_three = 1;
                
                /***
                 * Total energy of the simulation.
                 * Change from the previous total energy to the current one. (set in a perturbation or flip function)
                 * Some kind of energy unit
                 ***/
                double e_total;
                double e_change;
                double e_g;
                /***
                 * J1, J2, J3 parameters.
                 * 
                 ***/
                double j1 = 0.1;
                double j2 = 0.1;
                double j3 = 1;
                /***
                 * current temperature. Can be set through a setter function. (I require a setter function
                 *      so that it can be made clear that there hasn't been a calculation yet. This is done
                 *      through the calculated variable.)
                 ***/
                double beta; 
                bool calculated = false;
                
                /***
                 * When you change the temperature, the chances are (perhaps radically) altered. Therefore, you need to keep doing
                 *      random perturbations until it stabilises. The way this is done is looking at a number, e.g. 10, of perturbations.
                 * After 10 perturbations, you compare the previous and current energies. If the difference is less than accuracy percent,
                 *      it has stabilised.
                 ***/
                double accuracy;
                /***
                 * When looking at, for instance, the heat capacity, fluctuations are used. But how many fluctuations do you want to sample?
                 * Well, sample_amount, that is.
                 ***/
                int sample_amount;
                
                /***
                 * Yay, matrices. Let's go over these.
                 ***/
                double field_r[L3][9]                 = {{0}};
                double field_u_x[L3][9]                = {{0}};
                double field_u_y[L3][9]                = {{0}};
                double Uz[L3][9]                = {{0}};
                double U_bath[U_order][9]       ={{0}}; 
                double s[L3]                    = {0};
                double mpc_urx[L3][9]           = {{0}};
                double mpc_ury[L3][9]           = {{0}};
                double mpc_urz[L3][9]           = {{0}}; 
                double rmc_matrices[rmc_number_total][9] = {{0}};
                
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