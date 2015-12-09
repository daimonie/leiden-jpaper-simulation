#ifndef SIMULATION_H
#def SIMULATION_H

class simulation
{
        //everything is public, for simplicity.
        public:
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
                 * Energy unit or equivalently ground state energy
                 ***/
                double e_total;
                double e_change;
                double e_ground;
                /***
                 * J1, J2, J3 parameters.
                 * 
                 ***/
                double j_one = 0.1;
                double j_two = 0.1;
                double j_three = 1;
                /***
                 * current temperature.
                 ***/
                double beta;  
                
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
                 * For thermalisation, we first perturb _outer times and then _inner times, check for convergence, _inner times, ... until convergence.
                 ***/
                int thermalization_number_outer = 1000;
                int thermalization_number_inner = 10;
                /***
                 * Yay, matrices. Let's go over these.
                 * First, the r field; this is the rotation matrix R associated with a grid point.
                 * The u fields are the bonds. Bonds are associated with their 'left' neighbour; so the bond between i=1, j=2 is called 'bond 1'.
                 * The bath_field_u contains the possible values for the u fields.
                 * The mpc_ matrices are simply the product of the ising field with the u field with the r field. They are something used
                 *      in various calculations, so saving them saves up time (about twenty percent in the test).
                 * The Random Matrix Cache contains the possible values of the r field.
                 ***/
                double field_r[length_three][9]                 = {{0}};
                double field_u_x[length_three][9]               = {{0}};
                double field_u_y[length_three][9]               = {{0}};
                double field_u_z[length_three][9]               = {{0}};
                double bath_field_u[u_order][9]                 ={{0}}; 
                double field_s[length_three]                    = {0};
                double mpc_urx[length_three][9]                 = {{0}};
                double mpc_ury[length_three][9]                 = {{0}};
                double mpc_urz[length_three][9]                 = {{0}}; 
                double rmc_matrices[rmc_number_total][9]        = {{0}};
                
                //Function time. Their explanations are in the implementation file (simulation.cpp)
                simulation(int);
                void build_gauge_bath();
                void uniform_initialization();
                void random_initialization();
                double site_energy(int i);   
                void thermalization(); 
                double thermalization_times();
                void mpc_initialisation(); 
                void generate_rotation_matrices(); 
                void build_rotation_matrix(int i, double, double);
                void flipper(double, double, double, double, double); 
                void flip_r(int, double, double);
                void flip_u_x(int, double, double);
                void flip_u_y(int, double, double);
                void flip_u_z(int, double, double);
                double orderparameter_n(); 
                data estimate_beta_c(); 
};
#endif