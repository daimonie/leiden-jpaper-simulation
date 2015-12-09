#include "simulation.h"

/***
 * 
 * Constructor
 *      Doesn't do anything, mostly because the parameters aren't set yet.
 ***/
simulation::simulation ()
{
        printf("Welcome to simulation. This seems to work. \n");
}
    
/*** 
 * Build the bath of field_u .This contains all possible u_fields.
 * Note that this is a manual calculation/setting for the symmetry group under consideration.
 ***/
void simulation::build_gauge_bath()
{ 
        /***
         * Note:
         * bath_field_u[0][0] = 1;  could be written this->bath_field_u[0][0] = 1; 
         * the latter is more clear, but not required and often considered to be bad practise.
         * ***/
        /** 0, Identity **/
        bath_field_u[0][0] = 1; 
        bath_field_u[0][4] = 1; 
        bath_field_u[0][8] = 1; 
        
        /** 1, C2,z **/
        bath_field_u[1][0] = -1; 
        bath_field_u[1][4] = -1; 
        bath_field_u[1][8] = 1;
        
        /** 2, -C4+,z **/
        bath_field_u[2][1] = 1; 
        bath_field_u[2][3] = -1; 
        bath_field_u[2][8] = -1;
        
        /** 3, -C4-,z **/
        bath_field_u[3][1] = -1; 
        bath_field_u[3][3] = 1; 
        bath_field_u[3][8] = -1;              
        
        /** 4, C2,y **/
        bath_field_u[4][0] = -1; 
        bath_field_u[4][4] = 1; 
        bath_field_u[4][8] = -1;                      
        
        /** 5, C2,x **/
        bath_field_u[5][0] = 1; 
        bath_field_u[5][4] = -1; 
        bath_field_u[5][8] = -1;      
        
        /** 6, m x,-x. z **/
        bath_field_u[6][1] = -1; 
        bath_field_u[6][3] = -1; 
        bath_field_u[6][8] = 1;                       

        /** 7, m x,x,z**/
        bath_field_u[7][1] = 1; 
        bath_field_u[7][3] = 1; 
        bath_field_u[7][8] = 1;        
        
}


/*** 
 * Uniform initialisation. This is the fully polarised state at zero temperature, 
 *      used for a 'cooling' simulation.
 ***/
void simulation::uniform_initialization()
{
        #pragma omp parallel for
        for(int i = 0; i < L3; i++)
        {
                 field_r[i][0] = 1;
                 field_r[i][4] = 1;
                 field_r[i][8] = 1;
                                 
                 field_s[i] = 1;
                 
                 field_u_x[i][0] = 1;
                 field_u_x[i][4] = 1;
                 field_u_x[i][8] = 1;
         
                 field_u_y[i][0] = 1;
                 field_u_y[i][4] = 1;
                 field_u_y[i][8] = 1;
         
                 field_u_z[i][0] = 1;
                 field_u_z[i][4] = 1;
                 field_u_z[i][8] = 1;                           
        } 

}

/*** 
 * Random initialisation. This is a temperature state, which means that it is random,
 *       with the chances influenced by the temperature through the Boltzmann factor.
 ***/
void simulation::random_initialization()
{
        int i, x, y, z;
        for (i = 0; i < L3; i++)
        {
                build_rotation_matrix(i, dice(), dice());
         
                x = int(u_order * dice());
                y = int(u_order * dice()); 
                z = int(u_order * dice()); 
                
                copy(begin(bath_field_u[x]),end(bath_field_u[x]),begin(field_u_x[i]));     
                copy(begin(bath_field_u[y]),end(bath_field_u[y]),begin(field_u_y[i]));     
                copy(begin(bath_field_u[z]),end(bath_field_u[z]),begin(field_u_z[i]));     
        } 
} 


/*** 
 * Energy calculation. The expression would be:
 *  s_left s_right U^mu_right R_right R_left.transpose
 ***/
double simulation::site_energy(int i) 
{
        /****** find neighbour, checked*****/
        
        
        x_prev = i % length_one == 0 ? i - 1 + length_one : i - 1;
        x_next = (i + 1) % length_one == 0 ? i + 1 - length_one : i + 1;
        
        y_prev = i % length_two < length_one ? i - length_one + length_two : i - length_one;
        y_next = (i + length_one) % length_two < length_one ? i + length_one - length_two : i + length_one;
        
        z_prev = i < length_two ? i - length_two + length_three : i - length_two;
        z_next = i + length_two >= length_three ? i + length_two - length_three : i + length_two; 
        
        double result[6][9] = {{0}};   
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[x_prev],
        mpc_urx[i], 3, field_r[x_prev],3,
        0.0, result[0],3);
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[i],
        mpc_urx[x_next], 3, field_r[i],3,
        0.0, result[1],3); 
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[y_prev],
        mpc_ury[i], 3, field_r[y_prev],3,
        0.0, result[2],3);
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[i],
        mpc_ury[y_next], 3, field_r[i],3,
        0.0, result[3],3);
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[z_prev],
        mpc_urz[i], 3, field_r[z_prev],3,
        0.0, result[4],3);
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[i],
        mpc_urz[z_next], 3, field_r[i],3,
        0.0, result[5],3); 
        
        /******** total energy *****/ 
        double energy = 0;
         
        for(int j = 0; j < 6; j++)
        {  
                energy += j_one*result[j][0] + j_two*result[j][4] + j_three*result[j][8];
        }                               
        return energy;
}     


/*** 
 * This is just the loop that perturbs n times. The logic for the 'inner'
 * and 'outer' perturbation sweeps is the same, so this is a good way to change these around.
 ***/ 
double simulation::thermalization_times ( int n)
{
        double energy = 1.0 ;
        int i, j;   
        for (i = 0; i <  n; i++)
        {
                energy += e_total; 
                for (j = 0; j < L3*4; j++)
                {   
                        flipper (dice(), dice(), dice(), dice(), dice());
                }                                               
        }  
        return energy;
}
/*** 
 * When you change the temperature, the simulation has to 'run' for a bit before things settle down in a state that 
 *       has the correct temperature. This is the 'running'.
 * It is settled down when current_energy / previous_energy is about the same, where about is within accuracy percent.
 ***/
void simulation::thermalization()
{
        double energy_ratio = 0;         
        double energy_two = 0.0;
        double energy_one = thermalization_times (thermalization_number_outer);
        
        while (energy_ratio < 1-accuracy || energy_ratio > 1+accuracy)
        {
                energy_two = energy_one;  
                energy_one = thermalization_times( thermalization_number_inner );
                
                energy_ratio = energy_one/energy_two; 
        }       
}
/*** 
 * After initialising, either uniformly or randomly, you need to cache these matrix products.
 * The idea of the Matrix Product Cache is that it reduces the need for taking this products,
 * which was often done 3-4 times before the state of this point changed.
 ***/

void simulation::mpc_initialisation()
{
        for(int i = 0; i < L3; i++)
        {
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3,3,3,s[i],
                field_u_x[i], 3, field_r[i],3,
                0.0, mpc_urx[i],3);
                
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3,3,3,field_s[i],
                field_u_y[i], 3, field_r[i],3,
                0.0, mpc_ury[i],3);
                
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3,3,3,field_s[i],
                field_u_z[i], 3, field_r[i],3,
                0.0, mpc_urz[i],3);       
        }
        
}
/*** 
 * The Random Matrix Cache represents the discretisation of the 'random' rotation matrix.
 ***/ 
void simulation::generate_rotation_matrices ()
{
        //printf("Generating %d rotation matrices.. \n", rmc_number);

        int ii, jj, kk;
        
        int rmc_index = 0;
        
        double w,x,y,z, one, two, three;  
        
        for (ii = 0; ii < rmc_number; ii++)
        {
                for (jj = 0; jj < rmc_number; jj++)
                {
                        for (kk = 0; kk < rmc_number; kk++)
                        {
                                //what number am I
                                rmc_index = rmc_number*rmc_number * ii + rmc_number * jj + kk;
                                
                                // arbitrary quaternions require 4 parameters; but our length is fixed. 
                                one = 1.0 / rmc_number * ii;
                                two = 1.0 / rmc_number * jj;
                                three = 1.0 / rmc_number * kk;

                                // http://planning.cs.uiuc.edu/node198.html 
                                w = sqrt(1-one) * sin(2 * M_PI * two);
                                x = sqrt(1-one) * cos(2 * M_PI * two);
                                y = sqrt(one) * sin(2 * M_PI * three);
                                z = sqrt(one) * cos(2 * M_PI * three);


                                //https://en.wikipedia.org/wiki/Rotation_group_SO%283%29#Quaternions_of_unit_norm

                                rmc_matrices[rmc_index][0] = 1 - 2 * (y*y+z*z);
                                rmc_matrices[rmc_index][1] = 2 * (x*y-z*w);
                                rmc_matrices[rmc_index][2] = 2 * (x*z+y*w);
                                rmc_matrices[rmc_index][3] = 2 * (x*y+z*w);
                                rmc_matrices[rmc_index][4] = 1 - 2 * (x*x+z*z);
                                rmc_matrices[rmc_index][5] = 2 * (y*z-x*w);
                                rmc_matrices[rmc_index][6] = 2 * (x*z-y*w);
                                rmc_matrices[rmc_index][7] = 2 * (y*z+x*w);
                                rmc_matrices[rmc_index][8] = 1 - 2 * (x*x+y*y); 
                                 
                        }
                }
        }   
}

/*** 
 * Sets field_r[i] to a random rotation matrix. 
 * Note that jactus means 'a throw of the dice'. As in [alea jactus est]. 
 ***/

void simulation::build_rotation_matrix(int i, double jactus) 
{         
        int build_random = (int) (rmc_number_total * jactus);
        
        copy( begin(rmc_matrices[build_random]), end(rmc_matrices[build_random]), begin(field_r[i]));
        s[i] = jactus2 > 0.5 ? 1 : -1; 
}
/*** 
 * Flipper contains the logic for a random perturbation. Note that doing the switch manually can save some random numbers,
 *      but this way ensures the logic is encapsulated.
 * You can also clearly see that I wanted to keep the casts of the dice separate from the logic.
 ***/
void simulation::flipper (double jactus_one, double jactus_two, double jactus_three, double jactus_four, double jactus_five)
{ 
        int site = int(length_three*jactus_two);  
        switch(int(4 * jactus_two)) 
        {  
                case 0 : flip_r(site, jactus_three, jactus_four, jactus_five); 
                break;
                case 1 : flip_u_x(site, jactus_three, jactus_four);
                break;
                case 2 : flip_u_y(site, jactus_three, jactus_four); 
                break;
                case 3 : flip_u_z(site, jactus_three, jactus_four); 
                break;                                                                          
        }           
        
}

/*** 
 * Perturbation of the rotation matrix field_r[i]
 ***/
  
void simulation::flip_r(int i, double jactus_one, double jactus_two, double jactus_three) 
{
        double e_old;
        double s_save;
        double r_save[9]; 
         
        e_old = site_energy(i); 
        
        copy(begin(field_r[i]), end(field_r[i]),begin(r_save));
        s_save = field_s[i];
        
        
        double tmp_urx[9] = {0};
        double tmp_ury[9] = {0};
        double tmp_urz[9] = {0};
        
        copy(begin(mpc_urx[i]), end(mpc_urx[i]), begin(tmp_urx)); 
        copy(begin(mpc_ury[i]), end(mpc_ury[i]), begin(tmp_ury)); 
        copy(begin(mpc_urz[i]), end(mpc_urz[i]), begin(tmp_urz));
         
        
        build_rotation_matrix(i, jactus_one, jactus_two); //saves to R[i] and s[i]
        
        //assign new value to matrix product cache
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_x[i], 3, field_r[i],3,
        0.0, mpc_urx[i],3);
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_y[i], 3, field_r[i],3,
        0.0, mpc_ury[i],3);
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_z[i], 3, field_r[i],3,
        0.0, mpc_urz[i],3);  
        
        double e_new = site_energy(i);  
        e_change = e_new - e_old;
          
        
        bool flip_accepted = true;
        double change_chance = exp(-beta * e_change);
        
        if (e_change >= 0)
        { 
                flip_accepted = false;
        }     
        if (change_chance > jactus_three) 
        {
                flip_accepted = true;
        }     
        
        if(flip_accepted)
        {
                e_total += e_change; 
                
        }
        else
        { 
                 
                copy(begin(r_save), end(r_save), begin(field_r[i])); 
                
                
                copy(begin(tmp_urx), end(tmp_urx), begin(mpc_urx[i])); 
                copy(begin(tmp_ury), end(tmp_ury), begin(mpc_ury[i])); 
                copy(begin(tmp_urz), end(tmp_urz), begin(mpc_urz[i]));
                
                field_s[i] = s_save; 
        } 
}
/*** 
 * Perturbation of the field_u-x
 ***/

void flip_u_x(int i, double jactus_one, double jactus_two)
{
        int x_prev = i % L == 0 ? i - 1 + L : i - 1; 
        
        double bond[9] = {0}; 
         
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[x_prev],
        mpc_urx[i], 3, field_r[x_prev],3,
        0.0, bond,3);
         
        double e_old = j_one*bond[0] + j_two * bond[4] + j_three*bond[8];
         
        double u_save[9];
        copy(begin(field_u_x[i]),end(field_u_x[i]),begin(u_save));
         
        int j = int(u_order * jactus_one);
        
        double tmp_urx[9] = {0};
        copy(begin(mpc_urx[i]), end(mpc_urx[i]), begin(tmp_urx));
        
        copy(begin(bath_field_u[j]),end(bath_field_u[j]),begin(field_u_x[i]));
         
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_x[i], 3, field_r[i],3,
        0.0, mpc_urx[i],3); 
          
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[x_prev],
        mpc_urx[i], 3, field_r[x_prev],3,
        0.0, bond,3);
        
        double e_new = j_one*bond[0] + j_two * bond[4] + j_three*bond[8];
        
        e_change = e_new - e_old;
          
        bool flip_accepted = true;
        double change_chance = exp(-beta * e_change);
        
        if (e_change >= 0)
        { 
                flip_accepted = false;
        }     
        if (change_chance > jactus_two) 
        {
                flip_accepted = true;
        }        
        
        if(flip_accepted)
        {
                e_total += e_change; 
                
        }
        else
        { 
                copy(begin(u_save),end(u_save),begin(field_u_x[i]));
                copy(begin(tmp_urx), end(tmp_urx), begin(mpc_urx[i])); 
        } 
}

/*** 
 * Perturbation of field_u_y
 ***/

void flip_Uy(int i, double jactus_one, double jactus_two)
{  
        int y_prev = i % L2 < L ? i - L + L2 : i - L;
         
        
        double bond[9] = {0}; 
         
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[yp],
        mpc_ury[i], 3, field_r[yp],3,
        0.0, bond,3);
         
        double e_old = j_one*bond[0] + j_two * bond[4] + j_three*bond[8];
         
        double u_save[9];
        copy(begin(Uy[i]),end(Uy[i]),begin(U_save));
         
        int j = int(u_order * jactus_one);
        
        double tmp_ury[9] = {0};
        copy(begin(mpc_ury[i]), end(mpc_ury[i]), begin(tmp_ury));
        
        copy(begin(U_bath[j]),end(U_bath[j]),begin(Uy[i]));
         
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_y[i], 3, field_r[i],3,
        0.0, mpc_ury[i],3); 
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3,field_s[yp],
        mpc_ury[i], 3, field_r[yp],3,
        0.0, bond,3);
        
        double e_new = j_one*bond[0] + j_two * bond[4] + j_three*bond[8];
        
        e_change = e_new - e_old;
          
        bool flip_accepted = true;
        double change_chance = exp(-beta * e_change);
        if (e_change >= 0)
        { 
                flip_accepted = false;
        }     
        if (change_chance > jactus_two) 
        {
                flip_accepted = true;
        }     
        
        if(flip_accepted)
        {
                e_total += e_change; 
                
        }
        else
        { 
                copy(begin(u_save),end(u_save),begin(field_u_y[i]));
                copy(begin(tmp_ury), end(tmp_ury), begin(mpc_ury[i])); 
        }   
}
/*** 
 ***/
/*** 
 ***/
/*** 
 ***/
/*** 
 ***/
/*** 
 ***/
