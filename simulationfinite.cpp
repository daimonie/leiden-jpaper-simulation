#include "simulation.h"
#include "simulationfinite.h"
#include "data.h"
//includes from previous implementation
#include <iostream> 
#include <string> 
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <chrono> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <iomanip> 
#include <ctime>
#include <vector>  
#include <string>
#include "omp.h"

//random generation libraries
#include "dSFMT-src-2.2.3/dSFMT.h"
#include <random> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>

#include "order.h"
using namespace std;

simulation_finite::simulation_finite( int lattice_size ) : simulation (lattice_size)
{
    
}

double simulation_finite::energy_bond_x(int start, int end)
{ 
    double bond[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3, field_s[start]*field_s[end],
        mpc_urx[start], 3, field_r[end], 3,
        0.0, bond, 3);
    return j_one * bond[0] + j_two * bond[4] + j_three  * bond[8];  
}

double simulation_finite::energy_bond_y(int start, int end)
{ 
    double bond[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3, field_s[start] *field_s[end],
        mpc_ury[start], 3, field_r[end], 3,
        0.0, bond, 3);
    return j_one * bond[0] + j_two * bond[4] + j_three  * bond[8];  
}

double simulation_finite::energy_bond_z(int start, int end)
{ 
    double bond[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        3,3,3, field_s[start] * field_s[end],
        mpc_urz[start], 3, field_r[end], 3,
        0.0, bond, 3);
    return j_one * bond[0] + j_two * bond[4] + j_three  * bond[8];  
}

double simulation_finite::energy_plaquette_xy (int i)
{ 
    if( finite_k == 0 )
    {
        return 0.0;
    }
    int x_next = (i + 1) % length_one == 0 ? i + 1 - length_one : i + 1; 
    int y_next = (i + length_one) % length_two < length_one ? i + length_one - length_two : i + length_one; 
    
    double plaquette_part_one[9] = {0};
    double plaquette_part_two[9] = {0};
    double plaquette[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
    3,3,3,1,
    field_u_x[y_next], 3, field_u_y[i],3,
    0.0, plaquette_part_one,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_y[x_next], 3, plaquette_part_one,3,
    0.0, plaquette_part_two,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_x[i], 3, plaquette_part_two,3,
    0.0, plaquette,3);    
    
    return finite_k * ( plaquette[0] + plaquette[4] + plaquette[8]);
}
double simulation_finite::energy_plaquette_yz (int i)
{ 
    if( finite_k == 0 )
    {
        return 0.0;
    }
    int y_next = (i + length_one) % length_two < length_one ? i + length_one - length_two : i + length_one; 
    int z_next = i + length_two >= length_three ? i + length_two - length_three : i + length_two; 
    
    double plaquette_part_one[9] = {0};
    double plaquette_part_two[9] = {0};
    double plaquette[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
    3,3,3,1,
    field_u_y[z_next], 3, field_u_z[i],3,
    0.0, plaquette_part_one,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_z[y_next], 3, plaquette_part_one,3,
    0.0, plaquette_part_two,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_y[i], 3, plaquette_part_two,3,
    0.0, plaquette,3);    
    
    return finite_k * ( plaquette[0] + plaquette[4] + plaquette[8]);
}
double simulation_finite::energy_plaquette_zx (int i)
{
    
    if( finite_k == 0 )
    {
        return 0.0;
    }
    int x_next = (i + 1) % length_one == 0 ? i + 1 - length_one : i + 1; 
    int z_next = i + length_two >= length_three ? i + length_two - length_three : i + length_two; 
    
    double plaquette_part_one[9] = {0};
    double plaquette_part_two[9] = {0};
    double plaquette[9] = {0};
    
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
    3,3,3,1,
    field_u_z[x_next], 3, field_u_x[i],3,
    0.0, plaquette_part_one,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_x[z_next], 3, plaquette_part_one,3,
    0.0, plaquette_part_two,3);
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3,3,3,1,
    field_u_z[i], 3, plaquette_part_two,3,
    0.0, plaquette,3);    
     
    return finite_k * ( plaquette[0] + plaquette[4] + plaquette[8]);
}

double simulation_finite::total_bond_energy(int i)
{
        int x_prev = i % length_one == 0 ? i - 1 + length_one : i - 1;
        int x_next = (i + 1) % length_one == 0 ? i + 1 - length_one : i + 1;
        
        int y_prev = i % length_two < length_one ? i - length_one + length_two : i - length_one;
        int y_next = (i + length_one) % length_two < length_one ? i + length_one - length_two : i + length_one;
        
        int z_prev = i < length_two ? i - length_two + length_three : i - length_two;
        int z_next = i + length_two >= length_three ? i + length_two - length_three : i + length_two; 
        
        return energy_bond_x( x_prev, i) + energy_bond_x (i, x_next) + 
         energy_bond_y( y_prev, i) + energy_bond_y (i, y_next) + 
          energy_bond_z( z_prev, i) + energy_bond_z (i, z_next);
}
bool simulation_finite::flip_accepted(double e_change, double jactus_three)
{  
        bool result = true;
        double change_chance = exp(-beta * e_change);
        
        if (e_change >= 0)
        { 
                result = false;
        }     
        if (change_chance > jactus_three) 
        {
                result = true;
        }        
        
        return result;
}
void simulation_finite::flip_r(int i, double jactus_one, double jactus_two, double jactus_three)
{ 
        double e_old;
        double s_save;
        double r_save[9]; 
         
        e_old = total_bond_energy(i); 
        
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
        
        double e_new = total_bond_energy(i);  
        e_change = e_new - e_old;
         
        if(flip_accepted(e_change, jactus_three))
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

void simulation_finite::flip_u_x(int i, double jactus_one, double jactus_two)
{  
        int x_next = (i + 1) % length_one == 0 ? i + 1 - length_one : i + 1; 
        int y_prev = i % length_two < length_one ? i - length_one + length_two : i - length_one;  
        int z_prev = i < length_two ? i - length_two + length_three : i - length_two; 
         
        double e_old = energy_bond_x(i, x_next);
        
        e_old += energy_plaquette_zx(i);
        e_old += energy_plaquette_zx(z_prev);
        
        e_old += energy_plaquette_xy(i);
        e_old += energy_plaquette_xy(y_prev);
         
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
         
        double e_new = energy_bond_x(i, x_next);
        
        e_new += energy_plaquette_zx(i);
        e_new += energy_plaquette_zx(z_prev);
        
        e_new += energy_plaquette_xy(i);
        e_new += energy_plaquette_xy(y_prev);
        
        e_change = e_new - e_old;
           
        if(flip_accepted(e_change, jactus_two))
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

void simulation_finite::flip_u_y(int i, double jactus_one, double jactus_two)
{  
        int x_prev = i % length_one == 0 ? i - 1 + length_one : i - 1;  
        int y_next = (i + length_one) % length_two < length_one ? i + length_one - length_two : i + length_one; 
        int z_prev = i < length_two ? i - length_two + length_three : i - length_two; 
        
        double e_old = energy_bond_y(i, y_next);
         
        e_old += energy_plaquette_xy(i);
        e_old += energy_plaquette_xy(x_prev);
        
        e_old += energy_plaquette_yz(i);
        e_old += energy_plaquette_yz(z_prev);  
         
        double u_save[9];
        copy(begin(field_u_y[i]),end(field_u_y[i]),begin(u_save));
         
        int j = int(u_order * jactus_one);
        
        double tmp_ury[9] = {0};
        copy(begin(mpc_ury[i]), end(mpc_ury[i]), begin(tmp_ury));
        
        copy(begin(bath_field_u[j]),end(bath_field_u[j]),begin(field_u_y[i]));
         
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_y[i], 3, field_r[i],3,
        0.0, mpc_ury[i],3); 
         
        double e_new = energy_bond_y(i, y_next);
        
        e_new += energy_plaquette_xy(i);
        e_new += energy_plaquette_xy(x_prev);
        
        e_new += energy_plaquette_yz(i);
        e_new += energy_plaquette_yz(z_prev);   
        
        e_change = e_new - e_old; 
         
        if(flip_accepted(e_change, jactus_two))
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
 * Perturbation of field_u_z
 ***/

void simulation_finite::flip_u_z(int i, double jactus_one, double jactus_two)
{
        int x_prev = i % length_one == 0 ? i - 1 + length_one : i - 1;  
        int y_prev = i % length_two < length_one ? i - length_one + length_two : i - length_one;  
        int z_next = i + length_two >= length_three ? i + length_two - length_three : i + length_two; 
        
        
        double e_old = energy_bond_z(i, z_next);
        
        e_old += energy_plaquette_yz(i);
        e_old += energy_plaquette_yz(y_prev);
        
        e_old += energy_plaquette_zx(i);
        e_old += energy_plaquette_zx(x_prev);    
         
        double u_save[9];
        copy(begin(field_u_z[i]),end(field_u_z[i]),begin(u_save));
         
        int j = int(u_order * jactus_one);
        
        double tmp_urz[9] = {0};
        copy(begin(mpc_urz[i]), end(mpc_urz[i]), begin(tmp_urz));
        copy(begin(bath_field_u[j]),end(bath_field_u[j]),begin(field_u_z[i]));
         
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        3,3,3,field_s[i],
        field_u_z[i], 3, field_r[i],3,
        0.0, mpc_urz[i],3); 
        
        double e_new = energy_bond_z(i, z_next);
        
        e_new += energy_plaquette_yz(i);
        e_new += energy_plaquette_yz(y_prev);
        
        e_new += energy_plaquette_zx(i);
        e_new += energy_plaquette_zx(x_prev);    
        
        e_change = e_new - e_old;
        if(flip_accepted(e_change, jactus_two))
        { 
                e_total += e_change;  
        }
        else
        { 
                copy(begin(u_save),end(u_save),begin(field_u_z[i]));
                copy(begin(tmp_urz), end(tmp_urz), begin(mpc_urz[i]));
        }
        
} 
data simulation_finite::calculate()
{
    data results = simulation::calculate();
    results.finite_k = finite_k;
    
    return results;
} 

double simulation_finite::energy_total()
{ 
    return e_total;
}