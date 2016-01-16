#ifndef SIMULATION_FINITE_H
#define SIMULATION_FINITE_H

#include "simulation.h"
#include "data.h"

class simulation_finite : public simulation
{
    public:
        
        simulation_finite( int );
        double energy_bond_x(int, int);
        double energy_bond_y(int, int);
        double energy_bond_z(int, int);
        
        double energy_plaquette_xy (int);
        double energy_plaquette_yz (int);
        double energy_plaquette_zx (int);
        
        double total_bond_energy(int);
        
        double energy_total();
        
        void flip_r(int, double, double, double);
        void flip_u_x(int, double, double);
        void flip_u_y(int, double, double);
        void flip_u_z(int, double, double); 
        
        double finite_k = 0.0; 
        
        bool flip_accepted(double, double);
        
        data calculate();
};


#endif