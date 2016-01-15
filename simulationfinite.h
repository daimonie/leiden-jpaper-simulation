#ifndef SIMULATION_FINITE_H
#define SIMULATION_FINITE_H

#include "simulation.h"

class simulation_finite : public simulation
{
    public:
        double energy_plaquette_xy (int i);
        double energy_plaquette_yz (int i);
        double energy_plaquette_zx (int i);
        
        double total_bond_energy(int i);
        
        double energy_total();
        
        void flip_u_x(int, double, double);
        void flip_u_y(int, double, double);
        void flip_u_z(int, double, double); 
        
        double finite_k = 0.0;
};


#endif