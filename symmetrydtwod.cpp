#include "symmetry.h"
#include "symmetrydtwod.h"

void symmetry_d2d::bath()
{ 
        bath_field[0][0] = 1; 
        bath_field[0][4] = 1; 
        bath_field[0][8] = 1; 
         
        bath_field[1][0] = -1; 
        bath_field[1][4] = -1; 
        bath_field[1][8] = 1;
         
        bath_field[2][1] = 1; 
        bath_field[2][3] = -1; 
        bath_field[2][8] = -1;
         
        bath_field[3][1] = -1; 
        bath_field[3][3] = 1; 
        bath_field[3][8] = -1;              
         
        bath_field[4][0] = -1; 
        bath_field[4][4] = 1; 
        bath_field[4][8] = -1;                      
         
        bath_field[5][0] = 1; 
        bath_field[5][4] = -1; 
        bath_field[5][8] = -1;      
         
        bath_field[6][1] = -1; 
        bath_field[6][3] = -1; 
        bath_field[6][8] = 1;                       
 
        bath_field[7][1] = 1; 
        bath_field[7][3] = 1; 
        bath_field[7][8] = 1;   
}