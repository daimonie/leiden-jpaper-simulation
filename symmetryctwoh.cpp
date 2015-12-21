#include "symmetry.h"
#include "symmetryctwoh.h"

void symmetry_c2h::bath()
{   
	bath_field[0][0] = 1; 
	bath_field[0][4] = 1; 
	bath_field[0][8] = 1; 
	 
	bath_field[1][0] = -1; 
	bath_field[1][4] = -1; 
	bath_field[1][8] = 1;	
 
	bath_field[2][0] = -1; 
	bath_field[2][4] = -1; 
	bath_field[2][8] = -1; 
	 
	bath_field[3][0] = 1; 
	bath_field[3][4] = 1; 
	bath_field[3][8] = -1;
}