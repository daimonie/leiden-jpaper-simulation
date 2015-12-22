#include "symmetry.h"
#include "symmetryctwo.h"

void symmetry_c2::bath()
{   
	bath_field[0][0] = 1; 
	bath_field[0][4] = 1; 
	bath_field[0][8] = 1; 
	 
	bath_field[1][0] = -1; 
	bath_field[1][4] = -1; 
	bath_field[1][8] = 1;	 
}