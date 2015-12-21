#include "symmetry.h"
#include "symmetrysfour.h"

void symmetry_s4::bath()
{  
	
	bath_field[0][0] = 1; 
	 
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

}