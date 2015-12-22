#include "symmetry.h"
#include "symmetryctwov.h"

void symmetry_c2v::bath()
{    
	bath_field[0][0] = 1; 
	bath_field[0][4] = 1; 
	bath_field[0][8] = 1; 
	 
	bath_field[1][0] = -1; 
	bath_field[1][4] = -1; 
	bath_field[1][8] = 1;	
 
	bath_field[2][0] = 1; 
	bath_field[2][4] = -1; 
	bath_field[2][8] = 1; 
	 
	bath_field[3][0] = -1; 
	bath_field[3][4] = 1; 
	bath_field[3][8] = 1;
}
string symmetry_c2v::label()
{
	return string("C2v");
}
int symmetry_c2v::bath_size ()
{
	return 4;
	
}