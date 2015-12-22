#include "symmetry.h"
#include "symmetrystwo.h"

void symmetry_s2::bath()
{  
	
	bath_field[0][0] = 1; 
	bath_field[0][4] = 1; 
	bath_field[0][8] = 1; 
	 
	bath_field[1][0] = -1; 
	bath_field[1][4] = -1; 
	bath_field[1][8] = -1;	

}
string symmetry_s2::label()
{
	return string("S2");
}
int symmetry_s2::bath_size ()
{
	return 2;
}