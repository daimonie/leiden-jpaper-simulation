#include "symmetry.h"
#include "symmetrydtwoh.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>

void symmetry_d2h::bath()
{    
	bath_field[0][0] = 1;  
	
	 
	bath_field[0][0] = 1; 
	bath_field[0][4] = 1; 
	bath_field[0][8] = 1; 
	 
	bath_field[1][0] = -1; 
	bath_field[1][4] = -1; 
	bath_field[1][8] = 1;
	 
	bath_field[2][0] = -1; 
	bath_field[2][4] = 1; 
	bath_field[2][8] = -1;			
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1.0,
	bath_field[1], 3, bath_field[2],3,
	0.0, bath_field[3],3);
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,-1.0,
	bath_field[0], 3, bath_field[0],3,
	0.0, bath_field[4],3);
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,-1.0,
	bath_field[1], 3, bath_field[0],3,
	0.0, bath_field[5],3);
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,-1.0,
	bath_field[2], 3, bath_field[0],3,
	0.0, bath_field[6],3);	
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,-1.0,
	bath_field[3], 3, bath_field[0],3,
	0.0, bath_field[7],3);
}