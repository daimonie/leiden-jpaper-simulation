/****
 * argv = filename, C_or_H, 
 * beta_lower, beta_upper, beta_step_small, beta_step_big, beta_1, beta_2,
 * accurate, choice (E or S or Q Or F), output file
 * J1, J2, J3 , K sample_amount,
 * ****/

#include <iostream>
#include <fstream> 
#include <string> 
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <iomanip> 
#include "/Users/Ke/Coding/My_library/dSFMT-src-2.2.3/dSFMT.h"

#ifndef SIZE
	#define SIZE 6
#endif

using namespace std; 

void build_gauge_bath();
void uniform_initialization();
void random_initialization();
void build_rotation_matrix(int i);
void updating();
double energy_bond_x(int bo, int nd);  
double energy_bond_y(int bo, int nd);  
double energy_bond_z(int bo, int nd); 
double energy_pla_xy(int i);
double energy_pla_yz(int i);
double energy_pla_zx(int i); 
void flip_R(int i);
void flip_Ux(int i);
void flip_Uy(int i);
void flip_Uz(int i);
void thermalization();
void estimate_beta_c(char *output);
void measure_energy(char *output);
double orderparameter();
void measure_J_distribution (char *output);
void measure_effective_J (char *output);

/*************  Constants     ***************/
const int L = SIZE;
const int L2 = L*L;
const int L3 = L*L*L;
const int U_order = 24;

/*********** Define fields **************/
double R[L3][9] = {{0}};
double Ux[L3][9] = {{0}};
double Uy[L3][9] = {{0}};
double Uz[L3][9] = {{0}};
double U_bath[U_order][9]={{0}};

double s[L3] = {0}; // s for Ising field

dsfmt_t dsfmt; //defined a required variable

/********** Measure ************/
double E_total = 0, E_bonds = 0, E_plaques = 0, E_g = 0; //E_g for the ground state energy
double E_change = 0;
//int Racc =0, xacc =0 , yacc = 0, zacc =0, Rrej=0, xrej=0, yrej=0, zrej=0;

/**** parameters ****/
double J1 , J2 , J3 ;
double K ;
double beta, beta_lower, beta_upper, beta_step_small, beta_step_big, beta_1, beta_2;
double accurate;
int tau = 100;
int sample_amount;

int main(int argc, char **argv)
{
	dsfmt_init_gen_rand(&dsfmt, time(0)); //seed dsfmt
	
	build_gauge_bath();
	
	J1 = -atof(argv[11]);
	J2 = -atof(argv[12]);
	J3 = -atof(argv[13]);
	
	K = -atof(argv[14]);
	
	sample_amount = atof(argv[15]);

	char *C_or_H = argv[1];
	
	beta_lower = atof(argv[2]);
	beta_upper = atof(argv[3]);
	
	if (*C_or_H == 'C') // if cooling
		{			
			beta = beta_lower; // beta from small to big when cooling
			beta_step_small = atof(argv[4]);
			beta_step_big = atof(argv[5]);

			random_initialization();
			
			/**** compute the energy of all bonds, checked ****/
			int xn, yn, zn;
			for (int i=0; i < L3; i++)
				{ 
				xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
				yn = (i + L) % L2 < L ? i + L - L2 : i + L;
				zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
		
				E_bonds += energy_bond_x(i,xn) + energy_bond_y(i,yn) + energy_bond_z(i, zn);
				}

			/**** compute the energy of all plaquette, checked ****/
			for (int i = 0; i < L3; i++)
				{ E_plaques += energy_pla_xy(i) + energy_pla_yz(i) + energy_pla_zx(i);}	
	
				E_g = 3*L3*(J1 + J2 + J3) + 3*L3*3*K; //E_g =  ground state energy of bonds and of plaques
	
				E_total = E_bonds + E_plaques;			
			}
	
	if (*C_or_H == 'H') // if heating
		{			
			beta = beta_upper;
			beta_step_small = -atof(argv[4]); // the step is negative since beta is decreasing
			beta_step_big = -atof(argv[5]);
			
			uniform_initialization();
			
			/**** compute the energy of all bonds, checked ****/
			int xn, yn, zn;
			for (int i=0; i < L3; i++)
				{ 
				xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
				yn = (i + L) % L2 < L ? i + L - L2 : i + L;
				zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
		
				E_bonds += energy_bond_x(i,xn) + energy_bond_y(i,yn) + energy_bond_z(i, zn);
				}

			/**** compute the energy of all plaquette, checked ****/
			for (int i = 0; i < L3; i++)
				{ E_plaques += energy_pla_xy(i) + energy_pla_yz(i) + energy_pla_zx(i);}
	
				E_g = 3*L3*(J1 + J2 + J3) + 3*L3*3*K; //E_g =  ground state energy of bonds and of plaques
	
				E_total = E_bonds + E_plaques;						
			}

		
	/**** initialize beta and the determine the fine region****/		
	beta_1 = atof(argv[6]);
	beta_2 = atof(argv[7]);
	
	accurate = atof(argv[8]);
	
	char *choice = argv[9];

	if(*choice == 'E')
		{ estimate_beta_c(argv[10]);}
	
	if(*choice == 'S')
		{ measure_energy(argv[10]);}
		
	if(*choice == 'D')
		measure_J_distribution(argv[10]);		

	if(*choice == 'J')
		measure_effective_J(argv[10]);		
	return 0;
	}


/**** build gauge bath ****/
void build_gauge_bath()
{
	/**** O checked ****/
	/** 0, Identity **/
	U_bath[0][0] = 1; 
	U_bath[0][4] = 1; 
	U_bath[0][8] = 1; 
	
	/** 1, C2,z **/
	U_bath[1][0] = -1; 
	U_bath[1][4] = -1; 
	U_bath[1][8] = 1;
	
	/** 2, C2,y **/
	U_bath[2][0] = -1; 
	U_bath[2][4] = 1; 
	U_bath[2][8] = -1;			
	
	/** 3, C2,x **/
	U_bath[3][0] = 1; 
	U_bath[3][4] = -1; 
	U_bath[3][8] = -1;
	
	/** 4, 3+ x,x,x **/
	U_bath[4][2] = 1; 
	U_bath[4][3] = 1; 
	U_bath[4][7] = 1;
	
	/** 5, 3+ -x,x,-x **/
	U_bath[5][2] = 1; 
	U_bath[5][3] = -1; 
	U_bath[5][7] = -1;			
	
	/** 6, 3+ x,-x,-x **/
	U_bath[6][2] = -1; 
	U_bath[6][3] = -1; 
	U_bath[6][7] = 1;			

	/** 7, 3+ -x, -x, x**/
	U_bath[7][2] = -1; 
	U_bath[7][3] = 1; 
	U_bath[7][7] = -1;
	
	/** 8, 3- x,x,x **/
	U_bath[8][1] = 1; 
	U_bath[8][5] = 1; 
	U_bath[8][6] = 1;			
	
	/** 9, 3- x,-x,-x **/
	U_bath[9][1] = -1; 
	U_bath[9][5] = 1; 
	U_bath[9][6] = -1;
	
	/** 10, 3- -x,-x,x **/
	U_bath[10][1] = 1; 
	U_bath[10][5] = -1; 
	U_bath[10][6] = -1;			
	
	/** 11, 3- -x,x-x **/
	U_bath[11][1] = -1; 
	U_bath[11][5] = -1; 
	U_bath[11][6] = 1;
	
	/** 12, 2 x,y=x,z=0 **/
	U_bath[12][1] = 1; 
	U_bath[12][3] = 1; 
	U_bath[12][8] = -1;	
	
	/** 13, 2 x,-x,0  **/
	U_bath[13][1] = -1; 
	U_bath[13][3] = -1; 
	U_bath[13][8] = -1;
	
	/** 14, 4-, 0,0,z  **/
	U_bath[14][1] = 1; 
	U_bath[14][3] = -1; 
	U_bath[14][8] = 1;		
	
	/** 15, 4+ 0,0,z **/
	U_bath[15][1] = -1; 
	U_bath[15][3] = 1; 
	U_bath[15][8] = 1;	

	/** 16, 4-, x,0,0 **/
	U_bath[16][0] = 1; 
	U_bath[16][5] = 1; 
	U_bath[16][7] = -1;
	
	/** 17, 2, 0,y,y **/
	U_bath[17][0] = -1; 
	U_bath[17][5] = 1; 
	U_bath[17][7] = 1;	
	
	/** 18, 2 0,y-y **/
	U_bath[18][0] = -1; 
	U_bath[18][5] = -1; 
	U_bath[18][7] = -1;	
	
	/** 19, 4+ x,0,0 **/
	U_bath[19][0] = 1; 
	U_bath[19][5] = -1; 
	U_bath[19][7] = 1;	
	
	/** 20, 4+ 0,y,0 **/
	U_bath[20][2] = 1; 
	U_bath[20][4] = 1; 
	U_bath[20][6] = -1;
	
	/** 21, 2x x,0,x **/
	U_bath[21][2] = 1; 
	U_bath[21][4] = -1; 
	U_bath[21][6] = 1;
	
	/** 22, 4- 0,y,0 **/
	U_bath[22][2] = -1; 
	U_bath[22][4] = 1; 
	U_bath[22][6] = 1;			
	
	/** 23, 2 -x,0,-x **/
	U_bath[23][2] = -1; 
	U_bath[23][4] = -1; 
	U_bath[23][6] = -1;												

/****
	for (int i = 0; i < U_order; i++)
		{	cout << "U_bath \t" << i << endl;
			for(int j = 0; j<9;j++)
				if ((j+1)%3==0) {cout << U_bath[i][j] << endl;}
				 else {cout << U_bath[i][j] << '\t';}
			}
****/
	
	}

/**** initialize R, sigma, U ****/
void uniform_initialization()
{
	for(int i = 0; i < L3; i++)
		{
		 R[i][0] = 1;
		 R[i][4] = 1;
		 R[i][8] = 1;
		 		 
		 s[i] = 1;
		 
		 Ux[i][0] = 1;
		 Ux[i][4] = 1;
		 Ux[i][8] = 1;
	 
		 Uy[i][0] = 1;
		 Uy[i][4] = 1;
		 Uy[i][8] = 1;
	 
		 Uz[i][0] = 1;
		 Uz[i][4] = 1;
		 Uz[i][8] = 1;		  	 	 
			}
	/****** checked ****/	
	
	}

void random_initialization()
{
	int i, j;
	for (i = 0; i < L3; i++)
		{
			build_rotation_matrix(i);
		
			/** build Ux **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
			copy(U_bath[j],U_bath[j]+9,Ux[i]);	
			
			/** build Uy **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
			copy(U_bath[j],U_bath[j]+9,Uy[i]);	
			
			/** build Uz **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
			copy(U_bath[j],U_bath[j]+9,Uz[i]);
			
			
			}
		
	}


void updating()
{
	int site, i;
	/**** one sweep ****/
	for (i = 0; i < L3*4*tau ; i++)
			 {
				/**** choose a site ****/
				site = int(L3*dsfmt_genrand_close_open(&dsfmt));

				/**** randomly flip R, Ux, Uy, Uz ****/
				switch(int(4 * dsfmt_genrand_close_open(&dsfmt))) 
					{  
						case 0 : flip_R(site); 
									break;
						case 1 : flip_Ux(site);
									break;
						case 2 : flip_Uy(site); 
									break;
						case 3 : flip_Uz(site); 
									break;					 
							}
					}
	}


/**** Change value of R[i] = SO(3) and sigma[i] = 1/-1 ****/
void build_rotation_matrix(int i) 
{	
	
	/************* Generating u1,u2,u3,u4 belons to [-1,1)
	 * and u1^2 + u2^2 =< 1, u3^2 + u4^2 <=1 *************/
	
	double u1, u2, u3, u4;
	
	u1 = -1 + 2 * dsfmt_genrand_close_open(&dsfmt);
	
	do {
		u2 = -1 + 2 * dsfmt_genrand_close_open(&dsfmt);
		} while (u1*u1 + u2*u2 >= 1);
		
	u3 = -1 + 2* dsfmt_genrand_close_open(&dsfmt);
		
	do {
		u4 = -1 + 2 * dsfmt_genrand_close_open(&dsfmt);
		} while (u3*u3 + u4*u4 >= 1);	 
	
	/***** Checking has been done for this part!****/
	
	/******** calculate x1,x2,x3,x4 for a unit 4D vector*************/
	double x1, x2, x3, x4, x_foo;
	
	x_foo = sqrt((1 - u1*u1 - u2*u2)/(u3*u3 + u4*u4));
	
	x1 = u1;
	x2 = u2;
	x3 = u3*x_foo;
	x4 = u4*x_foo;

	/********Checking for coordinates and the length has been done!*******/
	
	
	/****** Define the rotational matrix ******/
	
	
	double x12, x22, x32, x1x2, x3x4, x1x3, x2x4, x2x3, x1x4;
	
	x12 = x1*x1; x22 = x2*x2; x32 = x3*x3;
	x1x2 = x1*x2; x3x4 = x3*x4;
	x1x3 = x1*x3; x2x4 = x2*x4;
	x2x3 = x2*x3; x1x4 = x1*x4;
	
	
	R[i][0] = 1 - 2*(x22+x32); R[i][1] = 2*(x1x2-x3x4); R[i][2] = 2*(x1x3+x2x4);
	R[i][3] = 2*(x1x2+x3x4); R[i][4] = 1-2*(x12+x32); R[i][5] = 2*(x2x3-x1x4);
	R[i][6] = 2*(x1x3-x2x4); R[i][7] = 2*(x2x3+x1x4); R[i][8] = 1-2*(x12+x22);
				  
	/**** Othorgonality is verified up to 10^(-16) ********/
	
	/**** change value of sigma[i]****/
	s[i] = dsfmt_genrand_close_open(&dsfmt) > 0.5 ? 1 : -1;	
	}
	
	
/**** i passed from main so not need to defined int again and again
 * s[i] also changed in build rotation
 * ****/

/**** bond energy of Bond = s[i] s[j] Tr[ R^T[i] Ux[i] R[j] ] ****/
double energy_bond_x(int bo, int nd)
{
	double bond_energy;
	double Bond[9] = {0}, foo[9] = {0};
	
	/** i = bo, the starting of the bond, j = nd, the end of the bond 
	 * Bond = s[bo] s[nd] R[bo]^T Ux[bo] R[nd]
	 * 		= s[bo] s[nd] Ux[bo] R[nd] R^T[bo] //this is for the convience of mutiplying J	
	 * foo = s[bo] s[nd] R[nd] R^T[bo]
	 * Bond = Ux[bo] foo
	 * **/ 
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[bo]*s[nd],
	R[nd], 3, R[bo],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[bo], 3, foo,3,
	0.0, Bond,3);

	bond_energy = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	return bond_energy;
	
	}


/**** bond energy of Bond = s[i] s[j] Tr[ R^T[i] Uy[i] R[j] ] ****/
double energy_bond_y(int bo, int nd)
{
	double bond_energy;
	double Bond[9] = {0}, foo[9] = {0};
	
	/** i = bo, the starting of the bond, j = nd, the end of the bond 
	 * Bond = s[bo] s[nd] R[bo]^T Uy[bo] R[nd]
	 * 		= s[bo] s[nd] Uy[bo] R[nd] R^T[bo] //this is for the convience of mutiplying J	
	 * foo = s[bo] s[nd] R[nd] R^T[bo]
	 * Bond = Uy[bo] foo
	 * **/ 
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[bo]*s[nd],
	R[nd], 3, R[bo],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[bo], 3, foo,3,
	0.0, Bond,3);

	bond_energy = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	return bond_energy;
	
	}
	
	
/**** bond energy of Bond = s[i] s[j] Tr[ R^T[i] Uz[i] R[j] ] ****/
double energy_bond_z(int bo, int nd)
{
	double bond_energy;
	double Bond[9] = {0}, foo[9] = {0};
	
	/** i = bo, the starting of the bond, j = nd, the end of the bond 
	 * Bond = s[bo] s[nd] R[bo]^T Uz[bo] R[nd]
	 * 		= s[bo] s[nd] Uz[bo] R[nd] R^T[bo] //this is for the convience of mutiplying J	
	 * foo = s[bo] s[nd] R[nd] R^T[bo]
	 * Bond = Uz[bo] foo
	 * **/ 
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[bo]*s[nd],
	R[nd], 3, R[bo],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[bo], 3, foo,3,
	0.0, Bond,3);

	bond_energy = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	return bond_energy;
	
	}	
	
/**** energy of a unit plaquette in the xy plane****/
double energy_pla_xy(int i)
{
	/** find the required two neighbours**/
	int xn, yn;
	xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
	yn = (i + L) % L2 < L ? i + L - L2 : i + L;
	
	/** compute the matrix of the plaquetee
	 * Pla = Ux[i] Uy[xn] Ux[yn]^T Uy[i]^T
	 * foo1 =  Ux[yn]^T Uy[i]^T
	 * foo2 = Uy[xn] foo1
	 * Pla = Ux[i] foo2
	 * **/
	double Pla[9] = {0}, foo1[9] = {0}, foo2[9] = {0};
	
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	3,3,3,1,
	Ux[yn], 3, Uy[i],3,
	0.0, foo1,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[xn], 3, foo1,3,
	0.0, foo2,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[i], 3, foo2,3,
	0.0, Pla,3);	
	
	/**** plaquett energy ****/
	double pla = Pla[0] + Pla[4] + Pla[8];
	
	return K*pla;
	
	}	

/**** energy of a unit plaquette in the yz plane****/
double energy_pla_yz(int i)
{
	/** find the required two neighbours**/
	int yn, zn;
	yn = (i + L) % L2 < L ? i + L - L2 : i + L;
	zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
	
	/** compute the matrix of the plaquetee
	 * Pla = Uy[i] Uz[yn] Uy[zn]^T Uz[i]^T
	 * foo1 =  Uy[zn]^T Uz[i]^T
	 * foo2 = Uz[yn] foo1
	 * Pla = Uy[i] foo2
	 * **/
	double Pla[9] = {0}, foo1[9] = {0}, foo2[9] = {0};
	
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	3,3,3,1,
	Uy[zn], 3, Uz[i],3,
	0.0, foo1,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[yn], 3, foo1,3,
	0.0, foo2,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[i], 3, foo2,3,
	0.0, Pla,3);	
	
	/**** plaquett energy ****/
	double pla = Pla[0] + Pla[4] + Pla[8];
	
	return K*pla;
	
	}	
	
/**** energy of a unit plaquette in the zx plane****/
double energy_pla_zx(int i)
{
	/** find the required two neighbours**/
	int zn, xn;
	zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
	xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
		
	/** compute the matrix of the plaquetee
	 * Pla = Uz[i] Ux[zn] Uz[xn]^T Ux[i]^T
	 * foo1 =  Uz[xn]^T Ux[i]^T
	 * foo2 = Ux[zn] foo1
	 * Pla = Uz[i] foo2
	 * **/
	double Pla[9] = {0}, foo1[9] = {0}, foo2[9] = {0};
	
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	3,3,3,1,
	Uz[xn], 3, Ux[i],3,
	0.0, foo1,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[zn], 3, foo1,3,
	0.0, foo2,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[i], 3, foo2,3,
	0.0, Pla,3);	
	
	/**** plaquett energy ****/
	double pla = Pla[0] + Pla[4] + Pla[8];
	
	return K*pla;
	
	}		

void flip_R(int i) 
{
	/****** find the 6 neighbours, checked*****/
	
	int xp, xn, yp, yn, zp, zn; 
	
	xp = i % L == 0 ? i - 1 + L : i - 1;
	xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
	
	yp = i % L2 < L ? i - L + L2 : i - L;
	yn = (i + L) % L2 < L ? i + L - L2 : i + L;
	
	zp = i < L2 ? i - L2 + L3 : i - L2;
	zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
	
	/**** energy of the six connected bonds****/
	double bond_xp_i, bond_i_xn;
	double bond_yp_i, bond_i_yn;
	double bond_zp_i, bond_i_zn;
	
	bond_xp_i = energy_bond_x(xp,i); bond_i_xn = energy_bond_x(i,xn);
	bond_yp_i = energy_bond_y(yp,i); bond_i_yn = energy_bond_y(i,yn);
	bond_zp_i = energy_bond_z(zp,i); bond_i_zn = energy_bond_z(i,zn);
	
	double E_old;
	E_old = bond_xp_i + bond_i_xn + bond_yp_i + bond_i_yn
			+ bond_zp_i + bond_i_zn;
	
	double s_save;
	double R_save[9];
	
	copy(R[i], R[i]+9,R_save); //save R[i] to R_save and s[i] to s_save
	s_save = s[i];
	
	build_rotation_matrix(i); //generate new R[i] and s[i]
	
	/**** recompute the energy of the six connected bonds, and save it to E_new****/
	bond_xp_i = energy_bond_x(xp,i); bond_i_xn = energy_bond_x(i,xn);
	bond_yp_i = energy_bond_y(yp,i); bond_i_yn = energy_bond_y(i,yn);
	bond_zp_i = energy_bond_z(zp,i); bond_i_zn = energy_bond_z(i,zn);
	
	double E_new;
	
	E_new = bond_xp_i + bond_i_xn + bond_yp_i + bond_i_yn
			+ bond_zp_i + bond_i_zn;
	
	E_change = E_new - E_old;
	
	/*******decide flip and change E_total********/
	if (E_change < 0)
		{E_total += E_change;}
	 else
		{ if (exp(-beta * E_change) > dsfmt_genrand_close_open(&dsfmt))
			{E_total += E_change;}
		   else	{copy(R_save, R_save+9, R[i]); 
				 s[i] = s_save;
			     } // change R[i] and s[i] back				
			}	
	}

void flip_Ux(int i) // the connected Bond is R^T[i] Ux[i] R[xn] 
{
	/** energy of the related bond **/
	int xn; // the end of the bond
	double E_b; // energy of the bond

	xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
	E_b = energy_bond_x(i, xn);
	
	/** energy of the related plaquetees**/
	int zp, yp; // orign of neighbouring plaquettes
	double E_pla; // energy of the four related plaquettes

	zp = i < L2 ? i - L2 + L3 : i - L2;	
	yp = i % L2 < L ? i - L + L2 : i - L;

	E_pla = energy_pla_zx(i) + energy_pla_zx(zp) + energy_pla_xy(i) + energy_pla_xy(yp);
	//E_pla /= beta; // cancel the beta-dependence of the defect penalty	
	
	double E_old;
	E_old = E_b + E_pla;
	
	/**** save Ux[i] to U_save ****/
	double U_save[9];
	copy(Ux[i],Ux[i]+9,U_save);
	
	/**** generate new Ux by choosing from U_bath ****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(U_bath[j],U_bath[j]+9,Ux[i]);
	
	/**** compute the bond enery again and store it to E_enw****/
	E_b = energy_bond_x(i, xn);
	E_pla = energy_pla_zx(i) + energy_pla_zx(zp) + energy_pla_xy(i) + energy_pla_xy(yp);
	//E_pla /= beta; //cancel the beta-dependence of the defect penalty		

	double E_new;
	E_new = E_b + E_pla;
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
	if (E_change < 0)
		{E_total += E_change;}
	 else { if (exp(-beta * E_change) > dsfmt_genrand_close_open(&dsfmt))
				{E_total += E_change;}
			 else {copy(U_save,U_save+9,Ux[i]);}	
			}		
	
	}

void flip_Uy(int i)
{
	/** energy of the related bond **/
	int yn; // the end of the bond
	double E_b; // energy of the bond
	
	yn = (i + L) % L2 < L ? i + L - L2 : i + L;
	
	E_b = energy_bond_y(i,yn);	
	
	/** energy of the related plaquetees**/
	int xp, zp; // orign of neighbouring plaquettes
	double E_pla; // energy of the four related plaquettes
	
	xp = i % L == 0 ? i - 1 + L : i - 1;
	zp = i < L2 ? i - L2 + L3 : i - L2;
	
	E_pla = energy_pla_xy(i) + energy_pla_xy(xp) + energy_pla_yz(i) + energy_pla_yz(zp);
	//E_pla /= beta;// cancel the beta-dependence of the defect penalty
	
	double E_old;
	E_old = E_b + E_pla;
	
	/**** save Uy[i] to U_save ****/
	double U_save[9];
	copy(Uy[i],Uy[i]+9,U_save);
	
	/**** generate choosing new Uy from U_bath****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(U_bath[j],U_bath[j]+9,Uy[i]);
	
	/**** compute the bond enery again and store it to E_new****/
	E_b = energy_bond_y(i,yn);
	E_pla = energy_pla_xy(i) + energy_pla_xy(xp) + energy_pla_yz(i) + energy_pla_yz(zp);
	//E_pla /= beta; // cancel the beta-dependence of the defect penalty
	
	double E_new;
	E_new = E_b + E_pla;
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
	if (E_change < 0)
		{E_total += E_change;}
	 else { if (exp(-beta * E_change) > dsfmt_genrand_close_open(&dsfmt))
				{E_total += E_change;}
			 else {copy(U_save,U_save+9,Uy[i]);}	
			}		
	
	}


void flip_Uz(int i)
{
	/** energy of the related bond **/
	int zn; // the end of the bond
	double E_b; // energy of the bond

	zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
	E_b = energy_bond_z(i, zn);
	
	/** energy of the related plaquetees**/
	int yp, xp; // orign of neighbouring plaquettes
	double E_pla; // energy of the four related plaquettes
	
	yp = i % L2 < L ? i - L + L2 : i - L;
	xp = i % L == 0 ? i - 1 + L : i - 1;

	E_pla = energy_pla_yz(i) + energy_pla_yz(yp) + energy_pla_zx(i) + energy_pla_zx(xp);
	//E_pla /= beta; // cancel the beta-dependence of the defect penalty	
	
	double E_old;
	E_old = E_b + E_pla;
	
	/**** save Uz[i] to U_save ****/
	double U_save[9];
	copy(Uz[i],Uz[i]+9,U_save);
	
	/**** generate new Uz by choosing new U from U_bath****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(U_bath[j],U_bath[j]+9,Uz[i]);
	
	/**** compute the bond enery again and store it to E_new****/
	E_b = energy_bond_z(i, zn);	
	E_pla = energy_pla_yz(i) + energy_pla_yz(yp) + energy_pla_zx(i) + energy_pla_zx(xp);
	//E_pla /= beta; // cancel the beta-dependence of the defect penalty
		
	double E_new;
	E_new = E_b + E_pla;
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
	if (E_change < 0)
		{E_total += E_change;}
	 else { if (exp(-beta * E_change) > dsfmt_genrand_close_open(&dsfmt))
				{E_total += E_change;}
			 else {copy(U_save,U_save+9,Uz[i]);}	
			}		
	}

void thermalization()
{
	double afoo = 0, s1 = 1, s2;

//int steps = 0;	
	
	while (afoo < 1-accurate || afoo > 1+accurate)
		{
			s2 = s1; // copy old value of s1 to s2.
			
			/**** updata configurations and calculate new s1 ****/
			for (int i = 0; i<1000; i++)
				{ s1 += E_total;
					updating();						
					}		
			afoo = s1/s2;
//steps++;		
			}
//cout << beta << '\t' << steps << "\t *1000 T" << endl;		
	}
	
void estimate_beta_c(char *output)
{
	double S1, S2,Cv;
//	double foo2, Q1, Q2, chi_Q; 
	double foo_s, s1, s2, chi_s;
	
	ofstream output_file;
	output_file.open(output);
	
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
		{ 
/** re-initialize quantites for the acception ratio **/
//Racc = 0; Rrej = 0; xacc = 0; xrej = 0;
//yacc = 0; yrej = 0; zacc = 0; zrej = 0;
			
		/**** re-thermalization and reset S1, S2****/
		  thermalization();
		  S1 = 0;
		  S2 = 0;
		  
		  s1 = 0;
		  s2 = 0;
		  
//		  Q1 = 0;		  Q2 = 0;
		  /**** measure ****/	  
		  for (int i = 0; i < sample_amount; i++)
		  { 
			updating();	
			S1 += E_total;	 //extensive here
			S2 += E_total * E_total;	// extensive^2
			
//			foo2 = orderparameter();
//			Q2 += foo2;
//			Q1 += sqrt(foo2);
			
			foo_s = 0;
			for(int k = 0; k < L3; k++) 
				{ foo_s += s[k];} //extensive
			foo_s /= L3; // intensive
			s1 += foo_s;
			s2 += foo_s*foo_s;
			 }	 
			
		S1 /= sample_amount;
		S2 /= sample_amount;
		Cv = (S2 - S1 * S1) * beta * beta / L3;	 
		
		S1 /= E_g; 
		
//		Q1/=sample_amount;
//		Q2 /= sample_amount;
//		chi_Q = (Q2 - Q1*Q1)*beta*L3;
		
		s1 /= sample_amount;
		s2 /= sample_amount;
		chi_s = (s2 -s1*s1)*beta*L3;

		output_file << 1/beta << '\t'<< S1 << '\t' << Cv << '\t' << s1<< '\t' 
								<< chi_s << endl;
								//<< Q1 << '\t' << chi_Q << endl;
		
		if ( (beta >= beta_1) && (beta <= beta_2) )
			{beta += beta_step_small;}
		else {beta += beta_step_big;}		
		
/** acception ratio**/		
//cout << "R, Ux, Uy, Uz" << endl;
//cout << beta << '\t' << Racc << '\t' << Rrej << '\t' << Racc*1.0/(Racc+Rrej) << endl; 
//cout << beta << '\t' << xacc << '\t' << xrej << '\t' << xacc*1.0/(xacc+xrej) << endl;
//cout << beta << '\t' << yacc << '\t' << yrej << '\t' << yacc*1.0/(yacc+yrej) << endl;  
//cout << beta << '\t' << zacc << '\t' << zrej << '\t' << zacc*1.0/(zacc+zrej) << endl; 		 		
			}
	output_file.close();
	}

void measure_energy(char *output)
{

	ofstream output_file;
	output_file.open(output);
	
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
		{ 

/** re-initialize quantites for the acception ratio **/
//Racc = 0; Rrej = 0; xacc = 0; xrej = 0;
//yacc = 0; yrej = 0; zacc = 0; zrej = 0;

			
		/**** re-thermalization ****/
		  thermalization();

		  
		 output_file << beta << '\t' << J1 << '\t' << J2 << '\t' << J3 << '\t'<<flush;

		  /**** measure ****/	  
		  for ( int i = 0; i < sample_amount; i++)
		  { 
			updating();
			output_file << E_total/L3 << '\t' << flush; 
			 }
			output_file << endl; 	 
			
					
		if ( (beta >= beta_1) && (beta <= beta_2) )
			{beta += beta_step_small;}
		else {beta += beta_step_big;}	
/** acception ratio**/		
//cout << "R, Ux, Uy, Uz" << endl;
//cout << beta << '\t' << Racc << '\t' << Rrej << '\t' << Racc*1.0/(Racc+Rrej) << endl; 
//cout << beta << '\t' << xacc << '\t' << xrej << '\t' << xacc*1.0/(xacc+xrej) << endl;
//cout << beta << '\t' << yacc << '\t' << yrej << '\t' << yacc*1.0/(yacc+yrej) << endl;  
//cout << beta << '\t' << zacc << '\t' << zrej << '\t' << zacc*1.0/(zacc+zrej) << endl; 		 	
			}
	output_file.close();
	}

/**** the order parameter for one sample****/

double orderparameter()
{
	/**** creat the T_abc tensor first; extensive at T, intensive at return Q2****/
	
	double T[3][3][3] = {{{0}}}; 
	
	/** compute every element of T, over the lattice
	 * T_abc = s[i] (l_a m_b n_c - 1/6 empsilon_abc)
	 * l_a = R[0], R[1], R[2]
	 * m_b = R[3], R[4], R[5]
	 * n_c = R[6], R[7], R[8]
	 * 27 terms in total
	 * **/
	 
	/**** l=0,m=0 ****/ 
	for(int i = 0; i < L3; i++)
		{ T[0][0][0] += s[i]*R[i][0]*R[i][3]*R[i][6];}	
	for(int i = 0; i < L3; i++)
		{ T[0][0][1] += s[i]*R[i][0]*R[i][3]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[0][0][2] += s[i]*R[i][0]*R[i][3]*R[i][8];}
	
	/**** l=0,m=1 ****/
	for(int i = 0; i < L3; i++)
		{ T[0][1][0] += s[i]*R[i][0]*R[i][4]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[0][1][1] += s[i]*R[i][0]*R[i][4]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[0][1][2] += s[i]*R[i][0]*R[i][4]*R[i][8] - s[i]*0.1666666666667;} // 1/6
	
	/**** l=0,m=2 ****/
	for(int i = 0; i < L3; i++)
		{ T[0][2][0] +=  s[i]*R[i][0]*R[i][5]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[0][2][1] +=  s[i]*R[i][0]*R[i][5]*R[i][7] + s[i]*0.1666666666667;}
	for(int i = 0; i < L3; i++)
		{ T[0][2][2] += s[i]*R[i][0]*R[i][5]*R[i][8];}
	
	/**** l=1,m=0 ****/ 
	for(int i = 0; i < L3; i++)
		{ T[1][0][0] +=  s[i]*R[i][1]*R[i][3]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[1][0][1] +=  s[i]*R[i][1]*R[i][3]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[1][0][2] +=  s[i]*R[i][1]*R[i][3]*R[i][8] + s[i]*0.1666666666667;}
	
	/**** l=1,m=1 ****/
	for(int i = 0; i < L3; i++)
		{ T[1][1][0] +=  s[i]*R[i][1]*R[i][4]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[1][1][1] +=  s[i]*R[i][1]*R[i][4]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[1][1][2] += s[i]*R[i][1]*R[i][4]*R[i][8];}
	
	/**** l=1,m=2 ****/
	for(int i = 0; i < L3; i++)
		{ T[1][2][0] +=  s[i]*R[i][1]*R[i][5]*R[i][6] - s[i]*0.1666666666667;}
	for(int i = 0; i < L3; i++)
		{ T[1][2][1] += s[i]*R[i][1]*R[i][5]*R[i][7] ;}
	for(int i = 0; i < L3; i++)
		{ T[1][2][2] += s[i]*R[i][1]*R[i][5]*R[i][8];}	
	
		/**** l=2,m=0 ****/ 
	for(int i = 0; i < L3; i++)
		{ T[2][0][0] +=  s[i]*R[i][2]*R[i][3]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[2][0][1] += s[i]*R[i][2]*R[i][3]*R[i][7] - s[i]*0.1666666666667;}	
	for(int i = 0; i < L3; i++)
		{ T[2][0][2] +=  s[i]*R[i][2]*R[i][3]*R[i][8];}
	
	/**** l=2,m=1 ****/
	for(int i = 0; i < L3; i++)
		{ T[2][1][0] += s[i]*R[i][2]*R[i][4]*R[i][6] + s[i]*0.1666666666667;}
	for(int i = 0; i < L3; i++)
		{ T[2][1][1] +=  s[i]*R[i][2]*R[i][4]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[2][1][2] +=  s[i]*R[i][2]*R[i][4]*R[i][8];}
	
	/**** l=2,m=2 ****/
	for(int i = 0; i < L3; i++)
		{ T[2][2][0] += s[i]*R[i][2]*R[i][5]*R[i][6];}
	for(int i = 0; i < L3; i++)
		{ T[2][2][1] +=  s[i]*R[i][2]*R[i][5]*R[i][7];}
	for(int i = 0; i < L3; i++)
		{ T[2][2][2] +=  s[i]*R[i][2]*R[i][5]*R[i][8];}
	 
	 /**** define the Q_Td tensor ****/
	 double Q[3][3][3] = {{{0}}};
	 
	 /**** with all same index ****/
	 Q[0][0][0] = 6*T[0][0][0];
	 Q[1][1][1] = 6*T[1][1][1];
	 Q[2][2][2] = 6*T[2][2][2];
	 
	 /**** with all different index ****/
	 Q[0][1][2] = 2*T[0][1][2] +  2*T[1][2][0] + 2*T[2][0][1] ;
	 Q[1][2][0] = Q[0][1][2];
	 Q[2][0][1] = Q[0][1][2];

	 Q[0][2][1] =  2*T[0][2][1]  + 2*T[1][0][2] + 2*T[2][1][0];
	 Q[1][0][2] = Q[0][2][1];
	 Q[2][1][0] = Q[0][2][1];

	 /**** with two same index ****/
	 Q[0][0][1] = 2*T[0][0][1] + 2*T[0][1][0] + 2*T[1][0][0];
	 Q[0][1][0] = Q[0][0][1];
	 Q[1][0][0] = Q[0][0][1];
	 
	 Q[0][0][2] = 2*T[0][0][2] + 2*T[0][2][0] + 2*T[2][0][0];
	 Q[0][2][0] = Q[0][0][2];
	 Q[2][0][0] = Q[0][0][2];
	 
	 Q[1][1][0] = 2*T[1][1][0] + 2*T[1][0][1] + 2*T[0][1][1];
	 Q[1][0][1] = Q[1][1][0];
	 Q[0][1][1] = Q[1][1][0];
	 
	 Q[1][1][2] = 2*T[1][1][2] + 2*T[1][2][1] + 2*T[2][1][1];
	 Q[1][2][1] = Q[1][1][2];
	 Q[2][1][1] = Q[1][1][2];
	 
	 Q[2][2][0] = 2*T[2][2][0] + 2*T[2][0][2] + 2*T[0][2][2];
	 Q[2][0][2] = Q[2][2][0];
	 Q[0][2][2] = Q[2][2][0];
	 
	 Q[2][2][1] = 2*T[2][2][1] + 2*T[2][1][2] + 2*T[1][2][2];
	 Q[2][1][2] = Q[2][2][1];
	 Q[1][2][2] = Q[2][2][1];		
	 
	 double Q2 = 0;

	 for(int a = 0; a < 3; a++)
		{for(int b = 0; b<3; b++)
			{for(int c=0; c<3; c++)
					{ Q2 += Q[a][b][c]*Q[a][b][c];}
			}
		}	 
		
	Q2 /= L3*L3;
		
	 return Q2;
	}

void measure_J_distribution(char *output)
{	
	/**
	 * Computing J_eff = Tr(R^T_i J U_{ij} R_j) = energy_bond(bo, nd) /s[bo]/s[nd]
	 *  R is defened no s, the function energy bond is defined with s
	 * J_eff[3*L3] are defined to save the value of J_eff of the 3N bonds,
	 * numbers of samples use the sample_amount in the header
	 * cout -J_eff in output, since J is defined with a minus sign
	 * **/
	int xn, yn, zn; 
	double bond; //bond energy
	int ss; // s[i]s[j]
	
	ofstream output_file;
	output_file.open(output);
	
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
		{ 

		/**** re-thermalization and reset J_eff****/
		  thermalization();
		  
		double J_eff[3*L3]={0}; // define J_eff at the current temperature 
		  
		  /**** measure ****/	  
		  for (int k = 0; k < sample_amount; k++)
		  { 
			updating();
					
			/** sweep the lattice **/
			for( int i = 0; i < L3; i++)
				{
					/** Jeff at x bond **/

					xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
					bond = energy_bond_x(i, xn);
					ss = s[i]*s[xn];					
					J_eff[3*i] += bond/ss;
					
					/** J_eff at y direction **/
					
					yn = (i + L) % L2 < L ? i + L - L2 : i + L;
					bond = energy_bond_y(i, yn);
					ss = s[i]*s[yn];					
					J_eff[3*i+1] += bond/ss;
					
					
					/** J_eff at z direction **/
					
					zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
					bond = energy_bond_z(i, zn);
					ss = s[i]*s[zn];					
					J_eff[3*i+2] += bond/ss;					
					
					}
					
			 }	 
		
		output_file << beta << '\t';	
		for (int j = 0; j < 3*L3; j++)
				output_file << -J_eff[j]/sample_amount << '\t';				

		output_file << endl;
		
		if ( (beta >= beta_1) && (beta <= beta_2) )
			{beta += beta_step_small;}
		else {beta += beta_step_big;}	
				 		
			}
	output_file.close();
	}

void measure_effective_J (char *output)
{	
	/**
	 * Computing J_eff = <Tr(R^T_i J U_{ij} R_j)> = <bond/(s[i]*s[j])>
	 *  R is defened no s, bond is computed with s
	 * cout -J_eff in output, since J is defined with a minus sign
	 * **/
	int xn, yn, zn; 
	int ss; // ss = s[i]s[j]
	double bond; //bond energy
	double J_eff; //different to J_eff in measure_J_distribution, J_eff here is a number, the average 
	
	ofstream output_file;
	output_file.open(output);
	
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
		{ 

		/**** re-thermalization and reset J_eff****/
		  thermalization();
		  
		J_eff = 0; // re-set J_eff at the current temperature 
		  
		  /**** measure ****/	  
		  for (int k = 0; k < sample_amount; k++)
		  { 
			updating();
					
			/** sweep the lattice **/
			for( int i = 0; i < L3; i++)
				{
					/**** J_eff at x direction ****/
					xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
					bond = energy_bond_x(i, xn);
					ss = s[i]*s[xn];		
					J_eff += bond/ss;
					
					/**** J_eff at y direction ****/					
					yn = (i + L) % L2 < L ? i + L - L2 : i + L;
					bond = energy_bond_y(i, yn);
					ss = s[i]*s[yn];
					J_eff += bond/ss;
					
					/**** J_eff at z direction ****/
					
					zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
					bond = energy_bond_z(i, zn);
					ss = s[i]*s[zn];
					J_eff += bond/ss;					
					
					}
					
			 }	 
		J_eff = J_eff/sample_amount/L3/3; // the additional 3 is due to there are 3L3 bonds
		
		output_file << beta << '\t' << -J_eff<< endl;

		if ( (beta >= beta_1) && (beta <= beta_2) )
			{beta += beta_step_small;}
		else {beta += beta_step_big;}	
				 		
			}
	output_file.close();
	}
