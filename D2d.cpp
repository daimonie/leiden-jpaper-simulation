/****
 * argv = filename, C_or_H, 
 * beta_lower, beta_upper, beta_step_small, beta_step_big, beta_1, beta_2,
 * accurate, choice (E or S or Q Or F), output file
 * J1, J2, J3 ,sample_amount
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
#include "dSFMT-src-2.2.3/dSFMT.h"
#include <ctime>
#include "omp.h"
#include <random> 

#ifndef SIZE
	#define SIZE 8
#endif

using namespace std; 

void build_gauge_bath();
void uniform_initialization();
void random_initialization();
void build_rotation_matrix(int i);
double site_energy(int i);  
void flip_R(int i);
void flip_Ux(int i);
void flip_Uy(int i);
void flip_Uz(int i);
void thermalization();
void estimate_beta_c();
void measure_energy(char *output);
double orderparameter_n();
void measure_J_distribution(char *output);
void measure_effective_J (char *output);
void generate_rotation_matrices();

void josko_diagnostics ();

/*************  Constants     ***************/
const int L = SIZE;
const int L2 = L*L;
const int L3 = L*L*L;
const int U_order = 8;

/*********** Define fields **************/
double R[L3][9] = {{0}};
double Ux[L3][9] = {{0}};
double Uy[L3][9] = {{0}};
double Uz[L3][9] = {{0}};
double U_bath[U_order][9]={{0}};

double s[L3] = {0}; // s for Ising field

dsfmt_t dsfmt; //defined a required variable

/********** Measure ************/
double E_total, E_change, E_g;
//int Racc =0, xacc =0 , yacc = 0, zacc =0, Rrej=0, xrej=0, yrej=0, zrej=0;

/**** parameters ****/
double J1 = 1, J2 = 1, J3 = 1;
double beta, beta_lower, beta_upper, beta_step_small, beta_step_big, beta_1, beta_2;
double accurate;
int tau = 100;
int sample_amount;
/**** rotation matrices cache ****/
const int rmc_number = 10; // Note: it generates rmc_number^3 matrices.
const int rmc_number_total = rmc_number*rmc_number*rmc_number;
double rmc_matrices[rmc_number_total][9] = {{0}};

/**** time to start the program ****/

int main(int argc, char **argv)
{
	time_t tstart, tend;
	tstart = time(0);
	dsfmt_init_gen_rand(&dsfmt, time(0)); //seed dsfmt 
        
        omp_set_num_threads(4); //I still want to use my computer :)
        generate_rotation_matrices(); 
        
	build_gauge_bath();
	
	J1 = -atof(argv[11]);
	J2 = -atof(argv[12]);
	J3 = -atof(argv[13]);
		

	sample_amount = atof(argv[14]);

	char *C_or_H = argv[1];
	
	beta_lower = atof(argv[2]);
	beta_upper = atof(argv[3]);
	
	if (*C_or_H == 'C') // if cooling
		{			
		beta = beta_lower; // beta from small to big when cooling
		beta_step_small = atof(argv[4]);
		beta_step_big = atof(argv[5]);

		random_initialization();

		for (int i=0; i < L3; i++)
		{
			E_total += site_energy(i);
		}

		E_total /= 2;
		E_g = L3*3*(J1+J2+J3);				
	}
	else if (*C_or_H == 'H') // if heating
	{			
		beta = beta_upper;
		beta_step_small = -atof(argv[4]); // the step is negative since beta is decreasing
		beta_step_big = -atof(argv[5]);

		uniform_initialization();

		for (int i=0; i < L3; i++)
		{
			E_total += site_energy(i);
		}

		E_total /= 2;
		E_g = L3*3*(J1+J2+J3);							
	}

		
	/**** initialize beta and the determine the fine region****/		
	beta_1 = atof(argv[6]);
	beta_2 = atof(argv[7]);
	
	accurate = atof(argv[8]);
	
	char *choice = argv[9];

        josko_diagnostics ();
	printf("#//Calculate from %2.3f to %2.3f, using {%2.3f, %2.3f, %2.3f}, accuracy %2.3f and %d samples \n", beta_lower, beta_upper, J1, J2, J3, accurate, sample_amount);
	printf("#//Maximum cores %d \n", omp_get_max_threads());
        if(*choice == 'E')
        {
               estimate_beta_c();
        }
        else if(*choice == 'S')
        {
                measure_energy(argv[10]);
                
        }
        else if(*choice == 'D')
        {
                measure_J_distribution(argv[10]);		
        }
        else if(*choice == 'J')
        {
                measure_effective_J(argv[10]);							
        }
        
	tend = time(0);
	
	printf("#// Time taken is %2.3f seconds", difftime(tend,tstart));

	return 0;
}
/**** Josko's Diagnostics ****/
void josko_diagnostics()
{
        int ii, jj, kk;
        int rmc_samples = rmc_number*rmc_number*rmc_number*10000;
        int build_random = 0;
        double rmc_matrix_average[9] = {0}; 
        for(ii = 0; ii < rmc_samples; ii++)
        {
                build_random = (int) (rmc_number_total * dsfmt_genrand_close_open(&dsfmt));
                
                for(jj = 0; jj < 9; jj++)
                {      
                        rmc_matrix_average[jj] += rmc_matrices[build_random][jj]/rmc_samples;
                }
                
        } 
        
        printf("Average of sampled matrices{ %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f) \n",  rmc_matrix_average[0], rmc_matrix_average[1], rmc_matrix_average[2],
                rmc_matrix_average[3], rmc_matrix_average[4], rmc_matrix_average[5],
                rmc_matrix_average[6], rmc_matrix_average[7], rmc_matrix_average[8]); 
	
        int rnd_samples = 100;
        double rnd_average = 0.0;
        double rnd_correlation = 0.0;
        double rnd_list[rnd_samples]; 
        for(ii = 0; ii < rnd_samples; ii++)
        {
                rnd_list[ii] = dsfmt_genrand_close_open(&dsfmt);
                rnd_average += rnd_list[ii] / rnd_samples;
        } 
        for(jj =0; jj < rnd_samples; jj++)
        {
                for(kk = 0; kk < rnd_samples; kk++)
                {   
                        rnd_correlation += rnd_list[jj] * rnd_list[ (jj+kk)%rnd_samples ] / rnd_samples;
                }
                rnd_correlation -= rnd_average; 
        } 
        rnd_correlation /= rnd_samples;
        printf("Test 4.1, average [%2.3f], correlation [%.3f]. \n",rnd_average, rnd_correlation);
	// Different engine
	rnd_average = 0.0;
	rnd_correlation = 0.0;
	
	printf("New random generator. \n");
	
	int rnd_samples2 = rnd_samples*rnd_samples;
        double rnd_list2[rnd_samples2]; 	  
	int j = 0;   
	
	std::mt19937_64 engine(ii);
	std::uniform_real_distribution<double> dice(0.0, 1.0);
		
		
	#pragma omp parallel for private(j, engine, dice)
        for(ii = 0; ii < rnd_samples; ii++)
        {
		engine = std::mt19937_64(ii);
		dice = uniform_real_distribution<double>(0.0, 1.0);
		for( jj = 0; jj < rnd_samples; jj++)
		{
			j = ii * rnd_samples + jj;
			
			rnd_list2[j] = dice(engine);
			rnd_average += rnd_list2[j] / rnd_samples2;
		}
        }    
	#pragma omp parallel for
        for(jj =0; jj < rnd_samples2; jj++)
        {
                for(kk = 0; kk < rnd_samples2; kk++)
                {   
                        rnd_correlation += rnd_list2[jj] * rnd_list2[ (jj+kk)%rnd_samples2 ] / rnd_samples2;
                }
                rnd_correlation -= rnd_average; 
        } 
        rnd_correlation /= rnd_samples2;
        printf("New engine, average [%2.3f], correlation [%.3f]. \n",rnd_average, rnd_correlation);
}
/**** generates rotation matrices ****/
void generate_rotation_matrices ()
{
        //printf("Generating %d rotation matrices.. \n", rmc_number);

        int ii, jj, kk;
        
        int rmc_index = 0;
        
        double w,x,y,z, u1, u2, u3;  
        
        for (ii = 0; ii < rmc_number; ii++)
        {
                for (jj = 0; jj < rmc_number; jj++)
                {
                        for (kk = 0; kk < rmc_number; kk++)
                        {
                                //what number am I
                                rmc_index = rmc_number*rmc_number * ii + rmc_number * jj + kk;
                                
                                // arbitrary quaternions require 4 parameters; but our length is fixed. 
                                u1 = 1.0 / rmc_number * ii;
                                u2 = 1.0 / rmc_number * jj;
                                u3 = 1.0 / rmc_number * kk;

                                // http://planning.cs.uiuc.edu/node198.html 
                                w = sqrt(1-u1) * sin(2 * M_PI * u2);
                                x = sqrt(1-u1) * cos(2 * M_PI * u2);
                                y = sqrt(u1) * sin(2 * M_PI * u3);
                                z = sqrt(u1) * cos(2 * M_PI * u3);


                                //https://en.wikipedia.org/wiki/Rotation_group_SO%283%29#Quaternions_of_unit_norm

                                rmc_matrices[rmc_index][0] = 1 - 2 * (y*y+z*z);
                                rmc_matrices[rmc_index][1] = 2 * (x*y-z*w);
                                rmc_matrices[rmc_index][2] = 2 * (x*z+y*w);
                                rmc_matrices[rmc_index][3] = 2 * (x*y+z*w);
                                rmc_matrices[rmc_index][4] = 1 - 2 * (x*x+z*z);
                                rmc_matrices[rmc_index][5] = 2 * (y*z-x*w);
                                rmc_matrices[rmc_index][6] = 2 * (x*z-y*w);
                                rmc_matrices[rmc_index][7] = 2 * (y*z+x*w);
                                rmc_matrices[rmc_index][8] = 1 - 2 * (x*x+y*y); 
                                 
                        }
                }
        }   
}

/**** Change value of R[i] = SO(3) and sigma[i] = 1/-1 ****/
void build_rotation_matrix(int i) 
{         
        int build_random = (int) (rmc_number_total * dsfmt_genrand_close_open(&dsfmt));
        
        copy( begin(rmc_matrices[build_random]), end(rmc_matrices[build_random]), begin(R[i]));
        s[i] = dsfmt_genrand_close_open(&dsfmt) > 0.5 ? 1 : -1; 
}
        
/**** build gauge bath ****/
void build_gauge_bath()
{
/**** D2d checked ****/
	/** 0, Identity **/
	U_bath[0][0] = 1; 
	U_bath[0][4] = 1; 
	U_bath[0][8] = 1; 
	
	/** 1, C2,z **/
	U_bath[1][0] = -1; 
	U_bath[1][4] = -1; 
	U_bath[1][8] = 1;
	
	/** 2, -C4+,z **/
	U_bath[2][1] = 1; 
	U_bath[2][3] = -1; 
	U_bath[2][8] = -1;
	
	/** 3, -C4-,z **/
	U_bath[3][1] = -1; 
	U_bath[3][3] = 1; 
	U_bath[3][8] = -1;		
	
	/** 4, C2,y **/
	U_bath[4][0] = -1; 
	U_bath[4][4] = 1; 
	U_bath[4][8] = -1;			
	
	/** 5, C2,x **/
	U_bath[5][0] = 1; 
	U_bath[5][4] = -1; 
	U_bath[5][8] = -1;	
	
	/** 6, m x,-x. z **/
	U_bath[6][1] = -1; 
	U_bath[6][3] = -1; 
	U_bath[6][8] = 1;			

	/** 7, m x,x,z**/
	U_bath[7][1] = 1; 
	U_bath[7][3] = 1; 
	U_bath[7][8] = 1;	 
	
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

}

void random_initialization()
{
	int i, j;
	for (i = 0; i < L3; i++)
        {
			build_rotation_matrix(i);
		
			/** build Ux **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
			copy(begin(U_bath[j]),end(U_bath[j]),begin(Ux[i]));	
			
			/** build Uy **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt)); 
                        copy(begin(U_bath[j]),end(U_bath[j]),begin(Uy[i]));     
			
			/** build Uz **/
			j = int(U_order * dsfmt_genrand_close_open(&dsfmt)); 
                        copy(begin(U_bath[j]),end(U_bath[j]),begin(Uz[i]));     
			
			
        }
		
}

	
double site_energy(int i)
{
	/****** find neighbour, checked*****/
	
	int xp, xn, yp, yn, zp, zn; 
	
	xp = i % L == 0 ? i - 1 + L : i - 1;
	xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
	
	yp = i % L2 < L ? i - L + L2 : i - L;
	yn = (i + L) % L2 < L ? i + L - L2 : i + L;
	
	zp = i < L2 ? i - L2 + L3 : i - L2;
	zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
	
	/**** matrix fields of six O(3)/H bonds, s[i]s[j] R^T_i U_ij R_j ***/
	
	/** Rfoo[i], i from 0 to 5 save bonds value of
	 * -x,x, -y, y, -z, z
	 * foo[9] an intermieda variable
	 *  **/
	
	double Rfoo[6][9] = {{0}};
	double foo[9] = {0};
	
	/** Rfoo[0] =  s[xp] s[i] (Ux[i] R[i] R[xp])^T) 
	 * foo = s[i]*s[xp] R[i] R[xp]^T
	 * Rfoo[0] = Ux[i] foo
	 * **/
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[xp],
	R[i], 3, R[xp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[i], 3, foo,3,
	0.0, Rfoo[0],3);
	
	
	/** Rfoo[1] = s[i] s[xn] (Ux[xn] R[xn] R[i]^T) 
	 * **/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[xn]*s[i],
	R[xn], 3, R[i],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[xn], 3, foo,3,
	0.0, Rfoo[1],3);

	/** Rfoo[2] = s[yp] s[i] (Uy[i] R[i]R[yp]^T) 
	 * **/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[yp],
	R[i], 3, R[yp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[i], 3, foo,3,
	0.0, Rfoo[2],3);
	
	/** Rfoo[3] = s[i] s[yn] (Uy[yn] R[yn] R[i]^T)
	 * **/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[yn]*s[i],
	R[yn], 3, R[i],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[yn], 3, foo,3,
	0.0, Rfoo[3],3);		

	/** Rfoo[4] = s[zp] s[i] (Uz[i] R[i] R[zp]^T) 
	 * **/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[zp],
	R[i], 3, R[zp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[i], 3, foo,3,
	0.0, Rfoo[4],3);
		
	/** Rfoo[5] = s[i] s[zn] (Uz[zn] R[zn] R[i]^T)
	 * **/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[zn]*s[i],
	R[zn], 3, R[i],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[zn], 3, foo,3,
	0.0, Rfoo[5],3);		
	
	/******** total energy *****/ 
	double Sfoo = 0;
	
	for(int j = 0; j < 6; j++)
        {
                Sfoo += J1*Rfoo[j][0] + J2*Rfoo[j][4] + J3*Rfoo[j][8];
        }				
	return Sfoo;
}	
	
/**** i passed from main so not need to defined int again and again
 * s[i] also changed in build rotation
 * ****/ 

void flip_R(int i) 
{
	double E_old;
        double s_save;
        double R_save[9];
        double E_new;
        
	E_old = site_energy(i); 
        
	copy(begin(R[i]), end(R[i]),begin(R_save)); //save R[i] to R_save
	s_save = s[i];
	
	build_rotation_matrix(i); //generate new R[i] and s[i]
        
	E_new = site_energy(i);
	
	E_change = E_new - E_old;
	
	/*******decide flip and change E_total********/
	if (E_change < 0)
        {
                E_total += E_change;
                
        }
	else
        {
                double change_chance = exp(-beta * E_new) /  (exp(-beta * E_new) + exp(-beta * E_old) );
                if (change_chance > dsfmt_genrand_close_open(&dsfmt))
                {
                        E_total += E_change;
                        
                }
                else
                {
                        copy(begin(R_save), end(R_save), begin(R[i])); 
                        s[i] = s_save;
                } // change R[i] and s[i] back				
        }	
}

void flip_Ux(int i)
{
	int xp = 0; 

	xp = i % L == 0 ? i - 1 + L : i - 1; 
	
	double Bond[9] = {0};
	double foo[9] = {0};

	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[xp],
	R[i], 3, R[xp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[i], 3, foo,3,
	0.0, Bond,3);
	
	/**** bond enery ****/
	double E_old;
	E_old = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	/**** save Ux[i] to U_save ****/
	double U_save[9];
	copy(begin(Ux[i]),end(Ux[i]),begin(U_save));
	
	/**** generate new Ux by choosing from U_bath ****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(begin(U_bath[j]),end(U_bath[j]),begin(Ux[i]));
	
	/**** compute the bond enery s[xp] s[i] Ux[i] R[i] R^T[xp] again****/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[xp],
	R[i], 3, R[xp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Ux[i], 3, foo,3,
	0.0, Bond,3);
	
	double E_new;
	E_new = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
	if (E_change < 0)
	{
		E_total += E_change;
	}
	else
	{
                double change_chance = exp(-beta * E_new) /  (exp(-beta * E_new) + exp(-beta * E_old) );
                if (change_chance > dsfmt_genrand_close_open(&dsfmt))
                {
			E_total += E_change;
		}
		else
		{
			copy(begin(U_save),end(U_save),begin(Ux[i]));
		}	
	}		
	
}

void flip_Uy(int i)
{
	/**** find neighbour, only one****/
	int yp;
	yp = i % L2 < L ? i - L + L2 : i - L;
	
	/**** bond energy  s[yp] s[i] Uy[i] R[i] R^T[yp] ****/

	double Bond[9] = {0};
	double foo[9] = {0}; 
	
	/** Bond = s[yp] s[i] (Uy[i] R[i] R[yp]^T)
	 * foo = s[i] s[yp] R[i] R[yp]^T
	 * Bond = Uy[i] foo
	 * **/ 
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[yp],
	R[i], 3, R[yp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[i], 3, foo,3,
	0.0, Bond,3);
	
	/**** bond enery ****/
	double E_old;
	E_old = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	/**** save Uy[i] to U_save ****/
	double U_save[9];
	copy(begin(Uy[i]),end(Uy[i]),begin(U_save));
	
	/**** generate choosing new Uy from U_bath****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(begin(U_bath[j]),end(U_bath[j]),begin(Uy[i]));
	
	/**** compute the bond enery s[yp] s[i] Uy[i] R[i] R^T[yp] again****/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[yp],
	R[i], 3, R[yp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uy[i], 3, foo,3,
	0.0, Bond,3);
	
	double E_new;
	E_new = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
	if (E_change < 0)
	{
                E_total += E_change;
                
        }
        else {
                double change_chance = exp(-beta * E_new) /  (exp(-beta * E_new) + exp(-beta * E_old) );
                if (change_chance > dsfmt_genrand_close_open(&dsfmt))
                {
                        E_total += E_change;
                        
                }
                else {
                        copy(begin(U_save),end(U_save),begin(Uy[i]));
                        
                }	
        }		
	
}


void flip_Uz(int i)
{
	/**** find neighbour, only one****/
	int zp;
	zp = i < L2 ? i - L2 + L3 : i - L2;
	
	/**** bond energy s[zp] s[i] Uz[i] R[i] R^T[zp] ****/

	double Bond[9] = {0};
	double foo[9] = {0}; 
	
	/** Bond = s[zp] s[i] Uz[i] R[i] R^T[zp]
	 * foo = s[i] s[zp] R[i] R^T[zp]
	 * Bond = Uz[i] foo
	 * **/ 
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[zp],
	R[i], 3, R[zp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[i], 3, foo,3,
	0.0, Bond,3);
	
	/**** bond enery ****/
	double E_old;
	E_old = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	/**** save Uz[i] to U_save ****/
	double U_save[9];
	copy(begin(Uz[i]),end(Uz[i]),begin(U_save));
	
	/**** generate new Uz by choosing new U from U_bath****/
	int j;
	j = int(U_order * dsfmt_genrand_close_open(&dsfmt));
	
	copy(begin(U_bath[j]),end(U_bath[j]),begin(Uz[i]));
	
	/**** compute the bond enery s[zp] s[i] Uz[i] R[i] R^T[zp] again****/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,s[i]*s[zp],
	R[i], 3, R[zp],3,
	0.0, foo,3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	Uz[i], 3, foo,3,
	0.0, Bond,3);
	
	double E_new;
	E_new = J1*Bond[0] + J2*Bond[4] + J3*Bond[8];
	
	E_change = E_new - E_old;
	
	/**** decide flip and change E_total ****/
        if (E_change < 0)
        {
                E_total += E_change;
                
        }
        else {
                double change_chance = exp(-beta * E_new) /  (exp(-beta * E_new) + exp(-beta * E_old) );
                if (change_chance > dsfmt_genrand_close_open(&dsfmt))
                {
                        E_total += E_change;

                }
                else {
                        copy(begin(U_save),end(U_save),begin(Uz[i]));

                }	
        }		
}

double thermalization_inner ( int N)
{
        double s1 = 1.0 ;
        int i, j, site;  
        for (i = 0; i <  N; i++)
        {
                s1 += E_total;
                /**** one sweep ****/
                for (j = 0; j < L3*4; j++)
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
        return s1;
}
void thermalization()
{
	double afoo = 0;

        double s1 = 0.0;
        double s2 = thermalization_inner ( 1000);
	
	while (afoo < 1-accurate || afoo > 1+accurate)
        {
                s2 = s1; // copy old value of s1 to s2.

                /**** updata configurations and calculate new s1 ****/ 
                s1 = thermalization_inner( 10 );
                
                afoo = s1/s2; 
        } 	
}
	
void estimate_beta_c()
{
	double S1, S2, Cv; 
	double foo2_n, Q1_n, Q2_n, chi_n;
	double foo_s, s1, s2, chi_s;
	int i, j, site;
  
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
	{ 
		//printf("#//Calculating for beta=%2.3f < %2.3f \n", beta, beta_upper);

		/** re-initialize quantites for the acception ratio **/
		//Racc = 0; Rrej = 0; xacc = 0; xrej = 0;
		//yacc = 0; yrej = 0; zacc = 0; zrej = 0;
		/**** re-thermalization and reset all quantities need for one temperature****/
		thermalization();
		S1 = 0; S2 = 0;
		s1 = 0; s2 = 0;
		Q1_n = 0; Q2_n = 0;
		/**** measure ****/ 	
                 
		for (j = 0; j < sample_amount; j++)
		{
			//printf("Num threads %d \n", omp_get_num_threads());
			//line used to check core number. 
                         
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
			S1 += E_total;	 
			S2 += E_total * E_total;	

			foo2_n = orderparameter_n(); //intensive
			Q2_n += foo2_n;
			Q1_n += sqrt(foo2_n);					

			foo_s = 0;
			for(int k = 0; k < L3; k++) 
			{
				foo_s += s[k];
			} //extensive
			foo_s /= L3; // intensive
			s1 += foo_s;
			s2 += foo_s*foo_s;
		}	 

		S1 /= sample_amount;
		S2 /= sample_amount;

		Cv = (S2 - S1 * S1) * beta * beta / L3;	 
		S1 /= E_g;	

		Q1_n/=sample_amount;
		Q2_n /= sample_amount;
		chi_n = (Q2_n - Q1_n*Q1_n)*beta*L3;									

		s1 /= sample_amount;
		s2 /= sample_amount;
		chi_s = (s2 -s1*s1)*beta*L3;

		printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n", beta, S1, Cv, s1, chi_s, Q1_n, chi_n); 
		if ( (beta >= beta_1) && (beta <= beta_2) )
		{
			beta += beta_step_small;
		}
		else
		{
			beta += beta_step_big;
		} 		 		
	} 
}


void measure_energy(char *output)
{
	int i, j, site;
	
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
		  for (j = 0; j < sample_amount; j++)
		  { for (i = 0; i < L3*4*tau ; i++)
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
			output_file << E_total << '\t' << flush; 
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
	

double orderparameter_n()
{
	double Q11_n = 0, Q22_n = 0, Q33_n = 0, Q12_n = 0, Q23_n = 0, Q13_n = 0, Q2_n = 0;
	
			for(int i = 0; i < L3; i++) 
				{Q11_n += 1.5*R[i][6] * R[i][6] - 0.5;}						
			for(int i = 0; i < L3; i++) 
				{Q22_n += 1.5*R[i][7] * R[i][7] - 0.5;}							
			for(int i = 0; i < L3; i++) 
				{Q33_n += 1.5*R[i][8] * R[i][8] - 0.5;}	
			for(int i = 0; i < L3; i++) 
				{Q12_n += 1.5*R[i][6] * R[i][7];}
			for(int i = 0; i < L3; i++) 
				{Q13_n += 1.5*R[i][6] * R[i][8];}
			for(int i = 0; i < L3; i++) 
				{Q23_n += 1.5*R[i][7] * R[i][8];}				
						
	Q2_n = (Q11_n*Q11_n + Q22_n*Q22_n + Q33_n *Q33_n
				 + 2*Q12_n*Q12_n + 2*Q13_n*Q13_n + 2*Q23_n*Q23_n)/L3/L3;	
	
	return Q2_n;
	}			
	
	
void measure_J_distribution(char *output)
{	
	/**
	 * Computing J_eff = Tr(R^T_i J U_{ij} R_j), R is defened no s.
	 * J_eff[10][3*L3] are defined to save the value of J_eff of the 3N bonds,
	 * beta will be couted by itself in the save outout file before J[i][0]
	 * 10 is an abitrary choose for the number of silce
	 * numbers of samples use the sample_amount in the header
	 * cout -J_eff in output, since J is defined with a minus sign
	 * **/
	 
	double foo[9] = {0}, xfoo[9]={0}, yfoo[9]={0}, zfoo[9]={0}; //xfoo = R^T[i] Ux[xn] R[xn]
	int site;
	int xn, yn, zn; 
	
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
			
			 for ( int j = 0; j < L3*4*tau ; j++) // to generate a new sample without autocorrealtion 
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
					
			/** sweep the lattice **/
			for( int i = 0; i < L3; i++)
				{
					/** J_eff at x direction 
					 * foo = R_j R^T_i = R[xn]R^T[i]
					 * xfoo = U[xn] foo
					 * as in Rfoo[1] but no s[i] fields
					 *  **/

					xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
					
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[xn], 3, R[i],3,
					0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Ux[xn], 3, foo,3,
					0.0, xfoo,3);
					
					J_eff[3*i] += J1 * xfoo[0] + J2 * xfoo[4] + J3 * xfoo[8];
					
					/** J_eff at y direction
					 * foo = R_j R^T_i = R[yn]R^T[i]
					 * yfoo = Uy[yn] foo
					 *  **/
					
					yn = (i + L) % L2 < L ? i + L - L2 : i + L;

					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[yn], 3, R[i],3,
					0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Uy[yn], 3, foo,3,
					0.0, yfoo,3);
					
					J_eff[3*i + 1] += J1 * yfoo[0] + J2 * yfoo[4] + J3 * yfoo[8];
					
					
					/** J_eff at z direction 
					 * foo = R_j R^T_i = R[zn]R^T[i]
					 * zfoo = Uz[zn] foo
					 * **/
					
					zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
					
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[zn], 3, R[i],3,
					0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Uz[zn], 3, foo,3,
					0.0, zfoo,3);
					
					J_eff[3*i + 2] += J1 * zfoo[0] + J2 * zfoo[4] + J3 * zfoo[8];					
					
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
	 * Computing J_eff = <Tr(R^T_i J U_{ij} R_j)>, R is defened no s.
	 * cout -J_eff in output, since J is defined with a minus sign
	 * **/
	 
	double foo[9] = {0}, xfoo[9]={0}, yfoo[9]={0}, zfoo[9]={0}; //xfoo = R^T[i] Ux[xn] R[xn]
	int site;
	int xn, yn, zn; 
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
			
			 for ( int j = 0; j < L3*4*tau ; j++) // to generate a new sample without autocorrealtion 
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
					
			/** sweep the lattice **/
			for( int i = 0; i < L3; i++)
				{
					/** J_eff at x direction 
					 * foo = R_j R^T_i = R[xn]R^T[i]
					 * xfoo = U[xn] foo
					 * as in Rfoo[1] but no s[i] fields
					 *  **/

					xn = (i + 1) % L == 0 ? i + 1 - L : i + 1;
					
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[xn], 3, R[i],3,
					0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Ux[xn], 3, foo,3,
					0.0, xfoo,3);
					
					J_eff += J1 * xfoo[0] + J2 * xfoo[4] + J3 * xfoo[8];
					
					/** J_eff at y direction
					 * foo = R_j R^T_i = R[yn]R^T[i]
					 * yfoo = Uy[yn] foo
					 *  **/
					
					yn = (i + L) % L2 < L ? i + L - L2 : i + L;

					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[yn], 3, R[i],3, 0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Uy[yn], 3, foo,3,
					0.0, yfoo,3);
					
					J_eff += J1 * yfoo[0] + J2 * yfoo[4] + J3 * yfoo[8];
					
					
					/** J_eff at z direction 
					 * foo = R_j R^T_i = R[zn]R^T[i]
					 * zfoo = Uz[zn] foo
					 * **/
					
					zn = i + L2 >= L3 ? i + L2 - L3 : i + L2;
					
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
					3,3,3,1,
					R[zn], 3, R[i],3,
					0.0, foo,3);
	
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					3,3,3,1,
					Uz[zn], 3, foo,3,
					0.0, zfoo,3);
					
					J_eff += J1 * zfoo[0] + J2 * zfoo[4] + J3 * zfoo[8];					
					
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
