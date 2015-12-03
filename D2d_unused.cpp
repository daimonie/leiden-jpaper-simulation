		
	


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
