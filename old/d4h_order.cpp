double dfc(int a, int b)
{
	if (a == b) return 1;
		else return 0;
	}

double orderparameter_D4h_1()
{
	/** define the order parameter
	 * checked at random limit (Q2 = 0) and uniform limit (Q2 = 30/11),
	 *  **/
	
	double Q[3][3][3][3];
	double foo_lm, foo_ml; // for storing llll, mmmm and nnnn
	int a, b, c, d, i;
	
	for(a = 0; a < 3; a++)
		{	for(b = 0; b < 3; b++)
				{	for(c = 0; c < 3; c++)
						{
							for(d = 0; d < 3; d++)
								{	
									/** reset foo_l.m.n **/
									foo_lm = 0; foo_ml = 0;
									
									/** average llmm, mmll first **/
									for(i = 0; i < L3; i++)
										{	foo_lm += R[i][a] * R[i][b] * R[i][3+c] * R[i][3+d];
											foo_ml += R[i][3+a] * R[i][3+b] * R[i][c] * R[i][d]; 
											}	
									foo_lm /= L3;
									foo_ml /= L3;
									
									/** Q = 2.5 (llmm+mmll) - delta_functions **/		
									Q[a][b][c][d] = 15.0/11.0 * (foo_lm + foo_ml) - 4.0/11.0 * dfc(a,b) * dfc(c,d) 
													+ 1.0/11.0 * ( dfc(a,c) * dfc(b,d) + dfc(a,d) * dfc(b,c) );
									}
							
							}
					}
			}
	
	/**** trace of Q ****/
	double Q2 = 0;
	
	for(a = 0; a < 3; a++)
		{	for(b = 0; b < 3; b++)
				{	for(c = 0; c < 3; c++)
						{	for(d = 0; d < 3; d++)
								{	Q2 += Q[a][b][c][d] * Q[a][b][c][d];
									}
							}
					}
			}
	
	return Q2;		
	
	}
	
double orderparameter_D4h_2()
{
	/** define the order parameter
	 * checked at random limit (Q2 = 0) and uniform limit (Q2 = 2),
	 *  **/
	
	double Q[3][3][3][3][3][3];
	double foo_lmn, foo_mln; // for storing llmmnn, mmllnn
	int a, b, c, d, e, f, i;
	
	for(a = 0; a < 3; a++)
	{	for(b = 0; b < 3; b++)
		{	for(c = 0; c < 3; c++)
			{	for(d = 0; d < 3; d++)
				{	for(e = 0; e < 3; e++)
					{	for(f = 0; f < 3; f++)
						{
							/** reset foo_l.m.n **/
							foo_lmn = 0; foo_mln = 0;
									
							/** average foo_l, foo_m, foo_n first **/
							for(i = 0; i < L3; i++)
							{	
								foo_lmn += R[i][a] * R[i][b] * R[i][3+c] * R[i][3+d] * R[i][6+e] * R[i][6+f];
								
								foo_mln += R[i][3+a] * R[i][3+b] * R[i][c] * R[i][d] * R[i][6+e] * R[i][6+f];	
							}	
							
							foo_lmn /= L3;
							foo_mln /= L3;
									
							/** Q =  (foo_lmn+ foo_mln)  **/		
							Q[a][b][c][d][e][f] =  foo_lmn + foo_mln; 	
						}
					}				
									
									
									
				}
							
			}
		}
	}
	
	/**** trace of Q ****/
	double Q2 = 0;
	
	for(a = 0; a < 3; a++)
	{	for(b = 0; b < 3; b++)
		{	for(c = 0; c < 3; c++)
			{	for(d = 0; d < 3; d++)
				{	for(e = 0; e < 3; e++)
					{	for(f = 0; f < 3; f++)
						{
							
							//cout << a << b << c << d << e << f << '\t' << Q[a][b][c][d][e][f] << '\t' << delta_function_Ih(a,b,c,d,e,f)<< endl;
						Q2 += Q[a][b][c][d][e][f] * Q[a][b][c][d][e][f];	
						}
					
					}
				}
			}
		}
	}
	
	return Q2;		
	
	}	

void measure_order_parameter_D4h(char *output)
{
	double foo2_1, Q1_1, Q2_1, chi_Q1;
	double foo2_2, Q1_2, Q2_2, chi_Q2;
	int i, j, site; 
	
	ofstream output_file;
	output_file.open(output);
	
	while ( ( (beta >= beta_lower) && (beta <= beta_upper) ))
		{ 
/** re-initialize quantites for the acception ratio **/
//Racc = 0; Rrej = 0; xacc = 0; xrej = 0;
//yacc = 0; yrej = 0; zacc = 0; zrej = 0;
			
		/**** re-thermalization and reset Q1, Q2****/
		  thermalization();
		  Q1_1 = 0; Q2_1 = 0;
		  Q1_2 = 0; Q2_2 = 0;
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
			
			foo2_1 = orderparameter_D4h_1();
			Q2_1 += foo2_1;
			Q1_1 += sqrt(foo2_1);
			
			foo2_2 = orderparameter_D4h_2();
			Q2_2 += foo2_2;
			Q1_2 += sqrt(foo2_2);
			 }	 
		
		Q1_1 /=sample_amount;
		Q2_1 /= sample_amount;
		chi_Q1 = (Q2_1 - Q1_1*Q1_1)*beta*L3;
		
		Q1_1 /=sqrt(30.0/11); // nomalize Q1, but be careful. This should be after chi_Q
		
		
		Q1_2 /=sample_amount;
		Q2_2 /= sample_amount;
		chi_Q2 = (Q2_2 - Q1_2*Q1_2)*beta*L3;
		
		Q1_2 /=sqrt(2); // nomalize Q1, but be careful. This should be after chi_Q

		output_file << 1/beta << '\t' << Q1_1 << '\t' << chi_Q1 << '\t' << Q1_2 << '\t' << chi_Q2<< endl;
		
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