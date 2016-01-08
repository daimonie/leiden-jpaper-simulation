#include "order.h" 
#include "order_d4h_2.h"
#include "simulation.h" 
double order_d4h_2::calculate(simulation * sweep)
{  
	double Q[3][3][3][3];
	double foo_lm, foo_ml; // for storing llll, mmmm and nnnn
	int a, b, c, d, i;
	double Q2 = 0;	
	for(a = 0; a < 3; a++)
	{	
		for(b = 0; b < 3; b++)
		{	
			for(c = 0; c < 3; c++)
			{
				for(d = 0; d < 3; d++)
				{	 
 
					foo_lm = 0; foo_ml = 0;
 
					for(i = 0; i < sweep->length_three; i++)
					{
						foo_lm += sweep->field_r[i][a] * sweep->field_r[i][b] * sweep->field_r[i][3+c] * sweep->field_r[i][3+d];
						foo_ml += sweep->field_r[i][3+a] * sweep->field_r[i][3+b] * sweep->field_r[i][c] * sweep->field_r[i][d]; 
					}	
					foo_lm /= sweep->length_three;
					foo_ml /= sweep->length_three;
 
					Q[a][b][c][d] = 15.0/11.0 * (foo_lm + foo_ml) - 4.0/11.0 * dfc(a,b) * dfc(c,d) 
					+ 1.0/11.0 * ( dfc(a,c) * dfc(b,d) + dfc(a,d) * dfc(b,c) );
					
					Q2 += Q[a][b][c][d] * Q[a][b][c][d];
				}
			}
		}
	}

	return Q2;	 
}
double order_d4h_2::dfc(int a, int b)
{
	if (a == b)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}