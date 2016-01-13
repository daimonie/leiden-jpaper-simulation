#include "order.h" 
#include "order_d4h_2.h"
#include "simulation.h" 
double order_d4h_2::calculate(simulation * sweep)
{  
	double Q[3][3][3][3][3][3];
	double foo_lmn, foo_mln; // for storing llmmnn, mmllnn
	int a, b, c, d, e, f, i;
	double Q2 = 0;	
	
        for(a = 0; a < 3; a++)
        {	
            for(b = 0; b < 3; b++)
            {	
                for(c = 0; c < 3; c++)
                {	
                    for(d = 0; d < 3; d++)
                    {	
                        for(e = 0; e < 3; e++)
                        {	
                            for(f = 0; f < 3; f++)
                            {
                                /** reset foo_l.m.n **/
                                foo_lmn = 0; foo_mln = 0;

                                /** average foo_l, foo_m, foo_n first **/
                                for(i = 0; i <  sweep->length_three; i++)
                                {	
                                    foo_lmn += sweep->field_r[i][a] * sweep->field_r[i][b] * sweep->field_r[i][3+c] * sweep->field_r[i][3+d] * sweep->field_r[i][6+e] * sweep->field_r[i][6+f];

                                    foo_mln += sweep->field_r[i][3+a] * sweep->field_r[i][3+b] * sweep->field_r[i][c] * sweep->field_r[i][d] * sweep->field_r[i][6+e] * sweep->field_r[i][6+f];	
                                }	

                                foo_lmn /=  sweep->length_three;
                                foo_mln /=  sweep->length_three;

                                /** Q =  (foo_lmn+ foo_mln)  **/		
                                Q[a][b][c][d][e][f] =  foo_lmn + foo_mln; 	 
                                Q2 += Q[a][b][c][d][e][f] * Q[a][b][c][d][e][f];	
                            } 
                        }
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