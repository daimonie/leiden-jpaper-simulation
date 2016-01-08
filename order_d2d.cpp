#include "order.h" 
#include "order_d2d.h"
#include "simulation.h" 
double order_d2d::calculate(simulation * sweep)
{ 
        int a, b, c, i; 
        
        double q_sum = 0.0;
        double order = 0;
        for(a = 0; a < 3; a++)
        {       for(b = 0; b < 3; b++)
                {       for(c = 0; c < 3; c++)
                        {       
                                q_sum = 0.0;
                                for(i = 0; i < sweep->length_three; i++)
                                {
                                        q_sum += sweep->field_s[i] * (sweep->field_r[i][a] * sweep->field_r[i][3+b] + sweep->field_r[i][b] * sweep->field_r[i][3+a]) * sweep->field_r[i][6+c];
                                }  
                                order += q_sum * q_sum;
                        }
                }
        } 
        return order / sweep->length_three / sweep->length_three;             
} 