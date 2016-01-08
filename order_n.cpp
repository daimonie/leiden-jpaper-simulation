#include "order.h" 
#include "order_n.h"
#include "simulation.h" 
double order_n::calculate(simulation * sweep)
{ 
        double one_one = 0, two_two = 0, three_three = 0, one_two = 0, two_three = 0, one_three = 0; 
        for(int i = 0; i < sweep->length_three; i++) 
        {
                one_one         += 1.5*sweep->field_r[i][6] * sweep->field_r[i][6] - 0.5; 
                two_two         += 1.5*sweep->field_r[i][7] * sweep->field_r[i][7] - 0.5; 
                three_three     += 1.5*sweep->field_r[i][8] * sweep->field_r[i][8] - 0.5; 
                one_two         += 1.5*sweep->field_r[i][6] * sweep->field_r[i][7]; 
                one_three       += 1.5*sweep->field_r[i][6] * sweep->field_r[i][8]; 
                two_three       += 1.5*sweep->field_r[i][7] * sweep->field_r[i][8];
                
        }                               
                                                
        double q = (one_one*one_one + two_two*two_two + three_three *three_three + 2*one_two*one_two + 2*one_three*one_three + 2*two_three*two_three) / sweep->length_three / sweep->length_three;        
        
        return q;
}  