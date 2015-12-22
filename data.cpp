#include "data.h"
#include <iostream>

void data::report ()
{
        printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%s\n",
	       beta, total_energy, heat_capacity, j_one, j_two, j_three, order_one, order_two, chi_order_one, chi_order_two,
	       point_group.c_str());
}