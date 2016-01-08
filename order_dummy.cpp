#include "order.h" 
#include "order_dummy.h"
#include "simulation.h" 
double order_dummy::calculate(simulation * sweep)
{ 
	(void)sweep;
	
	return 0.0;
}