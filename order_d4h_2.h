#ifndef ORDER_D4H_2_H
#define ORDER_D4H_2_H

#include "order.h"   
class order_d4h_2 : public order
{
	public:   
		double calculate(simulation *);
		double dfc(int,int);
};
#endif