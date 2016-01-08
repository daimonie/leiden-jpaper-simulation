#ifndef ORDER_D4H_H
#define ORDER_D4H_H

#include "order.h"   
class order_d4h : public order
{
	public:   
		double calculate(simulation *);
		double dfc(int,int);
};
#endif