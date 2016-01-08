#ifndef ORDER_DUMMY_H
#define ORDER_DUMMY_H

#include "order.h"   
class order_d4h_2 : public order
{
	public:   
		double calculate(simulation *);
		double dfc(int,int);
};
#endif