#ifndef ORDER_DUMMY_H
#define ORDER_DUMMY_H

#include "order.h"   
class order_dummy : public order
{
	public:   
		double calculate(simulation *);
};
#endif