#ifndef ORDER_H
#define ORDER_H

// #include "simulation.h"
class simulation; // forward declaration instead of circular dependency 
class order
{
	public: 
		virtual ~order(){}
		virtual double calculate(simulation *) = 0;
};
#endif