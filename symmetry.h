#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <string>
class symmetry
{
	public:
		string label 			= "default";
		int bath_size			= 8;
		virtual void bath() 		= 0; 
                double bath_field[100][9]	= {{0}};
};
#endif