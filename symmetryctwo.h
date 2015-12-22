#ifndef SYMMETRY_CTWO_H
#define SYMMETRY_CTWO_H

#include "symmetry.h"
class symmetry_c2 : public symmetry
{
	public:
		int bath_size	= 2;
		void bath(); 
		string label 	= "C2";
};


#endif