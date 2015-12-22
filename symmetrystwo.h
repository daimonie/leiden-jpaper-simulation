#ifndef SYMMETRY_STWO_H
#define SYMMETRY_STWO_H

#include "symmetry.h"
class symmetry_s2 : public symmetry
{
	public:
		int bath_size	= 2;
		void bath(); 
};


#endif