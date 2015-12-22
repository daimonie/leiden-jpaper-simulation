#ifndef SYMMETRY_DTWOD_H
#define SYMMETRY_DTWOD_H

#include "symmetry.h"
class symmetry_d2d : public symmetry
{
	public:
		int bath_size	= 8;
		void bath(); 
};


#endif