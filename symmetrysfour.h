#ifndef SYMMETRY_SFOUR_H
#define SYMMETRY_SFOUR_H

#include "symmetry.h"
class symmetry_s4 : public symmetry
{
	public:
		int bath_size	= 4;
		void bath(); 
};


#endif