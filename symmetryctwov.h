#ifndef SYMMETRY_CTWOV_H
#define SYMMETRY_CTWOV_H

#include "symmetry.h"
class symmetry_c2v : public symmetry
{
	public:
		int bath_size	= 4;
		void bath(); 
};


#endif