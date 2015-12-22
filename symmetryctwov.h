#ifndef SYMMETRY_CTWOV_H
#define SYMMETRY_CTWOV_H

#include "symmetry.h"
class symmetry_c2v : public symmetry
{
	public:
		string label 	= "C2v";
		int bath_size	= 4;
		void bath(); 
};


#endif