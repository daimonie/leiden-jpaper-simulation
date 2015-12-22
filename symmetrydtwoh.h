#ifndef SYMMETRY_DTWOH_H
#define SYMMETRY_DTWOH_H

#include "symmetry.h"
class symmetry_d2h : public symmetry
{
	public:
		int bath_size	= 8;
		void bath(); 
		string label 	= "D2h";
};


#endif