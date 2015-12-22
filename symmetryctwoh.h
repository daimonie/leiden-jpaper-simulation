#ifndef SYMMETRY_CTWOH_H
#define SYMMETRY_CTWOH_H

#include "symmetry.h"
class symmetry_c2h : public symmetry
{
	public:
		string label 	= "C2h";
		int bath_size	= 4;
		void bath(); 
};


#endif