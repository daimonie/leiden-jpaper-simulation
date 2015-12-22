#ifndef SYMMETRY_CTWO_H
#define SYMMETRY_CTWO_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_c2 : public symmetry
{
	public: 
		void bath();
		string label();
		int bath_size();
};


#endif