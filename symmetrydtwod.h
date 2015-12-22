#ifndef SYMMETRY_DTWOD_H
#define SYMMETRY_DTWOD_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_d2d : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif