#ifndef SYMMETRY_DTWOH_H
#define SYMMETRY_DTWOH_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_d2h : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif