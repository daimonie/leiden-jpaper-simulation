#ifndef SYMMETRY_STWO_H
#define SYMMETRY_STWO_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_s2 : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif