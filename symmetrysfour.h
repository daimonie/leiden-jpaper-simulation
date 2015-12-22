#ifndef SYMMETRY_SFOUR_H
#define SYMMETRY_SFOUR_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_s4 : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif