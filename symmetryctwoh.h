#ifndef SYMMETRY_CTWOH_H
#define SYMMETRY_CTWOH_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_c2h : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif