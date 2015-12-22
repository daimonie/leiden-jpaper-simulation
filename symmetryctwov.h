#ifndef SYMMETRY_CTWOV_H
#define SYMMETRY_CTWOV_H

#include <string>
#include "symmetry.h"
using namespace std;
class symmetry_c2v : public symmetry
{
	public:
		void bath();
		string label();
		int bath_size();
};


#endif