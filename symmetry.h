#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <string>
using namespace std;
class symmetry
{
	public: 
		virtual ~symmetry(){}
		virtual void bath()		= 0;
		virtual string label()		= 0;
		virtual int bath_size ()	= 0;
		
                double bath_field[100][9]	= {{0}}; 
};
#endif