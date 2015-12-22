#ifndef SYMMETRY_H
#define SYMMETRY_H
class symmetry
{
	public:
		int bath_size			= 8;
		virtual void bath() 		= 0; 
                double bath_field[100][9]	= {{0}};
};
#endif