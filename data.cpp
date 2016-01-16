#include "data.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <boost/format.hpp>
using namespace std;
void data::report ()
{ 
    stringstream format_order;
    for(int i = 0; i < (int) order.size(); i++)
    {
        format_order <<  boost::format("\t%.8f") % order[i];
        format_order << boost::format("\t%.8f") % chi_order[i];
    }
    printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f%s\t%s\n",
        beta, total_energy, heat_capacity,finite_k, j_one, j_two, j_three, format_order.str().c_str(),
        point_group.c_str());
}
void data::shout (FILE * file_handler)
{ 
    stringstream format_order;
    for(int i = 0; i < (int) order.size(); i++)
    {
        format_order <<  boost::format("\t%.8f") % order[i];
        format_order << boost::format("\t%.8f") % chi_order[i];
    }
    fprintf(file_handler, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f%s\t%s\n",
        beta, total_energy, heat_capacity,finite_k, j_one, j_two, j_three, format_order.str().c_str(),
        point_group.c_str());
}