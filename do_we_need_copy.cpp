#include <iostream>
#include <string>
using namespace std;
void print_matrix(string name, double matrix[])
{
        printf("%s = np.matrix([[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f])\n",
                name.c_str(), matrix[0], matrix[1], matrix[2],
                matrix[3], matrix[4], matrix[5],
                matrix[6], matrix[7], matrix[8]);
}  
int main(int argc, char **argv)
{
        double first[9] = {0};
        
        for(int j = 0; j < 10; j++)
        {
                first[j] = j+0.1;
        }
        
        auto second = first;
        
        double third[9] = {0};
        
        copy(begin(first), end(first), begin(third));
        
        for(int j = 0; j < 10; j++)
        {
                first[j] = 1/first[j];
        }
        
        printf("Do we need copy? Assign first, set second = first, change first. \n");
        
        print_matrix("first:", first);
        print_matrix("second:", second);
        print_matrix("third:", third);
        
       
        return 0;
}




