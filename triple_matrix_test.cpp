 
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h> 
#include <string>

using namespace std; 
void print_matrix (string , double[]);

int main()
{
	printf("Hello. \n");
	
	double A[9] = {0.1, 0, 0, 0, 0.1, 0, 0, 0, -0.1};	
	print_matrix("A", A);
	
	
	double B[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0, 0, 0, -0.1};	
	print_matrix("A", B);
	
	
	double C[9] = {0, 0, 0.1, 0, 0.1, 0, -0.1, 0, 0};	
	print_matrix("A", C);
	
	
	double foo[9] = {0};
	double foo2[9] = {0};
	double foo3[9] = {0};
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,1,
	B, 3, A,3,
	0.0, foo,3);
	
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,1,
	C, 3, foo,3,
	0.0, foo2,3);
	print_matrix("foo2", foo2);
	
	int index, i, j, k, l;
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			index = i*3 + j;
			for(l = 0; l < 3; l++)
			{		
				for(k = 0; k < 3; k++)
				{
					foo3[index] += C[i*3 + l] * B[l + k*3] * A[k*3+j];
				}
			}
		}
	}
	
	print_matrix("foo3", foo3);  
}
void print_matrix(string name, double matrix[])
{
	printf("%s = np.matrix([[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f])\n",
		name.c_str(), matrix[0], matrix[1], matrix[2],
		matrix[3], matrix[4], matrix[5],
		matrix[6], matrix[7], matrix[8]);
}