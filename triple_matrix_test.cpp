 
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h> 
#include <string>

using namespace std; 
void print_matrix (string , double[]);

int main()
{
	printf("Hello. \n");
	/***
uxi = np.matrix([[-1.000,0.000,0.000],[0.000,1.000,0.000],[0.000,0.000,-1.000])
ri = np.matrix([[-0.000,-0.951,0.309],[-0.588,0.250,0.769],[-0.809,-0.182,-0.559])
rj = np.matrix([[0.400,0.872,-0.283],[-0.000,-0.309,-0.951],[-0.917,0.380,-0.124])
***/
	double A[9] = {0.400,0.872,-0.283,-0.000,-0.309,-0.951,-0.917,0.380,-0.124};	
	double C[9] = {-0.000,-0.951,0.309,-0.588,0.250,0.769,-0.809,-0.182,-0.559};	
	double B[9] = {0.400,0.872,-0.283,-0.000,-0.309,-0.951,-0.917,0.380,-0.124};	
	
	print_matrix("A", A);
	print_matrix("B", B);
	print_matrix("C", C);
	
	
	double foo[9] = {0};
	double foo2[9] = {0};
	double foo3[9] = {0};
	 
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	3,3,3,1,
	B, 3, A,3,
	0.0, foo,3);
	
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	3,3,3,1,
	C, 3, foo,3,
	0.0, foo2,3);
	print_matrix("foo2", foo2);
	
	int index, i, j, k, l;
	
	double sum = 0.0;
	
	#pragma omp parallel for reduction(+ : sum, foo)
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			index = i*3 + j;
			for(l = 0; l < 3; l++)
			{		
				for(k = 0; k < 3; k++)
				{
 					foo3[index] += C[i*3 + l] * B[l*3 + k] * A[k+j*3];
					sum += C[i*3 + l] * B[l*3 + k] * A[k+j*3];
				}
			}
		}
	}
	print_matrix("foo3", foo3);  
	printf("Found sum %.3f .\n", sum);
}
void print_matrix(string name, double matrix[])
{
	printf("%s = np.matrix([[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f],[%.3f,%.3f,%.3f])\n",
		name.c_str(), matrix[0], matrix[1], matrix[2],
		matrix[3], matrix[4], matrix[5],
		matrix[6], matrix[7], matrix[8]);
}