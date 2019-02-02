#include<iostream>
#include"MyMatrix.h"
#include<time.h>
using namespace std;

int main()
{
	//double array1[] = {8, -3, 2, 4, 11, -1, 6, 3, 12};
	double array1[] = {0.269733, -0, -24.2 ,0.199247, -70.4, 0, -0.942096, -82, 0};
	double array2[] = {20, 33, 36};
	MySquareMatrix<double> A(3, array1);
	MyMatrix<double> B(3, 1, array2);
	MyMatrix<double> x(3, 1);
	x = A.SolveEqu(B);
	cout<<x;
	x = A.Jacobi(B, 20);
	cout<<x;
	x = A.Gauss_Seidel(B, 10);
	cout<<x;
	x = A.SOR(B, 20, 0.5);
	cout<<x;
	return 0;
}
