/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	6// enter your assignment number
#define eval		0		// set 0

//#include "myNM.h"

#include "C:\Users\82104\Documents\GitHub\tutorialNM\include/myNM.h"

int main(int argc, char* argv[])
{
	/*	 [※ DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	double x[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
	double t[] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0 };
	double dydx[sizeof(x)/sizeof(double)];
	Matrix mx= arr2Mat(x,21,1);
	Matrix mt= arr2Mat(t, 21, 1);
	Matrix gd = createMat(21, 1);
	Matrix gdf = createMat(21, 1);
	int m=sizeof(x)/sizeof(double);
	double x0 = 0.2;
	double tol = 0.00001;
	double NR_result;
	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	gd = gradient(mt, mx);
	gradient(t, x, dydx, m);
	gdf = gradientFunc(myFunc, mt);
	// 뉴턴 랩손은 출력란에 실행 함수를 넣음. 
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("\n*************************************************");
	printf("\n                  Input Matrix ");
	printf("\n*************************************************\n");
	printMat(mt, "t");
	printMat(mx, "x");
	printf("\n*************************************************");
	printf("\n           differentiation of Matrix ");
	printf("\n*************************************************\n");
	printMat(gd, "gradient");

	printf("\n*************************************************");
	printf("\n           differentiation of 1D array ");
	printf("\n*************************************************\n");
	printArr(dydx, "1D gradient",m);


	printf("\n*************************************************");
	printf("\n           differentiation of My Function ");
	printf("\n*************************************************\n");
	printMat(gdf, "Func gradient");


	printf("\n*************************************************");
	printf("\n           Newton-Rhapson method  Results     ");
	printf("\n*************************************************\n");


	printf("Newton-Rhapson method:\n");
	NR_result = newtonRaphsonfunc(func, dfunc, x0, tol);
	printf("Final NR Solution : %f \t", NR_result);
	printf("\n");

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(mt);		freeMat(mx);			freeMat(gd); freeMat(gdf);

	system("pause");
	return 0;
}