/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	5// enter your assignment number
#define eval		0		// set 0

#include "../../include/myNM.h"

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
	Matrix vecT = txt2Mat(path, "vecT");
	Matrix vecP = txt2Mat(path, "vecP");
	double fT = 0;
	double T = 150;
	double fTi = 0;
	double Ti = 75;
	double Tquery[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
	Matrix xq= arr2Mat(Tquery, 21, 1);
	Matrix yq = createMat(xq.rows, 1);
	for (int i = 0; i < xq.rows; i++) {

	}
	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	Matrix z = linearFit(vecT, vecP);
	yq = linearInterp(vecT, vecP, xq);
	fT = FindValueLCF(z, T);
	
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("\n*************************************************");
	printf("\n                  Input Matrix ");
	printf("\n*************************************************\n");
	printMat(vecT, "T");
	printMat(vecP, "P");

	printf("\n*************************************************");
	printf("\n Curvefitting  Coefficient z and Predict value for T ");
	printf("\n*************************************************\n");

	printMat(z, "z");
	printf("\n      Predict the pressure at T = %0.f",T);
	printf("\n               f(%f) = %f \n",T ,fT);


	printf("\n*************************************************");
	printf("\n		 Linear Spline interpolation");
	printf("\n*************************************************\n");

	printMat(yq, "yq");

	fTi = FindValueLIP(yq, xq, Ti);  // 오류 출력을 명확하게 보이기 위해 일부러 이곳에서 T에 따른 추정값을 찾는 함수를 놓음
	printf("\n      Predict the pressure at T = %0.f", Ti);
	printf("\n               f(%f) = %f \n\n", Ti, fTi);
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(vecT);		freeMat(vecP);		








	system("pause");
	return 0;
}