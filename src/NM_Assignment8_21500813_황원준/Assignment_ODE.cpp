/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	8// enter your assignment number
#define eval		0		// set 0

//#include "myNM.h"

#include "../../include/myNM.h"


int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
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
	double t0 = 0;
	double tf = 0.1;
	double h = 0.001;
	double yEu[100] ;
	double yEm[100] ;
	int m = sizeof(yEu)/sizeof(double);



	/*=========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	ode(func, yEu, t0, tf, h, Eu);
	ode(func, yEm, t0, tf, h, Em);
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("\n*************************************************");
	printf("\n            Output of Euler's method ");
	printf("\n*************************************************\n");
	printArr(yEu, "Eu",m);

	printf("\n*************************************************");
	printf("\n        Output of Modified Euler's method ");
	printf("\n*************************************************\n");
	printArr(yEm, "Em",m);


	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/	

	system("pause");
	return 0;
}