/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park, YKKIM
Created          : 2021-06-03
Modified         : 2021-06-03  by YKKIM
Language/ver     : C++ in MSVS2017

Description      : [Tutorial]ODE_IVP_student.c
-------------------------------------------------------------------------------*/


#include "../../include/myNM.h"
#include <stdio.h>
#include <math.h>




int main(int argc, char* argv[])
{

	/*-------------------------------------------------------------------*/
	// Single of 1st Order ODE
	/*-------------------------------------------------------------------*/

	//Parameter Definitions
	double a = 0;
	double b = 0.1;
	double h = 0.001;
	unsigned int N = (b - a) / h + 1;
	double y_EU[200] = { 0 };				//Cannot use y_EU[N]
	double y_EM[200] = { 0 };
	double y_RK2[200] = { 0 };
	double y_RK4[200] = { 0 };

	// Initial value
	double v0 = 0;

	// ODE solver
	ode(func, y_EU, a, b, h, v0,Eu);
	ode(func, y_EM, a, b, h, v0,Em);
	ode(func, y_RK2, a, b, h, v0,ODE_RK2);
	ode(func, y_RK4, a, b, h, v0, ODE_RK4);
	//odeEU(odeFunc_rc, y_EU, a, b, h, v0);
	//odeEM(odeFunc_rc, y_EM, a, b, h, v0);

	
	////////////////////////////////////////////////	
	// Exercise 1: Create a general form for RK2
	//odeRK2(odeFunc_rc, y_RK2, a, b, h, v0);

	//// Exercise 2: Create the standard form  for RK4
	//odeRK2(odeFunc_rc, y_RK2, a, b, h, v0);
	////////////////////////////////////////////////
		
	
	// Print outputs
	printf("/*-----------------------*/\n");
	printf("/ Single of 1st Order ODE /\n");
	printf("/*-----------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("t= %f\tyEU= %f\tyEM= %f\tyRK2= %f\tyRK4= %f\n", a + i * h, y_EU[i], y_EM[i], y_RK2[i], y_RK4[i]);
	printf("\n");



	///*-------------------------------------------------------------------*/
	//// 2nd Order ODE : MCK example
	///*-------------------------------------------------------------------*/

	//Parameter Definitions
	double t0 = 0;
	double tf = 1;
	h = 0.01;
	N = (tf - t0) / h + 1;
	double y[200] = { 0 };
	double v[200] = { 0 };

	// Initial values
	double y0 = 0;
	v0 = 0.2;

	// ODE solver: RK2
	sys2RK2(odeFunc_mck, y, v, t0, tf, h, y0, v0);

	///////////////////////////////////////////////////////////////
	// Exercise 3: Create the standard form  for RK4 for 2nd order	

	   


	// Print outputs
	printf("/*---------------------------*/\n");
	printf("/ 2nd Order ODE : MCK example /\n");
	printf("/*---------------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("t= %f\ty= %f\tv= %f\n", t0 + i * h, y[i], v[i]);
	printf("\n\n");

	sys2RK4(odeFunc_mck, y, v, t0, tf, h, y0, v0);
	printf("/*---------------------------*/\n");
	printf("/ 2nd Order ODE : RK4 /\n");
	printf("/*---------------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("t= %f\ty= %f\tv= %f\n", t0 + i * h, y[i], v[i]);
	printf("\n\n");

	// Copy and paste the output in MATLAB and PLOT
	for (int i = 0; i < N; i++)
		printf("%f\t%f\t%f\n", t0 + i * h, y[i], v[i]);


	system("pause");
	return 0;
}

// Gradient function for ODE - 1st order single eq.
//double odeFunc_rc(const double t, const double v) {
//
//	double tau = 1;
//	double T = 1 / tau;
//	double f = 10;
//	double Vm = 1;
//	double omega = 2 * PI * f;
//
//	return  -T * v + T * Vm * cos(omega * t);
//}
//





/*---------------------------------------------------------------------------------------------------------------------------*/


//void odeEU(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
//{
//	int N = (tf - t0) / h + 1;
//	double ti = t0;
//	double slope;
//
//	y[0] = y0;
//
//	for (int i = 0; i < N - 1; i++) {
//		slope = odeFunc(ti, y[i]);
//		y[i + 1] = y[i] + slope * h;
//		ti += h;
//	}
//}
//
//void odeEM(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
//{
//	int N = (tf - t0) / h + 1;
//	double ti = t0;
//	double y_EU;
//	double K1, K2;
//
//	y[0] = y0;
//	for (int i = 0; i < N - 1; i++) {
//
//		// First slope
//		K1 = odeFunc(ti, y[i]);
//
//		// Second slope
//		y_EU = y[i] + K1 * h;
//		K2 = odeFunc(ti + h, y_EU);
//
//		// Update 
//		y[i + 1] = y[i] + (K1 + K2) * 0.5 * h;
//		ti += h;
//	}
//}
//
//// EXERCISE


//void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
//{
//// EXERCISE
//	
//}
//
//
///*---------------------------------------------------------------------------------------------------------------------------*/
//
// ODE RK2:  one of 2nd order ODE <--> two of 1st order ODE
