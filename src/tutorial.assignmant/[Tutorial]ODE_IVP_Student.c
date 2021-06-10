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

//#define ODE_EU 0
//#define ODE_EM 1
//#define ODE_RK2 2
//#define ODE_RK4 3
//
////  PI is defined in  myNM.h
//#define PI 3.14159265368979323846264338327950288412
//
//// Problem1: Single equation of 1st order ODE
//double odeFunc_rc(const double t, const double v);
//
//// Problem2: Single equation of 2nd order ODE
//void odeFunc_mck(const double t, const double Y[], double dYdt[]);
//
//
//// Single Equation : odeEM, odeEU
//void odeEU(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeEM(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//
//// 2nd order Equations : sys2RK2, sys2RK4
//void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
//void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);


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
	ode(odeFunc_rc, y_EU, a, b, h, Eu);
	ode(odeFunc_rc, y_EU, a, b, h, Em);
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

	////Parameter Definitions
	//double t0 = 0;
	//double tf = 1;
	//h = 0.01;
	//N = (tf - t0) / h + 1;
	//double y[200] = { 0 };
	//double v[200] = { 0 };

	//// Initial values
	//double y0 = 0;
	//v0 = 0.2;

	//// ODE solver: RK2
	//sys2RK2(odeFunc_mck, y, v, t0, tf, h, y0, v0);

	/////////////////////////////////////////////////////////////////
	//// Exercise 3: Create the standard form  for RK4 for 2nd order	
	////sys2RK4(odeFunc_mck, y, v, t0, tf, h, y0, v0);
	//   


	//// Print outputs
	//printf("/*---------------------------*/\n");
	//printf("/ 2nd Order ODE : MCK example /\n");
	//printf("/*---------------------------*/\n");
	//printf(" - Total number of data N=%d \n", N);
	//for (int i = 0; i < N; i++)
	//	printf("t= %f\ty= %f\tv= %f\n", t0 + i * h, y[i], v[i]);
	//printf("\n\n");

	//	
	//// Copy and paste the output in MATLAB and PLOT
	//for (int i = 0; i < N; i++)
	//	printf("%f\t%f\t%f\n", t0 + i * h, y[i], v[i]);


	system("pause");
	return 0;
}

//// Gradient function for ODE - 1st order single eq.
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
//
//void odeFunc_mck(const double t, const double Y[], double dYdt[])
//{
//	double m = 1;
//	double c = 7;
//	double k = 6.9;
//	double f = 5;
//
//	double Fin = 2 * cos(2 * PI * f * t);
//
//	dYdt[0] = Y[1];	
//
//	// EXERCISE: MODIFY HERE
//	// HINT;   zdot= (-k*Y - c*Z + Fin)/m;
//	dYdt[1] = 0;
//}
//
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
//void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
//{
//	double C1 = 0.5;
//	double C2 = 0.5;
//	double alpha = 1;
//	double beta = alpha;  // alpha=beta
//		
//	int N = (tf - t0) / h + 1;
//	double ti = t0;
//	double y_EU;
//	double K1=0, K2=0;
//
//	// Initialization 
//	y[0] = y0;
//
//	for (int i = 0; i < N - 1; i++) 
//	{
//		// First slope  
//		// K1=?
//		
//
//		// Second slope  
//		// K2=? 
//		
//
//		// Update 
//		y[i + 1] = y[i] + (C1*K1 + C2 * K2) * h;
//		ti += h;
//	}
//}
//
//
//void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
//{
//// EXERCISE
//	
//}
//
//
///*---------------------------------------------------------------------------------------------------------------------------*/
//
//// ODE RK2:  one of 2nd order ODE <--> two of 1st order ODE
//void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
//{
//	int N = (tf - t0) / h + 1;
//	double ti = t0;
//
//	double K1[2] = { 0 };
//	double K2[2] = { 0 };
//	double Yin[2] = { 0 };
//	double K1_y1=0, K1_y2=0, K2_y1=0, K2_y2=0;
//
//
//	// Initial condition
//	y1[0] = y1_init;
//	y2[0] = y2_init;
//
//	for (int i = 0; i < N - 1; i++) {
//		
//		// Slope 1 : K1
//		Yin[0] = y1[i];		// z
//		Yin[1] = y2[i];		// dzdt		
//		odeFunc_sys2(ti, Yin, K1);
//		K1_y1 = K1[0];
//		K1_y2 = K1[1];
//
//		// Slope 2 : K2
//		// ADD CODE HERE
//
//		// Update
//		y1[i + 1] = y1[i] + (K1_y1 + K2_y1) * 0.5 * h;		
//		//y2[i + 1] = ? 		
//		ti += h;
//	}
//}
//
//
//// Classical RK4
//void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
//{
//
//	// EXERCISE
//
//}

