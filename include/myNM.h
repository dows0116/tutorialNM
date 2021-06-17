/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"
#define ODE_EU 0
#define ODE_EM 1
#define ODE_RK2 2
#define ODE_RK4 3
#define PI 3.14159265368979323846264338327950288412
extern double Func_q4(const double t, const double y);
extern void odePC(double Func_q4(const double x, const double y), double Ypc[], double t0, double tf, double h, double y0);
extern double integral13(double y[], double a, double b, int n);
extern void odeFunc_Q3(const double t, const double Y[], double dYdt[]);
extern double FindValueLCFX(Matrix z, double R);

extern void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
extern void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
extern void odeFunc_mck(const double t, const double Y[], double dYdt[]);

extern double odeFunc_rc(const double t, const double v);
//ode method = 0 : Euler,   1 : Modified Euler
extern void ode(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0, int method);

// add parameter y function
extern double func(const double x, const double y);

// integral trap
extern double trapz(double x[], double y[], int m);

// simpson 13
extern double integral(double func(const double x), double a, double b, int n);

// rectangular
extern double IntegrateRect(double _x[], double _y[], int _m);

// Modify newton 
extern double newtonRaphsonfunc(double func(const double x), double dfunc(const double x), double x0, double tol); //  newtonRaphson method 함수 헤더로 불러오기
extern double func(const double x);
extern double dfunc(const double x);
// gradient to array
extern void  gradient(double x[], double y[], double dydx[], int m);
//gradient of function 
extern double myFunc(const double x);
//gradientfunc
extern Matrix gradientFunc(double func(const double x), Matrix xin);

//gradient
extern Matrix	gradient(Matrix _x, Matrix _y);

// Linear interpolation 
extern Matrix linearInterp(Matrix x, Matrix y, Matrix xq);


// totalerror and residure
extern double totalerror(Matrix x, Matrix y, Matrix z);
//Find value Fx for T at Linear interpolation
double FindValueLIP(Matrix yq, Matrix Tq, double T);
// Curvefitting
extern Matrix linearFit(Matrix _x, Matrix _y);

//Find value Fx for T and coefficient of Z
extern double FindValueLCF(Matrix z, double T);

// Find Condition number of A
extern double Cond(Matrix A);

// Find Eigen vector use QR decomposition
extern void QReig(Matrix A, Matrix Eig);
// QR decomposition House hold
extern void QRdecomp(Matrix A, Matrix I, Matrix Q, Matrix R);

////myfunc FOR QR
//extern void myfunc(Matrix X, Matrix FX);

// norm
extern	double norm(Matrix vec, int k);

// inverse matrix
extern	double inv(Matrix _A, Matrix _Ainv);

// scaled partial pivoting
extern	void SPpivot(Matrix _A, Matrix _P, int _k);

// LU solve function
extern	void	solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix vectorX);
// LUdecomposition
extern	void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);



//Gauss Elimination
extern  void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	void	backSub(Matrix _A, Matrix _b, Matrix _x);

// Apply forward-substitution
extern	void	fwdSub(Matrix _A, Matrix _b, Matrix _x);
// 과제 1
extern double bisectionNL(double _a0, double _b0, double _tol);  //  bisection method 함수 헤더로 불러오기
extern double newtonRhapson(double x0, double tol); //  newtonRaphson method 함수 헤더로 불러오기

extern double NewtonBisection(double _a0, double _b0, double x0, double _tol); // hybrid 함수 헤더로 불러오기
#endif