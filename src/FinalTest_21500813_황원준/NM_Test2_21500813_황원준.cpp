/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 10-06-2021
Modified         : 10-06-2021
Language/ver     : C++ in MSVS2019

Description      : Test2_yourID.cpp
----------------------------------------------------------------*/

// Include your library "myNM.h"
#include "../../include/myNM.h"
	

int main(int argc, char* argv[])
{
	// Q1
	printf("\n**************************************************");
	printf("\n|                   Question 1.                  |");
	printf("\n**************************************************\n");
	double Q1_T_arr[] = { 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
	double Q1_R_arr[] = { 12390, 8610, 6120, 4440, 3282, 2467, 1885, 1459, 1147, 912 };


	double Q1_a = 0, Q1_b = 0, Q1_R = 0, Q1_T = 0;
	// Add your code here 
	// Add your code here 
	// Add your code here 
	Matrix Q1_Tmat = arr2Mat(Q1_T_arr, 10, 1);
	Matrix Q1_Rmat = arr2Mat(Q1_R_arr, 10, 1);
	Matrix Q1z = createMat(2, 1);
	Matrix Q1_lnR = createMat(Q1_Rmat.rows, 1);
	for (int i = 0; i < Q1_Rmat.rows; i++)
		Q1_lnR.at[i][0]= log(Q1_Rmat.at[i][0]);
	Q1z = linearFit(Q1_Tmat, Q1_lnR);
	Q1_a = exp(Q1z.at[0][0]);
	Q1_b = Q1z.at[1][0];
	Q1_R = exp(FindValueLCF(Q1z, 20));
	Q1_T = FindValueLCFX(Q1z, log(7000));
	printf("\nanswer of a) Parameter a = %f, b = %f\n", Q1_a, Q1_b);
	printf("\nanswer of b) When T=20, R = %f\n", Q1_R);
	printf("\nanswer of c) When R=7k, T = %f\n\n", Q1_T);



	// Q2
	printf("\n**************************************************");
	printf("\n|                   Question 2.                  |");
	printf("\n**************************************************\n");
	double Q2_x_arr[] = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	double Q2_T_arr[] = { 473, 446.3, 422.6, 401.2, 382, 364.3, 348.0, 332.7, 318.1, 304.0, 290.1 };
	double Q2_q0 = 0, Q2_qL = 0;
	// Add your code here 
	// Add your code here 
	// Add your code here 
	double Q2_dtdx[11] = { 0 };
	double K = 240;
	gradient(Q2_x_arr, Q2_T_arr, Q2_dtdx, 11);
	Q2_q0 = -K*Q2_dtdx[0];
	Q2_qL = -K*Q2_dtdx[10];
	printf("\n qx(x=0) = %f, qx(x=L) = %f\n\n", Q2_q0, Q2_qL);



	// Q3
	printf("\n**************************************************");
	printf("\n|                   Question 3.                  |");
	printf("\n**************************************************\n");

	double Q3_I = 0;
	// Add your code here 
	// Add your code here 
	// Add your code here 
	double Q3y1[21];
	double Q3y2[21];
	double th0 = PI / 4;
	double dthdt = 0;

	sys2RK4(odeFunc_Q3, Q3y1, Q3y2, 0, 2, 0.1, th0, dthdt);


	Q3_I=integral13(Q3y1, 0, 2, 21);

	// Print your output for a)
	printArr(Q3y1, "a",21);
	// Print your output for b)
	printArr(Q3y2, "b", 21);
	printf("\nanswer of c) theta(t=0to2) = %f\n\n", Q3_I);



	// Q4
	printf("\n**************************************************");
	printf("\n|                   Question 4.                  |");
	printf("\n**************************************************\n");

	double Y_PC[5] = { 0 };
	double Y_EM[5] = { 0 };


	// Add your code here 
	// Add your code here 
	// Add your code here 

	double t0 = 0;
	double tf = 4;
	double y0 = 3;
	double h = 1;

	ode(Func_q4, Y_EM, t0, tf, h, y0, ODE_EM);
	odePC(Func_q4, Y_PC, t0, tf, h, y0);
	// Print answer for a) Y_PC 
	printArr(Y_EM, "Y_EM", 5);
	// Print answer for b) Y_EM - Y_PC
	printArr(Y_PC, "Y_PC", 5);

	system("pause");
	return 0;
}