/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.h
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

typedef struct { 
	double** at;
	int rows, cols;
}Matrix;


extern void printArr2(double x[], double y[], const char* _name, int m);
// Print Array
extern void printArr(double x[], const char* _name,int m);

// Create a matrix from 1D-array
extern Matrix arr2Mat(double* _1Darray, int _rows, int _cols);

// Matrix multiply
extern Matrix MatMult(Matrix _A, Matrix _B);
// Matrix sum
extern Matrix Matsum(Matrix _A, Matrix _B);
// Matrix subtract
extern Matrix Matsub(Matrix _A, Matrix _B);
// initialization of Pivot elements
extern	void inipivot(Matrix _A);

// rank matrix
extern  double	rankMat(Matrix _A);
// paste Matrix fo Ax=b >> A|b Matrix
extern  Matrix	pasteMat(Matrix _A, Matrix _B);

//using namespace std;

// Create Matrix with specified size
extern	Matrix	createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void	freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix	txt2Mat(std::string _filePath, std::string _fileName);

//// Print matrix
extern	void	printMat(Matrix _A, const char* _name);


/// It is recommended to create the following functions.

// initialization of Matrix elements
extern	void	initMat(Matrix _A, double _val);

// Create matrix of all zeros
extern	Matrix	zeros(int _rows, int _cols);

// Create matrix of all ones
//extern	Matrix	ones(int _rows, int _cols);

// Create identity 
//extern	Matrix	eye(int _rows, int _cols);

// Create Transpose matrix
extern	Matrix	transpose(Matrix _A);

// Copy matrix
extern	Matrix	copyMat(Matrix _A); 

// Copy matrix Elements from A to B
extern	void	copyVal(Matrix _A, Matrix _B);


#endif