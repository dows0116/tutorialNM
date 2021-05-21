/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"


void printArr(double x[],const char* _name,int m) {
	printf("%s =\n", _name);
	for (int i = 0; i < m; i++)
		printf("%15.6f\n", x[i]);

	return;
}

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols) // ���� �ǽ��ð��� ���� ��̸� ��ķ� �ٲٴ� �Լ�. 
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}


Matrix	transpose(Matrix _A) {
	Matrix Out = createMat(_A.cols, _A.rows);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[j][i] = _A.at[i][j];
	
	return Out;

}

// Matrix Multiply
Matrix MatMult(Matrix _A, Matrix _B) { // ��� �� �Լ�. �� ����� ���� ������� �ప�� ���� ��쿡�� ���� ����.
	if (_A.cols != _B.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matrix Multiply' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}
	Matrix Out = createMat(_A.rows, _B.cols);
	initMat(Out, 0);


	for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++) 
			for (int k = 0; k < _A.cols; k++)
				Out.at[i][j] += _A.at[i][k] * _B.at[k][j];

	return Out;

}

// Matrix sum
Matrix Matsum(Matrix _A, Matrix _B) { // ��� �� �Լ�. ����İ� ������� ��� ���� ���ƾ����� ���� ����.
	if ((_A.cols != _B.cols) || (_A.rows != _B.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matrix Sum' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}
	Matrix Out = createMat(_A.rows, _B.cols);
	initMat(Out, 0);


	for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
				Out.at[i][j] += _A.at[i][j] + _B.at[i][j];

	return Out;
}


// Matrix subtract
Matrix Matsub(Matrix _A, Matrix _B) { // ��� �� �Լ� ����İ� ������� ��� ���� ���ƾ����� ���� ����
	if ((_A.cols != _B.cols) || (_A.rows != _B.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matrix Subtract' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}
	Matrix Out = createMat(_A.rows, _B.cols);
	initMat(Out, 0);


	for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
			Out.at[i][j] += _A.at[i][j] - _B.at[i][j];

	return Out;
}

// Copy matrix
Matrix	copyMat(Matrix _A) { //��� ���� ����. �ٸ� �ּҰ��� �ٷ� ����Ǵ� ���� �ƴ�.
	Matrix Out = createMat(_A.rows, _A.cols);
	initMat(Out, 0);
	for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
			Out.at[i][j] = _A.at[i][j];
	return Out;
}
// Copy matrix A to B
void	copyVal(Matrix _A, Matrix _B) { // ������� ���� ��Ŀ� ���� ������. �ּҰ��� �ٷ� ����ȴ�.
	if ((_A.cols != _B.cols) || (_A.rows != _B.rows)) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matrix Copy' function");
		printf("\n*************************************************\n");
		return;
	}
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_B.at[i][j] = _A.at[i][j];
	return ;
}

// initialization of Pivot elements

void inipivot(Matrix _A) { //  ���� ����� ��� P ( == I )�� ����� �Լ�. 
	if (_A.rows != _A.cols) //  A�� square�� �ƴҶ� ������ ������ �ȵǰ� ��.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return ;
	}
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++) {
			if (i == j)
				_A.at[i][j] = 1;
			else
				_A.at[i][j] = 0;
		}


	return ;

}


// create pasteMatrix
Matrix	pasteMat(Matrix _A, Matrix _B) // �� ��ĸ� ���̴� �Լ�. ���� ���̸� ������ A�ڿ� B�� ����������.
{
	if (_A.rows != _B.rows) { // �� ����� ����̰� �ٸ��� ���� ����� �������� �ʵ��� ��.
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'pasteMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_B.rows, _A.cols + _B.cols); //��µ� ��� ������ ����. 
	for (int i = 0; i < _A.rows; i++) // A����� �տ� ����.
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];

	for (int i = 0; i < _B.rows; i++) //B����� �ڿ� ����. 
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j + _A.cols] = _B.at[i][j];
	return Out;
}
// Rank Matrix
double	rankMat(Matrix _A) // A����� rank�� �Ǵ���.
{
	Matrix temp_A = createMat(1, _A.cols);
	double ranklevel = _A.rows; // ranklevel�� ��� �������� �ప���� �ʱⰪ�� ����.
	int ranktest = 0; 
	double divid = 0; 
	int zerorows = 0; // 0���θ� ������ ���� ������츦 ���� ����
	for (int r = 0; r < _A.rows; r++) 
	{
		zerorows = 0; 
		for (int k = 0; k < _A.cols; k++) // k�� ���� �� ���� ã�� ���� ��������.
		{
			if (_A.at[r][k] != 0) // r���� k��°���� 0�� �ƴҰ�� �ٸ� ��� ���Ҽ� �������� 0�� �ƴҶ��� ����.
			{
				for (int i = r + 1; i < _A.rows; i++) 
				{
					ranktest = 0;
					divid = 0;
					if (_A.at[i][k] != 0) // a.[r+1][k]�� 0�� �ƴҰ�쿣 a.[r][k]��°�� ������ ���̸� ���� �׸�ŭ�� ���ϰų� ������ �ٸ����鿡 ������. 
					{                     // a.[r+1][k]�� 0��  ��쿣 ������  a.[r][k]��°�� ���� �ٸ� ���� ���� rank level�� �������� ����. 
						if (fabs(_A.at[r][k]) >= fabs(_A.at[i][k])) // �ٸ����� �� ���� ���� �� ũ�Ŀ� ���� ���ϰų� ���� ���� ����.
							divid = _A.at[i][k] / _A.at[r][k];
						else
							divid = _A.at[i][k];
						for (int j = 0; j < _A.cols; j++)
						{
							if (_A.at[i][j] == 0 && _A.at[r][j] == 0)
								temp_A.at[0][j] = 0;
							else
								temp_A.at[0][j] = divid * _A.at[r][j] - _A.at[i][j];

							if (temp_A.at[0][j] != 0) // ���� ���� �ప�� ���̰� �ϳ��� ���� ��� �� ����� ���� �ٸ� ������ �Ǵ���.
								ranktest += 1;

						}
						if (ranktest == 0) // ���� ���� �ప�� ���̰� �ϳ��� ���� ��� �� ����� ���� ��ġ�ϴ� �������� ��ũ������ ����.
							ranklevel -= 1;
					}

				}
				break;
			}
			else // a.[r][k]��°�� 0�̸� 0���� ���ɼ��� �������� 0���� �Ǵ��ϴ� ������ +1
				zerorows += 1;
		}
		if (zerorows == _A.cols) // ���� ���� ���� 0�� ��� rank������ �����.
			ranklevel -= 1;
	}
	return ranklevel;
}


// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}		

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}

// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{    
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;

}

// Create matrix of all zeros
Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);


	return Out;
}