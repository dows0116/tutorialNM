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

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols) // 수업 실습시간에 사용된 어레이를 행렬로 바꾸는 함수. 
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
Matrix MatMult(Matrix _A, Matrix _B) { // 행렬 곱 함수. 앞 행렬의 열과 뒷행렬의 행값이 같은 경우에만 실행 가능.
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
Matrix Matsum(Matrix _A, Matrix _B) { // 행렬 합 함수. 앞행렬과 뒷행렬의 행과 열이 같아야지만 실행 가능.
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
Matrix Matsub(Matrix _A, Matrix _B) { // 행렬 차 함수 앞행렬과 뒷행렬의 행과 열이 같아야지만 실행 가능
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
Matrix	copyMat(Matrix _A) { //행렬 값을 넣음. 다만 주소값에 바로 저장되는 것이 아님.
	Matrix Out = createMat(_A.rows, _A.cols);
	initMat(Out, 0);
	for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
			Out.at[i][j] = _A.at[i][j];
	return Out;
}
// Copy matrix A to B
void	copyVal(Matrix _A, Matrix _B) { // 앞행렬을 뒤의 행렬에 값을 저장함. 주소값에 바로 저장된다.
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

void inipivot(Matrix _A) { //  넣은 행렬을 행렬 P ( == I )를 만드는 함수. 
	if (_A.rows != _A.cols) //  A가 square가 아닐때 오류로 실행이 안되게 함.
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
Matrix	pasteMat(Matrix _A, Matrix _B) // 두 행렬를 붙이는 함수. 행의 길이만 같으면 A뒤에 B를 붙히도록함.
{
	if (_A.rows != _B.rows) { // 두 행렬의 행길이가 다르면 오류 출력후 실행하지 않도록 함.
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'pasteMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_B.rows, _A.cols + _B.cols); //출력될 행렬 공간을 만듬. 
	for (int i = 0; i < _A.rows; i++) // A행렬을 앞에 복사.
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];

	for (int i = 0; i < _B.rows; i++) //B행렬을 뒤에 복사. 
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j + _A.cols] = _B.at[i][j];
	return Out;
}
// Rank Matrix
double	rankMat(Matrix _A) // A행렬의 rank를 판단함.
{
	Matrix temp_A = createMat(1, _A.cols);
	double ranklevel = _A.rows; // ranklevel은 행과 같음으로 행값으로 초기값을 설정.
	int ranktest = 0; 
	double divid = 0; 
	int zerorows = 0; // 0으로만 구성된 행이 있을경우를 위한 변수
	for (int r = 0; r < _A.rows; r++) 
	{
		zerorows = 0; 
		for (int k = 0; k < _A.cols; k++) // k는 행의 열 값을 찾기 위한 루프변수.
		{
			if (_A.at[r][k] != 0) // r행의 k번째열이 0이 아닐경우 다른 행과 비교할수 있음으로 0이 아닐때만 실행.
			{
				for (int i = r + 1; i < _A.rows; i++) 
				{
					ranktest = 0;
					divid = 0;
					if (_A.at[i][k] != 0) // a.[r+1][k]가 0이 아닐경우엔 a.[r][k]번째의 값과의 차이를 보고 그만큼을 곱하거나 나눠서 다른값들에 적용함. 
					{                     // a.[r+1][k]가 0일  경우엔 어차피  a.[r][k]번째의 값과 다른 값임 으로 rank level이 떨어지지 않음. 
						if (fabs(_A.at[r][k]) >= fabs(_A.at[i][k])) // 다른행의 두 값이 누가 더 크냐에 따라 곱하거나 나눌 값을 넣음.
							divid = _A.at[i][k] / _A.at[r][k];
						else
							divid = _A.at[i][k];
						for (int j = 0; j < _A.cols; j++)
						{
							if (_A.at[i][j] == 0 && _A.at[r][j] == 0)
								temp_A.at[0][j] = 0;
							else
								temp_A.at[0][j] = divid * _A.at[r][j] - _A.at[i][j];

							if (temp_A.at[0][j] != 0) // 같은 열의 행값이 차이가 하나라도 있을 경우 그 행들은 서로 다른 행임을 판단함.
								ranktest += 1;

						}
						if (ranktest == 0) // 같은 열의 행값이 차이가 하나라도 없을 경우 그 행들은 서로 일치하는 행임으로 랭크레벨을 낮춤.
							ranklevel -= 1;
					}

				}
				break;
			}
			else // a.[r][k]번째가 0이면 0행의 가능성이 있음으로 0행을 판단하는 변수를 +1
				zerorows += 1;
		}
		if (zerorows == _A.cols) // 행의 값이 전부 0일 경우 rank레벨을 낮춘다.
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