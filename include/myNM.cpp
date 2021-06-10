/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"

double odeFunc_rc(const double t, const double v) {

	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double omega = 2 * PI * f;

	return  -T * v + T * Vm * cos(omega * t);
}


double func(const double x, const double y) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * M_PI * f;
	double dvdt = -1*T * y + T * Vm * cos(w * x);
	return dvdt;
}


void ode(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0 , int method) {
	int m = (tf - t0) / h;
	double xt = t0;
	y[0] = 0;
	double slope1 = 0;
	double slope2 = 0;
	double yE = 0;
	double py = 0;
	if (method == Eu) { 	//Euler
		for (int i = 0; i < m - 1; i++) {
			xt = xt + h;
			y[i + 1] = y[i] + func(xt, y[i]) * h;
		}
	}
	else if (method == Em) { 	//Modified Euler 
		for (int i = 0; i < m - 1; i++) {
			slope1 = func(xt, y[i]); // 현재 값인 x[i]와 y[i]로 slope1 값을 구하면서 yE[i+1]를 구한다.
			yE = y[i] + func(xt, y[i]) * h;
			xt = xt + h; // i+1값 
			slope2 = func(xt, yE); // 기존 오일러로 구한 yE[i+1] 로 구한 slope 2값
			y[i + 1] = y[i] + (slope1 + slope2) * h / 2;
		}
	}
	else if (method == ODE_RK2) {
		double C1 = 0.5;
		double C2 = 0.5;
		double alpha = 1;
		double beta = alpha;  // alpha=beta

		int N = (tf - t0) / h + 1;
		double ti = t0;
		double y_EU;
		double K1 = 0, K2 = 0;

		// Initialization 
		y[0] = y0;

		for (int i = 0; i < N - 1; i++)
		{
			// First slope  
			// K1=?
			K1 = odeFunc_rc(ti, y[i]);

			// Second slope  
			// K2=? 
			y_EU = y[i] + beta * K1 * h;
			K2 = odeFunc_rc(ti + alpha * h, y_EU);

			// Update 
			y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;
			ti += h;
		}
	}
	else if (method == ODE_RK4) {

		double a = 0.5;

		int N = (tf - t0) / h + 1;
		double ti = t0;
		double y_EU;
		double K1 = 0, K2 = 0, K3 = 0, K4 = 0;

		// Initialization 
		y[0] = y0;

		for (int i = 0; i < N - 1; i++)
		{
			// First slope  
			// K1=?
			K1 = odeFunc_rc(ti, y[i]);

			// Second slope  
			// K2=? 
			y_EU = y[i] + a * K1 *h ;
			K2 = odeFunc_rc(ti+a*h, y_EU);

			y_EU = y[i] + a * K2 * h;
			K3 = odeFunc_rc(ti + a * h, y_EU);

			y_EU = y[i] + a * K3 * h;
			K4 = odeFunc_rc(ti + h, y_EU);
			// Update 
			y[i + 1] = y[i] + (K1 * 2*K2 + 2*K3 * K4) * h/6;
			ti += h;
		}
	}




	return;

}


double func(const double x) // f(x) 수식 값 입력
{
	return x*x*x;
}

double dfunc(const double x) // // f'(x) 수식 값 입력
{
	return 3*x*x;
}

double newtonRaphsonfunc(double func(const double x), double dfunc(const double x), double x0, double tol) // newtonRaphson method 함수 정의
{
	int i = 0; // iteration 값 정의
	double xn = 0; // xn 변수 정의 정확히는 xn+1
	double x = x0; // 초기 x 값을 넣는 변수, xn 
	int Nmax = 1000; // 반복 횟수는 1000번까지 iIteration 횟수 
	double ep = 0; // tolerance 비교 값을 위한 변수
	if (dfunc(x) != 0) { // f'(x)가 0이 되면 xn이 발산하기 때문에 0이 아닐때만 실행하도록 함
		for (i; i < Nmax; i++) {
			xn = -func(x) / dfunc(x) + x;
			ep = fabs((xn - x) / x);
			printf("K:%d \t", i);
			printf("X(n): %0.10f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			if (ep <= tol) // tolerance 이하가 되었을때 반복문을 정지
				break;
			else if (fabs(xn) > 10000 * fabs(x)) {  // 해를 찾지 못하고 급히 발산할 경우 발산 에러 출력
				printf("Error : Diverging error\n ");
				break;
			}
			x = xn;
		}
		if (i == Nmax && ep > tol)
			printf("Error : infinite loop error\n "); // 해를 찾지 못하고 MaxItr에 도달 했을 경우 무한루프 에러 출력
	}
	else  // f'(x) = 0일 경우 해를 찾을수 없기에 오류 출력 후 정지 
		printf("Error : f'(x) is 0 , you can't find solution by using Newton-Raphson Method\n");

	return xn;
}



void  gradient(double x[], double y[], double dydx[], int m) {
	double h = x[1] - x[0];
	if ( m < 2) { //  2보다 작아 함수를 잘 실행시킬수 없음으로 오류 출력
		printf("Error!! : 입력값의 길이 2보다 작음");
		return ;
	}
	else if (m == 2) {// 길이가 2일 경우에  2포인트 포워드 및 벡워드 사용
		dydx[0] = (y[1] - y[0]) / (h);
		dydx[1] = (y[1] - y[0]) / (h);
	}
	else {
		dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

		for (int i = 1; i < m-1; i++)
			dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);


		dydx[m-1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m-1]) / (2 * h);
	}
	return ;
	
}

double myFunc(const double x) {
	return  x * x * x;
}

Matrix gradientFunc(double myFunc(const double x), Matrix xin) {
	Matrix y = createMat(xin.rows, 1);
	Matrix fdx = createMat(xin.rows, 1);
	for (int i = 0; i < y.rows; i++)
		y.at[i][0] = myFunc(xin.at[i][0]);

	fdx = gradient(xin, y);

	return fdx;
}


Matrix	gradient(Matrix _x, Matrix _y) {
	Matrix fdx = createMat(_x.rows, 1);
	int mx = _x.rows;
	int m = mx - 1;
	int my = _y.rows;
	double h = _x.at[1][0] - _x.at[0][0];
	if ((mx != my) || (mx < 2)) { // x와 y의 길이가 안맞거나 2보다 작아 함수를 잘 실행시킬수 없음으로 오류 출력
		printf("Error!! : 입력값의 길이가 맞지 않음");
		return createMat(0, 0);
	}
	else if ((mx == my) && (mx == 2)) { // 길이가 2일 경우에  2포인트 포워드 및 벡워드 사용
		fdx.at[0][0] = (_y.at[1][0] - _y.at[0][0]) / (h);
		fdx.at[1][0] = (_y.at[1][0] - _y.at[0][0]) / (h);
	}
	else {
		fdx.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);

		for (int i = 1; i < m; i++)
			fdx.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);


		fdx.at[m][0] = (_y.at[m - 2][0] - 4 * _y.at[m - 1][0] + 3 * _y.at[m][0]) / (2 * h);
	}
	return fdx;
}


double FindValueLIP(Matrix yq, Matrix Tq, double T) {   // 특정 T에 대한 추정 값 yq를 찾는 함수
	double Fx = 0;
	int a = 0;
	int no = 0;
	for (int i = 0; i < Tq.rows; i++) { // T의 Tq=xq 위치를 찾고 없는 경우 no값을 누적시켜 오류를 출력하도록함.
		if (Tq.at[i][0] == T)
			a = i;
		else
			no++;
	}
	if (no == Tq.rows) {
		printf("\n*************************************************");  // 추정값으로 넣을  T가 Tquery 행렬 안에 없을 경우 오류 출력
		printf("\n ERROR!!: 추정시킬 값이 yq 안에 존재하지 않음 ");
		printf("\n*************************************************\n");
		return Fx = 0;
	}

	Fx = yq.at[a][0];

	return Fx;
}


Matrix linearInterp(Matrix x, Matrix y, Matrix xq) { // 리니어 인터폴레이션 함수 
	Matrix yq = createMat(xq.rows, 1);
	int mx = x.rows;
	int my = y.rows;
	if ((mx != my) || (mx < 2)) { // x와 y의 길이가 안맞거나 2보다 작아 함수를 잘 실행시킬수 없음으로 오류 출력
		printf("Error!! : 입력값의 길이가 맞지 않거나 2보다 작아 Linear Interpolation 을 실행 시킬 수 없음.");
		return zeros(yq.rows,yq.cols); 
	}
	int i = 0;
	for (int j = 0; j < xq.rows; j++) {


		if ((xq.at[j][0] >= x.at[i + 1][0]) || (xq.at[j][0] < x.at[i][0])) { // xq가 x에서 어디 위치하는지 찾고 그에 해당하는 yq의 직선을 찾기 위함.
			for (int k = i; k < x.rows-1; k++) {
				if ((xq.at[j][0] < x.at[k + 1][0]) && (xq.at[j][0] >= x.at[k][0])) {
					i = k;
					break;
				}
			}
		}

		// 위치를 찾은 리니어 함수에 값에 xq값을 넣어 yq값을 도출
		yq.at[j][0] = y.at[i][0] * (xq.at[j][0] - x.at[i + 1][0]) / (x.at[i][0] - x.at[i + 1][0]) + y.at[i + 1][0] * (xq.at[j][0] - x.at[i][0]) / (x.at[i + 1][0] - x.at[i][0]); 

	}

	return yq;

}

double FindValueLCF(Matrix z, double T) { // 커브 피팅에서 원하는 입력값 T의 추정값을 찾는 함수
	double Fx= 0;
	Fx = z.at[0][0] + z.at[1][0] * T;

	return Fx;
}

Matrix	linearFit(Matrix _x, Matrix _y) { // 리니어 커브 피팅 함수
	int mx = _x.rows;
	int my = _y.rows;
	double a1 = 0;
	double a0 = 0;
	if ((mx != my) || (mx < 2)) { //오류 출력
		printf("Error!! : 입력값의 길이가 맞지 않거나 2보다 작아 Linear curve fitting 을 실행 시킬 수 없음.");
		a1 = 0; a0 = 0;
	}
	else {

		double Sx = 0;
		double Sxx = 0;
		double Sxy = 0;
		double Sy = 0;
		int m = mx;
		for (int k = 0; k < m; k++) {
			Sx = Sx + _x.at[k][0];
			Sxx = Sxx + _x.at[k][0] * _x.at[k][0];;
			Sxy = Sxy + _x.at[k][0] * _y.at[k][0];
			Sy = Sy + _y.at[k][0];
		}
		a1 = (m * Sxy - Sx * Sy) / (m * Sxx - Sx * Sx);
		a0 = (Sxx * Sy - Sxy * Sx) / (m * Sxx - Sx * Sx);
	}
	Matrix z = createMat(2, 1);
	z.at[0][0] = a0;
	z.at[1][0] = a1;
	return z;
}

double Cond(Matrix A) { //컨디션 함수 A의 컨디션 넘버를 찾기 위한 함수


	Matrix At = createMat(A.cols, A.rows);     
	Matrix AtA = createMat(A.cols, A.cols);
	Matrix eig = createMat(AtA.rows, 1);
	Matrix Atemp = createMat(A.cols, A.cols);
	At = transpose(A);

	Atemp = MatMult(At, A);
	copyVal(Atemp, AtA);
	QReig(AtA,eig);   // AT*A의 컨디션 넘버를 QR factorization 으로 아이겐 벡터를 찾는다. 
	//printMat(eig, "eig");
	double normA = norm(eig, 3);  // 찾은 아이겐 벡터에서 최대값은 AT*A 에 저장한다.
	double normAinv = 1/ norm(eig, 4); //  AT*A의 역함수의 최대값은 아이겐 벡터의 0이 아닌  최소값을 1로 나눈 것이다. 
	double ans = sqrt(normA) * sqrt(normAinv);
	
	return ans;
}

void QReig(Matrix A, Matrix Eig) { // QR decopmsition 을 통해 QR을 구하고 A값이 similat matrix가 될 때 까지 반복하여 아이겐 벡터를 찾는다.
	Matrix Atemp = createMat(A.rows, A.cols);
	copyVal(A, Atemp);
	Matrix QTA = createMat(A.rows, A.cols);
	Matrix AQ = createMat(A.rows, A.cols);
	Matrix Q = createMat(A.rows, A.cols);
	Matrix R = createMat(A.rows, A.cols);
	Matrix Qt = createMat(Q.rows, Q.cols);
	Matrix I = createMat(Q.rows, Q.cols);
	copyVal(A, R);
	copyVal(I, Q);
	inipivot(I);
	int nmax = 100;
	for (int k = 0; k < nmax; k++) { // A가 similat matrix가 될 때 까지 반복하여 QRdecomp를 실행한다. 
		QRdecomp(Atemp,I,Q,R);
		Qt = transpose(Q);
		QTA = MatMult(Qt, Atemp);
		Atemp = MatMult(QTA, Q);
	}
	for (int i = 0; i < Eig.rows; i++)
		Eig.at[i][0] = Atemp.at[i][i];
	return;
}
// QRdecomp
void QRdecomp(Matrix A, Matrix I, Matrix Q, Matrix R) {

	// 조건식 스퀘어 , 랭크 a = n 이여야함.
	if (A.rows != A.cols) //  A가 square가 아닐때 오류로 실행이 안되게 함.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return ;
	}
	copyVal(A, R);
	copyVal(I, Q);
	Matrix c = createMat(R.rows, 1);
	Matrix v = createMat(R.rows, 1);
	Matrix vt = createMat(R.rows, 1);
	Matrix e = createMat(R.rows, 1);
	Matrix H = createMat(R.rows, R.cols);
	Matrix vvt = createMat(I.rows, I.cols);
	double normc = 0;
	double normv = 0;

	Matrix Rtemp = createMat(R.rows, R.cols);
	Matrix Qtemp = createMat(Q.rows, Q.cols);
	for (int k = 0; k < R.cols-1 ; k++) { 
		for (int i = 0; i < R.rows; i++)  
			c.at[i][0] = R.at[i][k];
		if(k>0)
		for (int i = 0; i < k; i++)
			c.at[i][0] = 0;

		normc = norm(c, 2);
		for (int i = 0; i < R.rows; i++)
			e.at[i][0] = 0;

		if(c.at[k][0]>=0)
			e.at[k][0] = normc;
		else 
			e.at[k][0] = -normc;

		v = Matsum(c, e);  // 여기까지 v = c + ||c||*e  구하는 알고리즘
		normv = norm(v, 2); // vtv는 norm2의 제곱이므로 v벡터의 norm2를 구함.
		vt = transpose(v);
		vvt = MatMult(v, vt);
		if (normv == 0) {
			printf("\n*************************************************"); // 0으로 나뉘어지는 경우 함수 중단.
			printf("\n ERROR!!: Division by zero : norm v = 0 ");
			printf("\n*************************************************\n");
			return;
		}
			
		for (int i = 0; i < vvt.rows; i++)
			for (int j = 0; j < vvt.cols; j++)
				vvt.at[i][j] = -2 * vvt.at[i][j]/(normv*normv);
		H = Matsum(I, vvt);
		Qtemp = MatMult(Q, H);
		Rtemp = MatMult(H, R); // 
		copyVal(Qtemp, Q);

		copyVal(Rtemp, R);

	}
	return;
}

// norm
double norm(Matrix vec, int k) { // norm 함수 1일때는 절대값을 모두 더한것, 2일때는 크기를 구하는것, 3은 최대값 , 4는 0이 아닌 최소값.
	double sum = 0;
	double max = 0;
	double min = 0;
	Matrix rsum = createMat(1, 1);
	initMat(rsum, 0);
	if (k == 1) {// norm 1 
		for (int i = 0; i < vec.rows; i++)
			sum += fabs(vec.at[i][0]);
		return sum;
	}
	else if (k == 2) {// norm 2
		Matrix vect = transpose(vec);

		rsum = MatMult(vect, vec);
		sum = sqrt(rsum.at[0][0]);
		return sum;
	}
	else if (k == 3) { // norm inf
		for (int i = 0; i < vec.rows; i++) {
			sum = fabs(vec.at[i][0]);
			if (sum > max)
				max = sum;
		}
		return max;
	}
	else if (k == 4) { // norm inf
		min = fabs(vec.at[0][0]);
		for (int i = 1; i < vec.rows; i++) {
			sum = fabs(vec.at[i][0]);
			if ((sum < min) && sum != 0)
				min = sum;
		}
		return min;
	}
}



//// myfunc
//void myfunc(Matrix X, Matrix FX) {
//	double F1 = 0;
//	for (int i = 0; i < X.rows; i++) {
//		FX.at[i][0] = 4 * X.at[i][0]; // 수식대입
//	}
//	// F1 = 수식 x,y 대입
//
//	return ;
//}





// inverse matrix
double inv(Matrix _A, Matrix _Ainv) { // 역행렬 함수, LU를 이용했고, Ainv의 한 열마다 LUsolve하여 값을 구해냄.

	Matrix _Atemp = createMat(_A.rows, _A.cols);
	Matrix _Ainvtemp = createMat(_Ainv.rows, _Ainv.cols);
	Matrix L = createMat(_A.rows, _A.cols);
	Matrix U = createMat(_A.rows, _A.cols);
	Matrix P = createMat(_A.rows, _A.cols);
	Matrix I = createMat(_A.rows, _A.cols);
	Matrix b = createMat(_A.rows, 1);
	Matrix X = createMat(_A.rows, 1);
	Matrix _Pvec = createMat(P.rows, 1);
	initMat(L, 0);
	initMat(U, 0);
	initMat(b, 0);
	inipivot(P);
	inipivot(I);
	copyVal(_A, _Atemp);
	int org = 0;
	// 조건식 스퀘어 , 랭크 a = n 이여야함.
	if (_A.rows != _A.cols) //  A가 square가 아닐때 오류로 실행이 안되게 함.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return 0;
	}
	double rankA = rankMat(_Atemp);
	double n = _A.rows;
	if ((rankA != n)) // A의 랭크가  equation n보다 낮을경우 해를 구할수없음을 오류로 출력
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f < n %0.f : infinite solution", rankA, n);
		printf("\n*************************************************\n");
		return 0;
	}
	printf("\n*************************************************");
	printf("\n            Inverse Matrix A  ");
	LUdecomp(_Atemp, L, U, P);
	for (int k = 0; k < _A.cols; k++) { // Ainv의 각 열마다를 LU solve함수를 써서 구했음.
		for (int i = 0; i < I.rows; i++)
			b.at[i][0] = I.at[i][k];
		solveLU(L, U, P, b, X);
		for (int i = 0; i < I.rows; i++)
			_Ainv.at[i][k] = X.at[i][0]; // solve로 구한 X값을 해당하는 열의 Ainv 임시값에 저장함.
	}
	return 0;
}

// LU solve function
void	solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix vectorX) {
	if ( ((_L.rows != _U.rows) || (_L.cols != _U.cols)) || ((_U.rows != _P.rows) || (_U.cols != _P.cols))|| (_P.rows != _b.rows) || ((_b.rows != vectorX.rows) || (_b.cols != vectorX.cols)))
	{ // 인풋과 아웃풋의 행렬이 잘못되었을때 오류 출력
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at L & U & P & b & x LU solve function");
		printf("\n*************************************************\n");
		return;
	}
	Matrix y = createMat(_U.rows, vectorX.cols); // y는 Ux 값
	Matrix _btemp = createMat(_b.rows, _b.cols);
	initMat(y, 0);
	_btemp = MatMult(_P, _b); // PA = LU = Pb 이므로 피벗을 b에도 적용함. 
	fwdSub(_L, _btemp, y); //  L(Ux=y)= Pb의 값을 포워드로 구함.
	backSub(_U, y, vectorX); // Ux = y이므로 백서브로 구함. 

	return ;
}



// LUdecomposition
void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P) {
	if (_A.rows != _A.cols) // 과제의 조건인 A가 square가 아닐때 오류로 실행이 안되게 함.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return;
	}
	if( ((_A.rows != _L.rows) || (_A.cols != _L.cols)) || ((_L.rows != _U.rows) || (_L.cols != _U.cols)) || ((_U.rows != _P.rows) || (_U.cols != _P.cols)) ) 
	{ // 인풋과 아웃풋의 행렬이 잘못되었을때 오류 출력
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A & L & U & P  LUdecomposition ");
		printf("\n*************************************************\n");
		return;
	}

	double rankA = rankMat(_A);
	double n = _A.rows;
	if ((rankA != n)) // A의 랭크가  equation n보다 낮을경우 해를 구할수없음을 오류로 출력
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f < n %0.f : infinite solution", rankA, n);
		printf("\n*************************************************\n");
		return;
	}
	printf("\n             LU decomposition");
	printf("\n*************************************************\n");

	Matrix _Atemp = createMat(_A.rows, _A.cols);
	Matrix _Ltemp = createMat(_L.rows, _L.cols);
	Matrix _Lsub = createMat(_L.rows, _L.cols);
	Matrix _LAtemp = createMat(_A.rows, _A.cols);
	Matrix _Lorig = createMat(_A.rows, _A.cols);
	Matrix _Pcopy = createMat(_P.rows, _P.cols);
	_Atemp = copyMat(_A);;
	_Ltemp = copyMat(_L);
	initMat(_Lsub, 0);
	inipivot(_Pcopy);
	Matrix _Ptemp = createMat(_P.rows, _P.cols);
	Matrix _Pvec = createMat(_P.rows, 1);
	int org = 0;
	for (int k = 0; k < _A.cols-1; k++) { 
	//	printf("\n At k = %d \n\n", k);
		inipivot(_Ptemp);
		SPpivot(_Atemp,_Ptemp, k); //Scaled pivot으로 피벗된 P를 구함, 
		_Ltemp = MatMult(_Ptemp, _Ltemp); // L P A 에 Scaled pivot 을 곱하여 피벗된 값을 저장한다.
		_Pcopy = MatMult(_Ptemp, _Pcopy);
	//	printMat(_Pcopy, "P");
		_Atemp = MatMult(_Ptemp, _Atemp);
		_Lsub = copyMat(_Ltemp);
		for (int i = k; i < _A.rows - 1; i++) {
			_Ltemp.at[i + 1][k] = _Atemp.at[i + 1][k] / _Atemp.at[k][k];
			if (_Atemp.at[k][k]==0) 
			{
				printf("\n*************************************************"); // 0으로 나뉘어지는 경우 함수 중단.
				printf("\n ERROR!!: Division by zero at '%d' LU decomposition ", k);
				printf("\n*************************************************\n");
				return;
			}
		}

		_Lsub = Matsub(_Ltemp, _Lsub); // 빼야되는 열의 L값만을 추출
		_LAtemp = MatMult(_Lsub, _Atemp); // LA를 저장
		_Atemp = Matsub(_Atemp, _LAtemp); // 기존 A값에 LA값을 빼서 U행렬을 만듬
		//printMat(_Atemp, "U");

		for (int i=0;i<_Pcopy.rows;i++) // 피벗 되기전의 L을 구하기 위한 피벗이 어떻게 바뀌었는 지를 순서를 숫자로 벡터에 저장.
			for (int j = 0; j < _Pcopy.cols; j++) {
				if (_Pcopy.at[i][j] == 1)
					_Pvec.at[i][0] = j;
			}
		
		for (int i = 0; i < _Lorig.rows; i++) // 순서 벡터를 이용하여 각 행을 알맞는 순서로 바꾸어서 L original 값을 구함.
			for (int j = 0; j < _Lorig.rows; j++) {
				org = _Pvec.at[i][0];
				_Lorig.at[org][j] = _Ltemp.at[i][j];
				if (i == j) // L = L + I 
					_Lorig.at[org][j] = 1;
			}

	//	printMat(_Lorig, "L_orig");
	}
	inipivot(_Ptemp); // I 행렬 만들어서 L = L + I 가 되도록 함. 
	_Ltemp = Matsum(_Ptemp, _Ltemp);
	copyVal(_Ltemp, _L); // 아웃풋 L U P에 값을 그대로 저장.
	copyVal(_Atemp, _U);
	copyVal(_Pcopy, _P);

	
	return ;
}

// scaled partial pivoting
void SPpivot(Matrix _A, Matrix _P, int _k) { // scaled partial pivoting 함수, 각 열은 k값에 따라 정해지고, 행마다의 max값을 찾아서 S.P 값을 찾아내는 함수
	int max = 0;							 // S.P 값이 만약 1이거나 동일할 경우 각 행의 A값중 큰값으로 결정.
	double maxj = 0;
	int maxrow = 0;
	int _kex = 0;
	double temp = 0;
	Matrix vecMax = createMat(_A.rows - _k, 1);

	for (int i = _k; i < _A.rows; i++) // 각 행의 MAX 값의 절대값을 찾아 vecMax에 sp값을 저장
	{
		max = _k;
		for (int j = _k + 1; j < _A.cols; j++)
		{
			if (fabs(_A.at[i][max]) < fabs(_A.at[i][j]))
				max = j;
		}
		maxj = fabs(_A.at[i][max]);
		vecMax.at[i - _k][0] = fabs(_A.at[i][_k]) / maxj;
	}

	for (int r = 1; r < vecMax.rows; r++) //가장 큰 SP값의 행을 찾아 k번째 행과 바꿈. 만약 k번째 행이 가장 클경우엔 피벗값은 그대로된다.
	{
		if (vecMax.at[maxrow][0] < vecMax.at[r][0])
			maxrow = r;
		else if (vecMax.at[maxrow][0] == vecMax.at[r][0])
		{
			for (int t = _k; t < _A.cols; t++)
			{
				if (_A.at[maxrow + _k][t] < _A.at[r + _k][t])
				{
					maxrow = r;
					break;
				}
				else if (_A.at[maxrow + _k][t] > _A.at[r + _k][t])
					break;
			}
		}
		else
			maxrow = maxrow;
	}

	if (_k != maxrow + _k) { // k번째 행과 가장 큰 행값을 바꿈.
		for (int w = 0; w < _P.cols; w++) {
			temp = _P.at[_k][w];
			_P.at[_k][w] = _P.at[maxrow + _k][w];
			_P.at[maxrow + _k][w] = temp;

		}
	}

	return;
}





//Gauss Elimination
void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d ) {
	if (_A.rows != _A.cols) // 과제의 조건인 A가 square가 아닐때 오류로 실행이 안되게 함.
	{ 
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return ;
	}
	if (_A.rows != _b.rows) // A와 벡터b의 행이 맞지 않을때 오류를 출력하고 실행이 안되게함.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'gaussElim' function");
		printf("\n*************************************************\n");
		return ;
	}

	Matrix Ab = createMat(_A.rows, _A.cols+_b.cols);
	initMat(Ab, 0);
		Ab = pasteMat(_A, _b); // paste 함수로 A와 b를 붙힌 A|b 를 만듬.
	
	double rankA = rankMat(_A);   // 오류를 찾기 위해 rank 함수를 실행하여 A와 A|B행렬의 rank를 구해 조건이 맞을경우에만 실행하도록함
	double rankAb = rankMat(Ab);
	double n = _A.rows; // 실행 조건을 위한 n equation 변수
	if (rankA < rankAb) // A의 랭크가 Ab보다 낮을경우 해를 구할수없음을 오류로 출력
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f < rank(A|b) %.0f : No solution",rankA,rankAb);
		printf("\n*************************************************\n");
		return ;
	}
	else if ((rankA == rankAb) && (rankA<n)) // A의 랭크가 Ab와 같고 equation n보다 낮을경우 해를 구할수없음을 오류로 출력
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f = rank(A|b) %.0f < n %0.f : infinite solution", rankA, rankAb, n);
		printf("\n*************************************************\n");
		return ;
	}

	printf("\n*************************************************"); // 정상 실행이 될경우 A와 B의 랭크를 보여주고 바로 실행하도록 함.
	printf("\n  Gauss Elimination : rank(A) %.0f = rank(A|b) %.0f ", rankA, rankAb);
	printf("\n*************************************************\n");

	for (int k = 0; k < _A.rows-1; k++) // gauss elimination으로 U|d 행렬을 구함. 
	{
		if(Ab.at[k][k]==0){ // a.11 , a.22와 같이 a.kk 의 값이 0일경우 0으로 나눠져 무한이 됨으로 오류 출력을 실행을 하지 않도록 함.
			printf("\n*************************************************");
			printf("\n ERROR!!: Division by zero at '%d' roof gauss Elimination ",k);
			printf("\n*************************************************\n");
			return ;
		}
		for (int i = k + 1; i < _A.rows; i++) { // 첫 열을 계산후 0으로 변함으로 temp에 a.ik번째 값을 저장해 계산하도록 함
			double temp = Ab.at[i][k];
			for (int j = k; j < Ab.cols; j++) {
				Ab.at[i][j] = Ab.at[i][j] - (temp / Ab.at[k][k]) * Ab.at[k][j];
			}
		}
	}
	for (int i = 0; i < _U.rows; i++) // A|b(U|d)의 A열 까지 U에 복사하도록함.
		for (int j = 0; j < _U.cols; j++)
			_U.at[i][j] = Ab.at[i][j];

	for (int i = 0; i < _d.rows; i++) // A|b(U|d)의 A열끝 값 이후부터는 d에 복사 하도록 함. 
		for (int j = 0; j < _d.cols; j++) 
			_d.at[i][j] = Ab.at[i][_A.cols+j];
	return ;
}
// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

void	backSub(Matrix _A, Matrix _b, Matrix _x) // 백 서브티션 수업때 배운 것을 토대로 함수를 void로 한다음 x에 출력값을 넣도록함
{
	if ((_A.rows != _b.rows) || (_b.rows != _x.rows) ) { // a와 b의 행이 다르거나 x의 행이 다를 경우 오류 출력후 실행이 안되도록함.
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'backSub' function");
		printf("\n*************************************************\n");
		return ;
	}
	for (int i = _A.rows-1; i >=0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * _x.at[j][0];
		if (_A.at[i][i] == 0) // _A.at[i][i] 가 0일경우 0을 나누게 됨으로 오류 출력후 실행이 안되도록 함. 
		{
			printf("\n*************************************************");
			printf("\n  ERROR!!: _A[%d][%d] is zero division at BackSub",i,i);
			printf("\n*************************************************\n");
			return;
		}

		_x.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
	return ;
}

void	fwdSub(Matrix _A, Matrix _b, Matrix _x) {
	if ((_A.rows != _b.rows) || (_b.rows != _x.rows)) { // a와 b의 행이 다르거나 x의 행이 다를 경우 오류 출력후 실행이 안되도록함.
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'ForwardSub' function");
		printf("\n*************************************************\n");
		return;
	}
	for (int i = 0; i <= _A.rows - 1; i++) 
	{
		double temp = 0;
		for (int j = 0; j < i+1; j++)
			temp += _A.at[i][j] * _x.at[j][0];
		if (_A.at[i][i] == 0) // _A.at[i][i] 가 0일경우 0을 나누게 됨으로 오류 출력후 실행이 안되도록 함. 
		{
			printf("\n*************************************************");
			printf("\n  ERROR!!:  A[%d][%d] is zero division at ForwardSub ", i, i);
			printf("\n*************************************************\n");
			return;
		}

		_x.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
}


// Q1

double Ffunc(double _x) // f(x) 수식 값 입력
{
	double F = 2.624 * (10 * 10 * 10 * 10) * exp(-0.03734 * _x) - 800;
	return F;
}

double Fdfunc(double _x) // // f'(x) 수식 값 입력
{
	double dF = 2.624 * (10 * 10 * 10 * 10) * (-0.03734) * exp(-0.03734 * _x);
	return dF;
}


double newtonRhapson(double x0, double tol) // newtonRaphson method 함수 정의
{
	int i = 0; // iteration 값 정의
	double xn = 0; // xn 변수 정의 정확히는 xn+1
	double x = x0; // 초기 x 값을 넣는 변수, xn 
	int Nmax = 1000; // 반복 횟수는 1000번까지 iIteration 횟수 
	double ep = 0; // tolerance 비교 값을 위한 변수
	if (Fdfunc(x) != 0) { // f'(x)가 0이 되면 xn이 발산하기 때문에 0이 아닐때만 실행하도록 함
		for (i; i < Nmax; i++) {
			xn = -func(x) / dfunc(x) + x;
			ep = fabs((xn - x) / x);
			printf("K:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			if (ep <= tol) // tolerance 이하가 되었을때 반복문을 정지
				break;
			else if (fabs(xn) > 10000 * fabs(x)) {  // 해를 찾지 못하고 급히 발산할 경우 발산 에러 출력
				printf("Error : Diverging error\n ");
				break;
			}
			x = xn;
		}
		if (i == Nmax && ep > tol)
			printf("Error : infinite loop error\n "); // 해를 찾지 못하고 MaxItr에 도달 했을 경우 무한루프 에러 출력
	}
	else  // f'(x) = 0일 경우 해를 찾을수 없기에 오류 출력 후 정지 
		printf("Error : f'(x) is 0 , you can't find solution by using Newton-Raphson Method\n");

	return xn;
}


double bisectionNL(double _a0, double _b0, double _tol) // bisection method 함수 정의 튜토리얼에 올라온 것 그대로 가져왔지만, 해를 못찾는 경우를 추가했습니다.
{
	int k = 0;
	int Nmax = 1000;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;
	double tst = 0;
	if (Ffunc(a) * Ffunc(b) < 0) {// f(a)*f(b) < 0 일 경우에만 바이섹션이 성립함으로 조건이 만족할때만 사용되게함.
		do {
			xn = (a + b) / 2;
			ep = fabs(Ffunc(xn));
			printf("Iteration:%d \t", k);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			if (Ffunc(a) * Ffunc(xn) < 0)
				b = xn;
			else
				a = xn;
			k++;
		} while (k<Nmax && ep>_tol);
		if (k<Nmax && ep>_tol) // a b 안에서 해를 찾지 못하고 MaxItr에 도달 했을 경우 무한루프 에러 출력
			printf("Error : infinite loop errors\n ");
	}
	else
		printf("Error : f(a)*f(b) >= 0 , you can't find solution by using Bisection Method\n");

	return xn;
}



double NewtonBisection(double a0, double b0, double x0, double tol) {
	double a = a0;
	double b = b0;
	double aa = a;
	double bb = b;
	double xT = x0;
	double xnn = 0;
	int i = 0; // iteration 값 정의
	double xn = 0; // xn 변수 정의 정확히는 xn+1
	int Nmax = 1000; // 반복 횟수는 1000번까지 iIteration 횟수 
	double ep = 0; // tolerance 비교 값을 위한 변수
	for (i; i < Nmax; i++) {
		if ((Fdfunc(xT) != 0 && (xT <= b && xT >= a))) { // f'(x)가 0이 되면 xn이 발산하기 때문에 0이 아닐때만 실행하도록 함
			xn = -Ffunc(xT) / Fdfunc(xT) + xT;        // xT가 a와 b 범위 안에 있을 경우 실행
			ep = fabs((xn - xT) / xT);
			printf("Iteration:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);
			xnn = xT;
			xT = xn;

			if (ep <= tol) // tolerance 이하가 되었을때 반복문을 정지
				break;
		}
		else if (Fdfunc(xT) == 0 || (xT > b || xT < a)) { // f'(x) = 0일 경우와, xT가 a 와 b 범위 밖에 있을때 바이섹션 실행
			if (xT > b) {
				bb = xT;
				aa = a;
			}
			else if (xT < a) {
				bb = b;
				aa = xT;
			}
			xn = (aa + bb) / 2;
			ep = fabs(Ffunc(xn));
			printf("Iteration:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f Used bisection\n", ep);
			xnn = xT;
			xT = xn;
			if (ep <= tol) // tolerance 이하가 되었을때 반복문을 정지
				break;
		}
		else if ((xT < b && xT > a) && (xn == xnn || (fabs(xn) > 10 * fabs(xnn) || (fabs(xn) > fabs(xnn) / 10) && i != 0))) { // xn == xnn이 되어 발산이 될 경우와, xn값이 급격히 커져 발산될경우 바이섹션으로 실행
			if (xT > b) {
				bb = xT;
				aa = a;
			}
			else if (xT < a) {
				bb = b;
				aa = xT;
			}
			xn = (aa + bb) / 2;
			ep = fabs(Ffunc(xn));
			printf("Iteration:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f Used bisection\n", ep);
			xnn = xT;
			xT = xn;
			if (ep <= tol) // tolerance 이하가 되었을때 반복문을 정지
				break;
		}
	}
	return xn;
}