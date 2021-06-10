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
			slope1 = func(xt, y[i]); // ���� ���� x[i]�� y[i]�� slope1 ���� ���ϸ鼭 yE[i+1]�� ���Ѵ�.
			yE = y[i] + func(xt, y[i]) * h;
			xt = xt + h; // i+1�� 
			slope2 = func(xt, yE); // ���� ���Ϸ��� ���� yE[i+1] �� ���� slope 2��
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


double func(const double x) // f(x) ���� �� �Է�
{
	return x*x*x;
}

double dfunc(const double x) // // f'(x) ���� �� �Է�
{
	return 3*x*x;
}

double newtonRaphsonfunc(double func(const double x), double dfunc(const double x), double x0, double tol) // newtonRaphson method �Լ� ����
{
	int i = 0; // iteration �� ����
	double xn = 0; // xn ���� ���� ��Ȯ���� xn+1
	double x = x0; // �ʱ� x ���� �ִ� ����, xn 
	int Nmax = 1000; // �ݺ� Ƚ���� 1000������ iIteration Ƚ�� 
	double ep = 0; // tolerance �� ���� ���� ����
	if (dfunc(x) != 0) { // f'(x)�� 0�� �Ǹ� xn�� �߻��ϱ� ������ 0�� �ƴҶ��� �����ϵ��� ��
		for (i; i < Nmax; i++) {
			xn = -func(x) / dfunc(x) + x;
			ep = fabs((xn - x) / x);
			printf("K:%d \t", i);
			printf("X(n): %0.10f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			if (ep <= tol) // tolerance ���ϰ� �Ǿ����� �ݺ����� ����
				break;
			else if (fabs(xn) > 10000 * fabs(x)) {  // �ظ� ã�� ���ϰ� ���� �߻��� ��� �߻� ���� ���
				printf("Error : Diverging error\n ");
				break;
			}
			x = xn;
		}
		if (i == Nmax && ep > tol)
			printf("Error : infinite loop error\n "); // �ظ� ã�� ���ϰ� MaxItr�� ���� ���� ��� ���ѷ��� ���� ���
	}
	else  // f'(x) = 0�� ��� �ظ� ã���� ���⿡ ���� ��� �� ���� 
		printf("Error : f'(x) is 0 , you can't find solution by using Newton-Raphson Method\n");

	return xn;
}



void  gradient(double x[], double y[], double dydx[], int m) {
	double h = x[1] - x[0];
	if ( m < 2) { //  2���� �۾� �Լ��� �� �����ų�� �������� ���� ���
		printf("Error!! : �Է°��� ���� 2���� ����");
		return ;
	}
	else if (m == 2) {// ���̰� 2�� ��쿡  2����Ʈ ������ �� ������ ���
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
	if ((mx != my) || (mx < 2)) { // x�� y�� ���̰� �ȸ°ų� 2���� �۾� �Լ��� �� �����ų�� �������� ���� ���
		printf("Error!! : �Է°��� ���̰� ���� ����");
		return createMat(0, 0);
	}
	else if ((mx == my) && (mx == 2)) { // ���̰� 2�� ��쿡  2����Ʈ ������ �� ������ ���
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


double FindValueLIP(Matrix yq, Matrix Tq, double T) {   // Ư�� T�� ���� ���� �� yq�� ã�� �Լ�
	double Fx = 0;
	int a = 0;
	int no = 0;
	for (int i = 0; i < Tq.rows; i++) { // T�� Tq=xq ��ġ�� ã�� ���� ��� no���� �������� ������ ����ϵ�����.
		if (Tq.at[i][0] == T)
			a = i;
		else
			no++;
	}
	if (no == Tq.rows) {
		printf("\n*************************************************");  // ���������� ����  T�� Tquery ��� �ȿ� ���� ��� ���� ���
		printf("\n ERROR!!: ������ų ���� yq �ȿ� �������� ���� ");
		printf("\n*************************************************\n");
		return Fx = 0;
	}

	Fx = yq.at[a][0];

	return Fx;
}


Matrix linearInterp(Matrix x, Matrix y, Matrix xq) { // ���Ͼ� ���������̼� �Լ� 
	Matrix yq = createMat(xq.rows, 1);
	int mx = x.rows;
	int my = y.rows;
	if ((mx != my) || (mx < 2)) { // x�� y�� ���̰� �ȸ°ų� 2���� �۾� �Լ��� �� �����ų�� �������� ���� ���
		printf("Error!! : �Է°��� ���̰� ���� �ʰų� 2���� �۾� Linear Interpolation �� ���� ��ų �� ����.");
		return zeros(yq.rows,yq.cols); 
	}
	int i = 0;
	for (int j = 0; j < xq.rows; j++) {


		if ((xq.at[j][0] >= x.at[i + 1][0]) || (xq.at[j][0] < x.at[i][0])) { // xq�� x���� ��� ��ġ�ϴ��� ã�� �׿� �ش��ϴ� yq�� ������ ã�� ����.
			for (int k = i; k < x.rows-1; k++) {
				if ((xq.at[j][0] < x.at[k + 1][0]) && (xq.at[j][0] >= x.at[k][0])) {
					i = k;
					break;
				}
			}
		}

		// ��ġ�� ã�� ���Ͼ� �Լ��� ���� xq���� �־� yq���� ����
		yq.at[j][0] = y.at[i][0] * (xq.at[j][0] - x.at[i + 1][0]) / (x.at[i][0] - x.at[i + 1][0]) + y.at[i + 1][0] * (xq.at[j][0] - x.at[i][0]) / (x.at[i + 1][0] - x.at[i][0]); 

	}

	return yq;

}

double FindValueLCF(Matrix z, double T) { // Ŀ�� ���ÿ��� ���ϴ� �Է°� T�� �������� ã�� �Լ�
	double Fx= 0;
	Fx = z.at[0][0] + z.at[1][0] * T;

	return Fx;
}

Matrix	linearFit(Matrix _x, Matrix _y) { // ���Ͼ� Ŀ�� ���� �Լ�
	int mx = _x.rows;
	int my = _y.rows;
	double a1 = 0;
	double a0 = 0;
	if ((mx != my) || (mx < 2)) { //���� ���
		printf("Error!! : �Է°��� ���̰� ���� �ʰų� 2���� �۾� Linear curve fitting �� ���� ��ų �� ����.");
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

double Cond(Matrix A) { //����� �Լ� A�� ����� �ѹ��� ã�� ���� �Լ�


	Matrix At = createMat(A.cols, A.rows);     
	Matrix AtA = createMat(A.cols, A.cols);
	Matrix eig = createMat(AtA.rows, 1);
	Matrix Atemp = createMat(A.cols, A.cols);
	At = transpose(A);

	Atemp = MatMult(At, A);
	copyVal(Atemp, AtA);
	QReig(AtA,eig);   // AT*A�� ����� �ѹ��� QR factorization ���� ���̰� ���͸� ã�´�. 
	//printMat(eig, "eig");
	double normA = norm(eig, 3);  // ã�� ���̰� ���Ϳ��� �ִ밪�� AT*A �� �����Ѵ�.
	double normAinv = 1/ norm(eig, 4); //  AT*A�� ���Լ��� �ִ밪�� ���̰� ������ 0�� �ƴ�  �ּҰ��� 1�� ���� ���̴�. 
	double ans = sqrt(normA) * sqrt(normAinv);
	
	return ans;
}

void QReig(Matrix A, Matrix Eig) { // QR decopmsition �� ���� QR�� ���ϰ� A���� similat matrix�� �� �� ���� �ݺ��Ͽ� ���̰� ���͸� ã�´�.
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
	for (int k = 0; k < nmax; k++) { // A�� similat matrix�� �� �� ���� �ݺ��Ͽ� QRdecomp�� �����Ѵ�. 
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

	// ���ǽ� ������ , ��ũ a = n �̿�����.
	if (A.rows != A.cols) //  A�� square�� �ƴҶ� ������ ������ �ȵǰ� ��.
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

		v = Matsum(c, e);  // ������� v = c + ||c||*e  ���ϴ� �˰���
		normv = norm(v, 2); // vtv�� norm2�� �����̹Ƿ� v������ norm2�� ����.
		vt = transpose(v);
		vvt = MatMult(v, vt);
		if (normv == 0) {
			printf("\n*************************************************"); // 0���� ���������� ��� �Լ� �ߴ�.
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
double norm(Matrix vec, int k) { // norm �Լ� 1�϶��� ���밪�� ��� ���Ѱ�, 2�϶��� ũ�⸦ ���ϴ°�, 3�� �ִ밪 , 4�� 0�� �ƴ� �ּҰ�.
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
//		FX.at[i][0] = 4 * X.at[i][0]; // ���Ĵ���
//	}
//	// F1 = ���� x,y ����
//
//	return ;
//}





// inverse matrix
double inv(Matrix _A, Matrix _Ainv) { // ����� �Լ�, LU�� �̿��߰�, Ainv�� �� ������ LUsolve�Ͽ� ���� ���س�.

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
	// ���ǽ� ������ , ��ũ a = n �̿�����.
	if (_A.rows != _A.cols) //  A�� square�� �ƴҶ� ������ ������ �ȵǰ� ��.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return 0;
	}
	double rankA = rankMat(_Atemp);
	double n = _A.rows;
	if ((rankA != n)) // A�� ��ũ��  equation n���� ������� �ظ� ���Ҽ������� ������ ���
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f < n %0.f : infinite solution", rankA, n);
		printf("\n*************************************************\n");
		return 0;
	}
	printf("\n*************************************************");
	printf("\n            Inverse Matrix A  ");
	LUdecomp(_Atemp, L, U, P);
	for (int k = 0; k < _A.cols; k++) { // Ainv�� �� �����ٸ� LU solve�Լ��� �Ἥ ������.
		for (int i = 0; i < I.rows; i++)
			b.at[i][0] = I.at[i][k];
		solveLU(L, U, P, b, X);
		for (int i = 0; i < I.rows; i++)
			_Ainv.at[i][k] = X.at[i][0]; // solve�� ���� X���� �ش��ϴ� ���� Ainv �ӽð��� ������.
	}
	return 0;
}

// LU solve function
void	solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix vectorX) {
	if ( ((_L.rows != _U.rows) || (_L.cols != _U.cols)) || ((_U.rows != _P.rows) || (_U.cols != _P.cols))|| (_P.rows != _b.rows) || ((_b.rows != vectorX.rows) || (_b.cols != vectorX.cols)))
	{ // ��ǲ�� �ƿ�ǲ�� ����� �߸��Ǿ����� ���� ���
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at L & U & P & b & x LU solve function");
		printf("\n*************************************************\n");
		return;
	}
	Matrix y = createMat(_U.rows, vectorX.cols); // y�� Ux ��
	Matrix _btemp = createMat(_b.rows, _b.cols);
	initMat(y, 0);
	_btemp = MatMult(_P, _b); // PA = LU = Pb �̹Ƿ� �ǹ��� b���� ������. 
	fwdSub(_L, _btemp, y); //  L(Ux=y)= Pb�� ���� ������� ����.
	backSub(_U, y, vectorX); // Ux = y�̹Ƿ� �鼭��� ����. 

	return ;
}



// LUdecomposition
void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P) {
	if (_A.rows != _A.cols) // ������ ������ A�� square�� �ƴҶ� ������ ������ �ȵǰ� ��.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return;
	}
	if( ((_A.rows != _L.rows) || (_A.cols != _L.cols)) || ((_L.rows != _U.rows) || (_L.cols != _U.cols)) || ((_U.rows != _P.rows) || (_U.cols != _P.cols)) ) 
	{ // ��ǲ�� �ƿ�ǲ�� ����� �߸��Ǿ����� ���� ���
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A & L & U & P  LUdecomposition ");
		printf("\n*************************************************\n");
		return;
	}

	double rankA = rankMat(_A);
	double n = _A.rows;
	if ((rankA != n)) // A�� ��ũ��  equation n���� ������� �ظ� ���Ҽ������� ������ ���
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
		SPpivot(_Atemp,_Ptemp, k); //Scaled pivot���� �ǹ��� P�� ����, 
		_Ltemp = MatMult(_Ptemp, _Ltemp); // L P A �� Scaled pivot �� ���Ͽ� �ǹ��� ���� �����Ѵ�.
		_Pcopy = MatMult(_Ptemp, _Pcopy);
	//	printMat(_Pcopy, "P");
		_Atemp = MatMult(_Ptemp, _Atemp);
		_Lsub = copyMat(_Ltemp);
		for (int i = k; i < _A.rows - 1; i++) {
			_Ltemp.at[i + 1][k] = _Atemp.at[i + 1][k] / _Atemp.at[k][k];
			if (_Atemp.at[k][k]==0) 
			{
				printf("\n*************************************************"); // 0���� ���������� ��� �Լ� �ߴ�.
				printf("\n ERROR!!: Division by zero at '%d' LU decomposition ", k);
				printf("\n*************************************************\n");
				return;
			}
		}

		_Lsub = Matsub(_Ltemp, _Lsub); // ���ߵǴ� ���� L������ ����
		_LAtemp = MatMult(_Lsub, _Atemp); // LA�� ����
		_Atemp = Matsub(_Atemp, _LAtemp); // ���� A���� LA���� ���� U����� ����
		//printMat(_Atemp, "U");

		for (int i=0;i<_Pcopy.rows;i++) // �ǹ� �Ǳ����� L�� ���ϱ� ���� �ǹ��� ��� �ٲ���� ���� ������ ���ڷ� ���Ϳ� ����.
			for (int j = 0; j < _Pcopy.cols; j++) {
				if (_Pcopy.at[i][j] == 1)
					_Pvec.at[i][0] = j;
			}
		
		for (int i = 0; i < _Lorig.rows; i++) // ���� ���͸� �̿��Ͽ� �� ���� �˸´� ������ �ٲپ L original ���� ����.
			for (int j = 0; j < _Lorig.rows; j++) {
				org = _Pvec.at[i][0];
				_Lorig.at[org][j] = _Ltemp.at[i][j];
				if (i == j) // L = L + I 
					_Lorig.at[org][j] = 1;
			}

	//	printMat(_Lorig, "L_orig");
	}
	inipivot(_Ptemp); // I ��� ���� L = L + I �� �ǵ��� ��. 
	_Ltemp = Matsum(_Ptemp, _Ltemp);
	copyVal(_Ltemp, _L); // �ƿ�ǲ L U P�� ���� �״�� ����.
	copyVal(_Atemp, _U);
	copyVal(_Pcopy, _P);

	
	return ;
}

// scaled partial pivoting
void SPpivot(Matrix _A, Matrix _P, int _k) { // scaled partial pivoting �Լ�, �� ���� k���� ���� ��������, �ึ���� max���� ã�Ƽ� S.P ���� ã�Ƴ��� �Լ�
	int max = 0;							 // S.P ���� ���� 1�̰ų� ������ ��� �� ���� A���� ū������ ����.
	double maxj = 0;
	int maxrow = 0;
	int _kex = 0;
	double temp = 0;
	Matrix vecMax = createMat(_A.rows - _k, 1);

	for (int i = _k; i < _A.rows; i++) // �� ���� MAX ���� ���밪�� ã�� vecMax�� sp���� ����
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

	for (int r = 1; r < vecMax.rows; r++) //���� ū SP���� ���� ã�� k��° ��� �ٲ�. ���� k��° ���� ���� Ŭ��쿣 �ǹ����� �״�εȴ�.
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

	if (_k != maxrow + _k) { // k��° ��� ���� ū �ప�� �ٲ�.
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
	if (_A.rows != _A.cols) // ������ ������ A�� square�� �ƴҶ� ������ ������ �ȵǰ� ��.
	{ 
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at A matrix isn't square");
		printf("\n*************************************************\n");
		return ;
	}
	if (_A.rows != _b.rows) // A�� ����b�� ���� ���� ������ ������ ����ϰ� ������ �ȵǰ���.
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'gaussElim' function");
		printf("\n*************************************************\n");
		return ;
	}

	Matrix Ab = createMat(_A.rows, _A.cols+_b.cols);
	initMat(Ab, 0);
		Ab = pasteMat(_A, _b); // paste �Լ��� A�� b�� ���� A|b �� ����.
	
	double rankA = rankMat(_A);   // ������ ã�� ���� rank �Լ��� �����Ͽ� A�� A|B����� rank�� ���� ������ ������쿡�� �����ϵ�����
	double rankAb = rankMat(Ab);
	double n = _A.rows; // ���� ������ ���� n equation ����
	if (rankA < rankAb) // A�� ��ũ�� Ab���� ������� �ظ� ���Ҽ������� ������ ���
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f < rank(A|b) %.0f : No solution",rankA,rankAb);
		printf("\n*************************************************\n");
		return ;
	}
	else if ((rankA == rankAb) && (rankA<n)) // A�� ��ũ�� Ab�� ���� equation n���� ������� �ظ� ���Ҽ������� ������ ���
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: rank(A) %.0f = rank(A|b) %.0f < n %0.f : infinite solution", rankA, rankAb, n);
		printf("\n*************************************************\n");
		return ;
	}

	printf("\n*************************************************"); // ���� ������ �ɰ�� A�� B�� ��ũ�� �����ְ� �ٷ� �����ϵ��� ��.
	printf("\n  Gauss Elimination : rank(A) %.0f = rank(A|b) %.0f ", rankA, rankAb);
	printf("\n*************************************************\n");

	for (int k = 0; k < _A.rows-1; k++) // gauss elimination���� U|d ����� ����. 
	{
		if(Ab.at[k][k]==0){ // a.11 , a.22�� ���� a.kk �� ���� 0�ϰ�� 0���� ������ ������ ������ ���� ����� ������ ���� �ʵ��� ��.
			printf("\n*************************************************");
			printf("\n ERROR!!: Division by zero at '%d' roof gauss Elimination ",k);
			printf("\n*************************************************\n");
			return ;
		}
		for (int i = k + 1; i < _A.rows; i++) { // ù ���� ����� 0���� �������� temp�� a.ik��° ���� ������ ����ϵ��� ��
			double temp = Ab.at[i][k];
			for (int j = k; j < Ab.cols; j++) {
				Ab.at[i][j] = Ab.at[i][j] - (temp / Ab.at[k][k]) * Ab.at[k][j];
			}
		}
	}
	for (int i = 0; i < _U.rows; i++) // A|b(U|d)�� A�� ���� U�� �����ϵ�����.
		for (int j = 0; j < _U.cols; j++)
			_U.at[i][j] = Ab.at[i][j];

	for (int i = 0; i < _d.rows; i++) // A|b(U|d)�� A���� �� ���ĺ��ʹ� d�� ���� �ϵ��� ��. 
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

void	backSub(Matrix _A, Matrix _b, Matrix _x) // �� ����Ƽ�� ������ ��� ���� ���� �Լ��� void�� �Ѵ��� x�� ��°��� �ֵ�����
{
	if ((_A.rows != _b.rows) || (_b.rows != _x.rows) ) { // a�� b�� ���� �ٸ��ų� x�� ���� �ٸ� ��� ���� ����� ������ �ȵǵ�����.
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'backSub' function");
		printf("\n*************************************************\n");
		return ;
	}
	for (int i = _A.rows-1; i >=0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * _x.at[j][0];
		if (_A.at[i][i] == 0) // _A.at[i][i] �� 0�ϰ�� 0�� ������ ������ ���� ����� ������ �ȵǵ��� ��. 
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
	if ((_A.rows != _b.rows) || (_b.rows != _x.rows)) { // a�� b�� ���� �ٸ��ų� x�� ���� �ٸ� ��� ���� ����� ������ �ȵǵ�����.
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
		if (_A.at[i][i] == 0) // _A.at[i][i] �� 0�ϰ�� 0�� ������ ������ ���� ����� ������ �ȵǵ��� ��. 
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

double Ffunc(double _x) // f(x) ���� �� �Է�
{
	double F = 2.624 * (10 * 10 * 10 * 10) * exp(-0.03734 * _x) - 800;
	return F;
}

double Fdfunc(double _x) // // f'(x) ���� �� �Է�
{
	double dF = 2.624 * (10 * 10 * 10 * 10) * (-0.03734) * exp(-0.03734 * _x);
	return dF;
}


double newtonRhapson(double x0, double tol) // newtonRaphson method �Լ� ����
{
	int i = 0; // iteration �� ����
	double xn = 0; // xn ���� ���� ��Ȯ���� xn+1
	double x = x0; // �ʱ� x ���� �ִ� ����, xn 
	int Nmax = 1000; // �ݺ� Ƚ���� 1000������ iIteration Ƚ�� 
	double ep = 0; // tolerance �� ���� ���� ����
	if (Fdfunc(x) != 0) { // f'(x)�� 0�� �Ǹ� xn�� �߻��ϱ� ������ 0�� �ƴҶ��� �����ϵ��� ��
		for (i; i < Nmax; i++) {
			xn = -func(x) / dfunc(x) + x;
			ep = fabs((xn - x) / x);
			printf("K:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			if (ep <= tol) // tolerance ���ϰ� �Ǿ����� �ݺ����� ����
				break;
			else if (fabs(xn) > 10000 * fabs(x)) {  // �ظ� ã�� ���ϰ� ���� �߻��� ��� �߻� ���� ���
				printf("Error : Diverging error\n ");
				break;
			}
			x = xn;
		}
		if (i == Nmax && ep > tol)
			printf("Error : infinite loop error\n "); // �ظ� ã�� ���ϰ� MaxItr�� ���� ���� ��� ���ѷ��� ���� ���
	}
	else  // f'(x) = 0�� ��� �ظ� ã���� ���⿡ ���� ��� �� ���� 
		printf("Error : f'(x) is 0 , you can't find solution by using Newton-Raphson Method\n");

	return xn;
}


double bisectionNL(double _a0, double _b0, double _tol) // bisection method �Լ� ���� Ʃ�丮�� �ö�� �� �״�� ����������, �ظ� ��ã�� ��츦 �߰��߽��ϴ�.
{
	int k = 0;
	int Nmax = 1000;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;
	double tst = 0;
	if (Ffunc(a) * Ffunc(b) < 0) {// f(a)*f(b) < 0 �� ��쿡�� ���̼����� ���������� ������ �����Ҷ��� ���ǰ���.
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
		if (k<Nmax && ep>_tol) // a b �ȿ��� �ظ� ã�� ���ϰ� MaxItr�� ���� ���� ��� ���ѷ��� ���� ���
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
	int i = 0; // iteration �� ����
	double xn = 0; // xn ���� ���� ��Ȯ���� xn+1
	int Nmax = 1000; // �ݺ� Ƚ���� 1000������ iIteration Ƚ�� 
	double ep = 0; // tolerance �� ���� ���� ����
	for (i; i < Nmax; i++) {
		if ((Fdfunc(xT) != 0 && (xT <= b && xT >= a))) { // f'(x)�� 0�� �Ǹ� xn�� �߻��ϱ� ������ 0�� �ƴҶ��� �����ϵ��� ��
			xn = -Ffunc(xT) / Fdfunc(xT) + xT;        // xT�� a�� b ���� �ȿ� ���� ��� ����
			ep = fabs((xn - xT) / xT);
			printf("Iteration:%d \t", i);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);
			xnn = xT;
			xT = xn;

			if (ep <= tol) // tolerance ���ϰ� �Ǿ����� �ݺ����� ����
				break;
		}
		else if (Fdfunc(xT) == 0 || (xT > b || xT < a)) { // f'(x) = 0�� ����, xT�� a �� b ���� �ۿ� ������ ���̼��� ����
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
			if (ep <= tol) // tolerance ���ϰ� �Ǿ����� �ݺ����� ����
				break;
		}
		else if ((xT < b && xT > a) && (xn == xnn || (fabs(xn) > 10 * fabs(xnn) || (fabs(xn) > fabs(xnn) / 10) && i != 0))) { // xn == xnn�� �Ǿ� �߻��� �� ����, xn���� �ް��� Ŀ�� �߻�ɰ�� ���̼������� ����
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
			if (ep <= tol) // tolerance ���ϰ� �Ǿ����� �ݺ����� ����
				break;
		}
	}
	return xn;
}