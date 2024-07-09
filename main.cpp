#include "stdafx.h"
#include <iostream>
using namespace std;

double *bm, **x, *Xk, **beta, **xi, **w, **matr, *p1, *p2, *p3, *pp1, *pp2, *pp3;
int n;
double uravn(double*, double, int);
double kasat(double(*f)(double*, double, int), double(*fpr1)(double), double(*fpr2)(double), double eps, double a, double b, double*cn, int n);
void gelg(double**, double*, int, double, int&);
void ymn(double**, int);
void BETA(double*, double*);
void XI(double**, double**);
void check(double*, double**, double**);

void gelg(double **x, double *bm, int n, double eps, int &kr)
{
	int i, j, k, imax, jmax, rint;
	double max, rdouble, xkk;
	long double dxkk;
	int *cod = new int[n];
	for (i = 0; i<n; i++) cod[i] = i;
	kr = 0;

	//Цикл нахождения главных элементов
	for (k = 0; k <= n - 2; k++)
	{

		//максимальный элемент матрицы
		max = fabs(x[k][k]); imax = k; jmax = k;	//Берем диагональный элемент
		for (i = k; i<n; i++)
			for (j = k; j < n; j++)
			{
				rdouble = fabs(x[i][j]);			//Берем поочередно каждый элемент
				if (rdouble > max)					//сравниваем диагональный элемент с текущим
				{
					max = rdouble; imax = i; jmax = j;	//если текущий больше диагонального, меняем их местами
				}
			}

		//перестановка строк и правых частей
		for (int p = 0; p < n; p++)
		{
			rdouble = x[imax][p];  //поочередно взятый элемент меняем местами с максимальным по строкам
			x[imax][p] = x[k][p];
			x[k][p] = rdouble;
		}
		rdouble = bm[imax];
		bm[imax] = bm[k];
		bm[k] = rdouble;
		for (int p = 0; p < n; p++)	//по столбцам
		{
			rdouble = x[p][k];
			x[p][k] = x[p][jmax];
			x[p][jmax] = rdouble;
		}
		rint = cod[k];
		cod[k] = cod[jmax];
		cod[jmax] = rint;

		// Проверка главного элемента на <eps;
		xkk = x[k][k];
		if (fabs(xkk) < eps)
		{
			kr = k + 1;
			delete[]cod;
			return;
		}

		//Деление текущей строки и правых частей на  найденный главный эелемент dxkk;
		dxkk = xkk;
		for (j = k; j < n; j++) x[k][j] = x[k][j] / dxkk;
		bm[k] = bm[k] / dxkk;
		for (i = k + 1; i <= n - 1; i++)
		{
			dxkk = x[i][k];
			for (j = k; j<n; j++) x[i][j] = x[i][j] - dxkk*x[k][j];
			bm[i] = bm[i] - dxkk*bm[k];
		}
	}

	// Нахождение решения системы

	if (fabs(x[n - 1][n - 1]) < eps*eps)
	{
		kr = n;
		delete[]cod;
		return;
	}
	bm[n - 1] /= x[n - 1][n - 1];
	for (k = n - 2; k >= 0; k--)
	{
		dxkk = 0;
		for (j = k + 1; j<n; j++) dxkk = dxkk + x[k][j] * bm[j];
		bm[k] = bm[k] - dxkk;
	}

	//Привести к нормальному виду

	for (i = 0; i <= n - 2; i++)
	{
		rint = cod[i];
		imax = i;
		for (j = i + 1; j <= n - 1; j++) if (cod[j]<rint) { rint = cod[j]; imax = j; }
		rint = cod[i];
		cod[i] = cod[imax];
		cod[imax] = rint;
		rdouble = bm[imax];
		bm[imax] = bm[i];
		bm[i] = rdouble;
	}
	delete[]cod;
}
void ymn(double **a, int n)
{
	cout << "c0 : " << endl;
	double *c0 = new double[n];
	c0[0] = 1; c0[1] = 0; c0[2] = 1;
	for (int i = 0; i<n; i++) cout << c0[i] << " ";
	cout << endl << endl; 

	cout << "bm[i] : " << endl;

	for (int i = 0; i<n; i++) bm[i] = c0[i];

	for (int i = 0; i<n; i++) cout << bm[i] << " ";
	cout << endl << endl; 

	/**************************************/

	x = new double *[n];
	for (int i = 0; i<n; i++) x[i] = new double[n];
	for (int i = 0; i<n; i++) for (int j = 0; j<n; j++) x[i][j] = 0;
	//----------------------------------------------------------
	for (int k = 0; k<n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j<n; j++)
				x[i][n - k - 1] += a[i][j] * c0[j];
		}
		for (int m = 0; m<n; m++) c0[m] = x[m][n - k - 1];
	}

	/**************************************/

	cout << "intermediate values :" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << x[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	/**************************************/

	for (int j = 0; j<n - 1; j++)
	{
		double l;
		for (int i = 0; i<n; i++)
		{
			l = x[i][j];
			x[i][j] = x[i][j + 1];
			x[i][j + 1] = l;
		}
	}

	double r;
	for (int i = 0; i<n; i++)
	{
		r = x[i][n - 1];
		x[i][n - 1] = bm[i];
		bm[i] = r;
	}
	/**************************************/
	cout << "Matrix : " << endl;;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << x[i][j] << " ";
		}
		cout << endl;
	}
}
double df(double x)
{
	return 3 * x*x - 12 * x + 11;
}
double ddf(double x)
{
	return 6 * x - 12;
}
double uravn(double*bm, double lamb, int n)
{
	return lamb*lamb*lamb - 6 * lamb*lamb + 11 * lamb - 6;
}
double kasat(double(*f)(double*, double, int), double(*fpr1)(double), double(*fpr2)(double), double eps, double a, double b, double*cn, int n)
{
	double x0 = b;
	if (f(bm, a, n)*fpr2(a)>0) x0 = a;
	else if (f(bm, b, n)*fpr2(b)>0) x0 = b;
	double m1 = a;
	double M2 = b;
	double xn1, xn = x0, xn2, xn3;
	xn1 = xn - f(bm, xn, n) / fpr1(xn);
	xn2 = xn1 - f(bm, xn1, n) / fpr1(xn1);
	xn3 = xn2 - f(bm, xn2, n) / fpr1(xn2);
	while (fabs(xn1 - xn) > sqrt(2 * m1*eps / M2)) return xn2;
	while (fabs(xn2 - xn1) > sqrt(2 * m1*eps / M2)) return xn3;
	return xn1;
}
void BETA(double*xk, double *bm)
{
	beta = new double*[n];
	for (int k = 0; k < n; k++) beta[k] = new double[n];
	for (int j = 0; j < n; ++j)
	{
		beta[0][j] = 1;
	}
	for (int j = 1; j < n; ++j)
		for (int i = 0; i < n; ++i)
		{
			beta[j][i] = Xk[i] * beta[j - 1][i] - bm[j - 1];
		}
}
void XI(double**beta, double **x)
{
	xi = new double *[n];
	for (int k = 0; k<n; k++) xi[k] = new double[n];
	for (int t = 1; t <= n; t++) {
		for (int i = 1; i <= n; i++) {
			double s = 0;
			for (int j = 1; j <= n; j++)
			{
				s += beta[j - 1][t - 1] * w[i - 1][j - 1];
			} xi[i - 1][t - 1] = s;
		}
	}
}
void check(double *Xk, double**xi, double **matr)
{
	p1 = new double[n];
	p2 = new double[n];
	p3 = new double[n];
	pp1 = new double[n];
	pp2 = new double[n];
	pp3 = new double[n];


	for (int k = 0; k < n; k++) {
		p1[k] = 0.;
		p2[k] = 0.;
		p3[k] = 0.;
		pp1[k] = 0.;
		pp2[k] = 0.;
		pp3[k] = 0.;
	}



	for (int i = 0; i < n; i++) {
		p1[i] = xi[0][i] * Xk[i];
		p2[i] = xi[1][i] * Xk[i];
		p3[i] = xi[2][i] * Xk[i];
	}

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++) {
			pp1[i] += xi[0][i] * matr[j][i];
			pp2[i] += xi[1][i] * matr[j][i];
			pp3[i] += xi[2][i] * matr[j][i];
		}
	}

}


int main()
{
	int l;
	cout << "Set equation number : " << endl;
	cin >> n;
	bm = new double[n];

	int z;
	cout << "matrix dimension:" << endl;
	//-----------------Matrix formation-------------------//

	cout << "Innocent matrix : " << endl;
	matr = new double*[n];
	for (int i = 0; i < n; i++)
	{
		matr[i] = new double[n];
	}

	matr[0][0] = 1; matr[0][1] = 0; matr[0][2] = -2;
	matr[1][0] = 2; matr[1][1] = 3; matr[1][2] = -1;
	matr[2][0] = 0; matr[2][1] = 0; matr[2][2] = 2;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << matr[i][j] << " ";
		}
		cout << endl;
	}


	ymn(matr, n);

	w = new double*[n];
	for (int i = 0; i < n; i++)
	{
		w[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			{
				w[i][j] = x[i][j];
			}
		}

	cout << "W : " << endl;;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << w[i][j] << " ";

		}
		cout << endl;
	}

	gelg(x, bm, n, 1e-5, l);



	cout << "Solutions : " << endl << endl;
	for (int i = 0; i<n; i++) cout << bm[i] << endl << endl;


	Xk = new double[n];
	Xk[0] = kasat(uravn, df, ddf, 0.01, 2.9, 3.1, bm, n);
	Xk[1] = kasat(uravn, df, ddf, 0.01, 1.9, 2.1, bm, n);
	Xk[2] = kasat(uravn, df, ddf, 0.01, 0.9, 1.1, bm, n);


	cout << "Equation solution: " << endl;
	for (int k = 0; k < n; k++)
	{
		cout << "x=" << Xk[k] << endl;
	}
	cout << endl;

	BETA(Xk, bm);
	cout << "Eigenvalues = " << endl;
	for (int j = 0; j<n; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			cout << beta[j][i] << " ";
		}
		cout << endl;
	}
	cout << endl;


	XI(beta, x);
	cout << "Eigenvector = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << xi[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	check(Xk, xi, matr);

	double**ch1 = new double*[n];
	double**ch2 = new double*[n];
	for (int k = 0; k < n; ++k) {
		ch1[k] = new double[n];
		ch2[k] = new double[n];
	}

	for (int k = 0; k < n; k++) {
		ch1[0][k] = p1[k];
		ch1[1][k] = p2[k];
		ch1[2][k] = p3[k];

		ch2[0][k] = pp1[k];
		ch2[1][k] = pp2[k];
		ch2[2][k] = pp3[k];
	}

	cout << "CHECK 1 = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << ch1[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "CHECK 2 = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << ch2[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	system("pause");
	return 0;
}
