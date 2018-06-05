#ifndef SVD_H_
#define SVD_H_

#include<cmath>
#include "matrix.h"
#include "Math.h"
using namespace std;

/**
 * We should turn it into vector<double>, 1D array to save memory and speed
 */
class svd{
public:
	svd();
//	svd(double **A, int m, int n);
	svd(vector<double> & A, const int mm, const int nn);
//	svd(double **u, double **v, double *ss, int mm, int nn): U(u), V(v), s(ss), m(mm), n(nn){};
	virtual ~svd();

	matrix getMatrixU () const { return matrix(U, m, iMin( m+1, n) );};
	matrix getMatrixV () const { return matrix(V, n, n);};
	matrix getMatrixS () const;
	vector<double>  getU ()const { return U;};
	vector<double>  getV ()const { return V;};
	vector<double>  getSingularValues() const {return s;};

private:
	vector<double> U;
	vector<double> V;
	vector<double> s;
//	double **U;
//	double **V;
//	double  *s;
	int m;
	int n;
};

/**  utility routines **/
    /** Check magnitude of difference of scalars. **/
static void check(double x, double y) {
	double eps = pow(2.0, -50.0);
	if (x == 0 & fabs(y) < eps) return;
	if (y == 0 & fabs(x) < eps) return;
	if (fabs(x-y) > eps * dMax(fabs(x),fabs(y))) {
		cout<<"The difference x-y is too large: x = "
				<< x << "  y = " << y<<endl;
		exit(1);
	}
}

    /** Check norm of difference of "vectors". **/
static void check(double *x, double *y, int m, int n) {
    	if (m == n ) {
    		for (int i=0; i < m; i++) {
    			check(x[i], y[i]);
    		}
    	} else {
    		cout<<"Attempt to compare vectors of different lengths\n";
    		exit(1);
    	}
    }

//    /** Check norm of difference of arrays. **/
//void check( double** x, double **y, int m, int n ) {
//	matrix *A = new matrix(x, m, n);
//	matrix *B = new matrix(y, m, n);
//	check(A,B);
//}
//
///** Check norm of difference of Matrices. **/
//void check(matrix *X, matrix *Y) {
//	double eps = pow(2.0, -50.0);
//	if (X->norm1() == 0. & Y->norm1() < eps) return;
//	if (X->minus(Y)->norm1() > 100.0 * eps * dMax(X->norm1(),Y->norm1())) {
//		cout<<"The norm of (X-Y) is too large: \n";
//		exit(1);
//	}
//}

#endif /*SVD_H_*/
