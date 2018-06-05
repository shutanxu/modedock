#include <cmath>
#include <iostream>

#include "svd.h"
//#include "Math.h"
using namespace std;

svd::svd(){}

//svd::svd(double **A, int mm, int nn){
//	// Derived from LINPACK code.
//
//	// Initialize.
//	m =mm;
//	n = nn;
//	cout<<"inside svd \n";
//	for( int j =0; j<m; j++){
//		for (int k = 0; k< n ; k++)
//			cout<<A[j][k]<<" ";
//		cout<<endl;
//	}
//	int nu = iMin ( m, n );
//	s = new double [iMin(m+1, n)];
//	U = new double* [m];
//	for(int i = 0; i < m; i++){
//		U[i] = new double [nu];
//		for(int j = 0; j<nu; j++)
//			U[i][j] = 0.0;
//	}
//	V = new double* [n];
//	for(int i = 0; i < n; i++)
//		V[i] = new double [n];
//	double *e = new double [n];
//	double *work = new double [m];
//	bool wantu = true;
//	bool wantv = true;
//	// Reduce A to bidiagonal form, storing the diagonal elements
//	// in s and the super-diagonal elements in e.
//	int nct = iMin(m-1, n);
//	int nrt = iMax(0, iMin(n-2, m) );
//	cout<<iMin(m+1, n)<<" ns "<<nct<<"  "<<nrt<<endl;
//
//	for (int k = 0; k < iMax(nct, nrt); k++) {
//		if (k < nct) {
//			// Compute the transformation for the k-th column and
//			// place the k-th diagonal in s[k].
//			// Compute 2-norm of k-th column without under/overflow.
//			s[k] = 0;
//			for (int i = k; i < m; i++) {
//				s[k] = hypot(s[k], A[i][k]);
//			}
//			if (s[k] != 0.0) {
//				if (A[k][k] < 0.0) {
//					s[k] = -s[k];
//				}
//				for (int i = k; i < m; i++) {
//					A[i][k] /= s[k];
//				}
//				A[k][k] += 1.0;
//			}
//			s[k] = -s[k];
//		}
//
//		for (int j = k+1; j < n; j++) {
//			if ((k < nct) & (s[k] != 0.0))  {
//				// Apply the transformation.
//				double t = 0;
//				for (int i = k; i < m; i++) {
//					t += A[i][k]*A[i][j];
//				}
//				t = -t/A[k][k];
//				for (int i = k; i < m; i++) {
//					A[i][j] += t*A[i][k];
//				}
//			}
//			// Place the k-th row of A into e for the
//			// subsequent calculation of the row transformation.
//			e[j] = A[k][j];
//		}
//		if (wantu & (k < nct)) {
//			// Place the transformation in U for subsequent back
//			// multiplication.
//			for (int i = k; i < m; i++) {
//				U[i][k] = A[i][k];
//			}
//		}
//		if (k < nrt) {
//			// Compute the k-th row transformation and place the
//			// k-th super-diagonal in e[k].
//			// Compute 2-norm without under/overflow.
//			e[k] = 0;
//			for (int i = k+1; i < n; i++) {
//				e[k] = hypot(e[k], e[i]);
//			}
//			if (e[k] != 0.0) {
//				if (e[k+1] < 0.0) {
//					e[k] = -e[k];
//				}
//				for (int i = k+1; i < n; i++) {
//					e[i] /= e[k];
//				}
//				e[k+1] += 1.0;
//			}
//			e[k] = -e[k];
//			if ((k+1 < m) & (e[k] != 0.0)) {
//				// Apply the transformation.
//				for (int i = k+1; i < m; i++) {
//					work[i] = 0.0;
//				}
//				for (int j = k+1; j < n; j++) {
//					for (int i = k+1; i < m; i++) {
//						work[i] += e[j] * A[i][j];
//					}
//				}
//				for (int j = k+1; j < n; j++) {
//					double t = -e[j]/e[k+1];
//					for (int i = k+1; i < m; i++) {
//						A[i][j] += t * work[i];
//					}
//				}
//			}
//			if (wantv) {
//				// Place the transformation in V for subsequent
//				// back multiplication.
//				for (int i = k+1; i < n; i++) {
//					V[i][k] = e[i];
//				}
//			}
//		}
//	}
//
//	cout<<"Here U\n";
//	for( int j =0; j<n; j++){
//		for (int k = 0; k< n ; k++)
//			cout<<U[j][k]<<", ";
//		cout<<endl;
//	}
//
//	// Set up the final bidiagonal matrix or order p.
//	int p = iMin(n, m+1);
//	if (nct < n) {
//		s[nct] = A[nct][nct];
//	}
//	if (m < p) {
//		s[p-1] = 0.0;
//	}
//	if (nrt+1 < p) {
//		e[nrt] = A[nrt][p-1];
//	}
//	e[p-1] = 0.0;
//
//	// If required, generate U.
//	if (wantu) {
//		for (int j = nct; j < nu; j++) {
//			for (int i = 0; i < m; i++) {
//				U[i][j] = 0.0;
//			}
//			U[j][j] = 1.0;
//		}
//		for (int k = nct-1; k >= 0; k--) {
//			if (s[k] != 0.0) {
//				for (int j = k+1; j < nu; j++) {
//					double t = 0;
//					for (int i = k; i < m; i++) {
//						t += U[i][k]*U[i][j];
//					}
//					t = -t / U[k][k];
//					for (int i = k; i < m; i++) {
//						U[i][j] += t*U[i][k];
//					}
//				}
//				for (int i = k; i < m; i++ ) {
//					U[i][k] = -U[i][k];
//				}
//				U[k][k] = 1.0 + U[k][k];
//				for (int i = 0; i < k-1; i++) {
//					U[i][k] = 0.0;
//				}
//			} else {
//				for (int i = 0; i < m; i++) {
//					U[i][k] = 0.0;
//				}
//				U[k][k] = 1.0;
//			}
//		}
//	}
//
//	// If required, generate V.
//	if (wantv) {
//		for (int k = n-1; k >= 0; k--) {
//			if ((k < nrt) & (e[k] != 0.0)) {
//				for (int j = k+1; j < nu; j++) {
//					double t = 0;
//					for (int i = k+1; i < n; i++) {
//						t += V[i][k] * V[i][j];
//					}
//					t = -t/V[k+1][k];
//					for (int i = k+1; i < n; i++) {
//						V[i][j] += t * V[i][k];
//					}
//				}
//			}
//			for (int i = 0; i < n; i++) {
//				V[i][k] = 0.0;
//			}
//			V[k][k] = 1.0;
//		}
//	}
//
//	// Main iteration loop for the singular values.
//	int pp = p-1;
//	int iter = 0;
//	const double epsL = pow(2.0, -52.0);
//	while (p > 0) {
//		int k,kase;
//
//		// Here is where a test for too many iterations would go.
//
//		// This section of the program inspects for
//		// negligible elements in the s and e arrays.  On
//		// completion the variables kase and k are set as follows.
//
//		// kase = 1     if s(p) and e[k-1] are negligible and k<p
//		// kase = 2     if s(k) is negligible and k<p
//		// kase = 3     if e[k-1] is negligible, k<p, and
//		//              s(k), ..., s(p) are not negligible (qr step).
//		// kase = 4     if e(p-1) is negligible (convergence).
//
//		for (k = p-2; k >= -1; k--) {
//			if (k == -1) {
//				break;
//			}
//			if (fabs(e[k]) <= eps*(fabs(s[k]) + fabs(s[k+1]))) {
//				e[k] = 0.0;
//				break;
//			}
//		}
//		if (k == p-2) {
//			kase = 4;
//		} else {
//			int ks;
//			for (ks = p-1; ks >= k; ks--) {
//				if (ks == k) {
//					break;
//				}
//				double t = (ks != p ?  fabs(e[ks]) : 0.0) + (ks != k+1 ? fabs(e[ks-1]) : 0.0);
//				if ( fabs(s[ks]) <= eps*t )  {
//					s[ks] = 0.0;
//					break;
//				}
//			}
//			if (ks == k) {
//				kase = 3;
//			} else if (ks == p-1) {
//				kase = 1;
//			} else {
//				kase = 2;
//				k = ks;
//			}
//		}
//		k++;
//
//		// Perform the task indicated by kase.
//		switch (kase) {
//		// Deflate negligible s(p).
//		case 1: {
//			double f = e[p-2];
//			e[p-2] = 0.0;
//			for (int j = p-2; j >= k; j--) {
//				double t = hypot(s[j], f );
//				double cs = s[j] / t;
//				double sn = f / t;
//				s[j] = t;
//				if (j != k) {
//					f = -sn*e[j-1];
//					e[j-1] = cs*e[j-1];
//				}
//				if (wantv) {
//					for (int i = 0; i < n; i++) {
//						t = cs*V[i][j] + sn*V[i][p-1];
//						V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
//						V[i][j] = t;
//					}
//				}
//			}
//		}
//		break;
//
//		// Split at negligible s(k).
//		case 2: {
//			double f = e[k-1];
//			e[k-1] = 0.0;
//			for (int j = k; j < p; j++) {
//				double t = hypot(s[j], f);
//				double cs = s[j] / t;
//				double sn = f / t;
//				s[j] = t;
//				f = -sn*e[j];
//				e[j] = cs*e[j];
//				if (wantu) {
//					for (int i = 0; i < m; i++) {
//						t = cs*U[i][j] + sn*U[i][k-1];
//						U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
//						U[i][j] = t;
//					}
//				}
//			}
//		}
//		break;
//
//		// Perform one qr step.
//		case 3: {
//			// Calculate the shift.
//			double scale = dMax( dMax( dMax( dMax (
//					fabs(s[p-1]), fabs(s[p-2])), fabs(e[p-2]) ),
//					fabs(s[k]) ), fabs(e[k]) );
//			double sp = s[p-1] / scale;
//			double spm1 = s[p-2] / scale;
//			double epm1 = e[p-2] / scale;
//			double sk = s[k] / scale;
//			double ek = e[k] / scale;
//			double b = ( (spm1 + sp)*(spm1 - sp) + epm1*epm1) / 2.0;
//			double c = (sp*epm1)*(sp*epm1);
//			double shift = 0.0;
//			if ((b != 0.0) | (c != 0.0)) {
//				shift = sqrt(b*b + c);
//				if (b < 0.0) {
//					shift = -shift;
//				}
//				shift = c/(b + shift);
//			}
//			double f = (sk + sp)*(sk - sp) + shift;
//			double g = sk*ek;
//
//			// Chase zeros.
//			for (int j = k; j < p-1; j++) {
//				double t = hypot(f, g);
//				double cs = f / t;
//				double sn = g / t;
//				if (j != k) {
//					e[j-1] = t;
//				}
//				f = cs*s[j] + sn*e[j];
//				e[j] = cs*e[j] - sn*s[j];
//				g = sn*s[j+1];
//				s[j+1] = cs*s[j+1];
//				if (wantv) {
//					for (int i = 0; i < n; i++) {
//						t = cs*V[i][j] + sn*V[i][j+1];
//						V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
//						V[i][j] = t;
//					}
//				}
//				t = hypot(f, g);
//				cs = f / t;
//				sn = g / t;
//				s[j] = t;
//				f = cs*e[j] + sn*s[j+1];
//				s[j+1] = -sn*e[j] + cs*s[j+1];
//				g = sn*e[j+1];
//				e[j+1] = cs*e[j+1];
//				if (wantu && (j < m-1)) {
//					for (int i = 0; i < m; i++) {
//						t = cs*U[i][j] + sn*U[i][j+1];
//						U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
//						U[i][j] = t;
//					}
//				}
//			}
//			e[p-2] = f;
//			iter = iter + 1;
//		}
//		break;
//
//		// Convergence.
//
//		case 4: {
//
//			// Make the singular values positive.
//			if (s[k] <= 0.0) {
//				s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
//				if (wantv) {
//					for (int i = 0; i <= pp; i++) {
//						V[i][k] = -V[i][k];
//					}
//				}
//			}
//
//			// Order the singular values.
//			while (k < pp) {
//				if (s[k] >= s[k+1]) {
//					break;
//				}
//				double t = s[k];
//				s[k] = s[k+1];
//				s[k+1] = t;
//				if (wantv && (k < n-1)) {
//					for (int i = 0; i < n; i++) {
//						t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
//					}
//				}
//				if (wantu && (k < m-1)) {
//					for (int i = 0; i < m; i++) {
//						t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
//					}
//				}
//				k++;
//			}
//			iter = 0;
//			p--;
//		}
//		break;
//		}
//	}
//}

// this routine will change A !!
svd::svd(vector<double>&  A, const int mm, const int nn){
	// Derived from LINPACK code.
	// Initialize.
	m =mm;
	n = nn;
	const int nu = iMin ( m, n );
	const int mu = iMin(m+1, n);
	s = vector<double> (mu, 0.0);
	U = vector<double> (m * nu, 0.0);
	double *e = new double [n];
	double *work = new double [m];
	for(int i = 0; i < m; i++){
		work[i] = 0.0;
//		U[i] = new double [nu];
//		for(int j = 0; j<nu; j++)
//			U[i][j] = 0.0;
	}
//	V = new double* [n];
	V = vector<double> (n * n, 0.0);
	for(int i = 0; i < n; i++){
		e[i] = 0.0;
//		V[i] = new double [n];
//		for(int j = 0; j<n; j++)
//			V[i][j] = 0.0;
	}

	bool wantu = true;
	bool wantv = true;
	// Reduce A to bidiagonal form, storing the diagonal elements
	// in s and the super-diagonal elements in e.
	const int nct = iMin(m-1, n);
	const int nrt = iMax(0, iMin(n-2, m) );

	for (int k = 0; k < iMax(nct, nrt); k++) {
		if (k < nct) {
			// Compute the transformation for the k-th column and
			// place the k-th diagonal in s.at(k).
			// Compute 2-norm of k-th column without under/overflow.
			s[k] = 0;
			for (int i = k; i < m; i++) {
				s[k] = hypot(s[k], A[ i * n + k ]);
			}
			if (s[k] != 0.0) {
				if (A[ k * n + k ] < 0.0) {
					s[k] = -s[k];
				}
				for (int i = k; i < m; i++) {
					A[ i * n + k ] /= s[k];
				}
				A[ k * n + k ] += 1.0;
			}
			s[k] = -s[k];
		}

		for (int j = k+1; j < n; j++) {
			if ((k < nct) & (s[k] != 0.0))  {
				// Apply the transformation.
				double t = 0;
				for (int i = k; i < m; i++) {
					t += A[ i * n + k ] * A[ i * n + j ];
				}
				t = -t/A[ k * n + k ];
				for (int i = k; i < m; i++) {
					A[ i * n + j ] += t * A[ i * n + k ];
				}
			}
			// Place the k-th row of A into e for the
			// subsequent calculation of the row transformation.
			e[j] = A[ k * n + j ];
		}
		if (wantu & (k < nct)) {
			// Place the transformation in U for subsequent back
			// multiplication.
			for (int i = k; i < m; i++) {
				U[ i * nu + k ] = A[ i * n + k ];
			}
		}
		if (k < nrt) {
			// Compute the k-th row transformation and place the
			// k-th super-diagonal in e[k].
			// Compute 2-norm without under/overflow.
			e[k] = 0;
			for (int i = k+1; i < n; i++) {
				e[k] = hypot(e[k], e[i]);
			}
			if (e[k] != 0.0) {
				if (e[k+1] < 0.0) {
					e[k] = -e[k];
				}
				for (int i = k+1; i < n; i++) {
					e[i] /= e[k];
				}
				e[k+1] += 1.0;
			}
			e[k] = -e[k];
			if ((k+1 < m) & (e[k] != 0.0)) {
				// Apply the transformation.
				for (int i = k+1; i < m; i++) {
					work[i] = 0.0;
				}
				for (int j = k+1; j < n; j++) {
					for (int i = k+1; i < m; i++) {
						work[i] += e[j] * A[ i * n + j ];
					}
				}
				for (int j = k+1; j < n; j++) {
					double t = -e[j] / e[k+1];
					for (int i = k+1; i < m; i++) {
						A[ i * n + j ] += t * work[i];
					}
				}
			}
			if (wantv) {
				// Place the transformation in V for subsequent
				// back multiplication.
				for (int i = k+1; i < n; i++) {
					V[ i * n + k ] = e[i];
				}
			}
		}
	}

	// Set up the final bidiagonal matrix or order p.
	int p = iMin(n, m+1);
	if (nct < n) {
		s[nct] = A[ nct * n + nct]; //A[nct][nct];
	}
	if (m < p) {
		s[p-1] = 0.0;
	}
	if (nrt+1 < p) {
		e[nrt] = A[ nrt * n + p -1]; // A[nrt][p-1];
	}
	e[p-1] = 0.0;

	// If required, generate U.
	if (wantu) {
		for (int j = nct; j < nu; j++) {
			for (int i = 0; i < m; i++) {
				U[ i * nu + j ]= 0.0; //U[i][j] = 0.0;
			}
			U[ j * nu + j ] = 1.0; //U[j][j] = 1.0;
		}
		for (int k = nct-1; k >= 0; k--) {
			if (s[k] != 0.0) {
				for (int j = k+1; j < nu; j++) {
					double t = 0;
					for (int i = k; i < m; i++) {
						t += U[ i * nu + k] * U[ i * nu + j ];
					}
					t = -t / U[ k * nu + k ];
					for (int i = k; i < m; i++) {
						U[ i * nu + j ] += t * U[ i * nu + k];
					}
				}
				for (int i = k; i < m; i++ ) {
					U[ i * nu + k] = -U[ i * nu + k];
				}
				U[ k * nu + k ] = 1.0 + U[ k * nu + k ];
				for (int i = 0; i < k-1; i++) {
					U[ i * nu + k] = 0.0;
				}
			} else {
				for (int i = 0; i < m; i++) {
					U[ i * nu + k] = 0.0;
				}
				U[ k * nu + k ] = 1.0;
			}
		}
	}

	// If required, generate V.
	if (wantv) {
		for (int k = n-1; k >= 0; k--) {
			if ((k < nrt) & (e[k] != 0.0)) {
				for (int j = k+1; j < nu; j++) {
					double t = 0;
					for (int i = k+1; i < n; i++) {
						t += V[ i * n + k ] * V[ i * n + j ];
					}
					t = -t / V[ k * n +n+ k ];
					for (int i = k+1; i < n; i++) {
						V[ i * n + j ] += t * V[ i * n + k ];
					}
				}
			}
			for (int i = 0; i < n; i++) {
				V[ i * n + k ] = 0.0;
			}
			V[ k * n + k ] = 1.0;
		}
	}

	// Main iteration loop for the singular values.
	int pp = p-1;
	int iter = 0;
	const double epsL = pow(2.0, -52.0);
	while (p > 0) {
		int k,kase;

		// Here is where a test for too many iterations would go.

		// This section of the program inspects for
		// negligible elements in the s and e arrays.  On
		// completion the variables kase and k are set as follows.

		// kase = 1     if s(p) and e[k-1] are negligible and k<p
		// kase = 2     if s(k) is negligible and k<p
		// kase = 3     if e[k-1] is negligible, k<p, and
		//              s(k), ..., s(p) are not negligible (qr step).
		// kase = 4     if e(p-1) is negligible (convergence).

		for (k = p-2; k >= -1; k--) {
			if (k == -1) {
				break;
			}
			if (fabs(e[k]) <= epsL*(fabs(s[k]) + fabs(s[ k+1] ))) {
				e[k] = 0.0;
				break;
			}
		}
		if (k == p-2) {
			kase = 4;
		} else {
			int ks;
			for (ks = p-1; ks >= k; ks--) {
				if (ks == k) {
					break;
				}
				double t = (ks != p ?  fabs(e[ks]) : 0.0) + (ks != k+1 ? fabs(e[ks-1]) : 0.0);
				if ( fabs(s[ks]) <= epsL*t )  {
					s[ks] = 0.0;
					break;
				}
			}
			if (ks == k) {
				kase = 3;
			} else if (ks == p-1) {
				kase = 1;
			} else {
				kase = 2;
				k = ks;
			}
		}
		k++;

		// Perform the task indicated by kase.
		switch (kase) {
		// Deflate negligible s(p).
		case 1: {
			double f = e[p-2];
			e[p-2] = 0.0;
			for (int j = p-2; j >= k; j--) {
				double t = hypot(s[ j ], f );
				double cs = s[ j ] / t;
				double sn = f / t;
				s[ j ] = t;
				if (j != k) {
					f = -sn*e[j-1];
					e[j-1] = cs*e[j-1];
				}
				if (wantv) {
					for (int i = 0; i < n; i++) {
						t = cs * V[ i * n + j ] + sn * U[ i * n + p - 1 ];
						U[ i * n + p - 1 ] = -sn*V[ i * n + j ] + cs*U[ i * n + p - 1 ];
						V[ i * n + j ] = t;
					}
				}
			}
		}
		break;

		// Split at negligible s(k).
		case 2: {
			double f = e[k-1];
			e[k-1] = 0.0;
			for (int j = k; j < p; j++) {
				double t = hypot(s[ j ], f);
				double cs = s[ j ] / t;
				double sn = f / t;
				s[ j ] = t;
				f = -sn*e[j];
				e[j] = cs*e[j];
				if (wantu) {
					for (int i = 0; i < m; i++) {
						t = cs*U[ i * nu + j ] + sn*U[ i * nu + k - 1 ];
						U[ i * nu + k - 1 ] = -sn*U[ i * nu + j ] + cs*U[ i * nu + k - 1 ];
						U[ i * nu + j ] = t;
					}
				}
			}
		}
		break;

		// Perform one qr step.
		case 3: {
			// Calculate the shift.
			double scale = dMax( dMax( dMax( dMax (
					fabs(s[p-1]), fabs(s[p-2])), fabs(e[p-2]) ),
					fabs(s[k]) ), fabs(e[k]) );
			double sp = s[p-1] / scale;
			double spm1 = s[p-2] / scale;
			double epm1 = e[p-2] / scale;
			double sk = s[k] / scale;
			double ek = e[k] / scale;
			double b = ( (spm1 + sp)*(spm1 - sp) + epm1*epm1) / 2.0;
			double c = (sp*epm1)*(sp*epm1);
			double shift = 0.0;
			if ((b != 0.0) | (c != 0.0)) {
				shift = sqrt(b*b + c);
				if (b < 0.0) {
					shift = -shift;
				}
				shift = c/(b + shift);
			}
			double f = (sk + sp)*(sk - sp) + shift;
			double g = sk*ek;

			// Chase zeros.
			for (int j = k; j < p-1; j++) {
				double t = hypot(f, g);
				double cs = f / t;
				double sn = g / t;
				if (j != k) {
					e[j-1] = t;
				}
				f = cs*s[ j ] + sn*e[j];
				e[j] = cs*e[j] - sn*s[ j ];
				g = sn*s[ j + 1 ];
				s[ j + 1 ] = cs*s[ j + 1 ];
				if (wantv) {
					for (int i = 0; i < n; i++) {
						t = cs*V[ i * n + j ] + sn*V[ i * n + j  + 1 ]; //[i][j+1];
						V[ i * n + j  + 1 ] = -sn*V[ i * n + j ] + cs* V[ i * n + j  + 1 ]; //V[i][j+1];
						V[ i * n + j ] = t;
					}
				}
				t = hypot(f, g);
				cs = f / t;
				sn = g / t;
				s[ j ] = t;
				f = cs*e[j] + sn*s[ j+1];
				s[ j+1] = -sn*e[j] + cs*s[ j+1];
				g = sn*e[j+1];
				e[j+1] = cs*e[j+1];
				if (wantu && (j < m-1)) {
					for (int i = 0; i < m; i++) {
						t = cs * U[ i * nu + j ] + sn * U[ i * nu + j + 1 ];
						U[ i * nu + j + 1 ] = -sn*U[ i * nu + j ] + cs*U[ i * nu + j + 1 ];
						U[ i * nu + j ] = t;
					}
				}
			}
			e[p-2] = f;
			iter = iter + 1;
		}
		break;

		// Convergence.

		case 4: {

			// Make the singular values positive.
			if (s[k] <= 0.0) {
				s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
				if (wantv) {
					for (int i = 0; i <= pp; i++) {
						V[ i * n + k ] = -V[ i * n + k ];
					}
				}
			}

			// Order the singular values.
			while (k < pp) {
				if (s[k] >= s[ k+1 ]) {
					break;
				}
				double t = s[k];
				s[k] = s[ k+1 ];
				s[ k+1 ] = t;
				if (wantv && (k < n-1)) {
					for (int i = 0; i < n; i++) {
						t = V[ i * n + k + 1 ]; //[i][k+1];
						V[ i * n + k + 1 ] = V[ i * n + k ];
						V[ i * n + k ] = t;
					}
				}
				if (wantu && (k < m-1)) {
					for (int i = 0; i < m; i++) {
						t = U[ i * nu + k + 1 ]; U[ i * nu + k + 1 ] = U[ i * nu + k]; U[ i * nu + k] = t;
					}
				}
				k++;
			}
			iter = 0;
			p--;
		}
		break;
		}
	}
	delete[] e;
	delete[] work;
}

matrix svd::getMatrixS () const {
	vector<double>A(n * m, 0.0);
	for (int i = 0; i < m; i++)
		A[i * n + i] = s [i];
	return matrix(A, m, n);
}

svd::~svd(){
	U.clear();
	V.clear();
//	s.clear();
//	for(int i =0; i < m; i++)
//		delete [] U[i];
//	delete [] U;
//	for(int i =0; i < n; i++)
//		delete [] V[i];
//	delete [] V;
//	delete [] s;
	;
}
