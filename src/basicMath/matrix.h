/**
 * Modified from a Java class previously written by Lincong Wang(Myself)
 * The Java class itself is modified from the Jama numeric package
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include<string>
#include<vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "Math.h"
using namespace std;

class matrix{
public:

	matrix();
	matrix (const size_t m, const size_t n) ;
	matrix (const size_t m, const size_t n, double s);
	matrix (double **A);
	matrix (float **A, const size_t m,  const size_t n) ;
	matrix (double **A, const size_t m, const size_t n) ;
	matrix (const vector<double>& a, const size_t mm=3, const size_t nn=3)
		: A(a),  m(mm), n(nn){;}
	matrix (double *vals, const size_t m, const size_t size) ;
	matrix (const matrix&);
	matrix &operator=(const matrix& );

	vector<double>  getArray () const {return A;}

	bool getARow(const size_t r, vector<double>& row ){
		if ( r >= m )
			return false;
		for (unsigned int j = 0; j < n ; j++)
			row[ j ] = A[ r * n + j ] ;
		return true;
	}

	size_t getRowDimension () const {return m;}
	size_t getColumnDimension () const {return n;}
	double getElement(const size_t i, const size_t j) const{ return A.at(n*i+j);}
	double** getArrayCopy() const;  // return should be changed to vector<double>
	matrix rotationMat(double, string);
	matrix rotationMatCS(double, string);
	matrix eulerMat (double alpha, double beta, double gamma);

	matrix eulerMat (double sinAlpha, double cosAlpha, double sinBeta, double cosBeta,
			double sinGamma, double cosGamma);

	matrix rotate2XYplane(double *normalUnitVector, const size_t dim, bool);

	void set (const int i, const int j, const double s) {A.at(n*i+j) = s;}

	matrix transpose() const;
	matrix plus (const matrix&) const ;
	matrix minus (const matrix&)  const;
	matrix times (const double) const;

	//	matrix times (matrix*)  const;
	matrix times (const matrix& B)  const {
		if (B.m != n) {
			cerr<<"matrix inner dimensions must agree B:."<<endl;
			exit(0);
		}
		unsigned int nb = B.n;
		vector<double>C(m * nb , 0.0) ;
		for (unsigned int i = 0; i < m; i++)
			for (unsigned int j = 0; j < nb; j++)
				for (unsigned int k = 0; k < n; k++)
					C[ i * nb + j ] += A [ i * n + k ] * B.A[ k * nb + j ];  //C_{ij} = A_{ik} * B_{kj}
		return matrix(C, m, nb);
	}

	//	matrix times (matrix*)  const;
	bool times (const matrix& B, matrix &mat)  const {
		if (B.m != n) {
			cerr<<"matrix inner dimensions must agree BB:."<<endl;
			return false;
		}
		unsigned int nb = B.n;
		vector<double>C(m * nb , 0.0) ;
		for (unsigned int i = 0; i < m; i++)
			for (unsigned int j = 0; j < nb; j++)
				for (unsigned int k = 0; k < n; k++)
					C[ i * nb + j ] += A [ i * n + k ] * B.A[ k * nb + j ];  //C_{ij} = A_{ik} * B_{kj}
		mat = matrix(C, m, nb);
		return true;
	}

	matrix changeHandness() const;
	void changeHandness2(){
		for(unsigned int i = 0; i < m; i++){ A [ n*i ] = -A [ n * i ]; A [ n*i  + 2 ] = A [ n*i + 2] ;} }

	vector<double> times(const double*, const size_t =dim3) const;

	double* times(double* const, const size_t ) const;

	bool times(const double*, const size_t, vector<double>& ) const;

	vector<float> times( const float*, const size_t =dim3) const;

	bool times(const float*, const size_t, vector<float>& ) const;

	vector<double> times (const vector<double>& B, const size_t len=dim3) const{
		if ( len != n) {
			cerr<<"matrix inner dimensions must agree."<<endl;
			exit(0);
		}
		vector<double> vec(m, 0.0);
		for (unsigned int j = 0; j < m; j++)
			for (unsigned int k = 0; k < n; k++)
				vec [ j ] += A [ j * n + k ] * B [ k ]; //A[j][k] * B[k];
		return vec;
	}

	bool times (const vector<double>& B, const size_t len, vector<double>& vec) const{
		if ( len != n || vec.size() != m ) {
			cerr<<"matrix inner dimensions must agree."<<endl;
			return false;
		}
		for (unsigned int j = 0; j < m; j++)
			for (unsigned int k = 0; k < n; k++)
				vec [ j ] += A [ j * n + k ] * B [ k ]; //A[j][k] * B[k];
		return true;
	}

	vector<float> times (const vector<float>& B, const size_t len=dim3) const{
		if ( len != n) {
			cerr<<"matrix inner dimensions must agree."<<endl;
			exit( 0 );
		}
		vector<float> vec(m, 0.0);
		for (unsigned int j = 0; j < m; j++)
			for (unsigned int k = 0; k < n; k++)
				vec [ j ]+= A [ j * n + k ] * B [ k ]; //A[j][k] * B[k];
		return vec;
	}

	bool times (const vector<float>& B, const size_t len, vector<float>& vec) const{
		if ( len != n || vec.size() != m ) {
			cerr<<"matrix inner dimensions must agree."<<endl;
			return false;
		}
		for (unsigned int j = 0; j < m; j++)
			for (unsigned int k = 0; k < n; k++)
				vec [ j ] += A [ j * n + k ] * B [ k ]; //A[j][k] * B[k];
		return true;
	}

	matrix uminus() const;
	double norm1();
	double trace();
	double det(const matrix& ) const;

	void print() const;

	virtual ~matrix();

private:
	vector<double> A;    //a m-by-n matrix, the index will be
	size_t m;     //row m and column n.
	size_t n;
};

//static matrix identity(unsigned int, unsigned int);
//
static matrix rotIdentity (const size_t m, const size_t n) {
	vector<double> A (m*n );
	for (unsigned int i = 0; i < m; i++) {
		for (unsigned int j = 0; j < n; j++) {
			A [ i *n + j ]= (i == j ? 1.0 : 0.0);
		}
	}
	return matrix(A, m, n);
}

static matrix rotationMatCS(double cosTheta, string axis) {
	vector<double> A(9, 0.0); // 3*3 - 1
	if (axis=="+x"){
		A [ 0] =  1.0;
		A [ 4 ]= cosTheta;
		A [ 5] =  sqrt(1.0 - cosTheta* cosTheta );
		A [ 7] =  -A[5];
		A [ 8] =  cosTheta;
	} else if (axis=="-x"){
		A [ 0] =  1.0;
		A [ 4] =  cosTheta;
		A [ 5] =  - sqrt(1.0 - cosTheta* cosTheta );
		A [ 7] =  -A[5];
		A [ 8] =  cosTheta;
	}else if (axis=="+y"){
		A [ 0] =  cosTheta;
		A [ 2] =  -sqrt(1.0 - cosTheta* cosTheta );
		A [ 4] =  1.0;
		A [ 6] =  - A [ 2];
		A [ 8] =  cosTheta;
	}else if (axis=="-y"){
		A [ 0] =  cosTheta;
		A [ 2] =  sqrt(1.0 - cosTheta* cosTheta );
		A [ 4] =  1.0;
		A [ 6] =  - A [ 2];
		A [ 8] =  cosTheta;
	}else if (axis=="+z"){
		A [ 0] =  cosTheta;
		A [ 1] =  sqrt(1.0 - cosTheta* cosTheta );
		A [ 3] =  -A[1];
		A [ 4] =  cosTheta;
		A [ 8] =  1.0;
	}else if (axis=="-z"){
		A [ 0] =  cosTheta;
		A [ 1] =  -sqrt(1.0 - cosTheta* cosTheta );
		A [ 3] =  -A[1];
		A [ 4] =  cosTheta;
		A [ 8] =  1.0;
	}
	return matrix(A, 3, 3);
}

inline static matrix rotationMatCS(const double sinTheta, const double cosTheta, string axis) {
	vector<double> A(9, 0.0); // 3*3 - 1
	if (axis=="+x"){  // +y --> +z
		A [ 0] =  1.0;
		A [ 4 ]= cosTheta;
		A [ 5] = sinTheta;
		A [ 7] =  -A[5];
		A [ 8] =  cosTheta;
	} else if (axis=="-x"){
		A [ 0] =  1.0;
		A [ 4] =  cosTheta;
		A [ 5] =  - sinTheta;
		A [ 7] =  -A[5];
		A [ 8] =  cosTheta;
	}else if (axis=="+y"){   // the angle is from +x --> +z
		A [ 0] =  cosTheta;
		A [ 2] =  -sinTheta;
		A [ 4] =  1.0;
		A [ 6] =  - A [ 2];
		A [ 8] =  cosTheta;
	}else if (axis=="-y"){
		A [ 0] =  cosTheta;
		A [ 2] =  sinTheta;
		A [ 4] =  1.0;
		A [ 6] =  - A [ 2];
		A [ 8] =  cosTheta;
	}else if (axis=="+z"){  //+x --> +y
		A [ 0] =  cosTheta;
		A [ 1] =  sinTheta;
		A [ 3] =  -A[1];
		A [ 4] =  A [ 0] ;
		A [ 8] =  1.0;
	}else if (axis=="-z"){
		A [ 0] =  cosTheta;
		A [ 1] =  -sinTheta;
		A [ 3] =  -A[1];
		A [ 4] =  cosTheta;
		A [ 8] =  1.0;
	}
	return matrix(A, 3, 3);
}

/**
For some applications, it is helpful to be able to make a rotation with a given axis.
 Given a unit vector u = (ux, uy, uz), where ux2 + uy2 + uz2 = 1,
 the matrix for a rotation by an angle of θ about an axis in the direction of u is

    R = \begin{bmatrix} \cos \theta +u_x^2 \left(1-\cos \theta\right) & u_x u_y
    \left(1-\cos \theta\right) - u_z \sin \theta & u_x u_z \left(1-\cos \theta\right)
    + u_y \sin \theta \\ u_y u_x \left(1-\cos \theta\right) + u_z \sin \theta & \cos \theta
    + u_y^2\left(1-\cos \theta\right) & u_y u_z \left(1-\cos \theta\right) - u_x \sin \theta
    \\ u_z u_x \left(1-\cos \theta\right) - u_y \sin \theta & u_z u_y \left(1-\cos \theta\right)
    + u_x \sin \theta & \cos \theta + u_z^2\left(1-\cos \theta\right) \end{bmatrix}

This can be written more concisely as

    R = \mathbf{u}\otimes\mathbf{u} + \cos\theta(I-\mathbf{u}\otimes\mathbf{u})
    + \sin\theta[\mathbf u]_{\times},

where [\mathbf u]_{\times} is the skew symmetric form of u, ⊗ is the tensor product
and I is the Identity matrix. This is a matrix form of Rodrigues' rotation formula, with

    \mathbf{u}\otimes\mathbf{u} = \begin{bmatrix} u_x^2 & u_x u_y & u_x u_z
    \\[3pt] u_x u_y & u_y^2 & u_y u_z \\[3pt] u_x u_z & u_y u_z & u_z^2 \end
    {bmatrix},\qquad [\mathbf u]_{\times} = \begin{bmatrix} 0 & -u_z & u_y
    \\[3pt] u_z & 0 & -u_x \\[3pt] -u_y & u_x & 0 \end{bmatrix}.

If the 3D space is right-handed, this rotation will be counterclockwise for an observer
placed so that the axis u goes in his or her direction (Right-hand rule).
*/
/**
 * c = cos(θ); s = sin(θ); C = 1-c
 * [ xxC+c   xyC-zs  xzC+ys ]
[ yxC+zs  yyC+c   yzC-xs ]
[ zxC-ys  zyC+xs  zzC+c  ]
 *
 */

inline static matrix rotMatAlongAnAxis(const double sinTheta, const double cosTheta,
		const vector<double> &coord) {
	vector<double> A(9, 0.0);
	double C = 1.0 - cosTheta;
	A [0] =  coord[0] * coord[0] * C + cosTheta;
	A [4] =  coord[1] * coord[1] * C + cosTheta;
	A [8] =  coord[2] * coord[2] * C + cosTheta;

	A [1] =  coord[0] * coord[1] * C - coord[2] * sinTheta;
	A [3] =  coord[1] * coord[0] * C + coord[2] * sinTheta;

	A [2] =  coord[0] * coord[2] * C + coord[1] * sinTheta;
	A [6] =  coord[2] * coord[0] * C  - coord[1] * sinTheta;

	A [5] =  coord[1] * coord[2] * C  - coord[0] * sinTheta;
	A [7] =  coord[2] * coord[1] * C + coord[0] * sinTheta;

	return matrix(A, 3, 3);
}

inline static matrix rotMatAlongAnAxis(const double sinTheta, const double cosTheta,
		double *const coord) {
	vector<double> A(9, 0.0);
	double C = 1.0 - cosTheta;
	A [0] =  coord[0] * coord[0] * C + cosTheta;
	A [4] =  coord[1] * coord[1] * C + cosTheta;
	A [8] =  coord[2] * coord[2] * C + cosTheta;

	A [1] =  coord[0] * coord[1] * C - coord[2] * sinTheta;
	A [3] =  coord[1] * coord[0] * C + coord[2] * sinTheta;

	A [2] =  coord[0] * coord[2] * C + coord[1] * sinTheta;
	A [6] =  coord[2] * coord[0] * C  - coord[1] * sinTheta;

	A [5] =  coord[1] * coord[2] * C  - coord[0] * sinTheta;
	A [7] =  coord[2] * coord[1] * C + coord[0] * sinTheta;

	return matrix(A, 3, 3);
}


/** matrix transpose.
@return    A', a new matrix
 */
inline matrix matrix::transpose () const {
	  vector<double> C ( m * n );
	  for (unsigned int i = 0; i < m; i++)
		  for (unsigned int j = 0; j < n; j++)
			  C[ j * m + i ] = A [ i * n + j ];
	  return matrix(C, n, m);
}

inline double* matrix::times(double* const v, const size_t len=dim3 ) const{
	if ( len != n) {
		cerr<<"matrix inner dimensions must agree."<<endl;
		exit(0);
	}
	double *vec = new double[m];
	for (unsigned int j = 0; j < m; j++) {
		vec[j] = 0.0;
		for (unsigned int k = 0; k < n; k++) {
			vec [ j ]+= A [ j * n + k ] * v[k]; //A[j][k] * B[k];
		}
	}
	return vec;
}

inline vector<float> matrix::times ( const float* B, const size_t len) const {
	 if ( len != n) {
		 cerr<<"matrix inner dimensions must agree."<<endl;
		 exit(0);
	 }
	 vector<float> vec(m, 0.0);
	 for (unsigned int j = 0; j < m; j++) {
		 for (unsigned int k = 0; k < n; k++) {
			 vec [ j ]+= A [ j * n + k ] * B[k]; //A[j][k] * B[k];
		 }
	 }
	 return vec;
 }

inline bool matrix::times ( const float* B, const size_t len, vector<float>& vec) const {
	 if ( len != n || vec.size() != m ) {
		 cerr<<"matrix inner dimensions must agree."<<endl;
		 return false;
	 }
	 for (unsigned int j = 0; j < m; j++) {
		 for (unsigned int k = 0; k < n; k++) {
			 vec [ j ]+= A [ j * n + k ] * B[k]; //A[j][k] * B[k];
		 }
	 }
	 return true;
 }

inline vector<double> matrix::times (const double* B, const size_t len) const {
	if ( len != n) {
		cerr<<"matrix inner dimensions must agree."<<endl;
		exit(0);
	}
	vector<double> vec(m, 0.0);
	for (unsigned int j = 0; j < m; j++) {
		for (unsigned int k = 0; k < n; k++) {
			vec [ j ] += A [ j * n + k ] * B[k]; //A[j][k] * B[k];
		}
	}
	return vec;
}

inline bool matrix::times (const double* B, const size_t len, vector<double>& vec) const {
	if ( len != n || vec.size() != m ) {
		cerr<<"matrix inner dimensions must agree."<<endl;
		return false;
	}
	for (unsigned int j = 0; j < m; j++) {
		vec [ j ]  = 0.0;
		for (unsigned int k = 0; k < n; k++)
			vec [ j ] += A [ j * n + k ] * B[k]; //A[j][k] * B[k];
	}
	return true;
}

inline matrix matrix::eulerMat (double sinAlpha, double cosAlpha, double sinBeta, double cosBeta,
		double sinGamma, double cosGamma){
	vector<double> mat(9);
	mat.at(0) = cosAlpha * cosBeta * cosGamma - sinAlpha * sinGamma;
	mat.at(1) = sinAlpha * cosBeta * cosGamma + cosAlpha * sinGamma;
	mat.at(2) = -sinBeta * cosGamma;
	mat.at(3) = -cosAlpha * cosBeta * sinGamma - sinAlpha * cosGamma;
	mat.at(4) = -sinAlpha * cosBeta * sinGamma + cosAlpha * cosGamma;
	mat.at(5) = sinGamma * sinBeta;
	mat.at(6) = sinBeta * cosAlpha;
	mat.at(7) = sinAlpha * sinBeta;
	mat.at(8) = cosBeta;
	return (matrix(mat, 3, 3));
}

inline matrix matrix::rotate2XYplane(double *const normalUnitVector, const size_t dim, bool reverse){

	if ( dim != 3 ){
		cerr<<"rotate2XYplane requires the dimension to be 3\n";
		exit(0);
	}
	double cosX = normalUnitVector[0];
	double cosY = normalUnitVector[1];
	double cosZ = normalUnitVector[2];

	double cosBeta = cosZ;
	double sinBeta = sqrt(1.0 - cosBeta * cosBeta);
	double sinGamma = 0.00;
	double cosGamma = 1.00;
	if( sinBeta > 10.0 * Tiny){
		sinGamma = cosY / sinBeta;
		cosGamma = - cosX / sinBeta;
	}
	double sinAlpha = 0.0;
	double cosAlpha = 1.0;
	matrix rot = eulerMat(sinAlpha, cosAlpha, sinBeta, cosBeta, sinGamma, cosGamma);
	if(reverse)
		return rot.transpose();
	else
		return rot;
}
//
//inline void translation( const vector<float>& originCoordVec, vector<float>& newCoordVec,
//		const vector<float>& center)































//static matrix rotIdentity (const size_t m, const size_t n) ;

#endif /*MATRIX_H_*/
