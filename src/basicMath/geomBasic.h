/**
 * Implementation of some simple methods for 2D and 3D geometry
 */

#ifndef GEOMBASIC_H_
#define GEOMBASIC_H_

//The following is ugly but I have no choice
#define eps1 1.0E-12
#define eps22 1.0E-11
#define eps4 1.0E-10
#define eps5 5.1E-5  // since the atom coordinates is originally stored using three digits
#define eps7 1.0E-7

#ifndef mathPI
#define mathPI 3.141592653589793
#endif

#ifndef CST
#define CST  180.00 / 3.141592653589793
#endif

#ifndef toRad
#define toRad  3.141592653589793 / 180.00
#endif

#include <cmath>
#include <iostream>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>

#include "Math.h"
#include "matrix.h"

using namespace std;

class geomPlane{
public:
	geomPlane(const int pid){ planeId = pid;};
	geomPlane(  const int pid,
			const vector<float>& c, bool onXY, const vector<double>& cc,
			const vector<double>& ee,	const double r,
			double a = 0.0, bool renorm=false )
	: planeId(pid), coords(c), onXYplane(onXY), center(cc), coef(ee), residual2Plane(r),
	  area(a){
		norm = vector<double> ( dim3 );
		if(renorm) {
			double length = coef[0] * coef[0] +coef[1]*coef[1] + coef[2]*coef[2];
			norm[0] = coef[0] / length;
			norm[1] = coef[1] / length;
			norm[2] = coef[2] / length;
		}else {
			norm[0] = coef[0] ;
			norm[1] = coef[1] ;
			norm[2] = coef[2];
		}
	};

	vector<double> getCoef() const { return coef; };
	vector<double> getNorm() const { return norm; };

	vector<double> normalize () {
		double length = coef[0] * coef[0] +coef[1]*coef[1] +coef[2]*coef[2];
		norm = vector<double> ( dim3 );
		norm[0] = coef[0] / length;
		norm[1] = coef[1] / length;
		norm[2] =  coef[2]  / length;
		return norm;
	}

	float getRelativeTiltAngle() const { return relativeTiltAngle; };
	void  setRelativeTiltAngle(const float angle) { relativeTiltAngle = angle; };

	~geomPlane(){;};

public:

	// an index, typically will be residue number for both peptide plane and NT
	int planeId;

	// input of real data
	vector<float> coords;

	// output
	bool onXYplane;
	vector<double> center;
	double residual2Plane;

	// an optional parameter that may not be computed
	double area;
	// a tilt angle of the plane relative to its previous one
	// the first one is given the value 0.0;
	float relativeTiltAngle;

private:
	vector<double> coef;
	vector<double> norm;
};

struct gplaneIdComparatorPtr: public binary_function<geomPlane* , geomPlane*, bool>{
	bool operator() ( geomPlane* const e1, geomPlane* const e2 ) const {
		return ( e1->planeId < e2->planeId? true: false);
	}
};

class geomBasic{
public:
	const static size_t dim3 = 3;
//	static const double eps = 1.0E-14;
//	static const double mathPI = 3.141592653589793;
//	static const double cst = 180.00 / 3.141592653589793;

	geomBasic();
	virtual ~geomBasic();

//	bool colinear(double* p1, double* p2, double *p3, const size_t d);
//
//	bool colinear( const vector<double>& p1, const vector<double>& p2,
//			const vector<double> & p3, const size_t d);
//
//	bool colinear( const vector<float>& p1, const vector<float>& p2,
//			const vector<float> & p3);

	bool coplanar( double *p1, double *p2, double *p3, double *p4, const size_t d);
	bool coplanar( double **p, int size, const size_t d);
	double *crossProd(double *p1, double *p2, const size_t d);
	//bool crossProd(double *p1, double *p2, double *c, const size_t d);
	bool triangleArea(double *v1, double *v2, double *v3, double *area, const size_t d);
//	double dotProd(double *p1, double *p2, const size_t d);
//	double dotProd(const vector<double>& p1, const vector<double>& p2, const size_t d);
//	bool planeAngleP2Mid(double *p1, double* p2, double* p3, double *angle);
//	bool planeAngleP2Mid( const vector<double>& p1, const vector<double>& p2,
//			const vector<double>& p3, double& angle);
	// a x + b y + c z +d = 0.0
	double* coefOfPlaneEqu(double* p1, double* p2, double* p3, const size_t d);
	double dihedralAngle(double *p1, double *p2, double *p3, double *p4, const size_t d);

private:
	double *p1;
	double *p2;
	double *p3;
	double *p4;
	int dim;
};

inline double *internucleiVec(double *p1, double *p2, size_t d){
	double *v = new double[d];
	for(unsigned int i = 0; i<d; i++)
		v[i] = p2[i] - p1[i];
	return v;
}

inline float det ( const vector<float> &mat){
	if ( mat.size() != 9 ){
		cerr<<"The function only computes det of 3D,  but the dimension d ="<<mat.size()<<endl;
		exit(0);
	}
	return (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]

	       - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ]) ;
}

inline double det ( const vector<double> &mat ){
	return (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]

	       - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ]) ;
}

inline double det (const double *const mat, const size_t d){
	if ( d!=3){
		cerr<<"The function only computes det of 3D,  but the dimension d ="<<d<<endl;
		exit(0);
	}
	return (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]

	       - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ]) ;
}

inline double det (double **arrV, int d){
	if ( d!=3){
		cerr<<"The function only computes det of 3D,  but the dimension d ="<<d<<endl;
		exit(0);
	}

	return ( arrV[0][0] * arrV[1][1] * arrV[2][2] + arrV[1][0] * arrV[2][1] * arrV[0][2]
	          + arrV[2][0] * arrV[1][2] * arrV[0][1] - arrV[2][0] * arrV[1][1] * arrV[0][2]
	           - arrV[1][2] * arrV[2][1] * arrV[0][0] - arrV[2][2] * arrV[1][0] * arrV[0][1] );
}

inline double dotProd(const double *const p1, const double *const p2, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests points in 3D: but the dimension d ="<<d<<endl;
		exit(0);
	}
	double c = 0.0;
	for(unsigned int i = 0; i < d; i++)
		c +=  p1[i] * p2[i];
	return c;
}

inline float dotProd( const vector<float>& p1, const vector<float>& p2, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests points in 3D: but the dimension d ="<<d<<endl;
		exit(0);
	}
	float c = 0.0;
	for(unsigned int i = 0; i < d; i++)
		c +=  p1[i] * p2[i];
	return c;
}

inline double dotProd( const vector<double>& p1, const vector<double>& p2, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests points in 3D: but the dimension d ="<<d<<endl;
		exit(0);
	}
	double c = 0.0;
	for(unsigned int i = 0; i < d; i++)
		c +=  p1[i] * p2[i];
	return c;
}

/**
 * We need to move the origin to the center in order to be robust
 */
inline bool coplanar(double *p1, double *p2, double *p3, double* p4,const size_t d){
	if ( d!=3){
		cerr<<"The function only tests points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
	double **mat = new double *[d] ;
	for(unsigned int i = 0; i<d; i++)
		mat[i] = new double[d];
	mat[0][0] = p4[0] - p1[0];
	mat[0][1] = p4[1] - p1[1];
	mat[0][2] = p4[2] - p1[2];
	mat[1][0] = p3[0] - p1[0];
	mat[1][1] = p3[1] - p1[1];
	mat[1][2] = p3[2] - p1[2];
	mat[2][0] = p2[0] - p1[0];
	mat[2][1] = p2[1] - p1[1];
	mat[2][2] = p2[2] - p1[2];
	double detVal = det(mat, d);
	for(unsigned int i = 0; i<d; i++)
		delete[] mat[i] ;
	delete [] mat;
	if (detVal < eps1 )
		return true;
	return false;
}

inline bool coplanar(double **p, int size, const size_t d){
	if ( d!=3 || size < 4){
		cerr<<"The function only tests points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
	double **mat = new double *[d] ;
	for(unsigned int i = 0; i<d; i++)
		mat[i] = new double[d];
	double *p1=0, *p2=0, *p3=0, *p4=0;
	const int loopSize = size - 3;
	int k = 0;
	double detVal  = 0.0;
	bool isAPlane = true;
	while( k < loopSize ){
		p1 = p[ k ];
		p2 = p[ k + 1 ] ;
		p3 = p[ k + 2 ];
		p4 = p[ k + 3 ];
		mat[0][0] = p4[0] - p1[0];
		mat[0][1] = p4[1] - p1[1];
		mat[0][2] = p4[2] - p1[2];
		mat[1][0] = p3[0] - p1[0];
		mat[1][1] = p3[1] - p1[1];
		mat[1][2] = p3[2] - p1[2];
		mat[2][0] = p2[0] - p1[0];
		mat[2][1] = p2[1] - p1[1];
		mat[2][2] = p2[2] - p1[2];
		detVal = det(mat, d);
		if ( fabs(detVal) > eps1 ) {
			isAPlane = false;
		}
		k++;
	}
	for(unsigned int i = 0; i<d; i++)
		delete[] mat[i] ;
	delete [] mat;
	delete [] p1;
	delete [] p2;
	delete [] p3;
	delete [] p4;
	return isAPlane;
}

/**
 * It has a bug if one points is the origin
 */
inline bool colinear( const double *const p1, const double *const p2,
		const double *const p3, const size_t dim ){
	if ( dim != 3 ){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<dim<<endl;
		return false;
	}

	double mat[dim * dim];
	mat[0] = p1[0] ;
	mat[1] = p1[1] ;
	mat[2] = p1[2] ;
	mat[3] = p2[0] ;
	mat[4] = p2[1] ;
	mat[5] = p2[2] ;
	mat[6] = p3[0] ;
	mat[7] = p3[1] ;
	mat[8] = p3[2] ;
	double detVal = det(mat, dim);
	if (fabs(detVal) > eps1 )
		return false;

	double v1[dim3] = {p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]};
	double v2[dim3] = {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]};

	double norm1 = 0.0, norm2 = 0.0;
	for (unsigned int i = 0; i < dim; ++i) {
		norm1 += v1[i] * v1[i];
		norm2 += v2[i] * v2[i];
	}
	if ( norm1 < eps1 || norm2 < eps1 )
		return true;


	double c = dotProd(v1, v2, dim) / sqrt( norm1 * norm2 );
	if (fabs(c-1.0) < eps1)
		return true;
	if ( fabs(c+1.0)<eps1)
		return true;
	return false;
}

/**
 * It has a bug if one points is the origin, so we should move them
 */
inline bool colinear( const vector<float>& p1, const vector<float>& p2,
		const vector<float> & p3){

	const size_t d = p1.size();
	if ( d != 3 ){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
//	vector<double> center(d, 0.0);
//	for (unsigned int i = 0; i < 3; ++i)  {
//		center[i] += p1[i] + p2[i]+p3[i];
//	}
//	for (unsigned int i = 0; i < 3; ++i)
//			center[i] /= 3;


	vector<double> mat( d * d , 0.0);
	mat[0] = p1[0] ;
	mat[1] = p1[1] ;
	mat[2] = p1[2] ;

	mat[3] = p2[0];
	mat[4] = p2[1] ;
	mat[5] = p2[2] ;

	mat[6] = p3[0] ;
	mat[7] = p3[1] ;
	mat[8] = p3[2] ;

	double detVal = det(mat);

	if (fabs(detVal) > eps1 )
		return false;

	double v1[dim3] = {p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]};
	double v2[dim3] = {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]};

	double norm1 = 0.0, norm2 = 0.0;
	for (unsigned int i = 0; i < d; ++i) {
		norm1 += v1[i] * v1[i];
		norm2 += v2[i] * v2[i];
	}
	if ( norm1 < eps1 || norm2 < eps1 )
		return true;


	double c = dotProd(v1, v2, d) / sqrt( norm1 * norm2 );
	if (fabs(c-1.0) < eps1)
		return true;
	if ( fabs(c+1.0)<eps1)
		return true;
	return false;

}


/**
 * It has a bug if one points is the origin
 */
inline bool colinear( const vector<double>& p1, const vector<double>& p2,
		const vector<double> & p3, const size_t d){
	if ( d != 3 ){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
//	vector<double> center(d, 0.0);
//	for (unsigned int i = 0; i < 3; ++i)
//		center[i] += p1[i] + p2[i]+p3[i];
//
//	for (unsigned int i = 0; i < 3; ++i)
//		center[i] /= 3;

	double mat [d*d] ;
	mat[0] = p1[0] ;
	mat[1] = p1[1] ;
	mat[2] = p1[2] ;
	mat[3] = p2[0] ;
	mat[4] = p2[1] ;
	mat[5] = p2[2] ;
	mat[6] = p3[0] ;
	mat[7] = p3[1] ;
	mat[8] = p3[2] ;
	double detVal = det(mat, d);

	if (fabs(detVal) > eps1 )
		return false;

	double v1[dim3] = {p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]};
	double v2[dim3] = {p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]};

	double norm1 = 0.0, norm2 = 0.0;
	for (unsigned int i = 0; i < d; ++i) {
		norm1 += v1[i] * v1[i];
		norm2 += v2[i] * v2[i];
	}
	if ( norm1 < eps1 || norm2 < eps1 )
		return true;


	double c = dotProd(v1, v2, d) / sqrt( norm1 * norm2 );
	if (fabs(c-1.0) < eps1)
		return true;
	if ( fabs(c+1.0)<eps1)
		return true;
	return false;

}

inline vector<double>internucleiVec(const vector<double>&p1, const vector<double>&p2, size_t d){
	vector<double> v(d);
	for(unsigned int i = 0; i<d; i++)
		v[i] = p2[i] - p1[i];
	return v;
}

inline vector<float>internucleiVec(const vector<float>&p1, const vector<float>&p2, size_t d){
	vector<float> v(d);
	for(unsigned int i = 0; i<d; i++)
		v[i] = p2[i] - p1[i];
	return v;
}

//static int sgn(double x){
//	if( x > 0.0) return 1;
//	else if(x < 0.0) return -1;
//	else return 0;
//}

inline double norm(double* v, size_t d){
	double len = 0.0;
	for(unsigned int i = 0; i<d; i++)
			len += v[i]  * v[i];
	return (sqrt(len));
}

inline double norm( const vector<double>& v, size_t d){
	double len = 0.0;
	for(unsigned int i = 0; i<d; i++)
			len += v[i]  * v[i];
	return (sqrt(len));
}

inline bool crossProd(double *p1, double *p2, double *c, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	return true;
}

inline bool crossProd(double *p1, double *p2, double *c, const size_t d,
		bool normalize ) {
	if ( d!=3){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	if ( normalize ) {
		double r = sqrt( c[0] * c[0] + c[1] * c[1] + c[2] * c[2] ) ;
		c[0] /= r;
		c[1] /= r;
		c[2] /= r;
	}
	return true;
}

template <class T> bool crossProd( const vector<T>& p1, const vector<T>& p2,
		vector<T>& c, bool normalize=false ){
	if ( p1.size() != 3 ) {
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<p1.size()<<endl;
		return false;
	}
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	if ( normalize ) {
		double r = sqrt( c[0] * c[0] + c[1] * c[1] + c[2] * c[2] ) ;
		c[0] /= r;
		c[1] /= r;
		c[2] /= r;
	}
	return true;
}

template <class T> bool crossDotProd( const vector<T>& p1, const vector<T>& p2,
		vector<T>& c, T& cosT, bool normalize ) {
	if ( p1.size() != 3 ) {
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<p1.size()<<endl;
		return false;
	}
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	if ( normalize ) {
		double p1Norm = p1[0] * p1[0] +  p1[1] * p1[1] + p1[2] * p1[2] ;
		double p12Norm = p1Norm * ( p2[0] * p2[0] +  p2[1] * p2[1] + p2[2] * p2[2] );
		p12Norm = sqrt ( p12Norm );
		cosT = ( p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2] ) / p12Norm;
		double r = sqrt( c[0] * c[0] + c[1] * c[1] + c[2] * c[2] ) ;
		c[0] /= r;
		c[1] /= r;
		c[2] /= r;
	}
	return true;
}


inline void mTimesV ( const vector<double>& m, const vector<float>& v,
		vector<float>&u ) {
	const size_t dim = v.size();
	for ( unsigned int i = 0; i < dim; ++i ) {
		u[i] = 0.0;
		for (unsigned int j = 0; j < dim; ++j)
			u[i] += m [dim * i + j ] * v[j];
	}
}


inline bool interAngleF( const vector<float> &v1, const vector<float> & v2,
		const  size_t dim, float& cosTheta){

	double v1Len = 0.0, v2Len = 0.0;
	for (unsigned int i = 0; i < dim; i++) {
		v1Len += v1[i] * v1[i];
		v2Len += v2[i] * v2[i];
	}
	if ( v1Len < Tiny || v2Len < Tiny){
		cout<<"divide by ZERO in interAngle \n";
		return false;
	}
	cosTheta = 0.0;
	for (unsigned int i = 0; i < dim; i++)
		cosTheta += v1[i] * v2[i];
	cosTheta /= sqrt (v1Len * v2Len) ;

	return true;
}

inline bool interAngleD( const vector<double> &v1, const vector<double> & v2,
		const  size_t dim, double& cosTheta){

	double v1Len = 0.0, v2Len = 0.0;
	for (unsigned int i = 0; i < dim; i++) {
		v1Len += v1[i] * v1[i];
		v2Len += v2[i] * v2[i];
	}
	if ( v1Len < Tiny || v2Len < Tiny){
		cout<<"divide by ZERO in interAngle \n";
		return false;
	}
	cosTheta = 0.0;
	for (unsigned int i = 0; i < dim; i++)
		cosTheta += v1[i] * v2[i];
	cosTheta /= sqrt (v1Len * v2Len) ;

	return true;
}

inline bool interAngleD( const double* const v1, const double * const v2,
		const  size_t dim, double& angle){

	double v1Len = 0.0, v2Len = 0.0;
	for ( unsigned int i = 0; i < dim; i++){
		v1Len += v1[i] * v1[i];
		v2Len += v2[i] * v2[i];
	}
	if ( v1Len < Tiny || v2Len < Tiny){
		cout<<"divide by ZERO in interAngle \n";
		return false;
	}
	double cosValue = 0.0;

	for (unsigned int i = 0; i < dim; i++)
		cosValue += v1[i] * v2[i];

	angle = cosValue / sqrt (v1Len * v2Len) ;

//	printf("Angle:  %12.8f%12.8f%12.8f\%12.8f\n", angle, cosValue, v1Len, v2Len);

	return true;
}

inline bool interAngleF( const float* const v1, const float * const v2,
		const  size_t dim, float & angle){

	double v1Len = 0.0, v2Len = 0.0;
	for ( unsigned int i = 0; i < dim; i++){
		v1Len += v1[i] * v1[i];
		v2Len += v2[i] * v2[i];
	}
	if ( v1Len < Tiny || v2Len < Tiny){
		cout<<"divide by ZERO in interAngle \n";
		return false;
	}
	double cosValue = 0.0;

	for (unsigned int i = 0; i < dim; i++)
		cosValue += v1[i] * v2[i];

	angle = cosValue / sqrt (v1Len * v2Len) ;

//	printf("Angle:  %12.8f%12.8f%12.8f\%12.8f\n", angle, cosValue, v1Len, v2Len);

	return true;
}
/**
 * The angle computed may differ from the others by 180-$\theta$
 */
inline bool interAngleF(const vector<float>& p1, const vector<float>& p2,
		const vector<float>& p3, float& cosTheta ) {
	const size_t dim3 = 3;
	vector<float> v1(dim3, 0.0), v2(dim3, 0.0);
	for (unsigned int i = 0; i < dim3; ++i ) {
		v1[i] = p3[i] - p2[i];
		v2[i] = p1[i] - p2[i];
	}
	return interAngleF( v1, v2, dim3, cosTheta );
}

//
//inline bool planeAngleP2Mid(double* p1, double* p2, double* p3, double *angle){
//	const size_t dim3 = 3;
//	vector<float> r21(dim3, 0.0), r23(dim3, 0.0);
//	double lenR21 = 0.0,  lenR23 = 0.0;
//	for (unsigned int i = 0; i < dim3; ++i) {
//		r21[i] = p1[i] - p2[i];
//		r23[i] = p3[i] - p2[i];
//		lenR21 += r21[i]  *  r21[i] ;
//		lenR23 += r23[i]  *  r23[i] ;
//	}
//	if( lenR21 < eps1 || lenR23 < eps1 ) {
//		return false;
//	}
//	if( colinear(p1, p2, p3, dim3) ){
//		if( dotProd(r21, r23, dim3) > 0.0 )
//			*angle = 0.1*eps1;
//		else
//			*angle =  180.0;
//		return true;
//	}
//	double length = sqrt ( lenR21 * lenR23 );
//	double cosAngle = dotProd(r21, r23, dim3) / length;
//	angle =RadToDegrees * (acos(cosAngle));
//	return true;
//}

inline bool planeAngleP2Mid( const vector<float>& p1, const vector<float>& p2,
		const vector<float> & p3, float& angle){
	const size_t dim3 = 3;
	vector<float> r21(dim3, 0.0), r23(dim3, 0.0);
	double lenR21 = 0.0,  lenR23 = 0.0;
	for (unsigned int i = 0; i < dim3; ++i) {
		r21[i] = p1[i] - p2[i];
		r23[i] = p3[i] - p2[i];
		lenR21 += r21[i]  *  r21[i] ;
		lenR23 += r23[i]  *  r23[i] ;
	}
	if( lenR21 < eps1 || lenR23 < eps1 ) {
		return false;
	}
//	if( colinear(p1, p2, p3) ){
//		if( dotProd(r21, r23, dim3) > 0.0 )
//			angle= 0.1*eps1;
//		else
//			angle =  180.0;
//		return true;
//	}
	double length = sqrt ( lenR21 * lenR23 );
	double cosAngle = dotProd(r21, r23, dim3) / length;
	angle =RadToDegrees * (acos(cosAngle));
	return true;
}

inline bool planeAngleP2Mid( const vector<double>& p1, const vector<double>& p2,
		const vector<double> & p3, double& angle ) {

	const size_t dim3 = 3;
	vector<float> r21(dim3, 0.0), r23(dim3, 0.0);
	double lenR21 = 0.0,  lenR23 = 0.0;
	for (unsigned int i = 0; i < dim3; ++i) {
		r21[i] = p1[i] - p2[i];
		r23[i] = p3[i] - p2[i];
		lenR21 += r21[i]  *  r21[i] ;
		lenR23 += r23[i]  *  r23[i] ;
	}
	if( lenR21 < eps1 || lenR23 < eps1 ) {
		return false;
	}
//	if( colinear(p1, p2, p3, dim3) ){
//		if( dotProd(r21, r23, dim3) > 0.0 )
//			angle= 0.1*eps1;
//		else
//			angle =  180.0;
//		return true;
//	}
	double length = sqrt ( lenR21 * lenR23 );
	double cosAngle = dotProd(r21, r23, dim3) / length;
	angle =RadToDegrees * (acos(cosAngle));
	return true;
}

/**
 *  compute dihedral angle from three vectors,
 */
inline double signedDihedralAngle( double *const b1, double *const b2, double *const b3,const size_t d ){
	double b12c [d];
	double b23c [d];
	crossProd(b2, b3, b23c, d);
	crossProd(b1, b2, b12c, d);
	double x = dotProd(b12c, b23c, d);
	double y = dotProd(b1, b23c, d);
	y *= norm(b2, d);
	double phi = -atan2(y, x);
//	if (phi > Pi)
//		phi -= TwoPi;
//	else
//		if (phi < -Pi)
//			phi += TwoPi;
	return phi;
}


/**
 * compute dihedral angle from four points
 * */
inline bool signedDihedralAngle( double *const p1, double *const p2, double *const p3,
		double *const p4, const size_t d, double& angle ){
	if( colinear (p1, p2, p3, d) ){
		cerr<<"p1, p2 and p3 colinear "<<endl;
		return false;
	}else if (colinear(p2,p3,p4,d) ){
		cerr<<"p2, p3, p4 colinear"<<endl;
		return false;
	}
	double *b1 = internucleiVec(p2, p1, d);
	double *b2 = internucleiVec(p3, p2, d);
	double *b3 = internucleiVec(p4, p3, d);
	double b12c [d];
	double b23c [d];
	crossProd(b2, b3, b23c, d);
	crossProd(b1, b2, b12c, d);
	double x = dotProd(b12c, b23c, d);
	double y = dotProd(b1, b23c, d);
	y *= norm(b2, d);
	angle  = -atan2(y, x);
	delete [] b1;
	delete [] b2;
	delete [] b3;
	return true;
}

/**
 * Compute the dihedral angle formed by four points $p1, p2, p3, p4$
 * or between two planes formed by $p1, p2, p3$ and $p2, p3, p4$ with
 * $p2--p3$ as the shared edge
 *	                                     p4
 *	                                     /
 *									   /  b3
 *						 b2		/
 *	            p2------------- p3
 *			   /
 *			 /  b1
 *		   /
 *		 p1
 *		 SEE http://en.wikipedia.org/wiki/Dihedral_angle
 *	$\varphi = \operatorname{atan2} \left( \left([\mathbf{b}_1 \times \mathbf{b}_2]\times [\mathbf{b}_2
 *	\times \mathbf{b}_3]\right) \cdot \frac{\mathbf{b}_2}{|\mathbf{b}_2|},
 *	 [\mathbf{b}_1 \times \mathbf{b}_2] \cdot [\mathbf{b}_2 \times \mathbf{b}_3] \right)$
 */
inline bool signedDihedralAngle( const vector<float>& p1,  const vector<float>& p2,
		const vector<float>& p3,  const vector<float>& p4, float& angle ){
	const size_t d = 3;
	if( colinear (p1, p2, p3) ){
		cerr<<"p1, p2 and p3 colinear "<<endl;
		return false;
	}else if (colinear(p2, p3, p4) ){
		cerr<<"p2, p3, p4 colinear"<<endl;
		return false;
	}
	vector<float> b1 = internucleiVec(p2, p1, d);
	vector<float> b2 = internucleiVec(p3, p2, d);
	vector<float> b3 = internucleiVec(p4, p3, d);
	vector<float>b12c(d, 0.0);
	vector<float>b23c(d, 0.0);
	crossProd(b2, b3, b23c );  // the norm of the plane composed of $p2, p3, p4$
	crossProd(b1, b2, b12c );   // the norm of the plane composed of $p1, p2, p3$

	float x = dotProd(b12c, b23c, d);
	float y = dotProd(b1, b23c, d);
	y *= norm(b2, d);
	angle  = - RadToDegrees * atan2(y, x);
	return true;
}


/**
 * compute the polar coordinate for a 3D point
 @param:  polarCoords return, respectively,
				 $radius,$\cos{\theta}, \sin{\theta}, cos{\phi}, sin{\phi} $
 */
inline int toPolarTrigonoF( double *const coords, const size_t dim, double *polarCoords ){
	if ( dim != 3){
		cerr<<"The function only computes polar coordinates for 3D,  but the dimension d ="<<dim<<endl;
		exit(0);
	}
	double dis = 0.0;
	for (unsigned int i = 0; i < dim; i++ )
		dis += coords[i]  * coords[i] ;

	const double r = polarCoords[0] = sqrt( dis );   // the radius
	if ( polarCoords[0] < Tiny  )
		return 0;
	polarCoords[1] = coords[2] / r;  //  $ \cos{\theta} $
	polarCoords[2] = sqrt( 1.0 - coords[2] * coords[2] / dis ); //  $ \sin{\theta} $
	if ( polarCoords[2]  < Tiny  ){
		polarCoords[1]  = 1.0;
		return 1;
	}
	const double RsinTheta = r * polarCoords[2];
	polarCoords[3] =coords[0] /  RsinTheta ; // $ \cos{\phi} $
	polarCoords[4] =coords[1] / RsinTheta ; // $ \sin{\phi} $
	return 2;
}

/**
 * compute the polar coordinate for a 3D point
 @param:  polarCoords return, respectively,
				 $radius,$\cos{\theta}, \sin{\theta}, cos{\phi}, sin{\phi} $
 */
inline int toPolarTrigonoF(const vector<double>&coords, const size_t dim, double *polarCoords ){
	if ( dim != 3){
		cerr<<"The function only computes polar coordinates for 3D,  but the dimension d ="<<dim<<endl;
		exit(0);
	}
	double dis = 0.0;
	for (unsigned int i = 0; i < dim; i++ )
		dis += coords[i]  * coords[i] ;

	const double r = polarCoords[0] = sqrt( dis );   // the radius
	if ( polarCoords[0] < Tiny  )
		return 0;
	polarCoords[1] = coords[2] / r;  //  $ \cos{\theta} = z / r $
	//  $ \sin{\theta} = \sqrt( 1.0 -  \cos{\theta}^2) $
	polarCoords[2] = sqrt( 1.0 - coords[2] * coords[2] / dis );
	if ( fabs( polarCoords[2] ) < Tiny  ){
		polarCoords[1]  = 1.0;
		return 1;
	}
	const double RsinTheta = r * polarCoords[2]; //  $ r \sin{\theta} $
	polarCoords[3] =coords[0] /  RsinTheta ; // $ \cos{\phi} $
	polarCoords[4] =coords[1] / RsinTheta ; // $ \sin{\phi} $
	return 2;
}

/**
 * compute the polar coordinate for an 3D point
 @param:  polarCoords for returning the $radius, $cos(\theta)$,  $sin( \phi )$ and  $cos( \phi )$
 */
inline int toPolarCoord(double *const coords, const size_t dim, double *polarCoords ){
	if ( dim != 3){
		cerr<<"The function only computes polar coordinates for 3D,  but the dimension d ="<<dim<<endl;
		exit(0);
	}
	double dis = 0.0;
	for (unsigned int i = 0; i < dim; i++ )
		dis += coords[i]  * coords[i] ;
	const double r = polarCoords[0] = sqrt( dis );   // the radius
	if ( polarCoords[0] < Tiny  )
		return 0;
	polarCoords[1] = coords[2] / r;  //  $ \cos{\theta} $
	if ( ( 1.0 - coords[2] * coords[2] / dis ) < 0.0 ) {
		polarCoords[1]  = 1.0;
		return 1;
	}
	polarCoords[2] = sqrt( 1.0 - coords[2] * coords[2] / dis ); // $ \sin{\theta} $
	if ( polarCoords[2]  < Tiny  ){
		polarCoords[1]  = 1.0;
		return 1;
	}
	polarCoords[3] = coords[0] / r / polarCoords[2]  ; // $ \cos{\phi} $
	polarCoords[4] = coords[1] / r / polarCoords[2]  ; // $ \sin{\phi} $
	// atan2 return [ -pi, +pi ], add $\pi$ if want to change to the range between $[0, 2\pi]$
	polarCoords[5] = atan2( coords[1], coords[0] );
//	if ( polarCoords[5] < 0.0 )
//		polarCoords[5] += TwoPi;
	return 2;
}

/**
 * compute the polar coordinate for an 3D point[
 * 		polarCoords[0] = r;
 * 		polarCoords[1] = $\theta$;
 * 		polarCoords[2] = $\phi$;
 @param:  polarCoords for returning the $radius, $\theta$, and  $\phi$
 The returned $\theta$ is always in the range of [ 0.0, $\pi$ ]
 How about the $\phi$, it must be returned as [0.0, $2\pi$ ]
 The issue of discontinuous point at the 0 and $2\pi$ could become a problem
 */
inline int toPolarAngle( double *const coords, const size_t dim, double *polarCoords ){
	if ( dim != 3){
		cerr<<"Can only computes polar coordinates for 3D,  but the dimension d ="<<dim<<endl;
		exit(0);
	}
	double dis = 0.0;
	for (unsigned int i = 0; i < dim; i++ )
		dis += coords[i]  * coords[i] ;

	const double r = sqrt( dis );
	if ( r < Tiny  )
		return 0;
	polarCoords[0] = r;
	const double sinT = fabs( 1.0 - coords[2] * coords[2] / dis);
	if ( sinT < Tiny  ){
		polarCoords[1]  = 0.0;
		return 1;
	}
	polarCoords[1]  = acos( coords[2] / r ); // $theta$
	// atan2 return [ -pi, +pi ], add $\pi$ if want to change to the range between $[0, 2\pi]$
	polarCoords[2] = atan2( coords[1], coords[0] );
//	if ( polarCoords[2] < 0.0 )
//		polarCoords[2] += TwoPi;
	return 2;
}

/**
 * compute the intersection points if the segment does intersect with the sphere
 * @param edge:        four elements representing two ends
 * @param dimID:      0--x, 1--y and 2--z
 * @param center:     the 3D coordinate for the center of the sphere
 * @param radius:     the radius of the sphere
 */
inline size_t gridSphereIntersection( double *const edge, const size_t dimID, double *const center,
		const double radius, double *s0, double *s1 ){

	double low, up, x, y, z;
	double a = 1.0, b = 0.0, c = 0.0; // the coefficients

	if (dimID == 0 ){ // along X-axis
		low = edge[0];
		up  = edge[1];
		y    = edge[2];
		z    = edge[3];
		b = - 2.0 * center[0];
		c = center[0] * center[0] + ( y - center[1] ) * ( y - center[1] ) +  ( z - center[2] ) * ( z - center[2] )
						 - radius * radius;
	}else if  (dimID == 1 ){
		x = edge[0];
		low = edge[1];
		up  = edge[2];
		z = edge[3];
		b = - 2.0 * center[1];
		c = center[1] * center[1] + ( x - center[0] ) * ( x - center[0] ) +  ( z - center[2] ) * ( z - center[2] )
									- radius * radius;

	}else if  (dimID == 2 ){
		x = edge[0];
		y = edge[1];
		low = edge[2];
		up  = edge[3];
		b = - 2.0 * center[2];
		c = center[2] * center[2] + ( x - center[0] ) * ( x - center[0] ) +  ( y - center[1] ) * ( y - center[1] )
								   - radius * radius;
	}

	size_t noOfSolution = gsl_poly_solve_quadratic(a, b, c, s0, s1);
	if ( noOfSolution == 0 )
		return 0;
	if (  *s0 >= low && *s0 < up  && *s1 >= low && *s1 < up )
		return 2;
	if (  *s0 >= low && *s0 < up )
		return 1;
	if (  *s1 >= low && *s1 < up ){
		*s0 = *s1;
		return 1;
	}
	return 0;
}

/** */
inline bool segmentSphereIntersection( double *const p1, double *const p2, double *const center,
		const double radius, double &s0, double &s1){

	double e0 = p2[0] - p1[0];
	double e1 = p2[1] - p1[1];
	double e2 = p2[2] - p1[2];
	double f0  = p1[0] - center[0];
	double f1  = p1[1] - center[1];
	double f2  = p1[2] - center[2];
	 // the coefficients
	double a = e0 * e0 + e1 * e1 + e2 * e2;
	double b = 2.0 * ( e0 * f0 + e1 * f1 + e2 * f2) ;
	double c = f0 * f0 + f1 * f1 + f2 * f2 - radius * radius;
	size_t noOfSolution = gsl_poly_solve_quadratic(a, b, c, &s0, &s1);
	if ( noOfSolution == 0 )
		return false;
	if (  ( s0 > 0.0 && s0 < 1.0 ) || ( s1  > 0.0 && s1 < 1.0 )  )
		return true;
	else
		return false;
}

/**
 * projection a 3D point $(x_0, y_0, z_0)$ into a plane represented by
 * $ax + by + cz +d = 0$ with a normal $(n_x, n_y, n_z) = (a, b, c)$
 * the perpendicular line through the point is
 * $\frac{x-x_0} {a} = \frac{y-y_0} {b} = \frac{z-z_0} {c}$ ;
 * so we can derive the following equation for $z$
 * $  (a^2 / c + b^2 / c + c) z = (a^2 / c + b^2 / c ) z_0 - ( a x_0 + b y_0 + d )$
 */
inline bool projection2Plane( const vector<double>& coef, const vector<float>& p,
		vector<float>& pj ) {
	double a = coef[0];
	double b = coef[1];
	double c = coef[2];
	double d = coef[3];
	double m = ( a * a + b * b);
	if(  m  < eps1 )  { // an XY plane
		pj[2] = 0.0;
		pj[1] = p[1];
		pj[0] = p[0] ;
		return true;
	}else if (  (a * a + c * c) < eps1  ) {  // an XZ plane
		pj[2] = p[2];
		pj[1] = 0.0;
		pj[0] = p[0] ;
		return true;
	}else if (  (b * b + c * c) < eps1 ) {  // an YZ plane
		pj[2] = p[2];
		pj[1] = p[1] ;
		pj[0] = 0.0;
	}else  if ( ( m + c * c ) < eps1 || fabs(c) < eps1 || fabs ( m / c + c ) < eps1 )
		return false;
	pj[2] = ( p[2] * m / c - a * p[0] - b * p[1] - d ) / ( m / c + c);
	pj[1] = p[1] + b * ( pj[2] - p[2]) / c ;
	pj[0] = p[0] + a * ( pj[2] - p[2]) / c ;
	return true;
}

inline bool coefGenerator (double *v1, double *v2, double v1x, double v2x, double *coef){

	double detVal = v1[0] * v2[1] - v1[1] * v2[0];
	if ( fabs(detVal) < eps22 )
		return false;
	double cx = ( v1[1] * v2[2] - v2[1] * v1[2]) / detVal;
	double dx = ( v1x * v2[1] - v2x * v1[1]) / detVal;
	double cy = ( v2[0] * v1[2] - v1[0] * v2[2]) / detVal;
	double dy = ( v1[0] * v2x - v2[0] * v1x) / detVal;
	coef[0] = cx * cx + cy * cy + 1.0;
	coef[1]  = 2.0 * ( cx * dx + cy * dy );
	coef[2]  = dx* dx + dy * dy - 1.0;
	coef[3] = cx;
	coef[4] = dx;
	coef[5] = cy;
	coef[6] = dy;
	return true;
}


inline bool coefGenerator (const vector<double>& v1, const vector<double>& v2,
		double v1x, double v2x, double *coef){

	double detVal = v1[0] * v2[1] - v1[1] * v2[0];
	if ( fabs(detVal) < eps1 )
		return false;
	double cx = ( v1[1] * v2[2] - v2[1] * v1[2]) / detVal;
	double dx = ( v1x * v2[1] - v2x * v1[1]) / detVal;
	double cy = ( v2[0] * v1[2] - v1[0] * v2[2]) / detVal;
	double dy = ( v1[0] * v2x - v2[0] * v1x) / detVal;
	coef[0] = cx * cx + cy * cy + 1.0;
	coef[1]  = 2.0 * ( cx * dx + cy * dy );
	coef[2]  = dx* dx + dy * dy - 1.0;
	coef[3] = cx;
	coef[4] = dx;
	coef[5] = cy;
	coef[6] = dy;
	return true;
}

/**
 * We only return a single matrix, it seems that the returned matrix is a left rotation
 * under certain conditions
 * @param mat: must be a 9 elements matrix,
 * @param dim: must be 3
 */

inline bool rotmat(double *v1, double *v2, double *v1N, double *v2N, const size_t dim, double *mat) {
	double coef[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	if ( ! coefGenerator(v1, v2, v1N[0], v2N[0], coef ) ) {
		cout<<"returned 01\n";
		return false;
	}
	double x0;
	double x1;
	double discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	int noOfSolution = 0;
	if (fabs (discriminant) < eps5 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a1 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] /coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 01a\n";
			return false;
		}
	}
	double mm1[6] = { coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1 };

	if ( ! coefGenerator(v1, v2, v1N[1], v2N[1], coef ) ) {
		cout<<"returned 02\n";
		return false;
	}
	//	cout<<coef[0]<<" - "<<coef[1]<<"  "<< coef[2]<<endl;
	//	printf("%3.12f%3.12f%3.12f\n", coef[0], coef[1], coef[2]);
	//	double d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	//	cout<<"d1 = "<<d<<endl;
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps4 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a2 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 022a\n";
			cout<<coef[0]<<" "<<coef[1]<<"  "<< coef[2]<<endl;
			//			d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
			//			cout<<"d2 = "<<d<<endl;
			return false;
		}
	}
	double mm2[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};


	if ( ! coefGenerator(v1, v2, v1N[2], v2N[2], coef ) ) {
		cout<<"returned 03\n";
		return false;
	}
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps4 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a3 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 03a\n";
			return false;
		}
	}
	double mm3[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};

	// double mat[ dim * dim ] ; // the 9 elements of the rotation matrix
	double v1b[dim];
	double v2b[dim];

	for (int i = 0; i < 2; i++ ){
		mat[0] = mm1[ dim * i ];
		mat[1] = mm1[ dim * i + 1];
		mat[2] = mm1[ dim * i + 2];
		for (int j = 0; j < 2; j++ ){
			mat[3] = mm2[ dim * j ];
			mat[4] = mm2[ dim * j + 1];
			mat[5] = mm2[ dim * j + 2];
			for (int k = 0; k < 2; k++){
				mat[6] = mm3[ dim * k ];
				mat[7] = mm3[ dim * k + 1];
				mat[8] = mm3[ dim * k + 2];
			}

			for ( unsigned int p = 0; p < dim; p++){
				v1b [ p ] = 0.0;
				v2b [ p ] = 0.0;
				for ( unsigned int q = 0; q < dim; q++) {
					v1b [ p ] += mat [ p * dim + q ] * v1[ q];
					v2b [ p ] += mat [ p * dim + q ] * v2[ q ];
				}
			}
			if (
					fabs(v1b[0] - v1N[0] ) < eps22 &&
					fabs(v1b[1] - v1N[1] ) < eps22 &&
					fabs(v1b[2] - v1N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs(v2b[2] - v2N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs ( det(mat, dim) - 1.0) < eps4
			) {
				//				matrix cbMat(mat, 3, 9);
				//				vector<double> v1bb = cbMat.transpose().times(v1N, 3);
				//				vector<double> v2bb = cbMat.transpose().times(v2N, 3);
				//				cout<<"inside:\n";
				//				for(unsigned int m = 0; m < dim; m++)
				//					cout<<v1bb[m]<<" ";
				//				cout<<endl;
				//				for(unsigned int m = 0; m < dim; m++)
				//					cout<<v2bb[m]<<" ";
				//				cout<<endl;
				//				cout<<endl;
				return true;
			}
			//			cout<<"\n********************************\n";
			//			cout<<det(mat, dim)<<endl;signedDihedralAngle
			//			cout<<v1b[0]<<"  "<<v1b[1]<<"  "<<v1b[2]<<endl;
			//			cout<<v2b[0]<<"  "<<v2b[1]<<"  "<<v2b[2]<<endl;
			//			cout<<"\n********************************\n";
		}
	}
	cout<<"can not compute back \n";
	return false;
}


/**
 * We only return a single matrix, it seems that the returned matrix is a left rotation
 * under certain conditions
 * the final yz plane is defined as follows:
 @param v1N  the rotated vector of v1 ( along +Y)
 @param v2N  the roated vector of v2  ( with $\theta$ angle from +Z )
 *

 v2N (C')     +Z
 \                 |
    \              |
       \           |
          \        |
             \     | $\theta$
                   N ------------------------------+Y ( v1N ) (CA)
 * @param mat: must be a 9 elements array (a 3x3 matrix),
 * @param dim: must be 3
 */

inline bool rotmat(const vector<double>& v1, const vector<double>& v2, const vector<double>& v1N,
		const vector<double>&v2N, const size_t dim, double *mat){

	double coef[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	if ( ! coefGenerator(v1, v2, v1N[0], v2N[0], coef ) ) {
		cout<<"returned 01\n";
		return false;
	}
	double x0;
	double x1;
	double discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	int noOfSolution = 0;
	if (fabs (discriminant) < eps5 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a1 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] /coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 01a\n";
			return false;
		}
	}
	double mm1[6] = { coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1 };

	if ( ! coefGenerator(v1, v2, v1N[1], v2N[1], coef ) ) {
		cout<<"returned 02\n";
		return false;
	}
//	cout<<coef[0]<<" - "<<coef[1]<<"  "<< coef[2]<<endl;
//	printf("%3.12f%3.12f%3.12f\n", coef[0], coef[1], coef[2]);
//	double d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
//	cout<<"d1 = "<<d<<endl;
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps4 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a2 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 022a\n";
			cout<<coef[0]<<" "<<coef[1]<<"  "<< coef[2]<<endl;
//			d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
//			cout<<"d2 = "<<d<<endl;
			return false;
		}
	}
	double mm2[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};

	if ( ! coefGenerator(v1, v2, v1N[2], v2N[2], coef ) ) {
		cout<<"returned 03\n";
		return false;
	}
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps4 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a3 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 03a\n";
			return false;
		}
	}
	double mm3[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};

	// double mat[ dim * dim ] ; // the 9 elements of the rotation matrix
	double v1b[dim];
	double v2b[dim];

	for (int i = 0; i < 2; i++ ){
		mat[0] = mm1[ dim * i ];
		mat[1] = mm1[ dim * i + 1];
		mat[2] = mm1[ dim * i + 2];
		for (int j = 0; j < 2; j++ ){
			mat[3] = mm2[ dim * j ];
			mat[4] = mm2[ dim * j + 1];
			mat[5] = mm2[ dim * j + 2];
			for (int k = 0; k < 2; k++){
				mat[6] = mm3[ dim * k ];
				mat[7] = mm3[ dim * k + 1];
				mat[8] = mm3[ dim * k + 2];
			}

			for ( unsigned int p = 0; p < dim; p++){
				v1b [ p ] = 0.0;
				v2b [ p ] = 0.0;
				for ( unsigned int q = 0; q < dim; q++) {
					v1b [ p ] += mat [ p * dim + q ] * v1[ q];
					v2b [ p ] += mat [ p * dim + q ] * v2[ q ];
				}
			}
			if (
					fabs(v1b[0] - v1N[0] ) < eps22 &&
					fabs(v1b[1] - v1N[1] ) < eps22 &&
					fabs(v1b[2] - v1N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs(v2b[2] - v2N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs ( det(mat, dim) - 1.0) < eps4
			) {
//				matrix cbMat(mat, 3, 9);
//				vector<double> v1bb = cbMat.transpose().times(v1N, 3);
//				vector<double> v2bb = cbMat.transpose().times(v2N, 3);
//				cout<<"inside:\n";
//				for(unsigned int m = 0; m < dim; m++)
//					cout<<v1bb[m]<<" ";
//				cout<<endl;
//				for(unsigned int m = 0; m < dim; m++)
//					cout<<v2bb[m]<<" ";
//				cout<<endl;
//				cout<<endl;
				return true;
			}
			//			cout<<"\n********************************\n";
			//			cout<<det(mat, dim)<<endl;signedDihedralAngle
			//			cout<<v1b[0]<<"  "<<v1b[1]<<"  "<<v1b[2]<<endl;
			//			cout<<v2b[0]<<"  "<<v2b[1]<<"  "<<v2b[2]<<endl;
			//			cout<<"\n********************************\n";
		}
	}
	cout<<"can not compute back \n";
	return false;
}

/**
 * This routine has a critical problem when the coef for the quadratic equation is very large that
 * results a relative large discriminant, and cause the equation unsolvable. What is the solution ?
 */
inline bool rotmat( const vector<double>& v1, const vector<double>& v2,
		const vector<double>& v1N, const vector<double>&v2N,
		const size_t dim, vector<double> &mat ) {

	double coef[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	if ( ! coefGenerator(v1, v2, v1N[0], v2N[0], coef ) ) {
		cout<<"returned 01\n";
		return false;
	}
	double x0;
	double x1;
	double discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	int noOfSolution = 0;
	if (fabs (discriminant) < eps5 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a1 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] /coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 01a\n";
			return false;
		}
	}
	double mm1[6] = { coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1 };

	if ( ! coefGenerator(v1, v2, v1N[1], v2N[1], coef ) ) {
		cout<<"returned 02\n";
		return false;
	}
//	cout<<coef[0]<<" - "<<coef[1]<<"  "<< coef[2]<<endl;
//	printf("%3.12f%3.12f%3.12f\n", coef[0], coef[1], coef[2]);
//	double d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
//	cout<<"d1 = "<<d<<endl;
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps5 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a2 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 022a\n";
			cout<<coef[0]<<" "<<coef[1]<<"  "<< coef[2]<<endl;
			double d = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
			cout<<"d2 = "<<d<<endl;
			return false;
		}
	}
	double mm2[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};


	if ( ! coefGenerator(v1, v2, v1N[2], v2N[2], coef ) ) {
		cout<<"returned 03\n";
		return false;
	}
	discriminant = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
	if (fabs (discriminant) < eps5 ) {
		if ( fabs ( coef[0]) < eps1 ) {
			cout<<"a3 = "<<coef[0]<<endl;
			return false;
		}
		x0 = x1 = -0.50 * coef[1] / coef[0];
	}else {
		noOfSolution = gsl_poly_solve_quadratic(coef[0], coef[1], coef[2], &x0, &x1);
		if ( noOfSolution == 0 ) {
			cout<<"returned 03a\n";
			return false;
		}
	}
	double mm3[6] = {coef[3] * x0 + coef[4], coef[5] * x0 + coef[6], x0,
			coef[3] * x1 + coef[4], coef[5] * x1 + coef[6], x1};

	// double mat[ dim * dim ] ; // the 9 elements of the rotation matrix
	double v1b[dim];
	double v2b[dim];

	for (int i = 0; i < 2; i++ ){
		mat[0] = mm1[ dim * i ];
		mat[1] = mm1[ dim * i + 1];
		mat[2] = mm1[ dim * i + 2];
		for (int j = 0; j < 2; j++ ){
			mat[3] = mm2[ dim * j ];
			mat[4] = mm2[ dim * j + 1];
			mat[5] = mm2[ dim * j + 2];
			for (int k = 0; k < 2; k++){
				mat[6] = mm3[ dim * k ];
				mat[7] = mm3[ dim * k + 1];
				mat[8] = mm3[ dim * k + 2];
			}

			for ( unsigned int p = 0; p < dim; p++){
				v1b [ p ] = 0.0;
				v2b [ p ] = 0.0;
				for ( unsigned int q = 0; q < dim; q++) {
					v1b [ p ] += mat [ p * dim + q ] * v1[ q];
					v2b [ p ] += mat [ p * dim + q ] * v2[ q ];
				}
			}
			if (
					fabs(v1b[0] - v1N[0] ) < eps22 &&
					fabs(v1b[1] - v1N[1] ) < eps22 &&
					fabs(v1b[2] - v1N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs(v2b[2] - v2N[2] ) < eps22 &&
					fabs(v2b[0] - v2N[0] ) < eps22 &&
					fabs ( det(mat) - 1.0) < eps4
			) {
//				matrix cbMat(mat, 3, 9);
//				vector<double> v1bb = cbMat.transpose().times(v1N, 3);
//				vector<double> v2bb = cbMat.transpose().times(v2N, 3);
//				cout<<"inside:\n";
//				for(unsigned int m = 0; m < dim; m++)
//					cout<<v1bb[m]<<" ";
//				cout<<endl;
//				for(unsigned int m = 0; m < dim; m++)
//					cout<<v2bb[m]<<" ";
//				cout<<endl;
//				cout<<endl;
				return true;
			}
			//			cout<<"\n********************************\n";
			//			cout<<det(mat, dim)<<endl;
			//			cout<<v1b[0]<<"  "<<v1b[1]<<"  "<<v1b[2]<<endl;
			//			cout<<v2b[0]<<"  "<<v2b[1]<<"  "<<v2b[2]<<endl;
			//			cout<<"\n********************************\n";
		}
	}
	cout<<"can not compute back \n";
	return false;
}

inline void mTimesV ( const vector<double>& m, const double* const v,
		vector<double>&u ) {
	const size_t dim = u.size();
	for ( unsigned int i = 0; i < dim; ++i ) {
		u[i] = 0.0;
		for (unsigned int j = 0; j < dim; ++j)
			u[i] += m [dim * i + j ] * v[j];
	}
}

inline void mTimesV ( const vector<double>& m, const vector<double>& v,
		vector<double>&u ) {
	const size_t dim = v.size();
	for ( unsigned int i = 0; i < dim; ++i ) {
		u[i] = 0.0;
		for (unsigned int j = 0; j < dim; ++j)
			u[i] += m [dim * i + j ] * v[j];
	}
}

inline void mTimesV ( const double *const  m, const vector<double>& v,
		vector<double>&u ) {
	const size_t dim = v.size();
	for ( unsigned int i = 0; i < dim; ++i ) {
		u[i] = 0.0;
		for (unsigned int j = 0; j < dim; ++j)
			u[i] += m [dim * i + j ] * v[j];
	}
}


inline void mTimesV ( const size_t dim, const double *const  m, const double *const v,
		double u[] ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		u[i] = 0.0;
		for (unsigned int j = 0; j < dim; ++j)
			u[i] += m [dim * i + j ] * v[j];
	}
}

inline void mTimesN ( const size_t dim, const double *const  m, const double *const n,
		double w[] ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

inline void mTimesN ( const size_t dim, const double *const  m, const double *const n,
		vector<double>& w ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

inline void mTimesN ( const size_t dim, const double *const  m, const vector<double>& n,
		vector<double>& w ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

inline void mTimesN ( const size_t dim, const vector<double>& m, const double *const n,
		double w[] ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

inline void mTimesN ( const size_t dim, const vector<double>& m, const double *const n,
		vector<double>& w ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

inline void mTimesN ( const size_t dim, const vector<double>& m, const vector<double>& n,
		vector<double>& w ) {

	for ( unsigned int i = 0; i < dim; ++i ) {
		for (unsigned int j = 0; j < dim; ++j){
			w[ dim * i + j ] = 0.0;
			for (unsigned int k =0; k < dim; ++k)
				w[ dim * i + j ]  += m [dim * i + k ] * n[ k * dim + j ];
		}
	}
}

/**
 *
 * Rotation matrix given an axis and an angle.  make a rotation along the given axis.
 * Given a unit vector $u = (u_x, u_y, u_z)$, where $u_x^2 + u_y^2 + u_z^2 = 1$,
 * the matrix for a rotation by an angle of $\theta$ along $u$ is
 *     R = \begin{bmatrix}
 *      u_x^2+(1-u_x^2)c & u_x u_y(1-c)-u_zs & u_x u_z(1-c)+u_ys \\
 *      u_x u_y(1-c)+u_zs & u_y^2+(1-u_y^2)c & u_y u_z(1-c)-u_xs \\
 *      u_x u_z(1-c)-u_ys & u_y u_z(1-c)+u_xs & u_z^2+(1-u_z^2)c
 *      \end{bmatrix}
*/
inline void rotationAlongAxisWithAngle(const double s, const double c, double *const u,
		 double mat[] ){
	const double cc = 1.0 - c;
	mat[0] = c + cc *  u[0] * u [0] ;
	mat[1] = cc *  u[0] * u [1] - u[2] * s ;
	mat[2] = cc *  u[0] * u [2] + u[1] * s  ;
	mat[3] = cc *  u[0] * u [1] + u[2] * s ;
	mat[4] = c + cc *  u[1] * u [1] ;
	mat[5] = cc *  u[1] * u [2] - u[0] * s ;
	mat[6] =  cc *  u[0] * u [2] - u[1] * s ;
	mat[7] =  cc *  u[1] * u [2] + u[0] * s ;
	mat[8] = c + cc *  u[2] * u [2] ;
}

inline void rotationAlongAxisWithAngle(const double s, const double c, double *const u,
		 vector<double> &mat ){
	const double cc = 1.0 - c;
	mat[0] = c + cc *  u[0] * u [0] ;
	mat[1] = cc *  u[0] * u [1] - u[2] * s ;
	mat[2] = cc *  u[0] * u [2] + u[1] * s  ;
	mat[3] = cc *  u[0] * u [1] + u[2] * s ;
	mat[4] = c + cc *  u[1] * u [1] ;
	mat[5] = cc *  u[1] * u [2] - u[0] * s ;
	mat[6] =  cc *  u[0] * u [2] - u[1] * s ;
	mat[7] =  cc *  u[1] * u [2] + u[0] * s ;
	mat[8] = c + cc *  u[2] * u [2] ;
}

inline void rotationAlongAxisWithAngle(const double s, const double c,
		const vector<float>&u, vector<double> &mat ){
	const double cc = 1.0 - c;
	mat[0] = c + cc *  u[0] * u [0] ;
	mat[1] = cc *  u[0] * u [1] - u[2] * s ;
	mat[2] = cc *  u[0] * u [2] + u[1] * s  ;
	mat[3] = cc *  u[0] * u [1] + u[2] * s ;
	mat[4] = c + cc *  u[1] * u [1] ;
	mat[5] = cc *  u[1] * u [2] - u[0] * s ;
	mat[6] =  cc *  u[0] * u [2] - u[1] * s ;
	mat[7] =  cc *  u[1] * u [2] + u[0] * s ;
	mat[8] = c + cc *  u[2] * u [2] ;
}


/**
 @param m: the size of the atoms in the structure
 */
inline matrix rotMatrixByRMSD(const vector<float>& structCoordA,
		const vector<float>& structCoordB, vector<float>& centerA,
		vector<float> &centerB, bool& isSingularZero, double& rmsVal ){

	const size_t m = structCoordA.size() / dim3;
	if( m == 0 || structCoordA.size() != structCoordB.size() ) {
		cerr<<"Errors in RMSD cal, there in NO atoms in the array \n";
		return matrix(3,3);
	}
//	cout<<"noAtom4Overlap = "<<m<<endl;

	for (unsigned int j = 0; j< dim3; j++)
		centerA[j] = centerB[j] = 0.0;
	for ( unsigned int i = 0; i < m; i++){
		for ( unsigned int j = 0; j< dim3 ; j++) {
			centerA[ j ] += structCoordA[ dim3 * i + j];
			centerB[ j ] += structCoordB[ dim3 * i + j];
		}
	}

	for ( unsigned  int j = 0; j<dim3; j++) {
		centerA[j] /= m;
		centerB[j]  /= m;
	}
	const size_t size =  structCoordA.size();

	double arrA [size];
	double arrB [size];

	for (unsigned int i=0; i<m; i++){  //compute the Sum of the square of the distance between the center and all the points
		for (unsigned int j = 0; j<dim3; j++){
			arrA[ dim3*i+j ] = structCoordA[ dim3 * i + j] - centerA[j];
			arrB[ dim3*i+j ] = structCoordB[ dim3 * i + j] - centerB[j];
		}
	}

	double  arrU[ dim3 * dim3 ];
	for (unsigned int i = 0; i< dim3 ; i++ ){
		for  (unsigned int j = 0; j< dim3; j++) {
			arrU[ i * dim3 + j ] = 0.0;
			for (unsigned int k = 0; k< m; k++)
				arrU[ i * dim3 + j ]  += arrB [k * dim3 + j ] * arrA[ k * dim3 + i ];
		}
	}
	gsl_matrix_view U	= gsl_matrix_view_array (arrU, dim3, dim3);
	gsl_matrix *V = gsl_matrix_alloc (dim3, dim3);
	gsl_vector *sValues = gsl_vector_alloc ( dim3);
	gsl_vector *work = gsl_vector_alloc ( dim3 );

	gsl_linalg_SV_decomp (&U.matrix, V, sValues, work);
	for (unsigned int k = 0; k < dim3; k++){
		if ( fabs ( gsl_vector_get(sValues, k) ) < 1.0E-12) {
			isSingularZero = true;
			gsl_matrix_free (V);
			gsl_vector_free (work);
			gsl_vector_free (sValues);
			return matrix(3,3);
		}
	}

	//matrix rot = U.times( V.transpose() );
	const size_t dim9 = 9;
	vector<double> rotArray ( dim9, 0.0);
	for (unsigned int i = 0; i < dim3; i++ )
		for  (unsigned int j = 0; j < dim3; j++)
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3+j ] +=gsl_matrix_get(&U.matrix, i, k) * gsl_matrix_get (V, j, k);

	double rms = 0.0;
	double arrB2[ dim3 ];
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ]; //A[j][k] * B[k];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}

	double detValue = det(rotArray);
	if( fabs( detValue + 1.0) > 1.0E-8)  {
//		printf("rms0=%8.4f\n", rms);
		rmsVal = rms;
		isSingularZero = false;
		gsl_matrix_free (V);
		gsl_vector_free (work);
		gsl_vector_free (sValues);
		return matrix(rotArray, dim3, dim3);
	}

//	cout<<"left handedness !\n";
	double arrV[ dim9 ];
//	vector<double> rmsArr(dim3, 0.0);
	for (unsigned int i = 0; i < dim3; i++ )
		for  (unsigned int j = 0; j < dim3; j++)
			arrV[ dim3 * i + j ] = gsl_matrix_get (V, i, j);
	gsl_matrix_view Vv = gsl_matrix_view_array (arrV, dim3, dim3);

	 // reverse the x-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 0,  - gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 1,  gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
		gsl_matrix_set ( V, i, 2,  gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}
//	printf("rms1=%8.4f\n", rms);
	double minRMS = rms;
	size_t indexOfMinRms = 0;
	// reverse the y-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 1,  - gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 0,  gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
		gsl_matrix_set ( V, i, 2,  gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}
	if( rms < minRMS ) {
		minRMS = rms;
		indexOfMinRms = 1;
	}
//	printf("rms2=%8.4f\n", rms);

	// reverse the z-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 2,  - gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 0,  gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
		gsl_matrix_set ( V, i, 1,  gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray[ p * dim3 + q ] * arrB [ dim3 * i + q ]; //A[j][k] * B[k];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}
	if( rms < minRMS ) {
		minRMS = rms;
		indexOfMinRms = 2;
	}
//	printf("rms3=%8.4f\n", rms);

	rmsVal = minRMS;

	for (unsigned int i = 0; i< dim3; i++ ) {
		for (unsigned int j = 0; j < dim3; ++j )  {
			if ( j != indexOfMinRms)
				gsl_matrix_set ( V, i, j,  gsl_matrix_get (&Vv.matrix,  i, j) ) ;
			else
				gsl_matrix_set ( V, i, j,  - gsl_matrix_get (&Vv.matrix,  i, j) ) ;
		}
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}

	isSingularZero = false;
	gsl_matrix_free (V);
	gsl_vector_free (work);
	gsl_vector_free (sValues);
	return matrix(rotArray, dim3, dim3);
}

inline matrix fit2PlaneByRotation( const vector<float>&coords,
		const vector<double>& coef, bool& isSingularZero, double& rmsVal,
		bool debug = false ) {

	vector<float> coordVec(dim3, 0.0);
	vector<float>  pAtCoordVec( dim3, 0.0 );
	vector<float>  projCoordVec( coords.size(), 0.0 );
	bool canBeProjected = true;
	const size_t n = (size_t) coords.size() / dim3;
	double residual = 0.0;

	for ( unsigned int i = 0; i < n; ++i ){
		pAtCoordVec = vector<float>(dim3, 0.0);
		for (unsigned int k = 0; k < dim3; k++)
			coordVec[k] =  coords[ dim3 * i + k];
		if ( projection2Plane( coef, coordVec, pAtCoordVec) ) {
			for (unsigned int k = 0; k < dim3; k++) {
				projCoordVec[dim3 * i +k] = pAtCoordVec[k];
				residual +=(pAtCoordVec[k] - coordVec[k]) * ( pAtCoordVec[k] - coordVec[k]);
			}
		}else if (debug){
			canBeProjected = false;
			cout<<"could not be projected into a plane !\n";
			//			return false; //			break;
		}
	}
	if(debug) {
		residual = sqrt ( residual / n );
		printf("Fit to a plane, residual=%10.5f\n", residual);
	}

	vector<float> centerA(dim3, 0.0);
	vector<float> centerB(dim3, 0.0);
	matrix rot = rotMatrixByRMSD( projCoordVec, coords, centerA, centerB,
			isSingularZero, rmsVal );
	if( ! isSingularZero) {
		if(debug)
			printf("Rot to a plane, rms=%10.5f\n", sqrt(rmsVal/n));
		return rot;
	}
	if(debug)
		cout<<"singular !\n";
	return rotIdentity(dim3, dim3);
}

/**
 * See the document of LeastSquareFit.pdf for the detail
 * z = coef[0] x + coef[1] y + coef[2], the plane equation
 * the normal vector to a plane $ax + by + cz +d = 0 $ is (a, b, c) using the gradient
 * so in this case, it is (coef[0], coef[1] , -1)
 */
inline bool fit2Plane( double *const coords, double *coef, bool *onXYplane,
		const size_t n, double *center ){
	const size_t dim = 3;
	if( n < 3 ){
		cout<<"Must have more than three points in order to define a plane\n";
		return false;
	}
	double x = 0.0, y = 0.0, z = 0.0;
	double xy = 0.0, xz= 0.0, yz = 0.0;
	double xx = 0.0, yy= 0.0, zz = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x += coords[ dim * i ];
		y += coords[ i * dim + 1];
		z += coords[ i * dim + 2];
		xx += coords[ i * dim ] * coords[ i * dim ];
		yy += coords[ i * dim + 1 ] * coords[ i * dim +1];
		zz += coords[ i * dim + 2] * coords[ i * dim + 2];
		xy += coords[ i * dim ] * coords[ i * dim + 1];
		xz += coords[ i * dim ] * coords[ i * dim + 2];
		yz += coords[ i * dim + 1] * coords[ i *dim + 2];
	}
	// in the xy plane
	center[ 0 ] = x / n;
	center[ 1 ] = y / n;
	center[ 2 ] = z / n;
	if ( fabs(z) < eps1 ){
		*onXYplane = true;
//		coef[0] = 0.0;
//		coef[1] = 0.0;
//		coef[2] = 1.0;
		return true;
	}
	*onXYplane = false;

	double a_data[] = { xx, xy, x,
					                 xy, yy, y,
					                 x,   y,  (double)n};
	double b_data[] = { xz, yz, z };
	gsl_matrix_view m  = gsl_matrix_view_array (a_data, dim, dim);
	gsl_vector_view b = gsl_vector_view_array (b_data, dim);
	gsl_vector *solution = gsl_vector_alloc (dim);
	int s;
	gsl_permutation * p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );

//	printf ("coef = \n");
//	gsl_vector_fprintf (stdout, solution, "%g");
	for ( unsigned int k = 0; k < dim; k++)
		coef[k] = gsl_vector_get(solution, k);

	gsl_permutation_free (p);
	gsl_vector_free ( solution );

	return true;
}

/**
 * See the document of LeastSquareFit.pdf for the detail
 * z = coef[0] x + coef[1] y + coef[2], the plane equation
 * the normal vector to a plane $ax + by + cz +d = 0 $ is (a, b, c) using the gradient
 * so in this case, it is (coef[0], coef[1] , -1)
 */

inline bool fit2Plane(const vector<float>&coords, double *coef, bool *onXYplane,
		double *center ){
	const size_t dim = 3;
	const size_t n = (size_t) coords.size() / dim;
	if( n < 3){
		cout<<"Must have more than three points for defining a plane\n";
		return false;
	}
	double x = 0.0, y = 0.0, z = 0.0;
	double xy = 0.0, xz= 0.0, yz = 0.0;
	double xx = 0.0, yy= 0.0, zz = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x += coords[ dim * i ];
		y += coords[ i * dim + 1];
		z += coords[ i * dim + 2];
		xx += coords[ i * dim ] * coords[ i * dim ];
		yy += coords[ i * dim + 1 ] * coords[ i * dim +1];
		zz += coords[ i * dim + 2] * coords[ i * dim + 2];
		xy += coords[ i * dim ] * coords[ i *dim + 1];
		xz += coords[ i * dim ] * coords[ i * dim + 2];
		yz += coords[ i * dim + 1] * coords[ i *dim + 2];
	}

	center[ 0 ] = x / n;
	center[ 1 ] = y / n;
	center[ 2 ] = z / n;

	// in the xy plane
	if ( fabs(z) < eps1 ){
		*onXYplane = true;
		coef[0] = coef[1] = 0.0;
		coef[2] = 1.0;
		return true;
	}
	*onXYplane = false;

	double a_data[] = { xx, xy, x,
					                 xy, yy, y,
					                 x,   y,  (double)n};
	double b_data[] = { xz, yz, z };
	gsl_matrix_view m  = gsl_matrix_view_array ( a_data, dim, dim );
	gsl_vector_view b    = gsl_vector_view_array ( b_data, dim );
	gsl_vector *solution = gsl_vector_alloc (dim);
	int s;
	gsl_permutation *p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );

//	printf ("coef = \n");
//	gsl_vector_fprintf (stdout, solution, "%g");

	for (unsigned int k = 0; k < dim; k++)
		coef[k] = gsl_vector_get(solution, k);

	gsl_permutation_free (p);
	gsl_vector_free ( solution );

	return true;
}

/**
 * See the document of LeastSquareFit.pdf for the detail
 * z = coef[0] x + coef[1] y + coef[2], the plane equation
 * the normal vector to a plane $ax + by + cz +d = 0 $ is (a, b, c) using the gradient
 * so in this case, it is (coef[0], coef[1] , -1)
 * We need to fit to three different planes: however, it may be somewhat trickier
 * since for the case of a sphere, there is no unique plane
 * However, when return the coefficients in $coef$
 * float[] coef ={a, b, c, d},
 *
 * Normal could be both righ-handed and left handed
 */
inline bool fit2Plane( const vector<float>&coords, bool &onXYplane,
		vector<double>& center, vector<double>& coef, double &residual,
		bool debug = false) {

	const size_t n = (size_t) coords.size() / dim3;
	if( n < dim3 ) {
		cout<<"Must have more than three points for defining a plane\n";
		return false;
	}
	double x = 0.0, y = 0.0, z = 0.0,  xt = 0.0, yt = 0.0, zt = 0.0;
	double xy = 0.0, xz= 0.0, yz = 0.0;
	double xx = 0.0, yy= 0.0, zz = 0.0;

	for(unsigned int i = 0; i < n; i++){
		x = coords[ dim3 * i ];
		y = coords[ i * dim3 + 1 ] ;
		z = coords[ i * dim3 + 2 ] ;
		xt += x ;
		yt += y ;
		zt += z ;
		xx += x* x ;
		yy += y * y;
		zz += z * z ;
		xy += x * y;
		xz += x * z;
		yz +=  y * z;
	}
	center[ 0 ] = xt / n;
	center[ 1 ] = yt / n;
	center[ 2 ] = zt / n;

	// in the xy plane
	if ( fabs(z) < eps1 ){
		onXYplane = true;
		residual = sqrt ( fabs(z) / n );
		coef[0] = coef[1] = 0.0;
		coef[2] = 1.0;
		return true;
	}
	onXYplane = false;

	// CASE I: FIT to place  0 = solution[0] x + solution[1] y - z + solution[2]
	// or  z = solution[0] x + solution[1] y + solution[2]
	double a_data[] = { xx, xy, xt,
					                 xy, yy, yt,
					                 xt,   yt,  (double)n};
	double b_data[] = { xz, yz, zt };

	gsl_matrix_view m  = gsl_matrix_view_array (a_data, dim3, dim3);
	gsl_vector_view b = gsl_vector_view_array (b_data, dim3);
	gsl_vector *solution = gsl_vector_alloc (dim3);
	gsl_permutation *p = gsl_permutation_alloc (dim3);

	int s;
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );

	coef[0] = gsl_vector_get(solution, 0) ;
	coef[1] = gsl_vector_get(solution, 1) ;
	coef[2] = -1.0;
	coef[3] = gsl_vector_get(solution, 2) ;

	residual = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x = coords[ dim3 * i ];
		y = coords[ dim3 * i + 1];
		z = coords[ dim3 * i + 2];
		residual+=(coef[0] * x+coef[1] * y+coef[3] - z ) * (coef[0] * x +coef[1]*y+coef[3]- z) ;
	}
	if (debug)
		printf("Fit to a XY plane, residual 2 =%10.5f\n", residual );

	// CASE 2: FIT to place  0 = -x + solution[0] y + solution[1] z + solution[2]
	// or x = solution[0] y + solution[1] z + solution[2]
	a_data[0] = yy;
	a_data[1] = yz;
	a_data[2] = yt;
	a_data[3] = yz;
	a_data[4] = zz;
	a_data[5] = zt;
	a_data[6] = yt;
	a_data[7] = zt;

	b_data[0] = xy;
	b_data[1] = xz;
	b_data[2] = xt;

	m  = gsl_matrix_view_array (a_data, dim3, dim3);
	b = gsl_vector_view_array (b_data, dim3);
	s = 0;
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );
	double abc[dim3];
	for (unsigned int k = 0; k < dim3; k++)
		abc[k] = gsl_vector_get(solution, k) ;
	double residual2 = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x = coords[ dim3 * i ];
		y = coords[ dim3 * i + 1];
		z = coords[ dim3 * i + 2];
		residual2 +=(abc[0] * y+abc[1] * z+abc[2] - x) * (abc[0] * y +abc[1]*z+abc[2]- x) ;
	}

	if( residual2 < residual ) {
		coef[0] =  -1.0;
		coef[1] = gsl_vector_get(solution, 0) ;
		coef[2] = gsl_vector_get(solution, 1) ;
		coef[3] = gsl_vector_get(solution, 2) ;
		residual = residual2;
	}
	if (debug)
		printf("Fit to a YZ plane, residual=%10.5f\n", residual2 );

	// CASE 3: FIT to place  0 =  solution[1] x - y + solution[0] z + solution[2]
	// or y =solution[0] z +  solution[1] x + solution[2]
	a_data[0] = zz;
	a_data[1] = xz;
	a_data[2] = zt;
	a_data[3] = xz;
	a_data[4] = xx;
	a_data[5] = xt;
	a_data[6] = zt;
	a_data[7] = xt;

	b_data[0] = yz;
	b_data[1] = xy;
	b_data[2] = yt;

	m  = gsl_matrix_view_array (a_data, dim3, dim3);
	b = gsl_vector_view_array (b_data, dim3);
	s = 0;
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );
	for (unsigned int k = 0; k < dim3; k++)
		abc[k] = gsl_vector_get(solution, k) ;

	residual2 = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x = coords[ dim3 * i ];
		y = coords[ dim3 * i + 1];
		z = coords[ dim3 * i + 2];
		residual2 +=(abc[0] *z+abc[1]* x+abc[2] - y) * (abc[0] *z +abc[1]*x+abc[2]- y) ;
	}
	if( residual2 < residual ) {
		coef[0] = gsl_vector_get(solution, 1) ;
		coef[1] =  -1.0;
		coef[2] = gsl_vector_get(solution, 0) ;
		coef[3] = gsl_vector_get(solution, 2) ;
		residual = residual2;
	}
	if (debug)
		printf("Fit to a XZ plane, residual=%10.5f\n",  residual2 );

	gsl_permutation_free (p);
	gsl_vector_free ( solution );

//	if ( canBeProjected ){
//		for ( unsigned int i = 0; i < n; ++i )
//			printf("%10.3f%10.3f%10.3f%5s%10.3f%10.3f%10.3f\n",
//					atomCoordVec[i][0], atomCoordVec[i][1], atomCoordVec[i][2],
//					"---", pAtCoordVec[ i][0],  pAtCoordVec[i][1], pAtCoordVec[ i][2] );
//	}
	return true;
}


/**
 * See the document of LeastSquareFit.pdf for the detail
 * z = coef[0] x + coef[1] y + coef[2], the plane equation
 * the normal vector to a plane $ax + by + cz +d = 0 $ is (a, b, c) using the gradient
 * so in this case, it is (coef[0], coef[1] , -1)
 */

inline bool fit2Plane(const vector<float>&coords, double *coef, bool &onXYplane,
		double *center, double &residual ){
	const size_t dim3 = 3;
	const size_t n = (size_t) coords.size() / dim3;
	if( n < 3){
		cout<<"Must have more than three points for defining a plane\n";
		return false;
	}
	double x = 0.0, y = 0.0, z = 0.0;
	double xy = 0.0, xz= 0.0, yz = 0.0;
	double xx = 0.0, yy= 0.0, zz = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x += coords[ dim3 * i ];
		y += coords[ i * dim3 + 1];
		z += coords[ i * dim3 + 2];
		xx += coords[ i * dim3 ] * coords[ i * dim3 ];
		yy += coords[ i * dim3 + 1 ] * coords[ i * dim3 +1];
		zz += coords[ i * dim3 + 2] * coords[ i * dim3 + 2];
		xy += coords[ i * dim3 ] * coords[ i *dim3 + 1];
		xz += coords[ i * dim3 ] * coords[ i * dim3 + 2];
		yz += coords[ i * dim3 + 1] * coords[ i *dim3 + 2];
	}

	center[ 0 ] = x / n;
	center[ 1 ] = y / n;
	center[ 2 ] = z / n;
	// in the xy plane
	if ( fabs(z) < eps1 ){
		onXYplane = true;
		residual = sqrt ( fabs(z) / n );
		coef[0] = coef[1] = 0.0;
		coef[2] = 1.0;
		return true;
	}

	onXYplane = false;

	double a_data[] = { xx, xy, x,
					                 xy, yy, y,
					                 x,   y,  (double)n};
	double b_data[] = { xz, yz, z };
	gsl_matrix_view m  = gsl_matrix_view_array (a_data, dim3, dim3);
	gsl_vector_view b = gsl_vector_view_array (b_data, dim3);
	gsl_vector *solution = gsl_vector_alloc (dim3);
	int s;
	gsl_permutation * p = gsl_permutation_alloc (dim3);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );

//	printf ("coef = \n");
//	gsl_vector_fprintf (stdout, solution, "%g");

	double normal2Plane [dim3];
	for (unsigned int k = 0; k < dim3; k++)
		normal2Plane[k] = coef[k] = gsl_vector_get(solution, k);
	normal2Plane[2]= -1.0;
	double v1LenSq = 0.0;
	for (unsigned int k = 0; k < dim3; k++)
		v1LenSq += normal2Plane[k] * normal2Plane[k] ;
	double norm = sqrt( v1LenSq );
	for (unsigned int k = 0; k < dim3; k++)
		normal2Plane[k] /= norm;

	residual = 0.0;
	double coord[dim3];
	double length2center = 0.0;
	double crossProduct = 0.0;
	for(unsigned int i = 0; i < n; i++) {
		coord[0] = coords[ dim3 * i ] - center[0];
		coord[1] = coords[ i * dim3 + 1 ] - center[1];
		coord[2] = coords[ i * dim3 + 2 ] - center[2];
		length2center = 0.0;
		for( unsigned int s = 0; s < dim3; s++)
			length2center += coord[s] * coord[s];
		if ( length2center < Tiny )
			residual += Tiny;
		else{
			length2center = sqrt( length2center );
			crossProduct = 0.0;
			for (unsigned int s = 0; s < dim3; s++)
				crossProduct += (normal2Plane [s] * coord[s]) / length2center;
			residual += fabs(crossProduct);
		}
	}
	residual /= n ;

	gsl_permutation_free (p);
	gsl_vector_free ( solution );

	return true;
}


/**
 * See the document of LeastSquareFit.pdf for the detail
 * z = coef[0] x + coef[1] y + coef[2], the plane equation
 * the normal vector to a plane $ax + by + cz +d = 0 $ is (a, b, c) using the gradient
 * so in this case, it is (coef[0], coef[1] , -1)
 */

inline bool fit2Plane( const vector<double>& coords, double *coef, bool *onXYplane,
		double *center ){
	const size_t dim = 3;
	const size_t n = (size_t) coords.size() / dim ;
	if( n < 3){
		cout<<"Must have more than three points for defining a plane\n";
		return false;
	}
	double x = 0.0, y = 0.0, z = 0.0;
	double xy = 0.0, xz= 0.0, yz = 0.0;
	double xx = 0.0, yy= 0.0, zz = 0.0;
	for(unsigned int i = 0; i < n; i++){
		x += coords[ dim * i ];
		y += coords[ i * dim + 1];
		z += coords[ i * dim + 2];
		xx += coords[ i * dim ] * coords[ i * dim ];
		yy += coords[ i * dim + 1 ] * coords[ i * dim +1];
		zz += coords[ i * dim + 2] * coords[ i * dim + 2];
		xy += coords[ i * dim ] * coords[ i *dim + 1];
		xz += coords[ i * dim ] * coords[ i * dim + 2];
		yz += coords[ i * dim + 1] * coords[ i *dim + 2];
	}


	center[ 0 ] = x / n;
	center[ 1 ] = y / n;
	center[ 2 ] = z / n;
	// in the xy plane
	if ( fabs(z) < eps1 ){
		*onXYplane = true;
		coef[0] = coef[1] = 0.0;
		coef[2] = 1.0;
		return true;
	}

	*onXYplane = false;
	double a_data[] = { xx, xy, x,
					                 xy, yy, y,
					                 x,   y,  (double)n};
	double b_data[] = { xz, yz, z };
	gsl_matrix_view m  = gsl_matrix_view_array (a_data, dim, dim);
	gsl_vector_view b = gsl_vector_view_array (b_data, dim);
	gsl_vector *solution = gsl_vector_alloc (dim);
	int s;
	gsl_permutation * p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, solution );

//	printf ("coef = \n");
//	gsl_vector_fprintf (stdout, solution, "%g");

	for (unsigned int k = 0; k < dim; k++)
		coef[k] = gsl_vector_get(solution, k);

	gsl_permutation_free (p);
	gsl_vector_free ( solution );

	return true;
}



// area2D_Polygon(): computes the area of a 2D polygon
//    Input:  int n = the number of vertices in the polygon
//            Point* V = an array of n+2 vertices
//                       with V[n]=V[0] and V[n+1]=V[1]
//    Return: the (float) area of the polygon
inline double area2DPolygon( double *V, const size_t n ) {
	double area = 0;
	unsigned int   i, j, k;     // indices

	for (i = 1, j = 2, k = 0; i <= n; i++, j++, k++)
		area += V[ 2 * i  ] * ( V[ 2 * j + 1 ] - V[ 2 * k + 1 ] );
	return area / 2.0;
}
//===================================================================

// area3D_Polygon(): computes the area of a 3D planar polygon
//    Input:  int n = the number of vertices in the polygon
//            Point* V = an array of n+2 vertices in a plane
//                       with V[n]=V[0] and V[n+1]=V[1]
//            Point N = unit normal vector of the polygon's plane
//    Return: the (float) area of the polygon
inline double  area3DPolygon( double *V, double *N, const size_t n ) {
	double area = 0;
	const size_t dim = 3;
	double  an, ax, ay, az;  // abs value of normal and its coords
	int   coord;           // coord to ignore: 1=x, 2=y, 3=z
	unsigned int   i, j, k;         // loop indices

	// select largest abs coordinate to ignore for projection
	ax = ( N[0] > 0 ? N[0] : -N[0] );     // abs x-coord
	ay = ( N[1] > 0 ? N[1] : -N[1] );     // abs y-coord
	az = ( N[2] > 0 ? N[2] : -N[2] );     // abs z-coord

	coord = 3;                     // ignore z-coord
	if (ax > ay) {
		if (ax > az)
			coord = 1;    // ignore x-coord
	} else if (ay > az)
		coord = 2;       // ignore y-coord

	// compute area of the 2D projection
	for (i=1, j=2, k=0; i<=n; i++, j++, k++) {
		switch (coord) {
		case 1:
			area += ( V[dim *i +1 ] * ( V[dim *j +2 ] - V[dim *k + 2 ] ) );
			continue;
		case 2:
			area += ( V[dim * i ] * ( V[dim *j +2 ] - V[dim *k + 2 ] ) );
			continue;
		case 3:
			area += ( V[dim * i ] * ( V[ dim * j + 1] - V[dim * k + 1] ) );
			continue;
		}
	}

	// scale to get area before projection
	an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
	switch (coord) {
	case 1:
		area *= (an / (2*ax));
		break;
	case 2:
		area *= (an / (2*ay));
		break;
	case 3:
		area *= (an / (2*az));
	}
	return area;
}


/**
 * This implementation is mainly for the computation of the volume of a tetrahedron
 * so the rank should be six ($d_{12}^2, d_{13}^2, d_{14}^2, d_{23}^2,d_{24}^2, d_{34}^2)$
 */
inline bool CayleyMengerDeterminant( double *const dist, const size_t rank, double &volume){
	const size_t expectedRank = 6;
	if (rank != expectedRank )
		return false;
	volume = 0.0;
	float d12 = dist[0], d13 =dist[1], d14 =dist[2], d23 =dist[3], d24 =dist[4], d34 =dist[5];
	const size_t  datDim = 9;
	double mat[ datDim ];
	mat[0] = - d12 * d12;
	mat[1] =   d23 * d23 - d13  * d13;
	mat[2] =   d24 * d24 - d14  * d14;
	mat[3] =   d23 * d23 - d12  * d12;
	mat[4] = - d13 * d13;
	mat[5] =   d34 * d34 - d14 * d14;
	mat[6] =   d24 * d24 - d12 * d12;
	mat[7] =   d34 * d34 - d13 * d13;
	mat[8] = - d14 * d14;
	volume = - ( mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]
	                           - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ] ) ;

	mat[0] =   d12 * d12;
	mat[1] =   d23 * d23 - d13  * d13;
	mat[2] =   d24 * d24 - d14  * d14;
	mat[3] =   d13 * d13;
	mat[4] = - d13 * d13;
	mat[5] =   d34 * d34 - d14 * d14;
	mat[6] =   d14 * d14;
	mat[7] =   d34 * d34 - d13 * d13;
	mat[8] = - d14 * d14;
	volume += (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]
		        - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ] );

	mat[0] =   d12 * d12;
	mat[1] =  - d12 * d12;
	mat[2] =   d24 * d24 - d14  * d14;
	mat[3] =   d13 * d13;
	mat[4] =   d23 * d23 - d12  * d12;
	mat[5] =   d34 * d34 - d14 * d14;
	mat[6] =   d14 * d14;
	mat[7] =   d24 * d24 - d12 * d12;
	mat[8] = - d14 * d14;
	volume -=  (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]
		        - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ] );

	mat[0] =   d12 * d12;
	mat[1] =  - d12 * d12;
	mat[2] =   d23 * d23 - d13 * d13;
	mat[3] =   d13 * d13;
	mat[4] =   d23 * d23 - d12  * d12;
	mat[5] = - d13 * d13;
	mat[6] =   d14 * d14;
	mat[7] =   d24 * d24 - d12 * d12;
	mat[8] =   d34 * d34 - d13 * d13;
	volume  += (mat[ 0 ] * mat[ 4 ] * mat[ 8 ] + mat[ 3] * mat[ 7] * mat[ 2]  + mat[ 6] * mat[ 5] * mat[ 1]
			        - mat[ 2]  * mat[ 4] *  mat[ 6] -mat[ 5]  * mat[ 7] *  mat[ 0] -mat[ 8]  * mat[ 3] *  mat[ 1 ] ) ;

	return true;
}

inline void fixPointApproximation ( double *radius, double *radiusArr, double  const *center,
		double *abc, const vector<float>& coords ){
//compute the average radius
	const size_t dim = 3;
	size_t m = (size_t) ( coords.size() / dim );
	if( m == 0 ){
		cout<<"no data: "<< m <<endl ;
		exit(1);
	}
	double a = 0.0, b = 0.0, c = 0.0;
	for (unsigned int i = 0; i < m; i++ ){
		a += ( abc[0] - coords[ dim *i ]) / radiusArr[ i ];
		b += ( abc[1] - coords[ dim *i +1 ] ) / radiusArr[ i ];
		c += ( abc[2] - coords[ dim *i +2 ] ) / radiusArr[ i ];
	}
	abc[0] = center[0] + (*radius) * a / m;
	abc[1] = center[1] + (*radius) * b / m;
	abc[2] = center[2] + (*radius) * c / m;
//	cout<<"Before ="<< (*radius) <<endl;
	double  r = 0.0, dis = 0.0 ;
	double coord[dim];
	for( unsigned int i = 0; i < m; i++){
		dis = 0.0;
		for ( unsigned int p = 0; p< dim; p++) {
			coord[ p ] = coords[ dim * i + p ] - abc[p] ;
			dis += coord[ p ] * coord[ p ];
		}
		radiusArr[ i ] = sqrt(dis);
		r += radiusArr[ i ] ;
	}
	*radius = r / m;
//	cout<<abc[0]<<"  "<<abc[1]<<"  "<<abc[2]<<" After ="<< (*radius) <<endl;
}

inline double fit2Sphere( const vector<float>& coords,  double  *abc){
	const size_t dim = 3;
	size_t m = (size_t) ( coords.size() / dim );
	if( m == 0 ){
		cout<<"no data: "<< m <<endl ;
		exit(1);
	}

	// Compute the center
	double  center[dim];
	size_t p = 0;
	for ( p = 0; p< dim; p++)
		center[ p ] = 0.0;
	for ( unsigned int i = 0; i < m; i++)
		for ( p = 0; p< dim; p++)
			center[p] += coords[ dim * i + p ];
	for ( p = 0; p< dim; p++)
		abc[p] = center[p] /= m;

//	unsigned int fpIndex = 0;
	float dis = 0.0;
	double *radiusArr = new double[m];
	double  r = 0.0;
	double coord[dim];
	for( unsigned int i = 0; i < m; i++){
		dis = 0.0;
		for ( p = 0; p< dim; p++) {
			coord[ p ] = coords[ dim * i + p ] - center[p] ;
			dis += coord[ p ] * coord[ p ];
		}
		radiusArr [ i ] = sqrt(dis);
		r += radiusArr[ i ] ;
	}
	double radius = r / m;
//	cout<<"Does radius improve ?\n ";
	const size_t nCycle = 32;
	size_t n = 0;
	while ( n < nCycle ){
		fixPointApproximation ( &radius, radiusArr, center, abc, coords );
		n++;
//		cout<<radius<<endl;
	}
	return radius;
}


/**
 * Both $arrA$ and $arrB$ must be centered first
 @param m: the size of the atoms in the structure
 */
inline matrix rotMatrixByRMSD4C(const vector<float>&arrA, const vector<float>&arrB,
		bool& isSingularZero, double& rmsVal ){

	const size_t m = (size_t) ( arrA.size() / dim3 );
	double  arrU[ dim3 * dim3 ];
	for (unsigned int i = 0; i < dim3 ; i++ ){
		for  (unsigned int j = 0; j < dim3; j++) {
			arrU[ i * dim3 + j ] = 0.0;
			for (unsigned int k = 0; k < m; k++)
				arrU[ i * dim3 + j ] += arrB [ k * dim3 + j ] * arrA[ k * dim3 + i ];
		}
	}

	gsl_matrix_view U	= gsl_matrix_view_array (arrU, dim3, dim3);
	gsl_matrix *V = gsl_matrix_alloc (dim3, dim3 );
	gsl_vector *sValues = gsl_vector_alloc ( dim3 );
	gsl_vector *work = gsl_vector_alloc ( dim3 );

	gsl_linalg_SV_decomp (&U.matrix, V, sValues, work);
	for (unsigned int k = 0; k < dim3; k++){
		if ( fabs ( gsl_vector_get( sValues, k ) ) < 1.0E-12) {
			isSingularZero = true;
			rmsVal= 0.0;
			gsl_matrix_free (V);
			gsl_vector_free (work);
			gsl_vector_free (sValues);
			return matrix(dim3, dim3);
		}
	}

	// matrix rot = U.times( V.transpose() );
	const size_t dim9 = 9;
	vector<double> rotArray ( dim9, 0.0);
	for (unsigned int i = 0; i < dim3; i++ )
		for  (unsigned int j = 0; j < dim3; j++)
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;

	double rms = 0.0;
	double arrB2[ dim3 ];
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ]; //A[j][k] * B[k];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}

	double detValue = det(rotArray);
	if( fabs( detValue + 1.0) > 1.0E-8)  {
		rmsVal = rms;
		isSingularZero = false;
		gsl_matrix_free (V);
		gsl_vector_free (work);
		gsl_vector_free (sValues);
		return matrix(rotArray, dim3, dim3);
	}

//	cout<<"left handedness !\n";
	double arrV[ dim9 ];
	for (unsigned int i = 0; i < dim3; i++ )
		for  (unsigned int j = 0; j < dim3; j++)
			arrV[ dim3 * i + j ] = gsl_matrix_get (V, i, j);
	gsl_matrix_view Vv = gsl_matrix_view_array (arrV, dim3, dim3);

	// reverse the x-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 0,  - gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 1,  gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
		gsl_matrix_set ( V, i, 2,  gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}

	double minRMS = rms;
	size_t indexOfMinRms = 0;
	// reverse the y-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 1,  - gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 0,  gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
		gsl_matrix_set ( V, i, 2,  gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray [ p * dim3 + q ] * arrB [ dim3 * i + q ];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}
	if( rms < minRMS ) {
		minRMS = rms;
		indexOfMinRms = 1;
	}

	// reverse the z-axis
	for (unsigned int i = 0; i< dim3; i++ )
		gsl_matrix_set ( V, i, 2,  - gsl_matrix_get (&Vv.matrix,  i, 2) ) ;
	for (unsigned int i = 0; i< dim3; i++ ) {
		gsl_matrix_set ( V, i, 0,  gsl_matrix_get (&Vv.matrix,  i, 0) ) ;
		gsl_matrix_set ( V, i, 1,  gsl_matrix_get (&Vv.matrix,  i, 1) ) ;
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	rms = 0.0;
	for (unsigned int i = 0; i < m; i++){
		for (unsigned int p = 0; p < dim3; p++){
			arrB2[ p ] = 0.0;
			for ( unsigned int q = 0; q < dim3; q++)
				arrB2 [ p ]+= rotArray[ p * dim3 + q ] * arrB [ dim3 * i + q ]; //A[j][k] * B[k];
			rms += (  (arrB2[p] - arrA[ dim3 * i + p ] ) * (arrB2[ p ] - arrA[ i * dim3 + p ] ) );
		}
	}
	if( rms < minRMS ) {
		minRMS = rms;
		indexOfMinRms = 2;
	}

	rmsVal = minRMS;
	for (unsigned int i = 0; i< dim3; i++ ) {
		for (unsigned int j = 0; j < dim3; ++j )  {
			if ( j != indexOfMinRms)
				gsl_matrix_set ( V, i, j,  gsl_matrix_get (&Vv.matrix,  i, j) ) ;
			else
				gsl_matrix_set ( V, i, j,  - gsl_matrix_get (&Vv.matrix,  i, j) ) ;
		}
	}
	for (unsigned int i = 0; i < dim3; i++ ) {
		for  (unsigned int j = 0; j < dim3; j++) {
			rotArray [ i * dim3 + j ]  = 0.0;
			for (unsigned int k = 0; k< dim3; k++)
				rotArray [ i * dim3 + j ] += gsl_matrix_get (&U.matrix, i, k) * gsl_matrix_get(V, j, k)  ;
		}
	}
	isSingularZero = false;
	gsl_matrix_free (V);
	gsl_vector_free (work);
	gsl_vector_free (sValues);
	return matrix(rotArray, dim3, dim3);
}

/**
 * BY the following figure, from the length of ca2ca ( $intervalS$ that is rather fixed,
 * about 3.82\AA) we can compute the extension angle $ta$ from $b$ and $a$
 * since the center of the helix circle must lies on the line that bisect P1 and P2j
 * and perpendicular to the line from P1P2j, so
 * we have $\sin{t_a/2} = \frac{ L_{P1P2j} }{2a}
 * The fitted origin $fOrigin$ could be back-computed from the standard helix origin
 * The fitted normal $fNorm$ is computed from the z-axis of the standard helix frame
 * The standard helix sits on XY plane with the helix axis along +Z axis
 * and the helix origin is the origin of its frame
 *                                                         P2
 *                                   center         /   |
 *                                       |        /         |
 *                                        | /               |
 *                                    /    |                |   b
 *                               /          | a            |
 *                         /___________|__________|
 *                     P1    $pLength$            P2j (the projection of P2 onto the helix plane )
 *
 * Typical value for protein helix
 * const float b0 = 1.10, deltaB = 0.02; // nb = 30, so  1.2 -- 1.2+0.6
 * const float a0 = 1.95, deltaA = 0.02; // 1.728-- 2.728
 * 	const size_t na = 40;
	const size_t nb = 40;
 * 	const size_t resolution = 20;
 */

inline bool fit4Points2Helix ( const vector<float>& helixVec,
		const float b0, const float deltaB, const size_t nb, // parameters for helix pitch
		const float a0, const float deltaA, const size_t na, // parameters for helix circle
		const size_t resolution,
		vector<float>& p12Vec, vector<float>& p23Vec, vector<float>& p34Vec,
		vector<float>& fOrigin, vector<float>& fNorm,
		vector<float>& fittedHelicalParam, bool debug=false) {

	vector<float> ghelixVec = helixVec;
	const size_t noOfPoints = (size_t) ( ghelixVec.size() / dim3 );
	if( noOfPoints < dim4 ) {
		cout<<"Need at least 4 points for fitting into helix !\n";
		return false;
	}
	// $p1, p2, p3$ and $p4$ save the original 4 points
	vector<float> p1(dim3, 0.0),  p2(dim3, 0.0), p3(dim3, 0.0), p4(dim3, 0.0);
	vector<float> gCenter(dim3, 0.0);
	float dis12 = 0.0, dis23= 0.0, dis34= 0.0;
	unsigned int s = 0;

	for ( s = 0; s < dim3; ++s ) {
		p1[s] = ghelixVec[ s];
		p2[s] = ghelixVec[s+dim3];
		p3[s] = ghelixVec[s+dim6];
		p4[s] = ghelixVec[s+dim9];
		gCenter[s] += 	p1[s]  + p2[s] +p3[s] + p4[s] ;
		dis12 += ( p2[s] - p1[s] ) * ( p2[s] - p1[s] );
		dis23 += ( p2[s] - p3[s] ) * ( p2[s] - p3[s] );
		dis34 += ( p3[s] - p4[s] ) * ( p3[s] - p4[s] );
	}
	const float interval = (sqrt(dis12) + sqrt(dis23) + sqrt(dis23) ) / 3.0;
//	if( debug )
//		printf("Length of P12= %10.6f\n", interval) ;

	for ( s = 0; s < dim3; s++)
		gCenter[s] /=  noOfPoints;

	for ( unsigned int i = 0; i < noOfPoints; ++i ) {
		for ( s = 0; s < dim3; s++)
			ghelixVec[ dim3 * i  + s ] -= gCenter[s];
	}

	float b = 0.0, a = 0.0;
	vector<float> sHelixVec ( dim3 * noOfPoints, 0.0);
	vector<float> sCenter(dim3, 0.0);
	vector<float> fCenter(dim3, 0.0);
	float ta=0.0, pLength=0.0;
	float st =0.0, ct=0.0, st2= 0.0, ct2= 0.0, st3= 0.0, ct3= 0.0, ss=0.0, cc=0.0;
	bool  isSingularZero;
	matrix rot, fRot;
	int minP = -1, minQ = -1;
	double minRMS  = 100.0 / smallEPS;
	double rmsVal = 0.0;
	const double intervalS =  interval * interval ;
	vector<float> fHelixVec;

	for (unsigned int p = 0; p < na; ++p) {
		a = a0 + p * deltaA;

		for ( unsigned int q = 0; q < nb; ++q ) {
			b = b0 + q * deltaB;

			pLength = sqrt( intervalS  - b * b );
			ta = 2.0 * asin( 0.5 * pLength / a );
			st = sin(ta);
			ct = cos(ta);
			ss = st * st;
			cc = ct * ct;
			st2 = 2.0 * st * ct;   // \sin{2x} = 2 \sin{x} \cos{x}
			ct2 = 2.0 * cc - 1.0; // \cos{2x} = 2 \cos{x} \cos{x} - 1
			st3 = 3.0 * st  - 4.0 * ss * st; // \six{3x} = 3 \sin{x} -4 \sin^3 {x}
			ct3 = 4.0 * ct * cc - 3.0 * ct; // \cos{3x} = 4 \cos^3 {x} - 3 \cos{x}

			sHelixVec[0] = a;
			sHelixVec[1] = 0.0;
			sHelixVec[2] = 0.0;

			sHelixVec[3] = a * ct;
			sHelixVec[4] = a * st;
			sHelixVec[5] = b;

			sHelixVec[6] = a * ct2;
			sHelixVec[7] = a * st2;
			sHelixVec[8] = 2.0 * b;

			sHelixVec[9]   = a * ct3;
			sHelixVec[10] = a * st3;
			sHelixVec[11] = 3.0 * b;

			sCenter[2] = sCenter[1] = sCenter[0] = 0.0;
			for ( s = 0; s < dim3; s++)
				sCenter[s] += sHelixVec[dim6+s]+sHelixVec[dim9+s]
				                                          +sHelixVec[dim3+s]+sHelixVec[s];

			for ( s = 0; s < dim3; s++)
				sCenter[s] /=  noOfPoints;
			for ( unsigned int i = 0; i < noOfPoints; ++i ) {
				for (  s = 0; s < dim3; s++)
					sHelixVec[ dim3 * i  + s ] -= sCenter[s];
			}
			rmsVal = 0.0;
			isSingularZero = true;
			rot = rotMatrixByRMSD4C(ghelixVec, sHelixVec, isSingularZero, rmsVal );
			if ( ! isSingularZero ) {
				if ( rmsVal < minRMS ) {
					minRMS = rmsVal;
					minP = p;
					minQ = q;
					fHelixVec = sHelixVec;
					fRot = rot;
					fCenter = sCenter;
				}
			}
		}
	}

	if( minP == -1) {
		cerr<<"can NOT fit to a helix !\n";
		return false;
	}
	fittedHelicalParam[0] = a = a0 + minP * deltaA;
	fittedHelicalParam[1] = b = b0 + minQ * deltaB;
	pLength = sqrt( intervalS - b * b );
	fittedHelicalParam[2] = ta = dim2 * asin( 0.5 * pLength / a );
	fittedHelicalParam[3] = sqrt ( minRMS / noOfPoints);
	float t = 0.0;
	float tta = ta / resolution;
	const float bb = b / resolution;
	vector<float> coord(dim3, 0.0), coordN(dim3, 0.0);

	p12Vec = vector<float> ( dim3 * resolution, 0.0 );
	// the back-computed $p1$
	for ( s = 0; s < dim3; s++)
		coord[s] = fHelixVec[s];
	coordN = fRot.times(coord);
	for( s = 0; s<dim3; ++s)
		p12Vec[s] = coordN[s] + gCenter[s];

	for ( unsigned int i = 1; i < resolution; ++i ) {
		t = i * tta;
		coord[0] = a *  cos(t) - fCenter[0];
		coord[1] = a *  sin(t)  - fCenter[1];
		coord[2] = bb * i  - fCenter[2];
		coordN = fRot.times(coord);
		for( s = 0; s<dim3; ++s )
			p12Vec[ dim3 * i + s] = coordN[s] + gCenter[s];
	}

	// the back-computed $p2$
	p23Vec = vector<float> (dim3 * resolution, 0.0);
	for ( s = 0; s < dim3; s++)
		coord[s] = fHelixVec[s+dim3];
	coordN = fRot.times(coord);
	for( s = 0; s<dim3; ++s)
		p23Vec[s] = coordN[s] + gCenter[s];

	for ( unsigned int i = 1; i < resolution; ++i ) {
		t = i * tta + ta;
		coord[0] = a *  cos(t) - fCenter[0];
		coord[1] = a *  sin(t)  - fCenter[1];
		coord[2] = b + bb * i  - fCenter[2];
		coordN = fRot.times(coord);
		for( s = 0; s<dim3; ++s)
			p23Vec[ dim3 * i + s] = coordN[s] + gCenter[s];
	}

	p34Vec = vector<float> (dim3 * resolution, 0.0);
	// the back-computed $p3$
	for ( s = 0; s < dim3; s++)
		coord[s] = fHelixVec[s+dim6];
	coordN = fRot.times(coord);
	for( s = 0; s<dim3; ++s )
		p34Vec[s] = coordN[s] + gCenter[s];
	ta *= 2.0;
	b *= 2.0;
	for ( unsigned int i = 1; i < resolution; ++i ) {
		t = i * tta + ta;
		coord[0] = a *  cos(t) - fCenter[0];
		coord[1] = a *  sin(t)  - fCenter[1];
		coord[2] = b + bb * i  - fCenter[2];
		coordN = fRot.times(coord);
		for( s = 0; s<dim3; ++s)
			p34Vec[ dim3 * i + s] = coordN[s] + gCenter[s];
	}

	// compute $fOrigin$
	for( s = 0; s<dim3; ++s)
		coord[s] = 0.0 - fCenter[s];
	coordN = fRot.times(coord);
	for( s = 0; s<dim3; ++s)
		fOrigin[s] = coordN[s] + gCenter[s];

	// compute $fNorm$
	vector<float> unitZ(dim3, 0.0);
	unitZ[2] = 1.0;
	fNorm = fRot.times(unitZ);

	if(debug)  {
		printf("a=%7.3f, b=%7.3f, ta=%7.3f\n", a, b, ta * rad2Deg );
//		fRot.print();
		printf("RMSD=%8.4f\n", sqrt ( minRMS / noOfPoints)) ;
		// back-computation
		float rms = 0.0;
		rmsVal = 0.0;
		for (unsigned int i = 0; i < noOfPoints; ++i) {
			for(unsigned int s = 0; s<dim3; ++s)  {
				printf("%10.6f ", helixVec[dim3 * i +s ] );
				coord[s] = fHelixVec[dim3 * i +s ];
			}
			cout<<endl;
			rms = 0.0;
			coordN = fRot.times(coord);
			for(unsigned int s = 0; s<dim3; ++s) {
				sHelixVec[dim3 * i + s ] = coordN[s] + gCenter[s];
				rms += ( helixVec[dim3 * i +s ] - sHelixVec[dim3 * i +s ]  )
						*  ( helixVec[dim3 * i +s ] - sHelixVec[dim3 * i +s ]  );
				printf("%10.6f ", sHelixVec[dim3 * i +s ] );
			}
			rmsVal += rms;
			cout<<"rms = "<<sqrt(rms)<<endl;
		}
		printf("RMSD2 =%8.4f\n", sqrt ( rmsVal / noOfPoints)) ;
		cout<<endl;
	}
	return true;
}

inline bool fitnPoints2Helix ( const vector<float>& helixVec,
		const float b0, const float deltaB, const size_t nb, // parameters for helix pitch
		const float a0, const float deltaA, const size_t na, // parameters for helix circle
		const size_t resolution,
		vector<float>& fPointVec,  //vector<float>& p23Vec, vector<float>& p34Vec,
		vector<float>& fOrigin, vector<float>& fNorm,
		vector<float>& fittedHelicalParam,
		bool  isLefthanded = false, bool debug=false ) {

	vector<float> ghelixVec = helixVec;
	const size_t noOfPoints = (size_t) ( ghelixVec.size() / dim3 );
	if( noOfPoints < dim4 ) {
		cout<<"Need at least 4 points for fitting into helix !\n";
		return false;
	}

	vector<float> p1(dim3, 0.0),  p2(dim3, 0.0);
	vector<float> gCenter(dim3, 0.0);
	unsigned int s = 0;

	for ( s = 0; s < dim3; ++s )
		gCenter[s] = p1[s] = ghelixVec[ s];

	// compute the interatom Ca-Ca or P-P distances
	vector<float> disVec(noOfPoints - 1, 0.0 );
	for ( unsigned int i = 1; i < noOfPoints; ++i ) {
		disVec[ i - 1 ] = 0.0;
		for ( s = 0; s < dim3; ++s ) {
			p2[s] = ghelixVec[ i * dim3 + s];
			disVec[ i - 1 ] += ( p2[s] - p1[s] ) * ( p2[s] - p1[s] );
			gCenter[s] += p2[s];
		}
		p1 = p2;
	}
	float interval = 0.0;
	for ( unsigned int i = 0; i < disVec.size(); ++i )
		interval += sqrt(disVec[i]);
	interval /= (noOfPoints - 1);
	if( debug )
		printf("Length of P12= %10.6f\n", interval) ;

	// Move to its own center
	for ( s = 0; s < dim3; s++)
		gCenter[s] /=  noOfPoints;
	for ( unsigned int i = 0; i < noOfPoints; ++i ) {
		for ( s = 0; s < dim3; s++)
			ghelixVec[ dim3 * i  + s ] -= gCenter[s];
	}

	// Compute standard aligned with +Z helix
	float b = 0.0, a = 0.0;
	vector<float> sHelixVec ( dim3 * noOfPoints, 0.0);
	vector<float> sCenter(dim3, 0.0);
	vector<float> fCenter(dim3, 0.0);
	float ta=0.0, pLength=0.0;
//	float st =0.0, ct=0.0, st2= 0.0, ct2= 0.0, st3= 0.0, ct3= 0.0, ss=0.0, cc=0.0;
	bool  isSingularZero;
	matrix rot, fRot;
	int minP = -1, minQ = -1;
	double minRMS  = 100.0 / smallEPS;
	double rmsVal = 0.0, aveTa = 0.0;
	const double intervalS =  interval * interval ;
	vector<float> fHelixVec;

	for (unsigned int p = 0; p < na; ++p) {
		a = a0 + p * deltaA;

		for ( unsigned int q = 0; q < nb; ++q ) {
			b = b0 + q * deltaB;
			ta = 0.0;

//			pLength = sqrt( intervalS  - b * b );
//			ta = 2.0 * asin( 0.5 * pLength / a );

			sCenter[0] = sHelixVec[0] = a;
			sCenter[1] = sHelixVec[1] = 0.0;
			sCenter[2] = sHelixVec[2] = 0.0;
			for (unsigned k = 1; k < noOfPoints; ++k) {
				pLength = sqrt( disVec[k-1]  - b * b );
				ta += 2.0 * asin( 0.5 * pLength / a );
				sHelixVec[ k * dim3 ]       = a * cos( ta);
				sHelixVec[ k * dim3 + 1 ] = a * sin( ta);
//				sHelixVec[ k * dim3 ]       = a * cos( k* ta);
//				sHelixVec[ k * dim3 + 1 ] = a * sin( k * ta);
				sHelixVec[ k * dim3 + 2 ] = b * k ;
				sCenter[0] += sHelixVec[ k * dim3] ;
				sCenter[1] += sHelixVec[ k * dim3 + 1] ;
				sCenter[2] += sHelixVec[ k * dim3 + 2] ;
			}

			// Move to its own center
			for ( s = 0; s < dim3; s++)
				sCenter[s] /=  noOfPoints;
			for ( unsigned int i = 0; i < noOfPoints; ++i ) {
				for ( s = 0; s < dim3; s++)
					sHelixVec[ dim3 * i  + s ] -= sCenter[s];
			}

			rmsVal = 0.0;
			isSingularZero = true;
			rot = rotMatrixByRMSD4C(ghelixVec, sHelixVec, isSingularZero, rmsVal );
			if ( ! isSingularZero ) {
				if ( rmsVal < minRMS ) {
					minRMS = rmsVal;
					minP = p;
					minQ = q;
					fHelixVec = sHelixVec;
					fRot = rot;
					fCenter = sCenter;
					aveTa = ta / ( noOfPoints - 1);
				}
			}
		}
	}

	if( minP == -1)
		return false;
	fittedHelicalParam[0] = a = a0 + minP * deltaA;
	fittedHelicalParam[1] = b = b0 + minQ * deltaB;
//	pLength = sqrt( intervalS - b * b );
	fittedHelicalParam[2] = aveTa; // = dim2 * asin( 0.5 * pLength / a );
	fittedHelicalParam[3] = sqrt ( minRMS / noOfPoints);
//	float t = 0.0;
//	float tta = ta / resolution;
//	const float bb = b / resolution;
	vector<float> coord(dim3, 0.0), coordN(dim3, 0.0);
//	fPointVec = vector<float> ( (noOfPoints - 1) *  dim3 * resolution, 0.0 );
//	for ( unsigned int i = 0; i < noOfPoints - 1; ++i ) {
//		// add the first point
//		for ( unsigned int s = 0; s < dim3; s++)
//			coord[s] = fHelixVec[ s + dim3 * i ];
//		coordN = fRot.times(coord);
//		for(unsigned int s = 0; s<dim3; ++s)
//			coordN[s] += gCenter[s];
//		fPointVec.insert( fPointVec.end(), coordN.begin(), coordN.end() );
//		// add the fitted intermediate points for that interval
//		for ( unsigned int j = 1; j < resolution; ++j ) {
//			t = j * tta + i * ta;
//			coord[0] = a *  cos(t) - fCenter[0];
//			coord[1] = a *  sin(t)  - fCenter[1];
//			coord[2] = i * b + j * bb - fCenter[2];
//			coordN = fRot.times(coord);
//			for(unsigned int s = 0; s<dim3; ++s)
//				coordN[s] += gCenter[s];
//			fPointVec.insert( fPointVec.end(), coordN.begin(), coordN.end() );
//		}
//	}
	// compute $fOrigin$
	coord[0] = 0.0 - fCenter[0];
	coord[1] = 0.0 - fCenter[1];
	coord[2] = 0.0 - fCenter[2];
	coordN = fRot.times(coord);
	for(unsigned int s = 0; s<dim3; ++s)
		fOrigin[s] =  coordN[s] + gCenter[s];

	// compute $fNorm$
	vector<float> unitZ(dim3, 0.0);
	unitZ[2] = 1.0;
	fNorm = fRot.times(unitZ);

//	if(debug)  {
//		printf("a=%7.3f, b=%7.3f, ta=%7.3f\n", a, b, ta * rad2Deg );
////		fRot.print();
//		printf("RMSD=%8.4f\n", sqrt ( minRMS / noOfPoints)) ;
////		// back-computation
////		for (unsigned int i = 0; i < noOfPoints; ++i) {
////			for(unsigned int s = 0; s<dim3; ++s)  {
////				printf("%10.6f ", helixVec[dim3 * i +s ] );
////				coord[s] = fHelixVec[dim3 * i +s ];
////			}
////			cout<<endl;
////			coordN = fRot.times(coord);
////			for(unsigned int s = 0; s<dim3; ++s) {
////				sHelixVec[dim3 * i + s ] = coordN[s] + gCenter[s];
////				printf("%10.6f ", sHelixVec[dim3 * i +s ] );
////			}
////			cout<<endl;
////		}
////		cout<<endl;
//	}
	return true;
}

/**
 *
 * GOOD for protein helix
 * 	const float b0 = 1.10, deltaB = 0.02; // nb = 30, so  1.2 -- 1.2+0.6
	const float a0 = 1.95, deltaA = 0.02; // 1.728-- 2.728
	const size_t na = 40;
	const size_t nb = 40;
	const size_t resolution = 20;
 */
inline void generateHelixPolyline( const vector<string>& idVec,
		const vector<float>& ncacoCoordVec,
		const float b0, const float deltaB, const size_t nb, // parameters for helix pitch
		const float a0, const float deltaA, const size_t na, // parameters for helix circle
		const size_t resolution,
		vector<float>& polyline,
		vector<float>& circleCenterVec,  // the center of the helix using the quadruple
		vector<float>& normalVec,   // the normal of for each quadruple
		vector<float>& fittedParVec,  // the helical parameters for each quadruple
		bool fitplane = false, bool debug = false ) {

	const size_t noOfAtom = (size_t) (ncacoCoordVec.size() / dim3 );
	if( noOfAtom < dim4 )
		return;
	vector<float> ghelixVec(dim12, 0.0);
	bool debugFit = false; // true;
	float ta = 0.0;
	//	if(debug)
	//		debugFit = true;

	//	double minRMS = 0.0;
	vector<float> fittedHelicalParam(dim4, 0.0);
	//  p12 means the points between CA1 and CA2,
	// the coordinate of CA1 itself is saved as the first 3 elements of p12Vec etc.
	// pVec for previously-generated and  ppVec for that generated two steps back
	vector<float> p12Vec, p23Vec, p34Vec, p12pVec,  p23pVec, p12ppVec;
	vector<float> fOrigin(dim3, 0.0), fNorm(dim3, 0.0);
	ghelixVec.assign( ncacoCoordVec.begin(), ncacoCoordVec.begin() + dim12);
	fit4Points2Helix ( ghelixVec, b0, deltaB, nb,  a0, deltaA, na, resolution,
			p12Vec, p23Vec, p34Vec, fOrigin, fNorm, fittedHelicalParam,
			debugFit );
	fittedParVec.insert(fittedParVec.end(), fittedHelicalParam.begin(),
			fittedHelicalParam.end() );
	normalVec.insert(normalVec.end(), fNorm.begin(), fNorm.end() );

	double residual = 1.0 / mediumEPS;
	bool onXY = false;
	vector<double> center(dim3, 0.0);
	vector<double> norm( dim3, 0.0 );
	vector<double> coef(dim4, 0.0);
	if(fitplane) {
		if( ! fit2Plane( ghelixVec, onXY, center, coef, residual ) )
			residual = 1.0 / mediumEPS;
//		else
//			residual = sqrt( residual / dim4 );
	}

	// fit to plane
	if(debug) {
		ta += fittedHelicalParam[2] ;
		if (fitplane) {
			printf("residID,   a       b        (a+b)   angle / ta       RMSD   residual2P \n");
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%7.2f %8.4f %8.3f\n",
					idVec.front().c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta * rad2Deg, fittedHelicalParam[3],
					sqrt( residual / dim4 ) );
		}else {
			printf("residID,   a       b      (a+b)   angle / ta      RMSD\n");
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f\n",
					idVec.front().c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3] );
		}
	}
	circleCenterVec.insert(circleCenterVec.end(), fOrigin.begin(), fOrigin.end() );
	// save $p12Vec$
	polyline.insert(polyline.end(), p12Vec.begin(), p12Vec.end() );

	// Only 4 CAs, save p23 and p34 and RETURN
	if ( noOfAtom == dim4 ) {
		polyline.insert(polyline.end(), p23Vec.begin(), p23Vec.end() );
		polyline.insert(polyline.end(), p34Vec.begin(), p34Vec.end() );
		polyline.insert(polyline.end(), ncacoCoordVec.end() - dim3, ncacoCoordVec.end());
		circleCenterVec.push_back( fOrigin[0] + 4.0*fNorm[0]);
		circleCenterVec.push_back( fOrigin[1] + 4.0*fNorm[1]);
		circleCenterVec.push_back( fOrigin[2] + 4.0*fNorm[2]);
//		for ( unsigned int j =1; j < idVec.size();  ++j )
//			cout<< idVec[j]<<endl;
		return;
	}

	// save the previously-generated
	p12pVec = p23Vec;
	p23pVec = p34Vec;
	// compute a helix using a new quadruple
	ghelixVec.assign( ncacoCoordVec.begin()+dim3, ncacoCoordVec.begin()+dim15);
	fit4Points2Helix ( ghelixVec, b0, deltaB, nb,  a0, deltaA, na, resolution,
			p12Vec, p23Vec, p34Vec, fOrigin, fNorm, fittedHelicalParam,
			debugFit );
	fittedParVec.insert(fittedParVec.end(),
			fittedHelicalParam.begin(), fittedHelicalParam.end() );
	normalVec.insert(normalVec.end(), fNorm.begin(), fNorm.end() );

	if(fitplane)
		if( ! fit2Plane( ghelixVec, onXY, center, coef,  residual ) )
			residual = 1.0 / mediumEPS;

	if(debug) {
		ta += fittedHelicalParam[2] ;
		if(ta > TwoPi )
			ta -= TwoPi;
		if (fitplane) {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f %8.3f\n",
					idVec[1].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3],
					sqrt( residual / dim4 ) );
		}else {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f\n",
					idVec[1].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3] );
		}
	}

	circleCenterVec.insert(circleCenterVec.end(), fOrigin.begin(), fOrigin.end() );
	// Compute  $p23Vec$ by averaging the current p23Vec and the saved p23Vec
	// the current p23Vec is $p12Vec$ and the previous one is $p12pVec$
	for ( unsigned int s = 0; s < p12Vec.size(); s++)
		p12Vec[s] = 0.5 * ( p12Vec[s] +p12pVec[s] ) ;
	// save $p23Vec$
	polyline.insert(polyline.end(), p12Vec.begin(), p12Vec.end() );

	// Only 5 CAs, save p34 and p45 and RETURN
	if ( noOfAtom == dim5 ) {
		// Compute and save $p34Vec$
		for ( unsigned int s = 0; s < p23Vec.size(); s++)
			p23Vec[s] = 0.5 * ( p23Vec[s] +p23pVec[s] ) ;
		polyline.insert(polyline.end(), p23Vec.begin(), p23Vec.end() );
		// save $p45Vec$
		polyline.insert(polyline.end(), p34Vec.begin(), p34Vec.end() );
		// save the original last CA as the end point
		polyline.insert(polyline.end(), ncacoCoordVec.end() - dim3, ncacoCoordVec.end());
		circleCenterVec.push_back( fOrigin[0] + 4.0*fNorm[0]);
		circleCenterVec.push_back( fOrigin[1] + 4.0*fNorm[1]);
		circleCenterVec.push_back( fOrigin[2] + 4.0*fNorm[2]);
		return;
	}
	// save the previously-generated and  previously-previously-generated $p34$
	p12ppVec = p23pVec; // it is the $p34$ generated in the 1st quadruple
	p12pVec = p23Vec;     // it is the $p34$ generated in the 2nd quadruple
	p23pVec = p34Vec;    // it is the $p45$ generated in the 2nd quadruple

	ghelixVec.assign(ncacoCoordVec.begin()+dim6, ncacoCoordVec.begin()+dim6+dim12);
	fit4Points2Helix ( ghelixVec, b0, deltaB, nb,  a0, deltaA, na, resolution,
			p12Vec, p23Vec, p34Vec, fOrigin, fNorm, fittedHelicalParam, debugFit );
	fittedParVec.insert(fittedParVec.end(),
			fittedHelicalParam.begin(), fittedHelicalParam.end() );
	normalVec.insert(normalVec.end(), fNorm.begin(), fNorm.end() );

	if(fitplane)
		if( ! fit2Plane( ghelixVec, onXY, center, coef, residual ) )
			residual = 1.0 / mediumEPS;

	if(debug) {
		ta += fittedHelicalParam[2] ;
		if(ta > TwoPi )
			ta -= TwoPi;
		if (fitplane) {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f %8.3f\n",
					idVec[2].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3],
					sqrt( residual / dim4 ) );
		}else {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f\n",
					idVec[2].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta* rad2Deg, fittedHelicalParam[3] );
		}
	}

	circleCenterVec.insert(circleCenterVec.end(), fOrigin.begin(), fOrigin.end() );
	// Compute and save $p34Vec$ by averaging over the  $p34$ computed
	// in the 1st, 2nd and the current ( the 3rd) quadruples
	for ( unsigned int s = 0; s < p12Vec.size(); s++)
		p12Vec[s] = ( p12Vec[s] +p12pVec[s] +p12ppVec[s] ) / 3.0 ;
	polyline.insert(polyline.end(), p12Vec.begin(), p12Vec.end() );

	// Only 6 CAs, save $p45Vec$ and $p56Vec$ and RETURN
	if ( noOfAtom == dim6 ) {
		// Compute and save $p45Vec$
		for ( unsigned int s = 0; s < p23Vec.size(); s++)
			p23Vec[s] = 0.5 * ( p23Vec[s] +p23pVec[s] ) ;
		polyline.insert(polyline.end(), p23Vec.begin(), p23Vec.end() );
		// save $p56Vec$
		polyline.insert(polyline.end(), p34Vec.begin(), p34Vec.end() );
		polyline.insert(polyline.end(), ncacoCoordVec.end() - dim3, ncacoCoordVec.end());
		circleCenterVec.push_back( fOrigin[0] + 4.0*fNorm[0]);
		circleCenterVec.push_back( fOrigin[1] + 4.0*fNorm[1]);
		circleCenterVec.push_back( fOrigin[2] + 4.0*fNorm[2]);
		return;
	}

	// save the previously-generated and  previously-previously-generated $p45$
	p12ppVec = p23pVec; // it is the $p45$ generated in the 2nd quadruple
	p12pVec = p23Vec;   // it is the $p45$ generated in the 3rd quadruple
	p23pVec = p34Vec;    // it is the $p56$ generated in the 3rd quadruple

	// Get into the loop only when there are > 7 CAs
	// So what happens if we have only 7 CAs, no loop, go directly to the next
	unsigned int i = 0;
	for ( i = 3; i < noOfAtom - dim4; ++i ) {
		ghelixVec.assign( ncacoCoordVec.begin() + dim3 * i,
				ncacoCoordVec.begin() + dim3 * i + dim12 );
		fit4Points2Helix ( ghelixVec, b0, deltaB, nb,  a0, deltaA, na, resolution,
				p12Vec, p23Vec, p34Vec, fOrigin, fNorm, fittedHelicalParam, debugFit );

		fittedParVec.insert(fittedParVec.end(),
				fittedHelicalParam.begin(),fittedHelicalParam.end() );

		normalVec.insert(normalVec.end(), fNorm.begin(), fNorm.end() );

		if(fitplane)
			if( ! fit2Plane( ghelixVec, onXY, center, coef, residual ) )
				residual = 1.0 / mediumEPS;
		if(debug) {
			ta += fittedHelicalParam[2] ;
			if(ta > TwoPi )
				ta -= TwoPi;
			if (fitplane) {
				printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f %8.3f\n",
						idVec[i].c_str(),
						fittedHelicalParam[0], fittedHelicalParam[1],
						(fittedHelicalParam[0]+fittedHelicalParam[1]),
						fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3],
						sqrt( residual / dim4 ) );
			}else {
				printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.3f\n",
						idVec[i].c_str(),
						fittedHelicalParam[0], fittedHelicalParam[1],
						(fittedHelicalParam[0]+fittedHelicalParam[1]),
						fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3] );
			}
		}

		circleCenterVec.insert(circleCenterVec.end(), fOrigin.begin(), fOrigin.end() );
		for ( unsigned int s = 0; s < p12Vec.size(); s++)
			p12Vec[s] = ( p12Vec[s] +p12pVec[s] +p12ppVec[s] ) / 3.0 ;
		polyline.insert(polyline.end(),  p12Vec.begin(), p12Vec.end() );
		p12ppVec = p23pVec;
		p12pVec = p23Vec;
		p23pVec = p34Vec;
	}

	// the last two:
	//	const size_t last12begin = ncacoCoordVec.size() - dim12;
	ghelixVec.assign( ncacoCoordVec.end() - dim12, ncacoCoordVec.end() );
	fit4Points2Helix ( ghelixVec, b0, deltaB, nb,  a0, deltaA, na, resolution,
			p12Vec, p23Vec, p34Vec, fOrigin, fNorm, fittedHelicalParam, debugFit );
	fittedParVec.insert( fittedParVec.end(),
			fittedHelicalParam.begin(), fittedHelicalParam.end() );
	normalVec.insert(normalVec.end(), fNorm.begin(), fNorm.end() );

	if(fitplane)
		if( ! fit2Plane( ghelixVec, onXY, center, coef, residual ) )
			residual = 1.0 / mediumEPS;

	if(debug) {
		ta += fittedHelicalParam[2] ;
		if( ta > TwoPi )
			ta -= TwoPi;
		if (fitplane) {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f %8.3f\n",
					idVec[i].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta*rad2Deg, fittedHelicalParam[3],
					sqrt( residual / dim4 ) );
		}else {
			printf("%5s%8.3f %8.3f %8.3f %8.3f/%8.3f %8.4f\n",
					idVec[i].c_str(),
					fittedHelicalParam[0], fittedHelicalParam[1],
					(fittedHelicalParam[0]+fittedHelicalParam[1]),
					fittedHelicalParam[2] * rad2Deg, ta * rad2Deg, fittedHelicalParam[3] );
		}
	}

	circleCenterVec.insert(circleCenterVec.end(), fOrigin.begin(), fOrigin.end() );
	// if only 7 CAs, this is p45
	for ( unsigned int s = 0; s < p12Vec.size(); s++)
		p12Vec[s] = ( p12Vec[s] +p12pVec[s] +p12ppVec[s] ) / 3.0 ;
	polyline.insert(polyline.end(),  p12Vec.begin(), p12Vec.end() );
	// if only 7 CAs, this is p56
	for ( unsigned int s = 0; s < p23Vec.size(); s++)
		p23Vec[s]  = 0.5 * ( p23Vec[s] +p23pVec[s] ) ;
	polyline.insert(polyline.end(), p23Vec.begin(), p23Vec .end() );
	// if only 7 CAs, this is p67
	polyline.insert(polyline.end(), p34Vec.begin(), p34Vec.end() );
	polyline.insert(polyline.end(), ncacoCoordVec.end() - dim3, ncacoCoordVec.end());
	circleCenterVec.push_back( fOrigin[0] + 4.0*fNorm[0]);
	circleCenterVec.push_back( fOrigin[1] + 4.0*fNorm[1]);
	circleCenterVec.push_back( fOrigin[2] + 4.0*fNorm[2]);

//	for ( unsigned int j = i+1; j < idVec.size();  ++j )
//		cout<< idVec[j]<<endl;

	return;
}

//
//inline double areaOfPolygon( double *v, const size_t n){
//	double area = 0;
//	int   i, j, k;     // indices
//
//	for (i=1, j=2, k=0; i<=n; i++, j++, k++) {
//		area += (v[i].x * (V[j].y - V[k].y);
//	}
//	return area / 2.0;
//}
//
///**
// * Codes for computing areas in 2D and 3D
// */
//// a Point (or vector) is defined by its coordinates
//typedef struct {int x, y, z;} Point;    // exclude z for 2D
//// a Triangle is given by three points: Point V0, V1, V2
//// a Polygon is given by:
////        int n = number of vertex points
////        Point* V[] = an array of points with V[n]=V[0], V[n+1]=V[1]
//
//// Note: for efficiency low-level functions are declared to be inline.
//
//// isLeft(): tests if a point is Left|On|Right of an infinite line.
////    Input:  three points P0, P1, and P2
////    Return: >0 for P2 left of the line through P0 and P1
////            =0 for P2 on the line
////            <0 for P2 right of the line
//inline int
//isLeft( Point P0, Point P1, Point P2 )
//{
//    return ( (P1.x - P0.x) * (P2.y - P0.y)
//            - (P2.x - P0.x) * (P1.y - P0.y) );
//}
////===================================================================
//
//// orientation2D_Triangle(): test the orientation of a triangle
////    Input:  three vertex points V0, V1, V2
////    Return: >0 for counterclockwise
////            =0 for none (degenerate)
////            <0 for clockwise
//inline int
//orientation2D_Triangle( Point V0, Point V1, Point V2 )
//{
//    return isLeft(V0, V1, V2);
//}
////===================================================================
//
//// area2D_Triangle(): compute the area of a triangle
////    Input:  three vertex points V0, V1, V2
////    Return: the (float) area of T
//inline float
//area2D_Triangle( Point V0, Point V1, Point V2 )
//{
//    return (float)isLeft(V0, V1, V2) / 2.0;
//}

////===================================================================
//
//// orientation2D_Polygon(): tests the orientation of a simple polygon
////    Input:  int n = the number of vertices in the polygon
////            Point* V = an array of n+1 vertices with V[n]=V[0]
////    Return: >0 for counterclockwise
////            =0 for none (degenerate)
////            <0 for clockwise
////    Note: this algorithm is faster than computing the signed area.
//int
//orientation2D_Polygon( int n, Point* V )
//{
//    // first find rightmost lowest vertex of the polygon
//    int rmin = 0;
//    int xmin = V[0].x;
//    int ymin = V[0].y;
//
//    for (int i=1; i<n; i++) {
//        if (V[i].y > ymin)
//            continue;
//        if (V[i].y == ymin) {    // just as low
//            if (V[i].x < xmin)   // and to left
//                continue;
//        }
//        rmin = i;          // a new rightmost lowest vertex
//        xmin = V[i].x;
//        ymin = V[i].y;
//    }
//
//    // test orientation at this rmin vertex
//    // ccw <=> the edge leaving is left of the entering edge
//    if (rmin == 0)
//        return isLeft( V[n-1], V[0], V[1] );
//    else
//        return isLeft( V[rmin-1], V[rmin], V[rmin+1] );
//}
//===================================================================

//static double norm(double*, int );
//static double det(double **, int);
//static double* interVector(double *p1, double *p2, int d);
//static int sgn(double x);

/**
 * Three functions to compute a xyz coordinate from a grid point
 */
inline void grid2Coord( const vector<size_t>& gridXYZ, const vector<double>& bioMinMax,
		const size_t gIndex, const double bioDelta, double &x, double &y, double &z ){
	size_t indexX = gIndex % gridXYZ[0];
	size_t indexT = ( gIndex - indexX ) / gridXYZ[0] ;
	size_t indexY = indexT % gridXYZ[1] ;
	size_t indexZ = ( indexT - indexY ) / gridXYZ[1] ;
	x = (indexX + 0.50) * bioDelta + bioMinMax[0]  ;
	y = ( indexY + 0.50) * bioDelta + bioMinMax[1] ;
	z = ( indexZ + 0.50) * bioDelta + bioMinMax[2] ;
}

inline void grid2Coord( const vector<size_t>& gridXYZ, const vector<double>& bioMinMax,
		const size_t gIndex, const double bioDelta, double *coords ){
	size_t indexX = gIndex % gridXYZ[0];
	size_t indexT = ( gIndex - indexX ) / gridXYZ[0] ;
	size_t indexY = indexT % gridXYZ[1] ;
	size_t indexZ = ( indexT - indexY ) / gridXYZ[1] ;
	coords[0] = (indexX+ 0.50)  * bioDelta + bioMinMax[0] ;
	coords[1] = (indexY+ 0.50)  * bioDelta + bioMinMax[1] ;
	coords[2] = (indexZ+ 0.50)  * bioDelta + bioMinMax[2] ;
}

inline bool grid2Coord( const vector<size_t>& gridXYZ, const vector<double>& bioMinMax,
		const matrix& mat, const size_t gIndex, const double bioDelta, double *coords ){
	size_t indexX = gIndex % gridXYZ[0];
	size_t indexT = ( gIndex - indexX ) / gridXYZ[0] ;
	size_t indexY = indexT % gridXYZ[1] ;
	size_t indexZ = ( indexT - indexY ) / gridXYZ[1] ;
	const size_t dim3 = 3;
	vector<double> rCoord(dim3, 0.0);
	coords[0] = (indexX + 0.50) * bioDelta + bioMinMax[0] ;
	coords[1] = (indexY + 0.50) * bioDelta + bioMinMax[1] ;
	coords[2] = (indexZ + 0.50) * bioDelta + bioMinMax[2] ;
	if ( mat.times( coords, dim3, rCoord ) ){
		for (unsigned int i = 0; i < dim3; ++i)
			coords[i] = rCoord[i];
		return true;
	}else{
		cerr<<"could not compute the coordinate\n";
		return false;
	}
}

inline int grid2PolarCoord( const vector<size_t>& gridXYZ, const vector<double>& bioMinMax,
		const size_t gIndex, const double bioDelta, double *polarCoords ){
	size_t indexX = gIndex % gridXYZ[0];
	size_t indexT = ( gIndex - indexX ) / gridXYZ[0] ;
	size_t indexY = indexT % gridXYZ[1] ;
	size_t indexZ = ( indexT - indexY ) / gridXYZ[1] ;
	const size_t dim3 = 3;
	double coords[dim3];
	coords[0] = (indexX+ 0.50)  * bioDelta + bioMinMax[0] ;
	coords[1] = (indexY+ 0.50)  * bioDelta + bioMinMax[1] ;
	coords[2] = (indexZ+ 0.50)  * bioDelta + bioMinMax[2] ;
	double dis = 0.0;
	for (unsigned int i = 0; i < dim3; i++ )
		dis += coords[i]  * coords[i] ;

	if ( dis  < Tiny  )
		return 0;
	const double sinT = fabs( 1.0 - coords[2] * coords[2] / dis);
	if ( sinT < Tiny  ){
		polarCoords[1]  = 0.0;
		return 1;
	}
	polarCoords[0] = sqrt( dis );
	polarCoords[1] = acos( coords[2] / polarCoords[0] ); // $theta$
	// atan2 return [ -pi, +pi ], add $\pi$ if want to change to the range between $[0, 2\pi]$
	polarCoords[2] = atan2( coords[1], coords[0] );
	return 2;
}
//
//inline void extropolation(double *p1, double *dirCos,
//		const double delta, ){
//
//}

/**
 * compute the rotation matrix to a plane defined by the following three atoms
 * 	( atomOfChi<---Top--->-center ) with
 * +Y axis along Top-->center direction,
 * +Z axis along Top-->atomOfChi direction with
 *      an angle $\theta$ from the Top-->atomOfChi vector
 *
 *  vector $v1N$ represents rotated Top--->center vector that lies on +Y
 *  vector $v2N$ represents rotated Top--->atomOfChi vector that lies
 *  with a $\theta$ angle from +Z
 *  while the angle formed by $v1N$ and $v2N$ > 90

 ( $v2N$ )
    atomOfChi (prev)       +Z
         \                             /\
            \                           |
               \       $\theta$    |
                  \                     |
                      \                 |
                          \             |
                        	  \         |
                       	   	      \     |
                                      Top -------------------->center ( +Y ) ( $v1N$ )
 *
 * Possible numerical error due to the limited accuracy of pdb coordinates and since it is stored as
 * float, it will lead to large error in the computation of _ rotmat_
 */
inline bool rot2YZplane( double *const atomOfChi, double *const top, double *const center,
		vector<double>& mat, bool transpose = false, bool debug = false) {
//	debug = true;
	double cosTheta = 0.0;
	vector<double> v1 = n1ToN2VecD( top, center, dim3 );
	vector<double> v2 = n1ToN2VecD( top, atomOfChi, dim3 );

	double norm = 0.0, sinT = 0.0, cosT = 0.0,  length = 0.0;
	double dirCos[dim3]={0.0, 0.0, 0.0};
	double rotDir[dim3] = {0.0, 0.0, 0.0};
	vector<double> v1ToYaxisMat ( dim9 );
	v1ToYaxisMat[0] = v1ToYaxisMat[4] = v1ToYaxisMat[8] = 1.0;  // identity matrix
	vector<double> v2N = v2;

	// rotate $v1$ to Y-axis
	if( fabs(v1[0]) > smallEPS || fabs(v1[2]) > smallEPS ) {
		norm = sqrt( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] );
		dirCos[0] = v1[0] / norm;
		dirCos[1] = v1[1] / norm;
		dirCos[2] = v1[2] / norm;
		// rotDir is perpendicular to the plane formed by $v1$ and
		// the Y-axis of the original frame because $rotDir[1] == 0.0$
		rotDir [0] = -dirCos[2];
		rotDir[2]  =  dirCos[0];
		// the angle needed to rotate $v1$ to Y-axis
		cosT = dirCos[1];
		sinT =  sqrt( 1.0 - cosT * cosT);;
//		if( cosT < 0 )
//			sinT = - sinT;
		length = 1.0  / sqrt( rotDir[0] * rotDir[0] + rotDir[2] * rotDir[2]);
		rotDir[0] *= length;
		rotDir[2] *= length;
		rotationAlongAxisWithAngle( sinT, cosT, rotDir, v1ToYaxisMat );
		if(debug) {
			vector<double> v11(dim3, 0.0);
			mTimesV( v1ToYaxisMat, v1, v11);
			printArray( v11);
		}
		mTimesV( v1ToYaxisMat, v2, v2N );
	}

	// HERE is a bug
	if(debug) {
		cout<<"V2 After the rotation !\n";
		printArray( v2N);
	}
	vector<double> v2ToYZplaneMat (dim9, 0.0);
	vector<double> mm (dim9, 0.0);
	v2ToYZplaneMat[0] = v2ToYZplaneMat[4] = v2ToYZplaneMat[8] = 1.0;  // identity matrix
	if( fabs(v2N[0] ) > smallEPS ){
		norm =  sqrt( v2N[0] * v2N[0] + v2N[2] * v2N[2]);
		// this is the angle needed to rotate $v2N$ to YZ plane
		cosT =   v2N[2] / norm; // cos(theta);
		sinT =  sqrt( 1.0 - cosT * cosT);;
		// rotate along -Y axis
		rotDir[0] = 0.0;
		rotDir[2] = 0.0;
		if( v2N[0] > 0.0)
			rotDir[1] = -1.0;
		else
			rotDir[1] = 1.0;
		rotationAlongAxisWithAngle( sinT, cosT, rotDir, v2ToYZplaneMat );

		mTimesN ( dim3, v2ToYZplaneMat, v1ToYaxisMat, mm );
		// check whether it has rotated correctly
		vector<double> v22( dim3, 0.0 );
		mTimesV( mm, v2, v22 );
		if(debug) {
			cout<<"V2 After the 2nd rotation !\n";
			printArray( v22);
		}
//		if( fabs(v22[0]) > smallEPS ){ // rotate with respect to +Y axis
//			if(debug)
//				cout<<"Rotate along +y\n";
//			rotDir[0] = 0.0;
//			rotDir[1] = -1.0;
//			rotDir[2] = 0.0;
//			rotationAlongAxisWithAngle(sinT, cosT, rotDir, v2ToYZplaneMat );
//			mTimesN( dim3, v2ToYZplaneMat, v1ToYaxisMat, mm);
////			mTimesV( mm, v2, v22 );
//			mTimesV( v2ToYZplaneMat, v2N, v22 );
//		}
//		if(debug) {
//			cout<<"V2 After the 22nd  rotation !\n";
//			printArray( v22);
//		}
	}

	if(transpose){
		for (unsigned int  s = 0; s < dim3; ++s)
			for (unsigned int t = 0; t < dim3; ++t )
				mat [ dim3 * s + t ] = mm[ dim3 * t + s ];
	}else
		mat = mm;

//	cout<<"V2 VBefore  the rotation !\n";
//	printArray( v2);
//	cout<<"V2 After the rotation ! !\n";
//	vector<double> v2Coord(dim3, 0.0);
//	mTimesV( mat, v2, v2Coord);
//	printArray( v2Coord);
//
//	cout<<"V1 VBefore  the rotation !\n";
//	printArray( v1);
//	cout<<"V1 After the rotation !\n";
//	vector<double> v1Coord(dim3, 0.0);
//	mTimesV( mat, v1, v1Coord);
//	printArray( v1Coord);
//
//	vector<double> mm2 (dim3 * dim3);
//
//	for ( unsigned int s = 0; s < dim3; ++s)
//		for (unsigned int t = 0; t < dim3; ++t )
//			mm2 [ dim3 * s + t ] = mat [ dim3 * t + s ];
//
//	vector<double> v11(dim3, 0.0);
//	mTimesV( mm2, v1Coord, v11);
//	cout<<"V1 recovered !\n";
//	printArray( v11);
//
//	vector<double> v22(dim3, 0.0);
//	mTimesV( mm2, v2Coord, v22);
//	cout<<"V2 recovered !\n";
//	printArray( v22);

	return true;
//	double rotDir [] = {0.0, -dirCos[2], dirCos[1] };
//		double cosT = dirCos[0];
//		double sinT =  sqrt( 1.0 - cosT * cosT);;
//		if( cosT < 0 )
//			sinT = - sinT;
//		double length = 1.0  / sqrt( rotDir[1] * rotDir[1] + rotDir[2] * rotDir[2]);
//		rotDir [1] *= length;
//		rotDir [2] *= length;
//		vector<double> mat2 (sp2::dim3 * sp2::dim3, 0.0);
//		rotationAlongAxisWithAngle( -sinT, cosT, rotDir, mat2 );
//
//		vector<double> v1N ( sp2::dim3, 0.0);
//		vector<double> v2N ( sp2::dim3, 0.0);
//		vector<double> v1Coord(sp2::dim3, 0.0);
//		mTimesV(mat2, v1, v1Coord);
//		cout<<"After the rotation !\n";
//		printArray( v1Coord);

//	cout<<"\n compute the rotation matrix\n";

//	if( interAngleD( v1, v2, sp2::dim3, cosTheta ) ){
//		printf("cosTheta=%8.3f\n", (acos(cosTheta )*RadToDegrees  ) );
//
//		// $v1N$ and $v2N$
//		v1N[1] = norm( v1, sp2::dim3 ) ; // the length of $v1$
//		double normOfV2 = norm( v2, sp2::dim3 ) ;
//		v2N[1] = normOfV2 * cosTheta ;
//		v2N[2] = normOfV2 * sqrt ( 1.0 - cosTheta * cosTheta ); //sinTheta ;
//		if( rotmat( v1, v2, v1N, v2N, sp2::dim3, mat ) ) {
//			return true;
//		}
//		cerr<<"can not compute rotmat\n";
//		return false;
//	}
//
//	cerr<<"can not compute rotmat2\n";
//	return false;
}

inline bool rot2YZplane( const vector<float>&atomOfChi, const vector<float>& top,
		const vector<float>& center, vector<double>& mat, bool transpose = false) {
	const size_t dim3 = 3;

	vector<double> v1(dim3, 0.0);
	vector<double> v2 (dim3, 0.0);
	for (unsigned int i = 0; i <dim3; ++i) {
		v1[i] =  center[i] -  top[i];
		v2[i] =  atomOfChi[i] -  top[i];
	}

	double norm = 0.0, sinT = 0.0, cosT = 0.0,  length = 0.0;
	double dirCos[dim3];
	double rotDir[dim3] = {0.0, 0.0, 0.0};
	vector<double> v1ToYaxisMat (dim3 * dim3);
	v1ToYaxisMat[0] = v1ToYaxisMat[4] = v1ToYaxisMat[8] = 1.0;  // identity matrix
	vector<double> v2N = v2;
	if( fabs(v1[0]) > smallEPS || fabs(v1[2]) > smallEPS ) {
		norm = sqrt( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] );
		dirCos[0] = v1[0] / norm;
		dirCos[1] = v1[1] / norm;
		dirCos[2] = v1[2] / norm;
		rotDir[0] = -dirCos[2];
		rotDir[2] =  dirCos[0];
		cosT = dirCos[1];
		sinT =  sqrt( 1.0 - cosT * cosT );;
//		if( cosT < 0 )
//			sinT = - sinT;
		length = 1.0  / sqrt( rotDir[0] * rotDir[0] + rotDir[2] * rotDir[2]);
		rotDir[0] *= length;
		rotDir[2] *= length;
		rotationAlongAxisWithAngle( sinT, cosT, rotDir, v1ToYaxisMat );

		vector<double> v11(dim3, 0.0);
		mTimesV( v1ToYaxisMat, v1, v11);
//		printArray( v11);
		mTimesV( v1ToYaxisMat, v2, v2N );
	}
//	cout<<"V2 After the rotation !\n";
//	printArray( v2N);

	vector<double> v2ToYZplaneMat (dim3 * dim3);
	vector<double> mm (dim3 * dim3);
	v2ToYZplaneMat[0] = v2ToYZplaneMat[4] = v2ToYZplaneMat[8] = 1.0;  // identity matrix
	if( fabs ( v2N[0] ) > smallEPS ){
		//		norm = sqrt( v2N[0] * v2N[0] + v2N[1] * v2N[1] + v2N[2] * v2N[2] );
		// make a dihedral angle
//		double p1[]= {v2N[0], v2N[1], v2N[2]};
//		double p2[] = {0.0, 0.0, 0.0};
//		double p3[] = {1.0, 0.0, 0.0};
//		double p4[] = {1.0, v2N[1], v2N[2]};
//
//		geomBasic *gb = new geomBasic();
//		double theta = 0.0;
//		if (! gb->signedDihedralAngle(p1, p2, p3, p4, dim3, theta) ){
//			cerr<<"Plane: can not compute $\\chi$ \n";
//			delete gb;
//			return false;
//		}
		norm =  sqrt( v2N[0] * v2N[0] + v2N[2] * v2N[2] );
		cosT =   v2N[2] / norm; //sin(theta);
		sinT =  sqrt( 1.0 - cosT * cosT);;
//		if( cosT < 0 )
//			sinT = - sinT;
		rotDir [0] = 0.0;
		rotDir [1] = -1.0;
		rotDir[2]  = 0.0;
		rotationAlongAxisWithAngle( sinT, cosT, rotDir, v2ToYZplaneMat );

		mTimesN( dim3, v2ToYZplaneMat, v1ToYaxisMat, mm);
		vector<double> v22(dim3, 0.0);
		mTimesV( mm, v2, v22);
		if( fabs( v22[0] ) > smallEPS ) {
			rotationAlongAxisWithAngle ( -sinT, cosT, rotDir, v2ToYZplaneMat );
			mTimesN( dim3, v2ToYZplaneMat, v1ToYaxisMat, mm );
		}
//		printArray( v22);
	}

	if(transpose){
		for (unsigned int  s = 0; s < dim3; ++s)
			for (unsigned int t = 0; t < dim3; ++t )
				mat [ dim3 * s + t ] = mm[ dim3 * t + s ];
	}else
		mat = mm;

//	cout<<"V2 VBefore  the rotation !\n";
//	printArray( v2);
//	cout<<"V2 After the rotation ! !\n";
//	vector<double> v2Coord(dim3, 0.0);
//	mTimesV( mat, v2, v2Coord);
//	printArray( v2Coord);
//
//	cout<<"V1 VBefore  the rotation !\n";
//	printArray( v1);
//	cout<<"V1 After the rotation !\n";
//	vector<double> v1Coord(dim3, 0.0);
//	mTimesV( mat, v1, v1Coord);
//	printArray( v1Coord);
//
//	vector<double> mm2 (dim3 * dim3);
//
//	for ( unsigned int s = 0; s < dim3; ++s)
//		for (unsigned int t = 0; t < dim3; ++t )
//			mm2 [ dim3 * s + t ] = mat [ dim3 * t + s ];
//
//	vector<double> v11(dim3, 0.0);
//	mTimesV( mm2, v1Coord, v11);
//	cout<<"V1 recovered !\n";
//	printArray( v11);
//
//	vector<double> v22(dim3, 0.0);
//	mTimesV( mm2, v2Coord, v22);
//	cout<<"V2 recovered !\n";
//	printArray( v22);

	return true;
//	double rotDir [] = {0.0, -dirCos[2], dirCos[1] };
//		double cosT = dirCos[0];
//		double sinT =  sqrt( 1.0 - cosT * cosT);;
//		if( cosT < 0 )
//			sinT = - sinT;
//		double length = 1.0  / sqrt( rotDir[1] * rotDir[1] + rotDir[2] * rotDir[2]);
//		rotDir [1] *= length;
//		rotDir [2] *= length;
//		vector<double> mat2 (sp2::dim3 * sp2::dim3, 0.0);
//		rotationAlongAxisWithAngle( -sinT, cosT, rotDir, mat2 );
//
//		vector<double> v1N ( sp2::dim3, 0.0);
//		vector<double> v2N ( sp2::dim3, 0.0);
//		vector<double> v1Coord(sp2::dim3, 0.0);
//		mTimesV(mat2, v1, v1Coord);
//		cout<<"After the rotation !\n";
//		printArray( v1Coord);

//	cout<<"\n compute the rotation matrix\n";

//	if( interAngleD( v1, v2, sp2::dim3, cosTheta ) ){
//		printf("cosTheta=%8.3f\n", (acos(cosTheta )*RadToDegrees  ) );
//
//		// $v1N$ and $v2N$
//		v1N[1] = norm( v1, sp2::dim3 ) ; // the length of $v1$
//		double normOfV2 = norm( v2, sp2::dim3 ) ;
//		v2N[1] = normOfV2 * cosTheta ;
//		v2N[2] = normOfV2 * sqrt ( 1.0 - cosTheta * cosTheta ); //sinTheta ;
//		if( rotmat( v1, v2, v1N, v2N, sp2::dim3, mat ) ) {
//			return true;
//		}
//		cerr<<"can not compute rotmat\n";
//		return false;
//	}
//
//	cerr<<"can not compute rotmat2\n";
//	return false;
}



#endif /*GEOMBASIC_H_*/




















