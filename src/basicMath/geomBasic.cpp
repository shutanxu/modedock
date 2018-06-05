#include <cmath>
#include <iostream>

#include "matrix.h"
#include "geomBasic.h"
//#include "../constant.h"
using namespace std;

geomBasic::geomBasic(){
	p1 = p2 = p3 = p4 = 0;
}

geomBasic::~geomBasic(){
	delete [] p1;
	delete [] p2;
	delete [] p3;
	delete [] p4;
}

/**
 The three-point form of the plane:
 the same equation can be used to detect whether
 the four points  coplanar, (x,y, z; x1, y1, z1; x2, y2, z2; x3, y3, z3;}
 |x     y    z     1  |
 |x1  y1   z1   1  |   =   |x-x1  y-y1  z-z1    |
 |x2   y2  z2   1  |         |x2-x1  y2-y1  z2-z1 |  = 0
 |x3    y3   z3 1  |         |x3-x1  y3-y1  z3-z1|
 The final equation is
  	   ax+by+cz+d = 0
 */

double* geomBasic::coefOfPlaneEqu(double *p1, double *p2, double *p3, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		exit(0);
	}
	double a2 = p2[0] - p1[0];
	double b2 = p2[1] - p1[1];
	double c2 = p2[2] - p1[2];
	double a3 = p3[0] - p1[0];
	double b3 = p3[1] - p1[1];
	double c3 = p3[2] - p1[2];
	double *coef = new double[4];
	coef[0] = b2 * c3 - c3 - b2;
	coef[1] = a3 *c2 - c3 * a2;
	coef[2] = a2 * b3 - a3 * b2;
	coef[3] = -( p1[0] * coef[0] +p1[1] * coef[1] +p1[2] * coef[2] );
	return coef;
}

double* geomBasic::crossProd(double *p1, double *p2, const size_t d){
	if ( d!=3){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		exit(0);
	}
	double *c = new double[d];
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	return c;
}

bool geomBasic::triangleArea(double *v1, double *v2, double *v3, double *area, const size_t d){
	if ( d != 3 ){
		cerr<<"The function only tests 3 points in 3D: but the dimension d ="<<d<<endl;
		return false;
	}
	double p1 [3] = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
	double p2 [3] = {v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]};
	double c[3];
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	*area = 0.50 * ( c[0] * c[0] + c[1] * c[1] + c[2] * c[2] );
	return true;
}


/**
 * The following method does return the NOT, correctly,
 *  positive and negative signs for the dihedral angles
 * since it is is trick to determine the sign here.
 * The dihedral angle between the two planes
 \begin{eqnarray*}
a1 x + b1 y + c1 z + d1	&=&	0	\\
a2 x + b2 y + c2 z + d2	&=&	0	\\
\end{eqnarray*}
with the respective  normal vectors $n1=(a1, b1, c1) , n2=(a2, b2, c2)$
is simply given via the dot product of the two normals,
$\cos{\theta}	=	n1 \dot n2$
* the computation is efficient but the angle returned may NOT conform to the
* standard, then the question is how?,
*  for the computation of energy, I do not think it matters since the different is
* either the returned value or the minus of the returned value
*/

double geomBasic::dihedralAngle(double *p1, double *p2, double *p3, double *p4, const size_t d ){
//	const double eps = 1.0E-14;
//	const double mathPI = 3.141592653589793;
	if( colinear (p1, p2, p3, d) ){
		cout<<"p1, p2 and p3 colinear "<<endl;
		return 0.0;
	}else if (colinear(p2, p3, p4, d) ){
		cout<<"p2, p3, p3 colinear"<<endl;
		return 0.0;
	}
	double *coef1 =  coefOfPlaneEqu(p1, p2, p3, d);
	double *coef2 =  coefOfPlaneEqu(p2, p3, p4, d);
	double len1 = coef1[0] * coef1[0] +  coef1[1] * coef1[1] + coef1[2] * coef1[2] ;
	double len2 = coef2[0] * coef2[0] +  coef2[1] * coef2[1] + coef2[2] * coef2[2] ;
	if( len1 < eps1 || len2 < eps1 ){
		cerr<<"The length of normal vector is 0 "<<endl;
		exit(0);
	}
	double cosfi = coef1[0] *coef2[0] + coef1[1] *coef2[1] + coef1[2] *coef2[2] / sqrt(len1*len2);
	if (cosfi >=  1.0 ) return 0.0;
	if (cosfi <= -1.0 ) return mathPI;
	delete [] coef1;
	delete [] coef2;
	return  acos(cosfi);
}

