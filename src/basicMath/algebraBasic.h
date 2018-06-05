#ifndef ALGEBRABASIC_H_
#define ALGEBRABASIC_H_

#include "Math.h"

#define epsMinus14 1.0E-14
#define epsMinus10 1.0E-10
#define epsMinus7 1.0E-7
//#define mathPI 3.141592653589793
//#define CST  180.00 / 3.141592653589793
//#define toRad  3.141592653589793 / 180.00

class algebraBasic {
public:
	algebraBasic();
	bool trigoSolver(double a, double b, double* sinCosValue);
	virtual ~algebraBasic();
};

/**
 * SOLVE the equation:
 *     $a \sin{x} + b \cos{x} = c$
 *   Compute $x_0$
 *           let $r = \sqrt( a * a + b * b ), \sin{x_0} = a / r, \cos {x_0} = b / r;
 *  IF $ c/ r > 1.0$
 *  return false
 we have then
    $\cos( x - x_0) = c / r$;
    $x - x_0 = \acos(c / r) $;
    SOVLE for $x$
 */

template<class T> bool trigoSolver ( const T& a, const T& b, const T& c,
		vector<T>& solutionVec ) {
	if ( fabs(a) < smallEPS && fabs(b) < smallEPS ) {
		std::cerr<<"Both $a$ and $b$ too small in _trigoSolver_ \n";
		return false;
	}
	T solution;
	if ( fabs(a) < smallEPS ) { //  b \cos{x} = c
		T cb = c / b;
		if( cb > 1.0) {
			if ( (cb - 1.0)  > mediumEPS )
				return false;
			solutionVec.push_back( 0.0 );
			return true;
		}
		else if ( cb < -1.0 ) {
			if ( ( - 1.0 - cb )  > mediumEPS )
				return false;
			solutionVec.push_back( Pi );
			return true;
		}
		solutionVec.push_back(acos( cb ) );
		solutionVec.push_back(-acos( cb ) );
		return true;
	}else if ( fabs(b) < smallEPS ) {  //  $a \sin{x} = c$
		T ca = c / a;
		if( ca > 1.0){
			if ( ( ca - 1.0 ) > mediumEPS )
				return false;
			solutionVec.push_back( PiOverTwo );
			return true;
		}
		else if ( ca < -1.0 ){
			if ( ( - 1.0 - ca )  > mediumEPS )
				return false;
			solutionVec.push_back( TwoPi - PiOverTwo );
			return true;
		}
		solutionVec.push_back(asin( ca ) );
		solutionVec.push_back(Pi - asin( ca ) );
		return true;
	}

	T r = sqrt( a * a + b * b );
	T cr = c / r;
	T y0;
	if( cr > 1.0) { //  cos( cr ) = 1.0;
		if ( (double) ( cr - 1.0 ) > mediumEPS ) {
			return false;
		}else
			y0 = 0.0;
	}
	else if ( cr < -1.0 ) {
		if ( ( - cr - 1.0 ) > mediumEPS )
			return false;
		else
			y0 = Pi;
	}else
		y0  = acos( cr );

	// 2 solutions: y0 and \pi - y0
	T y[2] = {y0, - y0 };
	// 2 solutions:  x0 and -x0
	T x0 = acos( b / r );
	if(  a < 0.0)
		x0 = - x0;

	//The first solution
	solution = y[0] + x0;
	T left= a * sin(solution) + b * cos(solution);
	if ( fabs(left -c) < mediumEPS )
		solutionVec.push_back(solution);
	// the second solution
	solution = y[1] +  x0;
	T left2= a * sin(solution) + b * cos(solution);
	if ( fabs(left2 -c) < mediumEPS )
		solutionVec.push_back(solution);
	if( solutionVec.empty() ) {
		std::cerr<<"No final solution  in _trigoSolver_ \n";
//		printf("y0=%12.8f, x0=%12.8f, s0=%12.8f\n ",
//				y0*RadToDegrees, x0*RadToDegrees, y[0] + x0 );
//		printf("left = %12.8f, eps =%12.8f\n ", left, fabs(left - c) );
//		printf("y0=%12.8f, x0=%12.8f, s0=%12.8f\n ",
//				y0*RadToDegrees, x0*RadToDegrees, y[1] + x0 );
//		printf("left2 = %12.8f, eps =%12.8f\n ", left2, fabs(left2 - c) );
		return false;
	}
	return true;
}

#endif /*ALGEBRABASIC_H_*/
