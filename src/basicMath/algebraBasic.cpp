#include <cmath>
#include <iostream>
#include "algebraBasic.h"
using namespace std;

algebraBasic::algebraBasic(){}

algebraBasic::~algebraBasic(){}

bool algebraBasic:: trigoSolver(double a, double b, double* sinCosValue){
	if ( fabs(a) < epsMinus14 )
		return false;
	else if ( fabs(b) < epsMinus14){
		sinCosValue[0] = 1.0;
		sinCosValue[1] = 0.0;
		sinCosValue[2] = 1.0;
		sinCosValue[3] = 0.0;
		return true;
	}else{
		double t = ( b + sqrt( a*a + b * b) ) / a;
		double denom = 1.0 + t * t;
		sinCosValue[0] = 2.0 * t / denom;
		sinCosValue[1] = ( 1.0 - t * t ) / denom;

		t = ( b - sqrt( a*a + b * b) ) / a;
		denom = 1.0 + t * t;
		sinCosValue[2] = 2.0 * t / denom;
		sinCosValue[3] = ( 1.0 - t * t ) / denom;
	}
	return true;
}
