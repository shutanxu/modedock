/*
 * overlay.h
 *
 *  Created on: Oct 27, 2013
 *      Author: stan
 */

#ifndef OVERLAY_H_
#define OVERLAY_H_

#include <vector>
#include <map>
#include <iostream>
#include<math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

vector< vector<double> >
 overlayMatrix(const vector<float>&structCoordA,
		 const vector<float>&structCoordB, const bool& debug = false );

//vector< vector<double> >
// overlayMatrix2(const vector<float>&structCoordA,
//		 const vector<float>&structCoordB, const bool& debug = false );

vector<float>
overlayCoord( const vector<float>&structCoordA,
               const vector<float>&structCoordB, const bool& debug=false );



float
overlayRMSD( const vector<float>&structCoordA,
const vector<float>&structCoordB, const bool& debug=false );

#endif /* OVERLAY_H_ */
