/*
 * matrixCalculate.h
 *
 *  Created on: Nov 26, 2013
 *      Author: stan
 */

#ifndef MATRIXCALCULATE_H_
#define MATRIXCALCULATE_H_

#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<fstream>
#include<time.h>
#include<sstream>

using namespace std;

template<typename _Tp>
vector<_Tp>
matrixMultiplyVector( const vector< vector<_Tp> >& matrix, const vector<_Tp>& vec ){
	vector<_Tp>		newVec( matrix.size(), 0 ) ;

	for( size_t i=0; i<matrix.size() ; i++ ){
		if( matrix[i].size() != vec.size() ){
			throw			"can not multiply matrix with vector ";
		}
		for( size_t j=0; j<matrix[i].size() ; j++ ){
			newVec[i] += matrix[i][j] * vec[j];
		}
	}
	if(0){
		cout.width(15);
		cout<<"old vector";
		for( size_t i=0; i<vec.size(); i++ ){
			cout.width(10);
			cout<<vec[i];
		}
		cout<<endl;
		cout.width(15);
		cout<<"new vector";
		for( size_t i=0; i<newVec.size(); i++ ){
			cout.width(10);
			cout<<newVec[i];
		}
		cout<<endl;
	}
	return					newVec;
}

#endif /* MATRIXCALCULATE_H_ */
