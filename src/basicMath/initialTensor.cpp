#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "initialTensor.h"
using namespace std;

initialTensor::initialTensor(){}

initialTensor::~initialTensor(){;}


/**
 * The dimension of second axis must be 3
 */
bool initialTensor::initialization(){
	float ixx = 0.0, iyy = 0.0 , izz = 0.0;
	float ixy = 0.0, ixz = 0.0, iyz = 0.0;
	size_t m = (size_t) ( coordVec.size() / dim );
	if( m == 0 ){
		cout<<"no data for the computation of tensor: "<< m <<endl ;
		return false;
	}

	// Compute the center
	size_t p = 0;
	center = vector<float>(dim, 0.0);
	for ( unsigned int i = 0; i < m; i++)
		for ( p = 0; p< dim; p++)
			center[p] += coordVec[ dim * i + p ];
	for ( p = 0; p< dim; p++)
		center[p] /= m;

	float dis = 0.0;
	farthestDistance= 0.0;
	float coord[dim];
	for( unsigned int i = 0; i < m; i++){
		for ( p = 0; p< dim; p++)
			coord[ p ] = coordVec[ dim * i + p ] - center[p] ;

		for(unsigned int p = 0; p < dim; p++){
			if( p !=0)
				ixx += coord[ p ] * coord[ p ];
			if ( p != 1)
				iyy += coord[ p ] * coord[ p ];
			if( p != 2)
				izz += coord[ p ] * coord[ p ];
		}
		ixy -= coord[ 0 ] * coord[1];
		ixz -= coord[ 0 ] * coord[ 2];
		iyz -= coord[ 1] * coord[ 2];
		dis = (coord[ 0 ] * coord[ 0 ] + coord[ 1] * coord[1] + coord[ 2] * coord [ 2 ] );
		if (dis > farthestDistance){
			farthestDistance = dis;
			fpIndex = i;
		}
	}
	farthestDistance = sqrt(farthestDistance);
	farthestPoint = vector<float>(dim, 0.0);
	for ( p = 0; p< dim; p++)
		farthestPoint[p] = coordVec[ dim * fpIndex + p ];

	tensor = vector<float>( dim * dim, 0.0);

	tensor[0] = ixx;
	tensor[1] = ixy;
	tensor[2] = ixz;
	tensor[3] = ixy;
	tensor[4] = iyy;
	tensor[5] = iyz;
	tensor[6] = ixz;
	tensor[7] = iyz;
	tensor[8] = izz;
	return true;
}
//
//initialTensor::initialTensor( const vector<Atom*>& atomVec){
//
//	const size_t m = atomVec.size();
//	if( m == 0 ) {
//		cout<<"no data in the atomVec \n";
//		exit(1);
//	}
//	 vector<float> coord( dim * m, 0.0 );
//	 for (unsigned int i = 0; i< m; i++){
//		 coord [ dim * i ] =  atomVec.at(i)->getX();
//		 coord [ dim * i + 1] = atomVec.at(i)->getY();
//		 coord [ dim * i + 2] = atomVec.at(i)->getZ();
//	 }
//	 initialTensor(coord);
//}


matrix  initialTensor::principalMoment( vector<float> &eigenValue, vector<double>&eigenVector ){


	matrix rot;
	return rot;
}



