/*
 * coord.cpp
 *
 *  Created on: Nov 15, 2014
 *      Author: stan
 */

#include"coord.h"

double
getCoordDis( const Coord& co1, const Coord& co2 ){
	return 	sqrt( (double(co1.x)-double(co2.x))*(double(co1.x)-double(co2.x)) +
			(double(co1.y)-double(co2.y))*(double(co1.y)-double(co2.y)) +
			(double(co1.z)-double(co2.z))*(double(co1.z)-double(co2.z)) );
}

vector<Coord>
inverseCoordVec( const vector<Coord>& cVec ){
	vector<Coord> coordVec;
	if( cVec.size() == 0 ){
		cout<<"#############!"<<endl;
		return cVec;
	}
	for( int i=cVec.size()-1; i>-1; i-- ){
		coordVec.push_back( cVec[i] );
	}
	return coordVec;
}

