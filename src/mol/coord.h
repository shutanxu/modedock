/*
 * coord.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef COORD_H_
#define COORD_H_

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<vector>
#include<math.h>

using namespace std;

struct Coord{
	Coord(){ x=y=z = 0; };
	Coord( const float& a, const float& b, const float& c): x(a),y(b),z(c){}


	float		getDis()const{ return sqrt( x*x + y*y + z* z ); }
	Coord		getNormalizeCoord( ){
		if( getDis() == 0 ){
			Coord		co( x, y, z );
			return		co;
		}else{
			Coord		co( x/getDis(), y/getDis(), z/getDis() );
			return		co;
		}
	}
	void			print()const{ cout<<"x:"<<x<<" y:"<<y<<" z:"<<z<<endl; }

	Coord&          operator=(const Coord& c){ x=c.x; y=c.y; z=c.z; return *this; }
	float	x,y,z;
};

double
getCoordDis( const Coord& co1, const Coord& co2 );

inline bool
operator==( const Coord& c1, const Coord& c2 ){
	return		( c1.x == c2.x && c1.y == c2.y && c1.z == c2.z );
}

inline bool
operator!=( const Coord& c1, const Coord& c2 ){
	return		!( c1 == c2 );
}

inline Coord
cross( const Coord& a, const Coord& b ){
	return Coord( a.y*b.z - a.z*b.y,
		       a.z*b.x - a.x*b.z,
		       a.x*b.y - a.y*b.x );
}

inline double
operator*(const Coord &a, const Coord &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Coord
operator*(const Coord &a, const int &b)
{
	return 	Coord( a.x*b, a.y*b, a.z*b );
}

inline Coord
operator/(const Coord &a, const float &b){
	if( b == 0 ){
		return a;
	}
	return Coord( a.x/b, a.y/b, a.z/b );
}

inline Coord
operator-( const Coord &a, const Coord &b ){
//	cout<<"a-b:"<<a.x-b.x<<endl;
	return		Coord( a.x-b.x, a.y-b.y, a.z-b.z );
}

inline Coord
operator+( const Coord &a, const Coord &b ){
//	cout<<"a-b:"<<a.x-b.x<<endl;
	return		Coord( a.x+b.x, a.y+b.y, a.z+b.z );
}

vector<Coord>
inverseCoordVec( const vector<Coord>& cVec );

#endif /* COORD_H_ */
