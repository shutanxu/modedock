/*
 * bond.cpp
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#include"bond.h"

using namespace std;

bool
operator==( const Bond& bd1, const Bond& bd2 ){
	Atom			a1 = bd1.get_firstAtom();
	Atom			a2 = bd1.get_secAtom();
	Atom			b1 = bd2.get_firstAtom();
	Atom			b2 = bd2.get_secAtom();
	if( ( a1 == b1 && a2 == b2 ) ||
			( a2 == b1 && a1 == b2 ) ){
		return true;
	}else{
		return false;
	}
}

Bond&
Bond::operator =(const Bond& bd){
	index=bd.get_index();
	firstAtom=bd.get_firstAtom();
	secondAtom=bd.get_secAtom();
	type = bd.get_type();
	return *this;
}

void
Bond::set_fistAtomCoord(const Coord co){
	firstAtom.set_coord(co.x, co.y, co.z);
}

void
Bond::set_secAtomCoord(const Coord co){
	secondAtom.set_coord(co.x, co.y, co.z);
}

void
Bond::print()const{
	cout.width(5);
	cout<<"bond: "<<right<<index<<" type:"<<type<<endl;;
	firstAtom.print();
	secondAtom.print();
//	cout<<endl;
}

//----------------------------------------------------------

