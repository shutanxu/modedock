/*
 * bond.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef BOND_H_
#define BOND_H_

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<vector>
#include<math.h>

#include"atom.h"

using namespace std;

class Bond{
public:
	Bond(){}
	Bond( const size_t& ind, const Atom& at1, const Atom& at2, const string& ty ):
		index(ind),firstAtom(at1),secondAtom(at2),type(ty){}

	size_t		get_index()const{ return index; }
	Atom		get_firstAtom()const{ return firstAtom; }
	size_t		get_firstAtomIndex()const{ return firstAtom.get_index(); }
	string		get_firstAtomName()const{ return firstAtom.get_name(); }
	size_t		get_firstAtomResidIndex()const{ return firstAtom.get_residueIndex(); }

	size_t		get_secAtomIndex()const{ return secondAtom.get_index() ; }
	Atom		get_secAtom()const{ return secondAtom; }
	string		get_secAtomName()const{ return secondAtom.get_name() ; }
	size_t		get_secAtomResidIndex()const{ return secondAtom.get_residueIndex() ; }
	string		get_type()const{ return	type; }
	friend	bool					operator == (const Bond& bd1, const Bond& bd2);
	Bond&       operator=(const Bond& bd);
	void        set_fistAtomCoord(const Coord co);
	void        set_secAtomCoord(const Coord co);
	void		print()const;
private:
	size_t		index;
	Atom		firstAtom;
	Atom		secondAtom;
	string		type;
};

#endif /* BOND_H_ */
