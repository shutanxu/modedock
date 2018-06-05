/*
 * atom.h
 *
 *  Created on: Sep 21, 2013
 *      Author: stan
 */

#ifndef ATOM_H_
#define ATOM_H_

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<vector>
#include<math.h>

#include"coord.h"
#include"../algorithm/tokenBySpace.h"

using namespace std;

static const char* ResidueType[]={
			"ALA",
			"ARG",
			"ASN",
			"ASP",
			"CYS",
			"GLU",
			"GLN",
			"GLY",
			"HIS",
			"ILE",
			"LEU",
			"LYS",
			"MET",
			"PHE",
			"PRO",
			"SER",
			"THR",
			"TRP",
			"TYR",
			"VAL",
			"MSE"
};

static const char* DNAtype[]={
	"DA", "DT","DG","DC","A","T","G","C",
};

static const char* RNAtype[]={
	"DA", "DG","DC","DU","A","G","C","U",
};

class Atom{
public:
	Atom();
	Atom( const vector<string>& strVec );
	Atom( const string& oneLine );
	Atom( const size_t& id, const string& na, const string& tp, const Coord& co ):
		index(id), name(na), type( tp ), coord( co ){set_vdwRadiusDefault();}
	Atom( const size_t& id, const string& na, const string& ty, const string& ch, const Coord& co,
			const string& resdNm ,const size_t& resdInd, const float& chg ):
			index(id), name(na), type(ty), chainName(ch), coord(co),
			residueName(resdNm), residueIndex(resdInd), charge( chg ){ vdwRadius = 0; set_vdwRadiusDefault(); }

	bool				 isResidue();
	bool				 isDNA();
	bool				 isRNA();
	bool				 isLigand();

	void                 readMol2( const vector<string>& strVec );

	size_t				 get_index()const{ return index; }
	string				 get_name()const{ return name; }
	string				 get_type()const{ return type; }
	size_t				 get_sequ()const{ return	sequ; }
	string				 get_residueName()const{ return residueName ; }
	size_t				 get_residueIndex()const{ return residueIndex ;  }
	string				 get_chainName()const{ return chainName; }
	Coord				 get_coord()const{ return coord; }
	float				 get_charge()const{ return	charge; }
	float				 get_vdwRadius()const{ return vdwRadius; }
	vector<float>		 get_color();
    vector<Atom>         get_attachedH()const{ return attachedH; }
    vector<Atom>         get_attachedC()const{ return attachedC; }

	void				 set_index( const size_t& val ){ index = val; }
	void				 set_name( const string& val ){ name = val; }
	void				 set_type( const string& val ){ type = val; }
	void				 set_coord( const float& x, const float& y, const float& z ){ coord.x = x; coord.y = y; coord.z = z; }
	void				 set_charge( const float& ch ){ charge = ch; }
	void				 set_vdwRadiusDefault();
	void				 set_vdwRadius( const float& vdwR ){ vdwRadius = vdwR ; }
    void                 set_attachedH( const vector<Atom> atVec ){ attachedH = atVec; }
    void                 set_attachedC( const vector<Atom> atVec ){ attachedC = atVec; }
	int					 print()const;
	friend bool			 operator == (const Atom& at1, const Atom& at2);
	Atom&                operator=(const Atom& );
private:
	size_t				 index;
	string				 name;
	string				 type;
	Coord				 coord;
	size_t				 sequ;
	string				 residueName;
	size_t				 residueIndex;
	string				 chainName;
	float				 charge;
	float				 vdwRadius;
    vector<Atom>         attachedH;
    vector<Atom>         attachedC;
};

class AtomEqual{
public:
	AtomEqual( const Atom& a ):at(a){}
	bool								operator()( const Atom a )const{ return a == at; }
private:
	Atom		at;
};

struct	SortByAtomIndex{
	inline	bool				operator()( const Atom& a1, const Atom& a2 ){
		return		a1.get_index() < a2.get_index();
	}
};

inline bool
operator!=(  const Atom& at1, const Atom& at2 ){
	return		!( at1==at2 );
}

Atom
searchAtom( const vector<Atom>& atomVec, const size_t& index );

Atom
searchAtom(const vector<Atom>& atomVec, const size_t& index,
		  const string& name, const size_t& residIndex, const string& chainName) ;

Atom
searchAtom( const vector<Atom>& atomVec, const string& atomName,
		 const size_t& atomIndex, const size_t& residIndex );

vector<Atom>
searchAtom( const vector<Atom>& atomVec, const string& atomName );

Coord
centerOfAtoms( const vector<Atom>& atomVec );

bool
atomsCrash(const Atom& at1, const Atom& at2 );

bool
atomsCrash( const vector<Atom>& atVec1, const vector<Atom>& atVec2, bool computeHydrogen=true );

bool
atomsCrash( const Atom& at, const vector<Atom>& atVec, bool computeHydrogen=true );

float
atomDis( const Atom& at1,  const Atom& at2 );

float
atomVecRMSD( const vector<Atom>& at1,  const vector<Atom>& at2 );

float
atomSquareDis( const Atom& at1,  const Atom& at2 );

float
atomVecSquareRMSD( const vector<Atom>& at1,  const vector<Atom>& at2 );

pair<Atom, Atom>
findClosestAtom( const vector<Atom>& atomVec1,
		const vector<Atom>& atomVec2 );

Atom
findClosestAtom( const Atom& atom1,
		const vector<Atom>& atomVec2 );

float
atomAngle( const Atom& dir1Start, const Atom& dir1End,
		const Atom& dir2Start, const Atom& dir2End );

float
atomAngle( const Atom& end1, const Atom& vertexAt, const Atom& end2 );

float
atomsMaxRange( const vector<Atom>& atomVec );

bool
setAttachedH( Atom& at, const vector<Atom> atVec );

bool
setAttachedC( Atom& at, const vector<Atom> atVec );

float
atomsMinX(  const vector<Atom>& atomVec );
float
atomsMaxX(  const vector<Atom>& atomVec );
float
atomsMinY(  const vector<Atom>& atomVec );
float
atomsMaxY(  const vector<Atom>& atomVec );
float
atomsMinZ(  const vector<Atom>& atomVec );
float
atomsMaxZ(  const vector<Atom>& atomVec );

#endif /* ATOM_H_ */
