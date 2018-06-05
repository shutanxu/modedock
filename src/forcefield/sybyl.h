/*
 * sybyl.h
 *
 *  Created on: Nov 15, 2013
 *      Author: stan
 */

#ifndef SYBYL_H_
#define SYBYL_H_

#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<vector>

#include"../mol/atom.h"
#include"../mol/molecular.h"

using namespace std;

struct SybylAtom{
public:
	SybylAtom(){}
	SybylAtom( const string& ty, const float& vd, const float& ta,
			const float& io, const float& po, const float& we, const float& fo,
			const string& ge, const size_t&nu, const bool& id, const bool& ia ):
				type(ty), vdwRadius(vd), taffK(ta), ionisationPotential(io), polarizability(po),
				weight(we), formalCharge(fo), geometry(ge), numOfNeighbors(nu),
				isDonor(id), isAcceptor(ia){ cutoffDis = 1.5; };
	SybylAtom( const vector<string>& vec );
	void 					print();
	string					type;
	float					vdwRadius;
	float					taffK;
	float					ionisationPotential;
	float					polarizability;
	float					weight;
	float					formalCharge;
	string					geometry;
	size_t					numOfNeighbors;
	bool					isDonor;
	bool					isAcceptor;
	float					cutoffDis;
};

/*
 * # Directive to indicate the start of the hydrogen bond table.
START_H_BOND_TABLE

# Each H-Bond entry has 4 fields for donors, 5 fields for acceptor or
# donor acceptor and 6 fields for metals.
# 1. The hydrogen bonding type.
# 2. The atom type associated with this hydrogen bonding type.
# 3. The entry for the atom type is either default, DEF, or elucidated,
# ELU.
# 4. The entry is either for a donor, D; acceptor, A; donor-acceptor
# DA or metal, M.
# 5. For an acceptor or donor-acceptor field 5 is NONE, LP_DIR or
# LP_PLANE, indicating that the acceptor does have any lone pair
# directionality when forming h-bonds, forms h-bonds in the direction
# of lone pairs, or in the plane of the lone pairs, respectively.
# 5. For a metal this field is a comma separated list of co-ordination
# numbers.  Currently only 4 and 6 are supported, though bifurcated
# coordination geometries are detected.
# 6. The coordination distance for metals.
 */
struct SybylHydrogenBondParam {
	SybylHydrogenBondParam(){}
	SybylHydrogenBondParam( const string& bt, const string& at, const string& atE,
			const string& da, const string& lp, const string& cn, const float& md):
			bondType(bt), atomType(at), atomTypeEntry(atE), donorOrAcceptor(da),
			lonePairDirection(lp), coordinationNumber( cn ), metalDistance(md){};
	SybylHydrogenBondParam( const vector<string>& vec );

	string					bondType;
	string					atomType;
	string					atomTypeEntry;
	string					donorOrAcceptor;
	string					lonePairDirection;
	string					coordinationNumber;
	float					metalDistance;
};

struct SybylHBondEnergyParam{

	SybylHBondEnergyParam( const string& dt, const string& at, const float& me ):
		donorType( dt ), acceptorType( at ), maxEnergy( me ){}
	SybylHBondEnergyParam( const vector<string>& strVec );
	string					donorType;
	string					acceptorType;
	float					maxEnergy;
};

struct	SybylHydrogenBond{
	SybylHydrogenBond(const SybylAtom& ac,  const SybylAtom& d, const Atom& hy,
			const Atom& lp, const float& hyDis, const float& hyAngle, const float& en ): acceptorAtom(ac), donorAtom( d ),
			hydrogenAtom(hy), lonePair(lp), hydrogenLpDistance( hyDis ), hydrogenLpAngle( hyAngle ), energy( en ){}
	~SybylHydrogenBond(){};
	SybylAtom					acceptorAtom;
	SybylAtom					donorAtom;
	Atom						hydrogenAtom;
	Atom						lonePair;
	float						hydrogenLpDistance;
	float						hydrogenLpAngle;
	float						energy;
};

SybylHydrogenBond
get_maxHB( const vector<SybylHydrogenBond>& hbVec );

struct VdwCrash{
	VdwCrash( const string s1, const string s2, const float d ):
		n1(s1),n2(s2),dis(d){}
	string n1;    // name of the first atom
	string n2;    // name of the first atom
	float  dis;   // minimum distance of atom n1 n2 have VDW crash
};

class Sybyl{
public:
	Sybyl(){}
	Sybyl( const string& fileName ):ffpFileName( fileName ){
		readFile( fileName );
		initCrashDis();
	}
	void			readFile(const string& fileName );
	void			readGoldParmFile( const string& fileName );
	vector<Atom>	getHydrogenBondDonorAtoms( const vector<Atom>& atVec );
	vector<Atom>	getHydrogenBondAcceptorAtoms( const vector<Atom>& atVec );
	float			computeGoldVDW( const Molecular& molA, const Molecular& molB );
	float			computeGoldVDW( const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
	float			computeGoldVDW( const Atom& at1, const Atom& at2 );
	float			computeVDW( const Atom& at1, const Atom& at2 );
	float			computeVDW( const Atom& at1, const Atom& at2, const float dis );
	float			computeVDW( const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
	float			computeVDW( const Molecular mol );
	float			computeCutVDW( const Atom& at1, const Atom& at2 );
    float			computeCutVDW( const vector<Atom>& atVec1, const vector<Atom>& atVec2, bool computeHydrogen=true );

	float           computeESP( const Atom at1, const Atom at2);
	float           computeESP( const vector<Atom> v1, const vector<Atom> v2 );

	float			computeHBondEnergy(  const Molecular& molA, const Molecular& molB );
	float			computeHBondEnergy2(  const Molecular& molA, const Molecular& molB );
	float			computeHBondEnergy(  const Atom& at1, const Atom& at2 );
	float			computeHBondEnergy( const Atom& acceptor,	const Atom	& lonePair,
									    const Atom& donor, const Atom& hydrogen);
	vector<SybylHydrogenBond>			getHydrogenBondVec()const{ return hbVec; }
	Coord           getIdealHBondAtom( const Atom& donor, const Atom& donorH,
			                           const string& acceptorName);
	bool            atomsCrash( const Atom at, const Atom at2 );
	bool            atomsCrash( const Atom at, const vector<Atom> atVec );
	bool            atomsCrash( const vector<Atom> atVec1, const vector<Atom> atVec2 );

private:
	SybylAtom		getSybylAtom( const Atom& at );
	string			getHbondType( const string& atomType, const bool& isDonor );
	float			getHbondEnergy( const string& bondType1, const string& bondType2 );
	vector<Atom>	getAcceptorLonePair( const Molecular& mol, const size_t& acceptorIndex );
	vector<Atom>	getDonorHydrogen( const Molecular& mol, const size_t& donorIndex );
	void            initCrashDis();
private:
	string								ffpFileName;
	vector<SybylAtom>					atomVec;
	vector<SybylHydrogenBondParam>		hydrogenBondParamVec;
	vector<SybylHBondEnergyParam>		hbondEnergyParamVec;
	vector<SybylHydrogenBond>			hbVec;
	vector<VdwCrash>                    crashDisVec;
};

#endif /* SYBYL_H_ */
