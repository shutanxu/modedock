/*
 * molecular.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef MOLECULAR_H_
#define MOLECULAR_H_

#include<iostream>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#include<string>

#include"chain.h"
#include"bond.h"
#include"ss.h"
#include"../algorithm/matrixCalculate.h"
#include"../algorithm/isomorphism.h"
//#include"../algorithm/convert.h"

using namespace std;

class	Molecular{
public:
	Molecular();
	Molecular( const Molecular& mol );
	Molecular( const string& file );
	Molecular( const bool& empty, const string& na, const string& fn,
			const Coord& ce,const vector<string>& head,
			const vector<Chain>& ch, const vector<Bond>& bd, const SStructure& s ):
				isEmpty(empty),name(na),fileName(fn),center(ce),fileHeader(head),
				chainVec(ch),bondVec(bd),ss(s){
		computeCenter();
	}

    void set_name( const string n ){ name = n; }

	vector<Atom>			get_atomVec()const;
	vector<Atom>            get_atomVec( const string& name );
	vector<Bond> 			get_bondVec()const{ return bondVec; };
	vector<Bond>            get_rotableBondVec()const;
	vector<Residue>			get_residueVec()const ;
	vector<Chain>			get_chainVec()const{ return chainVec; }
	string					get_fileName()const{ return fileName; }
	string					get_name()const{ return name; }
	Coord					get_center()const{ return	center; }
	Atom                    get_centerAtom()const;
	float                   get_centerAtomRadius()const;
	vector<string >			get_fileHeader()const{ return fileHeader; }
	vector<vector<size_t> >	get_ringAtomIndex()const;
	vector<size_t>			get_methyC();
	Atom                    get_atom(const size_t& index)const;

	SStructure				get_secondaryStructure()const{ return ss; }
	vector<Bond>			constrainedRotatableBond( );
	void					print()const;
	void					print( string& fileName );
	void					printPreResidue()const;
	void					set_atomVec( const vector<Atom>& val );
	bool					empty()const{ return isEmpty; }
	bool                    innerCrash( float tolDis=0)const;
	void                    update_center();
	void					update_atomCoord( const Atom& at );
	void					update_atomCoordVec( const vector<Atom>& atVec );
	void					computeSS(){calculateBetaSheets(); calculateAlphaHelices();calculateTurns(); }
	Molecular&				operator=( const Molecular& mol );
	std::string float2string(float number){
	    std::ostringstream buff;
	    buff<<number;
	    return buff.str();
	}

private:
	int						computeCenter();
	void					set_preResidue();

	void					calculateBetaSheets(   );
	void					calculateAlphaHelices(    );
	void					calculateTurns(  );

	bool                    isEmpty;
	string                  name;
	string                  fileName ;
	Coord                   center;
	vector<string >         fileHeader;
	vector<Chain>           chainVec;
	vector<Bond>            bondVec;
	SStructure              ss;						// secondary structure
};

vector<Molecular>
readMolecularFile( const string& fileName );

vector<Molecular>
readPdb( const string& fileName );

vector<Molecular>
readMol2( const string& fileName );

vector<Molecular>
readMol( const string& fileName );

vector<Molecular>
read_pdbqt( const string& fileName );

vector<Atom>
overlay( const Molecular& standMol, const Molecular& oldMol,
         const bool& debug );

vector<Atom>
overlay( const vector<Atom> standAtVec, const vector<Atom> oldAtVec,
         const bool& debug  );

vector<Atom>
overlay(  const Molecular& standMol, const Molecular& oldMol,
        const vector<size_t>& olAtomIndex, const bool& debug );

//---------------------------------------
//Molecular
//overlay( const Molecular& standMol,
//         const Molecular& oldMol,
//         const bool& debug= false );

//Molecular
//overlay( const Molecular& standMol,
//         const Molecular& oldMol,
//         const vector<size_t>& olAtomIndex,
//         const bool& debug= false );

Molecular
overlay( const Molecular& standMol,
         const Molecular& oldMol,
         const string& overlayAtomName,
         const bool& debug= false );
//---------------------------------------

float
overlayRMSD( const Molecular MolA, const Molecular MolB );

float
overlayRMSD_CA( const Molecular MolA, const Molecular MolB );

float
pairRMSD( const Molecular& MolA, const Molecular& MolB );

float
pairRMSD_CA( const Molecular& molA, const Molecular& molB );

float
maxPairRMSD( const vector<Molecular> molVec );

//float
//overlayRMSD( const Molecular& MolA, const Molecular& MolB );

float
isomorphismRMSD( const Molecular MolA, const Molecular MolB, bool showMatrix = false );

bool
molecularCrash( const Molecular& mol1, const Molecular& mol2);

Molecular
rotateAroundBond( const Molecular& rotableMol, const Bond& bd,
		const float& degreeAngle );

bool
rotateAroundBond( const Molecular& rotableMol, const Bond& bd,
		const float& degreeAngle, Molecular& output );

Molecular
rotateAroundBond( const Molecular& rotableMol,
                  const size_t& fixedAtomIndex,
                  const size_t& rotableAtomIndex,
                  const float& degreeAngle );

bool
rotateAroundBond( const Molecular& rotableMol, const size_t& fixedAtomIndex,
		const size_t& rotableAtomIndex, const float& degreeAngle, Molecular& output );

Molecular
rotateAroundAtom( const Molecular& rotableMol, const size_t& centerAtomIndex,
                  const float& roundXAngle, const float& roundYAngle,
                  const float& roundZAngle );

Molecular
rotateAroundCenter( const Molecular& rotableMol, const Coord& center,
                    const float& roundXAngle, const float& roundYAngle,
                    const float& roundZAngle );

Molecular
rotateAroundCoord( const Molecular& rotableMol, const Coord& rotateCenter,
                   const float& roundXAngle, const float& roundYAngle,
                   const float& roundZAngle );

Molecular
translateMol(  const Molecular& mol, const Coord& direction, const float& distance );

vector<Atom>
getBindingSiteAtoms( const Molecular& mol, const Coord& center, const float& diameter );

vector<Atom>
getSubAtomVec( const Molecular& mol, const vector<string>& atNameVec,
			   const vector<size_t>& atIndexVec, const vector<size_t>& atResidNameVec );

vector<Bond>
getRotableBonds( const Molecular& rotableMol );

vector<pair<Atom, Atom> >
getBondAtomPairs( const Molecular& mol, const vector<pair<size_t, string> >& atVec );

vector<pair<Atom, Atom> >
getBondAtomPairs( const Molecular& mol, const vector< Atom >& atVec );

void
getConformationalIsomers( const Molecular& inputMol, const float& bondRotateDegree,
		                  vector<Molecular>& newMol, int bondIndex=0 );

#endif /* MOLECULAR_H_ */
