/*
 * residue.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef RESIDUE_H_
#define RESIDUE_H_

#include"atom.h"

using namespace std;

enum MolType{
	MOL_PROTEIN,
	MOL_DNA,
	MOL_RNA,
	MOL_LIGAND,
	MOL_WATHER
};

// --------------------------------------------------------------------
// for computation of hydrogen bond energy by DSSP
const double
	kSSBridgeDistance = 3.0,
	kMinimalDistance = 0.5,
	kMinimalCADistance = 9.0,
	kMinHBondEnergy = -9.9,
	kMaxHBondEnergy = -0.5,
	kCouplingConstant = -27.888,	//	= -332 * 0.42 * 0.2
	kMaxPeptideBondLength = 2.5;

const double
	kRadiusN = 1.65,
	kRadiusCA = 1.87,
	kRadiusC = 1.76,
	kRadiusO = 1.4,
	kRadiusSideAtom = 1.8,
	kRadiusWater = 1.4;

// --------------------------------------------------------------------

class Residue{
public:
	Residue(){ name = ""; index = 0; atomVec.clear();  }
	Residue( const Residue& res );
	Residue( const string& na ){ name = na; index = 0; atomVec.clear(); preResidue = NULL; set_molType();}
	Residue( const string& n, const size_t& i, const vector<Atom>& atVec ):
		name(n), index(i), atomVec( atVec ){ preResidue = NULL; chainName = atVec[0].get_chainName(); set_molType();}
	Residue( const string& n, const size_t& i, const vector<Atom>& atVec,  Residue *pre ):
		name(n), index(i), atomVec( atVec ){ preResidue = pre;  chainName = atVec[0].get_chainName(); set_molType(); }
	~Residue(){ /*atomVec.clear(); delete preResidue;*/ }
	void					set_preResidue(  Residue* pre ){ preResidue = pre; }
	void					set_name( string& na ){ name = na; }

	vector<Atom>	        get_atomVec()const{ return atomVec ; }
	Atom					get_C()const;
	Atom					get_CA()const;
	Atom					get_N()const;
	Atom					get_O()const;
	string					get_name()const{ return name; }
	size_t					get_index()const{ return index; }
	string					get_chainName()const{ return chainName; }
	Residue*				get_preResidue()const{ return	 preResidue; }
	MolType				    get_molType()const{ return type; }

	void					clear(){  name = ""; index = 0; atomVec.clear(); preResidue = NULL; }
	void					print()const;
	bool					empty()const{ return atomVec.empty(); };
	static float			calculateHBondEnergyDSSP(  const Residue& res1, const Residue& pre1, const Residue& res2, const Residue& pre2 );
	static	bool		    isAminoAcid( const Residue& res );
	bool					isAminoAcid();
	static	Residue 	    getResidue( const vector<Residue>& resVec, size_t index );
	Residue&				operator=( const Residue& res );

	void					update_atomCoord( const Atom& at );
private:
	void					set_molType();

	string					name;
	size_t					index ;
	string					chainName ;
	vector<Atom>	        atomVec;
	MolType				    type;
	Residue*				preResidue;
};

class
HasIndex: public unary_function<Residue, bool>{
public:
	HasIndex( size_t ind ): index( ind ){}
	bool	operator()( const Residue& resd )const{
		return				( resd.get_index() == index );
	}
private:
	size_t			index;
};

#endif /* RESIDUE_H_ */
