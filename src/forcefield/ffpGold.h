/*
 * ffpGold.h
 *
 *  Created on: Nov 11, 2013
 *      Author: stan
 */

#ifndef FFPGOLD_H_
#define FFPGOLD_H_

#include<iostream>
#include<stdlib.h>
#include<vector>

#include"../algorithm/tokenBySpace.h"
#include"../mol/atom.h"
#include"../mol/molecular.h"

using namespace std;

//////////////////////////////////////  gold.params  /////////////////////////////////////////////

const	float		GOLD_VDW_MINDIS = 3.0;

struct GOLD_ATOM_PARAM{
	GOLD_ATOM_PARAM(){}
	GOLD_ATOM_PARAM( const string& ty, const float& vd, const float& ta,
			const float& io, const float& po, const float& we, const float& fo,
			const string& ge, const size_t&nu, const bool& id, const bool& ia ):
				type(ty), vdwRadius(vd), taffK(ta), ionisationPotential(io), polarizability(po),
				weight(we), formalCharge(fo), geometry(ge), numOfNeighbors(nu),
				isDonor(id), isAcceptor(ia){ /*minCutoffDis = 3.0; maxCutoffDis = 30.0;*/ };
	void						print()const;
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
//	float					minCutoffDis;
//	float					maxCutoffDis;
};

struct GOLD_H_BOND{
	string					type;
	string					atomType;
	string					defaultElucidated;
};

struct GOLD_HB_ENERGIES{
	string					donorType;
	string					acceptorType;
	int						maxHBEnergy;
};

struct GOLD_TORSION{

};

struct
GOLD_PARAMETER{
	vector<GOLD_ATOM_PARAM>		atomParVec;
	vector<GOLD_H_BOND>				hbParVec;
	vector<GOLD_HB_ENERGIES>		hbEnergyVec;
	vector<GOLD_TORSION>				torParVec;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

GOLD_PARAMETER
readGoldParamFile( const string& fileName );

float
goldVDW( const Atom& at1, const Atom& at2, const GOLD_PARAMETER& gp );

float
goldVDW( const Molecular& mol1, const Molecular& mol2, const GOLD_PARAMETER& gp );

#endif /* FFPGOLD_H_ */
