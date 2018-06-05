/*
 * ss.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef SS_H_
#define SS_H_

#include<algorithm>
#include"atom.h"
#include"residue.h"
//#include"molecular.h"

using namespace std;

struct SStructure{

//	vector< vector< Atom > >				getHelixAtomVec()const{ return helixAtomsVecs; }
//	vector< vector< Atom > >				getSheetAtomVec()const{ return sheetAtomsVecs; }
//	vector< vector< Atom > >				getTurnAtomVec()const{ return turnAtomsVecs; }

//	string												fileName;
//	vector< vector< size_t > >			helixResVecs;
//	vector< vector< size_t > >			sheetResVecs;
//	vector< vector< size_t > >			turnResVecs;

//	bool												isEmpty()const{  }

	vector< vector< Atom > >				helixAtomsVecs;
	vector< vector< Atom > >				sheetAtomsVecs;
	vector< vector< Atom > >				turnAtomsVecs;
};

//---------------------------------------------------------------------------------------------------------------------------

bool
dssp( const vector<Atom>& atVec, const Residue& res1, const Residue& res2 );
void
dssp( const vector<Atom>& atomVec, const vector<Residue>& residueVec );

#endif /* SS_H_ */
