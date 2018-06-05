/*
 * chain.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef CHAIN_H_
#define CHAIN_H_

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<vector>
#include<math.h>

#include"residue.h"

using namespace std;

// subchain indicates the residues, dna, rna, ligand or water
class SubChain{
public:
	SubChain(){residueVec.clear();};
	SubChain( const SubChain& sub ){ type = sub.get_molType(); residueVec = sub.get_ResidueVec(); }
	SubChain( const Residue& res ){ type = res.get_molType(); residueVec.push_back( res ); }
	vector<Residue>    get_ResidueVec()const{ return residueVec; }
	vector<Atom>       get_AtomVec()const;
	void               add_Residue( const Residue& res );
	MolType            get_molType()const{ return type; }
	void               update_atomCoord( const Atom& at );
	SubChain&          operator=( const SubChain& rf );
private:
	MolType            type;				// protein, DNA, RNA, ligand or water
	vector<Residue>    residueVec;
};

class Chain{
public:
	Chain(){ subChainVec.clear(); }
	Chain( const Residue& res ){ SubChain sub(res); name = res.get_chainName(); subChainVec.push_back( sub ); }
	Chain( const Chain& ch ){ name = ch.get_name(); subChainVec = ch.get_subChainVec(); }
	void								add_Residue( const Residue& res );

//	void								set_ResidueVec( const vector<Residue>& rsVec ){ residueVec = rsVec; }
	void								set_ResidueVec( const vector<Residue>& rsVec );
	void								set_name( const string& st ){ name = st; }
	string                              get_name()const{ return name; }
	vector<Atom>                        get_AtomVec()const;
	vector<Residue>                     get_ResidueVec()const;
	vector<SubChain>                    get_subChainVec()const{ return subChainVec; }
	void                                clear(){ name = ""; subChainVec.clear(); }
	bool                                empty();
	Chain&                              operator=( const Chain& rf );
	void                                update_atomCoord( const Atom& at );
	void                                print();
private:
	string                              name;
//	vector<Residue>		residueVec;
	vector<SubChain>                    subChainVec;
};

#endif /* CHAIN_H_ */
