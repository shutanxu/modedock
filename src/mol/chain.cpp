/*
 * chain.cpp
 *
 *  Created on: May 15, 2014
 *      Author: stan
 */

#include"chain.h"

using namespace std;

void
SubChain::add_Residue( const Residue& res ){
	if( ! residueVec.empty() ){
		if( type != res.get_molType() ){
//			cout<<"--- warnning!!! different residue type in chain, is below ligand ? "<<endl;
//			res.print();
		}
		residueVec.push_back( res );
		type = res.get_molType();
//		res.print();
	}else{
		residueVec.push_back( res );
		type = res.get_molType();
	}
}

vector<Atom>
SubChain::get_AtomVec()const{
	vector<Atom>			atVec;
	atVec.clear();

	for( size_t i=0; i<residueVec.size(); i++ ){
//		cout<<"--- residue:"<<i<<endl;
		if( ! residueVec[i].get_atomVec().empty() ){
//			atVec.insert( atVec.end(), residueVec[i].get_atomVec().begin(), residueVec[i].get_atomVec().end() );
			vector<Atom>		resAtVec = residueVec[i].get_atomVec();
			for( size_t j=0; j<resAtVec.size(); j++ ){
				atVec.push_back( resAtVec[j] );
			}
		}
	}
	return							atVec;
}

void
SubChain::update_atomCoord( const Atom& at ){
		for( size_t i=0; i<residueVec.size(); i++ ){
			residueVec[i].update_atomCoord( at );
		}
}

SubChain&
SubChain::operator =( const SubChain& rf ){
	type = rf.get_molType();
//	vector<Residue>		resVec( rf.get_ResidueVec().size() );
	vector<Residue>     tempRes( rf.get_ResidueVec().size() );
	residueVec = tempRes;
	for( size_t i=0; i<residueVec.size(); i++ ){
		residueVec[i] = rf.get_ResidueVec()[i];
	}

	return *this;

}

//----------------------------------------------------------------------------------------------------------------------------

bool
Chain::empty(){
	for( size_t i=0; i<subChainVec.size(); i++ ){
		if( subChainVec[i].get_ResidueVec().size() > 0 ){
			return false;
		}
	}
	return true;
}

void
Chain::set_ResidueVec( const vector<Residue>& rsVec ){
	subChainVec.clear();
	for( size_t i=0; i<rsVec.size(); i++ ){
		bool success = false;
		for( size_t j=0; j<subChainVec.size(); j++ ){
			if( subChainVec[j].get_molType() == rsVec[i].get_molType() ){
				subChainVec[j].add_Residue( rsVec[i] );
				success = true;
				break;
			}
		}
		if( ! success ){
			SubChain		sub( rsVec[i] );
			subChainVec.push_back( sub );
		}
	}
}

vector<Atom>
Chain::get_AtomVec()const{
	vector<Atom>			atVec;

	for( size_t i=0; i<subChainVec.size(); i++ ){
		if( ! subChainVec[i].get_AtomVec().empty() ){

			vector<Atom>		subAtVec = subChainVec[i].get_AtomVec();
			for( size_t j=0; j<subAtVec.size(); j++ ){
				atVec.push_back( subAtVec[j] );
			}
		}
	}
	return							atVec;
}

vector<Residue>
Chain::get_ResidueVec()const{
	vector<Residue>		resVec;
	for( size_t i=0; i<subChainVec.size(); i++ ){
		if( ! subChainVec[i].get_ResidueVec().empty() ){
            vector<Residue> tempResVec = subChainVec[i].get_ResidueVec();
            for( size_t j=0; j<tempResVec.size(); j++ ){
                resVec.push_back( tempResVec[j] );
			}
		}
	}
	return		resVec;
}

Chain&
Chain::operator =( const Chain& rf ){
	name = rf.name;
//	vector<SubChain>		subChainVec( rf.get_subChainVec().size() );
	vector<SubChain>        sub(  rf.get_subChainVec().size() );
	subChainVec = sub;

	for( size_t i=0; i<subChainVec.size(); i++ ){
		subChainVec[i] = rf.get_subChainVec()[i];
	}
//	cout<<"subChain:"<<subChainVec.size()<<endl;
	return *this;
//	residueVec = resVec;
}

void
Chain::update_atomCoord( const Atom& at ){
	for( size_t i=0; i<subChainVec.size(); i++ ){
		subChainVec[i].update_atomCoord( at );
	}
}

void
Chain::add_Residue( const Residue& res ){
	for( size_t i=0; i<subChainVec.size(); i++ ){
//		if( subChainVec[i].get_molType() == res.get_molType() ){
		if( (subChainVec[i].get_molType() == res.get_molType()) ||
			( subChainVec[i].get_molType() == MOL_PROTEIN &&
			  res.get_molType() == MOL_LIGAND ) ){
			subChainVec[i].add_Residue( res );
			return;
		}
	}
	SubChain			subCh( res );
	subChainVec.push_back( subCh );
}

void
Chain::print(){
	cout<<"chain: "<<name<<endl;
	for( size_t i=0; i<subChainVec.size(); i++ ){
//		cout<<MOL_PROTEIN<<endl;
		cout<<"  subchain:"<<subChainVec[i].get_molType()<<" "<<subChainVec[i].get_ResidueVec().size()<<endl;
		for( size_t j=0; j<subChainVec[i].get_ResidueVec().size(); j++ ){
			subChainVec[i].get_ResidueVec()[j].print();
		}
	}
}
