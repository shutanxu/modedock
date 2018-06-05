/*
 * residue.cpp
 *
 *  Created on: May 12, 2014
 *      Author: stan
 */

#include"residue.h"

using namespace std;

Residue::Residue( const Residue& res ){
	name = res.get_name();
	index = res.get_index();
	chainName = res.get_chainName();
	atomVec = res.get_atomVec();
	preResidue = res.get_preResidue();
	type = res.get_molType();
}

Atom
Residue::get_C()const{
	Atom	at;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_name() == "C" ){
			return	atomVec[i];
		}
	}
	return	at;
}

Atom
Residue::get_CA()const{
	Atom	at;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_name() == "CA" ){
			return	atomVec[i];
		}
	}
	return	at;
}

Atom
Residue::get_N()const{
	Atom	at;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_name() == "N" ){
			return	atomVec[i];
		}
	}
	return	at;
}

Atom
Residue::get_O()const{
	Atom	at;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_name() == "O" ){
			return	atomVec[i];
		}
	}
	return	at;
}

void
Residue::print()const{
	for( size_t i=0; i<atomVec.size(); i++ ){
		atomVec[i].print();
	}
}

float
Residue::calculateHBondEnergyDSSP( const Residue& res1, const Residue& pre1, const Residue& res2, const Residue& pre2 ){
	float 			energy1 = 0;

	//------------------------------------------ res1 donor, res2 acceptor ----------------------------------------------------
	Atom			h1 = res1.get_N();
	Atom			n1 = res1.get_N();
	Atom			o2 = res2.get_O();
	Atom			c2 = res2.get_C();
	if( res1.get_name() != "PRO" && res1.get_name() != "HOH" ){
//		cout<<"res1:"<<res1.get_index()<<","<<res1.get_name()<<" p:"<<res1.get_preResidue()<<endl;
//		cout<<"	pre:"<< res1.get_preResidue()->get_index()<<","<< res1.get_preResidue()->get_name()<<endl;
		Atom		pc = pre1.get_C();
		Atom		po = pre1.get_O();
//		cout<<"pc: ";
//		pc.print();
//		cout<<"po: ";
//		po.print();
		float		dist = atomDis( pc, po );
		h1.set_coord(  (	(pc.get_coord().x - po.get_coord().x ) / dist),
				( (pc.get_coord().y - po.get_coord().y ) / dist),
				( (pc.get_coord().z - po.get_coord().z ) / dist) );
	}

	float			distHO = atomDis(  h1, o2 );
	float			distHC = atomDis(  h1, c2 );
	float			distNC = atomDis(  n1, c2 );
	float			distNO = atomDis(  n1, o2 );

	if (distHO < kMinimalDistance || distHC < kMinimalDistance || distNC < kMinimalDistance || distNO < kMinimalDistance)
		energy1 = kMinHBondEnergy;
	else
		energy1 = kCouplingConstant / distHO - kCouplingConstant / distHC + kCouplingConstant / distNC - kCouplingConstant / distNO;

	//---------------------------------------------- res1 acceptor, res2 donor ------------------------------------------------
	float				energy2 = 0;
	h1 = res2.get_N();
	n1 = res2.get_N();
	o2 = res1.get_O();
	c2 = res1.get_C();
	if( res2.get_name() != "PRO" && res2.get_name() != "HOH"  ){
		Atom		pc = pre2.get_C();
		Atom		po = pre2.get_O();
		float		dist = atomDis( pc, po );
		h1.set_coord(  (	(pc.get_coord().x - po.get_coord().x ) / dist),
				( (pc.get_coord().y - po.get_coord().y ) / dist),
				( (pc.get_coord().z - po.get_coord().z ) / dist) );
	}

	distHO = atomDis(  h1, o2 );
	distHC = atomDis(  h1, c2 );
	distNC = atomDis(  n1, c2 );
	distNO = atomDis(  n1, o2 );

	if (distHO < kMinimalDistance || distHC < kMinimalDistance || distNC < kMinimalDistance || distNO < kMinimalDistance)
		energy2 = kMinHBondEnergy;
	else
		energy2 = kCouplingConstant / distHO - kCouplingConstant / distHC + kCouplingConstant / distNC - kCouplingConstant / distNO;

	//-------------------------------------------------------------------
	float			energy = 0;
	energy = energy1 < energy2 ? energy1 : energy2;

	if( energy < kMinHBondEnergy ){
		energy = kMinHBondEnergy;
	}
	return			energy;
}

bool
Residue::isAminoAcid( const Residue& res ){
	string		type = res.get_name();
	bool		flag = false;
	if( type == "ALA" ){
		flag = true;
	}else if( type == "ARG" ){
		flag = true;
	}else if( type == "ASN" ){
		flag = true;
	}else if( type == "ASP" ){
		flag = true;
	}else if( type == "CYS" ){
		flag = true;
	}else if( type == "GLU" ){
		flag = true;
	}else if( type == "GLY" ){
		flag = true;
	}else if( type == "GLN" ){
		flag = true;
	}else if( type == "HIS" ){
		flag = true;
	}else if( type == "ILE" ){
		flag = true;
	}else if( type == "LEU" ){
		flag = true;
	}else if( type == "LYS" ){
		flag = true;
	}else if( type == "MET" ){
		flag = true;
	}else if( type == "PHE" ){
		flag = true;
	}else if( type == "PRO" ){
		flag = true;
	}else if( type == "SER" ){
		flag = true;
	}else if( type == "THR" ){
		flag = true;
	}else if( type == "TRP" ){
		flag = true;
	}else if( type == "TYR" ){
		flag = true;
	}else if( type == "VAL" ){
		flag = true;
	}else if( type == "MSE" ){
		flag = true;
	}

	return flag;
}




Residue
Residue::getResidue( const vector<Residue>& resVec, size_t index ){
	Residue res;
	for( size_t i=0; i<resVec.size(); i++ ){
		if( resVec[i].get_index() == index ){
			return resVec[i];
		}
	}
	return res;
}

Residue&
Residue::operator=( const Residue& res ){
	name = res.get_name();
	index = res.get_index();
	chainName = res.get_chainName();
	atomVec = res.get_atomVec();
	preResidue = res.get_preResidue();
//	cout<<"atom:"<<atomVec.size()<<endl;
	return *this;
}

void
Residue::update_atomCoord( const Atom& at ){
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_index() == at.get_index() &&
				atomVec[i].get_name() == at.get_name() ){
//			cout<<"before:";
//			atomVec[i].print();
			atomVec[i].set_coord( at.get_coord().x, at.get_coord().y, at.get_coord().z );
//			cout<<"after:";
//			atomVec[i].print();
		}
	}
}

bool
Residue::isAminoAcid(){
	string		type = name;
	bool		flag = false;
	if( type == "ALA" ){
		flag = true;
	}else if( type == "ARG" ){
		flag = true;
	}else if( type == "ASN" ){
		flag = true;
	}else if( type == "ASP" ){
		flag = true;
	}else if( type == "CYS" ){
		flag = true;
	}else if( type == "GLU" ){
		flag = true;
	}else if( type == "GLY" ){
		flag = true;
	}else if( type == "GLN" ){
		flag = true;
	}else if( type == "HIS" ){
		flag = true;
	}else if( type == "ILE" ){
		flag = true;
	}else if( type == "LEU" ){
		flag = true;
	}else if( type == "LYS" ){
		flag = true;
	}else if( type == "MET" ){
		flag = true;
	}else if( type == "PHE" ){
		flag = true;
	}else if( type == "PRO" ){
		flag = true;
	}else if( type == "SER" ){
		flag = true;
	}else if( type == "THR" ){
		flag = true;
	}else if( type == "TRP" ){
		flag = true;
	}else if( type == "TYR" ){
		flag = true;
	}else if( type == "VAL" ){
		flag = true;
	}else if( type == "MSE" ){
		flag = true;
	}
	return flag;
}

void
Residue::set_molType(){
	if( name == "DA" || name == "DT" || name == "DG"  || name == "DC" ){
		type = MOL_DNA;
	}else if(  name == "A" || name == "U" || name == "G"  || name == "C" ){
		type = MOL_RNA;
	}else if( name == "HOH" ){
		type = MOL_WATHER;
	}else if( !isAminoAcid() ){
		type = MOL_LIGAND;
	}else{
		type = MOL_PROTEIN;
	}
}
