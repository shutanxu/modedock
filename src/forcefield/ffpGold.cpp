/*
 * ffpGold.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: stan
 */

#include<fstream>
#include<iostream>

#include"ffpGold.h"

using namespace std;

void
GOLD_ATOM_PARAM::print()const{
	cout.width(25);
	cout<<left<<"type:"<<type<<endl;
	cout.width(25);
	cout<<left<<"vdwRadius:"<<vdwRadius<<endl;
	cout.width(25);
	cout<<left<<"taffK:"<<taffK<<endl;
	cout.width(25);
	cout<<left<<"ionisationPotential:"<<ionisationPotential<<endl;
	cout.width(25);
	cout<<left<<"polarizability:"<<polarizability<<endl;
	cout.width(25);
	cout<<left<<"weight:"<<weight<<endl;
	cout.width(25);
	cout<<left<<"formalCharge:"<<formalCharge<<endl;
	cout.width(25);
	cout<<left<<"geometry:"<<geometry<<endl;
	cout.width(25);
	cout<<left<<"numOfNeighbors:"<<numOfNeighbors<<endl;
	cout.width(25);
	cout<<left<<"isDonor:"<<isDonor<<endl;
	cout.width(25);
	cout<<left<<"isAcceptor:"<<isAcceptor<<endl;
}

GOLD_PARAMETER
readGoldParamFile( const string& fileName ){
	ifstream		iff( fileName.c_str() );
	if( !iff ){
		cerr<<"---error can not open "<<fileName<<"---"<<endl;
	}
	string							oneLine;
	vector<string>			strVec;
	bool							startRecord = false ;

	string					type;
	float					vdwRadius;
	float					taffK;
	float					ionisationPotential;
	float					polarizability;
	float					weight;
	float					formalCharge;
	string					geometry;
	size_t					numOfNeighbors;
	bool					isDonor = false;
	bool					isAcceptor = false;
	vector<GOLD_ATOM_PARAM>		allGold;
	allGold.clear();

	while( getline( iff, oneLine ) ){
		if( oneLine == "START_ATOM_TABLE" ){
			startRecord = true;
		}else if( oneLine == "END_ATOM_TABLE" ){
			startRecord = false;
			break;
		}

		if( startRecord ){
			strVec = tokenBySpace( oneLine );
			if( strVec[0] != "#" && strVec.size() == 11 ){
				type 						=			strVec[0];
				vdwRadius 				= 			atof( strVec[1].c_str() );
				taffK							= 			atof( strVec[2].c_str() );
				ionisationPotential 	= 			atof( strVec[3].c_str() );
				polarizability  			= 			atof( strVec[4].c_str() );
				weight						= 			atof( strVec[5].c_str() );
				formalCharge			= 			atof( strVec[6].c_str() );
				geometry					=			strVec[7];
				numOfNeighbors		=			atoi( strVec[8].c_str() );
				if( strVec[9] == "N" ){
					isDonor = false;
				}else if( strVec[9] == "Y" ){
					isDonor = true;
				}
				if( strVec[10] == "N" ){
					isAcceptor = false;
				}else if( strVec[10] == "Y" ){
					isAcceptor = true;
				}

				GOLD_ATOM_PARAM 		goldAtom( type, vdwRadius, taffK, ionisationPotential, polarizability,
						weight, formalCharge, geometry,  numOfNeighbors, isDonor, isAcceptor );
				allGold.push_back( goldAtom );
			}
		}
	}
	GOLD_PARAMETER			goldP;
	goldP.atomParVec = allGold;
	return	goldP;
}

float
goldVDW( const Atom& at1, const Atom& at2, const GOLD_PARAMETER& gp ){
	GOLD_ATOM_PARAM	atPar1;
	GOLD_ATOM_PARAM	atPar2;
	for( size_t i=0; i<gp.atomParVec.size(); i++ ){
		if( at1.get_type() == gp.atomParVec[i].type ){
			atPar1 = gp.atomParVec[i];
		}
		if( at2.get_type() == gp.atomParVec[i].type ){
			atPar2 = gp.atomParVec[i];
		}
	}

	if( atPar1.vdwRadius == 0 || atPar2.vdwRadius == 0 ){
		return	0;
	}

	float		dis = atomDis( at1, at2 );
	float		a = dis / ( atPar1.vdwRadius + atPar2.vdwRadius );
	if( a == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float		vdwVal = 0;
	/**
	 * notice that the vdw crash is considered to be zero when they under 3 \AA
	 * this correlates to GOLD ext_vdw value
	 */
	if( a < 1.5 && dis > GOLD_VDW_MINDIS ){
		float			a12 = a*a*a*a*a*a*a*a*a*a*a*a;
		float			a6 = a*a*a*a*a*a;
		float			k = sqrt( atPar1.taffK * atPar2.taffK );
		vdwVal = k * ( 1.0/a12 - 2.0/a6 );
	}
	return				vdwVal;
}

float
goldVDW( const Molecular& mol1, const Molecular& mol2, const GOLD_PARAMETER& gp ){
	vector<Atom>			at1 = mol1.get_atomVec();
	vector<Atom>			at2 = mol2.get_atomVec();
	float							potential = 0;

	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += goldVDW ( at1[i], at2[j], gp );

			if(0){
				cout<<"atom:  ";
				cout.width(4);
				cout<<left<<i;
				cout.width(8);
				cout<<left<<j<<"  poteintial:"<<potential<<endl;
			}
		}
	}
//	cout<<"potential:"<<potential<<endl;
	return potential;
}
