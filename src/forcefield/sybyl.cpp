/*
 * sybyl.cpp
 *
 *  Created on: Nov 16, 2013
 *      Author: stan
 */

#include<iostream>
#include<fstream>
#include<vector>

#include"sybyl.h"
#include"../algorithm/tokenBySpace.h"
#include"../algorithm/graph.h"

using namespace std;

//////////////////////////////////////// SybylAtom  ////////////////////////////////////////////////////////////

SybylAtom::SybylAtom( const vector<string>& strVec ){
	if( strVec.size() <11 ){
		cout<<" !"<<strVec.size()<<endl;
		cerr<<" can not read atom parameter"<<endl;
		return;
	}
	type 					=			strVec[0];
	vdwRadius 				= 			atof( strVec[1].c_str() );
	taffK					= 			atof( strVec[2].c_str() );
	ionisationPotential 	= 			atof( strVec[3].c_str() );
	polarizability  		= 			atof( strVec[4].c_str() );
	weight					= 			atof( strVec[5].c_str() );
	formalCharge			= 			atof( strVec[6].c_str() );
	geometry				=			strVec[7];
	numOfNeighbors			=			atoi( strVec[8].c_str() );
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
}

void
SybylAtom::print(){
	cout.width(6);
	cout<<right<<type;
	cout.width(8);
	cout<<right<<vdwRadius;
	cout.width(8);
	cout<<right<<taffK;
	cout.width(8);
	cout<<right<<ionisationPotential;
	cout.width(8);
	cout<<right<<polarizability;
	cout.width(8);
	cout<<right<<weight;
	cout.width(8);
	cout<<right<<formalCharge;
	cout.width(8);
	cout<<right<<geometry;
	cout.width(8);
	cout<<right<<numOfNeighbors;
	cout.width(8);
	cout<<right<<isDonor;
	cout.width(8);
	cout<<right<<isAcceptor<<endl;;
}

////// SybylHydrogenBondParam /////////////////////////////////////////////

SybylHydrogenBondParam::SybylHydrogenBondParam( const vector<string>& strVec ){
	if( strVec.size() <6 ){
		cout<<strVec.size()<<endl;
		cerr<<" can not read Hydrogen Bond parameter"<<endl;
		return;
	}

	bondType 			= strVec[0];
	atomType 			= strVec[1];
	atomTypeEntry 	= strVec[2];
	donorOrAcceptor = strVec[3];

	if( strVec[3] == "M" ){
		coordinationNumber = strVec[4];
		metalDistance = atof( strVec[5].c_str() );
	}else{
		lonePairDirection = strVec[4];
	}
}

//////////////////////////////////////// SybylHBondEnergyParam //////////////////////////////////////////////////////

SybylHBondEnergyParam::SybylHBondEnergyParam( const vector<string>& strVec ){
	if( strVec.size() != 3 ){
		cerr<<" can not read Hydrogen Bond Energy parameter"<<endl;
		return;
	}

	donorType 			= strVec[0];
	acceptorType 		= strVec[1];
	maxEnergy 			= atof( strVec[2].c_str() );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SybylHydrogenBond
get_maxHB( const vector<SybylHydrogenBond>& hbVec ){
	if( hbVec.empty() ){
		throw 			" no hydrogen bond to get ";
	}

	// energy is minus, thus the smaller the better.
	float							maxEnergy = 0;
	SybylHydrogenBond	hb = hbVec[0];
	for( size_t i=0; i<hbVec.size(); i++ ){
		if( hbVec[i].energy < maxEnergy ){
			maxEnergy = hbVec[i].energy ;
			hb = hbVec[i] ;

		}
	}
	return							hb;
}

//////////////////////////////////////////////// Sybyl //////////////////////////////////////////////////////////////////

void
Sybyl::readGoldParmFile( const string& fileName ){
	ifstream		iff( fileName.c_str() );
	if( !iff ){
		cerr<<"---error can not open "<<fileName<<"---"<<endl;
	}
	string							oneLine;
	vector<string>			strVec;
	bool							startRecordAtom 	= false ;
	bool							startHBond 				= false;
	bool							startHBondEnergy	= false;

	vector<SybylAtom>		allAtom;
	allAtom.clear();

	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		if( strVec[0] == "START_ATOM_TABLE" ){
			startRecordAtom = true;
		}else if( strVec[0] == "END_ATOM_TABLE" ){
			startRecordAtom = false;
		}

		if( strVec[0] == "START_H_BOND_TABLE" ){
			startHBond = true;
		}else if(  strVec[0] =="END_H_BOND_TABLE" ){
			startHBond = false;
		}

		if( strVec[0] == "START_H_BOND_ENERGIES" ){
			startHBondEnergy 	= true;
		}else if( strVec[0] == "END_H_BOND_ENERGIES" ){
			startHBondEnergy	= false;
		}

		if( startRecordAtom ){
			if( strVec[0] != "#" && strVec.size() == 11 ){
				SybylAtom			atom( strVec );
				atomVec.push_back( atom );
			}
		}
		if( startHBond ){
			if( strVec[0] != "#" && strVec.size() > 5 ){
				SybylHydrogenBondParam 			hb( strVec );
				hydrogenBondParamVec.push_back( hb );
			}
		}
		if( startHBondEnergy ){
			if( strVec[0] != "#" && strVec.size() == 3 ){
				SybylHBondEnergyParam				hbEnergy( strVec );
				hbondEnergyParamVec.push_back( hbEnergy );
			}
		}
	}

	if(0){
		cout<<"atomVec:"<<atomVec.size()<<endl;
		cout<<"hydrogenBondParamVec:"<<hydrogenBondParamVec.size()<<endl;
		cout<<"hbondEnergyParamVec:"<<hbondEnergyParamVec.size()<<endl;
	}
	if(0){
		for( size_t i=0; i<atomVec.size(); i++  ){
			atomVec[i].print();
		}
	}
}

void
Sybyl::readFile(const string& fileName){
	readGoldParmFile( fileName );
}

vector<Atom>
Sybyl::getHydrogenBondDonorAtoms( const vector<Atom>& atVec ){
    vector<Atom>			donorVec;
	for( size_t i=0; i<atVec.size(); i++ ){
        Atom at=atVec[i];
        SybylAtom			sy1 = getSybylAtom( at );
        if( sy1.isDonor && setAttachedH(at, atVec) ){
            donorVec.push_back( at );
		}
	}
	return		donorVec;
}

vector<Atom>
Sybyl::getHydrogenBondAcceptorAtoms( const vector<Atom>& atVec ){
	vector<Atom>			acceptorVec;
	for( size_t i=0; i<atVec.size(); i++ ){
        Atom at=atVec[i];
        SybylAtom			sy1 = getSybylAtom( at );
        if( sy1.isAcceptor && setAttachedC(at, atVec) ){
            acceptorVec.push_back( at );
		}
	}
	return		acceptorVec;
}

float
Sybyl::computeGoldVDW(  const Molecular& molA, const Molecular& molB ){
	vector<Atom>			at1 = molA.get_atomVec();
	vector<Atom>			at2 = molB.get_atomVec();
	float							potential = 0;
//	cout<<"at1:"<<at1.size()<<" at2:"<<at2.size()<<endl;
	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += computeGoldVDW( at1[i], at2[j] );
		}
	}

	return potential;
}

float
Sybyl::computeGoldVDW( const vector<Atom>& at1, const vector<Atom>& at2 ){
	float							potential = 0;
//	cout<<"at1:"<<at1.size()<<" at2:"<<at2.size()<<endl;
	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += computeGoldVDW( at1[i], at2[j] );
		}
	}

	return potential;
}


float
Sybyl::computeGoldVDW(  const Atom& at1, const Atom& at2 ){
	SybylAtom			sy1 = getSybylAtom( at1 );
	SybylAtom			sy2 = getSybylAtom( at2 );

	if( sy1.vdwRadius == 0 ){
		return 0;
	}
	if( sy2.vdwRadius == 0 ){
		return 0;
	}

	float					dis = atomDis( at1, at2 );
	float					a = dis/( sy1.vdwRadius + sy2.vdwRadius );
	if( a == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float					vdwVal = 0;

	/**
	 * notice that the vdw crash is considered to be zero when they under 3 \AA
	 * this correlates to GOLD ext_vdw value
	 */
	if( a < 1.5 && dis > 3 ){
		float			a12 = a*a*a*a*a*a*a*a*a*a*a*a;
		float			a6 = a*a*a*a*a*a;
		float			k = sqrt( sy1.taffK * sy2.taffK );
		vdwVal = k * ( 1.0/a12 - 2.0/a6 );
	}

	return vdwVal;
}

float
Sybyl::computeVDW(  const Atom& at1, const Atom& at2, const float dis1 ){
	SybylAtom			sy1 = getSybylAtom( at1 );
	SybylAtom			sy2 = getSybylAtom( at2 );

	if( sy1.vdwRadius == 0 ){
		return 0;
	}
	if( sy2.vdwRadius == 0 ){
		return 0;
	}

	sy1.vdwRadius=1.52;
	sy2.vdwRadius=1.52;

//	float					dis = atomDis( at1, at2 );
	float					dis = dis1;
	float					a = dis/( sy1.vdwRadius + sy2.vdwRadius );
	if( a == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float					vdwVal = 0;

	float			a12 = a*a*a*a*a*a*a*a*a*a*a*a;
	float			a6 = a*a*a*a*a*a;
	float			k = sqrt( sy1.taffK * sy2.taffK );
	vdwVal = k * ( 1.0/a12 - 2.0/a6 );
	return vdwVal;
}

float
Sybyl::computeVDW(  const Atom& at1, const Atom& at2 ){
	SybylAtom			sy1 = getSybylAtom( at1 );
	SybylAtom			sy2 = getSybylAtom( at2 );

	if( sy1.vdwRadius == 0 ){
		return 0;
	}
	if( sy2.vdwRadius == 0 ){
		return 0;
	}

	float					dis = atomDis( at1, at2 );
	float					a = dis/( sy1.vdwRadius + sy2.vdwRadius );
	if( a == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float					vdwVal = 0;

	float			a12 = a*a*a*a*a*a*a*a*a*a*a*a;
	float			a6 = a*a*a*a*a*a;
	float			k = sqrt( sy1.taffK * sy2.taffK );
	vdwVal = k * ( 1.0/a12 - 2.0/a6 );
	return vdwVal;
}

float
Sybyl::computeVDW( const vector<Atom>& at1, const vector<Atom>& at2 ){
	float							potential = 0;
	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += computeVDW( at1[i], at2[j] );
		}
	}
	return potential;
}

//float
//Sybyl::computeVDW( const Molecular mol ){
//	vector<Atom> atVec=mol.get_atomVec();
//	float sum=0;
//	for( size_t i=0; i<atVec.size()-1; i++ ){
//		for( size_t j=i+1; j<atVec.size(); j++ ){
//			Atom a1=atVec[i];
//			Atom a2=atVec[j];
//			float dis=atomDis(a1,a2);
//			if( dis>5.0 ){
//				float val=computeVDW(a1, a2);
//				sum+=val;
//			}
//		}
//	}
//	return sum;
//}

float
Sybyl::computeVDW( const Molecular mol ){
    vector<Atom> atVec = mol.get_atomVec();
    vector<Bond> bdVec = mol.get_bondVec();
    float sum=0;

    vector<vector<size_t> > matrix = connectMatrix(atVec, bdVec);

    for( size_t i=0; i<atVec.size()-1; i++ ){
        for( size_t j=i+1; j<atVec.size(); j++ ){
            Atom a1=atVec[i];
            Atom a2=atVec[j];
            if( matrix[i][j]>3 ){
                float val=computeVDW(a1, a2);
                sum+=val;
            }else if( matrix[i][j] == 3 ){
                float val=computeVDW(a1, a2);
                sum+=val/2.0 ;
            }
        }
    }
    return sum;
}

float
Sybyl::computeCutVDW(  const Atom& at1, const Atom& at2 ){
	SybylAtom			sy1 = getSybylAtom( at1 );
	SybylAtom			sy2 = getSybylAtom( at2 );

	if( sy1.vdwRadius == 0 ){
		return 0;
	}
	if( sy2.vdwRadius == 0 ){
		return 0;
	}

	float					dis = atomDis( at1, at2 );

	float					a = dis/( sy1.vdwRadius + sy2.vdwRadius );
	if( a == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		at1.print();
		at2.print();
		return 0;
	}

	float					vdwVal = 0;

	float			a12 = a*a*a*a*a*a*a*a*a*a*a*a;
	float			a6 = a*a*a*a*a*a;
	float			k = sqrt( sy1.taffK * sy2.taffK );
	vdwVal = k * ( 1.0/a12 - 2.0/a6 );

	//cut maximum energy
//	vdwVal = vdwVal>5 ? 5:vdwVal;
	vdwVal = vdwVal>1 ? 1:vdwVal;

	return vdwVal;
}

float
Sybyl::computeCutVDW( const vector<Atom>& at1, const vector<Atom>& at2, bool computeHydrogen ){
	float							potential = 0;
	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
            if( computeHydrogen ){
                potential += computeCutVDW( at1[i], at2[j] );
            }else if( at1[i]. get_type() != "H" && at2[j]. get_type() != "H" ){
                potential += computeCutVDW( at1[i], at2[j] );
            }
		}
	}
	return potential;
}



/**
 * electrostatic potential
 */
float
Sybyl::computeESP( const Atom a1, const Atom a2){
	float val=0;
	float p1=a1.get_charge();
	float p2=a2.get_charge();
	if( p1>1 || p1<-1 || p2 >1 || p2<-1 ){
		val=0;
		return val;
	}
	float eps=1;  //ε0=8.85*10^(-12)F/m
	float dis=atomDis(a1,a2);
	if( dis==0 ){
		cout<<"can not compute electrostatic potential for the same atom"<<endl;
		val=0;
		return val;
	}
	val= p1*p2/dis;
//	a1.print();
//	cout<<p1<<endl;
//	a2.print();
//	cout<<p2<<endl;
//	cout<<"  val:"<<val<<endl<<endl;
	return val;
}

float
Sybyl::computeESP( const vector<Atom> v1, const vector<Atom> v2 ){
	float val=0;
	for( size_t i=0; i<v1.size(); i++ ){
		for( size_t j=0; j<v2.size(); j++ ){
			float d=computeESP( v1[i], v2[j] );
			val+=d;
		}
	}
	return val;
}

SybylAtom
Sybyl::getSybylAtom( const Atom& at ){
	if( atomVec.empty() ){
		cerr<<"can not find sybyl force filed "<<endl;
	}
	string		residueName = at.get_residueName().substr( 0, 3 );
	if(  residueName == "ARG" && at.get_type() == "N.pl3" ){
		for( size_t i=0; i<atomVec.size(); i++ ){
			if( atomVec[i].type == "N.plc" ){
				return			atomVec[i];
			}
		}
	}else{
		for( size_t i=0; i<atomVec.size(); i++ ){
			if( atomVec[i].type == at.get_type() ){
				return atomVec[i];
			}
		}
	}

    // S.O2
    for( size_t i=0; i<atomVec.size(); i++ ){
        string tp1=atomVec[i].type.substr(0, atomVec[i].type.find(".") );
        string tp2=at.get_type().substr(0, at.get_type().find(".") );
        if( tp1 == tp2 ){
            return atomVec[i];
        }
    }

    for( size_t i=0; i<atomVec.size(); i++ ){
        string tp1=atomVec[i].type.substr(0, atomVec[i].type.find(".") );
        string tp2=at.get_type().substr(0, 1 );
        if( tp1 == tp2 ){
            return atomVec[i];
        }
    }

    // in case no atom matches
    return atomVec[0];
}

float
Sybyl::computeHBondEnergy( const Atom& at1, const Atom& at2 ){
	SybylAtom			syAt1 = getSybylAtom( at1 );
	SybylAtom			syAt2 = getSybylAtom( at2 );
	if( ( !syAt1.isAcceptor && ! syAt1.isDonor ) || ( !syAt2.isAcceptor && ! syAt2.isDonor ) ){
		return 0;
	}
	if( atomDis( at1, at2 ) > 3.5 ){
		return 0;
	}

	float			hbondEnergy = 0;
	string	bondType1 = "non";
	string	bondType2 = "non";
	if( syAt1.isAcceptor && syAt2.isDonor ){
		bondType1 = getHbondType( syAt1.type, false );
		bondType2 = getHbondType( syAt2.type, true );

		hbondEnergy = getHbondEnergy( bondType1, bondType2 );
	}else if( syAt1.isDonor && syAt2.isAcceptor ){
		bondType1 = getHbondType( syAt1.type, true );
		bondType2 = getHbondType( syAt2.type, false );

		hbondEnergy = getHbondEnergy( bondType1, bondType2 );
	}
	if( 1 /*hbondEnergy != 0 */){
		cout.width(5);
		cout<<left<<at1.get_index();
		cout.width(5);
		cout<<left<<at1.get_name();
		cout.width(5);
		cout<<left<<bondType1;
		cout<<" : ";
		cout.width(5);
		cout<<left<<at2.get_index();
		cout.width(5);
		cout<<left<<at2.get_name();
		cout.width(5);
		cout<<left<<bondType2;
		cout<<" ->";
		cout<<hbondEnergy<<endl;
	}

	return			hbondEnergy;
}

float
Sybyl::computeHBondEnergy(  const Atom& acceptor, const Atom & lonePair,
							const Atom& donor, const Atom& hydrogen	){
	float			maxDisH_LP = 2.5;
	float			disH_LP = atomDis( lonePair, hydrogen );
//	cout<<"disH_LP:"<<disH_LP<<endl;

	float			dis = atomDis( acceptor, hydrogen );
	float			disAcceptorDonor = atomDis( acceptor, donor );

	float			distance_wt = 0;

	SybylAtom			syAt1 = getSybylAtom( acceptor );
	SybylAtom			syAt2 = getSybylAtom( donor );

	//----------------------------------------------------------distance weight-------------------------------------------------------------------------------
	const	float		minDis 		= 0;							// default 0.5 in file gold.param, 0.2 in original Gold Paper, 1.2
	const	float		maxDis 		= 2 ;							// default 2.0 in file gold.param, 3.5 in original Gold paper, 2.5
	const	float		minDis2 		= 2.376 ;
	const	float		maxDis2		= 3.89 ;

	const 	float		maxWeight 	= 0;
	const	float		wt_multiply	= 1.0;

	if( disAcceptorDonor < 1.0 || disH_LP > maxDisH_LP ){
		distance_wt = 0;
	}else	if( dis > maxDis2 || dis < minDis  ){
		distance_wt = 0;
	}else if( dis > maxDis && dis < minDis2 ){
		distance_wt = 1;
	}else if ( dis > minDis && dis < maxDis ){
//		distance_wt =  2 * ( dis - minDis )*( dis - minDis ) / ( ( maxDis - minDis )* ( maxDis - minDis ) ) ;
		distance_wt = -0.0155 + 0.1173 * dis + 0.4216 * dis * dis ;
	}else if( dis > minDis2 && dis < maxDis2 ){
//		distance_wt =  ( (float)( dis - maxDis2 )*( dis - maxDis2 ) ) / (float)( ( maxDis2 - minDis2 )* ( maxDis2 - minDis2 ) );
		distance_wt = 6.5864 - 3.3809 * dis + 0.4338 * dis * dis ;
		if( distance_wt < 0 ){
			distance_wt = -distance_wt;
		}
	}

	if( distance_wt < 0.01 ){
		distance_wt = 0;
	}else if( distance_wt > 1.0 ){
		distance_wt = 1.0;
	}

	//-------------------------------------------------------------angle weight---------------------------------------------------------------------------------

	float			angle = atomAngle( lonePair, acceptor,  hydrogen, donor );
	float			angle_wt = 0;

	float			minAngle = 60;							// 60 in Gold file
	float			maxAngle = 160;							// 160 in Gold file

	if( angle > maxAngle ){
		angle_wt = 1;
	}else if( angle < minAngle ){
		angle_wt = 0;
	}else{
//		angle_wt = ( ( ( angle - minAngle )*( angle - minAngle ) ) / ( ( maxAngle - minAngle )*( maxAngle - minAngle ) ))  ;
		angle_wt = ( ( angle - minAngle ) / ( maxAngle - minAngle )  )  ;
	}

	if( angle_wt < 0.01 ){
		angle_wt = 0;
	}

	//------------------------------------------------------------------------  angle donor-hydrogen-lonePair ----------------------------------------------------------------------------------
	float				angleDonorHydrogenLP = atomAngle( donor, hydrogen, lonePair );
	float				angleWT_donorHydrogenLP = 0;
	float				minAngle_donorHydrogenLP = 120 ;
	float				maxAngle_donorHydrogenLP = 140 ;

	if( angleDonorHydrogenLP > maxAngle_donorHydrogenLP  ){
		angleWT_donorHydrogenLP = 1;
	}else if( angleDonorHydrogenLP  <  minAngle_donorHydrogenLP ){
		angleWT_donorHydrogenLP = 0;
	}else {
		angleWT_donorHydrogenLP = ( ( angleDonorHydrogenLP -  minAngle_donorHydrogenLP ) * ( angleDonorHydrogenLP - minAngle_donorHydrogenLP ) ) /
				(  ( maxAngle_donorHydrogenLP - minAngle_donorHydrogenLP )* ( maxAngle_donorHydrogenLP - minAngle_donorHydrogenLP )  );
	}

	if( angleWT_donorHydrogenLP < 0.01 ){
		angleWT_donorHydrogenLP = 0;
	}

//	cout.width( 25 );
//	cout<<left<<"angleDonorHydrogenLP: ";
//	cout.width( 10 );
//	cout<<left<<angleDonorHydrogenLP<<"    angleWT_donorHydrogenLP : "<<angleWT_donorHydrogenLP <<endl;

	//--------------------------------------------------------------------  angle donor-hydrogen-acceptor  --------------------------------------------------------------------------------------
	float				angle_donorHydrogenAcceptor = atomAngle( donor, hydrogen, acceptor );
	float				angleWT_donorHydrogenAcceptor = 0;
	float				minAngle_donorHydrogenAcceptor = 140 ;
	float				maxAngle_donorHydrogenAcceptor = 170 ;

	if( angle_donorHydrogenAcceptor > maxAngle_donorHydrogenAcceptor ){
		angleWT_donorHydrogenAcceptor = 1 ;
	}else if( angle_donorHydrogenAcceptor < minAngle_donorHydrogenAcceptor ){
		angleWT_donorHydrogenAcceptor = 0 ;
	}else {
//		angleWT_donorHydrogenAcceptor = (  ( angle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor ) *  ( angle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor )  ) /
//				(  ( maxAngle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor ) *  ( maxAngle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor ) );

		angleWT_donorHydrogenAcceptor = (  ( angle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor ) ) / (  ( maxAngle_donorHydrogenAcceptor - minAngle_donorHydrogenAcceptor )   ) ;
		if( angleWT_donorHydrogenAcceptor < 0 ){
			angleWT_donorHydrogenAcceptor = - angleWT_donorHydrogenAcceptor ;
		}
	}

	if( angleWT_donorHydrogenAcceptor < 0.01 ){
		angleWT_donorHydrogenAcceptor = 0;
	}

	//--------------------------------------------------------------------  angle hydrogen-lonePair-acceptor ---------------------------------------------------------------------------------

	float				angle_hydrogenLpAcceptor = atomAngle( hydrogen, lonePair, acceptor ) ;
	float				angleWT_hydrogenLpAcceptor = 0;
	float				minAngle_hydrogenLpAcceptor = 40;
	float				maxAngle_hydrogenLpAcceptor = 170;

	if( angle_hydrogenLpAcceptor >  maxAngle_hydrogenLpAcceptor ){
		angleWT_hydrogenLpAcceptor = 1 ;
	}else if( angle_hydrogenLpAcceptor <  minAngle_hydrogenLpAcceptor ){
		angleWT_hydrogenLpAcceptor = 0 ;
	}else{
		angleWT_hydrogenLpAcceptor = ( ( angle_hydrogenLpAcceptor - minAngle_hydrogenLpAcceptor ) *  ( angle_hydrogenLpAcceptor - minAngle_hydrogenLpAcceptor ) )/
				(  (maxAngle_hydrogenLpAcceptor - minAngle_hydrogenLpAcceptor  )*(maxAngle_hydrogenLpAcceptor - minAngle_hydrogenLpAcceptor ) );
	}

//	if( angleWT_hydrogenLpAcceptor < 0.01 ){
//		 angleWT_hydrogenLpAcceptor  = 0 ;
//	}

//	cout.width(25);
//	cout<<"angle_hydrogenLpAcceptor: ";
//	cout.width(10);
//	cout<<angle_hydrogenLpAcceptor<<"   angleWT_hydrogenLpAcceptor"<<angleWT_hydrogenLpAcceptor <<endl;

	//------------------------------------------------------------------------  angle donor-lonePair-acceptor  ---------------------------------------------------------------------------------
	float				angle_donorLpAcceptor = atomAngle( donor, lonePair, acceptor ) ;
	float				angleWT_donorLpAcceptor = 0;
	float				minAngle_donorLpAcceptor = 40;
	float				maxAngle_donorLpAcceptor = 170;

	if( angle_donorLpAcceptor >  maxAngle_donorLpAcceptor ){
		angleWT_donorLpAcceptor = 1 ;
	}else if( angle_donorLpAcceptor <  minAngle_donorLpAcceptor ){
		angleWT_donorLpAcceptor = 0 ;
	}else{
		angleWT_donorLpAcceptor = ( ( angle_donorLpAcceptor - minAngle_donorLpAcceptor ) *  ( angle_donorLpAcceptor - minAngle_donorLpAcceptor ) )/
				(  (maxAngle_donorLpAcceptor - minAngle_donorLpAcceptor  )*(maxAngle_donorLpAcceptor - minAngle_donorLpAcceptor ) );
	}

	//----------------------------------------------------------------------- overall weight -------------------------------------------------------------------------
//	float			wt = ( distance_wt * angleWT_donorHydrogenLP * angleWT_hydrogenLpAcceptor );
//	float			wt = ( distance_wt * angle_wt * angleWT_hydrogenLpAcceptor );
	float			wt = ( distance_wt * angle_wt * angleWT_donorHydrogenAcceptor );

//	cout<<"dis_wt:"<<distance_wt<<" angle_wt:"<<angle_wt<<" angleWT_d:"<<angleWT_donorHydrogenLP<<endl;
	if( wt < maxWeight ){
		return 0;
	}

	if( ( !syAt1.isAcceptor && ! syAt1.isDonor ) || ( !syAt2.isAcceptor && ! syAt2.isDonor ) ){
		cout<<" no hbond formed"<<endl;
		return 0;
	}

	float			hbondEnergy = 0;
	string			bondType1 = "non";
	string			bondType2 = "non";
	float			stdBondEnergy = 0;

	if( syAt1.isAcceptor && syAt2.isDonor ){
		bondType1 = getHbondType( syAt1.type, false );
		bondType2 = getHbondType( syAt2.type, true );
		hbondEnergy = getHbondEnergy( bondType1, bondType2 ) * wt ;
		stdBondEnergy = getHbondEnergy( bondType1, bondType2 );
	}else if( syAt1.isDonor && syAt2.isAcceptor ){
		bondType1 = getHbondType( syAt1.type, true );
		bondType2 = getHbondType( syAt2.type, false );
		hbondEnergy = getHbondEnergy( bondType1, bondType2 ) * wt ;
		stdBondEnergy = getHbondEnergy( bondType1, bondType2 );
	}

	SybylHydrogenBond				sHB( syAt1, syAt2, hydrogen, lonePair,  dis, angle,  hbondEnergy );
	hbVec.push_back( sHB );

	if( wt != 0 ){
		cout.width(50);
		cout<<"";
		cout<<"angle:";
		cout.width(10);
		cout<<angle;
		cout<<"angle_donorHyLp:";
		cout.width(10);
		cout<<angleDonorHydrogenLP;
		cout<<"angle_hyLpAccept:";
		cout.width(10 );
		cout<<angle_hydrogenLpAcceptor;
		cout<<"angle_donorHydrogenAcceptor:";
		cout.width( 10 );
		cout<<angle_donorHydrogenAcceptor;
		cout<<"dis:";
		cout.width(10);
		cout<<dis;
		cout<<endl;
	}

	if( 0 ){
		cout.width(45);
		cout<<"";
		cout<<"Acceptor:";
		cout.width(5);
		cout<<left<<acceptor.get_index();
		cout.width(4);
		cout<<left<<acceptor.get_name();
		cout.width(5);
		cout<<left<<acceptor.get_type()<<" ";
		cout<<"LonePair:";
		cout.width(5);
		cout<<left<<lonePair.get_index()<<" ";
		cout.width(5);
		cout<<left<<lonePair.get_name();
		cout<<"Donor:";
		cout.width(5);
		cout<<left<<donor.get_index()<<" ";
		cout.width(4);
		cout<<left<<donor.get_name() ;
		cout.width(5);
		cout<<left<<donor.get_type()<<" " ;
//		cout.width(10);
//		cout<<left<<donor.get_residueName();
		cout<<"Hydrogen:";
		cout.width(5);
		cout<<left<<hydrogen.get_index()<<" ";
		cout.width(5);
		cout<<left<<hydrogen.get_name();
		cout<<"stdE: ";
		cout.width(3);
		cout<<stdBondEnergy;
		cout<<"energy:";
		cout.width(10);
		cout<<hbondEnergy;
		cout<<" dis:";
		cout.width(10);
		cout<<left<<dis;
		cout<<"dis_wt:";
		cout.width(10);
		cout<<left<<distance_wt;
		cout<<"angle:";
		cout.width(10);
		cout<<left<<angle;
		cout<<"ang_wt:";
		cout.width(10);
		cout<<left<<angle_wt;
		cout<<"wt:"<<wt<<endl;
//		cout<<"wt:"<<wt<<" = "<<distance_wt<<" + "<<angle_wt<<endl;
	}
	return			hbondEnergy;
}

float
Sybyl::computeHBondEnergy2(  const Molecular& molA, const Molecular& molB ){
	vector<Atom>			at1 			= molA.get_atomVec();
	vector<Atom>			at2 			= molB.get_atomVec();
	float							potential 	= 0;
	float							dis 			= 0;
	const	float				maxDonorAcceptorDis = 6.0;
	hbVec.clear();

	for( size_t i=0; i<at1.size(); i++ ){
		vector<float>						allHBforOneAtom;
		allHBforOneAtom.clear();
		for( size_t j=0; j<at2.size(); j++ ){
			SybylAtom			syAt1 = getSybylAtom( at1[i] );
			SybylAtom			syAt2 = getSybylAtom( at2[j] );
			dis = atomDis( at1[i], at2[j] );

			if( syAt1.isAcceptor && syAt2.isDonor && dis < maxDonorAcceptorDis ){

				vector<Atom>			acceptorLP = getAcceptorLonePair( molA, at1[i].get_index() );
				vector<Atom>			donorHydrogen = getDonorHydrogen( molB, at2[j].get_index() );
				if( !acceptorLP.empty() && !donorHydrogen.empty() ){

					Atom							closestHydrogen = findClosestAtom( at1[i], donorHydrogen );
					Atom							closestLonePair 	= findClosestAtom( closestHydrogen, acceptorLP );
					potential += computeHBondEnergy( at1[i],  closestLonePair, at2[j], closestHydrogen );
				}
			}else if( syAt2.isAcceptor && syAt1.isDonor  && dis < maxDonorAcceptorDis ){

				vector<Atom>			acceptorLP = getAcceptorLonePair( molB, at2[j].get_index() );
				vector<Atom>			donorHydrogen = getDonorHydrogen( molA, at1[i].get_index() );
				if( !acceptorLP.empty() && !donorHydrogen.empty() ){

					Atom							closestHydrogen = findClosestAtom( at2[j], donorHydrogen );
					Atom							closestLonePair 	= findClosestAtom( closestHydrogen, acceptorLP );
					potential += computeHBondEnergy( at2[j], closestLonePair, at1[i], closestHydrogen ) ;
				}
			}
		}
	}
	return potential;
}

float
Sybyl::computeHBondEnergy(  const Molecular& molA, const Molecular& molB ){
	vector<Atom>			at1 			= molA.get_atomVec();
	vector<Atom>			at2 			= molB.get_atomVec();
	float							potential 	= 0;
	float							dis 			= 0;
	const	float				maxDonorAcceptorDis = 6.0;
	hbVec.clear();

	for( size_t i=0; i<at1.size(); i++ ){
		float										maxHbEnergyOneAtom = 0;
		for( size_t j=0; j<at2.size(); j++ ){
			SybylAtom			syAt1 = getSybylAtom( at1[i] );
			SybylAtom			syAt2 = getSybylAtom( at2[j] );
			dis = atomDis( at1[i], at2[j] );
			if( syAt1.isAcceptor && syAt2.isDonor && dis < maxDonorAcceptorDis ){
				vector<Atom>   acceptorLP = getAcceptorLonePair( molA, at1[i].get_index() );
				vector<Atom>   donorHydrogen = getDonorHydrogen( molB, at2[j].get_index() );
				if( !acceptorLP.empty() && !donorHydrogen.empty() ){
					Atom       closestHydrogen = findClosestAtom( at1[i], donorHydrogen );
					Atom       closestLonePair = findClosestAtom( closestHydrogen, acceptorLP );
					//notice that the energy is minus than 0, such that -6 will be accepted than -3;
					float      energy = computeHBondEnergy( at1[i],  closestLonePair, at2[j], closestHydrogen );
					maxHbEnergyOneAtom = energy < maxHbEnergyOneAtom ? energy : maxHbEnergyOneAtom;
				}
			}else if( syAt2.isAcceptor && syAt1.isDonor  && dis < maxDonorAcceptorDis ){
				vector<Atom>   acceptorLP = getAcceptorLonePair( molB, at2[j].get_index() );
				vector<Atom>   donorHydrogen = getDonorHydrogen( molA, at1[i].get_index() );
				if( !acceptorLP.empty() && !donorHydrogen.empty() ){
					Atom       closestHydrogen = findClosestAtom( at2[j], donorHydrogen );
					Atom       closestLonePair = findClosestAtom( closestHydrogen, acceptorLP );
					float      energy = computeHBondEnergy( at2[j], closestLonePair, at1[i], closestHydrogen ) ;
					//notice that the energy is minus than 0, such that -6 will be accepted than -3;
					maxHbEnergyOneAtom = energy < maxHbEnergyOneAtom ? energy : maxHbEnergyOneAtom;
				}
			}
		}
		potential += maxHbEnergyOneAtom ;
	}
	return potential;
}

Coord
Sybyl::getIdealHBondAtom( const Atom& donor, const Atom& donorH, const string& acceptorName ){
	Coord dCoord = donor.get_coord();
	Coord hCoord = donorH.get_coord();

	float dhDis = getCoordDis( dCoord, hCoord );
	Coord acceptor;
	string type = acceptorName.substr(0,1);

	if( type == "N" ){
		float ha = 2.8;
		float dx = dCoord.x-hCoord.x;
		float dy = dCoord.y-hCoord.y;
		float dz = dCoord.z-hCoord.z;

		float x = hCoord.x+(ha/dhDis)*(hCoord.x-dCoord.x);
		float y = hCoord.y+(ha/dhDis)*(hCoord.y-dCoord.y);
		float z = hCoord.z+(ha/dhDis)*(hCoord.z-dCoord.z);

		Coord c(x, y, z);
		acceptor = c;
	}else if( type == "O" ){
		float ha = 2.8;
		float dx = dCoord.x-hCoord.x;
		float dy = dCoord.y-hCoord.y;
		float dz = dCoord.z-hCoord.z;

		float x = hCoord.x+(ha/dhDis)*(hCoord.x-dCoord.x);
		float y = hCoord.y+(ha/dhDis)*(hCoord.y-dCoord.y);
		float z = hCoord.z+(ha/dhDis)*(hCoord.z-dCoord.z);

		Coord c(x, y, z);
		acceptor = c;
	}
	return acceptor;
}

/////////////////////////////////////////////////////////////////////////////////////

vector<Atom>
Sybyl::getAcceptorLonePair(  const Molecular& mol, const size_t& acceptorIndex ){
//	SybylAtom			at = getSybylAtom( mol.get_atomVec()[ acceptorIndex ] );
	SybylAtom			at = getSybylAtom( searchAtom( mol.get_atomVec(), acceptorIndex ) );
	if( !at.isAcceptor ){
		cout<<"error in getting acceptor lone pair"<<endl;;
		throw ;
	}
	vector<Atom>	lpVec;
	lpVec.clear();
	vector<Bond>		allBond = mol.get_bondVec() ;
	vector<Atom>	allAtom = mol.get_atomVec();

	//check if atom $index$ has already lone pairs denoted as $LP$
	for( size_t i=0; i<allBond.size(); i++ ){
		if( allBond[i].get_firstAtomIndex() == acceptorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_secAtomIndex() );
			if( connectedAt.get_type() == "LP" ){
				lpVec.push_back( connectedAt );
			}
		}else	if( allBond[i].get_secAtomIndex() == acceptorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_firstAtomIndex() );
			if( connectedAt.get_type() == "LP" ){
				lpVec.push_back( connectedAt );
			}
		}
	}
	if( ! lpVec.empty() ){
		return lpVec;
	}
//	cerr<<"no original lone pair from acceptor atom:"<<acceptorIndex<<endl;
}

vector<Atom>
Sybyl::getDonorHydrogen( const Molecular& mol, const size_t& donorIndex ){

	SybylAtom				at = getSybylAtom( searchAtom( mol.get_atomVec(), donorIndex ) );
	if( ! at.isDonor ){
		Atom			tempAt =  searchAtom( mol.get_atomVec(), donorIndex ) ;
		cout<<"atom:"<<tempAt.get_index()<<" "<<tempAt.get_name()<<endl;
		cout<<"!!! error in getting donor hydrogen atom"<<donorIndex<<endl;;
		throw ;
	}

	vector<Atom>	hydrogenVec;
	hydrogenVec.clear();
	vector<Bond>		allBond = mol.get_bondVec() ;
	vector<Atom>	allAtom = mol.get_atomVec();

	for( size_t i=0; i<allBond.size(); i++ ){
		if( allBond[i].get_firstAtomIndex() == donorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_secAtomIndex() );
			if( connectedAt.get_type() == "H" ){
				hydrogenVec.push_back( connectedAt );
			}
		}else	if( allBond[i].get_secAtomIndex() == donorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_firstAtomIndex() );
			if( connectedAt.get_type() == "H" ){
				hydrogenVec.push_back( connectedAt );
			}
		}
	}

	//---------------------------------------------------take lone pair in TYR  as hydrogen, as they are rotable in geometry--------------------------------------------------------------
	for( size_t i=0; i<allBond.size(); i++ ){
		if( allBond[i].get_firstAtomIndex() == donorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_secAtomIndex() );
			if( connectedAt.get_type() == "LP" ){
				if(connectedAt.get_residueName().substr(0, 3) == "TYR" ||
						connectedAt.get_residueName().substr(0, 3) == "THR" ||
						connectedAt.get_residueName().substr(0, 3) == "SER" ){
					hydrogenVec.push_back( connectedAt );
				}
			}
		}else	if( allBond[i].get_secAtomIndex() == donorIndex ){
			Atom  connectedAt = searchAtom( allAtom, allBond[i].get_firstAtomIndex() );
			if( connectedAt.get_type() == "LP"){
				if( connectedAt.get_residueName().substr(0, 3) == "TYR" ||
						connectedAt.get_residueName().substr(0, 3) == "THR" ||
						connectedAt.get_residueName().substr(0, 3) == "SER"			){
					hydrogenVec.push_back( connectedAt );
				}
			}
		}
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------

	if( ! hydrogenVec.empty() ){
		return hydrogenVec;
	}

//	cerr<<"no original hydrogen atom to donor atom:"<<donorIndex<<endl;
}

string
Sybyl::getHbondType( const string& atType, const bool& isDonor ){

	string			donor;
	if( isDonor ){
		donor = "D";
	}else{
		donor = "A";
	}

	for( size_t i=0; i<hydrogenBondParamVec.size(); i++ ){
		if( hydrogenBondParamVec[i].atomType == atType &&
			( hydrogenBondParamVec[i].donorOrAcceptor == "DA" ||  hydrogenBondParamVec[i].donorOrAcceptor == donor  )	){
			return hydrogenBondParamVec[i].bondType;
		}
	}
	cerr<<" can not find hydrogen bond type "<<hydrogenBondParamVec.size()<<endl;
}

float
Sybyl::getHbondEnergy( const string& bondTypeA, const string& bondTypeB ){
	float			energy = 0;
	string			bondType1 = bondTypeA;
	string			bondType2 = bondTypeB;
//	if( bondType1 == "NARA" ){
//		bondType1 = "N2A";
//	}else if( bondType2 == "NARA" ){
//		bondType2 = "N2A";
//	}
	bondType1 = ( bondTypeA == "NARA" ? "N2A" : bondTypeA );
	bondType2 = ( bondTypeB == "NARA" ? "N2A" : bondTypeB );

	for( size_t i=0; i<hbondEnergyParamVec.size(); i++ ){
		if( ( bondType1 == hbondEnergyParamVec[i].donorType && bondType2 == hbondEnergyParamVec[i].acceptorType ) ||
				( bondType2 == hbondEnergyParamVec[i].donorType && bondType1 == hbondEnergyParamVec[i].acceptorType  )){
			return hbondEnergyParamVec[i].maxEnergy;
		}
	}
	return 0;
}

void
Sybyl::initCrashDis(){
	crashDisVec.clear();
	string disVec[][3]={
			{"N","N","2.76179"},
			{"N","C","2.89542"},
			{"N","O","2.73506"},
			{"N","S1","2.89542"}, // S.o S.o2
			{"N","S2","2.98452"}, // S.3 S.2 S.m S.a
			{"N","P","2.98452"},
			{"N","H","2.71724"},
			{"N","F","2.69052"},
			{"N","Cl","2.93997"},
			{"N","Br","3.02907"},
			{"N","I","3.14487"},
			{"N","Si","2.44998"},
			{"N","Na","2.44998"},
			{"N","K","2.44998"},
			{"N","Ca","2.36979"},
			{"N","Li","2.44998"},
			{"N","Al","2.44998"},
			{"N","Mg","2.04016"},
			{"N","Zn","2.12034"},
			{"N","Fe","2.12034"},
			{"N","Mn","2.17380"},
			{"C","C","3.02907"},
			{"C","O","2.86870"},
			{"C","S1","3.02907"},
			{"C","S2","3.11815"},
			{"C","P","3.11815"},
			{"C","H","2.85088"},
			{"C","F","2.82416"},
			{"C","Cl","3.07361"},
			{"C","Br","3.16270"},
			{"C","I","3.27851"},
			{"C","Si","2.58361"},
			{"C","Na","2.58361"},
			{"C","K","2.58361"},
			{"C","Ca","2.50343"},
			{"C","Li","2.58361"},
			{"C","Al","2.58361"},
			{"C","Mg","2.17380"},
			{"C","Zn","2.25397"},
			{"C","Fe","2.25397"},
			{"C","Mn","2.30744"},
			{"O","O","2.70834"},
			{"O","S1","2.86870"},
			{"O","S2","2.95779"},
			{"O","P","2.95779"},
			{"O","H","2.69051"},
			{"O","F","2.66379"},
			{"O","Cl","2.91325"},
			{"O","Br","3.00233"},
			{"O","I","3.11815"},
			{"O","Si","2.42325"},
			{"O","Na","2.42325"},
			{"O","K","2.42325"},
			{"O","Ca","2.34307"},
			{"O","Li","2.42325"},
			{"O","Al","2.42325"},
			{"O","Mg","2.01344"},
			{"O","Zn","2.09362"},
			{"O","Fe","2.09362"},
			{"O","Mn","2.14707"},
			{"S1","S1","3.02907"},
			{"S1","S2","3.11815"},
			{"S1","P","3.11815"},
			{"S1","H","2.85088"},
			{"S1","F","2.82416"},
			{"S1","Cl","3.07361"},
			{"S1","Br","3.16270"},
			{"S1","I","3.27851"},
			{"S1","Si","2.58361"},
			{"S1","Na","2.58361"},
			{"S1","K","2.58361"},
			{"S1","Ca","2.50343"},
			{"S1","Li","2.58361"},
			{"S1","Al","2.58361"},
			{"S1","Mg","2.17380"},
			{"S1","Zn","2.25397"},
			{"S1","Fe","2.25397"},
			{"S1","Mn","2.30744"},
			{"S2","S2","3.20724"},
			{"S2","P","3.20724"},
			{"S2","H","2.93997"},
			{"S2","F","2.91324"},
			{"S2","Cl","3.16270"},
			{"S2","Br","3.25179"},
			{"S2","I","3.36761"},
			{"S2","Si","2.67270"},
			{"S2","Na","2.67270"},
			{"S2","K","2.67270"},
			{"S2","Ca","2.59252"},
			{"S2","Li","2.67270"},
			{"S2","Al","2.67270"},
			{"S2","Mg","2.26289"},
			{"S2","Zn","2.34307"},
			{"S2","Fe","2.34307"},
			{"S2","Mn","2.39652"},
			{"P","P","3.20724"},
			{"P","H","2.93997"},
			{"P","F","2.91324"},
			{"P","Cl","3.16270"},
			{"P","Br","3.25179"},
			{"P","I","3.36761"},
			{"P","Si","2.67270"},
			{"P","Na","2.67270"},
			{"P","K","2.67270"},
			{"P","Ca","2.59252"},
			{"P","Li","2.67270"},
			{"P","Al","2.67270"},
			{"P","Mg","2.26289"},
			{"P","Zn","2.34307"},
			{"P","Fe","2.34307"},
			{"P","Mn","2.39652"},
			{"H","H","2.67270"},
			{"H","F","2.64598"},
			{"H","Cl","2.89542"},
			{"H","Br","2.98452"},
			{"H","I","3.10033"},
			{"H","Si","2.40543"},
			{"H","Na","2.40543"},
			{"H","K","2.40543"},
			{"H","Ca","2.32525"},
			{"H","Li","2.40543"},
			{"H","Al","2.40543"},
			{"H","Mg","1.99562"},
			{"H","Zn","2.07580"},
			{"H","Fe","2.07580"},
			{"H","Mn","2.12925"},
			{"F","F","2.61925"},
			{"F","Cl","2.86870"},
			{"F","Br","2.95779"},
			{"F","I","3.07361"},
			{"F","Si","2.37870"},
			{"F","Na","2.37870"},
			{"F","K","2.37870"},
			{"F","Ca","2.29852"},
			{"F","Li","2.37870"},
			{"F","Al","2.37870"},
			{"F","Mg","1.96890"},
			{"F","Zn","2.04907"},
			{"F","Fe","2.04907"},
			{"F","Mn","2.10253"},
			{"Cl","Cl","3.11815"},
			{"Cl","Br","3.20724"},
			{"Cl","I","3.32306"},
			{"Cl","Si","2.62815"},
			{"Cl","Na","2.62815"},
			{"Cl","K","2.62815"},
			{"Cl","Ca","2.54797"},
			{"Cl","Li","2.62815"},
			{"Cl","Al","2.62815"},
			{"Cl","Mg","2.21835"},
			{"Cl","Zn","2.29852"},
			{"Cl","Fe","2.29852"},
			{"Cl","Mn","2.35198"},
			{"Br","Br","3.29633"},
			{"Br","I","3.41215"},
			{"Br","Si","2.71724"},
			{"Br","Na","2.71724"},
			{"Br","K","2.71724"},
			{"Br","Ca","2.63706"},
			{"Br","Li","2.71724"},
			{"Br","Al","2.71724"},
			{"Br","Mg","2.30744"},
			{"Br","Zn","2.38762"},
			{"Br","Fe","2.38762"},
			{"Br","Mn","2.44107"},
			{"I","I","3.52797"},
			{"I","Si","2.83306"},
			{"I","Na","2.83306"},
			{"I","K","2.83306"},
			{"I","Ca","2.75288"},
			{"I","Li","2.83306"},
			{"I","Al","2.83306"},
			{"I","Mg","2.42325"},
			{"I","Zn","2.50343"},
			{"I","Fe","2.50343"},
			{"I","Mn","2.55689"},
			{"Si","Si","2.13816"},
			{"Si","Na","2.13816"},
			{"Si","K","2.13816"},
			{"Si","Ca","2.05798"},
			{"Si","Li","2.13816"},
			{"Si","Al","2.13816"},
			{"Si","Mg","1.72835"},
			{"Si","Zn","1.80853"},
			{"Si","Fe","1.80853"},
			{"Si","Mn","1.86198"},
			{"Ca","Ca","1.97780"},
			{"Ca","Li","2.13816"},
			{"Ca","Al","2.13816"},
			{"Ca","Mg","1.64817"},
			{"Ca","Zn","1.72835"},
			{"Ca","Fe","1.72835"},
			{"Ca","Mn","1.78180"},
			{"Mg","Mg","1.31853"},
			{"Mg","Zn","1.39872"},
			{"Mg","Fe","1.39872"},
			{"Mg","Mn","1.45217"},
			{"Zn","Zn","1.47890"},
			{"Zn","Fe","1.47890"},
			{"Zn","Mn","1.53235"},
			{"Mn","Mn","1.58580"}
	};
	crashDisVec.clear();
	for( int i=0; i<201; i++ ){
		string n1=disVec[i][0];
		string n2=disVec[i][1];
		float  d=atof(disVec[i][2].c_str());
		VdwCrash vc(n1, n2, d);
		crashDisVec.push_back(vc);
	}

	if(0){
		for( size_t i=0; i<crashDisVec.size(); i++ ){
			cout<<crashDisVec[i].n1<<" "<<crashDisVec[i].n2<<" "<<crashDisVec[i].dis<<endl;
		}
	}
}

/**
 * atoms crash if their vdw energy is larger than zero
 */
bool
Sybyl::atomsCrash( const vector<Atom> atVec1, const vector<Atom> atVec2 ){
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			if( atomsCrash(atVec1[i], atVec2[j]) ){
				return true;
			}
		}
	}
	return false;
}

bool
Sybyl::atomsCrash( const  Atom at, const vector<Atom> atVec ){

	for( size_t i=0; i<atVec.size(); i++ ){
		if( atomsCrash(at, atVec[i]) ){
			return true;
		}
	}

	return false;
}

bool
Sybyl::atomsCrash( const Atom at1, const Atom at2 ){
	SybylAtom syA1=getSybylAtom(at1);
	SybylAtom syA2=getSybylAtom(at2);

	string t1=at1.get_type();
	string t2=at2.get_type();
	if( t1==string("S.o") || t1==string("S.o2") ){
		t1=string("S1");
	}
	if( t1==string("S.3") || t1==string("S.2") ||
			t1==string("S.m")||t1==string("S.a") ){
		t1=string("S2");
	}
	if( t2==string("S.o") || t2==string("S.o2") ){
		t2=string("S1");
	}
	if( t2==string("S.3") || t2==string("S.2") ||
			t2==string("S.m")||t2==string("S.a") ){
		t2=string("S2");
	}
	if( t1!=string("S1") && t1!=string("S2") ){
		t1=t1.substr( 0, t1.find(".") );
	}
	if( t2!=string("S1") && t2!=string("S2") ){
		t2=t2.substr( 0, t2.find(".") );
	}

	float cutDis=0;

	for( size_t i=0; i<crashDisVec.size(); i++ ){
		if( ( t1==crashDisVec[i].n1 && t2==crashDisVec[i].n2 )||
				( t1==crashDisVec[i].n2 && t2==crashDisVec[i].n1 )){
			cutDis= crashDisVec[i].dis;
			break;
		}
	}
	float dis = atomDis(at1, at2);

	if( dis<cutDis ){
		return true;
	}

	return false;
}
