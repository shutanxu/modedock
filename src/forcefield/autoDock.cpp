/*
 * autoDock.cpp
 *
 *  Created on: Jan 3, 2014
 *      Author: stan
 */


#include"autoDock.h"
#include"../algorithm/tokenBySpace.h"

//static AD_ATOM_PARM AD_ATOM_PARM_NULL;

using namespace std;

static AD_ATOM_PARM AD_ATOM_PARM_NULL;

AD_ATOM_PARM::AD_ATOM_PARM( const vector<string>& strVec ){
	atomType 	= strVec[1];
	Rii 				=	atof( strVec[2].c_str() );
	espii				=	atof( strVec[3].c_str() );
	vol				=	atof( strVec[4].c_str() );
	solpar			=	atof( strVec[5].c_str() );
	Rij_hb			=	atof( strVec[6].c_str() );
	espij_hb		=	atof( strVec[7].c_str() );
	hbond			=	atof( strVec[8].c_str() );
	rec_index		=	atoi( strVec[9].c_str() );
	map_index	=	atoi( strVec[10].c_str() );
	bond_index	=	atoi( strVec[11].c_str() );
}

void AD_PARM::readParamFile(  const string& fileName ){
	ifstream			iif( fileName.c_str() );
	if( ! iif ){
		cerr<<"can not open param file "<<fileName<<endl;
		throw ;
	}

	string					oneLine;
	vector<string>	strVec;
	while( getline( iif, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		if( strVec[0] != "#" ){
			if( strVec[0] == "FE_coeff_vdW" ){
				FE_coeff_vdw = atof( strVec[1].c_str() );
			}else if( strVec[0] == "FE_coeff_hbond" ){
				FE_coeff_hbond = atof( strVec[1].c_str() );
			}else if( strVec[0] == "FE_coeff_estat" ){
				FE_coeff_estat = atof( strVec[1].c_str() );
			}else if( strVec[0] == "FE_coeff_desolv" ){
				FE_coeff_desolv = atof( strVec[1].c_str() );
			}else if( strVec[0] == "FE_coeff_tors" ){
				FE_coeff_tors = atof( strVec[1].c_str() );
			}else if( strVec[0] == "atom_par" ){
				AD_ATOM_PARM				atParm( strVec );
				atomParmVec.push_back( atParm );
			}
		}
	}
}

//------------------------------------------ AutoDock -------------------------------------------------------------

AD_ATOM_PARM
AutoDock::getAtomParm( const Atom& at ){
	for( size_t i=0; i<atom_parm.atomParmVec.size(); i++ ){
		if( atom_parm.atomParmVec[i].atomType == at.get_type() ){
			return			atom_parm.atomParmVec[i];
		}
	}
	return AD_ATOM_PARM_NULL;
}


float	AutoDock::computeVDW_cut( const Atom& at1, const Atom& at2 ){

	AD_ATOM_PARM				autoAtom1 = getAtomParm( at1 );
	AD_ATOM_PARM				autoAtom2 = getAtomParm( at2 );

	if( &autoAtom1 == &AD_ATOM_PARM_NULL  ||    &autoAtom2 == &AD_ATOM_PARM_NULL){
		return 			0 ;
	}

	if( autoAtom1.Rii == 0 ){
//		cerr<<"no appropriate force field for atom "<<at1.get_index()<<" "<<at1.get_name()<<" "<<at1.get_type()<<endl;
		return				0;
	}
	if( autoAtom2.Rii == 0 ){
//		cerr<<"no appropriate force field for atom "<<at2.get_index()<<" "<<at2.get_name()<<" "<<at2.get_type()<<endl;
		return				0;
	}

	float								dis = atomDis( at1, at2 );
	if( dis == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float			maxDis = 8.0 ;
	if( dis >= maxDis  ){
		return 0 ;
	}

	//cut --------------------------
	if( dis <( ( autoAtom1.Rii + autoAtom2.Rii ) / 2.0 ) ){
//		cout<<"dis:"<<dis<<endl;
//		cout<<"vdwR1:"<<autoAtom1.Rii <<"  vdwR2:"<<autoAtom2.Rii<<endl;
		return 0;
	}

	float								r12 = pow( dis, 12 );
	float								r6 = pow( dis, 6 );
//	cout<<"dis:"<<dis<<" r12:"<<r12<<" r6:"<<r6<<endl;

//	float								epsi_12 = sqrt( autoAtom1.espii * autoAtom2.espii );
//	epsi_12 = epsi_12 * atom_parm.FE_coeff_vdw ;
	float									epsi_12 = sqrt( autoAtom1.espii * atom_parm.FE_coeff_vdw * autoAtom2.espii * atom_parm.FE_coeff_vdw );


	float								r_eqmXY = ( autoAtom1.Rii + autoAtom2.Rii ) / 2.0 ;

	float								A_ij = 		epsi_12 * pow( r_eqmXY, 12 );
	float								B_ij = 2 * epsi_12 * pow( r_eqmXY, 6 );

	if( 0 ){
		cout<<"at1:"<<at1.get_name()<<" at2:"<<at2.get_name();
		dis = 4.0;
		cout<<" espi_12:"<<epsi_12<<endl<<endl;
	}

	float								energy = 0 ;
//	energy = ( A_ij / r12 ) - ( B_ij / r6 ) ;

//from AutoDock tutorial Version 3.0.3 p16 V(r)
	energy = ( 0.5 * epsi_12 * pow(r_eqmXY, 12) / r12 ) - ( 2.0 * epsi_12 * pow(r_eqmXY, 6) / r6 ) ;

	if(0){
		cout.width(4);
		cout<<left<<autoAtom1.atomType;
		cout.width(4);
		cout<<autoAtom2.atomType;
		cout.width(18);
		cout<< 1.0 * epsi_12 * pow(r_eqmXY, 12);
		cout<<"r_eqmXY12:";
		cout.width(18);
		cout<<pow(r_eqmXY, 12);;
		cout.width(18);
		cout<<2.0 * epsi_12 * pow(r_eqmXY, 6) ;
		cout<<"r_eqmXY6:";
		cout.width(18);
		cout<<pow(r_eqmXY, 6);
		cout<<endl;
	}
	return								energy ;
}

float	AutoDock::computeVDW( const Atom& at1, const Atom& at2 ){

	AD_ATOM_PARM				autoAtom1 = getAtomParm( at1 );
	AD_ATOM_PARM				autoAtom2 = getAtomParm( at2 );

	if( &autoAtom1 == &AD_ATOM_PARM_NULL  ||    &autoAtom2 == &AD_ATOM_PARM_NULL){
		return 			0 ;
	}

	if( autoAtom1.Rii == 0 ){
//		cerr<<"no appropriate force field for atom "<<at1.get_index()<<" "<<at1.get_name()<<" "<<at1.get_type()<<endl;
		return				0;
	}
	if( autoAtom2.Rii == 0 ){
//		cerr<<"no appropriate force field for atom "<<at2.get_index()<<" "<<at2.get_name()<<" "<<at2.get_type()<<endl;
		return				0;
	}

	float								dis = atomDis( at1, at2 );
	if( dis == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float			maxDis = 8.0 ;
	if( dis >= maxDis  ){
		return 0 ;
	}

	float								r12 = pow( dis, 12 );
	float								r6 = pow( dis, 6 );
//	cout<<"dis:"<<dis<<" r12:"<<r12<<" r6:"<<r6<<endl;

//	float								epsi_12 = sqrt( autoAtom1.espii * autoAtom2.espii );
//	epsi_12 = epsi_12 * atom_parm.FE_coeff_vdw ;
	float									epsi_12 = sqrt( autoAtom1.espii * atom_parm.FE_coeff_vdw * autoAtom2.espii * atom_parm.FE_coeff_vdw );


	float								r_eqmXY = ( autoAtom1.Rii + autoAtom2.Rii ) / 2.0 ;

	float								A_ij = 		epsi_12 * pow( r_eqmXY, 12 );
	float								B_ij = 2 * epsi_12 * pow( r_eqmXY, 6 );

	if( 0 ){
		cout<<"at1:"<<at1.get_name()<<" at2:"<<at2.get_name();
		dis = 4.0;
		cout<<" espi_12:"<<epsi_12<<endl<<endl;
	}

	float								energy = 0 ;
//	energy = ( A_ij / r12 ) - ( B_ij / r6 ) ;

//from AutoDock tutorial Version 3.0.3 p16 V(r)
	energy = ( 0.5 * epsi_12 * pow(r_eqmXY, 12) / r12 ) - ( 2.0 * epsi_12 * pow(r_eqmXY, 6) / r6 ) ;

	if(0){
		cout.width(4);
		cout<<left<<autoAtom1.atomType;
		cout.width(4);
		cout<<autoAtom2.atomType;
		cout.width(18);
		cout<< 1.0 * epsi_12 * pow(r_eqmXY, 12);
		cout<<"r_eqmXY12:";
		cout.width(18);
		cout<<pow(r_eqmXY, 12);;
		cout.width(18);
		cout<<2.0 * epsi_12 * pow(r_eqmXY, 6) ;
		cout<<"r_eqmXY6:";
		cout.width(18);
		cout<<pow(r_eqmXY, 6);
		cout<<endl;
	}

	return								energy ;
}

float
AutoDock::computeVDW( const Molecular& molA, const Molecular& molB ){
	vector<Atom>			at1 = molA.get_atomVec();
	vector<Atom>			at2 = molB.get_atomVec();
	float							potential = 0;

	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += computeVDW( at1[i], at2[j] );
		}
	}

	return potential;
}

float
AutoDock::computeVDW(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 ){
	float				overall = 0;
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			overall += computeVDW( atVec1[i], atVec2[j] );
		}
	}
	return				overall;
}

float
AutoDock::computeVDW_cut(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 ){
	float				overall = 0;
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			overall += computeVDW_cut( atVec1[i], atVec2[j] );
		}
	}
	return				overall;
}

float
AutoDock::computeHBondEnergy(  const Atom& at1, const Atom& at2 ){
	AD_ATOM_PARM				autoAtom1 = getAtomParm( at1 );
	AD_ATOM_PARM				autoAtom2 = getAtomParm( at2 );

	if( autoAtom1.Rij_hb == 0 || autoAtom2.Rij_hb == 0 ){
		return				0;
	}

	float								dis = atomDis( at1, at2 );
	if( dis == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float								r12 = pow( dis, 12 );
	float								r10 = pow( dis, 10 );

	float								epsi_12 = sqrt( autoAtom1.espij_hb * autoAtom2.espij_hb );
	float								r_eqmXY = ( autoAtom1.Rii + autoAtom2.Rii ) / 2.0 ;

	float								C_ij = 5 * epsi_12 * pow( r_eqmXY, 12 );
	float								D_ij = 6 * epsi_12 * pow( r_eqmXY, 10 );

	float								energy = ( C_ij / r12 ) - ( D_ij / r10 ) ;
	return								energy ;
}

float
AutoDock::computeHBondEnergy( const Molecular& molA, const Molecular& molB ){
	vector<Atom>			at1 = molA.get_atomVec();
	vector<Atom>			at2 = molB.get_atomVec();
	float							potential = 0;

	for( size_t i=0; i<at1.size(); i++ ){
		for( size_t j=0; j<at2.size(); j++ ){
			potential += computeHBondEnergy( at1[i], at2[j] );
		}
	}

	return potential;
}

float
AutoDock::computeHBondEnergy(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 ){
	float				overall = 0;
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			overall += computeHBondEnergy( atVec1[i], atVec2[j] );
		}
	}
	return				overall;
}

float
AutoDock::computeElectrostatic(  const Atom& at1, const Atom& at2 ){
	AD_ATOM_PARM				autoAtom1 = getAtomParm( at1 );
	AD_ATOM_PARM				autoAtom2 = getAtomParm( at2 );
//	cout<<"at1:"<<at1.get_charge()<<" " << at2.get_charge() <<endl;
	if( at1.get_charge() == 0 || at2.get_charge() == 0 ){
		return 			0;
	}

	float								dis = atomDis( at1, at2 );
	if( dis == 0 ){
		cerr<<"can not compute potential between atom "
				<<at1.get_index()<<" "<<at2.get_index()<<"  distance is zero!!!"<<endl;
		return 0;
	}

	float								epsi = 0;
	float								A = -8.5525 ;
	float								B = 78.4 - A ;
	float								k = 7.7839 ;
	float								lamd = 0.003627 ;
	epsi = A + B / ( 1 + k * exp( - lamd * B * dis ) ) ;

	float								potential = 0;
	if( epsi != 0 ){
		potential = ( at1.get_charge() * at2.get_charge() ) / ( epsi * dis ) ;
	}
	return 	potential ;
}

float
AutoDock::computeElectrostatic(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 ){
	float				overall = 0;
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			overall += computeElectrostatic( atVec1[i], atVec2[j] );
		}
	}
	return				overall;
}

//------------------------------------------------------------------------------------------------------------------------
vector<AutoDock_Result>
readFile_dlg(  const string& fileName ){
	ifstream			iif( fileName.c_str() );
	if( ! iif ){
		cerr<<"can not open param file "<<fileName<<endl;
		throw ;
	}

	vector<AutoDock_Result>			ligandResultVec;

	string					oneLine;
	vector<string>	strVec;
	size_t					count = 0;

	size_t							ligand_index =0 ;
	Atom							atom;
	float							bindEnergy = 0;
	float							final_inter = 0;
	float							vdw_hbond_desolv = 0;
	float							electrostatic = 0;
	float							total_internal = 0;
	float							torsional = 0;
	float							unbound = 0;
	vector<AutoDock_ATOM>			atVec ;

	while( getline( iif, oneLine ) ){
		strVec = tokenBySpace( oneLine );

		if( strVec[0] == "DOCKED:" ){
			if( strVec[1] == "MODEL" ){
				ligand_index = atoi( strVec[2].c_str() ) ;
				ligand_index =0 ;
				bindEnergy = 0;
				final_inter = 0;
				vdw_hbond_desolv = 0;
				electrostatic = 0;
				total_internal = 0;
				torsional = 0;
				unbound = 0;
				atVec.clear();
			}
			if( strVec.size() > 2 && strVec[1] == "USER" ){
				if( strVec[3] == "Free" ){
					bindEnergy = atof( strVec[8].c_str() ) ;
				}else if( strVec[2] == "(1)" ){
					final_inter = atof( strVec[7].c_str() );
				}else	if( strVec[2] == "vdW" ){
					vdw_hbond_desolv = atof( strVec[9].c_str() );
				}else	if( strVec[2] == "Electrostatic" ){
					electrostatic = atof( strVec[5].c_str() );
				}else if( strVec[2] == "(2)" ){
					total_internal = atof( strVec[8].c_str() );
				}else if( strVec[2] == "(3)" ){
					torsional = atof( strVec[7].c_str() );
				}else if( strVec[2] == "(4)" ){
					unbound = atof( strVec[8].c_str() );
				}
			}
			if(  strVec[1] == "ATOM" ){
				size_t				index = atoi( strVec[2].c_str() );
				string				name = strVec[3];
				Coord				coord;
				coord.x = atof( strVec[ strVec.size() - 7].c_str() );
				coord.y = atof( strVec[ strVec.size() - 6].c_str() );
				coord.z = atof( strVec[ strVec.size() - 5 ].c_str() );
				float				vdw = atof( strVec[ strVec.size() - 4 ].c_str() );
				float				elec = atof( strVec[  strVec.size() - 3 ].c_str() );
				float				q = atof( strVec[ strVec.size() - 2 ].c_str() );
				string				type =  strVec[ strVec.size() - 1 ].c_str() ;

				AutoDock_ATOM	atom( index, name, type, coord, vdw, elec, q );
				atVec.push_back( atom );
			}
			if( strVec[1] == "ENDMDL" ){
				AutoDock_Result				ad( ligand_index, bindEnergy, final_inter, vdw_hbond_desolv, electrostatic, total_internal,
						torsional, unbound, atVec );
				ligandResultVec.push_back( ad );

			}
		}
	}
	return			 ligandResultVec ;
}
