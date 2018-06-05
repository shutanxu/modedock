/*
 * goldLigand.cpp
 *
 *  Created on: Sep 22, 2013
 *      Author: stan
 */

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include"goldLigand.h"
#include"../algorithm/clusterdata.h"
#include"../algorithm/tokenBySpace.h"
//#include"../clusterdata.hpp"

using namespace std;

GoldLigand::~GoldLigand(){;}

int
GoldLigand::readGoldScore( const string& fileName ){
	ifstream		inStream( fileName.c_str() );
	string			oneLine;
	while( getline( inStream, oneLine ) ){
		if( oneLine == "> <Gold.Goldscore.Fitness>" ){
			getline( inStream, oneLine );
			fitness = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.Vdw>" ){
			getline( inStream, oneLine );
			inter_vdw = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.Torsion>" ){
			getline( inStream, oneLine );
			 inter_tors= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.External.Vdw>" ){
			getline( inStream, oneLine );
			 ext_vdw = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.HBond>" ){
			getline( inStream, oneLine );
			 inter_hb = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.External.HBond>" ){
			getline( inStream, oneLine );
			 ext_hb= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.Correction>" ){
			getline( inStream, oneLine );
			 inter_corre= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Reference.RMSD>" ){
			getline( inStream, oneLine );
			 rmsd_ref= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Score>" ){
			getline( inStream, oneLine );
			getline( inStream, oneLine );
			vector<string>		tokenVec = tokenBySpace( oneLine );
			inter_energy = atof( tokenVec[4].c_str() );
		}
	}
	inStream.close();
	return 1;
}

//------------------------------------------------------------------------------------------------

int
printGoldResult( const vector<GoldLigand>& goldLigandVec ){
	ofstream		oof( "bestranking.lst" );
	oof<<"#     Fitness  S(hb_ext) S(vdw_ext)  S(hb_int)    S(int)      intcor RMS(heavy) time                               File name                Ligand name"<<endl;
	for( size_t i=0; i<goldLigandVec.size(); i++ ){
		oof.width(10);
		oof<<right<<goldLigandVec[i].get_fitness();
		oof.width(11);
		oof<<goldLigandVec[i].get_ext_hb() ;
		oof.width(11);
		oof<<goldLigandVec[i].get_ext_vdw() ;
		oof.width(11);
		oof<<goldLigandVec[i].get_inter_hb() ;
		oof.width(11);
		oof<<goldLigandVec[i].get_inter_energy();
		oof.width(11);
		oof<<goldLigandVec[i].get_inter_corre() ;
		oof.width(9);
		oof<<goldLigandVec[i].get_rmsd_ref() ;
		oof.width(11);
		oof<<0;
		oof<<"    ";
//		oof.width(40);
		oof<<"'"<<goldLigandVec[i].get_fileName()<<"'" ;
		oof<<"   ";
		oof<<"****";
		oof<<endl;
	}
	oof.close();
}

int
printGoldHeader(){
	cout.width(5);
	cout<<"Index";
	cout.width(15);
	cout<<"fitness";
	cout.width(15);
	cout<<"inter_vdw";
	cout.width(15);
	cout<<"inter_tors";
	cout.width(15);
	cout<<"inter_hb";
	cout.width(15);
	cout<<"inter_coore";
	cout.width(15);
	cout<<"ext_hb";
	cout.width(15);
	cout<<"ext_vdw";
	cout.width(15);
	cout<<"rmsd_ref";
	cout.width(15);
	cout<<"pairDis";

	cout<<endl;
}

int
print( const vector<GoldLigand>& goldLigandVec,
		const 	vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix ){

	for( size_t i=0; i<clusterVec.size(); i++ ){
		printGoldHeader();
		vector<float>	minScoreVec(8, 1e10 );
		vector<float>	maxScoreVec(8, -1e10);
		vector<float>	scoreVec( 8, 0 );
		cout<<"---cluster "<<i+1<<"------ number of structures: "<<clusterVec[i].size()<<endl;
		for( size_t j=0; j<clusterVec[i].size(); j++ ){
			cout.width(5);
			cout<<clusterVec[i][j];

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_fitness();
			minScoreVec[0] = minScoreVec[0] < goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					minScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			maxScoreVec[0] = maxScoreVec[0] > goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					maxScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			scoreVec[0] += goldLigandVec[ clusterVec[i][j] ].get_fitness();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			minScoreVec[1] = minScoreVec[1] < goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					minScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			maxScoreVec[1] = maxScoreVec[1] > goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					maxScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			scoreVec[1] += goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			minScoreVec[2] = minScoreVec[2] < goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					minScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			maxScoreVec[2] = maxScoreVec[2] > goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					maxScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			scoreVec[2] += goldLigandVec[ clusterVec[i][j] ].get_inter_tors();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			minScoreVec[3] = minScoreVec[3] < goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					minScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			maxScoreVec[3] = maxScoreVec[3] > goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					maxScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			scoreVec[3] += goldLigandVec[ clusterVec[i][j] ].get_inter_hb();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			minScoreVec[4] = minScoreVec[4] < goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					minScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			maxScoreVec[4] = maxScoreVec[4] > goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					maxScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			scoreVec[4] += goldLigandVec[ clusterVec[i][j] ].get_inter_corre();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			minScoreVec[5] = minScoreVec[5] < goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					minScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			maxScoreVec[5] = maxScoreVec[5] > goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					maxScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			scoreVec[5] += goldLigandVec[ clusterVec[i][j] ].get_ext_hb();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			minScoreVec[6] = minScoreVec[6] < goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					minScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			maxScoreVec[6] = maxScoreVec[6] > goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					maxScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			scoreVec[6] += goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();

			cout.width(15);
			cout<< goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			minScoreVec[7] = minScoreVec[7] < goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					minScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			maxScoreVec[7] = maxScoreVec[7] > goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					maxScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			scoreVec[7] += goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();

			cout<<endl;
		}
		cout<<endl;
		cout.width(5);
		cout<<"min";
		for( size_t k=0; k<8; k++ ){
			cout.width(15);
			cout<<minScoreVec[k];
		}
		cout.width(15);
		cout<<getMinPair( distMatrix, clusterVec[i] );
		cout<<endl;

		cout.width(5);
		cout<<"max";
		for( size_t k=0; k<8; k++ ){
			cout.width(15);
			cout<<maxScoreVec[k];
		}
		cout.width(15);
		float	maxPair = getMaxPair( distMatrix, clusterVec[i] );
		cout<<maxPair;
		cout<<endl;

		cout.width(5);
		cout<<"aver";
		for( size_t k=0; k<8; k++ ){
			cout.width(15);
			if( clusterVec[i].size() != 0 ){
				cout<<scoreVec[k]/ clusterVec[i].size();
			}else{
				cout<<scoreVec[k] ;
			}
		}
		cout.width(15);
		cout<<getAveragePair( distMatrix, clusterVec[i] );
		cout<<endl<<endl;
	}
}

int
printGoldHeader(  ofstream& oof  ){
//	ofstream	oof = off;
	oof.width(5);
	oof<<"Index";
	oof.width(15);
	oof<<"fitness";
	oof.width(15);
	oof<<"inter_vdw";
	oof.width(15);
	oof<<"inter_tors";
	oof.width(15);
	oof<<"inter_hb";
	oof.width(15);
	oof<<"inter_coore";
	oof.width(15);
	oof<<"ext_hb";
	oof.width(15);
	oof<<"ext_vdw";
	oof.width(15);
	oof<<"rmsd_ref";
	oof.width(15);
	oof<<"overlay_ref";
	oof.width(15);
	oof<<"pairDis";

	oof<<endl;
}

int
print( const vector<GoldLigand>& goldLigandVec,
		const 	vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix,
		const string& fileName ){

	ofstream		off( fileName.c_str() );

	for( size_t i=0; i<clusterVec.size(); i++ ){

		vector<float>	minScoreVec(8, 1e10 );
		vector<float>	maxScoreVec(8, -1e10);
		vector<float>	scoreVec( 8, 0 );
		off<<"---cluster "<<i+1<<"------ number of structures: "<<clusterVec[i].size()<<endl;
		printGoldHeader(off);
		for( size_t j=0; j<clusterVec[i].size(); j++ ){
			off.width(5);
			off<<clusterVec[i][j];

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_fitness();
			minScoreVec[0] = minScoreVec[0] < goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					minScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			maxScoreVec[0] = maxScoreVec[0] > goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					maxScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			scoreVec[0] += goldLigandVec[ clusterVec[i][j] ].get_fitness();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			minScoreVec[1] = minScoreVec[1] < goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					minScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			maxScoreVec[1] = maxScoreVec[1] > goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					maxScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			scoreVec[1] += goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			minScoreVec[2] = minScoreVec[2] < goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					minScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			maxScoreVec[2] = maxScoreVec[2] > goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					maxScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			scoreVec[2] += goldLigandVec[ clusterVec[i][j] ].get_inter_tors();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			minScoreVec[3] = minScoreVec[3] < goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					minScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			maxScoreVec[3] = maxScoreVec[3] > goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					maxScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			scoreVec[3] += goldLigandVec[ clusterVec[i][j] ].get_inter_hb();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			minScoreVec[4] = minScoreVec[4] < goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					minScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			maxScoreVec[4] = maxScoreVec[4] > goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					maxScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			scoreVec[4] += goldLigandVec[ clusterVec[i][j] ].get_inter_corre();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			minScoreVec[5] = minScoreVec[5] < goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					minScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			maxScoreVec[5] = maxScoreVec[5] > goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					maxScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			scoreVec[5] += goldLigandVec[ clusterVec[i][j] ].get_ext_hb();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			minScoreVec[6] = minScoreVec[6] < goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					minScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			maxScoreVec[6] = maxScoreVec[6] > goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					maxScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			scoreVec[6] += goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			minScoreVec[7] = minScoreVec[7] < goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					minScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			maxScoreVec[7] = maxScoreVec[7] > goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					maxScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			scoreVec[7] += goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();

			off<<endl;
		}
		off<<endl;
		off.width(5);
		off<<"min";
		for( size_t k=0; k<8; k++ ){
			off.width(15);
			off<<minScoreVec[k];
		}
		off.width(15);
		off<<getMinPair( distMatrix, clusterVec[i] );
		off<<endl;

		off.width(5);
		off<<"max";
		for( size_t k=0; k<8; k++ ){
			off.width(15);
			off<<maxScoreVec[k];
		}
		off.width(15);
		float	maxPair = getMaxPair( distMatrix, clusterVec[i] );
		off<<maxPair;
		off<<endl;

		off.width(5);
		off<<"aver";
		for( size_t k=0; k<8; k++ ){
			off.width(15);
			if( clusterVec[i].size() != 0 ){
				off<<scoreVec[k]/ clusterVec[i].size();
			}else{
				off<<scoreVec[k] ;
			}
		}
		off.width(15);
		off<<getAveragePair( distMatrix, clusterVec[i] );
		off<<endl<<endl;
	}
	off.close();

}

int
print( const vector<GoldLigand>& goldLigandVec,
		const 	vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix,
		const	vector<float>& overlayDisVec,
		const string& fileName ){

	ofstream		off( fileName.c_str() );

	for( size_t i=0; i<clusterVec.size(); i++ ){

		vector<float>	minScoreVec(9, 1e10 );
		vector<float>	maxScoreVec(9, -1e10);
		vector<float>	scoreVec( 9, 0 );
		off<<"---cluster "<<i+1<<"------ number of structures: "<<clusterVec[i].size()<<endl;
		printGoldHeader(off);
		for( size_t j=0; j<clusterVec[i].size(); j++ ){
			off.width(5);
			off<<clusterVec[i][j];

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_fitness();
			minScoreVec[0] = minScoreVec[0] < goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					minScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			maxScoreVec[0] = maxScoreVec[0] > goldLigandVec[ clusterVec[i][j] ].get_fitness() ?
					maxScoreVec[0] : goldLigandVec[ clusterVec[i][j] ].get_fitness();
			scoreVec[0] += goldLigandVec[ clusterVec[i][j] ].get_fitness();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			minScoreVec[1] = minScoreVec[1] < goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					minScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			maxScoreVec[1] = maxScoreVec[1] > goldLigandVec[ clusterVec[i][j] ].get_inter_vdw() ?
					maxScoreVec[1] : goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();
			scoreVec[1] += goldLigandVec[ clusterVec[i][j] ].get_inter_vdw();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			minScoreVec[2] = minScoreVec[2] < goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					minScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			maxScoreVec[2] = maxScoreVec[2] > goldLigandVec[ clusterVec[i][j] ].get_inter_tors() ?
					maxScoreVec[2] : goldLigandVec[ clusterVec[i][j] ].get_inter_tors();
			scoreVec[2] += goldLigandVec[ clusterVec[i][j] ].get_inter_tors();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			minScoreVec[3] = minScoreVec[3] < goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					minScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			maxScoreVec[3] = maxScoreVec[3] > goldLigandVec[ clusterVec[i][j] ].get_inter_hb() ?
					maxScoreVec[3] : goldLigandVec[ clusterVec[i][j] ].get_inter_hb();
			scoreVec[3] += goldLigandVec[ clusterVec[i][j] ].get_inter_hb();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			minScoreVec[4] = minScoreVec[4] < goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					minScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			maxScoreVec[4] = maxScoreVec[4] > goldLigandVec[ clusterVec[i][j] ].get_inter_corre() ?
					maxScoreVec[4] : goldLigandVec[ clusterVec[i][j] ].get_inter_corre();
			scoreVec[4] += goldLigandVec[ clusterVec[i][j] ].get_inter_corre();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			minScoreVec[5] = minScoreVec[5] < goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					minScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			maxScoreVec[5] = maxScoreVec[5] > goldLigandVec[ clusterVec[i][j] ].get_ext_hb() ?
					maxScoreVec[5] : goldLigandVec[ clusterVec[i][j] ].get_ext_hb();
			scoreVec[5] += goldLigandVec[ clusterVec[i][j] ].get_ext_hb();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			minScoreVec[6] = minScoreVec[6] < goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					minScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			maxScoreVec[6] = maxScoreVec[6] > goldLigandVec[ clusterVec[i][j] ].get_ext_vdw() ?
					maxScoreVec[6] : goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();
			scoreVec[6] += goldLigandVec[ clusterVec[i][j] ].get_ext_vdw();

			off.width(15);
			off<< goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			minScoreVec[7] = minScoreVec[7] < goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					minScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			maxScoreVec[7] = maxScoreVec[7] > goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref() ?
					maxScoreVec[7] : goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();
			scoreVec[7] += goldLigandVec[ clusterVec[i][j] ].get_rmsd_ref();

			off.width(15);
			off<< overlayDisVec[ clusterVec[i][j] ] ;
			minScoreVec[8] = minScoreVec[8] < overlayDisVec[ clusterVec[i][j] ] ?
					minScoreVec[8] :  overlayDisVec[ clusterVec[i][j] ];
			maxScoreVec[8] = maxScoreVec[8] >  overlayDisVec[ clusterVec[i][j] ] ?
					maxScoreVec[8] :  overlayDisVec[ clusterVec[i][j] ] ;
			scoreVec[8] += overlayDisVec[ clusterVec[i][j] ] ;

			off<<endl;
		}
		off<<endl;
		off.width(5);
		off<<"min";
		for( size_t k=0; k<9; k++ ){
			off.width(15);
			off<<minScoreVec[k];
		}
		off.width(15);
		off<<getMinPair( distMatrix, clusterVec[i] );
		off<<endl;

		off.width(5);
		off<<"max";
		for( size_t k=0; k<9; k++ ){
			off.width(15);
			off<<maxScoreVec[k];
		}
		off.width(15);
		float	maxPair = getMaxPair( distMatrix, clusterVec[i] );
		off<<maxPair;
		off<<endl;

		off.width(5);
		off<<"aver";
		for( size_t k=0; k<9; k++ ){
			off.width(15);
			if( clusterVec[i].size() != 0 ){
				off<<scoreVec[k]/ clusterVec[i].size();
			}else{
				off<<scoreVec[k] ;
			}
		}
		off.width(15);
		off<<getAveragePair( distMatrix, clusterVec[i] );
		off<<endl<<endl;
	}
	off.close();
}

vector<GoldLigand>
read_goldRescoreFile( const string& fileName ){
	ifstream		inStream( fileName.c_str() );
	string			oneLine;
	float		fitness					=	0;
	float		inter_vdw				=0;
	float		inter_tors				=0;
	float		inter_hb				=0;
	float		inter_corre			=0;
	float		inter_energy		=0;
	float		ext_hb					=0;
	float		ext_vdw				=0;
	float		rmsd_ref				=0;
	vector<GoldLigand>				ligVec;
	while( getline( inStream, oneLine ) ){

		if( oneLine == "> <Gold.Goldscore.Fitness>" ){
			getline( inStream, oneLine );
			fitness = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.Vdw>" ){
			getline( inStream, oneLine );
			inter_vdw = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.Torsion>" ){
			getline( inStream, oneLine );
			 inter_tors= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.External.Vdw>" ){
			getline( inStream, oneLine );
			 ext_vdw = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Internal.HBond>" ){
			getline( inStream, oneLine );
			 inter_hb = atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.External.HBond>" ){
			getline( inStream, oneLine );
			 ext_hb= atof( oneLine.c_str() );
			 GoldLigand						newLig( fitness, inter_vdw, inter_tors, inter_hb, inter_corre,
					 inter_energy, ext_hb, ext_vdw, rmsd_ref );
			 ligVec.push_back( newLig );
				fitness					=	0;
				inter_vdw				=0;
				inter_tors				=0;
				inter_hb				=0;
				inter_corre			=0;
				inter_energy		=0;
				ext_hb					=0;
				ext_vdw				=0;
				rmsd_ref				=0;
		}		if( oneLine == "> <Gold.Goldscore.Internal.Correction>" ){
			getline( inStream, oneLine );
			 inter_corre= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Goldscore.Reference.RMSD>" ){
			getline( inStream, oneLine );
			 rmsd_ref= atof( oneLine.c_str() );
		}		if( oneLine == "> <Gold.Score>" ){
			getline( inStream, oneLine );
			getline( inStream, oneLine );
			vector<string>		tokenVec = tokenBySpace( oneLine );
			inter_energy = atof( tokenVec[4].c_str() );
		}
	}
	inStream.close();
	return				ligVec;
}
