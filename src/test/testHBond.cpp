/*
 * testHBond.cpp
 *
 *  Created on: Oct 4, 2014
 *      Author: stan
 */



#include"testHBond.h"

int
testHBond( int argc, char** argv ){
    return 0;
}

int
testHBond2( int argc, char** argv ){
//	QApplication application(argc,argv);
	vector<pair<float, float> >			hdVec;
	size_t								count = 0;

	vector<float>				xPlotVec  ;
	vector<vector<float> >		yPlotVec;
	vector<float>				yPlotVec1;
	vector<float>				yPlotVec2;
	ofstream					goldHbValueFile( "goldHbValueFile.dat" );
	ofstream					hydrogenAcceptorDisFile( "hydrogenAcceptorDisFile.dat" );
	for( size_t i= 1 ; i<2 ; i= i+1 ){

		stringstream		ss;
		ss<<i;
		string		ind = ss.str();

		string		fileNameB = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/dockResult/gold_soln_ligand_corina_m1_"
								+ ind + ".mol2" ;
		string		fileNameA = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/1ROB_protein.mol2"  ;

		string							goldParFile = "gold.params";
		Sybyl							goldPar( goldParFile );

		Molecular						molA( fileNameA );
		molA = readMolecularFile( fileNameA ).front();
		Molecular						molB( fileNameB );
		molB = readMolecularFile( fileNameB ).front();
        GoldLigand                  	lig( fileNameB );

		float							p = goldPar.computeHBondEnergy( molA, molB );

		vector<Atom> hDonorAtomVecA = goldPar.getHydrogenBondDonorAtoms( molA.get_atomVec() );
		vector<Atom> hDonorAtomVecB = goldPar.getHydrogenBondDonorAtoms( molB.get_atomVec() );
		vector<Atom> hAcceptorAtomVecA = goldPar.getHydrogenBondAcceptorAtoms( molA.get_atomVec() );
		vector<Atom> hAcceptorAtomVecB = goldPar.getHydrogenBondAcceptorAtoms( molB.get_atomVec() );
		for( size_t i=0; i<hDonorAtomVecA.size(); i++ ){
			for( size_t j=0; j<hAcceptorAtomVecB.size(); j++ ){
				float dis = atomDis( hDonorAtomVecA[i], hAcceptorAtomVecB[j] );
				if( dis < 3.2 ){
					hDonorAtomVecA[i].print();
					hAcceptorAtomVecB[j].print();
					cout<<endl;
				}
			}
		}
		for( size_t i=0; i<hAcceptorAtomVecA.size(); i++ ){
			for( size_t j=0; j<hDonorAtomVecB.size(); j++ ){
				float dis = atomDis( hAcceptorAtomVecA[i], hDonorAtomVecB[j] );
				if( dis < 3.2 ){
					hAcceptorAtomVecA[i].print();
					hDonorAtomVecB[j].print();
					cout<<endl;
				}
			}
		}

		cout<<i<<" p:"<<p<<endl;
	}

//	return application.exec();
	return 0;
}


int
testHBond4( int argc, char** argv ){
//	QApplication application(argc,argv);
	vector<pair<float, float> >			hdVec;
	size_t								count = 0;

	vector<float>				xPlotVec  ;
	vector<vector<float> >		yPlotVec;
	vector<float>				yPlotVec1;
	vector<float>				yPlotVec2;
	ofstream					goldHbValueFile( "goldHbValueFile.dat" );
	ofstream					hydrogenAcceptorDisFile( "hydrogenAcceptorDisFile.dat" );
	for( size_t i= 1 ; i<2 ; i= i+1 ){

		stringstream		ss;
		ss<<i;
		string		ind = ss.str();

		string		fileNameB = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/dockResult/gold_soln_ligand_corina_m1_"
								+ ind + ".mol2" ;
		string		fileNameA = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/1ROB_protein.mol2"  ;

		string							goldParFile = "gold.params";
		Sybyl							goldPar( goldParFile );

		Molecular						molA( fileNameA );
		molA = readMolecularFile( fileNameA ).front();
		Molecular						molB( fileNameB );
		molB = readMolecularFile( fileNameB ).front();
		GoldLigand					lig( fileNameB );

		float							p = goldPar.computeHBondEnergy( molA, molB );

		vector<Atom> hDonorAtomVecA = goldPar.getHydrogenBondDonorAtoms( molA.get_atomVec() );
		vector<Atom> hDonorAtomVecB = goldPar.getHydrogenBondDonorAtoms( molB.get_atomVec() );
		vector<Atom> hAcceptorAtomVecA = goldPar.getHydrogenBondAcceptorAtoms( molA.get_atomVec() );
		vector<Atom> hAcceptorAtomVecB = goldPar.getHydrogenBondAcceptorAtoms( molB.get_atomVec() );
		for( size_t i=0; i<hDonorAtomVecA.size(); i++ ){
			for( size_t j=0; j<hAcceptorAtomVecB.size(); j++ ){
				float dis = atomDis( hDonorAtomVecA[i], hAcceptorAtomVecB[j] );
				if( dis < 3.2 ){
					hDonorAtomVecA[i].print();
					hAcceptorAtomVecB[j].print();
					cout<<endl;
				}
			}
		}
		for( size_t i=0; i<hAcceptorAtomVecA.size(); i++ ){
			for( size_t j=0; j<hDonorAtomVecB.size(); j++ ){
				float dis = atomDis( hAcceptorAtomVecA[i], hDonorAtomVecB[j] );
				if( dis < 3.2 ){
					hAcceptorAtomVecA[i].print();
					hDonorAtomVecB[j].print();
					cout<<endl;
				}
			}
		}
		cout<<i<<" p:"<<p<<endl;
	}

//	return application.exec();
	return 0;
}
