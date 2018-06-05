/*
 * testGold.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: stan
 */




#include"testGold.h"

using namespace std;

int
testGoldVDW( int argc, char** argv ){
	QApplication application(argc,argv);
	vector<pair<float, float> >			vdwVec;
	size_t								count = 0;

	for( size_t i=1; i<=500 ; i++ ){
		stringstream		ss;
		ss<<i;
		string		ind = ss.str();

		string		fileNameB = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/dockResult/gold_soln_ligand_corina_m1_" + ind + ".mol2" ;
		string		fileNameA = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/1ROB_protein.mol2"  ;

		string							goldParFile = "gold.params";
		Sybyl							goldPar( goldParFile );

		Molecular						molA( fileNameA );
		molA = readMolecularFile( fileNameA ).front();
		Molecular						molB( fileNameB );
		molB = readMolecularFile( fileNameB ).front();
		GoldLigand					lig( fileNameB );

		float							p = goldPar.computeGoldVDW( molA, molB );

		if(p<100){
			vdwVec.push_back( make_pair( p, lig.get_ext_vdw() ) );
			cout.width(5);
			cout<<left<<i;
			cout<<"poteintial:";
			cout.width(15);
			cout<<left<<p;
			cout<<"gold value:"<<lig.get_ext_vdw()<<endl;
			count++;
		}
	}
	cout<<"count:"<<count<<endl;

//	plotData		pd( vdwVec );
//	pd.show();

	return application.exec();
}
