/*
 * testPDBBind.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: stan
 */

#include<utility>
#include"testPDBBind.h"

using namespace std;

int
testPDBbind(   int argc, char** argv   ){
	QApplication application(argc,argv);
	ifstream			iif1( "/mnt/2t/xushutan/CASF-2013/data-scoring/dSAS.TXT" );
	ifstream			iif2( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-Jain.TXT" );
	ifstream			iif3( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-LigScore1.TXT" );
	ifstream			iif4( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-LigScore2.TXT" );
	ifstream			iif5( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-LUDI1.TXT" );
	ifstream			iif6( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-LUDI2.TXT" );
	ifstream			iif7( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-LUDI3.TXT" );
	ifstream			iif8( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-PLP1.TXT" );
	ifstream			iif9( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-PLP2.TXT" );
	ifstream			iif10( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-PMF04.TXT" );
	ifstream			iif11( "/mnt/2t/xushutan/CASF-2013/data-scoring/DS-PMF.TXT" );
	ifstream			iif12( "/mnt/2t/xushutan/CASF-2013/data-scoring/GlideScore-SP.TXT" );
	ifstream			iif13( "/mnt/2t/xushutan/CASF-2013/data-scoring/GlideScore-XP.TXT" );
	ifstream			iif14( "/mnt/2t/xushutan/CASF-2013/data-scoring/GOLD-ASP.TXT" );
	ifstream			iif15( "/mnt/2t/xushutan/CASF-2013/data-scoring/GOLD-ChemPLP.TXT" );
	ifstream			iif16( "/mnt/2t/xushutan/CASF-2013/data-scoring/GOLD-ChemScore.TXT" );
	ifstream			iif17( "/mnt/2t/xushutan/CASF-2013/data-scoring/GOLD-GoldScore.TXT" );
	ifstream			iif18( "/mnt/2t/xushutan/CASF-2013/data-scoring/MOE-Affinity_dG.TXT" );
	ifstream			iif19( "/mnt/2t/xushutan/CASF-2013/data-scoring/MOE-Alpha_HB.TXT" );
	ifstream			iif20( "/mnt/2t/xushutan/CASF-2013/data-scoring/MOE-ASE.TXT" );
	ifstream			iif21( "/mnt/2t/xushutan/CASF-2013/data-scoring/MOE-London_dG.TXT" );
	ifstream			iif22( "/mnt/2t/xushutan/CASF-2013/data-scoring/SYBYL-ChemScore.TXT" );
	ifstream			iif23( "/mnt/2t/xushutan/CASF-2013/data-scoring/SYBYL-DScore.TXT" );
	ifstream			iif24( "/mnt/2t/xushutan/CASF-2013/data-scoring/SYBYL-GScore.TXT" );
	ifstream			iif25( "/mnt/2t/xushutan/CASF-2013/data-scoring/SYBYL-PMF.TXT" );
	ifstream			iif26( "/mnt/2t/xushutan/CASF-2013/data-scoring/XScore-Average.TXT" );
	ifstream			iif27( "/mnt/2t/xushutan/CASF-2013/data-scoring/XScore-HMScore.TXT" );
	ifstream			iif28( "/mnt/2t/xushutan/CASF-2013/data-scoring/XScore-HPScore.TXT" );
	ifstream			iif29( "/mnt/2t/xushutan/CASF-2013/data-scoring/XScore-HSScore.TXT" );

	string				sline;
	vector<pair<float, float> >			energyVec;
	getline( iif16, sline );
	while( getline( iif13, sline ) ){
		vector<string>			strVec = tokenBySpace( sline );
		float							kd = atof( strVec[1].c_str() );
		float							score_cry = atof( strVec[2].c_str() );
		if( kd<1000 && score_cry < 1000 ){
//			pair<float, float>		p = std::make_pair<float, float>( kd, score_cry );
//			energyVec.push_back( p );
		}
	}
//	plotData							pd( energyVec );
//	pd.show();
	return application.exec();
}
