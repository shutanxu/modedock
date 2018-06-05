/*
 * ss.cpp
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#include"ss.h"
#include"../algorithm/tokenBySpace.h"

//-------------------------------------------------------------------------------------------------------------------------

//---http://en.wikipedia.org/wiki/DSSP_%28protein%29
bool
dssp( const vector<Atom>& atVec, const Residue& res1, const Residue& res2 ){
	Atom							n1, c1, o1, h1 ;
	Atom							n2, c2, o2, h2 ;

	if( res1.get_atomVec().empty() || res2.get_atomVec().empty() ){
		return 					false;
	}

	for( size_t i=0; i<res1.get_atomVec().size(); i++ ){
		Atom						at = res1.get_atomVec()[i] ;

		if( at.get_name() == "N" ){
			n1 = at ;
		}else if( at.get_name() == "C" ){
			c1 = at ;
		}else if( at.get_name() == "O" ){
			o1 = at;
		}else if( at.get_name() == "H" ){
//			cout<<"?"<<endl;
			h1 = at;
		}
	}

	for( size_t i=0; i<res2.get_atomVec().size(); i++ ){
		Atom						at =  res2.get_atomVec()[i] ;
		if( at.get_name() == "N" ){
			n2 = at ;
		}else if( at.get_name() == "C" ){
			c2 = at ;
		}else if( at.get_name() == "O" ){
			o2 = at;
		}else if( at.get_name() == "H" ){
			cout<<"?"<<endl;
			h2 = at;
		}
	}

	float							dis_n1o2 = atomDis( n1, o2 );
	float							dis_o1n2 = atomDis( o1, n2 );
	float							dis_c1n2 = atomDis( c1, n2 ) ;
	float							dis_c2n1 = atomDis( c2, n1 ) ;
	float							energy = 0 ;

	if( dis_c1n2 < 2  || dis_c2n1<2  ){
		return 					false;
	}

	if( dis_n1o2 < dis_o1n2 ){
		float						r_on = atomDis( n1, o2 );
		float						r_ch = atomDis( h1, c2 );
		float						r_oh = atomDis( h1, o2 );
		float						r_cn = atomDis( n1, c2 );

		if( r_on == 0 || r_ch == 0 || r_oh== 0 || r_cn == 0 ){
			energy = 0;
		}else{
			energy = 0.084 * ( 1/r_on + 1/r_ch - 1/r_oh - 1/r_cn ) * 332 ;
		}
	}else {
		float						r_on = atomDis( n2, o1 );
		float						r_ch = atomDis( h2, c1 );
		float						r_oh = atomDis( h2, o1 );
		float						r_cn = atomDis( n2, c1 );

		if( r_on == 0 || r_ch == 0 || r_oh== 0 || r_cn == 0 ){
			energy = 0;
		}else{
			energy = 0.084 * ( 1/r_on + 1/r_ch - 1/r_oh - 1/r_cn ) * 332 ;
		}
	}

	if( energy < -0.5 ){
		//----hydrogen bond exist
		return 					true;
	}else{
		return						false;
	}
}

void
dssp( const vector<Atom>& atomVec, const vector<Residue>& residueVec ){

	vector< pair<size_t, size_t> >		helixVec;
	for( size_t i=0; i<residueVec.size(); i++ ){
		for( size_t j=i+ 1 ; j<residueVec.size(); j++ ){

			if(  dssp( atomVec, residueVec[i], residueVec[j] ) ){
//				cout.width(10);
//				cout<<left<<residueVec[i].get_index();
//				cout<<residueVec[j].get_index()<<endl;

				size_t				indexRange = 0;
				if( residueVec[i].get_index() >  residueVec[j].get_index() ){
					indexRange =  residueVec[i].get_index()  - residueVec[j].get_index() ;
				}else{
					indexRange =  residueVec[j].get_index()  - residueVec[i].get_index() ;
				}

//				indexRange = fabs(  residueVec[i].get_index()  - residueVec[j].get_index() );

				if( indexRange <= 4 ){

					helixVec.push_back( make_pair( residueVec[i].get_index(), residueVec[j].get_index() ) );
				}

			}
		}
	}

	cout<<"helix:"<<endl;
	for( size_t i=0; i<helixVec.size(); i++ ){
		cout<<helixVec[i].first<<" "<<helixVec[i].second<<endl;
	}
}
