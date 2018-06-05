/*
 * graph.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: stan
 */

#include<algorithm>
#include<list>

#include"graph.h"

using namespace std;

bool
visited( const vector<pair<size_t, bool> >& 	visitedVec, const size_t& index ){
	for( size_t i=0; i<visitedVec.size(); i++ ){
		if( visitedVec[i].first == index ){
			return 		visitedVec[i].second;
		}
	}
	cerr<<"visited error, can not find atom "<<index<<endl;
	return					false;
}

void
setVisited( vector<pair<size_t, bool> >& 	visitedVec, const size_t& index ){
	for( size_t i=0; i<visitedVec.size(); i++ ){
		if( visitedVec[i].first == index ){
			visitedVec[i].second = true;
			return;
		}
	}
	cerr<<" set visited error, can not find atom "<<index<<endl;
}

bool
allVisited(  const vector<pair<size_t, bool> >& 	visitedVec ){
	for( size_t i=0; i<visitedVec.size(); i++ ){
		if( visitedVec[i].second == false ){
			return false;
		}
	}
	return true;
}

void
dfsMol( const vector<Bond>&  				bdVec,
		const size_t&										index,
		vector<pair<size_t, bool> >&         	visitedVec,
		vector<size_t>&									historyVec )	{

	if( !visited( visitedVec, index ) ){
		setVisited( visitedVec, index );
		if( find( historyVec.begin(), historyVec.end(), index ) == historyVec.end() ){
			historyVec.push_back( index );
		}
	}
	vector<size_t>										neighborOfIndex;
	for( size_t i=0; i<bdVec.size(); i++ ){
		if( bdVec[i].get_firstAtomIndex() == index ){
			neighborOfIndex.push_back( bdVec[i].get_secAtomIndex() );
		}else if( bdVec[i].get_secAtomIndex() == index ){
			neighborOfIndex.push_back( bdVec[i].get_firstAtomIndex() );
		}
	}

	for( size_t i=0; i<neighborOfIndex.size(); i++ ){
		if( !visited( visitedVec, neighborOfIndex[i] ) ){
			dfsMol( bdVec, neighborOfIndex[i], visitedVec, historyVec );
		}
	}
}

/**
 * deep first search from direction: directionStartIndex -> directionNextIndex
 * @param bdVec
 * @param index
 * @param fatherIndex
 * @return
 */
vector<Atom>
dfsMol(	const vector<Atom>&          atomVec,
		const vector<Bond>&          bdVec,
		const	size_t&              directionStartIndex,
		const	size_t&              directionNextIndex ){

	vector< pair<size_t, bool> >				visitedVec;
	for( size_t i=0; i<atomVec.size(); i++ ){
		pair<size_t, bool>					visite;
		if( atomVec[i].get_index() == directionStartIndex ){
			pair<size_t, bool>				p( atomVec[i].get_index(), true );
			visite = p;
		}else{
			pair<size_t, bool>				p( atomVec[i].get_index(), false );
			visite = p;
		}
		visitedVec.push_back( visite );
	}
	vector<size_t>										historyVec;
	historyVec.push_back( directionStartIndex );

	vector< vector<size_t> >						allRingVec;

	//	dfsMol( bdVec, directionNextIndex, directionStartIndex, visitedVec, historyVec, allRingVec );
	dfsMol( bdVec, directionNextIndex, visitedVec, historyVec );

	vector<Atom>										searchedAtomVec;
	for( size_t i=0; i<historyVec.size(); i++ ){
		Atom			at = searchAtom( atomVec,  historyVec[i] );
		searchedAtomVec.push_back( at );
	}

	return			searchedAtomVec;
}

void
dfsMol( const vector<Bond>&             bdVec,
		const size_t&                   index,
		const size_t&                   fatherIndex,
		vector<pair<size_t, bool> >&    visitedVec,
		vector<size_t>&                 historyVec,
		vector<vector<size_t> >&        ringVec ){

	vector<size_t>						neighborVec;
	for( size_t i=0; i<bdVec.size(); i++ ){
		if( bdVec[i].get_firstAtomIndex( ) == index && bdVec[i].get_secAtomIndex() != fatherIndex  ){
			neighborVec.push_back( bdVec[i].get_secAtomIndex() );
		}else if( bdVec[i].get_secAtomIndex() == index && bdVec[i].get_firstAtomIndex() != fatherIndex ){
			neighborVec.push_back( bdVec[i].get_firstAtomIndex() );
		}
	}

	setVisited( visitedVec, index );

	for( size_t i=0; i<neighborVec.size(); i++ ){
		if( !visited( visitedVec, neighborVec[i] ) ){

			historyVec.erase( find( historyVec.begin(), historyVec.end(), index )  , historyVec.end() );
			//			if( find( historyVec.begin(), historyVec.end(), index ) != historyVec.end() ){
			//				historyVec.erase( find( historyVec.begin(), historyVec.end(), index )  , historyVec.end() );
			//			}
			historyVec.push_back( index );
			if(0){
				cout<<"history: ";
				for( size_t i=0; i<historyVec.size(); i++ ){
					cout<<historyVec[i]<<" ";
				}
				cout<<endl;
			}
			dfsMol( bdVec, neighborVec[i], index, visitedVec,  historyVec, ringVec  );
		}else if ( visited( visitedVec, neighborVec[i] )  &&
				find( historyVec.begin(), historyVec.end(), index ) == historyVec.end()   ) {
			historyVec.push_back( index );

			if(0){
				cout<<"ring: " ;
				vector<size_t>::iterator		it =  find( historyVec.begin(), historyVec.end(), neighborVec[i] ) ;
				for( it =  find( historyVec.begin(), historyVec.end(), neighborVec[i] ); it<historyVec.end(); it++ ){
					cout<<*it<<" ";
				}
				cout<<" "<<neighborVec[i];
				cout<<endl;
			}
			vector<size_t>					ring;
			vector<size_t>::iterator		it =  find( historyVec.begin(), historyVec.end(), neighborVec[i] ) ;
			for( it =  find( historyVec.begin(), historyVec.end(), neighborVec[i] ); it<historyVec.end(); it++ ){
				ring.push_back( *it );
			}
			ringVec.push_back( ring );
		}
	}
}

vector<Atom>
bfsMol( const vector<Atom>&             atomVec,
		const vector<Bond>&             bdVec,
		const size_t&                   startAtomIndex ){

	vector< pair<size_t, bool> >				visitedVec ;
	for( size_t i=0; i<atomVec.size(); i++ ){
		pair<size_t, bool>  p = make_pair< size_t, bool >( atomVec[i].get_index(), false );
		visitedVec.push_back( p );
	}

	Atom						startAtom = searchAtom( atomVec, startAtomIndex );
	list<Atom>					atomSeq;
	vector<Atom>				bfsAtomVec;
	atomSeq.push_back( startAtom );

	while( ! atomSeq.empty() ){

		Atom  at = atomSeq.front();
		bool  visited = false;
		for( size_t i=0; i<visitedVec.size(); i++ ){
			if( visitedVec[i].first  == at.get_index() ){
				if( visitedVec[i].second == true ){
//					cout<<"!!! atom visited ??? ";
//					at.print();
					visited = true;
				}
				visitedVec[i].second = true;
			}
		}

		if( visited == true ){
			atomSeq.pop_front();
		}else{
			atomSeq.pop_front();
			bfsAtomVec.push_back( at );
			for( size_t i=0; i<bdVec.size(); i++ ){
				if( bdVec[i].get_firstAtomIndex() == at.get_index() ){
					bool 		visitedFlag = false ;
					for( size_t j=0; j<visitedVec.size(); j++ ){
						if( visitedVec[j].first == bdVec[i].get_secAtomIndex() ){
							visitedFlag = visitedVec[j].second;
						}
					}
					if( !visitedFlag ){
						bdVec[i].get_secAtom().print();
						atomSeq.push_back( bdVec[i].get_secAtom() );
					}
				}
				if( bdVec[i].get_secAtomIndex() == at.get_index() ){
					bool 		visitedFlag = false ;
					for( size_t j=0; j<visitedVec.size(); j++ ){
						if( visitedVec[j].first == bdVec[i].get_firstAtomIndex() ){
							visitedFlag = visitedVec[j].second;
						}
					}
					if( !visitedFlag ){
						bdVec[i].get_firstAtom().print();
						atomSeq.push_back( bdVec[i].get_firstAtom() );
					}
				}
			}
		}
	}

	if( bfsAtomVec.size() != atomVec.size() ){
		cerr<<"!!! bfs error"<<endl;
	}

	return  bfsAtomVec;
}

/**
 * broad first search, the input bdVec must be connected
 */
vector<Bond>
getBfsBonds( const vector<Bond>& bdVec, const size_t& startAtomIndex ){

    vector<Atom> atomVec;
	for( size_t i=0; i<bdVec.size(); i++ ){
		Atom a1 = bdVec[i].get_firstAtom();
		Atom a2 = bdVec[i].get_secAtom();
		if( find( atomVec.begin(), atomVec.end(), a1 ) == atomVec.end() ){
			atomVec.push_back(a1);
		}
		if( find( atomVec.begin(), atomVec.end(), a2 ) == atomVec.end() ){
			atomVec.push_back(a2);
		}
	}

	vector<Bond> bfsBondVec;

	vector< pair<size_t, bool> >				visitedVec ;
	for( size_t i=0; i<atomVec.size(); i++ ){
		pair<size_t, bool>  p = make_pair< size_t, bool >( atomVec[i].get_index(), false );
		visitedVec.push_back( p );
	}

	Atom						startAtom = searchAtom( atomVec, startAtomIndex );
	list<Atom>					atomSeq;
	vector<Atom>				bfsAtomVec;
	atomSeq.push_back( startAtom );

	while( ! atomSeq.empty() ){

		Atom  at = atomSeq.front();
		bool  visited = false;
		for( size_t i=0; i<visitedVec.size(); i++ ){
			if( visitedVec[i].first  == at.get_index() ){
				if( visitedVec[i].second == true ){
//					cout<<"!!! atom visited ??? ";
//					at.print();
					visited = true;
				}
				visitedVec[i].second = true;
			}
		}

		if( visited == true ){
			atomSeq.pop_front();
		}else{
			atomSeq.pop_front();
			bfsAtomVec.push_back( at );
			for( size_t i=0; i<bdVec.size(); i++ ){
				if( bdVec[i].get_firstAtomIndex() == at.get_index() ){
					bool 		visitedFlag = false ;
					for( size_t j=0; j<visitedVec.size(); j++ ){
						if( visitedVec[j].first == bdVec[i].get_secAtomIndex() ){
							visitedFlag = visitedVec[j].second;
						}
					}
					if( !visitedFlag ){
//						bdVec[i].get_secAtom().print();
						atomSeq.push_back( bdVec[i].get_secAtom() );
						if( find( bfsBondVec.begin(), bfsBondVec.end(), bdVec[i] ) == bfsBondVec.end() ){
							bfsBondVec.push_back( bdVec[i] );
						}
					}
				}
				if( bdVec[i].get_secAtomIndex() == at.get_index() ){
					bool 		visitedFlag = false ;
					for( size_t j=0; j<visitedVec.size(); j++ ){
						if( visitedVec[j].first == bdVec[i].get_firstAtomIndex() ){
							visitedFlag = visitedVec[j].second;
						}
					}
					if( !visitedFlag ){
//						bdVec[i].get_firstAtom().print();
						atomSeq.push_back( bdVec[i].get_firstAtom() );
						if( find( bfsBondVec.begin(), bfsBondVec.end(), bdVec[i] ) == bfsBondVec.end() ){
							bfsBondVec.push_back( bdVec[i] );
						}
					}
				}
			}
		}
	}

	if( bfsAtomVec.size() != atomVec.size() ){
		cerr<<"!!! bfs error"<<endl;
	}

	return  bfsBondVec;
}

vector<vector<size_t> >
connectMatrix(const vector<Atom> atomVec, const vector<Bond> bdVec){

    size_t atNum = atomVec.size();
    vector<vector<size_t> > matrix( atNum, vector<size_t>(atNum, 0) );

    for( size_t atomInd=0; atomInd<atNum; atomInd++ ){
        vector< pair<size_t, bool> >				visitedVec ;
        for( size_t i=0; i<atomVec.size(); i++ ){
            pair<size_t, bool>  p = make_pair< size_t, bool >( atomVec[i].get_index(), false );
            visitedVec.push_back( p );
        }

        Node						startNode( atomVec[atomInd], 0);
        list<Node>					nodeSeq;
        vector<Atom>				bfsAtomVec;
        nodeSeq.push_back( startNode );

        while( ! nodeSeq.empty() ){

            Atom  at = nodeSeq.front().at;
            bool  visited = false;
            for( size_t i=0; i<visitedVec.size(); i++ ){
                if( visitedVec[i].first  == at.get_index() ){
                    if( visitedVec[i].second == true ){
                        //					cout<<"!!! atom visited ??? ";
                        //					at.print();
                        visited = true;
                    }
                    visitedVec[i].second = true;
                }
            }

            size_t ind = find( atomVec.begin(), atomVec.end(), at ) - atomVec.begin();

            int dep = nodeSeq.front().depth;
            matrix[atomInd][ind] = dep;

            if( visited == true ){
                nodeSeq.pop_front();
            }else{
                nodeSeq.pop_front();
                bfsAtomVec.push_back( at );
                for( size_t i=0; i<bdVec.size(); i++ ){
                    if( bdVec[i].get_firstAtomIndex() == at.get_index() ){
                        bool 		visitedFlag = false ;
                        for( size_t j=0; j<visitedVec.size(); j++ ){
                            if( visitedVec[j].first == bdVec[i].get_secAtomIndex() ){
                                visitedFlag = visitedVec[j].second;
                            }
                        }
                        if( !visitedFlag ){
//                            bdVec[i].get_secAtom().print();
                            Node n( bdVec[i].get_secAtom(), dep+1 );
                            nodeSeq.push_back( n );
                        }
                    }
                    if( bdVec[i].get_secAtomIndex() == at.get_index() ){
                        bool 		visitedFlag = false ;
                        for( size_t j=0; j<visitedVec.size(); j++ ){
                            if( visitedVec[j].first == bdVec[i].get_firstAtomIndex() ){
                                visitedFlag = visitedVec[j].second;
                            }
                        }
                        if( !visitedFlag ){
//                            bdVec[i].get_firstAtom().print();
                            Node n( bdVec[i].get_firstAtom(), dep+1 );
                            nodeSeq.push_back( n );
                        }
                    }
                }
            }
        }
    }

    if(0){
        cout.width(6);
        cout<<1;
        for( size_t i=1; i<matrix.size(); i++ ){
            cout.width(3);
            cout<<i+1;
        }
        cout<<endl;
        for( size_t i=0; i<matrix.size(); i++ ){
            cout.width(5);
            cout<<left<<i+1;
            for( size_t j=0; j<matrix[i].size(); j++ ){
                cout.width(3);
                cout<<matrix[i][j];
            }
            cout<<endl;
        }
    }

    return matrix;
}
