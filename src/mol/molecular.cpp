/*
 * molecular.cpp
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<algorithm>
#include<math.h>
#include<utility>
#include<functional>
#include<string>
#include<time.h>

#include"molecular.h"
#include"../algorithm/graph.h"
#include"../algorithm/tokenBySpace.h"
#include"../algorithm/overlay.h"
#include"../basicMath/geomBasic.h"

using namespace std;

const float maxBondLength = 1.70;

Molecular::Molecular(){
	isEmpty = false;
	name = string("");
	fileName = string("");
	center=Coord(0,0,0);
}

Molecular::Molecular(const string& file){
	fileHeader.clear();
	chainVec.clear();
	//	readFile( file );
}

Molecular::Molecular( const Molecular& mol ){
	isEmpty = mol.empty();
	name = mol.get_name();
	fileName = mol.get_fileName();
	center = mol.get_center();

	//	fileHeader = vector<string>( mol.get_fileHeader().size() );
	fileHeader = mol.get_fileHeader();

	//	chainVec = vector<Chain>( mol.get_chainVec().size() );
	chainVec = mol.get_chainVec();

	//	bondVec = vector<Bond>( mol.get_bondVec().size() );
	bondVec = mol.get_bondVec();

	ss = mol.get_secondaryStructure();
	//	cout<<"="<<endl;
}

vector<Molecular>
readMolecularFile( const string& file ){
	//	isEmpty = true;
	if( file.empty() ){
		Molecular	mol;
		vector<Molecular>		molVec;
		molVec.clear();
		return		molVec;
	}

	vector<Molecular>			molVec;
	molVec.clear();
	string						mol2( "mol2" );
	string						mol( "mol" );
	string						pdb( "pdb" );
	string						pdbqt( "pdbqt");

	size_t			typeIndex = file.find( '.' );
	string			type = file.substr( typeIndex+1 );
	string			fileName( file );
	for (std::string::reverse_iterator rit=fileName.rbegin(); rit!=fileName.rend(); ++rit){
		if( *rit == '.' ){
			size_t 		ind = rit - fileName.rbegin();
			type = fileName.substr(   fileName.length() - ind   );
			break;
		}
	}

	if( type == mol2 ){
		molVec = readMol2( file);
	}else if( type == mol ){
		molVec = readMol( file );
	}else if( type == pdb ){
		molVec = readPdb( file );
	}else if( type == "pdbqt" ){
		molVec = read_pdbqt( file );
	}else{
		cerr<<"--- warnning!!! can not read file "<<file<<"---"<<endl;
		return molVec ;
	}
	return  molVec ;
}

vector<Molecular>
readMol2( const string& file ){
    time_t t_start=clock();
	ifstream iff;
	iff.open( file.c_str() );
	if( !iff ){
		cerr<<"---error can not open "<<file<<"---"<<endl;
		throw;
	}
	string					oneLine;
	vector<string>			strVec;
	bool 					startAtom = false;
	bool 					startBond = false;
	bool					readHead = true;
	bool					recordName = false;
	vector<Atom>			residueAtomsVec;
	vector<Atom>			atomsOneMolVec;
	vector<Residue> 		residueVec;
	vector<Molecular>		molVec;
	size_t					preResidIndex 	= 10e10;
	size_t					newResidIndex = 10e10;
	Chain					chain;
	Residue*				preResidue;
	Atom					preAtom;
	Atom					tempAt;
	bool					isEmpty;
	string					name;
	string					fileName ;
	Coord					center;
	vector<string >			fileHeader;
	vector<Chain>			chainVec;
	vector<Bond>			bondVec;
	SStructure				ss;						// secondary structure

	while( getline( iff, oneLine ) ){

		strVec = tokenBySpace( oneLine );

		if( recordName && strVec.size() == 1 ){
			name = strVec[0];
			recordName = false;
		}
		if( strVec[0] == "@<TRIPOS>MOLECULE" ){
			recordName = true;
			startAtom = false;
			startBond = false;

			if( ! chain.empty() ){

				if( name == "" ){
					name = string( "no name" );
				}
				chainVec.push_back( chain );
				Molecular		mol( isEmpty, name, fileName, center,
						fileHeader, chainVec, bondVec, ss );
				if( ! mol.empty() ){
					cout<<"molecular: "<<name<<endl;
					molVec.push_back( mol );
				}else{
					//				cout<<"empty mol"<<endl;
				}

				chain.clear();
				chainVec.clear();
				name.clear();
				fileName.clear();
				fileHeader.clear();
				bondVec.clear();
				startAtom = false;
				startBond = false;
			}
		}
		if( strVec[0] == "@<TRIPOS>ATOM" ){
			readHead = false;
			startAtom = true;
			startBond = false;
		}else if( strVec[0] == "@<TRIPOS>BOND" ){
			startAtom = false;
			startBond = true;
        }else if( strVec[0] == "@<TRIPOS>SUBSTRUCTURE" ){

			if( name == "" ){
				name = string( "no name" );
			}
			chainVec.push_back( chain );
//			cout<<"bond:"<<bondVec.size()<<endl;
			Molecular		mol( isEmpty, name, fileName, center,
					fileHeader, chainVec, bondVec, ss );
			if( ! mol.empty() ){
				cout<<"molecular: "<<name<<endl;
				molVec.push_back( mol );
			}else{
				//				cout<<"empty mol"<<endl;
			}

			chain.clear();
			chainVec.clear();
			name.clear();
			fileName.clear();
			fileHeader.clear();
			bondVec.clear();
			startAtom = false;
			startBond = false;
		}
		if( startAtom && strVec.size() > 8 ){
			isEmpty = false;
            Atom at( strVec );
			atomsOneMolVec.push_back( at );
			tempAt = at;
			newResidIndex = at.get_residueIndex();
//			cout<<"res Ind:"<<newResidIndex<<endl;
//			at.print();
			if( newResidIndex == preResidIndex ) {
				residueAtomsVec.push_back( at );
			}else if( newResidIndex == preResidIndex + 1 &&  !residueAtomsVec.empty() ){
				Residue				resid;
				Residue				tempResid( preAtom.get_residueName(),
						preAtom.get_residueIndex(), residueAtomsVec );
				resid = tempResid;

				chain.add_Residue( resid );
				residueAtomsVec.clear();
				residueAtomsVec.push_back( at );
			}else if( newResidIndex != preResidIndex + 1 &&
					newResidIndex !=  preResidIndex && !residueAtomsVec.empty() ){
				Residue				resid;
				Residue				tempResid( preAtom.get_residueName(),
						preAtom.get_residueIndex(), residueAtomsVec );
				resid = tempResid;

				chain.add_Residue( resid );
				residueAtomsVec.clear();
				residueAtomsVec.push_back( at );
			}else{
				residueAtomsVec.push_back( at );
			}
			preResidIndex = newResidIndex;
			preAtom	= at ;
		}else if( startBond && strVec.size() > 3 ){
//            cout<<oneLine<<endl;
			if( ! residueAtomsVec.empty()  ){
				Residue				resid;
				Residue				tempResid( preAtom.get_residueName(),
						preAtom.get_residueIndex(), residueAtomsVec );
				resid = tempResid;
				chain.add_Residue( resid );
				residueAtomsVec.clear();
			}
			size_t	index = atoi( strVec[0].c_str() );
			size_t 	firstIndex = atoi( strVec[1].c_str() );
			size_t 	secondIndex = atoi( strVec[2].c_str() );

			Atom	firstAt =  searchAtom( atomsOneMolVec, firstIndex );
			Atom	secondAt =  searchAtom( atomsOneMolVec, secondIndex );
			string	bondType = strVec[3];

			Bond	bond( index,  firstAt, secondAt, bondType );
			bondVec.push_back( bond );
//            bond.print();
		}
		if( !startAtom && readHead ){
			fileHeader.push_back( oneLine );
		}
	}
	iff.close();

    if( ! chain.empty()  ){
		if( name == "" ){
			name = string( "no name" );
		}
		chainVec.push_back( chain );
		if(1){ //update bonds
			vector<Atom>					allAtomVec;
			for( size_t i=0; i<chainVec.size(); i++ ){
				vector<Atom>			chainAtVec = chainVec[i].get_AtomVec();
				for( size_t j=0; j<chainAtVec.size(); j++ ){
					allAtomVec.push_back( chainAtVec[j] );
				}
			}

            if( bondVec.empty() ){
                for( size_t i=0; i<allAtomVec.size(); i++  ){
                    for( size_t j=i+1; ( j<allAtomVec.size()  ); j++ ){
                        Atom a1 = allAtomVec[i];
                        Atom a2 = allAtomVec[j];
                        float	dist12 =  atomDis( a1, a2 ) ;
                        if( dist12 < maxBondLength && !(a1.get_type() == "H" && a2.get_type() == "H") ){
                            size_t 	firstIndex = a1.get_index();
                            size_t 	secondIndex = a2.get_index();
                            string	bondType = "!!";
                            Bond	bond( 0, a1, a2, bondType );
                            bondVec.push_back( bond );
                        }
                    }
                }
            }
		}

		Molecular		mol( isEmpty, name, fileName,
				center, fileHeader, chainVec, bondVec, ss );

		if( ! mol.empty() ){
			cout<<"molecular: "<<name<<endl;
			molVec.push_back( mol );
		}else{
			cout<<"empty mol"<<endl;
		}
		chain.clear();
		chainVec.clear();
		name.clear();
		fileName.clear();
		fileHeader.clear();
		bondVec.clear();
		startAtom = false;
		startBond = false;
	}
	if( molVec.size() > 4 ){
		cout<<"mol num:"<<molVec.size()<<endl;
	}
    cout<<"time:"<<float( clock() - t_start )/ CLOCKS_PER_SEC<<endl;
	return molVec;
}

vector<Molecular>
readMol( const string& file ){}

vector<Molecular>
readPdb( const string& file ){
    time_t t_start=clock();
	ifstream iff( file.c_str() );
	if( !iff ){
		cerr<<"---error can not open "<<file<<"---"<<endl;
		throw;
	}
	vector<Molecular>   molVec;
	string              oneLine;
	vector<string>      strVec;
	bool                startAtom = false;
	bool                startBond = false;
	bool                readHead = true;
	bool                recordName = false;
	vector<Atom>        residueAtomsVec;
	vector<Residue>     residueVec;
	size_t              preResidIndex 	= 10e10;
	size_t              newResidIndex = 10e10;
	Chain               chain;
	Residue*            preResidue;
	Atom                preAtom;
	Atom                tempAt;

	vector< pair< size_t, size_t > >        helixVec;			// pair< helix start residue index, helix end residue index >
	vector< pair< size_t, size_t > >        sheetVec;		// pair< sheet start residue index, sheet end residue index >

	bool                isEmpty;
	string              name;
	string              fileName  = file ;
	Coord               center;
	vector<string >     fileHeader;
	vector<Chain>       chainVec;
	vector<Bond>        bondVec;
	vector<Atom>        atomsOneMolVec;
	SStructure          ss;						// secondary structure

	Chain               residueChain;
	residueChain.clear();

	clock_t             t1 = clock();
	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		string			head = tokenBySpace( oneLine.substr(0, 6) )[0];
		//		cout<<"head:"<<head<<endl;

		if( strVec[0] == "MODEL" ){
			name = strVec.back();
		}
		if( strVec[0] == "HEADER" ){
			name = strVec.back();
		}
		if( strVec[0] == "ATOM" || strVec[0] == "HETATM" ){
			startAtom = true;
		}

		if( strVec[0] == "HELIX" ){
			size_t		startResid = atoi( strVec[5].c_str() ) ;
			size_t		endResid = atoi( strVec[8].c_str() ) ;
			pair<size_t, size_t>	helix = make_pair( startResid, endResid );
			helixVec.push_back( helix );
		}else if( strVec[0] == "SHEET" ){
			size_t		startResid = atoi( strVec[6].c_str() ) ;
			size_t		endResid = atoi( strVec[9].c_str() ) ;
			pair<size_t, size_t>	sheet = make_pair( startResid, endResid );
			sheetVec.push_back( sheet );
		}else	if( ( head == "ATOM" || head == "HETATM" ) && oneLine.length() >= 70 ){
			isEmpty = false;
			Atom at( oneLine );
			tempAt = at;
			atomsOneMolVec.push_back( at );
			//			at.print();

			if( at.get_residueIndex() != preAtom.get_residueIndex() && !residueAtomsVec.empty() ){
				Residue				resid( preAtom.get_residueName(),
						preAtom.get_residueIndex(), residueAtomsVec );

				if( chainVec.empty() ){
					Chain				ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain				ch( resid );
						chainVec.push_back( ch );
					}
				}
				residueAtomsVec.clear();
				residueAtomsVec.push_back( at );
			}else{
				residueAtomsVec.push_back( at );
			}

			preResidIndex = newResidIndex;
			preAtom	= at ;
		}else if( strVec[0] == "CONECT" ){
			for( size_t i=2; i<strVec.size()-2; i++ ){
				size_t 	firstIndex = atoi( strVec[1].c_str() );
				size_t 	secondIndex = atoi( strVec[i].c_str() );
				Atom	a1 = searchAtom( atomsOneMolVec, firstIndex );
				Atom	a2 = searchAtom( atomsOneMolVec, secondIndex );
				string	bondType = "";
				Bond	bond( 0, a1, a2, bondType );
				if( find( bondVec.begin(), bondVec.end(), bond ) == bondVec.end() ){
					bondVec.push_back( bond );
				}
			}
		}else if( strVec[0] == "TER" ){
			if( ! residueAtomsVec.empty() ){
				// save the last atom
				Residue				resid( tempAt.get_residueName(),
						                   tempAt.get_residueIndex(),
						                   residueAtomsVec );
				if( chainVec.empty() ){
					Chain ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain ch( resid );
						chainVec.push_back( ch );
					}
				}
				chain.clear();
				residueAtomsVec.clear();
				//				residueAtomsVec.push_back( at );
			}
		}else if (strVec[0]=="ENDMDL" || strVec[0] == "END"  ){
			atomsOneMolVec.clear();
			if( ! residueAtomsVec.empty() ){
				// save the last atom
				Residue				resid( tempAt.get_residueName(),
						tempAt.get_residueIndex(), residueAtomsVec );
				if( chainVec.empty() ){
					Chain				ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain				ch( resid );
						chainVec.push_back( ch );
					}
				}
				chain.clear();
				residueAtomsVec.clear();
			}

			// update bond information
			vector<Atom>					allAtomVec;
			for( size_t i=0; i<chainVec.size(); i++ ){
				vector<Atom>			chainAtVec = chainVec[i].get_AtomVec();
				allAtomVec.insert( allAtomVec.end(), chainAtVec.begin(), chainAtVec.end() );
			}

			for( size_t i=0; i<allAtomVec.size(); i++  ){
                for( size_t j=i+1; ( j<allAtomVec.size()   && j < i+50  ); j++ ){
					Atom a1 = allAtomVec[i];
					Atom a2 = allAtomVec[j];
					float	dist12 =  atomDis( a1, a2 ) ;
                    if( dist12 < maxBondLength && !(a1.get_type() == "H" && a2.get_type() == "H") ){
						size_t 	firstIndex = a1.get_index();
						size_t 	secondIndex = a2.get_index();
						string	bondType = "";
						Bond	bond( 0, a1, a2, bondType );
						bondVec.push_back( bond );
//						cout<<"$$$$$$$$$$$$$$$$$$$444"<<endl;
					}
				}
			}
			// set the previous residue pointer for each residue
			if( name == "" ){
				//				name = string( "no name" );
				name = fileName;
			}
			Molecular		mol( isEmpty, name, fileName, center,
					fileHeader, chainVec, bondVec, ss );
			if( ! mol.empty() ){
				cout<<"molecular: "<<name<<endl;
				molVec.push_back( mol );
			}else{
			}
			isEmpty = true;
			name  = "";
			center = Coord(0, 0, 0);
			fileHeader.clear();
			chainVec.clear();
			bondVec.clear();
			allAtomVec.clear();
		}
		if( !startAtom && readHead ){
			fileHeader.push_back( oneLine );
		}
	}

	//-------------------------------------------------------------------------
	if(  !chainVec.empty() ){
		atomsOneMolVec.clear();
		if( ! residueAtomsVec.empty() ){
			// save the last atom
			Residue				resid( tempAt.get_residueName(),
					tempAt.get_residueIndex(), residueAtomsVec );
			if( chainVec.empty() ){
				Chain				ch( resid );
				chainVec.push_back( ch );
			}else{
				bool				success = false;
				for( size_t i=0; i<chainVec.size(); i++ ){
					if( chainVec[i].get_name() == resid.get_chainName() ){
						chainVec[i].add_Residue( resid );

						success = true;
						break;
					}
				}
				if( !success ){
					Chain				ch( resid );
					chainVec.push_back( ch );

				}
			}
			chain.clear();
			residueAtomsVec.clear();
		}

		// update bond information
		vector<Atom>					allAtomVec;
		for( size_t i=0; i<chainVec.size(); i++ ){
			vector<Atom>			chainAtVec = chainVec[i].get_AtomVec();
			for( size_t j=0; j<chainAtVec.size(); j++ ){
				allAtomVec.push_back( chainAtVec[j] );
			}
		}

		for( size_t i=0; i<allAtomVec.size(); i++  ){
            for( size_t j=i+1; ( j<allAtomVec.size()   && j < i+50  ); j++ ){
				Atom a1 = allAtomVec[i];
				Atom a2 = allAtomVec[j];
				float	dist12 =  atomDis( a1, a2 ) ;
				if( dist12 < maxBondLength && !(a1.get_type() == "H" && a2.get_type() == "H") ){
//					size_t 	firstIndex = a1.get_index();
//					size_t 	secondIndex = a2.get_index();
					string	bondType = "";
					Bond	bond( 0, a1, a2, bondType );
					bondVec.push_back( bond );
				}
			}
		}
		// set the previous residue pointer for each residue
		if( name == "" ){
			//			name = string( "no name" );
			char c('/');
			int p = fileName.find_last_of(c);
			name = fileName.substr(p+1);
		}
		Molecular		mol( isEmpty, name, fileName, center,
				fileHeader, chainVec, bondVec, ss );
		if( ! mol.empty() ){
			cout<<"molecular: "<<name<<endl;
			molVec.push_back( mol );
		}else{
		}
		isEmpty = true;
		name  = "";
		center = Coord(0, 0, 0);
		fileHeader.clear();
		chainVec.clear();
		bondVec.clear();
		allAtomVec.clear();
	}
	//-------------------------------------------------------------------------

	if(  !chainVec.empty() ){
		cout<<"no data"<<endl;
	}
	if( molVec.empty() ){
		cerr<<"---warnning!!! not appropriate pdb file---"<<file<<endl;
	}
    cout<<"time:"<<float( clock() - t_start )/ CLOCKS_PER_SEC<<endl;
	return molVec;
}

vector<Molecular>
read_pdbqt(  const string& file ){
	ifstream iff( file.c_str() );
	if( !iff ){
		cerr<<"---error can not open "<<file<<"---"<<endl;
		throw;
	}
	vector<Molecular>							molVec;
	string												oneLine;
	vector<string>								strVec;
	bool 												startAtom = false;
	bool 												startBond = false;
	bool												readHead = true;
	bool												recordName = false;
	vector<Atom>								residueAtomsVec;
	vector<Residue> 							residueVec;
	size_t												preResidIndex 	= 10e10;
	size_t												newResidIndex = 10e10;
	Chain												chain;
	Residue*											preResidue;
	Atom												preAtom;
	Atom												tempAt;

	vector< pair< size_t, size_t > >		helixVec;			// pair< helix start residue index, helix end residue index >
	vector< pair< size_t, size_t > > 	sheetVec;		// pair< sheet start residue index, sheet end residue index >

	bool												isEmpty;
	string												name;
	string												fileName  = file ;
	Coord												center;
	vector<string >								fileHeader;
	vector<Chain>								chainVec;
	vector<Bond>									bondVec;
	vector<Atom>								atomsOneMolVec;
	SStructure										ss;						// secondary structure

	Chain												residueChain;
	residueChain.clear();

	clock_t												t1 = clock();
	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		string			head = tokenBySpace( oneLine.substr(0, 6) )[0];
		//		cout<<"head:"<<head<<endl;

		if( strVec[0] == "MODEL" ){
			name = strVec.back();
		}
		if( strVec[0] == "HEADER" ){
			name = strVec.back();
		}
		if( strVec[0] == "ATOM" ){
			startAtom = true;
		}

		if( strVec[0] == "HELIX" ){
			size_t		startResid = atoi( strVec[5].c_str() ) ;
			size_t		endResid = atoi( strVec[8].c_str() ) ;
			pair<size_t, size_t>	helix = make_pair( startResid, endResid );
			helixVec.push_back( helix );
		}else if( strVec[0] == "SHEET" ){
			size_t		startResid = atoi( strVec[6].c_str() ) ;
			size_t		endResid = atoi( strVec[9].c_str() ) ;
			pair<size_t, size_t>	sheet = make_pair( startResid, endResid );
			sheetVec.push_back( sheet );
		}else	if( ( head == "ATOM" || head == "HETATM" ) && strVec.size() >= 10 ){
			isEmpty = false;
			Atom at( oneLine );
			tempAt = at;
			atomsOneMolVec.push_back( at );
			//			at.print();

			if( at.get_residueIndex() != preAtom.get_residueIndex() && !residueAtomsVec.empty() ){
				Residue				resid( preAtom.get_residueName(), preAtom.get_residueIndex(), residueAtomsVec );

				if( chainVec.empty() ){
					Chain				ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain				ch( resid );
						chainVec.push_back( ch );
					}
					//					resid.print();
				}
				residueAtomsVec.clear();
				residueAtomsVec.push_back( at );
			}else{
				residueAtomsVec.push_back( at );
			}

			preResidIndex = newResidIndex;
			preAtom	= at ;
		}else if( strVec[0] == "CONECT" ){
			for( size_t i=2; i<strVec.size()-2; i++ ){
				size_t 	firstIndex = atoi( strVec[1].c_str() );
				size_t 	secondIndex = atoi( strVec[i].c_str() );
				Atom	a1 = searchAtom( atomsOneMolVec, firstIndex );
				Atom	a2 = searchAtom( atomsOneMolVec, secondIndex );
				string	bondType = "";
				Bond	bond( 0, a1, a2, bondType );
				if( find( bondVec.begin(), bondVec.end(), bond ) == bondVec.end() ){
					bondVec.push_back( bond );
				}
			}
		}else if( strVec[0] == "TER" ){
			if( ! residueAtomsVec.empty() ){
				// save the last atom
				Residue				resid( tempAt.get_residueName(), tempAt.get_residueIndex(), residueAtomsVec );
				if( chainVec.empty() ){
					Chain				ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain				ch( resid );
						chainVec.push_back( ch );
					}
				}
				chain.clear();
				residueAtomsVec.clear();
				//				residueAtomsVec.push_back( at );
			}
		}else if (strVec[0]=="ENDMDL" || strVec[0] == "END"  ){
			atomsOneMolVec.clear();
			if( ! residueAtomsVec.empty() ){
				// save the last atom
				Residue				resid( tempAt.get_residueName(), tempAt.get_residueIndex(), residueAtomsVec );
				if( chainVec.empty() ){
					Chain				ch( resid );
					chainVec.push_back( ch );
				}else{
					bool				success = false;
					for( size_t i=0; i<chainVec.size(); i++ ){
						if( chainVec[i].get_name() == resid.get_chainName() ){
							chainVec[i].add_Residue( resid );
							success = true;
							break;
						}
					}
					if( !success ){
						Chain				ch( resid );
						chainVec.push_back( ch );
					}
				}
				chain.clear();
				residueAtomsVec.clear();
			}

			// update bond information
			vector<Atom>					allAtomVec;
			for( size_t i=0; i<chainVec.size(); i++ ){
				vector<Atom>			chainAtVec = chainVec[i].get_AtomVec();
				for( size_t j=0; j<chainAtVec.size(); j++ ){
					allAtomVec.push_back( chainAtVec[j] );
				}
			}

			for( size_t i=0; i<allAtomVec.size(); i++  ){
				for( size_t j=i+1; ( j<allAtomVec.size()   && j < i+50  ); j++ ){
					Atom a1 = allAtomVec[i];
					Atom a2 = allAtomVec[j];
					float	dist12 =  atomDis( a1, a2 ) ;

					if( dist12 < maxBondLength && !(a1.get_type() == "H" && a2.get_type() == "H") ){
						size_t 	firstIndex = a1.get_index();
						size_t 	secondIndex = a2.get_index();
						string	bondType = "";
						Bond	bond( 0, a1, a2, bondType );
						bondVec.push_back( bond );
					}
				}
			}
			if( name == "" ){
				name = string( "no name" );
			}
			Molecular		mol( isEmpty, name, fileName, center, fileHeader, chainVec, bondVec, ss );

			if( ! mol.empty() ){
				//				cout<<"molecular: "<<name<<endl;
				molVec.push_back( mol );
			}else{
				//				cout<<"empty mol"<<endl;
			}
			isEmpty = true;
			name  = "";
			center = Coord(0, 0, 0);
			fileHeader.clear();
			chainVec.clear();
			bondVec.clear();
			allAtomVec.clear();
		}
		if( !startAtom && readHead ){
			fileHeader.push_back( oneLine );
		}
	}

	if( ! chainVec.empty() ){
		// update bond information
		vector<Atom>					allAtomVec;
		for( size_t i=0; i<chainVec.size(); i++ ){
			vector<Atom>			chainAtVec = chainVec[i].get_AtomVec();
			for( size_t j=0; j<chainAtVec.size(); j++ ){
				allAtomVec.push_back( chainAtVec[j] );
			}
		}

		for( size_t i=0; i<allAtomVec.size(); i++  ){
			for( size_t j=i+1; ( j<allAtomVec.size()   && j < i+50  ); j++ ){
				Atom a1 = allAtomVec[i];
				Atom a2 = allAtomVec[j];
				float	dist12 =  atomDis( a1, a2 ) ;

				if( dist12 < maxBondLength && !(a1.get_type() == "H" && a2.get_type() == "H") ){
					size_t 	firstIndex = a1.get_index();
					size_t 	secondIndex = a2.get_index();
					string	bondType = "";
					Bond	bond( 0, a1, a2, bondType );
					bondVec.push_back( bond );
				}
			}
		}
		if( name == "" ){
			name = string( "no name" );
		}
		Molecular		mol( isEmpty, name, fileName, center, fileHeader, chainVec, bondVec, ss );
		if( ! mol.empty() ){
			cout<<"molecular: "<<name<<endl;
			molVec.push_back( mol );
		}
	}

	if(  !chainVec.empty() ){
		cout<<"no data"<<endl;
	}

	if( molVec.empty() ){
		cerr<<"---warnning!!! not appropriate pdb file---"<<endl;
	}

	return molVec;
}

vector<Atom>
Molecular::get_atomVec()const{
	vector<Atom>				allAtomVec;
	for( size_t i=0; i<chainVec.size(); i++ ){
		if( ! chainVec[i].get_AtomVec().empty() ){
			vector<Atom>		chAtVec = chainVec[i].get_AtomVec();
			for( size_t j=0; j<chAtVec.size(); j++ ){
				allAtomVec.push_back( chAtVec[j] );
			}
		}
	}

	return allAtomVec;
}

vector<Atom>
Molecular::get_atomVec(  const string& name ){
	vector<Atom>				allAtomVec;

	vector<Atom> atVec;
	allAtomVec = get_atomVec();
	for( size_t i=0; i<allAtomVec.size(); i++ ){
		if( allAtomVec[i].get_name() == name ){
			atVec.push_back( allAtomVec[i] );
		}
	}

	return atVec;
}

Atom
Molecular::get_atom(const size_t& index)const{

	vector<Atom> allAtomVec = get_atomVec();
	for( size_t i=0; i<allAtomVec.size(); i++ ){
		if( allAtomVec[i].get_index() == index ){
			return allAtomVec[i];
		}
	}
	throw "can not find atom";
}

vector<Bond>
Molecular::get_rotableBondVec( )const{
	vector<Bond>						subBondVec;
//	vector<Bond>						bdVec = this->get_bondVec();
	vector< vector<size_t> >            ringAtIndexVec = this->get_ringAtomIndex();

	for( size_t i=0; i<bondVec.size(); i++ ){
		if( bondVec[i].get_type() == "1" && bondVec[i].get_firstAtom().get_type() != "H" &&
				bondVec[i].get_secAtom().get_type() != "H" ){

			int count = 0;
			//erase rotable bonds in ring
			for( size_t m=0; m<ringAtIndexVec.size(); m++ ){
				count = 0;
				for( size_t n=0; n<ringAtIndexVec[m].size(); n++ ){
					if( bondVec[i].get_firstAtom().get_index() == ringAtIndexVec[m][n] ){
						count ++;
					}
					if( bondVec[i].get_secAtom().get_index() == ringAtIndexVec[m][n] ){
						count ++;
					}
				}
				if( count == 2 ){
					break;
				}
			}

			//erase terminate rotable bond
			int countA1=0;
			int countA2=0;
			Atom a1=bondVec[i].get_firstAtom();
			Atom a2=bondVec[i].get_secAtom();
			for( size_t m=0; m<bondVec.size(); m++ ){
				if( m !=i ){
					if( a1==bondVec[m].get_firstAtom() || a1==bondVec[m].get_secAtom() ){
						countA1++;
					}
					if( a2==bondVec[m].get_firstAtom() || a2==bondVec[m].get_secAtom() ){
						countA2++;
					}
				}
			}

			//erase rotable bond with H as terminal
			int c1=0;
			int c2=0;
			for( size_t m=0; m<bondVec.size(); m++ ){
				if( m !=i ){
					if( ( a1==bondVec[m].get_firstAtom() && bondVec[m].get_secAtom().get_type() != string("H") ) ||
						( a1==bondVec[m].get_secAtom() && bondVec[m].get_firstAtom().get_type() != string("H") ) ){
						c1++;
					}
					if( ( a2==bondVec[m].get_firstAtom() && bondVec[m].get_secAtom().get_type() != string("H") ) ||
						( a2==bondVec[m].get_secAtom() && bondVec[m].get_firstAtom().get_type() != string("H") ) ){
						c2++;
					}
				}
			}

			if( count != 2 && countA1 != 0 && countA2 !=0 && c1 != 0 && c2 != 0  ){
				subBondVec.push_back( bondVec[i] );
			}
		}
	}

	return											subBondVec ;
}

vector<Residue>
Molecular::get_residueVec()const{
	vector<Residue>			resVec;
	resVec.clear();
	for( size_t i=0; i<chainVec.size(); i++ ){
//		if( ! chainVec[i].empty() ){
			vector<Residue>	tempRes = chainVec[i].get_ResidueVec() ;
            resVec.insert( resVec.end(), tempRes.begin(), tempRes.end() );
//		}
	}
	return								resVec;
}

/**
 * get the atom closet to the 'center'
 */
Atom
Molecular::get_centerAtom()const{
	vector<Atom> atVec=get_atomVec();
	float minDis = 10e10;
	int   index=0;
	for( size_t i=0; i<atVec.size(); i++ ){
		float dis = getCoordDis( center, atVec[i].get_coord() );
		if( dis < minDis ){
			minDis = dis;
			index=i;
		}
	}
	return atVec[index];
}

/**
 *  get the largest distance between all atoms and the center atom
 */
float
Molecular::get_centerAtomRadius()const{
	Atom centAt=get_centerAtom();
	vector<Atom> atVec=get_atomVec();
	float maxDis = -10e10;
	for( size_t i=0; i<atVec.size(); i++ ){
		float dis = getCoordDis( centAt.get_coord(), atVec[i].get_coord() );
		maxDis = maxDis>dis?maxDis:dis;
	}
	return maxDis;
}

void
Molecular::print()const{
//	cout<<"REMARK Molecular: "<<fileName<<endl;
	for( size_t i=0; i<chainVec.size(); i++ ){
		vector<Residue>		resVec = chainVec[i].get_ResidueVec();
		for( size_t j=0; j<resVec.size(); j++ ){
			vector<Atom>		atomVec = resVec[j].get_atomVec();
			for( size_t k=0; k<atomVec.size(); k++ ){
				atomVec[k].print();
			}
		}
	}
	cout<<"ENDMDL"<<endl;
}

void
Molecular::print( string& fileName ){
	ofstream off( fileName.c_str() );
	off<<"REMARK Molecular: "<<name<<endl;
	for( size_t i=0; i<chainVec.size(); i++ ){
		vector<Residue>		resVec = chainVec[i].get_ResidueVec();
		for( size_t j=0; j<resVec.size(); j++ ){
			vector<Atom>		atomVec = resVec[j].get_atomVec();
			for( size_t k=0; k<atomVec.size(); k++ ){
//				atomVec[k].print();
				Coord coord = atomVec[k].get_coord();
				string type = atomVec[k].get_type();
				off.width(6);
				off<<left<<"ATOM";
				off.width(5);
				off <<right<< atomVec[k].get_index()<<"  ";
				off.width(4);
				off <<left<< atomVec[k].get_name();
				off.width(3);
				off<<right<<atomVec[k].get_residueName()<<" ";
				off.width(1);
				off<<left<<atomVec[k].get_chainName() ;
				off.width(4);
				off<<right<<atomVec[k].get_residueIndex();
				off.width(3);
				off<<"";
				off.width(9);
				off <<right<< float2string(coord.x).substr( 0, 7 );
				off.width(8);
				off <<right<< float2string(coord.y).substr( 0, 7 );
				off.width(8);
				off <<right<< float2string(coord.z).substr( 0, 7 );
				off.width(24);
				off <<right<<type.substr(0, type.find('.'))<<endl;
			}
		}
	}
	off<<"ENDMDL"<<endl;
	off.close();
}

void
Molecular::printPreResidue()const{
	cout<<"Molecular:"<<fileName<<" preResidues"<<endl;
	for( size_t i=0; i<chainVec.size(); i++ ){
		vector<Residue>		resVec=chainVec[i].get_ResidueVec();
		for( size_t j=0; j<resVec.size(); j++ ){
			cout<<"i:"<<i<<" j:"<<j<<endl;
			Residue*					pRes = resVec[j].get_preResidue();
			if( pRes != NULL ){
				vector<Atom>		atomVec = pRes->get_atomVec();
				for( size_t k=0; k<atomVec.size(); k++ ){
					atomVec[k].print();
				}
			}else{
				cout<<"NULL"<<endl;
			}
		}
	}
}

int
Molecular::computeCenter(){
	center.x = center.y = center.z = 0;
	int count = 0;
	vector<Atom>		atomVec;

	for( size_t i=0; i<atomVec.size(); i++ ){
		center.x += atomVec[i].get_coord().x;
		center.y += atomVec[i].get_coord().y;
		center.z += atomVec[i].get_coord().z;
		count ++;
	}
	if( atomVec.empty() ){
		count = 1;
	}
	center.x /= count;
	center.y /= count;
	center.z /= count;
}

void
Molecular::set_preResidue(){
	for( size_t i=0; i<chainVec.size(); i++ ){
		vector<Residue>		rsVec = chainVec[i].get_ResidueVec();
		if( rsVec.size() == 1 ){
			rsVec[0].set_preResidue( NULL );
		}
		for( size_t j=1; j<rsVec.size(); j++ ){
			rsVec[j].set_preResidue( &rsVec[j-1] );
			//			rsVec[j].get_preResidue()->get_C().print();
		}
		chainVec[i].set_ResidueVec( rsVec );
	}
}

bool
Molecular::innerCrash(float tolDis )const{
	vector<Bond> bdVec = get_bondVec();
	vector<Atom> atVec = get_atomVec();

	for( size_t i=0; i<atVec.size()-1; i++ ){
		for( size_t j=i+1; j<atVec.size(); j++ ){
			Atom a1=atVec[i];
			Atom a2=atVec[j];
			bool flag = false;
			for( size_t k=0; k<bdVec.size(); k++ ){
				Atom b1=bdVec[k].get_firstAtom();
				Atom b2=bdVec[k].get_secAtom();
				if( (a1==b1 && a2==b2) || (a1==b2 && a2==b1) ){
					flag=true;
				}
			}
			if( !flag && ( getCoordDis(a1.get_coord(),a2.get_coord())
					< (a1.get_vdwRadius()+a2.get_vdwRadius() + tolDis) ) ){
				a1.print();
				a2.print();
				cout<<"coordDis:"<<getCoordDis(a1.get_coord(),a2.get_coord())<<" ";
				cout<<" sum:"<<(a1.get_vdwRadius()+a2.get_vdwRadius())<<endl;
				return true;
			}
		}
	}
	return false;
}

void
Molecular::update_center(){
	vector<Atom> atVec = get_atomVec();
	Coord c(0,0,0);
	for( size_t i=0; i<atVec.size(); i++ ){
		c.x=c.x+atVec[i].get_coord().x;
		c.y=c.y+atVec[i].get_coord().y;
		c.z=c.z+atVec[i].get_coord().z;
	}
	if( atVec.size() != 0 ){
		c=c/(float)atVec.size();
	}
	center = c;
}

void
Molecular::update_atomCoord(  const Atom& at   ){
	for( size_t i=0; i<chainVec.size(); i++ ){
		chainVec[i].update_atomCoord( at );
	}
	for( size_t i=0; i<bondVec.size(); i++ ){
		if( at.get_index()==bondVec[i].get_firstAtomIndex() ){
			bondVec[i].set_fistAtomCoord(at.get_coord());
//			bondVec[i].get_secAtom().print();
//			cout<<"    update:";
//			at.print();
		}
		if( at.get_index()==bondVec[i].get_secAtomIndex() ){
			bondVec[i].set_secAtomCoord(at.get_coord());
//			cout<<"    update:";
//			bondVec[i].get_firstAtom().print();
//			at.print();
		}
	}
}

void
Molecular::update_atomCoordVec(  const vector<Atom>& atVec   ){
	for( size_t j=0; j<atVec.size(); j++ ){
		for( size_t i=0; i<chainVec.size(); i++ ){
			chainVec[i].update_atomCoord( atVec[j] );
		}
	}

	size_t bdSize=bondVec.size();
	for( size_t i=0; i<bdSize; i++ ){
		size_t atSize = atVec.size();
		for( size_t j=0; j<atSize; j++ ){
			if( bondVec[i].get_firstAtom().get_index() == atVec[j].get_index() ){
				bondVec[i].set_fistAtomCoord(atVec[j].get_coord());
			}
			if( bondVec[i].get_secAtom().get_index() == atVec[j].get_index() ){
				bondVec[i].set_secAtomCoord(atVec[j].get_coord());
			}
		}
	}
}

void
Molecular::calculateAlphaHelices(   ){
	ss.helixAtomsVecs.clear();
	Chain							chain;
	vector<Residue>		resVec;
	ss.helixAtomsVecs.clear();
	for( size_t i=0; i<chainVec.size(); i++ ){
		chain 	= chainVec[i];
		resVec = chain.get_ResidueVec();
		vector<Atom>		helixAtVec;
		Atom						preAtom;
		Atom						newAtom;
		size_t 						helixAtomCount = 1;
		for( size_t j=1; j<resVec.size(); j++ ){
			for( size_t k=j+2; k<=j+5&& k<resVec.size(); k++ ){
				if( resVec[j].get_name() != "HOH" && resVec[k].get_name() != "HOH" ){
					//					cout<<"helices j:"<<j<<" k:"<<k<<endl;
					float		energy = Residue::calculateHBondEnergyDSSP( resVec[j], resVec[j-1], resVec[k], resVec[k-1] );
					if( energy < kMaxHBondEnergy ){
						//						cout<<resVec[j].get_index()<<","<<resVec[j].get_name()<<" "<<resVec[k].get_index()<<","<<resVec[k].get_name()<<": "<<energy<<endl;
						preAtom = newAtom;
						newAtom = resVec[j].get_C();
						if( newAtom.get_residueIndex() == preAtom.get_residueIndex() || newAtom.get_residueIndex() == ( preAtom.get_residueIndex()+1) ){
							if( newAtom.get_residueIndex() == ( preAtom.get_residueIndex()+1) ){
								helixAtomCount++;
							}
						}else if(  ! helixAtVec.empty() && helixAtomCount > 4 ){
							ss.helixAtomsVecs.push_back( helixAtVec );
							helixAtVec.clear();
							helixAtomCount = 1;
						}else{
							helixAtVec.clear();
							helixAtomCount = 1;
						}
						for( size_t m=j+1 ; m<=k; m++ ){
							if( find_if( helixAtVec.begin(), helixAtVec.end(), AtomEqual(resVec[m].get_C()) ) == helixAtVec.end() ){
								helixAtVec.push_back( resVec[m].get_C() );
							}
						}
					}
				}
			}
			sort( helixAtVec.begin(), helixAtVec.end(), SortByAtomIndex() );
			if( ! helixAtVec.empty() ){	}
		}
		if(  ! helixAtVec.empty() && helixAtomCount > 4 ){
			ss.helixAtomsVecs.push_back( helixAtVec );
		}
	}

	//	for( size_t i=0; i<ss.helixAtomsVecs.size(); i++ ){
	//		cout<<"helix:"<<i<<endl;
	//		for( size_t j=0; j<ss.helixAtomsVecs[i].size(); j++ ){
	//			ss.helixAtomsVecs[i][j].print();
	//		}
	//	}
}

void
Molecular::calculateBetaSheets(){
	Chain							chain;
	vector<Residue>		resVec;
	ss.sheetAtomsVecs.clear();
	for( size_t i=0; i<chainVec.size(); i++ ){
		chain 	= chainVec[i];
		resVec = chain.get_ResidueVec();
		vector<Atom>		sheetAtVec;

		vector<size_t>		hbResVec;
		Residue					hbRes;
		vector<size_t>	 	preHbResVec;
		Residue					preHbRes;
		bool						oneStep = false;
		for( int j=1; j<resVec.size(); j++ ){
			preHbRes = hbRes;
			preHbResVec = hbResVec;
			oneStep = false;
			hbResVec.clear();
			for( int k=1; k<resVec.size(); k++ ){
				if( abs( j - k ) > 3 ){
					if( resVec[j].get_name() != "HOH" && resVec[k].get_name() != "HOH" ){
						//						cout<<"sheet j:"<<j<<" k:"<<k<<endl;
						float		energy = Residue::calculateHBondEnergyDSSP( resVec[j], resVec[j-1], resVec[k], resVec[k-1] );
						if( energy < kMaxHBondEnergy ){
							//							cout<<resVec[j].get_index()<<","<<resVec[j].get_name()<<" "<<resVec[k].get_index()<<","<<resVec[k].get_name()<<": "<<energy<<endl;
							hbRes = resVec[j];
							hbResVec.push_back( resVec[k].get_index() );
						}
					}
				}
			}
			if( ( !hbResVec.empty() && !preHbResVec.empty() ) && ( hbRes.get_index() == ( preHbRes.get_index()+1 ) ) ){
				for( size_t m=0; m<hbResVec.size(); m++ ){
					for( size_t n=0; n<preHbResVec.size(); n++ ){
						if( hbResVec[m] == preHbResVec[n]+1 || hbResVec[m] == preHbResVec[n]-1 ){
							oneStep = true;
						}
					}
				}

				if( oneStep ){
					sheetAtVec.push_back( resVec[j].get_C()  );
				}
			}
			if( !sheetAtVec.empty()  && ( hbRes.get_index() != preHbRes.get_index()+1  || !oneStep ) ){
				if( sheetAtVec.size() > 3  ){
					ss.sheetAtomsVecs.push_back( sheetAtVec );
				}
				sheetAtVec.clear();
				hbResVec.clear();
			}
		}
	}

}

void
Molecular::calculateTurns( ){
	ss.turnAtomsVecs.clear();
    Chain				chain;
	vector<Residue>		resVec;

	bool							isHelix = false ;
	bool							isSheet = false ;
	Atom							preAtom;

	for( size_t i=0; i<chainVec.size(); i++ ){
		chain = chainVec[i];
		resVec = chain.get_ResidueVec();
		vector<Atom>			turnAtVec;

		for( size_t j=0; j<resVec.size(); j++ ){
			isHelix = false ;
			isSheet = false ;

			Atom		at = resVec[j].get_C();
			//			at.print();
			if( Residue::isAminoAcid( resVec[j] ) ){
				for( size_t m=0; m<ss.sheetAtomsVecs.size(); m++ ){
					for( size_t n=0; n<ss.sheetAtomsVecs[m].size(); n++ ){
						if( at == ss.sheetAtomsVecs[m][n] ){
							isSheet = true;

						}
					}
				}
				for( size_t m=0; m<ss.helixAtomsVecs.size(); m++ ){
					for( size_t n=0; n<ss.helixAtomsVecs[m].size(); n++ ){
						if( at == ss.helixAtomsVecs[m][n] ){
							isHelix = true;

						}
					}
				}
				if( !isSheet && ! isHelix ){
					if( at.get_residueIndex() != preAtom.get_residueIndex() + 1 && ! turnAtVec.empty()  ){

						Residue	reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
						Residue	reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
						if( !reFirst.empty() ){
							turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
						}
						if( !reEnd.empty() ){
							turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
						}

						reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
						reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
						if( !reFirst.empty() ){
							turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
						}
						if( !reEnd.empty() ){
							turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
						}

						reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
						reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
						if( !reFirst.empty() ){
							turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
						}
						if( !reEnd.empty() ){
							turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
						}

						ss.turnAtomsVecs.push_back( turnAtVec );
						turnAtVec.clear();
						turnAtVec.push_back( at );
						preAtom = at;
					}else{
						turnAtVec.push_back( at );
						preAtom = at;
					}
				}else{

				}
			}
		}
		if( ! turnAtVec.empty() ){
			Residue	reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
			Residue	reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
			if( !reFirst.empty() ){
				turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
			}
			if( !reEnd.empty() ){
				turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
			}

			reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
			reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
			if( !reFirst.empty() ){
				turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
			}
			if( !reEnd.empty() ){
				turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
			}

			reFirst= Residue::getResidue( resVec,  turnAtVec.front().get_residueIndex() - 1 );
			reEnd= Residue::getResidue( resVec,  turnAtVec.back().get_residueIndex() + 1 );
			if( !reFirst.empty() ){
				turnAtVec.insert( turnAtVec.begin(), reFirst.get_C() );
			}
			if( !reEnd.empty() ){
				turnAtVec.insert( turnAtVec.end(), reEnd.get_C() );
			}
			ss.turnAtomsVecs.push_back( turnAtVec );
		}
	}
	//	cout<<"turn:"<<ss.turnAtomsVecs.size()<<endl;
	//	for( size_t i=0; i<ss.turnAtomsVecs.size(); i++ ){
	//		cout<<"turn:"<<i<<endl;
	//		for( size_t j=0; j<ss.turnAtomsVecs[i].size(); j++ ){
	//			ss.turnAtomsVecs[i][j].print();
	//		}
	//	}
}

vector<vector<size_t> >
Molecular::get_ringAtomIndex()const{
	vector<Atom>								atomVec = get_atomVec();
	vector<pair<size_t, bool> >			visitedVec( atomVec.size() );
	vector<size_t>								historyVec;
	vector<vector<size_t> >				atomInRing;

	for( size_t i=0; i<atomVec.size(); i++ ){
		pair<size_t, bool>				visit;
		visit = make_pair( atomVec[i].get_index(), false );
		visitedVec[i] = visit;
	}

	vector<vector<size_t> >				ringIndexVec;

	size_t			startIndex = 0;
	for( size_t i=0; i<visitedVec.size(); i++ ){

		if( visitedVec[i].second == false ){ startIndex = visitedVec[i].first; break; }
	}
	dfsMol( bondVec,  startIndex, 10e10,visitedVec, historyVec, ringIndexVec );

	return	ringIndexVec;
}

vector<size_t>
Molecular::get_methyC(){
	vector<Atom>							atomVec = get_atomVec();
	vector<size_t>							methyVec;
	for( size_t i=0; i<atomVec.size(); i++ ){
		size_t 										count = 0;
		for( size_t m=0; m<bondVec.size(); m++ ){
			if( bondVec[m].get_firstAtomIndex() == atomVec[i].get_index()  ){
				Atom								atSec = searchAtom( atomVec, bondVec[m].get_secAtomIndex() );
				if( atSec.get_type() == "H" ){
					count++;
				}
			}
			if( bondVec[m].get_secAtomIndex() == atomVec[i].get_index()  ){
				Atom								atSec = searchAtom( atomVec, bondVec[m].get_firstAtomIndex() );
				if( atSec.get_type() == "H" ){
					count++;
				}
			}
		}
		if( count == 3 ){
			methyVec.push_back( atomVec[i].get_index() );
			//			cout<<"methy:"<<atomVec[i].get_index()<<endl;
		}
	}
	return		methyVec;
}

/**
 * bonds types 1 and that contains no leaf atom and not contained in a ring, not methy C
 * @return
 */
vector<Bond>
Molecular::constrainedRotatableBond(){
	vector<Bond>							rotBondVec;
	vector<vector<size_t> >		ringVec = get_ringAtomIndex() ;
	vector<size_t>						methyC = get_methyC();
	//	for( size_t i=0; i<ringVec.size(); i++ ){
	//		for( size_t j=0; j<ringVec[i].size(); j++ ){
	//			cout<<ringVec[i][j]<<" ";
	//		}
	//		cout<<endl;
	//	}
	for( size_t i=0; i<bondVec.size(); i++ ){
		if( bondVec[i].get_type() == "1" ){
			size_t count1 = 0;
			size_t count2 = 0;
			for( size_t j=0; j<bondVec.size(); j++ ){
				if( bondVec[i].get_firstAtomIndex() == bondVec[j].get_firstAtomIndex() ||
						bondVec[i].get_firstAtomIndex() == bondVec[j].get_secAtomIndex()	){
					count1 ++;
				}
				if( bondVec[i].get_secAtomIndex() == bondVec[j].get_firstAtomIndex() ||
						bondVec[i].get_secAtomIndex() == bondVec[j].get_secAtomIndex()	){
					count2 ++;
				}
			}
			bool		inRing = false;
			for( size_t m=0; m<ringVec.size(); m++ ){
				size_t 	count = 0;
				for( size_t n=0; n<ringVec[m].size(); n++ ){
					if( bondVec[i].get_firstAtomIndex() == ringVec[m][n] ){
						count++;
					}
					if( bondVec[i].get_secAtomIndex() == ringVec[m][n] ){
						count++;
					}
				}
				if( count == 2 ){
					inRing = true;
				}
			}
			bool		notMethyC = true;
			for( size_t k=0; k<methyC.size(); k++ ){
				if( methyC[k] == bondVec[i].get_firstAtomIndex() ||
						methyC[k] == bondVec[i].get_secAtomIndex() ){
					notMethyC = false;
				}
			}
			if( count1 > 1 && count2 > 1 && ! inRing && notMethyC ){
				rotBondVec.push_back( bondVec[i] );
			}
		}
	}
	return rotBondVec;
}

Molecular&
Molecular::operator =( const Molecular& mol ){
	if( this != &mol ){
		isEmpty = mol.empty();
		name = mol.get_name();
		fileName = mol.get_fileName();
		center = mol.get_center();

		vector<string> head( mol.get_fileHeader().size() );
		fileHeader = head;
		for( size_t i=0; i<mol.get_fileHeader().size(); i++ ){
            fileHeader[i] = mol.get_fileHeader()[i];
		}

        vector<Chain> mol_chainVec = mol.get_chainVec();
        vector<Chain> ch( mol_chainVec.size() );
		chainVec = ch;
		for( size_t i=0; i<chainVec.size();i++ ){
            chainVec[i] = mol_chainVec[i];
		}

        vector<Bond> mol_bdVec = mol.get_bondVec();
        vector<Bond> vb( mol_bdVec.size() );
		bondVec = vb;
		for( size_t i=0; i<bondVec.size(); i++ ){
            bondVec[i] = mol_bdVec[i];
		}

		ss = mol.get_secondaryStructure();

		//	cout<<"="<<endl;
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Atom>
overlay( const Molecular& standMol, const Molecular& oldMol, const bool& debug ){

    vector< vector<double> >	rotMatrix( 3, vector<double>(3, 0) );
    vector<Atom>						atomStd = standMol.get_atomVec();
    vector<Atom>						atomOld = oldMol.get_atomVec();
    vector<Atom>						atomNew = atomOld;
    vector<float>						coordStd;
    vector<float>						coordOld;
    vector<float>						translateVec( 3, 0 );
    vector<size_t> 						olAtomIndex;
    vector<Atom>						overlayAtomStd, overlayAtomOld;

    for( size_t i=0; i<atomStd.size(); i++ ){
        olAtomIndex.push_back( atomStd[i].get_index() );
    }

    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        overlayAtomStd.push_back( searchAtom( atomStd, olAtomIndex[i] ) );
        overlayAtomOld.push_back( searchAtom( atomOld, olAtomIndex[i] ) );
    }
    Coord						cen1 = centerOfAtoms( overlayAtomStd );
    Coord						cen2 = centerOfAtoms( overlayAtomOld );
    translateVec[0] = cen1.x - cen2.x;
    translateVec[1] = cen1.y - cen2.y;
    translateVec[2] = cen1.z - cen2.z;

    if(debug){
        cout<<"translate vector: "<<translateVec[0]<<" "
                <<translateVec[1]<<" "<<translateVec[2]<<endl;
    }

    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomStd , olAtomIndex[i] );
        coordStd.push_back( at.get_coord().x );
        coordStd.push_back( at.get_coord().y );
        coordStd.push_back( at.get_coord().z );
    }
    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomOld , olAtomIndex[i] );
        coordOld.push_back( at.get_coord().x  );
        coordOld.push_back( at.get_coord().y  );
        coordOld.push_back( at.get_coord().z  );
    }

    rotMatrix = overlayMatrix( coordStd, coordOld, true );

    if(debug){
        cout<<"rotate matrix2:"<<endl;
        for( size_t i=0; i<rotMatrix.size(); i++ ){
            for( size_t j=0; j<rotMatrix[i].size(); j++ ){
                cout.width(15);
                cout<<rotMatrix[i][j];
            }
            cout<<endl;
        }
    }

    for( size_t i=0; i<atomOld.size(); i++ ){
        vector<float>			coordNewAtom(3, 0);
        vector<float>			coordOldAtom;
        coordOldAtom.push_back( atomOld[i].get_coord().x   - cen2.x );
        coordOldAtom.push_back( atomOld[i].get_coord().y   - cen2.y );
        coordOldAtom.push_back( atomOld[i].get_coord().z   - cen2.z );
        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNewAtom[m] += rotMatrix[m][n] * coordOldAtom[n];
            }
        }
        atomNew[i].set_coord( coordNewAtom[0]  + translateVec[0] + cen2.x,
                coordNewAtom[1]  + translateVec[1]  + cen2.y,
                coordNewAtom[2]  + translateVec[2]  + cen2.z);
    }
//    newMol.set_atomVec( atomNew );
    return			atomNew;
}

vector<Atom>
overlay( const vector<Atom> standAtVec, const vector<Atom> oldAtVec,
         const bool& debug ){
	cout<<0<<endl;
    vector< vector<double> >	rotMatrix( 3, vector<double>(3, 0) );
    vector<Atom>						atomStd = standAtVec;
    vector<Atom>						atomOld = oldAtVec;
    vector<Atom>						atomNew = atomOld;
    vector<float>						coordStd;
    vector<float>						coordOld;
    vector<float>						translateVec( 3, 0 );
    vector<size_t> 						olAtomIndex;
    vector<Atom>						overlayAtomStd, overlayAtomOld;

    for( size_t i=0; i<atomStd.size(); i++ ){
        olAtomIndex.push_back( atomStd[i].get_index() );
    }

    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        overlayAtomStd.push_back( searchAtom( atomStd, olAtomIndex[i] ) );
        overlayAtomOld.push_back( searchAtom( atomOld, olAtomIndex[i] ) );
    }
    Coord						cen1 = centerOfAtoms( overlayAtomStd );
    Coord						cen2 = centerOfAtoms( overlayAtomOld );
    translateVec[0] = cen1.x - cen2.x;
    translateVec[1] = cen1.y - cen2.y;
    translateVec[2] = cen1.z - cen2.z;

    if(debug){
        cout<<"translate vector: "<<translateVec[0]<<" "
                <<translateVec[1]<<" "<<translateVec[2]<<endl;
    }

    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomStd , olAtomIndex[i] );
        coordStd.push_back( at.get_coord().x );
        coordStd.push_back( at.get_coord().y );
        coordStd.push_back( at.get_coord().z );
    }
    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomOld , olAtomIndex[i] );
        coordOld.push_back( at.get_coord().x  );
        coordOld.push_back( at.get_coord().y  );
        coordOld.push_back( at.get_coord().z  );
    }

    rotMatrix = overlayMatrix( coordStd, coordOld, false );

    if(debug){
        cout<<"rotate matrix2:"<<endl;
        for( size_t i=0; i<rotMatrix.size(); i++ ){
            for( size_t j=0; j<rotMatrix[i].size(); j++ ){
                cout.width(15);
                cout<<rotMatrix[i][j];
            }
            cout<<endl;
        }
    }

    for( size_t i=0; i<atomOld.size(); i++ ){
        vector<float>			coordNewAtom(3, 0);
        vector<float>			coordOldAtom;
        coordOldAtom.push_back( atomOld[i].get_coord().x - cen2.x );
        coordOldAtom.push_back( atomOld[i].get_coord().y - cen2.y );
        coordOldAtom.push_back( atomOld[i].get_coord().z - cen2.z );
        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNewAtom[m] += rotMatrix[m][n] * coordOldAtom[n];
            }
        }
        atomNew[i].set_coord( coordNewAtom[0]  + translateVec[0] + cen2.x,
                coordNewAtom[1]  + translateVec[1]  + cen2.y,
                coordNewAtom[2]  + translateVec[2]  + cen2.z);
    }

    return			atomNew;
}

vector<Atom>
overlay(  const Molecular& standMol, const Molecular& oldMol,
        const vector<size_t>& olAtomIndex, const bool& debug ){

    vector< vector<double> >	rotMatrix( 3, vector<double>(3, 0) );
    vector<Atom>						atomStd = standMol.get_atomVec();
    vector<Atom>						atomOld = oldMol.get_atomVec();
    vector<Atom>						atomNew = atomOld;
    vector<float>						coordStd;
    vector<float>						coordOld;
    Molecular							newMol = oldMol;
    vector<float>						translateVec( 3, 0 );

    vector<Atom>						overlayAtomStd, overlayAtomOld;
    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        overlayAtomStd.push_back( searchAtom( atomStd, olAtomIndex[i] ) );
        overlayAtomOld.push_back( searchAtom( atomOld, olAtomIndex[i] ) );
    }
    Coord						cen1 = centerOfAtoms( overlayAtomStd );
    Coord						cen2 = centerOfAtoms( overlayAtomOld );
    translateVec[0] = cen1.x - cen2.x;
    translateVec[1] = cen1.y - cen2.y;
    translateVec[2] = cen1.z - cen2.z;



    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomStd , olAtomIndex[i] );
        coordStd.push_back( at.get_coord().x );
        coordStd.push_back( at.get_coord().y );
        coordStd.push_back( at.get_coord().z );
    }
    for( size_t i=0; i<olAtomIndex.size(); i++ ){
        Atom			at = searchAtom( atomOld , olAtomIndex[i] );
        coordOld.push_back( at.get_coord().x  );
        coordOld.push_back( at.get_coord().y  );
        coordOld.push_back( at.get_coord().z  );
    }

    rotMatrix = overlayMatrix( coordStd, coordOld, true );

    if(debug){
        cout<<"translate vector: "<<translateVec[0]<<" "<<translateVec[1]<<" "<<translateVec[2]<<endl;
        cen1.print();
        cen2.print();
        cout<<"rotate matrix2:"<<endl;
        for( size_t i=0; i<rotMatrix.size(); i++ ){
            for( size_t j=0; j<rotMatrix[i].size(); j++ ){
                cout.width(15);
                cout<<rotMatrix[i][j];
            }
            cout<<endl;
        }
    }

    /*
     * move to the center, rotate, then move back,
     * then translate to the center of standard Molecular
     */
    for( size_t i=0; i<atomOld.size(); i++ ){
        vector<float>			coordNewAtom(3, 0);
        vector<float>			coordOldAtom;
        coordOldAtom.push_back( atomOld[i].get_coord().x   - cen2.x );
        coordOldAtom.push_back( atomOld[i].get_coord().y   - cen2.y );
        coordOldAtom.push_back( atomOld[i].get_coord().z   - cen2.z );
        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNewAtom[m] += rotMatrix[m][n] * coordOldAtom[n];
            }
        }
        atomNew[i].set_coord( coordNewAtom[0]  + translateVec[0] + cen2.x,
                coordNewAtom[1]  + translateVec[1]  + cen2.y,
                coordNewAtom[2]  + translateVec[2]  + cen2.z);
    }

    return			atomNew;
}

Molecular
overlay(  const Molecular& standMol,
        const Molecular& oldMol,
        const string& overlayAtomName,
        const bool& debug ){

    vector< vector<double> >	rotMatrix( 3, vector<double>(3, 0) );
    vector<Atom>						atomStd = standMol.get_atomVec();
    vector<Atom>						atomOld = oldMol.get_atomVec();
    vector<Atom>						atomNew = atomOld;
    vector<float>						coordStd;
    vector<float>						coordOld;
    Molecular									newMol = oldMol;
    vector<float>						translateVec( 3, 0 );

    vector<Atom>						overlayAtomStd, overlayAtomOld;

    overlayAtomStd = searchAtom( atomStd, overlayAtomName );
    overlayAtomOld = searchAtom( atomOld, overlayAtomName );

    //------------------------------------------------------------------

    Coord						cen1 = centerOfAtoms( overlayAtomStd );
    Coord						cen2 = centerOfAtoms( overlayAtomOld );
    translateVec[0] = cen1.x - cen2.x;
    translateVec[1] = cen1.y - cen2.y;
    translateVec[2] = cen1.z - cen2.z;

    if(debug){
    	cout<<"overlayAtomStd:"<<overlayAtomStd.size()<<endl;
    	cout<<"overlayAtomOld:"<<overlayAtomOld.size()<<endl;
        cout<<"translate vector: "<<translateVec[0]<<" "
                <<translateVec[1]<<" "<<translateVec[2]<<endl;
    }

    for( size_t i=0; i<overlayAtomStd.size(); i++ ){
        Atom			at = overlayAtomStd[i];
        coordStd.push_back( at.get_coord().x );
        coordStd.push_back( at.get_coord().y );
        coordStd.push_back( at.get_coord().z );
    }
    for( size_t i=0; i<overlayAtomOld.size(); i++ ){
        Atom			at = overlayAtomOld[i];
        coordOld.push_back( at.get_coord().x  );
        coordOld.push_back( at.get_coord().y  );
        coordOld.push_back( at.get_coord().z  );
    }

//	rotMatrix = overlayMatrix( coordStd, coordOld, debug );

    if(1){
        coordStd.clear();
        coordOld.clear();
        for( size_t i=0; i<overlayAtomStd.size(); i++ ){
            Atom			at = overlayAtomStd[i];
            coordStd.push_back( at.get_coord().x - cen1.x );
            coordStd.push_back( at.get_coord().y - cen1.y );
            coordStd.push_back( at.get_coord().z - cen1.z );
        }
        for( size_t i=0; i<overlayAtomOld.size(); i++ ){
            Atom			at = overlayAtomOld[i];
            coordOld.push_back( at.get_coord().x - cen2.x );
            coordOld.push_back( at.get_coord().y - cen2.y );
            coordOld.push_back( at.get_coord().z - cen2.z );
        }

        matrix rotM;
        double rmsd;
        bool isz = false;
        rotM = rotMatrixByRMSD4C( coordStd, coordOld, isz, rmsd);
        size_t c = rotM.getColumnDimension();
        size_t r = rotM.getRowDimension();
        vector<double> array = rotM.getArray();

        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                rotMatrix[i][j] = array[i*3+j];
            }
        }

    }

    if(debug){
        cout<<"translate vector: "<<translateVec[0]<<" "<<translateVec[1]<<" "<<translateVec[2]<<endl;
        cen1.print();
        cen2.print();
        cout<<"rotate matrix2:"<<endl;
        for( size_t i=0; i<rotMatrix.size(); i++ ){
            for( size_t j=0; j<rotMatrix[i].size(); j++ ){
                cout.width(15);
                cout<<rotMatrix[i][j];
            }
            cout<<endl;
        }
    }

    /*
     * move to the center, rotate, then move back,
     * then translate to the center of standard Molecular
     */
    for( size_t i=0; i<atomOld.size(); i++ ){
        vector<float>			coordNewAtom(3, 0);
        vector<float>			coordOldAtom;
        coordOldAtom.push_back( atomOld[i].get_coord().x   - cen2.x );
        coordOldAtom.push_back( atomOld[i].get_coord().y   - cen2.y );
        coordOldAtom.push_back( atomOld[i].get_coord().z   - cen2.z );
        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNewAtom[m] += rotMatrix[m][n] * coordOldAtom[n];
            }
        }
        atomNew[i].set_coord( coordNewAtom[0]  + translateVec[0] + cen2.x,
                coordNewAtom[1]  + translateVec[1]  + cen2.y,
                coordNewAtom[2]  + translateVec[2]  + cen2.z);
//        cout<<(coordNewAtom[0]  + translateVec[0] + cen2.x)<<endl;
    }

    newMol.update_atomCoordVec( atomNew );
    return			newMol;
}

float
overlayRMSD( const Molecular MolA, const Molecular MolB ){
    vector<Atom> atVec1 = MolA.get_atomVec();
    vector<Bond> bdVec1 = MolA.get_bondVec();
    vector<Atom> atVecNoH1;
    atVecNoH1.clear();
    vector<string> atNameVec1;
    vector<pair<int, int> > connectVecNoH1;
    connectVecNoH1.clear();

    for( size_t i=0; i<atVec1.size(); i++ ){
        string type = atVec1[i].get_type();
        if( type.substr(0, type.find(".") ) != "H" ){
            atVecNoH1.push_back( atVec1[i] );
            atNameVec1.push_back( atVec1[i].get_name().substr(0,1) );
        }
    }

    for( size_t i=0; i<bdVec1.size(); i++ ){
        Atom a1 = bdVec1[i].get_firstAtom();
        Atom a2 = bdVec1[i].get_secAtom();
        int count=0;
        int index1,index2;

        for( size_t j=0; j<atVecNoH1.size(); j++ ){
            Atom a=atVecNoH1[j];
            if( a==a1 ){
                index1=j;
                count++;
            }
            if( a==a2 ){
                index2=j;
                count++;
            }
        }
        if( count == 2 ){
            connectVecNoH1.push_back( make_pair(index1, index2) );
        }
    }

    //--------------------------------------------------------------
    vector<Atom> atVec2 = MolB.get_atomVec();
    vector<Bond> bdVec2 = MolB.get_bondVec();
    vector<Atom> atVecNoH2;
    atVecNoH2.clear();
    vector<string> atNameVec2;
    vector<pair<int, int> > connectVecNoH2;
    connectVecNoH2.clear();

    for( size_t i=0; i<atVec2.size(); i++ ){
        string type = atVec2[i].get_type();
        if( type.substr(0, type.find(".") ) != "H" ){
            atVecNoH2.push_back( atVec2[i] );
            atNameVec2.push_back( atVec2[i].get_name().substr(0,1) );
        }
    }

    for( size_t i=0; i<bdVec2.size(); i++ ){
        Atom a1 = bdVec2[i].get_firstAtom();
        Atom a2 = bdVec2[i].get_secAtom();
        int count=0;
        int index1,index2;
        for( size_t j=0; j<atVecNoH2.size(); j++ ){
            Atom a=atVecNoH2[j];
            if( a==a1 ){
                index1=j;
                count++;
            }
            if( a==a2 ){
                index2=j;
                count++;
            }
        }
        if( count == 2 ){
            connectVecNoH2.push_back( make_pair(index1, index2) );
        }
    }

    //--------------------------------------------------------------
    vector< vector<int> > matrix = isomorphism( atNameVec1, connectVecNoH1,
                                                atNameVec2, connectVecNoH2);

    // get atom with new sequences
    vector<Atom> atVec2NewSeq;
    for( size_t i=0; i<matrix.size(); i++ ){
        for( size_t j=0; j<matrix[i].size(); j++ ){
            if( matrix[i][j] == 1 ){
                atVec2NewSeq.push_back( atVecNoH2[j] );
                break;
            }
        }
    }

    if(1){
        for( size_t i=0; i<atVecNoH1.size(); i++ ){
            atVecNoH1[i].print();
            atVec2NewSeq[i].print();
            cout<<endl;
        }
    }

    vector<Atom> atNewVec2 = overlay( atVecNoH1, atVec2NewSeq, false );

    float dis = atomVecRMSD( atVecNoH1, atNewVec2 );
    return dis;
}

float
overlayRMSD_CA( const Molecular molA, const Molecular molB ){

	vector<Residue> resVecA=molA.get_residueVec();
	vector<Residue> resVecB=molB.get_residueVec();

	vector<Atom>	atomVecA;
	vector<Atom>	atomVecB;

	for( size_t i=0; i<resVecA.size(); i++ ){
		if( i>resVecB.size() ){
			return -1;
		}
		Atom	caA = resVecA[i].get_CA();
		Atom	caB = resVecB[i].get_CA();

		atomVecA.push_back(caA);
		atomVecB.push_back(caB);
	}

    vector<Atom> atNewVecB = overlay( atomVecA, atomVecB, true ) ;

    float dis = atomVecRMSD( atomVecA, atNewVecB );
    return dis;
}

float
pairRMSD( const Molecular& molA,
		const Molecular& molB ){

	vector<Atom>	atomA = molA.get_atomVec();
	vector<Atom>	atomB = molB.get_atomVec();

	vector < size_t > atomIndexVec(atomA.size(), 0);
	for (size_t i = 0; i < atomA.size(); i++) {
		atomIndexVec[i] = atomA[i].get_index();
	}
	vector<vector<size_t> > allVec;
	allVec.push_back( atomIndexVec );

	vector<float>		distVec;
	float				dist = 0;
	size_t				heavyAtomNum = 0;
	for( size_t i=0; i<atomB.size(); i++ ){
		if( atomB[i].get_name().at(0) != 'H' && atomB[i].get_name() != "****"  ){
			heavyAtomNum ++;
		}
	}

	for( size_t i=0; i<allVec.size(); i++ ){
		dist = 0;
		for( size_t j=0; j<heavyAtomNum; j++ ){
			if( atomA[ allVec[i][j]-1 ].get_name().at(0) != atomB[j].get_name().at(0) ){
				cerr<<"---- can not compute distance "<<molA.get_fileName()<<" : "<<molB.get_fileName()<<endl;
				atomA[ allVec[i][j] ].print();
				atomB[j].print();
				cerr<<atomA[ allVec[i][j] ].get_name().at(0) <<" : "<< atomB[j].get_name().at(0)<<endl;
				return  -1;
			}
			float		x1 =  atomA[ allVec[i][j]-1 ].get_coord().x;
			float		y1 =  atomA[ allVec[i][j]-1 ].get_coord().y;
			float		z1 =  atomA[ allVec[i][j]-1 ].get_coord().z;
			float		x2	 = atomB[j].get_coord().x;
			float		y2	 = atomB[j].get_coord().y;
			float		z2	 = atomB[j].get_coord().z;

			dist += ( x1-x2 )*( x1-x2 ) + (y1-y2)*(y1-y2) + (z1-z2)*( z1-z2 );
		}
		if( heavyAtomNum == 0 ){
			heavyAtomNum = 1;
		}
		dist /= heavyAtomNum;
		distVec.push_back( dist );
	}

	dist = *min_element( distVec.begin(), distVec.end() );
	dist = sqrt( dist );
	return		dist;
}

float
pairRMSD_CA( const Molecular& molA,
		const Molecular& molB ){

	vector<Residue> resVecA=molA.get_residueVec();
	vector<Residue> resVecB=molB.get_residueVec();

	vector<Atom>	atomVecA;
	vector<Atom>	atomVecB;

	for( size_t i=0; i<resVecA.size(); i++ ){
		if( i>resVecB.size() ){
			return -1;
		}
		Atom	caA = resVecA[i].get_CA();
		Atom	caB = resVecB[i].get_CA();

		atomVecA.push_back(caA);
		atomVecB.push_back(caB);
	}
	float dist = atomVecRMSD( atomVecA, atomVecB );
	cout<<"pair CA RMSD:"<<dist<<endl;
	return		dist;
}

float
maxPairRMSD( const vector<Molecular> molVec ){
	float maxDis = -10e10;
	for( size_t i=0; i<molVec.size()-1; i++ ){
		for( size_t j=i+1; j<molVec.size(); j++ ){
			float dis = pairRMSD(molVec[i], molVec[j]);
			maxDis = maxDis<dis? dis:maxDis;
		}
	}
	return maxDis;
}

float
isomorphismRMSD( const Molecular MolA, const Molecular MolB, bool showMatrix ){
    vector<Atom> atVec1 = MolA.get_atomVec();
    vector<Bond> bdVec1 = MolA.get_bondVec();
    vector<Atom> atVecNoH1;
    atVecNoH1.clear();
    vector<string> atNameVec1;
    vector<pair<int, int> > connectVecNoH1;
    connectVecNoH1.clear();

    for( size_t i=0; i<atVec1.size(); i++ ){
        string type = atVec1[i].get_type();
        if( type.substr(0, type.find(".") ) != "H" ){
            atVecNoH1.push_back( atVec1[i] );
            atNameVec1.push_back( atVec1[i].get_name().substr(0,1) );
        }
    }
    cout<<1<<endl;

    for( size_t i=0; i<bdVec1.size(); i++ ){
        Atom a1 = bdVec1[i].get_firstAtom();
        Atom a2 = bdVec1[i].get_secAtom();
        int count=0;
        int index1,index2;

        for( size_t j=0; j<atVecNoH1.size(); j++ ){
            Atom a=atVecNoH1[j];
            if( a==a1 ){
                index1=j;
                count++;
            }
            if( a==a2 ){
                index2=j;
                count++;
            }
        }
        if( count == 2 ){
            connectVecNoH1.push_back( make_pair(index1, index2) );
        }
    }

    cout<<2<<endl;
    //--------------------------------------------------------------
    vector<Atom> atVec2 = MolB.get_atomVec();
    vector<Bond> bdVec2 = MolB.get_bondVec();
    vector<Atom> atVecNoH2;
    atVecNoH2.clear();
    vector<string> atNameVec2;
    vector<pair<int, int> > connectVecNoH2;
    connectVecNoH2.clear();

    for( size_t i=0; i<atVec2.size(); i++ ){
        string type = atVec2[i].get_type();
        if( type.substr(0, type.find(".") ) != "H" ){
            atVecNoH2.push_back( atVec2[i] );
            atNameVec2.push_back( atVec2[i].get_name().substr(0,1) );
        }
    }

    for( size_t i=0; i<bdVec2.size(); i++ ){
        Atom a1 = bdVec2[i].get_firstAtom();
        Atom a2 = bdVec2[i].get_secAtom();
        int count=0;
        int index1,index2;
        for( size_t j=0; j<atVecNoH2.size(); j++ ){
            Atom a=atVecNoH2[j];
            if( a==a1 ){
                index1=j;
                count++;
            }
            if( a==a2 ){
                index2=j;
                count++;
            }
        }
        if( count == 2 ){
            connectVecNoH2.push_back( make_pair(index1, index2) );
        }
    }
    cout<<3<<endl;
    //--------------------------------------------------------------
    vector< vector<int> > matrix = isomorphism( atNameVec1, connectVecNoH1,
                                                atNameVec2, connectVecNoH2);

    // RMSD
    float allSquareDis = 0;
    size_t matrixSize=matrix.size();
    for( size_t i=0; i<matrixSize; i++ ){

        Atom a1=atVecNoH1[i];
        float minDis = 10e10;
        size_t ind = 0;
        size_t	matrixSizeI = matrix[i].size();
        for( size_t j=0; j<matrixSizeI; j++ ){
            if( matrix[i][j] == 1 ){
                Atom a2 = atVecNoH2[j];
                float dis = atomSquareDis(a1, a2);
                if( dis < minDis ){
                    minDis = dis;
                    ind = j;
                }
            }
        }

        for( size_t k=0; k<matrixSize; k++ ){
            if( k != i ){
                matrix[k][ind] = 0;
            }
        }
        allSquareDis += minDis;
    }
    cout<<4<<endl;
    if(showMatrix){
        cout.width(8);
        cout<<atVecNoH2[0].get_name();
        for( size_t i=1; i<matrix.size(); i++ ){
            cout.width(4);
            cout<<atVecNoH2[i].get_name();
        }
        cout<<endl;
        for( size_t i=0; i<matrix.size(); i++ ){
            cout.width(4);
            cout<<atVecNoH1[i].get_name();
            for( size_t j=0; j<matrix.size(); j++ ){
                cout.width(4);
                cout<<matrix[i][j];
            }
            cout<<endl;
        }
    }
    cout<<5<<endl;
    return sqrt( allSquareDis / (float)matrix.size() );
}

bool
molecularCrash( 	const Molecular& mol1,
		const Molecular& mol2){
	vector<Atom>			atVec1 = mol1.get_atomVec() ;
	vector<Atom>			atVec2 = mol2.get_atomVec() ;

	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
			float dis = atomDis( atVec1[i], atVec2[j] );
			if( dis < (atVec1[i].get_vdwRadius() + atVec2[j].get_vdwRadius() ) ){
				return						true;
			}
		}
	}
	return false;
}

Molecular
rotateAroundBond( const Molecular& rotableMol,
		const Bond& bd,
		const float& degreeAngle ){

	Molecular				mol = rotableMol ;
	Atom					fixedAtom, rotableAtom;
	fixedAtom = bd.get_firstAtom();
	rotableAtom = bd.get_secAtom();

	Coord			rotateAxis ;
	rotateAxis.x = rotableAtom.get_coord().x - fixedAtom.get_coord().x ;
	rotateAxis.y = rotableAtom.get_coord().y - fixedAtom.get_coord().y ;
	rotateAxis.z = rotableAtom.get_coord().z - fixedAtom.get_coord().z ;
	double			length = sqrt( rotateAxis.x * rotateAxis.x +
			rotateAxis.y* rotateAxis.y + rotateAxis.z * rotateAxis.z );
	if( length == 0 ){
		cout<<"rotate around a point ?"<<endl;
		return mol;
	}

	Coord				unitAxis;
	unitAxis.x = rotateAxis.x / length;
	unitAxis.y = rotateAxis.y / length;
	unitAxis.z = rotateAxis.z / length;
	double				angle = degreeAngle * 3.1415926 / 180.0 ;
	double				sinAngle 		= sin( angle );
	double				cosAngle 		= cos( angle );
	double				vsAngle 		= 1 - cos( angle );

	vector< vector<double> >		rotateMatrix( 3, vector<double>(3, 0) );
	rotateMatrix[0][0]		= unitAxis.x * unitAxis.x * vsAngle + cosAngle;
	rotateMatrix[0][1]		= unitAxis.x * unitAxis.y * vsAngle  - unitAxis.z * sinAngle ;
	rotateMatrix[0][2]		= unitAxis.x * unitAxis.z * vsAngle + unitAxis.y * sinAngle ;
	rotateMatrix[1][0]		= unitAxis.x * unitAxis.y * vsAngle + unitAxis.z * sinAngle ;
	rotateMatrix[1][1]		= unitAxis.y * unitAxis.y * vsAngle + cosAngle;
	rotateMatrix[1][2]		= unitAxis.y * unitAxis.z * vsAngle  - unitAxis.x * sinAngle ;
	rotateMatrix[2][0]		= unitAxis.x * unitAxis.z * vsAngle - unitAxis.y * sinAngle;
	rotateMatrix[2][1]		= unitAxis.y * unitAxis.z * vsAngle + unitAxis.x * sinAngle;
	rotateMatrix[2][2]		= unitAxis.z * unitAxis.z * vsAngle + cosAngle;

	vector<Atom>			rotableAtomAll =
			dfsMol( rotableMol.get_atomVec(), rotableMol.get_bondVec() ,
					fixedAtom.get_index() , rotableAtom.get_index() );

	//delete fixedAtomIndex and rotableAtomIndex from rotableAtomVec
	vector<Atom>			rotableAtomVec;
	for( size_t i=0; i<rotableAtomAll.size(); i++ ){
		if( rotableAtomAll[i].get_index() != fixedAtom.get_index() &&
				rotableAtomAll[i].get_index() != rotableAtom.get_index() ){
			rotableAtomVec.push_back( rotableAtomAll[i] );
		}
	}

	for( size_t i=0; i<rotableAtomVec.size(); i++ ){
		vector<double>		coord(3, 0);
		coord[0]  = rotableAtomVec[i].get_coord().x  -  fixedAtom.get_coord().x;
		coord[1]  = rotableAtomVec[i].get_coord().y  -  fixedAtom.get_coord().y;
		coord[2]  = rotableAtomVec[i].get_coord().z  -  fixedAtom.get_coord().z;

		vector<double>		newCoord(3, 0);
		newCoord = matrixMultiplyVector( rotateMatrix, coord );
		newCoord[0] += fixedAtom.get_coord().x;
		newCoord[1] += fixedAtom.get_coord().y;
		newCoord[2] += fixedAtom.get_coord().z;
		rotableAtomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );

		//update
        vector<Atom> tempAtVec = rotableMol.get_atomVec();
        for( size_t j=0; j<tempAtVec.size(); j++ ){
            if( tempAtVec[j].get_index() == rotableAtomVec[i].get_index() ){
                Atom 		newAt = tempAtVec[j];
                newAt.set_coord( newCoord[0], newCoord[1], newCoord[2] );
                mol.update_atomCoord( newAt );
            }
        }
	}
	return 				mol;
}

bool
rotateAroundBond( const Molecular& rotableMol, const Bond& bd,
		const float& degreeAngle, Molecular& output ){

	Molecular				mol;
	mol = rotableMol ;
	Atom					fixedAtom, rotableAtom;
    vector<Atom>            rotableMol_atVec = rotableMol.get_atomVec();
//	fixedAtom = bd.get_firstAtom();
	fixedAtom = mol.get_atom(bd.get_firstAtomIndex() );
//	rotableAtom = bd.get_secAtom();
	rotableAtom = mol.get_atom( bd.get_secAtomIndex() );

	Coord			rotateAxis ;
	rotateAxis.x = rotableAtom.get_coord().x - fixedAtom.get_coord().x ;
	rotateAxis.y = rotableAtom.get_coord().y - fixedAtom.get_coord().y ;
	rotateAxis.z = rotableAtom.get_coord().z - fixedAtom.get_coord().z ;
	double			length = sqrt( rotateAxis.x * rotateAxis.x +
			rotateAxis.y* rotateAxis.y + rotateAxis.z * rotateAxis.z );
	if( length == 0 ){
		cout<<"rotate around a point ?"<<endl;
		return false;
	}

	Coord				unitAxis;
	unitAxis.x = rotateAxis.x / length;
	unitAxis.y = rotateAxis.y / length;
	unitAxis.z = rotateAxis.z / length;
	double				angle = degreeAngle * 3.1415926 / 180.0 ;
	double				sinAngle 		= sin( angle );
	double				cosAngle 		= cos( angle );
	double				vsAngle 		= 1 - cos( angle );

	vector< vector<double> >		rotateMatrix( 3, vector<double>(3, 0) );
	rotateMatrix[0][0]		= unitAxis.x * unitAxis.x * vsAngle + cosAngle;
	rotateMatrix[0][1]		= unitAxis.x * unitAxis.y * vsAngle  - unitAxis.z * sinAngle ;
	rotateMatrix[0][2]		= unitAxis.x * unitAxis.z * vsAngle + unitAxis.y * sinAngle ;
	rotateMatrix[1][0]		= unitAxis.x * unitAxis.y * vsAngle + unitAxis.z * sinAngle ;
	rotateMatrix[1][1]		= unitAxis.y * unitAxis.y * vsAngle + cosAngle;
	rotateMatrix[1][2]		= unitAxis.y * unitAxis.z * vsAngle  - unitAxis.x * sinAngle ;
	rotateMatrix[2][0]		= unitAxis.x * unitAxis.z * vsAngle - unitAxis.y * sinAngle;
	rotateMatrix[2][1]		= unitAxis.y * unitAxis.z * vsAngle + unitAxis.x * sinAngle;
	rotateMatrix[2][2]		= unitAxis.z * unitAxis.z * vsAngle + cosAngle;

	vector<Atom>			rotableAtomAll =
			dfsMol( rotableMol.get_atomVec(), rotableMol.get_bondVec() ,
					fixedAtom.get_index() , rotableAtom.get_index() );

	//delete fixedAtomIndex and rotableAtomIndex from rotableAtomVec
	vector<Atom>			rotableAtomVec;
	for( size_t i=0; i<rotableAtomAll.size(); i++ ){
		if( rotableAtomAll[i].get_index() != fixedAtom.get_index() &&
				rotableAtomAll[i].get_index() != rotableAtom.get_index() ){
			rotableAtomVec.push_back( rotableAtomAll[i] );
		}
	}

	for( size_t i=0; i<rotableAtomVec.size(); i++ ){
		vector<double>		coord(3, 0);
		coord[0]  = rotableAtomVec[i].get_coord().x  -  fixedAtom.get_coord().x;
		coord[1]  = rotableAtomVec[i].get_coord().y  -  fixedAtom.get_coord().y;
		coord[2]  = rotableAtomVec[i].get_coord().z  -  fixedAtom.get_coord().z;

		vector<double>		newCoord(3, 0);
		newCoord = matrixMultiplyVector( rotateMatrix, coord );
		newCoord[0] += fixedAtom.get_coord().x;
		newCoord[1] += fixedAtom.get_coord().y;
		newCoord[2] += fixedAtom.get_coord().z;
		rotableAtomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );

		//update
        for( size_t j=0; j<rotableMol_atVec.size(); j++ ){
            if( rotableMol_atVec[j].get_index() == rotableAtomVec[i].get_index() ){
                Atom 		newAt = rotableMol_atVec[j];
				newAt.set_coord( newCoord[0], newCoord[1], newCoord[2] );
				mol.update_atomCoord( newAt );
			}
		}
	}

	//check if new molecular has atoms crash
	vector<Bond> bdVec = mol.get_bondVec();
	vector<Atom> atVec = mol.get_atomVec();
	for( size_t i=0; i<atVec.size()-1; i++ ){
		for( size_t j=i+1; j<atVec.size(); j++ ){
			Atom a1 = atVec[i];
			Atom a2 = atVec[j];
			float dis = getCoordDis( a1.get_coord(), a2.get_coord() );
			if( dis < 1.65 ){
				bool flag = false;
				for( size_t k=0; k<bdVec.size(); k++ ){
					Atom b1 = bdVec[k].get_firstAtom();
					Atom b2 = bdVec[k].get_secAtom();
					if( ( a1.get_index()==b1.get_index() && a2.get_index()==b2.get_index() ) ||
							( a1.get_index()==b2.get_index()&&a2.get_index()==b1.get_index() ) ){
						flag = true;
						//						a1.print();
						//						a2.print();
						//						b1.print();
						//						b2.print();
						//						cout<<endl;
						break;
					}
				}
				if( !flag ){
					output=rotableMol;
					return flag;
				}
			}
		}
	}
	mol.update_center();
	output = mol;
	return true;
}

Molecular
rotateAroundBond( const Molecular& rotableMol, const size_t& fixedAtomIndex,
		const size_t& rotableAtomIndex, const float& degreeAngle ){
	Molecular				mol( rotableMol );

	Atom					fixedAtom, rotableAtom;
    vector<Atom>            rotableMol_atVec = rotableMol.get_atomVec();
    for( size_t i=0; i<rotableMol_atVec.size(); i++ ){
        if( rotableMol_atVec[i].get_index() == fixedAtomIndex ){
            fixedAtom = rotableMol_atVec[i];
		}
        if( rotableMol_atVec[i].get_index() == rotableAtomIndex ){
            rotableAtom = rotableMol_atVec[i];
		}
	}

	Coord			rotateAxis ;
	rotateAxis.x = rotableAtom.get_coord().x - fixedAtom.get_coord().x ;
	rotateAxis.y = rotableAtom.get_coord().y - fixedAtom.get_coord().y ;
	rotateAxis.z = rotableAtom.get_coord().z - fixedAtom.get_coord().z ;
	float			length = sqrt( rotateAxis.x * rotateAxis.x +
			rotateAxis.y* rotateAxis.y + rotateAxis.z * rotateAxis.z );
	if( length == 0 ){
		cout<<"rotate around a point ?"<<endl;
		return mol;
	}

	Coord				unitAxis;
	unitAxis.x = rotateAxis.x / length;
	unitAxis.y = rotateAxis.y / length;
	unitAxis.z = rotateAxis.z / length;
	float				angle = degreeAngle * 3.1415926 / 180.0 ;
	float				sinAngle 		= sin( angle );
	float				cosAngle 		= cos( angle );
	float				vsAngle 		= 1 - cos( angle );

	vector< vector<float> >		rotateMatrix( 3, vector<float>(3, 0) );
	rotateMatrix[0][0]		= unitAxis.x * unitAxis.x * vsAngle + cosAngle;
	rotateMatrix[0][1]		= unitAxis.x * unitAxis.y * vsAngle - unitAxis.z * sinAngle ;
	rotateMatrix[0][2]		= unitAxis.x * unitAxis.z * vsAngle + unitAxis.y * sinAngle ;
	rotateMatrix[1][0]		= unitAxis.x * unitAxis.y * vsAngle + unitAxis.z * sinAngle ;
	rotateMatrix[1][1]		= unitAxis.y * unitAxis.y * vsAngle + cosAngle;
	rotateMatrix[1][2]		= unitAxis.y * unitAxis.z * vsAngle  - unitAxis.x * sinAngle ;
	rotateMatrix[2][0]		= unitAxis.x * unitAxis.z * vsAngle - unitAxis.y * sinAngle;
	rotateMatrix[2][1]		= unitAxis.y * unitAxis.z * vsAngle + unitAxis.x * sinAngle;
	rotateMatrix[2][2]		= unitAxis.z * unitAxis.z * vsAngle + cosAngle;


	vector<Atom>			rotableAtomAll =
			dfsMol( rotableMol.get_atomVec(), rotableMol.get_bondVec() , fixedAtomIndex , rotableAtomIndex );

	//delete fixedAtomIndex and rotableAtomIndex from rotableAtomVec
	vector<Atom>			rotableAtomVec;
	for( size_t i=0; i<rotableAtomAll.size(); i++ ){
		if( rotableAtomAll[i].get_index() != fixedAtomIndex && rotableAtomAll[i].get_index() != rotableAtomIndex ){
			rotableAtomVec.push_back( rotableAtomAll[i] );
		}
	}

	for( size_t i=0; i<rotableAtomVec.size(); i++ ){
		vector<float>		coord(3, 0);
		coord[0]  = rotableAtomVec[i].get_coord().x  -  fixedAtom.get_coord().x;
		coord[1]  = rotableAtomVec[i].get_coord().y  -  fixedAtom.get_coord().y;
		coord[2]  = rotableAtomVec[i].get_coord().z  -  fixedAtom.get_coord().z;

		vector<float>		newCoord(3, 0);
		newCoord = matrixMultiplyVector( rotateMatrix, coord );
		newCoord[0] += fixedAtom.get_coord().x;
		newCoord[1] += fixedAtom.get_coord().y;
		newCoord[2] += fixedAtom.get_coord().z;
		rotableAtomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );

		//update
        for( size_t j=0; j<rotableMol_atVec.size(); j++ ){
            if( rotableMol_atVec[j].get_index() == rotableAtomVec[i].get_index() ){
                Atom 		newAt = rotableMol_atVec[j];
				newAt.set_coord( newCoord[0], newCoord[1], newCoord[2] );
				mol.update_atomCoord( newAt );
			}
		}
	}
	return 				mol;
}

bool
rotateAroundBond( const Molecular& rotableMol, const size_t& fixedAtomIndex,
		const size_t& rotableAtomIndex, const float& degreeAngle, Molecular& output ){
	Molecular				mol;
	mol =  rotableMol ;

	Atom					fixedAtom, rotableAtom;
    vector<Atom>            rotableMol_atVec = rotableMol.get_atomVec();
    for( size_t i=0; i<rotableMol_atVec.size(); i++ ){
        if( rotableMol_atVec[i].get_index() == fixedAtomIndex ){
            fixedAtom = rotableMol_atVec[i];
		}
        if( rotableMol_atVec[i].get_index() == rotableAtomIndex ){
            rotableAtom = rotableMol_atVec[i];
		}
	}

	Coord			rotateAxis ;
	rotateAxis.x = rotableAtom.get_coord().x - fixedAtom.get_coord().x ;
	rotateAxis.y = rotableAtom.get_coord().y - fixedAtom.get_coord().y ;
	rotateAxis.z = rotableAtom.get_coord().z - fixedAtom.get_coord().z ;
	float			length = sqrt( rotateAxis.x * rotateAxis.x +
			rotateAxis.y* rotateAxis.y + rotateAxis.z * rotateAxis.z );
	if( length == 0 ){
		cout<<"rotate around a point ???"<<endl;
		cout<<"fixed:"<<fixedAtomIndex<<" rotable:"<<rotableAtomIndex<<endl;
		fixedAtom.print();
		rotableAtom.print();
		rotableMol.print();
		cout<<endl;
		return false;
	}

	Coord				unitAxis;
	unitAxis.x = rotateAxis.x / length;
	unitAxis.y = rotateAxis.y / length;
	unitAxis.z = rotateAxis.z / length;
	float				angle = degreeAngle * 3.1415926 / 180.0 ;
	float				sinAngle 		= sin( angle );
	float				cosAngle 		= cos( angle );
	float				vsAngle 		= 1 - cos( angle );

	vector< vector<float> >		rotateMatrix( 3, vector<float>(3, 0) );
	rotateMatrix[0][0]		= unitAxis.x * unitAxis.x * vsAngle + cosAngle;
	rotateMatrix[0][1]		= unitAxis.x * unitAxis.y * vsAngle  - unitAxis.z * sinAngle ;
	rotateMatrix[0][2]		= unitAxis.x * unitAxis.z * vsAngle + unitAxis.y * sinAngle ;
	rotateMatrix[1][0]		= unitAxis.x * unitAxis.y * vsAngle + unitAxis.z * sinAngle ;
	rotateMatrix[1][1]		= unitAxis.y * unitAxis.y * vsAngle + cosAngle;
	rotateMatrix[1][2]		= unitAxis.y * unitAxis.z * vsAngle  - unitAxis.x * sinAngle ;
	rotateMatrix[2][0]		= unitAxis.x * unitAxis.z * vsAngle - unitAxis.y * sinAngle;
	rotateMatrix[2][1]		= unitAxis.y * unitAxis.z * vsAngle + unitAxis.x * sinAngle;
	rotateMatrix[2][2]		= unitAxis.z * unitAxis.z * vsAngle + cosAngle;


	vector<Atom>			rotableAtomAll =
			dfsMol( rotableMol.get_atomVec(), rotableMol.get_bondVec() , fixedAtomIndex , rotableAtomIndex );

	//delete fixedAtomIndex and rotableAtomIndex from rotableAtomVec
	vector<Atom>			rotableAtomVec;
	for( size_t i=0; i<rotableAtomAll.size(); i++ ){
		if( rotableAtomAll[i].get_index() != fixedAtomIndex && rotableAtomAll[i].get_index() != rotableAtomIndex ){
			rotableAtomVec.push_back( rotableAtomAll[i] );
		}
	}

	for( size_t i=0; i<rotableAtomVec.size(); i++ ){
		vector<float>		coord(3, 0);
		coord[0]  = rotableAtomVec[i].get_coord().x  -  fixedAtom.get_coord().x;
		coord[1]  = rotableAtomVec[i].get_coord().y  -  fixedAtom.get_coord().y;
		coord[2]  = rotableAtomVec[i].get_coord().z  -  fixedAtom.get_coord().z;

		vector<float>		newCoord(3, 0);
		newCoord = matrixMultiplyVector( rotateMatrix, coord );
		newCoord[0] += fixedAtom.get_coord().x;
		newCoord[1] += fixedAtom.get_coord().y;
		newCoord[2] += fixedAtom.get_coord().z;
		rotableAtomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );

		//update
        for( size_t j=0; j<rotableMol_atVec.size(); j++ ){
            if( rotableMol_atVec[j].get_index() == rotableAtomVec[i].get_index() ){
                Atom 		newAt = rotableMol_atVec[j];
				newAt.set_coord( newCoord[0], newCoord[1], newCoord[2] );
				mol.update_atomCoord( newAt );
			}
		}
	}

	//check if new molecular has atoms crash
	vector<Atom> fixedAtomVec;
	vector<Atom> allAtomVec=rotableMol.get_atomVec();
	for( size_t i=0; i<allAtomVec.size(); i++ ){
		Atom a = allAtomVec[i];
		if( find( rotableAtomAll.begin(), rotableAtomAll.end(), a )==rotableAtomAll.end() &&
			a.get_index()!=fixedAtomIndex && a.get_index()!= rotableAtomIndex ){
			fixedAtomVec.push_back(a);
			for( size_t j=0; j<rotableAtomVec.size(); j++ ){
				float dis = getCoordDis( a.get_coord(), rotableAtomVec[j].get_coord() );
				if( dis<1.65 ){
					output = rotableMol;
					return false;
				}
			}
		}
	}
	/*
	vector<Bond> bdVec = mol.get_bondVec();
	vector<Atom> atVec = mol.get_atomVec();
	for( size_t i=0; i<atVec.size(); i++ ){
		for( size_t j=i+1; j<atVec.size(); j++ ){
			Atom a1 = atVec[i];
			Atom a2 = atVec[j];
			float dis = getCoordDis( a1.get_coord(), a2.get_coord() );
			if( dis < 1.65 ){
				bool flag = false;
				for( size_t k=0; k<bdVec.size(); k++ ){
					Atom b1 = bdVec[k].get_firstAtom();
					Atom b2 = bdVec[k].get_secAtom();
					if( ( a1.get_index()==b1.get_index() && a2.get_index()==b2.get_index() ) ||
						( a1.get_index()==b2.get_index()&&a2.get_index()==b1.get_index() ) ){
						flag = true;
						break;
					}
				}
				if( !flag ){
					output = rotableMol;
					return flag;
				}
			}
		}
	}
	*/

	output = mol;
	return true;
}

Molecular
rotateAroundAtom( const Molecular& rotableMol,
		const size_t& centerAtomIndex,
		const float& roundXAngle,
		const float& roundYAngle,
		const float& roundZAngle ){

	Molecular					newMol = rotableMol ;
	vector<Atom>				atomVec = newMol.get_atomVec();
	Atom						at = searchAtom( atomVec, centerAtomIndex );
	vector<float>				center(3, 0);
	center[0] = at.get_coord().x ;
	center[1] = at.get_coord().y ;
	center[2] = at.get_coord().z ;

	float							cos_x = cos( roundXAngle*3.1415926 / 180.0 );
	float							sin_x = sin( roundXAngle*3.1415926 / 180.0 );
	vector< vector<float> >			xMatrix(3, vector<float>(3, 0));
	xMatrix[1][1] = cos_x;
	xMatrix[1][2] = -sin_x;
	xMatrix[2][1] = sin_x;
	xMatrix[2][2] = cos_x;

	float						    cos_y = cos( roundYAngle*3.1415926 / 180.0 );
	float							sin_y = sin( roundYAngle*3.1415926 / 180.0 );
	vector< vector<float> >			yMatrix(3, vector<float>(3, 0));
	yMatrix[0][0] = cos_y;
	yMatrix[0][2] = -sin_y;
	yMatrix[2][0] = sin_y;
	yMatrix[2][2] = cos_y;

	float							cos_z = cos( roundZAngle*3.1415926 / 180.0 );
	float							sin_z = sin( roundZAngle*3.1415926 / 180.0 );
	vector< vector<float> >			zMatrix(3, vector<float>(3, 0));
	zMatrix[0][0] = cos_z;
	zMatrix[0][1] = -sin_z;
	zMatrix[1][0] = sin_z;
	zMatrix[1][1] = cos_z;

	for( size_t i=0; i<atomVec.size(); i++ ){
		vector<float>			coord(3, 0);
		coord[0] = atomVec[i].get_coord().x - center[0] ;
		coord[1] = atomVec[i].get_coord().y - center[1] ;
		coord[2] = atomVec[i].get_coord().z - center[2] ;
		vector<float>			newCoordX(3, 0);
		vector<float>			newCoordY(3, 0);
		vector<float>			newCoordZ(3, 0);
		newCoordX = matrixMultiplyVector( xMatrix, coord );
		newCoordX[0] = coord[0];

		newCoordY = matrixMultiplyVector( yMatrix, newCoordX );
		newCoordY[1] = newCoordX[1];

		newCoordZ = matrixMultiplyVector( zMatrix, newCoordY );
		newCoordZ[2] = newCoordY[2];

		vector<float>			newCoord(3, 0);
		newCoord[0] = newCoordZ[0] + center[0];
		newCoord[1] = newCoordZ[1] + center[1];
		newCoord[2] = newCoordZ[2] + center[2];

		atomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );
	}
	//	newLig.set_atomVec( atomVec );
	for( size_t i=0; i<atomVec.size(); i++ ){
		newMol.update_atomCoord( atomVec[i] );
	}
	return				newMol ;
}

Molecular
rotateAroundCenter( const Molecular& rotableMol,
		const Coord& rotateCenter,
		const float& roundXAngle,
		const float& roundYAngle,
		const float& roundZAngle ){

	Molecular							newMol = rotableMol ;
	vector<Atom>				atomVec = newMol.get_atomVec();

	vector<float>				center(3, 0);
	center[0] = rotateCenter.x ;
	center[1] = rotateCenter.y ;
	center[2] = rotateCenter.z ;

	float											cos_x = cos( roundXAngle*3.1415926 / 180.0 );
	float											sin_x = sin( roundXAngle*3.1415926 / 180.0 );
	vector< vector<float> >			xMatrix(3, vector<float>(3, 0));
	xMatrix[1][1] = cos_x;
	xMatrix[1][2] = -sin_x;
	xMatrix[2][1] = sin_x;
	xMatrix[2][2] = cos_x;

	float											cos_y = cos( roundYAngle*3.1415926 / 180.0 );
	float											sin_y = sin( roundYAngle*3.1415926 / 180.0 );
	vector< vector<float> >			yMatrix(3, vector<float>(3, 0));
	yMatrix[0][0] = cos_y;
	yMatrix[0][2] = -sin_y;
	yMatrix[2][0] = sin_y;
	yMatrix[2][2] = cos_y;

	float											cos_z = cos( roundZAngle*3.1415926 / 180.0 );
	float											sin_z = sin( roundZAngle*3.1415926 / 180.0 );
	vector< vector<float> >			zMatrix(3, vector<float>(3, 0));
	zMatrix[0][0] = cos_z;
	zMatrix[0][1] = -sin_z;
	zMatrix[1][0] = sin_z;
	zMatrix[1][1] = cos_z;

	for( size_t i=0; i<atomVec.size(); i++ ){
		vector<float>			coord(3, 0);
		coord[0] = atomVec[i].get_coord().x - center[0] ;
		coord[1] = atomVec[i].get_coord().y - center[1] ;
		coord[2] = atomVec[i].get_coord().z - center[2] ;
		vector<float>			newCoordX(3, 0);
		vector<float>			newCoordY(3, 0);
		vector<float>			newCoordZ(3, 0);
		newCoordX = matrixMultiplyVector( xMatrix, coord );
		newCoordX[0] = coord[0];

		newCoordY = matrixMultiplyVector( yMatrix, newCoordX );
		newCoordY[1] = newCoordX[1];

		newCoordZ = matrixMultiplyVector( zMatrix, newCoordY );
		newCoordZ[2] = newCoordY[2];

		vector<float>			newCoord(3, 0);
		newCoord[0] = newCoordZ[0] + center[0];
		newCoord[1] = newCoordZ[1] + center[1];
		newCoord[2] = newCoordZ[2] + center[2];

		atomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );
	}
	//	newLig.set_atomVec( atomVec );
	for( size_t i=0; i<atomVec.size(); i++ ){
		newMol.update_atomCoord( atomVec[i] );
	}
	return				newMol ;
}

Molecular
rotateAroundCoord( const Molecular& rotableMol,
		const Coord& rotateCenter,
		const float& roundXAngle,
		const float& roundYAngle,
		const float& roundZAngle ){


	Molecular							newMol = rotableMol ;
	vector<Atom>				atomVec = newMol.get_atomVec();
	vector<float>				center(3, 0);
	center[0] = rotateCenter.x;
	center[1] = rotateCenter.y;
	center[2] = rotateCenter.z;

	float											cos_x = cos( roundXAngle*3.1415926 / 180.0 );
	float											sin_x = sin( roundXAngle*3.1415926 / 180.0 );
	vector< vector<float> >			xMatrix(3, vector<float>(3, 0));
	xMatrix[1][1] = cos_x;
	xMatrix[1][2] = -sin_x;
	xMatrix[2][1] = sin_x;
	xMatrix[2][2] = cos_x;

	float											cos_y = cos( roundYAngle*3.1415926 / 180.0 );
	float											sin_y = sin( roundYAngle*3.1415926 / 180.0 );
	vector< vector<float> >			yMatrix(3, vector<float>(3, 0));
	yMatrix[0][0] = cos_y;
	yMatrix[0][2] = -sin_y;
	yMatrix[2][0] = sin_y;
	yMatrix[2][2] = cos_y;

	float											cos_z = cos( roundZAngle*3.1415926 / 180.0 );
	float											sin_z = sin( roundZAngle*3.1415926 / 180.0 );
	vector< vector<float> >			zMatrix(3, vector<float>(3, 0));
	zMatrix[0][0] = cos_z;
	zMatrix[0][1] = -sin_z;
	zMatrix[1][0] = sin_z;
	zMatrix[1][1] = cos_z;

	for( size_t i=0; i<atomVec.size(); i++ ){
		vector<float>			coord(3, 0);
		coord[0] = atomVec[i].get_coord().x - center[0] ;
		coord[1] = atomVec[i].get_coord().y - center[1] ;
		coord[2] = atomVec[i].get_coord().z - center[2] ;
		vector<float>			newCoordX(3, 0);
		vector<float>			newCoordY(3, 0);
		vector<float>			newCoordZ(3, 0);
		newCoordX = matrixMultiplyVector( xMatrix, coord );
		newCoordX[0] = coord[0];

		newCoordY = matrixMultiplyVector( yMatrix, newCoordX );
		newCoordY[1] = newCoordX[1];

		newCoordZ = matrixMultiplyVector( zMatrix, newCoordY );
		newCoordZ[2] = newCoordY[2];

		vector<float>			newCoord(3, 0);
		newCoord[0] = newCoordZ[0] + center[0];
		newCoord[1] = newCoordZ[1] + center[1];
		newCoord[2] = newCoordZ[2] + center[2];

		atomVec[i].set_coord( newCoord[0], newCoord[1], newCoord[2] );
	}

	for( size_t i=0; i<atomVec.size(); i++ ){
		newMol.update_atomCoord( atomVec[i] );
	}
	return				newMol ;

}

Molecular
translateMol( const Molecular& mol,
		const Coord& direction,
		const float& distance ){

	vector<Atom>    atomVec = mol.get_atomVec();

	Molecular       newMol;
	newMol = mol ;

	float length = sqrt( direction.x*direction.x +
			direction.y*direction.y +
			direction.z*direction.z );
	if( length == 0 ){
		cerr<<" translate ligand by zero direction"<<endl;
		return newMol;
	}
	Coord  unitDirection ;
	unitDirection.x = direction.x / length;
	unitDirection.y = direction.y / length;
	unitDirection.z = direction.z / length;

	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord newCoord;

		newCoord.x 	= atomVec[i].get_coord().x + unitDirection.x * distance ;
		newCoord.y 	= atomVec[i].get_coord().y + unitDirection.y * distance ;
		newCoord.z 	= atomVec[i].get_coord().z + unitDirection.z * distance ;
		atomVec[i].set_coord( newCoord.x, newCoord.y, newCoord.z );
	}

	for( size_t j=0; j<atomVec.size(); j++ ){
		newMol.update_atomCoord( atomVec[j] );
	}
	newMol.update_center();

	return newMol ;
}

vector<Atom>
getBindingSiteAtoms( 	const Molecular& mol,
		const Coord& center,
		const float& diameter ){
	vector<Atom>			atVec = mol.get_atomVec();
	vector<Atom>			bindsiteAtVec;
	for( size_t i=0; i<atVec.size(); i++ ){
		float						dis = getCoordDis( atVec[i].get_coord(), center );
		if( dis < diameter ){
			bindsiteAtVec.push_back( atVec[i] );
		}
	}
	return bindsiteAtVec;
}

vector<Atom>
getSubAtomVec( const Molecular& mol,  const vector<string>& atNameVec,
		const vector<size_t>& atIndexVec, const vector<size_t>& atResidNameVec ){
	vector<Atom>		atVec;
	vector<Atom>		molAtomVec = mol.get_atomVec();
	vector<string>		nameVec = atNameVec;
	vector<size_t>		indexVec = atIndexVec;
	vector<size_t>		resNameVec = atResidNameVec;
	for( size_t i=0; i<molAtomVec.size(); i++ ){
		for( size_t j=0; j<nameVec.size(); j++ ){
			if( 	molAtomVec[i].get_index() == indexVec[j] &&
					molAtomVec[i].get_name() == nameVec[j] &&
					molAtomVec[i].get_residueIndex() == resNameVec[j] ){
				atVec.push_back( molAtomVec[i] );
				nameVec.erase( nameVec.begin() + j  );
				indexVec.erase( indexVec.begin() + j  );
				resNameVec.erase( resNameVec.begin() + j );
			}
		}
		if( nameVec.empty() ){
			break;
		}
	}
	return						atVec;
}

vector<Bond>
getRotableBonds( const Molecular& rotableMol ){
	vector<Bond>						subBondVec;
	vector<Bond>						bdVec = rotableMol.get_bondVec();
	vector< vector<size_t> >			ringAtIndexVec = rotableMol.get_ringAtomIndex();

	for( size_t i=0; i<bdVec.size(); i++ ){
		if( bdVec[i].get_type() == "1" && bdVec[i].get_firstAtom().get_type() != "H" &&
				bdVec[i].get_secAtom().get_type() != "H" ){

			int count = 0;
			for( size_t m=0; m<ringAtIndexVec.size(); m++ ){
				count = 0;
				for( size_t n=0; n<ringAtIndexVec[m].size(); n++ ){
					if( bdVec[i].get_firstAtom().get_index() == ringAtIndexVec[m][n] ){
						count ++;
					}
					if( bdVec[i].get_secAtom().get_index() == ringAtIndexVec[m][n] ){
						count ++;
					}
				}
				if( count == 2 ){
					break;
				}
			}
			if( count != 2 ){
				subBondVec.push_back( bdVec[i] );
			}
		}
	}
	return											subBondVec ;
}


vector<pair<Atom, Atom> >
getBondAtomPairs( const Molecular& mol, const vector<pair<size_t, string> >& atVec ){
	vector<Bond>								molBondVec = mol.get_bondVec();
	vector<Atom>							molAtomVec = mol.get_atomVec();
	vector<pair<Atom, Atom> >		subBondAtomVec;
	for( size_t i=0; i<molBondVec.size(); i++ ){
		string										firstAtName = molBondVec[i].get_firstAtomName();
		size_t										firstAtIndex = molBondVec[i].get_firstAtomIndex();
		size_t										firstResidIndex = molBondVec[i].get_firstAtomResidIndex();

		string										secAtName = molBondVec[i].get_secAtomName();
		size_t										secAtIndex = molBondVec[i].get_secAtomIndex();
		size_t										secResidIndex = molBondVec[i].get_secAtomResidIndex();

		size_t										count = 0;
		for( size_t i=0; i<atVec.size(); i++ ){
			if( atVec[i].first == firstResidIndex && atVec[i].second == firstAtName ){
				count++;
				break;
			}
		}
		for( size_t i=0; i<atVec.size(); i++ ){
			if( atVec[i].first == secResidIndex && atVec[i].second == secAtName ){
				count++;
				break;
			}
		}
		if( count == 2 ){
			Atom									at1 = searchAtom( molAtomVec, firstAtName, firstAtIndex, firstResidIndex );
			Atom									at2 = searchAtom( molAtomVec, secAtName, secAtIndex, secResidIndex );
//			pair< Atom, Atom >			bdAtom = std::make_pair<Atom, Atom>( at1, at2 );
			pair< Atom, Atom >			bdAtom = std::make_pair( at1, at2 );
			subBondAtomVec.push_back( bdAtom );
		}
	}
	return 					subBondAtomVec;
}

vector<pair<Atom, Atom> >
getBondAtomPairs( const Molecular& mol, const vector< Atom >& atVec ){
	vector<Bond>								molBondVec = mol.get_bondVec();
	vector<pair<Atom, Atom> >		subBondAtomVec;
	//	cout<<"bond size: "<<molBondVec.size()<<endl;
	for( size_t i=0; i<molBondVec.size(); i++ ){
		Atom			a1 = molBondVec[i].get_firstAtom();
		Atom			a2 = molBondVec[i].get_secAtom();
		size_t			count = 0;
		for( size_t j=0; j<atVec.size(); j++ ){
			if( atVec[j] == a1 ){
				count++;
			}
			if( atVec[j] == a2 ){
				count ++;
			}
		}
		if( count >1 ){
//			pair< Atom, Atom >			bdAtom = std::make_pair<Atom, Atom>( a1, a2 );
			pair< Atom, Atom >			bdAtom = std::make_pair( a1, a2 );
			subBondAtomVec.push_back( bdAtom );
		}
	}

	return				subBondAtomVec;
}

void
getConformationalIsomers( const Molecular& inputMol, const float& bondRotateDegree,
		   vector<Molecular>& newMol, int bondIndex ){

	vector<Bond> bdVec = inputMol.get_rotableBondVec();
	cout<<endl<<"bondIndex:"<<bondIndex<<endl;
	cout<<"bdVec:"<<bdVec.size()<<endl;
	cout<<"vec:                  "<<newMol.size()<<endl;

	float degree =0;
	while( degree < 360 ){
		Molecular Mol = rotateAroundBond(inputMol, bdVec[bondIndex], degree);
		newMol.push_back(Mol);
		cout<<"degree:"<<degree<<endl;

		if( bondIndex+1 < bdVec.size() ){
			getConformationalIsomers(Mol, bondRotateDegree, newMol, bondIndex+1 );
			cout<<"?"<<endl;
		}
		degree=degree+bondRotateDegree;
		cout<<"        degree:"<<degree<<endl;
	}
}
