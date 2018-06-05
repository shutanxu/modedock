/*
 * clusterHelix.cpp
 *
 *  Created on: Sep 27, 2014
 *      Author: stan
 */

#include<iostream>
#include<fstream>
#include<stdlib.h>

#include"../mol/molecular.h"
#include"../algorithm/clusterdata.h"
#include"../algorithm/overlay.h"

#include"clusterHelix.h"

using namespace std;

void
clusterHelix(){
	/*
    string dir("PdbTest");
//	string dir("/media/stan/5387-1850/1009/DSSP12/");

	string input = dir + string("fileName");
    ifstream iff( input.c_str() );

    string oneLine;
    vector<Molecular> allMol;
    while( getline(iff, oneLine) ){
    	string fileName = dir + oneLine;
    	vector<Molecular> molVec = readPdb( fileName );
    	allMol.push_back( molVec.front() );
    	cout<<""<<allMol.size()<<endl;
    }

    vector< vector<Atom> > newMolAtomVec;
    vector<Molecular> newMolVec;
//    for( size_t i=0; i<allMol.size(); i++ ){
//    	string overlayAtomName("CA");
//        Molecular newM = overlay( allMol[0], allMol[i], overlayAtomName );
//    	newMolVec.push_back( newM );
//    	vector<Atom> atVec = newM.get_atomVec( overlayAtomName );
//    	newMolAtomVec.push_back( atVec );
//    	cout<<"overlay: "<<i<<endl;
//    }
    cout<<"---------"<<endl;

    for(size_t i=0; i<newMolVec.size(); i++){
    	string file = dir + newMolVec[i].get_name() + string(".pdb");
    	newMolVec[i].print( file );
    }
    cout<<"!!!!!!!!!!!!!"<<endl;


    vector<Molecular>             tempMolVec;
    vector< vector<Atom> >        tempAtomVec;
    for( size_t k=0; k<newMolVec.size(); k++ ){
    	string name("CA");
    	vector<Atom>    atVec1 = newMolVec[k].get_atomVec(name);
        if( atVec1.size() == 5 ){
    		tempMolVec.push_back( newMolVec[k] );
    		tempAtomVec.push_back( newMolAtomVec[k] );
    	}else{
    		cout<<"warnning!!! delete molecular:"<<newMolVec[k].get_name()<<endl;
    	}
    }
    newMolVec = tempMolVec;
    newMolAtomVec = tempAtomVec;


    cout<<"done"<<endl;

    vector< vector<float> > disMatrix(newMolVec.size(),
    		vector<float>(newMolVec.size(), 0) );

    for( size_t i=0; i<newMolAtomVec.size(); i++ ){
    	disMatrix[i][i] = 0;
    	for( size_t j=i+1; j<newMolAtomVec.size(); j++ ){
    		float rmsd = atomVecRMSD( newMolAtomVec[i], newMolAtomVec[j]  );
//    		float rmsd = atomVecSquareRMSD( atVec1, atVec2 );
    		disMatrix[i][j] = disMatrix[j][i] = rmsd;
    	}
    	cout<<"rmsd i: "<<i<<endl;
    }

    float rmsdCriteria = 1.5;

    vector< vector<size_t> > clusterResult = iterative_means( disMatrix,
    		rmsdCriteria );
//    vector< vector<size_t> > clusterResult = completeLink( disMatrix, rmsdCriteria );

    string resultFile = dir + string("result.txt");
    ofstream off( resultFile.c_str() );
    for( size_t i=0; i<clusterResult.size(); i++ ){
    	size_t centroid = getCentroid( disMatrix, clusterResult[i] );
    	off<<"$cluster:";
    	off.width(6);
    	off<<right<<i;
    	off<<" num:";
    	off.width(10);
    	off<<right<<clusterResult[i].size();
    	off<<" repre:";
    	off<<right<<newMolVec[centroid].get_name()<<".pdb"<<endl;
    	for( size_t j=0; j<clusterResult[i].size(); j++ ){
    		off<<newMolVec[ clusterResult[i][j] ].get_name()<<".pdb "<<endl;
    	}
    	off<<endl;
    }
    off.close();
    */
}

void
clusterHelix2(){
	/*
//	string dir("/home/stan/Desktop/CCCluster/We/14/");
    string dir("/home/stan/workspace/aDock/PdbTest/");
	string input = dir + string("fileName");
    ifstream iff( input.c_str() );

    string oneLine;
    vector<Molecular> allMol;
    while( getline(iff, oneLine) ){
    	string fileName = dir + oneLine;
    	vector<Molecular> molVec = readPdb( fileName );
    	allMol.push_back( molVec.front() );
    	cout<<""<<allMol.size()<<endl;
    }

    vector<Molecular> newMolVec;
    newMolVec.push_back( allMol[0] );
    for( size_t i=1; i<allMol.size(); i++ ){
        string overlayAtomName("CA");
        Molecular newM = overlay( allMol[0], allMol[i], overlayAtomName );
//        Molecular newM = overlay( allMol[0], allMol[i] );
        newMolVec.push_back( newM );
        cout<<"overlay: "<<i<<endl;
    }

    for(size_t i=0; i<newMolVec.size(); i++){
//    	newMolVec[i].print( );
    	string file = dir + newMolVec[i].get_name() + string(".pdb");
    	newMolVec[i].print( file );
    }

    vector< vector<float> > disMatrix(newMolVec.size(), vector<float>(newMolVec.size(), 0) );

    for( size_t i=0; i<newMolVec.size(); i++ ){
    	string name("CA");
    	vector<Atom>    atVec1 = newMolVec[i].get_atomVec(name);
        if( atVec1.size() != 5 ){
    		newMolVec.erase( newMolVec.begin() + i );
    		i--;
    		cout<<"# warnning!!! delete molecular "<<newMolVec[i].get_name()<<endl;
    	}
    }
    cout<<"done"<<endl;

    clock_t t0 = clock();
    for( size_t i=0; i<newMolVec.size(); i++ ){

    	string name("CA");
    	vector<Atom>    atVec1 = newMolVec[i].get_atomVec(name);
    	disMatrix[i][i] = 0;
    	for( size_t j=i+1; j<newMolVec.size(); j++){

    		vector<Atom>    atVec2 = newMolVec[j].get_atomVec(name);
//    		cout<<newMolVec[i].get_name()<<endl;
//    		cout<<newMolVec[j].get_name()<<endl;
    		float rmsd = atomVecRMSD( atVec1, atVec2 );
//    		float rmsd = atomVecSquareRMSD( atVec1, atVec2 );
    		disMatrix[i][j] = disMatrix[j][i] = rmsd;
//    		cout.width(10);
//    		cout<<rmsd;
    	}
    	cout<<"i:"<<i<<" "<<newMolVec.size()<<endl;
    	cout<<"time:"<<float( clock() - t0 )/ CLOCKS_PER_SEC<<endl;
    }

    float rmsdCriteria = 1.5;

    vector< vector<size_t> > clusterResult = iterative_means( disMatrix, rmsdCriteria );
//    vector< vector<size_t> > clusterResult = completeLink( disMatrix, rmsdCriteria );

    for( size_t i=0; i<clusterResult.size(); i++ ){
    	size_t centroid = getCentroid( disMatrix, clusterResult[i] );
    	cout<<"$cluster:";
    	cout.width(6);
    	cout<<right<<i;
    	cout<<" num:";
    	cout.width(10);
    	cout<<right<<clusterResult[i].size();
    	cout<<" repre:";
    	cout<<right<<newMolVec[centroid].get_name()<<".pdb"<<endl;
    	for( size_t j=0; j<clusterResult[i].size(); j++ ){
    		cout<<newMolVec[ clusterResult[i][j] ].get_name()<<".pdb "<<endl;
    	}
    	cout<<endl;
    }
}

void clusterFragments(){
//    string fileName("test");
    string fileName("ClusterCoord");
    ifstream iif( fileName.c_str() );
    string oneLine;
    vector<string> strVec;
    vector<float> cVec;
    vector< vector<float> > coordVec;

    while( getline(iif, oneLine) ){
        strVec = tokenBySpace(oneLine);

        if( strVec.size() < 3 ){
            coordVec.push_back(cVec);
            cVec.clear();
        }else{
            vector<float> coord(3,0);
            float x = atof( strVec[0].c_str() );
            float y = atof( strVec[1].c_str() );
            float z = atof( strVec[2].c_str() );
            cVec.push_back(x);
            cVec.push_back(y);
            cVec.push_back(z);
        }
    }
    coordVec.push_back(cVec);
    cVec.clear();

    vector<vector<float> > distMatrix(coordVec.size(),
                                     vector<float>(coordVec.size(),0));
    for( size_t i=0; i<coordVec.size(); i++ ){
        for( size_t j=i; j<coordVec.size(); j++ ){
            distMatrix[j][i]=distMatrix[i][j] =
            overlayRMSD(coordVec[i], coordVec[j]);
        }
//        cout<<i<<" "<<coordVec.size()<<endl;
    }

//    for( size_t i=0; i<distMatrix.size(); i++ ){
//        for( size_t j=0; j<distMatrix.size(); j++ ){
//            cout.width(12);
//            cout<<left<<distMatrix[i][j];
//        }
//        cout<<endl;
//    }

    float dis=0.5;
    vector<vector<size_t> > clusterResult = iterative_means(distMatrix, dis);
    for( size_t i=0; i<clusterResult.size(); i++ ){
        cout<<"cluster:"<<i<<" size:"<<clusterResult[i].size()<<endl;
        for( size_t j=0; j<clusterResult[i].size(); j++ ){
            cout<<clusterResult[i][j]<<" ";
        }
        cout<<endl<<endl;
    }
    */
}
