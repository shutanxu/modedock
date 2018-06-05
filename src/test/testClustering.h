/*
 * testClustering.h
 *
 *  Created on: Jun 24, 2014
 *      Author: stan
 */

#ifndef TESTCLUSTERING_H_
#define TESTCLUSTERING_H_

#include"../algorithm/clusterdata.h"

vector<pair<float, float> >
readDataFile( const string& inputFile ){
	ifstream			iff( inputFile.c_str() );
	if( !iff ){
		cerr<<"--- error can not open "<<inputFile<<"---"<<endl;
		throw;
	}

	vector< pair<float, float> >	inputCoordVec;
	string						oneLine;
	vector<string>		strVec;
	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		if( strVec.size() > 1 ){
			float					a = atof(strVec[1].c_str() );
			float					b = atof(strVec[2].c_str() );
			pair<float, float>	p = make_pair<float, float>( a, b );
			inputCoordVec.push_back( p );
		}
	}
	cout<<"input:"<<inputCoordVec.size()<<endl;
	return inputCoordVec;
}

vector< vector<float> >
read3dDataFile( const string& inputFile ){
//	string				inputFile( "bioData/tetra.lrn" );
	ifstream			iff( inputFile.c_str() );
	if( !iff ){
		cerr<<"--- error can not open "<<inputFile<<"---"<<endl;
		throw;
	}

	vector< vector<float> >	inputCoordVec;
	string						oneLine;
	vector<string>		strVec;
	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
			float					a = atof(strVec[1].c_str() );
			float					b = atof(strVec[2].c_str() );
			float					c = atof(strVec[3].c_str() );
			vector<float>	p;
			p.push_back( a );
			p.push_back( b );
			p.push_back( c );
			inputCoordVec.push_back( p );
	}
	cout<<"input:"<<inputCoordVec.size()<<endl;

	for( size_t i=0; i<inputCoordVec.size(); i++ ){
		cout<<i<<" :";
		for( size_t j=0; j<inputCoordVec[i].size(); j++ ){
			cout<<inputCoordVec[i][j]<<" ";
		}
		cout<<endl;
	}

	return inputCoordVec;
}

vector< vector<float> >
readFCPSdataFile( const string& inputFile ){
//	string				inputFile( "bioData/tetra.lrn" );
	ifstream			iff( inputFile.c_str() );
	if( !iff ){
		cerr<<"--- error can not open "<<inputFile<<"---"<<endl;
		throw;
	}

	vector< vector<float> >	inputCoordVec;
	string						oneLine;
	vector<string>		strVec;
	while( getline( iff, oneLine ) ){
		strVec = tokenBySpace( oneLine );
		if( strVec.size() == 4 ){
			float					a = atof(strVec[1].c_str() );
			float					b = atof(strVec[2].c_str() );
			float					c = atof(strVec[3].c_str() );
			vector<float>	p;
			p.push_back( a );
			p.push_back( b );
			p.push_back( c );
			inputCoordVec.push_back( p );
		}
	}
	cout<<"input:"<<inputCoordVec.size()<<endl;

	for( size_t i=0; i<inputCoordVec.size(); i++ ){
		cout<<i<<" :";
		for( size_t j=0; j<inputCoordVec[i].size(); j++ ){
			cout<<inputCoordVec[i][j]<<" ";
		}
		cout<<endl;
	}

	return inputCoordVec;
}

vector< vector<float> >
getDistMatrix( const vector< vector<float> >& inputData ){
	vector< vector<float> >		disMatrix( inputData.size(), vector<float>( inputData.size(), 0 ) );
	for( size_t i=0; i<inputData.size(); i++ ){
		disMatrix[i][i] = 0;
		for( size_t j=i+1; j<inputData.size(); j++ ){
			float 	dist = 0;
			if( inputData[i].size() == 3 ){
				dist = sqrt( ( inputData[i][0] - inputData[j][0] )*(  inputData[i][0] - inputData[j][0] ) +
						( inputData[i][1] - inputData[j][1] )*(  inputData[i][1] - inputData[j][1] ) +
						( inputData[i][2] - inputData[j][2] )*(  inputData[i][2] - inputData[j][2] ) );
			}
			disMatrix[i][j] = disMatrix[j][i] =dist;
		}
	}
	printDisMatrix( disMatrix );
	return disMatrix;
}

vector< vector<float> >
getDistMatrix( const vector<pair<float, float> >& inputData ){
	vector< vector<float> >		disMatrix( inputData.size(), vector<float>( inputData.size(), 0 ) );
	for( size_t i=0; i<inputData.size(); i++ ){
		disMatrix[i][i] = 0;
		for( size_t j=i+1; j<inputData.size(); j++ ){
			disMatrix[i][j] = disMatrix[j][i] =
					sqrt( ( inputData[i].first - inputData[j].first )*(  inputData[i].first - inputData[j].first ) +
							( inputData[i].second - inputData[j].second )*( inputData[i].second - inputData[j].second ) );
		}
	}
	printDisMatrix( disMatrix );
	return disMatrix;
}

void
testCircleClustering( const float& dia ){
	string			dataFile( "bioData/mlbench_circle.dat" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	cout<<"disM:"<<disMatrix.size()<<endl;
	float											diameter = dia;
	vector< vector<size_t> >			clusterResult = circleClustering( disMatrix, diameter, false );
	printCluster( clusterResult );
}

void
testCompleteLinkClustering( const float& dia ){
	string			dataFile( "bioData/mlbench_circle.dat" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	cout<<"disM:"<<disMatrix.size()<<endl;
	float											diameter = dia;
	vector< vector<size_t> >			clusterResult = completeLink( disMatrix, diameter, false );
	printCluster( clusterResult );
}

void
testBiDiversiveClustering( const float& dia ){
	string			dataFile( "bioData/mlbench_circle.dat" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	cout<<"disM:"<<disMatrix.size()<<endl;
	float											diameter = dia;
	vector< vector<size_t> >			clusterResult = biDiversive( disMatrix, diameter, false, false );
	printCluster( clusterResult );
}

void
testGeometricClustering( const float& dia ){
	string			dataFile( "bioData/mlbench_circle.dat" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	cout<<"disM:"<<disMatrix.size()<<endl;
	float											diameter = dia;
	vector< vector<size_t> >			clusterResult = iterative( disMatrix, diameter, false );
	printCluster( clusterResult );
}

void
testClusteringCircleData2D(){
	string			dataFile( "bioData/mlbench/circle_2d.txt" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 3.0; diameter= diameter+0.1 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = Hubert( disMatrix, diameter );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){

		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}


void
testClusteringCircleData3D(){
	string			dataFile( "bioData/mlbench/circle_3d.txt" );
	vector< vector<float>  >			data = read3dDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 3.5; diameter= diameter+0.1 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){
		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testClusteringhypecube3D(){
	string			dataFile( "bioData/mlbench/hypecube.txt" );
	vector< vector<float>  >			data = read3dDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 3.5; diameter= diameter+0.1 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){
		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testClusteringSmile(){
	string			dataFile( "bioData/mlbench/smiley.txt" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 3.0; diameter= diameter+0.1 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){

		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testClusteringThreeNorm2D(){
	string			dataFile( "bioData/mlbench/threenorm_2d.txt" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 11.0; diameter= diameter+0.2 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){
		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testClusteringThreeNorm3D(){
	string			dataFile( "bioData/mlbench/threenorm_3d.txt" );
	vector< vector<float> >	data = read3dDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.1; diameter < 11.0; diameter= diameter+0.2 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){
		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testClusteringTetra(){
	string			dataFile( "bioData/tetra.lrn" );
	vector< vector<float> >			data = readFCPSdataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
	for( float diameter = 0.9; diameter < 5.5; diameter= diameter+0.1 ){
		vector<float>						result;
		result.push_back( diameter );
		clusterResult = circleClustering( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = completeLink( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = biDiversive( disMatrix, diameter, false, false );
		cout<<clusterResult.size()<<" ";
		result.push_back( clusterResult.size() );
		clusterResult = iterative( disMatrix, diameter, false );
		cout<<clusterResult.size()<<" "<<endl;;
		result.push_back( clusterResult.size() );
		allCluster.push_back( result );
	}
	cout<<"end"<<endl;
	for( size_t i=0; i<allCluster.size(); i++ ){

		for( size_t j=0; j<allCluster[i].size(); j++ ){
			cout<<allCluster[i][j]<<" ";
		}
		cout<<endl;
	}
}

void
testHubert(){
	string			dataFile( "test.txt" );
	vector< pair<float, float> >	data = readDataFile( dataFile );
	vector< vector<float> >			disMatrix = getDistMatrix( data );
	vector< vector<size_t> >			clusterResult;
	vector< vector<float> >			allCluster;
//	for( float diameter = 2.8; diameter < 3.0; diameter= diameter+0.1 ){
//		vector<float>						result;
//		result.push_back( diameter );
//		clusterResult = HubertAlgorithmIterative( disMatrix, diameter );
//	}

//	clusterResult = Hubert( disMatrix,  (float)1.2 );
//	for( float dia = 0.8; dia<1.0; dia=dia+0.01 ){
		clusterResult = Hubert( disMatrix,  (float)1.5 );
//	}
}
#endif /* TESTCLUSTERING_H_ */
