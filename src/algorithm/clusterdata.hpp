/*
 * clusterdata.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: stan
 */

#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<fstream>
#include<time.h>
#include<sstream>

using namespace std;

template< typename _Tp >
struct clusterPair{
	size_t cluster_i;
	size_t cluster_j;

	//the distance between cluster_i and cluster_j
	_Tp		distance;
};

template <class T>
T GetMax (T a, T b) {
  return (a>b?a:b);
}

template< typename  _Tp >
vector<_Tp>
tokenizerBySpace( const _Tp& str ){
	string buf ="";  // Have a buffer string
//	trimSpace(str);
	stringstream ss(str);  // Insert the string into a stream
	vector<_Tp> tokens;  // Create vector to hold our words
	while (ss >> buf){
		tokens.push_back(buf);
		//     cout<<"buf :"<<buf<<endl;
	}
	//   cout<<"Size: "<<tokens.size()<<endl;
	buf ="";
	buf.clear();
	return tokens;
}

template<typename _Tp>
vector< vector<_Tp> >
getDisMatrix( const vector< _Tp >& inputDis ){
	size_t structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
	}

	vector< vector<_Tp> > inputMatrix(structNum, vector<_Tp>( structNum ));

	// initialize input distance matrix
	inputMatrix[0][0] = 0;
	for( size_t i=1; i<structNum; i++ ){
		for( size_t j=0; j<i; j++ ){
			inputMatrix[i][j] = inputMatrix[j][i] = inputDis[ i*(i-1)/2 + j  ];
		}
		inputMatrix[i][i] = 0;
	}

	return inputMatrix;
}

template<typename _Tp>
bool vectorComparatorFirst( const vector<_Tp>& v1,
		const vector<_Tp>& v2){
	size_t d1 = v1.size();
	size_t d2 = v2.size();
	if( d1 == d2 ){
		return ( v1.front() < v2.front() ? true: false );
	}else{
		return (d1 > d2? true: false);
	}
}

template<typename _Tp>
bool  vectorComparator(const vector<_Tp>& v1,
		 const vector<_Tp>& v2){
	size_t d1 = v1.size();
	size_t d2 = v2.size();
	return (d1 < d2? true: false);
};

template< typename _Tp >
bool clusterPairComparator(const clusterPair<_Tp>& clu1,
		const clusterPair<_Tp>& clu2){
	return (clu1.distance) < (clu2.distance) ? true: false ;
}

template< typename _Tp >
struct clusterPairComparator2{
	bool 	operator()( const clusterPair<_Tp>& clu1, const clusterPair<_Tp>& clu2 ){
//		return	clu1.distance < clu2.distance;
	}
};

template<typename _Tp>
vector< vector<_Tp> >
sortCluster( vector< vector<_Tp> >& input ){
	if( input.empty() ){
		cerr<<"---no input data in sortCluster---"<<endl;
	}
	if(0){
		cout<<"1"<<endl;
	}
	//sort cluster innder index sequence
	for(size_t i=0; i<input.size(); i++ ){
		sort(input[i].begin(), input[i].end() );
	}
	if(0){
		cout<<"2"<<endl;
	}

	//erase the empty cluster
	for( size_t i=0; i<input.size(); i++ ){
		if( 0 == input[i].size() ){
			input.erase( input.begin() + i );
			i--;
		}
	}
	sort(input.begin(), input.end(), vectorComparatorFirst<_Tp> );
	if(0){
		cout<<"3"<<endl;
	}
	//remove NULL cluster
	while( 0 == input.front().size()  ){
		input.erase( input.begin() );
	}
	if(0){
		cout<<"4"<<endl;
	}
	//sort cluster according to their size
	if(0){
		bool flag = true;
		while(flag ){
			flag = false;
			for(size_t i=0; i<input.size() - 1; i++ ){
				vector<_Tp> temp;
				if( input[i].size() < input[i+1].size() ){
					temp = input[i];
					input[i] = input[i+1];
					input[i+1] = temp;
					flag = true;
				}
			}
		}
	}
	return input;
}

template<typename _Tp>
_Tp
getMinPair( const vector<vector<_Tp > >& distMatrix , vector<size_t> groupVec ){
	_Tp		minVal = 1e10;
	for( size_t i=0; i<groupVec.size(); i++ ){
		for( size_t j=i+1; j<groupVec.size(); j++ ){
			minVal = minVal < distMatrix[ groupVec[i] ][ groupVec[j] ] ? minVal : distMatrix[ groupVec[i] ][ groupVec[j] ];
		}
	}
	return minVal;
}

template<typename _Tp>
_Tp
getMaxPair( const vector<vector<_Tp > >& distMatrix , vector<size_t> groupVec ){
	_Tp		maxVal = -1e10;
	for( size_t i=0; i<groupVec.size(); i++ ){
		for( size_t j=i+1; j<groupVec.size(); j++ ){
			maxVal = maxVal > distMatrix[ groupVec[i] ][ groupVec[j] ] ? maxVal : distMatrix[ groupVec[i] ][ groupVec[j] ];
		}
	}
	return maxVal;
}

template<typename _Tp>
_Tp
getAveragePair( const vector<vector<_Tp > >& distMatrix , vector<size_t> groupVec ){
	_Tp		aveVal = 0;
	size_t	count = 0;
	for( size_t i=0; i<groupVec.size(); i++ ){
		for( size_t j=i+1; j<groupVec.size(); j++ ){
			aveVal  += distMatrix[ groupVec[i] ][ groupVec[j] ];
			count ++;
		}
	}
	if( count != 0 ){
		aveVal /= count;
	}
	return aveVal;
}

template<typename _Tp>
size_t
getCentroid(  const vector<vector<_Tp > >& distMatrix , vector<size_t> indexesInCluster ){
	vector<_Tp> disVec( indexesInCluster.size(), 0 );
	for( size_t i=0; i<indexesInCluster.size(); i++ ){
		_Tp dis = 0;
		for( size_t j=0; j<indexesInCluster.size(); j++ ){
			size_t a = indexesInCluster[i];
			size_t b = indexesInCluster[j];
			dis += distMatrix[a][b];
		}
		disVec[i] = dis;
	}
	size_t minIndex = min_element( disVec.begin(), disVec.end() ) - disVec.begin() ;
	return indexesInCluster[minIndex];
}

/**
 * @brief		sort the cluster according to the cluster number
 * @param	input The input cluster to be sorted
 */
template<typename _Tp>
vector<vector<_Tp> >
sortCluster2(vector<vector<_Tp > >& input ){

	//sort cluster innder index sequence
	for(size_t i=0; i<input.size(); i++ ){
		sort(input[i].begin(), input[i].end() );
	}

	sort(input.begin(), input.end(), vectorComparator<_Tp> );

	//remove NULL cluster
	while( 0 == input.front().size()  ){
		input.erase( input.begin() );
	}

	//sort cluster according to their size
	if(0){
		bool flag = true;
		while(flag ){
			flag = false;
			for(size_t i=0; i<input.size() - 1; i++ ){
				vector<_Tp> temp;
				if( input[i].size() < input[i+1].size() ){
					temp = input[i];
					input[i] = input[i+1];
					input[i+1] = temp;
					flag = true;
				}
			}
		}
	}
	return input;
}

template<typename _Tp>
vector< vector<_Tp> >
readDisMatrix( const string& fileName ){
	ifstream iif;
	iif.open( fileName.c_str() );

	string sl;
	vector< vector<_Tp> > disMatrix;
	vector< _Tp > line;
	vector< string > tokenVec;

	disMatrix.clear();
	while( iif && !iif.eof() ){
		getline( iif, sl );
		tokenVec = tokenizerBySpace( sl );
		if( ! tokenVec.empty() ){
			line.clear();
			for( size_t i=0; i<tokenVec.size(); i++ ){
				line.push_back( (_Tp)atof( tokenVec[i].c_str() ) );
			}
			disMatrix.push_back( line );
		}
	}
	return disMatrix;
}

template<typename _Tp>
void
printClusterAverageRMSD( const vector< vector< _Tp > >& disMatrix,
		const vector< vector< size_t > >& clusterVec  ){
	if( disMatrix.empty() ){
		cerr<<"---no input data in printClusterAverageRMSD---"<<endl;
		throw ;
	}
	if( clusterVec.empty() ){
		cerr<<"---no input data in printClusterAverageRMSD---"<<endl;
		throw ;
	}

	cout<<"---cluster average RMSD---"<<endl;

	for( size_t clusterIndex = 0; clusterIndex < clusterVec.size(); clusterIndex++ ){
		_Tp averageDis = 0;
		size_t pairwiseNum = 0;
		pairwiseNum = ( clusterVec[ clusterIndex ].size() ) * ( clusterVec[clusterIndex].size() - 1 ) / 2;

		for( size_t structA = 0; structA < clusterVec[ clusterIndex ].size(); structA++ ){
			for( size_t structB = 0; structB < structA; structB++ ){
				averageDis += disMatrix[ clusterVec[ clusterIndex ][ structA ] ][ clusterVec[ clusterIndex ][ structB ] ];
				if(0){
					cout<<" "<<disMatrix[ clusterVec[ clusterIndex ][ structA ] ][ clusterVec[ clusterIndex ][ structB ] ]<<endl;
				}
			}
		}

		if( pairwiseNum != 0 ){
			averageDis /= pairwiseNum;
		}else{
			averageDis = 0;
		}

		cout<<"cluster: ";
		cout.width(5);
		cout<<left<<clusterIndex<<" average RMSD:"<<averageDis<<endl;
	}
}

template<typename _Tp>
void
printClusterCenterIndex( const vector< vector< _Tp > >& disMatrix,
		const vector< vector< size_t > >& clusterVec  ){
	if( disMatrix.empty() ){
		cerr<<"---no input data in printClusterCenterIndex---"<<endl;
		throw ;
	}
	if( clusterVec.empty() ){
		cerr<<"---no input data in printClusterCenterIndex---"<<endl;
		throw ;
	}

	cout<<"---cluster center data---"<<endl;

	for( size_t clusterIndex = 0; clusterIndex < clusterVec.size(); clusterIndex++ ){
		_Tp  maxDis = 1e10;
		size_t centerIndex;

		for( size_t i=0; i<clusterVec[ clusterIndex ].size(); i++ ){
			_Tp sumDis = 0;
			for( size_t j=0; j<clusterVec[ clusterIndex ].size(); j++ ){
				sumDis += disMatrix[ clusterVec[ clusterIndex ][ i ] ][ clusterVec[clusterIndex][ j ] ];
			}
			if( sumDis < maxDis && sumDis >= 0 ){
				maxDis = sumDis;
				centerIndex = clusterVec[ clusterIndex ][ i ];
			}
		}

		cout<<"cluster: ";
		cout.width(5);
		cout<<clusterIndex<<"  center data:"<<centerIndex<<endl;
	}
}


template<typename _Tp>
void
printDisMatrixAB(  const vector< vector< _Tp > >& disMatrix, string  fileName =  "matrix.txt"  ){
	ofstream oof;
	oof.open( fileName.c_str() );
	for( size_t i=0; i<disMatrix.size(); i++ ){
		for( size_t j=0; j<disMatrix[i].size(); j++ ){
			oof<<i<<":"<<j<<"->";
			oof.width(13);
			oof<<left<<disMatrix[i][j];
		}
		oof<<endl;
	}
	oof.close();
}

template<typename _Tp>
void
printDisMatrix( const vector< vector<_Tp> >& disMatrix ){
	ofstream oof;
	oof.open( "disMatrix.txt" );
	for( size_t i=0; i<disMatrix.size(); i++ ){
		oof.width(5);
		oof<<left<<i;
		for( size_t j=0; j<disMatrix[i].size(); j++ ){
			oof.width(3);
			oof<<right<<j<<":";
			oof.width(13);
			oof<<left<<disMatrix[i][j];
		}
		oof<<endl;
	}
	oof.close();
}

/**
 * @brief cout cluster result
 * @param resultVec
 * @return
 */
template<typename _Tp>
void
printCluster( const vector< vector< _Tp > >& resultVec ){
	size_t count = 0;
	if(1){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"# cluster"<<i<<": "<<endl;
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
				count++;
			}
			cout<<endl;
		}
	}
	cout<<"count:"<<count<<endl;
}

/**
 * @brief	print out cluster result
 * @param resultVec
 * @param disMatrix
 */
template<typename _Tp>
void
printCluster( const vector< vector<size_t> >& resultVec, const vector< vector<_Tp> >& disMatrix ){
	size_t count = 0;
	if(1){
		for( size_t i=0; i<resultVec.size(); i++ ){
			//compute largest distance in each cluster and print it out
			_Tp largestDis = largestInnerDis( disMatrix, resultVec[i] );
			cout<<"# cluster"<<i<<" num:"<<resultVec[i].size()<<" largest distance:"<<largestDis<<endl;
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
				count++;
			}
			cout<<endl;
		}
	}
	cout<<"count:"<<count<<endl;
}

/**
 * @brief		print out the cluster result
 * @param	clusterVec The cluster result
 */
template<typename _Tp>
void
printfCluster( const vector< vector<_Tp> > clusterVec, const string fileName ){
	cout<<"---cout clsuter result---"<<endl;
	if( clusterVec.empty() ){
		cerr<<" ---error no data after clustering---"<<endl;
		return;
	}

	ofstream of;
	of.open(  fileName.c_str()  );
	if( !of ){
		cerr<<"cant open printCluster.txt"<<endl;
		return;
	}

	//number of all ligands
	size_t		ligandNum = 0;
	for( size_t i=0; i<clusterVec.size(); i++ ){
		ligandNum += clusterVec[i].size();
	}
	of<<"# num:"<<ligandNum<<endl;

	for( size_t i=0; i<clusterVec.size(); i++ ){
		if( 0 != clusterVec[i].size() ){
			of<<"#  cluster:"<<i<<"  num:"<<clusterVec[i].size()<<endl;
			cout<<"#  cluster:"<<i<<"  num:"<<clusterVec[i].size()<<endl;
			for( size_t j=0; j<clusterVec[i].size(); j++ ){
				of<<clusterVec[i][j]<<"  ";
				cout<<clusterVec[i][j]<<"  ";
			}
			of<<endl;
			cout<<endl;
		}
	}
	of.close();
}

/**
 * @brief	  compare two vector, if all element are the same, then return true, else return false
 * @param vec1
 * @param vec2
 * @return
 */
template< typename _Tp >
bool
compareTwoVec( const vector<_Tp>& vec1, const vector<_Tp>& vec2 ){
	if( vec1.empty() || vec2.empty() ){
		return false;
	}
	if( vec1.size() != vec2.size()  ){
		return false;
	}

	vector< _Tp > v1( vec1 );
	vector< _Tp > v2( vec2 );
	sort( v1.begin(), v1.end() );
	sort( v2.begin(), v2.end() );

	bool flag = true;
	for( size_t i=0; i<v1.size(); i++ ){
		if( v1[ i ] != v2[ i ] ){
			flag = false;
			break;
		}
	}
	return flag;
}

/**
 * @brief		check whether a cluster if mergeable to other clusters
 * @param input
 * @return
 */
template<typename _Tp>
bool
checkMergeable( const vector< vector< _Tp > >& input ){
	bool flag = true;
	for( size_t i=0; i<input.size(); i++ ){
		if( input[i].empty() ){
			flag = false;
			break;
		}
	}
	return flag;
}

/**
 * @brief 		get the mergeable cluster index that have the smallest contradict with other mergeable clusters
 * @param splitVec
 * @param limiteDis
 * @return
 */
template< typename _Tp >
size_t
getMostLikelySplit( const vector< vector< vector<size_t> > >& splitVec, const _Tp& limiteDis ){
	if( splitVec.empty() ){
		cerr<<"---error no input data for getMostlikelySplit---"<<endl;
	}

	//the below is similar as $getConnectedMergable$
	//cluster connection matrix
	vector< vector<bool> > connectedVec( splitVec.size(), vector<bool>( splitVec.size(), false ) );

	for( size_t i=0; i<splitVec.size(); i++ ){
		for( size_t j=0; j<splitVec[i].size(); j++ ){
			//note the count indicate the number of clusters that have data share same clusters
			//ie cluster i has data j could be merged into s1 or s2, and cluster m has data n also
			//could be merged into s1 or s2, then count equals 2, note that
			size_t count = 0;
			for( size_t m =0; m<splitVec.size(); m++ ){
				if( m != i && checkMergeable( splitVec[i] ) && checkMergeable( splitVec[m] ) ){
					for( size_t n=0; n<splitVec[m].size(); n++ ){
						if( compareTwoVec( splitVec[i][j], splitVec[m][n] ) ){
							count++;
							if(0){
								cout<<"i:"<<i<<" j:"<<j<<" m:"<<m<<" n:"<<n<<" count:"<<count<<endl;
							}
						}
					}
				}
			}

			//note only when $count$ is larger than inputVec[i][j].size() could we say that cluster i and m combact
			//when merging.
			if( count >= splitVec[i][j].size() ){
				for( size_t m=0; m<splitVec.size(); m++ ){
					if( m != i ){
						for( size_t n=0; n<splitVec[m].size(); n++ ){
							if( compareTwoVec( splitVec[i][j], splitVec[m][n] ) ){
								connectedVec[i][m] = true;
							}
						}
					}
				}
			}
		}
	}

	//$countVec[i]$ indicate the number of contradict mergable clusters that one mergable cluster $i$ has
	vector<size_t> countVec( splitVec.size(), 0 );
	for( size_t i=0; i<connectedVec.size(); i++ ){
		size_t count = 0;
		for( size_t j=0; j<connectedVec[i].size(); j++ ){
			if( connectedVec[i][j] ){
				count++;
			}
		}
		countVec[i] = count;
	}

	//get the mergable cluster index that have the smallest contradict with other mergable clusters
	size_t min = 10e10;
	size_t minIndex = 0;
	for( size_t i=0; i<countVec.size(); i++ ){
		if( min > countVec[i] && checkMergeable( splitVec[i] ) ){
			if(0){
				cout<<"split i:"<<i<<endl;
			}
			min = countVec[i];
			minIndex = i;
			if( 0 == min ){
				break;
			}
		}
	}

	if(0){
		cout<<"countVec"<<endl;
		for( size_t i=0; i<countVec.size(); i++ ){
			cout<<" "<<i<<":"<<countVec[i];
		}
		cout<<endl;
	}

	return minIndex;
}

/**
 * @brief		return A[i][j][k] means the jth data in ith cluster could be merged into kth cluster
 * @param clusterVec
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector< vector<size_t> > >
getScratch( const vector< vector<size_t> >& clusterVec, const vector< vector<_Tp> >& inputMatrix,
		const _Tp& limiteDis  ){
	if( clusterVec.empty() || inputMatrix.empty() ){
		cerr<<"can't check if clusters could be scratched"<<endl;
	}

	//initialize scratch
	vector< vector< vector<size_t> > > scratch;
	for( size_t i=0; i<clusterVec.size(); i++ ){
		vector< vector< size_t > > vec1;
		for( size_t j=0; j<clusterVec[i].size(); j++ ){
			vector<size_t> vec2;
			vec1.push_back( vec2 );
		}
		scratch.push_back( vec1 );
	}

	for( size_t i=0; i<clusterVec.size()-1; i++ ){
		for( size_t j=0; j<clusterVec[i].size(); j++ ){
			bool check;
			size_t index = 0;

			//compare the ith data in ith cluster with all other clusters to see if it could be merged
			for( size_t m=0; m<clusterVec.size(); m++ ){
				if( m != i ){
					bool check = true;
					for( size_t n=0; n<clusterVec[m].size(); n++ ){
						if( inputMatrix[ clusterVec[i][j] ][ clusterVec[m][n] ] > limiteDis ){
							check = false;
						}
					}
					if( check ){
						scratch[ i ][ j ].push_back( m );
					}
				}
			}
		}
	}

	if(0){
		cout<<"---scratch all---"<<endl;
		for( size_t i=0; i<scratch.size(); i++ ){
			cout<<"cluster"<<i<<": "<<endl;
			for( size_t j=0; j<scratch[i].size(); j++ ){
				cout<<" j";
				cout.width(4);
				cout<<left<<j<<"->";
				cout.width(5);
				cout<<left<<clusterVec[i][j]<<": ";
				for( size_t k=0; k<scratch[i][j].size(); k++ ){
					cout.width(5);
					cout<<scratch[i][j][k]<<" ";
				}
				cout<<endl;
			}
//			cout<<endl;
		}
	}
	//print out only the mergeable clusters
	if(0){
		cout<<"---scratch o---"<<endl;
		size_t count = 0;
		for( size_t i=0; i<scratch.size(); i++ ){
			bool flag = true;
			for( size_t j=0; j<scratch[i].size(); j++ ){
				if( scratch[i][j].empty() ){
					flag = false;
				}
			}
			if( flag ){
				cout<<count<<" cluster"<<i<<": "<<endl;
				count++;
				for( size_t j=0; j<scratch[i].size(); j++ ){
					cout<<" j";
					cout.width(4);
					cout<<left<<j<<"->";
					cout.width(5);
					cout<<left<<clusterVec[i][j]<<": ";
					for( size_t k=0; k<scratch[i][j].size(); k++ ){
						cout.width(5);
						cout<<scratch[i][j][k]<<" ";
					}
					cout<<endl;
				}
				//			cout<<endl;
			}
		}
	}

	return scratch;
}

/**
 * @brief     at each step merge the mergeable cluster that share least contradict with other mergeable clusters
 * @param clusterInput
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
splitCluster_least( const vector< vector< size_t > >& clusterInput,
		const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis){

	if( clusterInput.empty() || inputMatrix.empty()){
		cerr<<"can't check if clusters could be scratched"<<endl;
	}

	vector< vector< size_t > > clusterVec( clusterInput );
	size_t count = 0;
	bool split = true;
	vector< vector< vector<size_t> > > scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );

	size_t splitIndex = 0;
	//get the first mergeable cluster index, which share least contradict with other mergeable clusters
	splitIndex = getMostLikelySplit( scratchVec, limiteDis );

	if(0){
		cout<<"------"<<endl;
		cout<<"the first mergeable cluster index:"<<splitIndex<<endl;
	}
	split = false;
	for( size_t i=0; i<scratchVec.size(); i++ ){
		if(checkMergeable( scratchVec[i] ) ){
			split = true;
		}
	}

	while( split ){
		if(0){
			cout<<"splitIndex:"<<splitIndex<<endl;
		}

		for( size_t j=0; j<scratchVec[ splitIndex ].size(); j++ ){
			clusterVec[ scratchVec[ splitIndex ][j].front() ].push_back( clusterVec[ splitIndex ][ j ] );
		}

		if(0){
			cout<<"removed cluster"<<splitIndex<<" :";
			for( size_t j=0; j<scratchVec[ splitIndex ].size(); j++ ){
				cout<<" "<<clusterVec[ splitIndex ][j];
			}
			cout<<endl;
		}
		//delete the splited cluster
		clusterVec.erase( clusterVec.begin() + splitIndex );

		//recompute the scratch vector
		scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );

		splitIndex = getMostLikelySplit( scratchVec, limiteDis );

		split = false;
		for( size_t i=0; i<scratchVec.size(); i++ ){
			if(checkMergeable( scratchVec[i] ) ){
				split = true;
			}
		}
	}

	return clusterVec;
}

/**
 * @brief		split cluster resutlt according to the scratch vector obtained from $getkScratch$
 * @param scratchVec
 * @param clusterVec
 * @param inputMatrix
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
splitCluster_rude( const vector< vector< size_t > >& clusterInput,
		const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis){

	if( clusterInput.empty() || inputMatrix.empty()){
		cerr<<"can't check if clusters could be scratched"<<endl;
	}

	vector< vector< size_t > > clusterVec( clusterInput );
	size_t count = 0;
	bool split = true;
	vector< vector< vector<size_t> > > scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );

	while( split ){
		//$split$ indicate whether cluster i could be splited into other clusters
		split = true;
		size_t clusterIndex = 0;
		for( size_t i=0; i<scratchVec.size(); i++ ){
			split = true;
			for( size_t j=0; j<scratchVec[i].size(); j++ ){
				if( scratchVec[i][j].empty() ){
					split = false;
					break;
				}
			}
			if( split ){
				clusterIndex = i;
				break;
			}
		}
		if( split ){
			for( size_t j=0; j<scratchVec[ clusterIndex ].size(); j++ ){
				//scratchVec[i][j].front() contains the cluster num that data cluster[i][j] could be merged into
				clusterVec[ scratchVec[ clusterIndex ][j].front() ].push_back( clusterVec[ clusterIndex ][j] );
			}
			if(1){
				cout<<"removed cluster"<<clusterIndex<<" :";
				for( size_t j=0; j<scratchVec[ clusterIndex ].size(); j++ ){
					cout<<" "<<clusterVec[ clusterIndex ][j];
				}
				cout<<endl;
			}
			//delete the splited cluster
			clusterVec.erase( clusterVec.begin() + clusterIndex );

			//recompute the scratch vector
			scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );
		}
//		split = false;
	}

	return clusterVec;
}

/**
 * @brief similar as $splitCluster$, the difference is this function cares more about
 * 		  the exhasut search way
 * @param clusterInput
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
splitClusterExhaust( const vector< vector< size_t > >& clusterInput,
		const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis){

	if( clusterInput.empty() || inputMatrix.empty()){
		cerr<<"can't check if clusters could be scratched"<<endl;
	}

	vector< vector< size_t > > clusterVec( clusterInput );
	size_t count = 0;

	bool split = true;
	vector< vector< vector<size_t> > > scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );

	if(1){
		vector< vector< vector<size_t> > > temp;
		temp = getIsolateMergable( scratchVec, limiteDis );
	}
	while( split ){
		//$split$ indicate whether cluster i could be splited into other clusters
		split = true;
		size_t clusterIndex = 0;
		for( size_t i=0; i<scratchVec.size(); i++ ){
			split = true;
			for( size_t j=0; j<scratchVec[i].size(); j++ ){
				if( scratchVec[i][j].empty() ){
					split = false;
					break;
				}
			}
			if( split ){
				clusterIndex = i;
				break;
			}
		}
		if( split ){
			for( size_t j=0; j<scratchVec[ clusterIndex ].size(); j++ ){
				//scratchVec[i][j].front() contains the cluster num that data cluster[i][j] could be merged into
				clusterVec[ scratchVec[ clusterIndex ][j].front() ].push_back( clusterVec[ clusterIndex ][j] );
			}
			if(1){
				cout<<"removed cluster"<<clusterIndex<<" :";
				for( size_t j=0; j<scratchVec[ clusterIndex ].size(); j++ ){
					cout<<" "<<clusterVec[ clusterIndex ][j];
				}
				cout<<endl;
			}
			//delete the splited cluster
			clusterVec.erase( clusterVec.begin() + clusterIndex );

			//recompute the scratch vector
			scratchVec = getScratch( clusterVec, inputMatrix, limiteDis );
		}
//		split = false;
	}

	return clusterVec;
}

/**
 * @brief 		get the connected mergable cluster matrix
 * @param inputVec
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector<bool> >
getConnectedMergable( vector< vector< vector<size_t> > >& inputVec, _Tp& limiteDis ){
	if( inputVec.empty() ){
		cerr<<"---error no input for checkIsoMergable---"<<endl;
	}
	//cluster connection matrix
	vector< vector<bool> > connectedVec( inputVec.size(), vector<bool>( inputVec.size(), false ) );

	for( size_t i=0; i<inputVec.size(); i++ ){
		for( size_t j=0; j<inputVec[i].size(); j++ ){
			//note the count indicate the number of clusters that have data share same clusters
			//ie cluster i has data j could be merged into s1 or s2, and cluster m has data n also
			//could be merged into s1 or s2, then count equals 2, note that
			size_t count = 0;
			for( size_t m =0; m<inputVec.size(); m++ ){
				if( m != i ){
					for( size_t n=0; n<inputVec[m].size(); n++ ){
						if( compareTwoVec( inputVec[i][j], inputVec[m][n] ) ){
							count++;
							if(0){
								cout<<"i:"<<i<<" j:"<<j<<" m:"<<m<<" n:"<<n<<" count:"<<count<<endl;
							}
						}
					}
				}
			}

			//note only when $count$ is larger than inputVec[i][j].size() could we say that cluster i and m combact
			//when merging.
			if( count >= inputVec[i][j].size() ){
				for( size_t m=0; m<inputVec.size(); m++ ){
					if( m != i ){
						for( size_t n=0; n<inputVec[m].size(); n++ ){
							if( compareTwoVec( inputVec[i][j], inputVec[m][n] ) ){
								connectedVec[i][m] = true;
							}
						}
					}
				}
			}

		}
	}

	if(1){
		cout<<"---connectedVec---"<<endl;
		cout<<"i:";
		cout<<"   ";
		for( size_t i=0; i<connectedVec.size(); i++ ){
			cout.width(3);
			cout<<left<<i;
		}
		cout<<endl;
		vector<size_t> countVec( connectedVec.size(), 0 );
		for( size_t i=0; i<connectedVec.size(); i++ ){
			cout<<"i:";
			cout.width(3);
			cout<<left<<i;

			size_t count = 0;
			for( size_t j=0; j<connectedVec[i].size(); j++ ){
				cout.width(3);
				if( connectedVec[i][j] ){
					cout<<left<<connectedVec[i][j];
				}else{
					cout<<left<<" ";
				}
				if( connectedVec[i][j] ){
					count++;
				}
			}
			countVec[i] = count;

			cout<<endl;
		}

		cout<<"count";
		for( size_t i=0; i<connectedVec.size(); i++ ){
			cout.width(3);
			cout<<left<<countVec[i];
		}
		cout<<endl;
	}
	return connectedVec;
}

/**
 * @brief get those mergable clusters which share no connection with other mergable clusters
 * 		  ie, if both cluster 1 and 2 could be merged into cluster A, and either first after
 * 		  merging cluster 1 to A will affects clustering 2 merging into A, or after merging
 * 		  cluster 2 to A will affect cluster 1 merging into to A.
 * @param clusterInput
 * @return
 * @note the input $inputVec$ is probably changed during this function
 */
template<typename _Tp>
vector< vector< vector<size_t> > >
getIsolateMergable( vector< vector< vector<size_t> > >& inputVec, _Tp& limiteDis ){
	if( inputVec.empty() ){
		cerr<<"---error no input for checkIsoMergable---"<<endl;
	}

	//obtain those that could be merged by other clusters
	vector< vector< vector<size_t> > > clusterInput;
	for( size_t i=0; i<inputVec.size(); i++ ){
		bool flag = true;
		for( size_t j=0; j<inputVec[i].size(); j++ ){
			if( inputVec[i][j].empty() ){
				flag = false;
			}
		}
		if( flag ){
			clusterInput.push_back( inputVec[i] );
		}
	}

	if(1){
		vector< vector<bool> > s = getConnectedMergable( clusterInput, limiteDis );
	}

	if(0){
		cout<<"---clusterInput---"<<endl;
		for( size_t i=0; i<clusterInput.size(); i++ ){
			cout<<"cluster "<<i<<endl;
			for( size_t j=0; j<clusterInput[i].size(); j++ ){
				cout<<" j"<<j<<": ";
				for( size_t k=0; k<clusterInput[i][j].size(); k++ ){
					cout<<" "<<clusterInput[i][j][k];
				}
				cout<<endl;
			}
			cout<<endl;
		}
	}

	bool flag = true;
	vector< vector< vector< size_t > > > resultVec;
	resultVec.clear();
	//cluster index
	for( size_t i=0; i<clusterInput.size(); i++ ){
		if(0){
			cout<<"i:"<<i<<endl;
		}
		//check whether cluster i share common destination clusters with other clusters, false if do share
		flag = true;
		for( size_t j=0; j<clusterInput[i].size(); j++ ){
			if(0){
				cout<<" j:"<<j<<endl;
			}
			for( size_t k=0; k<clusterInput[i][j].size(); k++ ){

				if(0){
					cout<<"i:"<<i<<" j:"<<j<<" k:"<<k<<endl;
					cout<<" clusterInput[i][j][k]:"<<clusterInput[i][j][k]<<endl;
				}
				//get cluster that data in cluster[i][j] belongs to
				//compare it with other data in other clusters
				for( size_t m=0; m<clusterInput.size(); m++ ){
					if( m != i ){
						for( size_t n=0; n<clusterInput[m].size(); n++ ){
							if( find( clusterInput[m][n].begin(), clusterInput[m][n].end(), clusterInput[i][j][k] ) !=
									clusterInput[m][n].end() ){
								flag = false;
								if(0){
									cout<<"m:"<<m<<" n:"<<n<<" i:"<<i<<" j:"<<j<<" k:"<<k
											<<" clusterInput"<<clusterInput[i][j][k]<<endl;
								}
							}
						}
					}
				}
				if(0){
					cout<<" flag:"<<flag<<endl;
				}
			}
		}
		if( flag ){
			resultVec.push_back( clusterInput[i] );
//			clusterInput.erase( clusterInput.begin() + i );
			if(0){
				cout<<"---isolated mergable---"<<endl;
				cout<<"cluster "<<i<<endl;
				for( size_t j=0; j<clusterInput[i].size(); j++ ){
					cout<<" j"<<j;
					for( size_t k=0; k<clusterInput[i][j].size(); k++ ){
						cout<<"  "<<clusterInput[i][j][k];
					}
					cout<<endl;
				}
				cout<<endl;
			}
		}
	}
	if(1){
		cout<<"---check isolated mergable clusters---"<<endl;
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"i:"<<i<<endl;
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" j:"<<j<<" ";
				for( size_t k=0; k<resultVec[i][j].size(); k++ ){
					cout<<" "<<resultVec[i][j][k];
				}
				cout<<endl;
			}
			cout<<endl;
		}
	}
//	cout<<"---check ios---"<<endl;
	return resultVec;
}

/**
 * @brief check the largest distance for each clusters
 * @param inputMatrix
 * @param clusterVec
 * @param disCriteria
 * @return
 */
template<typename _Tp>
bool
checkLargestInnerDis( const vector< vector<_Tp> >& inputMatrix,
		const vector< vector< size_t > >& clusterVec,
		const _Tp disCriteria = 10e10 ){
	if(1){
		cout<<"---get largest inner distance of each clusters---"<<endl;
	}
	for( size_t i=0; i<clusterVec.size(); i++ ){
		cout<<"cluster"<<i<<" ";
		_Tp maxDis = 0;

		//empty cluster
		if( clusterVec[i].empty() ){
			maxDis = 0;
		}else{
			for( size_t j=0; j<clusterVec[i].size()-1; j++ ){
				for( size_t k=j+1; k<clusterVec[i].size(); k++ ){
					if( maxDis < inputMatrix[ clusterVec[i][j] ][ clusterVec[i][k] ] ){
						maxDis = inputMatrix[ clusterVec[i][j] ][ clusterVec[i][k] ];
					}
				}
			}
		}
		cout<<" maxInnerDis:"<<maxDis;
		if( maxDis > disCriteria ){
			cout<<" warnning disCriteria:"<<disCriteria;
		}
		cout<<endl;
	}
	return true;
}

/**
 * @biref		Compute the largest distance between any two structures in clusterI and clusterJ
 * @param	inputMatrix pairwise distance matrix
 * @param	clusterI First cluster
 * @param	clusterJ Second cluster
 */
template<typename _Tp>
_Tp
largestDis(const vector<vector<_Tp> >& inputMatrix, const vector<size_t> clusterI,
		const vector<size_t> clusterJ){
	_Tp maxDis = 0;

	for(size_t i=0; i<clusterI.size(); i++ ){
		for(size_t j=0; j<clusterJ.size(); j++ ){
			if(maxDis < inputMatrix[clusterI[i] ][clusterJ[j] ] ){
				maxDis = inputMatrix[clusterI[i] ][clusterJ[j]];
			}
		}
	}
	return maxDis;
}

template<typename _Tp>
_Tp
largestDis(const vector<vector<_Tp> >& inputMatrix){
	_Tp maxDis = 0;

	for(size_t i=0; i<inputMatrix.size(); i++ ){
		for(size_t j=0; j<inputMatrix[i].size(); j++ ){
			if(maxDis < inputMatrix[i ][j] ){
				maxDis = inputMatrix[i][j];
			}
		}
	}
	return maxDis;
}


/**
 * @brief	count the number of largest pairwise distance
 * @param inputMatrix The input pairwise distance matrix
 * @return number of count
 */
template<typename _Tp>
size_t
computeNumOfLargestDis(const vector<vector<_Tp> >& inputMatrix ){
	size_t num = 0;
	_Tp		dis = 0;
	_Tp		maxDis = 0;
	for(size_t i=0; i<inputMatrix.size(); i++){
		for(size_t j=0; j<inputMatrix[i].size(); j++ ){
			dis = inputMatrix[i][j];
			if(dis > maxDis ){
				maxDis = dis;
			}
		}
	}
	for(size_t i=0; i<inputMatrix.size(); i++){
		for(size_t j=0; j<inputMatrix[i].size(); j++ ){
			if(maxDis == inputMatrix[i][j] ){
				cout<<i<<":"<<j<<" ";
				num ++;
			}
		}
	}
	cout<<"maxDis:"<<maxDis<<endl;
	return num;
}

/**
 * @brief		get the number of common data in two clusters
 * @param cluster1
 * @param cluster2
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
size_t
getCommonNum( const vector<size_t>& cluster1, const vector<size_t>& cluster2,
		const vector< vector<_Tp> >&inputMatrix, const _Tp& limiteDis ){
	if( cluster1.empty() || cluster2.empty() || inputMatrix.empty() ){
		cerr<<"---warnning!!! cant check common data in different clusters---"<<endl;
	}

	//count1 indicate the number of data in cluster1 that could be merged into cluster2
	size_t count1 = 0;
	for( size_t i=0; i<cluster1.size(); i++ ){
		bool flag = true;
		for( size_t j=0; j<cluster2.size(); j++ ){
			if( inputMatrix[ cluster1[ i ] ][ cluster2[ j ] ] > limiteDis ){
				flag = false;
			}
		}

		if( flag ){
			if(1){
				cout<<" clusteri[i]"<<cluster1[i]<<" limiteDis:"<<limiteDis<<endl;
			}
			count1++;
		}
	}

	//$count2$ is similar as $count1$
	size_t count2= 0;
	for( size_t i=0; i<cluster2.size(); i++ ){
		bool flag = true;
		for( size_t j=0; j<cluster1.size(); j++ ){
			if( inputMatrix[ cluster2[i] ][ cluster1[j] ] > limiteDis ){
				flag = false;
			}
		}

		if( flag ){
			count2++;
		}
	}

	if( count1 ){
		cout<<"count1:"<<count1<<endl;
	}
	if( count2 ){
		cout<<"count2:"<<count2<<endl;
	}
	size_t count = count1+count2;

	return count;
}

/**
 * @brief    deep first search
 * @param connectedVec
 * @param visitedVec
 * @param index
 * @return
 */
template<typename _Tp>
vector< _Tp >
dfs( const vector< vector<bool> >& connectedVec, vector< bool >& visitedVec, const _Tp& index){
	vector< size_t > vec;
	vector< size_t > temp;

	//data $index$ as first data
	if( !visitedVec[ index ] ){
		vec.push_back( index );
		visitedVec[ index ] = true;
	}

	for( size_t i=0; i<connectedVec[ index ].size(); i++ ){
		if( connectedVec[index][i] && !visitedVec[i] ){
			visitedVec[i] = true;
			vec.push_back( i );
			temp = dfs( connectedVec, visitedVec, i );
			vec.insert( vec.end(), temp.begin(), temp.end() );
		}
	}
	return vec;
}

/**
 * @brief    merge clusters once they share a common data
 * @param inputClusterVec
 * @param inputMatrix
 * @param criteriaDis
 * @param ratio   the ratio of common data num to cluster size,
 * 							only merge those clusters that the common data ratio is higher than that
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
mergeCluster_rude( const vector< vector< size_t > >& inputClusterVec,
		const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis,
		float ratioCri = 0.001 ){
	if( inputClusterVec.empty() || inputMatrix.empty() ){
		cerr<<"---warnning no input data for final merge---"<<endl;
	}

	//contain the merged cluster group vector
	vector< vector<size_t> > groupVec;

	vector< vector<size_t> > resultVec;

	if(0){
		vector< vector< vector<size_t> > > s = getScratch( inputClusterVec, inputMatrix, limiteDis );
	}

	//connection matrix
	vector< vector< bool > > connectedVec( inputClusterVec.size(), vector< bool >( inputClusterVec.size(), 0 ) );

	if(0){
		cout<<"i inputClusterVec.size():"<< inputClusterVec.size()<<endl;
	}

	//get the mergeable matrix, each element indicate that two clusters share common data
	for( size_t i=0; i<inputClusterVec.size() ; i++ ){
		for( size_t j=0; j<inputClusterVec.size(); j++ ){
			if( i != j ){
				if( inputClusterVec[i].empty() || inputClusterVec[j].empty() ){
					cerr<<"---warrning!!! no data in cluster "<<i<<" or "<<j<<"---"<<endl;
				}else{
					float ratio = 0;
					if( inputClusterVec[i].empty() || inputClusterVec[j].empty() ){
						ratio = 0;
					}else{
						size_t count = getCommonNum( inputClusterVec[ i ], inputClusterVec[ j ], inputMatrix, limiteDis );
						ratio = ( float )count/( inputClusterVec[i].size() + inputClusterVec[j].size() );
						if(0){
							if( count != 0 ){
								cout<<" i:"<<i<<" j:"<<j<<"  count:"<<count<<" ratio:"<<ratio;
								cout<<" clusteri.size()"<<inputClusterVec[i].size()<<" clusterj.size():"<<inputClusterVec[j].size()<<endl;
							}
						}
					}
					if( ratio > ratioCri ){
						connectedVec[ i ][ j ] = connectedVec[j][i] = 1;
					}
				}
			}
		}
	}

	if(0){
		cout<<"---connectedVec---"<<endl;
		for( size_t i=0; i<connectedVec.size(); i++ ){
			cout<<"i:";
			cout.width(3);
			cout<<left<<i<<" ";
			for( size_t j=0; j<connectedVec[i].size(); j++ ){
				cout<<" ";
				cout<<j<<":";
				cout.width(3);
				cout<<left<<connectedVec[i][j];
			}
			cout<<endl;
		}
	}
	if(0){
		cout<<"---connectedVec 1---"<<endl;
		for( size_t i=0; i<connectedVec.size(); i++ ){
			cout<<"i:";
			cout.width(3);
			cout<<left<<i<<"->";
			for( size_t j=0; j<connectedVec[i].size(); j++ ){
				if( connectedVec[i][j] ){
					cout<<" ";
					cout<<j<<":";
					cout.width(3);
					cout<<left<<connectedVec[i][j];
				}
			}
			cout<<endl;
		}
	}

	vector<bool> visitedVec( inputClusterVec.size(), false );
	for( size_t i=0; i<connectedVec.size(); i++ ){
		if( !visitedVec[i] ){
			vector<size_t> vec;
			vec = dfs( connectedVec, visitedVec, i );
			groupVec.push_back( vec );
		}
	}

	if(0){
		size_t count = 0;
		cout<<"---groupVec 0---"<<endl;
		for( size_t i=0; i<groupVec.size(); i++ ){
			cout<<i<<endl;
			for( size_t j=0; j<groupVec[i].size(); j++ ){
				cout<<"  "<<groupVec[i][j];
				count++;
			}
			cout<<endl;
		}
		cout<<"groupVec count:"<<count<<endl;
	}

	for( size_t i=0; i<groupVec.size(); i++ ){
		vector< size_t > largeCluster;
		largeCluster.clear();
		for( size_t j=0; j<groupVec[i].size(); j++ ){
			largeCluster.insert( largeCluster.end(), inputClusterVec[ groupVec[i][j] ].begin(), inputClusterVec[ groupVec[i][j] ].end() );
			if(1){
				cout<<" "<<groupVec[i][j]<<" largeCluster.size():"<<largeCluster.size()<<endl;
			}
//			largeCluster.insert( largeCluster.end(), inputClusterVec[ groupVec[i][j] ] );
		}
		resultVec.push_back( largeCluster );
	}
	if(1){
		cout<<"-----"<<endl;
		sortCluster( resultVec );
		printCluster( resultVec );
	}
	return resultVec;
}

/**
 * @brief merge cluster, somehow like the hierarchical
 * @param clusterVec
 * @param inputMatrix
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
mergeCluster( const vector< vector< size_t > >& inputClusterVec,
		const vector< vector< _Tp > >& inputMatrix, const _Tp& criteriaDis){
	if(1){
		cout<<"---merge cluster---"<<endl;
	}

	vector< vector< size_t > > clusterVec( inputClusterVec );

	if( clusterVec.size() <= 1 ){
		cerr<<"no need to merge cluster, empty or only one cluster"<<endl;
	}

	//identify whether continue to merge
	bool merge = false;
	size_t clusterI, clusterJ;
	for( size_t i=0; i<clusterVec.size(); i++ ){
		if( !merge ){
			for( size_t j=i+1; j<clusterVec.size(); j++  ){
				if( !merge ){
					if( largestDis( inputMatrix, clusterVec[i], clusterVec[j] ) < criteriaDis &&
							largestDis( inputMatrix, clusterVec[i], clusterVec[j] ) != 0 ){
						merge = true;
						clusterI = i;
						clusterJ = j;
					}
				}
			}
		}
	}
	while( merge ){
		//merge cluster
		if( clusterI < clusterJ ){
			clusterVec[ clusterI ].insert( clusterVec[ clusterI ].end(),
					clusterVec[ clusterJ ].begin(), clusterVec[ clusterJ ].end() );
			clusterVec[ clusterJ ].clear();
		}else{
			clusterVec[ clusterJ ].insert( clusterVec[ clusterJ ].end(),
					clusterVec[ clusterI ].begin(), clusterVec[ clusterI ].end() );
			clusterVec[ clusterI ].clear();
		}
		merge = false;

		//find two clusters that could be merged
		for( size_t i=0; i<clusterVec.size() ; i++ ){
			if( !merge ){
				for( size_t j=i+1; j<clusterVec.size() ; j++  ){
					if( !merge ){
						if( largestDis( inputMatrix, clusterVec[i], clusterVec[j] ) < criteriaDis &&
								largestDis( inputMatrix, clusterVec[i], clusterVec[j] ) != 0 ){
							merge = true;
							clusterI = i;
							clusterJ = j;
						}
					}
				}
			}
		}
	}

	for( size_t i=0; i<clusterVec.size(); i++ ){
		if( clusterVec[i].empty() ){
			clusterVec.erase( clusterVec.begin() + i );
			i--;
		}
	}
	clusterVec = sortCluster( clusterVec );

	if(1){
		checkClusterDis( clusterVec, inputMatrix, criteriaDis );
	}

	return clusterVec;
}

/**
 * @brief check the cluster result, check whether there exist
 *        two clusters' distance samller than  a criteria, which means
 *        these two clusters could be merged
 * @param clusterVec
 * @param inputMatrix
 * @param disCriteria
 * @return
 */
template<typename _Tp>
bool
checkClusterDis(  const vector< vector<size_t> >& clusterVec,
		const vector< vector<_Tp> >& inputMatrix, const _Tp disCriteria = 0){
//	cout<<"---check cluster distance---"<<endl;
	_Tp dis;
	_Tp maxDis = 0;
	size_t index1, index2;
	for( size_t i=0; i<clusterVec.size(); i++ ){
		for( size_t j=0; j<clusterVec.size() ; j++ ){
			if( i != j ){
				if(0){
					cout<<"i:"<<i<<" j:"<<j<<endl;
				}
				dis = 0;
				maxDis = 0;
				//get the largest distance between cluster i and cluster j
				for( size_t it=0; it<clusterVec[i].size(); it++ ){
					for( size_t jt=0; jt<clusterVec[j].size(); jt++ ){
						if(0){
							cout<<"it:"<<it<<" jt:"<<jt<<" inputMatrix.size():"<<inputMatrix.size()<<endl;
						}
						dis = inputMatrix[ clusterVec[i][it] ][ clusterVec[j][jt] ];
						if( dis > maxDis ){
							maxDis = dis;
							index1 = i;
							index2 = j;
						}
					}
				}
				if( maxDis < disCriteria ){
					cout<<"warrning cluster:"<<i<<"-"<<j<<" maxDis:"<<maxDis<<" limiteDis:"<<disCriteria<<endl;
				}
			}
		}
	}
	return true;
}

/**
 * @brief 	check if the input is matric
 * @param inputMatrc	The pairwise distance matrix
 * @brief matric means the triangle property: two edges' overall length
 * is further than the third edge's length
 */
template<typename _Tp>
bool
checkMatric(const vector<vector<_Tp> >& inputMatric, bool debug = false ){
	if(0){
		cout<<"---checkMatric---"<<endl;
	}
	if(inputMatric.size() < 3 ){
		cerr<<"cant check matric, must be more than 3 points!!!"<<endl;
	}
	bool flag = true;
	for(size_t i=2; i<inputMatric.size(); i++){
		for(size_t j=1; j<inputMatric[i].size(); j++ ){
			for(size_t k=1; k<j; k++ ){
				if(inputMatric[i][j] + inputMatric[i][k] < inputMatric[j][k] ||
				   inputMatric[i][j] + inputMatric[j][k] < inputMatric[i][k] ||
				   inputMatric[i][k] + inputMatric[j][k] < inputMatric[i][j]){
					flag = false;
					if(debug){
						cerr<<"i:"<<i<<" j:"<<j<<" k:"<<k<<endl;
						cerr<<i<<":"<<j<<"->"<<inputMatric[i][j]<<" ";
						cerr<<i<<":"<<k<<"->"<<inputMatric[i][k]<<" ";
						cerr<<j<<":"<<k<<"->"<<inputMatric[j][k]<<endl;;
					}
				}
			}
		}
	}
	if( !flag ){
		cerr<<"---warnning check distance matric failed!!!---"<<endl;
	}
	return flag;
}

/**
 * @brief		K-means clustering algorithm
 * @param	inputDis		The input distance matrix deposited in one dimension,
 * @param	clusterNum	The cluster number
 * @param	cluResultVec	The clustering result
 * @return 	true if success false if failed
 * @note		The 'inputDis' contains the pairwise distance of all structures in one dimension vector,
 * 					where $ ( i * ( i + 1 ) / 2) + j $th ( i>j ) item corresponding to  the [i][j]th item  in two dimension vectors.
 *
 *The correspondance between the one dimension 'inputDis' and two dimension distance matrix is as below
 * 	|    0                 0                  0                  0                0     |
 * 	| inputDis[0]     0                  0                  0                0     |
 * 	| inputDis[1] inputDis[2]       0                  0                0     |
 * 	| inputDis[3] inputDis[4] inputDis[5]         0                0     |
 * 	| inputDis[6] inputDis[7] inputDis[8]  inputDis[9]        0     |
 */
template< typename _Tp, typename _Size >
bool
k_means3( const vector< _Tp > inputDis,
		const _Size clusterNum , vector< vector< _Size > >& clusterResult ){
	cout<<"---kmeans1---"<<endl;

	_Size structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
		return false;
	}

	vector< vector<_Tp> > inputMatrix(structNum, vector<_Tp>( structNum ));

	if(0){
		cout<<"---inputDis---"<<endl;
		for( size_t i=0; i<inputDis.size(); i++ ){
			cout<<i<<":"<<inputDis.at(i)<<endl;
		}
	}

	// initialize input distance matrix
	inputMatrix[0][0] = 0;
	for( size_t i=1; i<structNum; i++ ){
		for( size_t j=0; j<i; j++ ){
			inputMatrix[i][j] = inputMatrix[j][i] = inputDis[ i*(i-1)/2 + j  ];
		}
		inputMatrix[i][i] = 0;
	}

	if(0){
		ofstream of( "disMatrix2.txt" );
		if( !of ){
			cout<<"cant open disMatrix.txt!"<<endl;
		}else {
			for( size_t i=0; i<inputMatrix.size(); i++ ){
				of.width(5);
				of<<i<<" ";
				for( size_t j=0; j<inputMatrix[i].size(); j++ ){
					of.width(5);
					of<<right<<j<<":";
					of.width(15);
					of<<left<<inputMatrix[i][j];
				}
				of<<endl;
			}
		}
		of.close();
	}

	clusterResult.clear();
	clusterResult = k_means( inputMatrix, clusterNum );
}

/**
 * @brief		K-means clustering algorithm
 * @param	inputDis		The input distance matrix deposited in one dimension,
 * @param	clusterNum	The cluster number
 * @param	cluResultVec	The clustering result
 * @return 	true if success false if failed
 * @note		The 'inputDis' contains the pairwise distance of all structures in one dimension vector,
 * 					where $ ( i * ( i + 1 ) / 2) + j $th item' corresponding to  the [i][j]th item  in two dimension vectors.
 */
template< typename _Tp, typename _Size >
bool
k_means2( const vector< _Tp > inputDis,
		const _Size clusterNum , vector< vector< _Size > >& clusterResult, _Size maxStep = 1024){
	cout<<"---kmeans---"<<endl;

	//error
	if( inputDis.empty() ){
		cerr<<"---error in k_menas, no pairwise distance data---"<<endl;
		return false;
	}
	if( clusterNum < 0 ){
		cerr<<"---error in k-menas, please input correct cluster number for k-menas---"<<endl;
		return false;
	}

	_Size structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
		return false;
	}

	if(0){
		cout<<"inputDis.size():"<<inputDis.size()<<"  structNum:"<<structNum<<endl;
		cout<<"---inputDis:---"<<endl;
		for( size_t i=0; i<inputDis.size(); i++ ){
			cout<<i<<":"<<inputDis.at(i)<<endl;
		}
		cout<<endl;
	}

	//contains the old cluster center vector
	vector<_Tp>				oldCenterVec;

	// all center vector, including all history centers
	vector< vector<_Tp> >				allCenterVec;

	//contains the new cluster center vector
	vector<_Tp>				newCenterVec;
	_Size					centerIndex = 0;

	//boolean flag whether clustering center change
	bool 					changed = true;

	//initialize resultVec
	vector< vector< _Tp > > resultVec( clusterNum );

	oldCenterVec.clear();
	for( _Size i=0; i<clusterNum; i++ ){
//		unsigned seed = 1000;
//		srand( seed );
		centerIndex = (_Tp)(rand()%structNum);
		if(1){
			cout<<"i:"<<i<<"  centerIndex:"<<centerIndex<<endl;
		}
		oldCenterVec.push_back( centerIndex );
	}

	if(1){
		oldCenterVec.at(0) = 20;
		oldCenterVec.at(1) = 7;
		oldCenterVec.at(2) = 6;
		oldCenterVec.at(3) = 11;
	}

	allCenterVec.push_back( oldCenterVec );

	//compute group for each data
	for( size_t i=0; i<(size_t)(structNum); i++ ){

		//distance to center
		vector< _Tp > centerDisVec;

		for( size_t j=0; j<oldCenterVec.size(); j++ ){

			_Tp	value;
			// compute the distance between ith and jth center
			if( i < oldCenterVec[j] ){
				value =  inputDis[ i*(i+1)/2 + oldCenterVec[j]  ];
			}else{
				value =  inputDis[ oldCenterVec[j] * ( oldCenterVec[j] +1 )/2 + i ] ;
			}
			centerDisVec.push_back( value );

			if(0){
				cout<<"structNum:"<<structNum<<"   i:"<<i<<"  oldCenterVec:"
						<<oldCenterVec[j]<<"  centerDisVec.back():"<<centerDisVec.back()<<"  value:"<<value<<endl;
			}
		}

		//compute the closest center for ith data
		_Size  minClusterIndex = ( _Size )( min_element( centerDisVec.begin(), centerDisVec.end() )
																- centerDisVec.begin() );

		// add ith data to its closest center
		// add ith data to tis closest center
		if( find( resultVec[ minClusterIndex ].begin(), resultVec[ minClusterIndex ].end(), i ) ==
				resultVec[ minClusterIndex ].end()  ){
			resultVec[ minClusterIndex ].push_back( i );
		}
	}

	if(1){
		cout<<"---oldCenterVec---"<<endl;
		for( size_t i=0; i<oldCenterVec.size(); i++ ){
			cout<<"i:"<<oldCenterVec[i]<<"  ";
		}
		cout<<endl;
		cout<<"---resultVec---"<<endl;
		for( size_t i=0; i<resultVec.size(); i++ ){
			for( size_t j=0; j<resultVec[i].size(); j++ ){
			cout<<" "<<resultVec[i][j];
			}
			cout<<endl;
		}
	}

	//maximum step
	size_t 	step  =0;

	while( changed && step < (size_t)maxStep ){

		step ++;
		if(1){
			cout<<"step:"<<step<<endl;
		}

		changed = false;
		newCenterVec.clear();

		// compute cluster center
		for( size_t i=0; i<resultVec.size(); i++ ){
			_Tp 	dis = 0;
			_Tp	maxDis = 10e10;
			_Size	centerIndex = 0;
			if( 1 == resultVec[i].size() ){
				centerIndex = resultVec[i].front();
			}else{
				for( _Size j=1; j<resultVec[i].size(); j++ ){
					dis = 0;
					for( _Size k=0; k<resultVec[i].size(); k++ ){
						if( j<k ){
							dis += inputDis[ j * (j+1)/2 + k ];
						}else if( k<j ){
							dis += inputDis[ k * (k+1)/2 + j ];
						}
					}
					if( maxDis > dis ){
						maxDis = dis;
						centerIndex = resultVec[i][j];
					}
				}
			}
			newCenterVec.push_back( centerIndex );
		}

		if(1){
			cout<<"---newCenterVec---"<<endl;
			for( size_t i=0; i<newCenterVec.size();  i++){
				cout<<" "<<newCenterVec[i];
			}
			cout<<endl;
		}

		//check if the new center have shown before
		for( size_t i=0; i<allCenterVec.size(); i++ ){
			bool fla = true;
			for( size_t j=0; j<newCenterVec.size(); j++ ){
				if(0){
					cout<<j<<" allCenterVec[i][j]:"<<allCenterVec[i][j]<<"  newCenterVec[j]:"<<newCenterVec[j];
				}
				if( allCenterVec[i][j] != newCenterVec[j] ){
					fla = false;
				}
			}
			if( fla == true ){
				changed = true;
			}
		}

		if( changed ){
			//store center
			oldCenterVec = newCenterVec;
			allCenterVec.push_back( oldCenterVec );
			for( size_t i=0; i<resultVec.size(); i++ ){
				resultVec[i].clear();
			}

			if(0){
				cout<<"---inputDis--"<<endl;
				for( size_t k=0; k<inputDis.size(); k++ ){
					cout<<k<<":"<<inputDis[k]<<"  ";
				}
				cout<<endl;
			}

			//cluster data according to their distance with centers
			for( size_t i=0; i<structNum; i++ ){
				vector< _Tp >  centerDisVec;
				for( size_t j=0; j<oldCenterVec.size(); j++ ){
					if(0){
						cout<<"i:"<<i<<"  oldCenterVec.ati(j):"<<oldCenterVec[j]<<endl;
					}
					// compute the distance between ith data and jth center
					if( i > oldCenterVec[j] ){
						centerDisVec.push_back( inputDis[ i * ( i - 1 )/2 + oldCenterVec[j]  ] );
					}else if( i < oldCenterVec[ j ] ){
						centerDisVec.push_back( inputDis[ (oldCenterVec[j] - 1 ) * ( oldCenterVec[j] )/2 + i  ]   );
					}else if( i == oldCenterVec[ j ] ){
						centerDisVec.push_back( (_Tp)0 );
					}
				}

				//compute the closest center for ith data
				_Size	minClusterIndex = ( _Size )( min_element( centerDisVec.begin(), centerDisVec.end() ) -
																centerDisVec.begin() );

				if(0){
					cout<<"i:"<<i<<"  minClusterIndex:"<<minClusterIndex;
					for( size_t k=0; k<centerDisVec.size();  k++){
						cout<<"  "<<k<<":"<<centerDisVec.at(k)<<" ";
					}
					cout<<endl;
				}

				//add ith data to its neareast center
				if( find( resultVec[ minClusterIndex ].begin(), resultVec[ minClusterIndex ].end(), i ) ==
						resultVec[ minClusterIndex ].end() ){
					resultVec[ minClusterIndex ].push_back( i );
				}
			}
		}

		if(0){
			cout<<"---oldCenterVec---"<<endl;
			for( size_t i=0; i<oldCenterVec.size(); i++ ){
				cout<<"i:"<<oldCenterVec[i]<<"  ";
			}
			cout<<endl;
		}
		if(1){
			cout<<"---resultVec---"<<endl;
			for( size_t i=0; i<resultVec.size(); i++ ){
				for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" "<<resultVec[i][j];
				}
				cout<<endl;
			}
		}
		if(0){
			return false;
		}
	}

	resultVec = sortCluster(resultVec );

	if(0){
		cout<<"resultVec.size():"<<resultVec.size()<<endl;
	}
	if(1){
		cout<<"---k-menas result---"<<endl;
		for(size_t i=0; i<resultVec.size(); i++){
			cout<<i<<" cluster : ";
			for(size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" "<<resultVec[i][j];
			}
			cout<<endl;
		}
	}

	string file = "kmeans.txt";
	printfCluster( resultVec, file );
}

/**
 * @biref 	K-means clustering algorithm using a pairwise distance
 * 			matrix and initial cluster number for clustering
 * @param	inputMatrix	Pairwise distance matrix
 * @param 	clusterNum	Initial cluster number for K-means algorithm
 * @return	classified two dimension vector
 * @note	The inputMatrix must be real rmsd distance, not rmsd square
 * 			distance or anything else.
 */
template<typename _Tp, typename _Size>
vector<vector<_Size> >
k_means( const vector<vector<_Tp> >& inputMatrix, const _Size clusterNum ){
	cout<<endl;
	cout<<"---kmeans---"<<endl;
	clock_t t1, t2;
	t1 = clock();

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for kmeans clustering!!!---"<<endl;
	}
	if(clusterNum == 0 || clusterNum < 0 ){
		cerr<<"---must initialize cluster number for kmeans!!!---"<<endl;
	}

	//temp vector
	vector<vector<_Size> > 	resultVec(clusterNum );

	//contains the old cluster center vector
	vector<_Tp>				oldCenterVec;

	//contains the new cluster center vector
	vector<_Tp>				newCenterVec;

	_Size					centerIndex = 0;
	_Size					matrixSize = inputMatrix.size();

	//Boolean flag whether clustering center changes
	bool 					changed = true;

	//initialize center randomly, put the center index in centerVec
	oldCenterVec.clear();
	srand( (unsigned)time(0) );
	for(_Size i=0; i<clusterNum; i++){
		centerIndex = rand()%inputMatrix.size();
		oldCenterVec.push_back(centerIndex );

		//the first index of resultVec[i] contains the center data index
		resultVec[i].push_back(centerIndex );
		cout<<i<<" center:"<<centerIndex<<endl;
	}

	//compute group for each data
	for(_Size i=0; i< matrixSize; i++ ){
		vector<_Size > centerDisVec;
		for(_Size j=0; j<oldCenterVec.size(); j++ ){

			//compute the distance between ith data and jth center
			centerDisVec.push_back( inputMatrix[i][oldCenterVec[j] ] );
		}
		//compute the closet center for ith data
		_Size minClusterIndex = (_Size)(min_element(centerDisVec.begin(),
				centerDisVec.end()) - centerDisVec.begin() );

		//add ith data to its closest center
		//add ith data to its closest center
		if(find(resultVec[minClusterIndex ].begin(), resultVec[minClusterIndex ].end(),
				i) == resultVec[minClusterIndex ].end() ){
			resultVec[minClusterIndex ].push_back(i);
		}
	}

	if(0){
		cout<<"1"<<endl;
	}
	while(changed ){
		changed = false;
		newCenterVec.clear();

		if(0){
//			resultVec = sortCluster(resultVec );
			cout<<"changed:"<<changed<<endl;
			for(size_t i=0; i<resultVec.size(); i++){
				cout<<i<<" ";
				for(size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<" "<<resultVec[i][j];
				}
				cout<<endl;
			}
			cout<<"old cluster center:"<<endl;
			for(size_t i=0; i<oldCenterVec.size(); i++ ){
				cout<<i<<" "<<oldCenterVec[i]<<endl;
			}
		}

		//compute cluster center
		for(_Size i=0;  i<resultVec.size(); i++){
			//for each data inner cluster compute its distance with other data
			//the one have smallest distance is center, it is also called as 'medoid'
			_Tp		dis = 0;
			_Tp 	maxDis = 10e10;
			_Size	centeIndex = 0;
			for(_Size j=0; j<resultVec[i].size(); j++ ){
				dis = 0;
				for(_Size k=0; k<resultVec[i].size(); k++ ){
					dis += inputMatrix[ resultVec[i][j] ][ resultVec[i][k] ];
				}
				if(maxDis > dis ){
					maxDis = dis;
					centerIndex = resultVec[i][j];
				}
			}
			newCenterVec.push_back(centerIndex );
		}

		if(0){
			cout<<"new center!!!"<<endl;
			for(size_t i=0; i<newCenterVec.size(); i++ ){
				cout<<i<<" "<<newCenterVec[i]<<endl;
			}
		}

		//check whether the center has moved
		for(_Size i=0; i<oldCenterVec.size(); i++ ){
			if(oldCenterVec[i] != newCenterVec[i] ){
				changed = true;
			}
		}

		if( changed ){
			if(0){
				cout<<"changed:"<<changed<<endl;
			}
			oldCenterVec = newCenterVec;

			for( size_t i=0; i<resultVec.size(); i++ ){
				resultVec[i].clear();
			}

			//cluster data according to their distance to centers
			for(_Size i=0; i< matrixSize; i++ ){
				if(0){
					cout<<"i:"<<i<<endl;
				}
				vector<_Size > centerDisVec;
				for(_Size j=0; j<oldCenterVec.size(); j++ ){

					//compute the distance between ith data and jth center
					centerDisVec.push_back( inputMatrix[i][oldCenterVec[j] ] );
				}
				//compute the closet center for ith data
				_Size minClusterIndex = (_Size)(min_element(centerDisVec.begin(),
						centerDisVec.end()) - centerDisVec.begin() );

				//add ith data to its closest center
				if(find(resultVec[minClusterIndex ].begin(), resultVec[minClusterIndex ].end(),
						i) == resultVec[minClusterIndex ].end() ){
					resultVec[minClusterIndex ].push_back(i);
					if(0){
						cout<<i<<"->"<<minClusterIndex<<endl;
					}
				}
			}
		}
	}
	if(1){
		cout<<"k-means running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	resultVec = sortCluster(resultVec );

	if(0){
		cout<<"resultVec.size():"<<resultVec.size()<<endl;
	}
	if(0){
		cout<<"---k-menas result---"<<endl;
		for(size_t i=0; i<resultVec.size(); i++){
			cout<<i<<" cluster : ";
			for(size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" "<<resultVec[i][j];
			}
			cout<<endl;
		}
	}

	printCluster( resultVec, inputMatrix );
	checkClusterDis( resultVec, inputMatrix );
	checkLargestInnerDis( inputMatrix, resultVec );

	if(0){
		string file = "kmeans.txt";
		printfCluster( resultVec, file );
	}

	return resultVec;
}

/**
 *	@biref 			Extended $K$-means algorithm, the initial cluster center is not randomly selected.
 *	@inputMatrix	The pairwise distance matrix
 *	@clusterNum 	The initial cluster number
 *	@return 		The cluster result
 *	@see kmeans
 *
 *	given neighborhood definiton, the data with most neighbors are substracted
 */
template<typename _Tp>
vector<vector<size_t> >
density(vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	cout<<"---density---"<<endl;

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for density !!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//neighbor number of a data
	size_t 					neighborNum = 0;

	//neighborhoods for all data
	vector<vector<size_t> > allNeighborVec(matrixSize );

	//neighborhoods for a data
	vector<size_t >			neighborVec;

	//unclassified data index
	vector<size_t>			unclassifiedVec;

	//initial the unclassified data index
	for(size_t i=0; i<matrixSize; i++){
		unclassifiedVec.push_back(i);
	}

	while( unclassifiedVec.size() != 0 ){
		if(0){
			cout<<"matrixSize:"<<matrixSize<<endl;
			cout<<"unclassifiedVec.size():"<<unclassifiedVec.size()<<endl;
		}
		allNeighborVec.clear();
		neighborVec.clear();

		//find the data with the lagest amount of
		//neighbors and take them as a cluster
		for(size_t i=0; i<unclassifiedVec.size(); i++ ){
			neighborVec.clear();
			for(size_t j=0; j<unclassifiedVec.size(); j++ ){
				if( inputMatrix[ unclassifiedVec[i] ][ unclassifiedVec[j] ] < limiteDis ){
					neighborVec.push_back( unclassifiedVec[j] );
				}
			}
			allNeighborVec.push_back( neighborVec );
		}

		if(0){
			cout<<"allNeighborVec:"<<endl;
			cout<<"allNeighborVec.size():"<<allNeighborVec.size()<<endl;
			for(size_t i=0; i<allNeighborVec.size(); i++ ){
				cout<<i<<": ";
				for(size_t j=0; j<allNeighborVec[i].size(); j++ ){
					cout<<allNeighborVec[i][j]<<" ";
				}
				cout<<endl;
			}
		}

		//find the largest cluster
		size_t maxClu = 0;
		size_t maxCluIndex = 0;
		for(size_t i=0; i<allNeighborVec.size(); i++ ){
			if( maxClu < allNeighborVec[i].size() ){
				maxClu = allNeighborVec[i].size();
				maxCluIndex = i;
			}
		}

		//add the largest cluster to classified
		resultVec.push_back( allNeighborVec[maxCluIndex] );

		if(0){
			cout<<"allNeighborVec:"<<endl;
			cout<<"maxCluIndex:"<<maxCluIndex<<endl;
			for(size_t i=0; i<allNeighborVec[maxCluIndex ].size(); i++ ){
				cout<<allNeighborVec[maxCluIndex ][i]<<" ";
			}
		}

		//remove the classified data from unclassifiedVec
		for(size_t i=0; i<allNeighborVec[maxCluIndex ].size(); i++ ){
			unclassifiedVec.erase( find(unclassifiedVec.begin(),
					unclassifiedVec.end(), allNeighborVec[maxCluIndex][i] ) );
		}
	}
	resultVec = sortCluster(resultVec );
	if(0){
		cout<<"2"<<endl;
		for(size_t i=0; i<resultVec.size(); i++){
			cout<<i<<" ";
			for(size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" "<<resultVec[i][j];
			}
			cout<<endl;
		}
	}
	return resultVec;
}

/**
 * @brief	Average distance between two clusters
 * @param inputMatrix 	The input pairwise distance matrix
 * @param clusterI		The first cluster vector
 * @param clusterJ		The second cluster vector
 *
 * compute the average distance between two clusters clusterI and clusterJ
 */
template<typename _Tp>
_Tp
averageLink(const vector<vector<_Tp> >& inputMatrix, const vector<size_t> clusterI,
		const vector<size_t> clusterJ){

	if( clusterI.size() > inputMatrix.size() ||
		clusterJ.size() > inputMatrix.size() ){
		cerr<<"clusterdata.h averageLink() error!!!"<<endl;
		return 0;
	}
	if( clusterI.size() == 0 || clusterJ.size() == 0 ){
		cerr<<"no data in cluster vector to compute !!!"<<endl;
		return 0;
	}

	_Tp averageDis = 0;


	for(size_t i=0; i<clusterI.size(); i++ ){
		for(size_t j=0; j<clusterJ.size(); j++ ){
			averageDis += inputMatrix[clusterI[i] ][clusterJ[j] ];
		}
	}

	if(clusterI.size() != 0 && clusterJ.size() != 0 ){
		averageDis /= (clusterI.size() * clusterJ.size() );
	}else{
		averageDis = 10e10;
	}

	return averageDis;
}

/**
 * @brief		Compute the largest distance between any two structures in a group
 * @param	inputMatrix The pariwise distance of all structures
 * @param	cluster 	Structure index
 */
template< typename _Tp >
_Tp
largestInnerDis( const vector< vector< _Tp > >& inputMatrix, const vector< size_t > cluster  ){
	_Tp maxDis = 0;
	if( cluster.empty() ){
		cout<<"---error largestInnerDis---"<<endl;
	}

	for( size_t i=0; i<cluster.size(); i++ ){
		for( size_t j=0; j<i; j++ ){
			if( cluster.at(i) > inputMatrix.size() || cluster.at(j) > inputMatrix.size() ){
				cerr<<"result invalid, cluster index "<<cluster.at(i)<<":"<<cluster.at(j)
						<<" exceed the size of distance matrix"<<inputMatrix.size()<<endl;
				return 0;
			}
			if( maxDis < inputMatrix[ cluster.at(i) ][ cluster.at(j) ] ){
				maxDis = inputMatrix[ cluster.at(i) ][ cluster.at(j) ];
			}
		}
	}

	return  maxDis;
}

/**
 * @brief		Hierarchical clustering algorithm using average-link
 * @param	inputDis 	The input distance matrix deposited in one dimension
 * @param	limiteDis	The criteria distance for hierarchical algorithm
 * @param	clusterResult	The cluster result, deposited in two dimension
 * @return	true if success false if failed
 * @note		The 'inputDis' contains the pairwise distance of all structures in one dimension
 *              	vector, where $ i * (i + 1) + j $th (i>j) item corresponding to the [i][j]th item of
 *              	distance matrix of two dimension
 *
  *The correspondance between the one dimension 'inputDis' and two dimension distance matrix is as below
 * 	|    0                 0                  0                  0                0     |
 * 	| inputDis[0]     0                  0                  0                0     |
 * 	| inputDis[1] inputDis[2]       0                  0                0     |
 * 	| inputDis[3] inputDis[4] inputDis[5]         0                0     |
 * 	| inputDis[6] inputDis[7] inputDis[8]  inputDis[9]        0     |
 */
template<typename _Tp>
bool
Hierarchical( const vector<_Tp> inputDis, const _Tp limiteDis,
		vector< vector< size_t > >& clusterResult  ){
	cout<<"---Hierarchical---"<<endl;

	size_t 	structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	// make sure the structNum is correct
	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
		return false;
	}

	vector< vector<_Tp> > inputMatrix(structNum, vector<_Tp>( structNum ));

	// initialize input distance matrix
	inputMatrix[0][0] = 0;
	for( size_t i=1; i<structNum; i++ ){
		for( size_t j=0; j<i; j++ ){
			inputMatrix[i][j] = inputMatrix[j][i] = inputDis[ i*(i-1)/2 + j  ];
		}
		inputMatrix[i][i] = 0;
	}

	clusterResult.clear();
	clusterResult = Hierarchical( inputMatrix, limiteDis );

	string	fileName = "hierarchical_one.txt";
	printfCluster( clusterResult, fileName );
	return true;
}

/**
 * @biref hierarchical clustering algorithm using average link
 * @param inputMatrix The input pairwise distance matrix
 * @param limiteDis The cutoff distance
 * @return cluster result
 */
template<typename _Tp>
vector<vector<size_t> >
Hierarchical2(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	cout<<"---hierarchical2---"<<endl;

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for Hierarchical!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	//the two clusters that have the smallest average distance
	size_t					clusterI, clusterJ;
	clusterI = clusterJ = 0;
	vector<size_t >			clusterIJVec;

	// sorted cluster distance vector, where the first pair of cluster in
	// 'sortedPairCluster' has the smallest distance, and the last pair
	// has the largest distance
	vector< clusterPair<_Tp> >  sortedPairCluster;

	//initialization
	for(size_t i=0; i<matrixSize; i++){
		resultVec.push_back( vector<size_t>(1, i) );
	}
	for(size_t i=0; i<resultVec.size(); i++ ){
		for(size_t j=0; j<i; j++ ){
			averageDis = averageLink(inputMatrix, resultVec[i], resultVec[j]);
			if( minDis > averageDis ){
				minDis = averageDis;

				//the smallest group index
				clusterI = i;
				clusterJ = j;

				//store the cluster pairs
				clusterPair<_Tp>	 	cluP;
				cluP.cluster_i = clusterI;
				cluP.cluster_j = clusterJ;
				cluP.distance = averageDis;
				sortedPairCluster.push_back(  cluP  );
			}
		}
	}
	sort( sortedPairCluster.begin(),  sortedPairCluster.end(), clusterPairComparator<_Tp> );

	while( minDis < limiteDis ){
		if(0){
			cout<<"minDis:"<<minDis<<endl;
			cout<<"clusterI:"<<clusterI<<endl;
			cout<<"clusterJ:"<<clusterJ<<endl;
			cout<<"resultVec.size():"<<resultVec.size()<<endl;
		}
		if(0){
			for( size_t i=0; i<resultVec.size(); i++ ){
				cout<<i<<": ";
				for( size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
		}

		//store the merged two clusters before they are deleted
		clusterIJVec.clear();
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterI ].begin(), resultVec[clusterI ].end());
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterJ ].begin(), resultVec[clusterJ ].end() );

		//delete clusteri clusterj in sortedPairCluster
		for( size_t i=0; i<sortedPairCluster.size(); i++ ){
			if( sortedPairCluster.at(i).cluster_i == clusterI || sortedPairCluster.at(i).cluster_i == clusterJ ||
				 sortedPairCluster.at(i).cluster_j == clusterI || sortedPairCluster.at(i).cluster_j == clusterJ){
				sortedPairCluster.erase( sortedPairCluster.begin() + i );
			}
		}

		if(0){
			cout<<"1"<<endl;
		}

		//erase the two merged cluster
		if(clusterI < clusterJ ){
			resultVec.erase( resultVec.begin() + clusterJ );
			resultVec.erase( resultVec.begin() + clusterI );
		}else{
			resultVec.erase( resultVec.begin() + clusterI );
			resultVec.erase( resultVec.begin() + clusterJ );
		}

		//combine the tow merged clusters into resultVec
		resultVec.push_back( clusterIJVec );

		if(0){
			cout<<"2"<<endl;
		}

		if(0){
			cout<<"---resultVec:---"<<endl;
			for( size_t i=0; i<resultVec.size(); i++ ){
				cout<<i<<": ";
				for( size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
			cout<<"------"<<endl;
		}

		if(0){
			cout<<"3"<<endl;
			cout<<"resultVec.size():"<<resultVec.size()<<endl;
		}

		//update sortedPairCluster
		for( size_t i=0; i<resultVec.size()-1; i++ ){
			if(0){
				cout<<i<<":  ";
			}

			_Tp clustDis = averageLink( inputMatrix,  resultVec[i], resultVec.back()  );

			//store the cluster pairs
			clusterPair<_Tp>	 	cluP;
			cluP.cluster_i = i;
			cluP.cluster_j = resultVec.size()  - 1;
			cluP.distance = clustDis;
			sortedPairCluster.push_back(  cluP  );
		}

		if(0){
			cout<<"4"<<endl;
		}

		sort( sortedPairCluster.begin(), sortedPairCluster.end(), clusterPairComparator<_Tp> );

		clusterI = sortedPairCluster.front().cluster_i;
		clusterJ = sortedPairCluster.front().cluster_j;

		if(0){
			cout<<"----sortedPairCluster---"<<endl;
			for( size_t i=0; i<sortedPairCluster.size(); i++ ){
				cout<<i<<"  ";
				cout<<sortedPairCluster.at(i).cluster_i<<":"<<sortedPairCluster.at(i).cluster_j
						<<"->"<<sortedPairCluster.at(i).distance<<"  ";
			}
		}


	}
	resultVec = sortCluster(resultVec );
	if(1){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"cluster"<<i<<": ";
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
			}
			cout<<endl;
		}
	}

	if(0){
		string	fileName = "hierarchical.txt";
		printfCluster( resultVec, fileName );
	}
	return resultVec;
}

/*
Hierarchical_fast( const vector< vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	cout<<"---hierarchical_fast---"<<endl;

}
*/

/**
 * @biref hierarchical clustering algorithm using average link
 * @param inputMatrix The input pairwise distance matrix
 * @param limiteDis The cutoff distance
 * @return cluster result
 */
template<typename _Tp>
vector<vector<size_t> >
Hierarchical(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	cout<<"---hierarchical---"<<endl;

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for Hierarchical!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	//the two clusters that have the smallest average distance
	size_t					clusterI, clusterJ;
	vector<size_t >			clusterIJVec;

	//initialization
	for(size_t i=0; i<matrixSize; i++){
		resultVec.push_back( vector<size_t>(1, i) );
	}
	for(size_t i=0; i<resultVec.size(); i++ ){
		for(size_t j=0; j<i; j++ ){
			averageDis = averageLink(inputMatrix, resultVec[i], resultVec[j]);
			if( minDis > averageDis ){
				minDis = averageDis;
				clusterI = i;
				clusterJ = j;
			}
		}
	}

	while( minDis < limiteDis ){
		if(0){
			cout<<"minDis:"<<minDis<<endl;
			cout<<"clusterI:"<<clusterI<<endl;
			cout<<"clusterJ:"<<clusterJ<<endl;
			cout<<"resultVec.size():"<<resultVec.size()<<endl;
		}
		if(0){
			for( size_t i=0; i<resultVec.size(); i++ ){
				cout<<i<<": ";
				for( size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
		}

		//store the merged two clusters before they are deleted
		clusterIJVec.clear();
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterI ].begin(), resultVec[clusterI ].end());
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterJ ].begin(), resultVec[clusterJ ].end() );

		//erase the two merged cluster
		if(clusterI < clusterJ ){
			resultVec.erase( resultVec.begin() + clusterJ );
			resultVec.erase( resultVec.begin() + clusterI );
		}else{
			resultVec.erase( resultVec.begin() + clusterI );
			resultVec.erase( resultVec.begin() + clusterJ );
		}

		//combine the tow merged clusters into resultVec
		resultVec.push_back( clusterIJVec );

		minDis = 10e10;
		for( size_t i=1; i<resultVec.size(); i++ ){
			for(size_t j=0; j<i; j++ ){
				if( i!= j ){
					if( 1 == resultVec[i].size() && 1 == resultVec[j].size() ){
						averageDis = inputMatrix[ resultVec[i].front() ][ resultVec[j].front() ];
					}else{
						averageDis = averageLink(inputMatrix, resultVec[i], resultVec[j] );
					}
					if(minDis > averageDis){
						minDis = averageDis;
						clusterI = i;
						clusterJ = j;
					}
				}
			}
		}
	}
	resultVec = sortCluster(resultVec );
	if(1){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"cluster"<<i<<": ";
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
			}
			cout<<endl;
		}
	}

	if(0){
		string	fileName = "hierarchical.txt";
		printfCluster( resultVec, fileName );
	}
	return resultVec;
}

/**
 * @brief	This clustering algorithm is similar as hierarchical algorithm, while it ensures the largest
 * 				inner distance of each cluster is smaller than a radius
 * 	@param	inputDis The input distance matrix deposited in one dimension
 * 	@param	limiteDis	The distance criteria for this algorithm
 * 	@param	radius		The largest inner distance of each cluster, typically the same value as limiteDis
 * 	@param	clusterResult		The clustering result
 * 	@return	true if success, false if failed
 *  					where $ ( i * ( i + 1 ) / 2) + j $th ( i>j ) item corresponding to  the [i][j]th item  in two dimension vectors.
 *
 * The correspondance between the one dimension 'inputDis' and two dimension distance matrix is as below
 * 	|    0                 0                  0                  0                0     |
 * 	| inputDis[0]     0                  0                  0                0     |
 * 	| inputDis[1] inputDis[2]       0                  0                0     |
 * 	| inputDis[3] inputDis[4] inputDis[5]         0                0     |
 * 	| inputDis[6] inputDis[7] inputDis[8]  inputDis[9]        0     |
 */
template<typename _Tp>
bool
Hierarchical_limit( const vector< _Tp > inputDis, const _Tp& limiteDis, const _Tp& radius,
								vector< vector< size_t > >& clusterResult ){
	cout<<"----Hierarchical_limite--"<<endl;

	size_t structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
		return false;
	}

	vector< vector<_Tp> > inputMatrix(structNum, vector<_Tp>( structNum ));

	// initialize input distance matrix
	inputMatrix[0][0] = 0;
	for( size_t i=1; i<structNum; i++ ){
		for( size_t j=0; j<i; j++ ){
			inputMatrix[i][j] = inputMatrix[j][i] = inputDis[ i*(i-1)/2 + j  ];
		}
		inputMatrix[i][i] = 0;
	}

	clusterResult.clear();
	clusterResult = Hierarchical_limit( inputMatrix, limiteDis, radius );

	return true;
}

/**
 * @biref hierarchical clustering algorithm using average link
 * @param inputMatrix The input pairwise distance matrix
 * @param limiteDis The cutoff distance
 * @param radius The largest distance inner distance
 * @return cluster result
 *
 * different with $Hierarchical$, this function merge two clusters not only when
 * they have the smallest average-link distance but also distance of any two are
 * uner a radius.
 */
template<typename _Tp>
vector<vector<size_t> >
Hierarchical_limit(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis,
		const _Tp& radius){
	clock_t t1, t2;
	t1 = clock();
	cout<<"---hierarchical_limit---"<<endl;

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for Hierarchical_limit!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	//the two clusters that have the smallest average distance
	size_t					clusterI, clusterJ;
	vector<size_t >			clusterIJVec;

	//initialization
	for(size_t i=0; i<matrixSize; i++){
		resultVec.push_back( vector<size_t>(1, i) );
	}
	for(size_t i=0; i<resultVec.size(); i++ ){
		for(size_t j=0; j<i; j++ ){
			averageDis = averageLink(inputMatrix, resultVec[i], resultVec[j]);
			if( minDis > averageDis ){
				minDis = averageDis;
				clusterI = i;
				clusterJ = j;
			}
		}
	}
	if(0){
		cout<<"---inputMatrix---"<<endl;
		ofstream of( "inputMatrix.txt" );
		for( size_t i=0; i<inputMatrix.size(); i++ ){
			of.width(5);
			of<<i<<" ";
			for( size_t j=0; j<inputMatrix[i].size(); j++ ){
				of.width(15);
				of<<inputMatrix[i][j]<<" ";
			}
			of<<endl;
		}
		of.close();
	}

	while( minDis < limiteDis/* && largestDis(inputMatrix,
			resultVec[clusterI ], resultVec[clusterJ ]) < radius*/ ){
		if(0){
			cout<<"minDis:"<<minDis<<endl;
			cout<<"clusterI:"<<clusterI<<endl;
			cout<<"clusterJ:"<<clusterJ<<endl;
			cout<<"resultVec.size():"<<resultVec.size()<<endl;
			cout<<"largestDis:"<<largestDis(inputMatrix,
					resultVec[clusterI ], resultVec[clusterJ ]);
			cout<<"limiteDis:"<<limiteDis<<endl;
		}
		if(0){
			for( size_t i=0; i<resultVec.size(); i++ ){
				cout<<i<<": ";
				for( size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
		}

		//store the merged two clusters before they are deleted
		clusterIJVec.clear();
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterI ].begin(), resultVec[clusterI ].end());
		clusterIJVec.insert(clusterIJVec.end(),
				resultVec[clusterJ ].begin(), resultVec[clusterJ ].end() );

		//erase the two merged cluster
		if(clusterI < clusterJ ){
			resultVec.erase( resultVec.begin() + clusterJ );
			resultVec.erase( resultVec.begin() + clusterI );
		}else{
			resultVec.erase( resultVec.begin() + clusterI );
			resultVec.erase( resultVec.begin() + clusterJ );
		}

		//combine the tow merged clusters into resultVec
		resultVec.push_back( clusterIJVec );

		minDis = 10e10;
		for( size_t i=0; i<resultVec.size(); i++ ){
			for(size_t j=0; j<resultVec.size(); j++ ){
				if( i!= j ){
					averageDis = averageLink(inputMatrix, resultVec[i], resultVec[j] );
					if(minDis > averageDis &&  largestDis(inputMatrix,
							resultVec[ i ], resultVec[ j ]) < radius  ){
						minDis = averageDis;
						clusterI = i;
						clusterJ = j;
					}
				}
			}
		}
		if(0){
			cout<<"minDis:"<<minDis<<endl;
		}
	}
	if(1){
		cout<<"t1:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	resultVec = sortCluster(resultVec );
	if(0){
		cout<<"2"<<endl;
		for(size_t i=0; i<resultVec.size(); i++){
			cout<<i<<" ";
			for(size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<" "<<resultVec[i][j];
			}
			cout<<endl;
		}
	}
	if(0){
		string	fileName = "hierarchical_limit.txt";
		printfCluster( resultVec, fileName );
	}

	if(0){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"cluster"<<i<<": ";
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
			}
			cout<<endl;
		}
	}

//	printCluster( resultVec );
	printCluster( resultVec, inputMatrix );
	return resultVec;
}

/**
 * @note this struct is particular for singleLink
 *  */
template<typename _Tp>
struct similarity{
public:
	//similarity
	_Tp		sim;
	size_t 	index;
public:
	bool operator()( const similarity & input )const{
		return (( input.sim == sim) && (input.index == index) );
	}
};

/**
 * @note	$similarity$ comparator based on sim
 */
template<typename _Tp>
bool similarityComp( const similarity<_Tp> & sim1,
		const similarity<_Tp>& sim2 ){
	return sim1.sim<sim2.sim ? true:false;
};

template<typename _Tp>
bool similarityNoZeroComp( const similarity<_Tp> & sim1,
		const similarity<_Tp>& sim2 ){
	if( sim1.sim ==0 || sim2.sim == 0 ){
		return false;
	}
	return sim1.sim<sim2.sim ? true:false;
};

/**
 * @brief	singlinkClustering
 * @param	inputMatrix distance input matrix
 * @param	limiteDis distance criteria
 * @note	the value in inputMatrix is distance between points
 */
template<typename _Tp>
vector<vector<size_t> >
singleLink(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	clock_t t1, t2;
	t1 = clock();
	cout<<endl<<"---singleLink Clustering---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for singleLink clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	// I
	vector< size_t >	I( matrixSize, 0 );

	// next-best-merge vector
	vector< similarity<_Tp>  >  nbm;
	for( size_t i=0; i<matrixSize; i++ ){
		nbm.push_back( similarity< _Tp >()  );
	}

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back(  c  );
	}

	if(0){
		cout<<"t2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}


	//initialization
	for( size_t n=0; n< matrixSize; n++ ){
		// determin the struct in C that have the smallest value sim, and give it to nbm[n]
		_Tp minSim = 10e10;
		size_t minIndex;

		for( size_t i=0; i< matrixSize; i++ ){
			//initialize C
			C[n][i].sim = inputMatrix[n][i];
			C[n][i].index = i;

			//get nbm[n]
			if( C[n][i].sim < minSim && C[n][i].sim != 0 ){
				minSim = C[n][i].sim;
				minIndex = i;
			}
		}
		I[ n ] = n;
		nbm[ n ] = C[n][ minIndex ];
	}

	//printf
	if(0){
		cout<<"nbm:"<<endl;
		for( size_t i=0; i<nbm.size(); i++ ){
			cout<<nbm[i].sim<<" ";
		}
		cout<<endl;
	}

	//A
	vector< vector< size_t > > A;
	_Tp minV = -10e10;

	for( size_t n=0;  n<matrixSize-1 ; n++ ){
		if(1){
			char ss[20] = "";
			for( size_t i =0; i<size_t(20*((float) n/matrixSize ) ); i++ ){
				ss[i] = '-';
			}
			if(size_t(20*((float) n/matrixSize ) )  <19  ){
				ss[ size_t(20*((float) n/matrixSize ) ) ] = '>';
			}
			for( size_t i=size_t(20*((float) n/matrixSize ) ) + 1; i<19; i++ ){
				ss[i] = ' ';
			}
			cerr<<"\r|"<<ss<<"|"<<" running time: "<< double( clock() - t1 )/ CLOCKS_PER_SEC<<" ";
		}
		if( minV < limiteDis ){

			size_t  i1, i2;

			// find i1, which is the index of smallest value in nbm[]
			minV = 10e10;
			size_t minI;
			for( size_t i=0; i<matrixSize; i++ ){
				if( nbm[i].sim < minV && I[i] == i  ){
					minV = nbm[i].sim;
					minI = i;
				}
			}
			i1 = minI;
			i2 = I[ nbm[i1].index ];

			if(0){
				cout<<endl<<endl;
				cout<<"i1:"<<i1<<"  i2:"<<i2<<endl;
			}

			vector<size_t> i12;
			i12.push_back(i1);
			i12.push_back(i2);
			A.push_back( i12 );

			if(0){
				cout<<"I:  ";
				for( size_t index = 0; index< I.size() ; index++ ){
					cout<<I[index]<<" ";
				}
				cout<<endl;
			}

			//update matrix C
			for( size_t i=0; i<matrixSize; i++ ){
				if( I[ i ] == i && i != i1 && i != i2  ){
					C[ i1 ][i].sim = C[ i ][ i1 ].sim = min( C[i1][i].sim, C[i2][i].sim );
					if(0){
						cout<<"i1:"<<i1<<"  i:"<<i<<endl;
					}
				}

				if( I[i] == i2 ){
					I[ i ] = i1;
				}
			}

			if(0){
				cout<<"I: ";
				for( size_t index = 0; index<I.size(); index++ ){
					cout<<I[ index ]<<" ";
				}
				cout<<endl;
			}

			if(0){
				cout<<"C:"<<endl;
				for( size_t index1 = 0; index1<C.size(); index1++ ){
					for( size_t index2 = 0; index2<C[index1].size(); index2++ ){
						cout.width(8 );
						cout<<C[index1][index2].sim<<" ";
					}
					cout<<endl;
				}
			}

			_Tp minSim = 10e10;
			for( size_t i=0; i<matrixSize; i++ ){
				if( I[i] == i && i != i1 && C[i1][i].sim < minSim ){
					minSim = C[i1][i].sim;
					minI = i;
				}
			}
			nbm[ i1 ] = C[i1][ minI ];

			if(0){
				cout<<"nbm:";
				for( size_t index=0; index<nbm.size(); index++ ){
					//				cout<<nbm[ index ].index<<","<<nbm[index].sim<<"   ";
					cout<<nbm[index].sim<<"   ";
				}
				cout<<endl;
			}

			//for comparasion with limiteDis to terminate the for loop
			minV = 10e10;
			for( size_t i=0; i<matrixSize; i++ ){
				if( nbm[i].sim < minV && I[i] == i  ){
					minV = nbm[i].sim;
				}
			}
		}
	}

	//get clustering result from vector I
	vector<size_t > tempV( I );
	vector<size_t > oneCluster;

	//flag == true means the value in tempV[i] has not shown before
	bool flag = true;

	resultVec.clear();
	for( size_t i=0; i<tempV.size(); i++ ){
		oneCluster.clear();
		flag = true;
		for( size_t j=0; j<i; j++ ){
			if( tempV[ j ]  == tempV[i]  ){
				flag = false;
			}
		}
		for( size_t j=i; j<tempV.size(); j++ ){
			if( flag && tempV[ i ] == tempV[ j ] ){
				oneCluster.push_back( j );
			}
		}
		if( ! oneCluster.empty() ){
			resultVec.push_back( oneCluster );
		}
	}

	if(1){
		cout<<"singleLinke running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	resultVec = sortCluster( resultVec );

	if(0){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"cluster"<<i<<": ";
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
			}
			cout<<endl;
		}
	}
	printCluster( resultVec, inputMatrix );

	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * @brief		similiar as single-link, however this algorithms ensures the inner distance is under $limiteDis$
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector<vector<size_t> >
singleLink_limit(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	clock_t t1, t2;
	t1 = clock();
	cout<<endl<<"---singleLink Clustering limite---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for singleLink limite clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	// I
	vector< size_t >	I( matrixSize, 0 );

	// next-best-merge vector
	vector< similarity<_Tp>  >  nbm;
	for( size_t i=0; i<matrixSize; i++ ){
		nbm.push_back( similarity< _Tp >()  );
	}

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back(  c  );
	}

	if(0){
		cout<<"t2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}


	//initialization
	for( size_t n=0; n< matrixSize; n++ ){
		// determin the struct in C that have the smallest value sim, and give it to nbm[n]
		_Tp minSim = 10e10;
		size_t minIndex;

		for( size_t i=0; i< matrixSize; i++ ){
			//initialize C
			C[n][i].sim = inputMatrix[n][i];
			C[n][i].index = i;

			//get nbm[n]
			if( C[n][i].sim < minSim && C[n][i].sim != 0 ){
				minSim = C[n][i].sim;
				minIndex = i;
			}
		}
		I[ n ] = n;
		nbm[ n ] = C[n][ minIndex ];
	}

	//printf
	if(0){
		cout<<"nbm:"<<endl;
		for( size_t i=0; i<nbm.size(); i++ ){
			cout<<nbm[i].sim<<" ";
		}
		cout<<endl;
	}

	//A
	vector< vector< size_t > > A;
	_Tp minV = -10e10;

	for( size_t n=0;  n<matrixSize-1  ; n++ ){
		if(1){
			char ss[20] = "";
			for( size_t i =0; i<size_t(20*((float) n/matrixSize ) ); i++ ){
				ss[i] = '-';
			}
			if(size_t(20*((float) n/matrixSize ) )  <19  ){
				ss[ size_t(20*((float) n/matrixSize ) ) ] = '>';
			}
			for( size_t i=size_t(20*((float) n/matrixSize ) ) + 1; i<19; i++ ){
				ss[i] = ' ';
			}
			cerr<<"\r|"<<ss<<"|"<<" running time: "<< double( clock() - t1 )/ CLOCKS_PER_SEC<<" ";
		}
		if( minV < limiteDis ){

			size_t  i1, i2;

			// find i1, which is the index of smallest value in nbm[]
			minV = 10e10;
			size_t minI;
			for( size_t i=0; i<matrixSize; i++ ){
				if( nbm[i].sim < minV && I[i] == i  ){
					minV = nbm[i].sim;
					minI = i;
				}
			}
			i1 = minI;
			i2 = I[ nbm[i1].index ];

			if(0){
				cout<<endl<<endl;
				cout<<"i1:"<<i1<<"  i2:"<<i2<<endl;
			}

			vector<size_t> i12;
			i12.push_back(i1);
			i12.push_back(i2);
			A.push_back( i12 );

			if(0){
				cout<<"I:  ";
				for( size_t index = 0; index< I.size() ; index++ ){
					cout<<I[index]<<" ";
				}
				cout<<endl;
			}

			//update matrix C
			for( size_t i=0; i<matrixSize; i++ ){
				if( I[ i ] == i && i != i1 && i != i2  ){
					//				C[ i1 ][i].sim = C[ i ][ i1 ].sim = min( C[i1][i].sim, C[i2][i].sim );
					//note: here is different from the $singleLink$
					C[ i1 ][i].sim = C[ i ][ i1 ].sim = max( C[i1][i].sim, C[i2][i].sim );

					//note: here is different from the $singleLink$, which has no this step
					//update the distance of all the rest data to the new cluster by the largest distance
					C[i2][i].sim = C[i][i2].sim = max( C[i1][i].sim, C[i2][i].sim );

					if(0){
						cout<<"i:"<<i<<" í1:"<<i1<<"  i2:"<<i2<<" C[i][i1]:"<<C[i1][i].sim<<" C[i2][i]:"<<C[i2][i].sim<<endl;
					}
				}

				if( I[i] == i2 ){
					I[ i ] = i1;
				}
			}

			if(0){
				cout<<"I: ";
				for( size_t index = 0; index<I.size(); index++ ){
					cout<<I[ index ]<<" ";
				}
				cout<<endl;
			}

			if(0){
				cout<<"C:"<<endl;
				for( size_t index1 = 0; index1<C.size(); index1++ ){
					for( size_t index2 = 0; index2<C[index1].size(); index2++ ){
						cout.width(8 );
						cout<<C[index1][index2].sim<<" ";
					}
					cout<<endl;
				}
			}

			_Tp minSim = 10e10;
			size_t minJ = 0;
			for( size_t i=0; i<matrixSize; i++ ){
				minSim = 10e10;
				for( size_t j=0; j<matrixSize; j++ ){
					if( I[j] == j && j != i && C[i][ j ].sim < minSim ){
						minSim = C[i][j].sim;
						minJ = j;
					}
				}
				nbm[ i ] = C[i][ minJ ];
			}

			if(0){
				cout<<"nbm:";
				for( size_t index=0; index<nbm.size(); index++ ){
					//				cout<<nbm[ index ].index<<","<<nbm[index].sim<<"   ";
					cout<<nbm[index].sim<<"   ";
				}
				cout<<endl;
			}

			//for comparasion with limiteDis to terminate the for loop
			minV = 10e10;
			for( size_t i=0; i<matrixSize; i++ ){
				if( nbm[i].sim < minV && I[i] == i  ){
					minV = nbm[i].sim;
				}
			}
			if(0){
				cout<<"minV:"<<minV<<endl;
			}
		}
	}
	cerr<<endl;
	//get clustering result from vector I
	vector<size_t > tempV( I );
	vector<size_t > oneCluster;

	//flag == true means the value in tempV[i] has not shown before
	bool flag = true;

	if(0){
		cout<<"I:"<<endl;
		for( size_t i=0; i<I.size(); i++ ){
			cout<<"I"<<i<<": "<<I[i]<<endl;
		}
	}

	//get clustering result from vector I
	resultVec.clear();
	for( size_t i=0; i<tempV.size(); i++ ){
		oneCluster.clear();
		flag = true;
		for( size_t j=0; j<i; j++ ){
			if( tempV[ j ]  == tempV[i]  ){
				flag = false;
			}
		}
		for( size_t j=i; j<tempV.size(); j++ ){
			if( flag && tempV[ i ] == tempV[ j ] ){
				oneCluster.push_back( j );
			}
		}
		if( ! oneCluster.empty() ){
			resultVec.push_back( oneCluster );
		}
	}

	if(1){
		cout<<"singleLink_limit running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	resultVec = sortCluster( resultVec );

	if(0){
		for( size_t i=0; i<resultVec.size(); i++ ){
			cout<<"cluster"<<i<<": ";
			for( size_t j=0; j<resultVec[i].size(); j++ ){
				cout<<resultVec[i][j]<<"  ";
			}
			cout<<endl;
		}
	}
	printCluster( resultVec, inputMatrix );

	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );
	vector< vector< vector<size_t> > > s = getScratch( resultVec, inputMatrix, limiteDis );

	return resultVec;
}

/**
 * @brief	singlinkClustering
 * @note  the inputMarix contains similarity value, which means higher value are merge first
 */
template<typename _Tp>
vector<vector<size_t> >
singleLink2(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis ){
	cout<<"---singleLink Clustering---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for singleLink clustering!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 10e10;

	// I
	vector< size_t >	I( matrixSize, 0 );

	// next-best-merge vector
	vector< similarity<_Tp>  >  nbm;
	for( size_t i=0; i<matrixSize; i++ ){
		nbm.push_back( similarity< _Tp >()  );
	}

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back(  c  );
	}

	//initialization
	for( size_t n=0; n< matrixSize; n++ ){
		for( size_t i=0; i< matrixSize; i++ ){
			C[n][i].sim = inputMatrix[n][i];
			C[n][i].index = i;
		}
		I[ n ] = n;
		nbm[ n ] = *max_element( C[n].begin(), C[n].end(), similarityComp<_Tp>  );

	}

	//printf
	if(1){
		cout<<"nbm:"<<endl;
		for( size_t i=0; i<nbm.size(); i++ ){
			cout<<nbm[i].sim<<" ";
		}
		cout<<endl;
	}

	//A
	vector< vector< size_t > > A;

	for( size_t n=0; n<matrixSize-1; n++ ){
		size_t  i1, i2;

		// find i1
		_Tp maxV = -10e10;
		size_t maxI;
		for( size_t i=0; i<matrixSize; i++ ){
			if( nbm[i].sim > maxV && I[i] == i  ){
				maxV = nbm[i].sim;
				maxI = i;
			}
		}
		i1 = maxI;
		i2 = I[ nbm[i1].index ];

		if(1){
			cout<<endl<<endl;
			cout<<"i1:"<<i1<<"  i2:"<<i2<<endl;
		}

		vector<size_t> i12;
		i12.push_back(i1);
		i12.push_back(i2);
		A.push_back( i12 );

		if(1){
			cout<<"I:  ";
			for( size_t index = 0; index< I.size() ; index++ ){
				cout<<I[index]<<" ";
			}
			cout<<endl;
		}

		//update matrix C
		for( size_t i=0; i<matrixSize; i++ ){
			if( I[ i ] == i && i != i1 && i != i2  ){
				C[ i1 ][i].sim = C[ i ][ i1 ].sim = max( C[i1][i].sim, C[i2][i].sim );
				if(1){
					cout<<"i1:"<<i1<<"  i:"<<i<<endl;
				}
			}

			if( I[i] == i2 ){
				I[ i ] = i1;
			}
		}

		if(1){
			cout<<"I: ";
			for( size_t index = 0; index<I.size(); index++ ){
				cout<<I[ index ]<<" ";
			}
			cout<<endl;
		}

		if(1){
			cout<<"C:"<<endl;
			for( size_t index1 = 0; index1<C.size(); index1++ ){
				for( size_t index2 = 0; index2<C[index1].size(); index2++ ){
					cout.width(8 );
					cout<<C[index1][index2].sim<<" ";
				}
				cout<<endl;
			}
		}

		_Tp maxSim = -10e10;
		for( size_t i=0; i<matrixSize; i++ ){
			if( I[i] == i && i != i1 && C[i1][i].sim > maxSim ){
				maxSim = C[i1][i].sim;
				maxI = i;
			}
		}
		nbm[ i1 ] = C[i1][ maxI ];

		if(1){
			cout<<"nbm:";
			for( size_t index=0; index<nbm.size(); index++ ){
				cout<<nbm[ index ].index<<","<<nbm[index].sim<<"   ";
			}
			cout<<endl;
		}
	}

	return A;
}

/**
 * @param i
 * @param k1
 * @param k2
 * @param inputMatrix
 * @return  max( sim(i, k1), sim(i, k2) )
 */
template<typename _Tp>
_Tp
completelinkDis( const size_t i, const size_t k1, const size_t k2, const vector< vector< _Tp > >& inputMatrix ){
	return inputMatrix[ i ][ k1 ] > inputMatrix[ i ][ k2 ]? inputMatrix[ i ][ k1 ] : inputMatrix[ i ][ k2 ];
}

template<typename _Tp>
_Tp
completelinkDis( const size_t i, const size_t k1, const size_t k2, const vector< vector<similarity<_Tp>  > > piVec ){
	if( piVec.empty() || i > piVec.size() ){
		cerr<<"---error to compute complete link distance---"<<endl;
	}

	if(0){
		cout<<"completelinkDis: "<<i<<endl;
	}

	_Tp dis1 = 0;
	_Tp dis2 = 0;
	bool flag1 = true;
	bool flag2 = true;
	for( size_t m=0; m<piVec[i].size() && ( flag1 || flag2 ); m++ ){
		if( piVec[i][m].index == k1 ){
			dis1 = piVec[i][m].sim;
			flag1 = false;
		}else if( piVec[i][m].index == k2 ){
			dis2 = piVec[i][m].sim;
			flag2 = false;
		}
	}

	if( dis1 == 0 || dis2 ==0  ){
		cerr<<"----warnning distance is 0 in "<<i<<" "<<k1<<" "<<k2<<"---"<<endl;
	}

	return dis1 > dis2 ? dis1:dis2;
}

template< typename _Tp >
_Tp
completelinkDis( const size_t k1, const size_t k2, const vector< similarity<_Tp> >& pVec ){
	if( pVec.empty() ){
		cerr<<"---error to compute complete link distance---"<<endl;
	}

	_Tp dis1 = 0;
	_Tp dis2 = 0;
	size_t count = 0;
	size_t i=0;
	while( count != 2 && i<pVec.size() ){
		if( pVec[i].index == k1 ){
			dis1 = pVec[i].sim;
			count++;
			i++;
		}else if( pVec[i].index == k2 ){
			dis2 = pVec[i].sim;
			count++;
			i++;
		}else{
			i++;
		}
	}
	return dis1 > dis2 ? dis1:dis2;
}

/**
 * @param i
 * @param k1
 * @param k2
 * @param inputMatrix
 * @return min( sim(i, k1), sim(i, k2) )
 */
template<typename _Tp>
_Tp
singlelinkDis( const size_t i, const size_t k1, const size_t k2, const vector< vector< _Tp > >& inputMatrix ){
	return inputMatrix[ i ][ k1 ] < inputMatrix[ i ][ k2 ]? inputMatrix[ i ][ k1 ] : inputMatrix[ i ][ k2 ];
}

/**
 * @brief		the fast complete linkage algorith( CLINK )
 * 					http://nlp.stanford.edu/IR-book/html/htmledition/time-complexity-of-hac-1.html
 * @param inputMatrix
 * @param limiteDis
  * @param splitFlag   indicate whether to split the result until no cluster is mergeable to other clusters
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
completeLink( const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis, const bool& splitFlag ){
	//check running time
	clock_t t1, t2;
	t1 = clock();

	cout<<endl<<"---completeLink algorithm---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for complete Link clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;

	//$minDis$ indicate the minmum distance, $minDis2$ indicate the second minimum distance
	//used for stopping hierarchial procedure
	_Tp						minDis = 0;
	_Tp						minDis2 = 0;

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//priority-queue
	vector< vector< similarity<_Tp> > > P;
	vector<_Tp> p;
	vector< size_t > I;

	if(0){
		cout<<"t1:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back( c );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		for( size_t j=0; j<matrixSize; j++ ){
			C[i][j].sim = inputMatrix[i][j];
			C[i][j].index = j;
		}
		I.push_back( 1 );
	}

	if(0){
		cout<<"t2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//initialize P
	for( size_t i=0; i<matrixSize; i++ ){
		c = C[i];
		sort( c.begin(), c.end(), similarityComp<_Tp> );
		c.erase( c.begin() );
		P.push_back( c );
	}

	vector< vector< size_t > > A;

	if(0){
		cout<<"t3:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//main loop
	for( size_t k=1; (k<matrixSize)  ; k++ ){
		if(1){
			char ss[20] = "";
			for( size_t i =0; i<size_t(20*((float) k/matrixSize ) ); i++ ){
				ss[i] = '-';
			}
			if(size_t(20*((float) k/matrixSize ) )  <19  ){
				ss[ size_t(20*((float) k/matrixSize ) ) ] = '>';
			}
			for( size_t i=size_t(20*((float) k/matrixSize ) ) + 1; i<19; i++ ){
				ss[i] = ' ';
			}
			cerr<<"\r|"<<ss<<"|"<<" complete link running time: "<< double( clock() - t1 )/ CLOCKS_PER_SEC;
		}
		//
		c.clear();
		for( size_t i=0; i<matrixSize; i++ ){
			c.push_back( P[i].front() );
		}

		//find out the index k1 and k2 which have the smallest distance
		size_t k1, k2, tempk;

		//the smallest and the second smallest distance
		minDis = 10e10;
		minDis2 = 10e10;

		//get the least distance and its index
		for( size_t i=0; i<matrixSize; i++ ){
			if(0){
				cout<<"c[i].sim:"<<c[i].sim<<endl;
			}
			if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
				minDis2 = minDis;
				minDis = c[i].sim;
				tempk = i;
			}else if( minDis2 > c[i].sim && minDis != c[i].sim &&  ( c[i].sim != 0 ) & ( I[i] != 0 )  ){
				minDis2 = c[i].sim;
			}

			if(0){
				cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
			}
		}
		k1 = tempk;
		k2 = P[ k1 ].front().index;

		if(0){
			cout<<"  t3.1:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
		}
		if(0){
			cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<" limiteDis:"<<limiteDis<<endl;
		}

		if(0){
			cout<<"P"<<endl;
			for( size_t i=0; i<P.size(); i++ ){
				cout.width(5);
				cout<<i<<" ";
				for( size_t j=0; j<P[i].size(); j++ ){
					cout.width(5);
					cout<<right<<P[i][j].index<<":";
					cout.width(10);
					cout<<left<<P[i][j].sim<<" ";
				}
				cout<<endl;
			}
		}

		if(0){
			cout<<"k1:"<<k1<<" k2:"<<k2<<endl;
		}
		if(0){
			cout<<"  t3.3:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
		}

		//index k2 labled removed
		I[k2] = 0;

		//clear P[k1], then update it
		P[k1].clear();

		//update P
		for( size_t i=0; i<matrixSize  ; i++ ){
			if(  1 == I[i] && i != k1 ){
//				_Tp dis = completelinkDis( i, k1, k2, P );
				_Tp dis = 0;
				dis = completelinkDis( k1, k2, P[i] );

				//delete C[i][k1] in P
				typename vector< similarity<_Tp> >::iterator it ;
				it = find_if( P[i].begin(), P[i].end(), C[i][k1] );
				P[i].erase( it );

				//delete C[i][k2] in P
				it =  find_if( P[i].begin(), P[i].end(), C[i][k2] );
				if( it != P[i].end() ){
					P[i].erase( it );
				}

				//to way to insert, different in running time
				bool test = false;
				//update the distance between the ith data and k1, which is the label of cluster that contains k1 and k2
				//					C[ i ][ k1 ].sim = completelinkDis( i, k1, k2, inputMatrix );
				C[ i ][ k1 ].sim = dis;

				if(0){
					cout<<"k2:"<<k2<<endl;
					cout<<"k1:"<<1<<" i:"<<i<<"  C[ i ][ k1 ]:"<<C[ i ][ k1 ].sim<<endl;
				}

				if( test ){
					P[i].push_back( C[i][k1] );
					sort( P[i].begin(), P[i].end(), similarityComp<_Tp> );
				}
				if( !test ){
					it = P[i].begin();
					//					cout<<" (*it).sim: "<<(*it).sim<<endl;
					while( (*it).sim < C[ i ][ k1 ].sim && ( it != P[ i ].end() ) ){
						it ++;
						//						cout<<" (*it).sim: "<<(*it).sim<<endl;
					}
					if( it != P[ i ].end() ){
						P[ i ].insert( it, C[ i ][ k1 ] );
					}else{
						P[ i ].push_back( C[ i ][ k1 ] );
					}
				}
				if(0){
					cout<<"C[ i ][ k1 ]:"<<C[ i ][ k1 ].sim<<endl;
					cout<<"P"<<i;
					for( size_t j=0; j<P[i].size(); j++ ){
						cout<<" "<<P[i][j].sim;
					}
					cout<<endl;
				}

				//update
				C[ k1 ][ i ].sim = C[ i ][ k1 ].sim ;
				if( test ){
					P[ k1 ].push_back( C[k1][ i ] );
					sort( P[ k1 ].begin(), P[ k1 ].end(), similarityComp<_Tp> );
				}
				if( !test ){
					it = P[k1].begin();
					while( (*it).sim < C[k1][i].sim && ( it != P[ k1 ].end() ) ){
						it ++;
						//						cout<<" (*it).sim: "<<(*it).sim<<endl;
					}
					if( it != P[k1].end() ){
						P[ k1 ].insert( it, C[k1][i] );
					}else{
						P[ k1 ].push_back( C[k1][i] );
					}
				}
				if(0){
					cout<<"k1:"<<1<<" i:"<<i<<"  C[ k1 ][ i ]:"<<C[ k1 ][ i ].sim<<endl;
					cout<<"P"<<k1;
					for( size_t j=0; j<P[k1].size(); j++ ){
						cout<<" "<<P[k1][ j ].sim;
					}
					cout<<endl;
				}
			}
		}

		//the smallest and the second smallest distance
		minDis = 10e10;
		minDis2 = 10e10;

		//get the least distance and its index
		for( size_t i=0; i<matrixSize; i++ ){
			if(0){
				cout<<"c[i].sim:"<<c[i].sim<<endl;
			}
			if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
				minDis2 = minDis;
				minDis = c[i].sim;
				tempk = i;
			}else if( minDis2 > c[i].sim && minDis != c[i].sim &&  ( c[i].sim != 0 ) & ( I[i] != 0 )  ){
				minDis2 = c[i].sim;
			}

			if(0){
				cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
			}
		}
		if(0){
			cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
		}
		if( minDis > limiteDis ){
//			cout<<endl;
//			cout<<" minDis: "<<minDis<<" limiteDis:"<<limiteDis<<endl;
			break;
		}

		//add new cluster
		vector< size_t > tempV;
		tempV.push_back( k1 );
		tempV.push_back( k2 );
		A.push_back( tempV );
		if(0){
			cout<<"  t3.2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
		}

		//update cluster result vector
		if(1){
			//indicate the way to update cluster result vector
			bool flagk1 = true;
			bool flagk2 = true;
			size_t i1, i2;
			i1 = i2 = 0;

			for( size_t i=0; i<resultVec.size(); i++ ){
				if( find( resultVec[i].begin(), resultVec[i].end(), k1 ) != resultVec[i].end() ){
					flagk1 = false;
					i1 = i;
				}else if( find( resultVec[i].begin(), resultVec[i].end(), k2 ) != resultVec[i].end() ){
					flagk2 = false;
					i2 = i;
				}
			}

			if(0){
				cout<<"i1:"<<i1<<" i2:"<<i2<<endl;
				cout<<"flagk1:"<<flagk1<<" flagk2:"<<flagk2<<endl;
				for( size_t i=0; i<resultVec.size(); i++ ){
					cout<<"# cluster"<<i<<endl;
					for( size_t j=0; j<resultVec[i].size(); j++ ){
						cout<<resultVec[i][j]<<" ";
					}
					cout<<endl;
				}
			}

			if( flagk1 == true &&  flagk2 == true ){
				resultVec.push_back( tempV );
			}else if( flagk1 == true && flagk2 == false ){
				resultVec[ i2 ].push_back( k1 );
			}else if( flagk2 == true && flagk1 == false ){
				resultVec[ i1 ].push_back( k2 );
			}else if( flagk1 == false && flagk2 == false ){
				resultVec[ i1 ].insert( resultVec[i1].end(), resultVec[ i2 ].begin(), resultVec[ i2 ].end() );
				resultVec[ i2 ].clear();
			}

			if(0){
				cout<<endl;
				cout<<"--------------------------------------------------------------------------------"<<endl;
				cout<<"k1:"<<k1<<" k2:"<<k2<<" minDis:"<<minDis<<endl;
				printCluster( resultVec );
				checkLargestInnerDis( inputMatrix, resultVec, limiteDis );
				cout<<"--------------------------------------------------------------------------------"<<endl;
				cout<<endl;
			}
		}
	}
	cerr<<endl;

	if(1){
		bool flag = true;
		for( size_t i=0; i<matrixSize; i++ ){
			flag = true;
			for( size_t j=0; j<resultVec.size(); j++ ){
				if( find( resultVec[j].begin(), resultVec[j].end(), i ) != resultVec[j].end() ){
					flag = false;
				}
			}

			if( flag == true ){
				vector< size_t > temp;
				temp.push_back(i);
				resultVec.push_back( temp );
			}
		}
	}

	if(0){
		cout<<"A:"<<endl;
		for( size_t i=0; i<A.size(); i++ ){
			cout<<A[i][0]<<" "<<A[i][1]<<endl;
		}
	}

	if(1){
		t2 = clock();
		float time = ( (float)t2 - ( float )t1 ) / CLOCKS_PER_SEC ;
		cout<<"completeLink running time:"<<time<<endl;
	}

	resultVec = sortCluster( resultVec );
	if(splitFlag){
		cout<<"cluster number before split:"<<resultVec.size()<<endl;
		resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
		resultVec = sortCluster( resultVec );
		cout<<"result after split:"<<endl;
	}
	printCluster( resultVec );

	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}


/**
 * @brief http://nlp.stanford.edu/IR-book/html/htmledition/time-complexity-of-hac-1.html
 * @param inputMatrix
 * @param limiteDis
 * @return
 * @note	complete-link, not average-link
 */
template<typename _Tp>
vector< vector< size_t > >
fastHierarchical( const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis ){

	clock_t t1, t2;
	t1 = clock();

	cout<<endl<<"---fast Hierarchical algorithm---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for fast Hierarchical clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 0;
	_Tp						minDis2 = 0;

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//priority-queue
	vector< vector< similarity<_Tp> > > P;
	vector<_Tp> p;
	vector< size_t > I;

	if(0){
		cout<<"t1:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back( c );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		for( size_t j=0; j<matrixSize; j++ ){
			C[i][j].sim = inputMatrix[i][j];
			C[i][j].index = j;
		}
		I.push_back( 1 );
	}

	if(0){
		cout<<"t2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//initialize P
	for( size_t i=0; i<matrixSize; i++ ){
		c = C[i];
		sort( c.begin(), c.end(), similarityComp<_Tp> );
		c.erase( c.begin() );
		P.push_back( c );
	}

	vector< vector< size_t > > A;

	if(0){
		cout<<"t3:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//main loop
	for( size_t k=1; (k<matrixSize)  ; k++ ){
		if( minDis	 < limiteDis ){
			if(0){
				cout<<"t"<<k<<":"<< double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			//
			c.clear();
			for( size_t i=0; i<matrixSize; i++ ){
				c.push_back( P[i].front() );
			}

			//find out the index k1 and k2 which have the smallest distance
			size_t k1, k2, tempk;

			//the second smallest distance
			minDis = 10e10;
			minDis2 = 10e10;

			for( size_t i=0; i<matrixSize; i++ ){
				if(0){
					cout<<"c[i].sim:"<<c[i].sim<<endl;
				}
				if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
					minDis2 = minDis;
					minDis = c[i].sim;
					tempk = i;
				}/*else if( (minDis2 > c[i].sim) && ( c[i].sim != 0 ) && ( I[i] != 0 ) ){
				minDis2 = c[i].sim;
			}*/
				if(0){
					cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
				}
			}
			k1 = tempk;
			k2 = P[ k1 ].front().index;
			if(0){
				cout<<"  t3.1:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			if(0){
				cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
			}

			if(0){
				cout<<"P"<<endl;
				for( size_t i=0; i<P.size(); i++ ){
					cout.width(5);
					cout<<i<<" ";
					for( size_t j=0; j<P[i].size(); j++ ){
						cout.width(5);
						cout<<right<<P[i][j].index<<":";
						cout.width(10);
						cout<<left<<P[i][j].sim<<" ";
					}
					cout<<endl;
				}
			}

			if(0){
				cout<<"k1:"<<k1<<" k2:"<<k2<<endl;
			}

			vector< size_t > tempV;
			tempV.push_back( k1 );
			tempV.push_back( k2 );
			A.push_back( tempV );
			if(0){
				cout<<"  t3.2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			if(1){
				bool flagk1 = true;
				bool flagk2 = true;
				size_t i1, i2;

				for( size_t i=0; i<resultVec.size(); i++ ){
					if( find( resultVec[i].begin(), resultVec[i].end(), k1 ) != resultVec[i].end() ){
						//					resultVec[i].push_back( k2 );
						flagk1 = false;
						i1 = i;
					}else if( find( resultVec[i].begin(), resultVec[i].end(), k2 ) != resultVec[i].end() ){
						//					resultVec[i].push_back( k1 );
						flagk2 = false;
						i2 = i;
					}
				}

				if(0){
					cout<<"i1:"<<i1<<" i2:"<<i2<<endl;
					cout<<"flagk1:"<<flagk1<<" flagk2:"<<flagk2<<endl;
					for( size_t i=0; i<resultVec.size(); i++ ){
						cout<<"# cluster"<<i<<endl;
						for( size_t j=0; j<resultVec[i].size(); j++ ){
							cout<<resultVec[i][j]<<" ";
						}
						cout<<endl;
					}
				}

				if( flagk1 == true &&  flagk2 == true ){
					resultVec.push_back( tempV );
				}else if( flagk1 == true && flagk2 == false ){
					resultVec[ i2 ].push_back( k1 );
				}else if( flagk2 == true && flagk1 == false ){
					resultVec[ i1 ].push_back( k2 );
				}else if( flagk1 == false && flagk2 == false ){
					resultVec[ i1 ].insert( resultVec[i1].end(), resultVec[ i2 ].begin(), resultVec[ i2 ].end() );
					resultVec[ i2 ].clear();
				}
			}
			if(0){
				cout<<"  t3.3:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			//index k2 labled removed
			I[k2] = 0;

			//clear P[k1]
			P[k1].clear();

			//update P
			for( size_t i=0; i<matrixSize  ; i++ ){

				if( I[i] == 1 && i != k1 ){

					//delete C[i][k1] in P
					typename vector< similarity<_Tp> >::iterator it ;
					it = find_if( P[i].begin(), P[i].end(), C[i][k1] );
					P[i].erase( it );

					//delete C[i][k2] in P
					it =  find_if( P[i].begin(), P[i].end(), C[i][k2] );
					if( it != P[i].end() ){
						P[i].erase( it );
					}

					bool test = false;
					//update the distance between the ith data and k1, which is the label of cluster that contains k1 and k2
					C[ i ][ k1 ].sim = completelinkDis( i, k1, k2, inputMatrix );
					if( test ){
						P[i].push_back( C[i][k1] );
						sort( P[i].begin(), P[i].end(), similarityComp<_Tp> );
					}
					if( !test ){
						it = P[i].begin();
						//					cout<<" (*it).sim: "<<(*it).sim<<endl;
						while( (*it).sim < C[ i ][ k1 ].sim && ( it != P[ i ].end() ) ){
							it ++;
							//						cout<<" (*it).sim: "<<(*it).sim<<endl;
						}
						if( it != P[ i ].end() ){
							P[ i ].insert( it, C[ i ][ k1 ] );
						}else{
							P[ i ].push_back( C[ i ][ k1 ] );
						}
					}
					if(0){
						cout<<"C[ i ][ k1 ]:"<<C[ i ][ k1 ].sim<<endl;
						cout<<"P"<<i;
						for( size_t j=0; j<P[i].size(); j++ ){
							cout<<" "<<P[i][j].sim;
						}
						cout<<endl;
					}

					//update
					C[ k1 ][ i ].sim = C[ i ][ k1 ].sim ;
					if( test ){
						P[ k1 ].push_back( C[k1][ i ] );
						sort( P[ k1 ].begin(), P[ k1 ].end(), similarityComp<_Tp> );
					}
					if( !test ){
						it = P[k1].begin();
						while( (*it).sim < C[k1][i].sim && ( it != P[ k1 ].end() ) ){
							it ++;
							//						cout<<" (*it).sim: "<<(*it).sim<<endl;
						}
						if( it != P[k1].end() ){
							P[ k1 ].insert( it, C[k1][i] );
						}else{
							P[ k1 ].push_back( C[k1][i] );
						}
					}
					if(0){
						cout<<"C[ k1 ][ i ]:"<<C[ k1 ][ i ].sim<<endl;
						cout<<"P"<<k1;
						for( size_t j=0; j<P[k1].size(); j++ ){
							cout<<" "<<P[k1][ j ].sim;
						}
						cout<<endl;
					}
				}

				for( size_t i=0; i<matrixSize; i++ ){
					if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
						minDis2 = minDis;
						minDis = c[i].sim;
						tempk = i;
					}
				}
			}
			if(0){
				cout<<"  t3.4:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
		}
	}

	if( 1 ){
		bool flag = true;
		for( size_t i=0; i<matrixSize; i++ ){
			flag = true;
			for( size_t j=0; j<resultVec.size(); j++ ){
				if( find( resultVec[j].begin(), resultVec[j].end(), i ) != resultVec[j].end() ){
					flag = false;
				}
			}

			if( flag == true ){
				vector< size_t > temp;
				temp.push_back(i);
				resultVec.push_back( temp );
			}
		}
	}

	if(0){
		cout<<"A:"<<endl;
		for( size_t i=0; i<A.size(); i++ ){
			cout<<A[i][0]<<" "<<A[i][1]<<endl;
		}
	}

	t2 = clock();
	float time = ( (float)t2 - ( float )t1 ) / CLOCKS_PER_SEC ;
	cout<<"fastHierarchical running time:"<<time<<endl;

	resultVec = sortCluster( resultVec );
	printCluster( resultVec );

	return resultVec;
}

/**
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
fastHierarchical_limit_alpha( const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis ){

	clock_t t1, t2;
	t1 = clock();

	cout<<"---fast Hierarchical limit algorithm---"<<endl;
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for fastHierarchical clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t					matrixSize = inputMatrix.size();

	//average link distance between two clusters
	_Tp						averageDis = 0;
	_Tp						minDis = 0;
	_Tp						minDis2 = 0;

	//struct matrix
	vector< vector<similarity<_Tp> > > C;
	vector< similarity<_Tp>  > c;

	//priority-queue
	vector< vector< similarity<_Tp> > > P;
	vector<_Tp> p;
	vector< size_t > I;

	//initialize C
	for( size_t i=0; i<matrixSize; i++ ){
		c.push_back( similarity<_Tp>()  );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		C.push_back( c );
	}
	for( size_t i=0; i<matrixSize; i++ ){
		for( size_t j=0; j<matrixSize; j++ ){
			C[i][j].sim = inputMatrix[i][j];
			C[i][j].index = j;
		}
		I.push_back( 1 );
	}

	if(0){
		cout<<"t2:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//initialize P
	for( size_t i=0; i<matrixSize; i++ ){
		c = C[i];
		sort( c.begin(), c.end(), similarityComp<_Tp> );
		c.erase( c.begin() );
		P.push_back( c );
	}

	vector< vector< size_t > > A;

	if(0){
		cout<<"t3:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	//main loop
	for( size_t k=1; (k<matrixSize)  ; k++ ){
		if( minDis	 < limiteDis ){
			if(0){
				cout<<"t"<<k<<":"<< double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			//
			c.clear();
			for( size_t i=0; i<matrixSize; i++ ){
				c.push_back( P[i].front() );
			}

			//find out the index k1 and k2 which have the smallest distance
			size_t k1, k2, tempk;

			//the second smallest distance
			minDis = 10e10;
			minDis2 = 10e10;

			for( size_t i=0; i<matrixSize; i++ ){
				if(0){
					cout<<"c[i].sim:"<<c[i].sim<<endl;
				}
				if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
					minDis2 = minDis;
					minDis = c[i].sim;
					tempk = i;
				}/*else if( (minDis2 > c[i].sim) && ( c[i].sim != 0 ) && ( I[i] != 0 ) ){
				minDis2 = c[i].sim;
			}*/
				if(0){
					cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
				}
			}
			k1 = tempk;
			k2 = P[ k1 ].front().index;

			if(0){
				cout<<"minDis:"<<minDis<<" minDis2:"<<minDis2<<endl;
			}

			if(0){
				cout<<"P"<<endl;
				for( size_t i=0; i<P.size(); i++ ){
					cout.width(5);
					cout<<i<<" ";
					for( size_t j=0; j<P[i].size(); j++ ){
						cout.width(5);
						cout<<right<<P[i][j].index<<":";
						cout.width(10);
						cout<<left<<P[i][j].sim<<" ";
					}
					cout<<endl;
				}
			}

			if(0){
				cout<<"k1:"<<k1<<" k2:"<<k2<<endl;
			}

			vector< size_t > tempV;
			tempV.push_back( k1 );
			tempV.push_back( k2 );
			A.push_back( tempV );

			if(0){
				cout<<"t"<<1<<":"<< double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}

			if(1){
				bool flagk1 = true;
				bool flagk2 = true;
				size_t i1, i2;

				for( size_t i=0; i<resultVec.size(); i++ ){
					if( find( resultVec[i].begin(), resultVec[i].end(), k1 ) != resultVec[i].end() ){
						//					resultVec[i].push_back( k2 );
						flagk1 = false;
						i1 = i;
					}else if( find( resultVec[i].begin(), resultVec[i].end(), k2 ) != resultVec[i].end() ){
						//					resultVec[i].push_back( k1 );
						flagk2 = false;
						i2 = i;
					}
				}

				if(0){
					cout<<"i1:"<<i1<<" i2:"<<i2<<endl;
					cout<<"flagk1:"<<flagk1<<" flagk2:"<<flagk2<<endl;
					for( size_t i=0; i<resultVec.size(); i++ ){
						cout<<"# cluster"<<i<<endl;
						for( size_t j=0; j<resultVec[i].size(); j++ ){
							cout<<resultVec[i][j]<<" ";
						}
						cout<<endl;
					}
				}

				if( flagk1 == true &&  flagk2 == true ){
					resultVec.push_back( tempV );
				}else if( flagk1 == true && flagk2 == false ){
					resultVec[ i2 ].push_back( k1 );
				}else if( flagk2 == true && flagk1 == false ){
					resultVec[ i1 ].push_back( k2 );
				}else if( flagk1 == false && flagk2 == false ){
					resultVec[ i1 ].insert( resultVec[i1].end(), resultVec[ i2 ].begin(), resultVec[ i2 ].end() );
					resultVec[ i2 ].clear();
				}
			}
			if(0){
				cout<<"t"<<1<<":"<< double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
			}
			//index k2 labled removed
			I[k2] = 0;

			//clear P[k1]
			P[k1].clear();

			//update P
			for( size_t i=0; i<matrixSize  ; i++ ){

				if( I[i] == 1 && i != k1 ){

					//delete C[i][k1] in P
					typename vector< similarity<_Tp> >::iterator it ;
					it = find_if( P[i].begin(), P[i].end(), C[i][k1] );
					P[i].erase( it );

					//delete C[i][k2] in P
					it =  find_if( P[i].begin(), P[i].end(), C[i][k2] );
					if( it != P[i].end() ){
						P[i].erase( it );
					}

					bool test = false;
					//update the distance between the ith data and k1, which is the label of cluster that contains k1 and k2
					C[ i ][ k1 ].sim = completelinkDis( i, k1, k2, inputMatrix );
					if( test ){
						P[i].push_back( C[i][k1] );
						sort( P[i].begin(), P[i].end(), similarityComp<_Tp> );
					}
					if( !test ){
						it = P[i].begin();
						//					cout<<" (*it).sim: "<<(*it).sim<<endl;
						while( (*it).sim < C[ i ][ k1 ].sim && ( it != P[ i ].end() ) ){
							it ++;
							//						cout<<" (*it).sim: "<<(*it).sim<<endl;
						}
						if( it != P[ i ].end() ){
							P[ i ].insert( it, C[ i ][ k1 ] );
						}else{
							P[ i ].push_back( C[ i ][ k1 ] );
						}
					}
					if(0){
						cout<<"C[ i ][ k1 ]:"<<C[ i ][ k1 ].sim<<endl;
						cout<<"P"<<i;
						for( size_t j=0; j<P[i].size(); j++ ){
							cout<<" "<<P[i][j].sim;
						}
						cout<<endl;
					}

					//update
					C[ k1 ][ i ].sim = C[ i ][ k1 ].sim ;
					if( test ){
						P[ k1 ].push_back( C[k1][ i ] );
						sort( P[ k1 ].begin(), P[ k1 ].end(), similarityComp<_Tp> );
					}
					if( !test ){
						it = P[k1].begin();
						while( (*it).sim < C[k1][i].sim && ( it != P[ k1 ].end() ) ){
							it ++;
							//						cout<<" (*it).sim: "<<(*it).sim<<endl;
						}
						if( it != P[k1].end() ){
							P[ k1 ].insert( it, C[k1][i] );
						}else{
							P[ k1 ].push_back( C[k1][i] );
						}
					}
					if(0){
						cout<<"C[ k1 ][ i ]:"<<C[ k1 ][ i ].sim<<endl;
						cout<<"P"<<k1;
						for( size_t j=0; j<P[k1].size(); j++ ){
							cout<<" "<<P[k1][ j ].sim;
						}
						cout<<endl;
					}
				}

				for( size_t i=0; i<matrixSize; i++ ){
					if( ( minDis > c[i].sim ) && ( c[i].sim != 0 ) & ( I[i] != 0 ) ){
						minDis2 = minDis;
						minDis = c[i].sim;
						tempk = i;
					}
				}
			}
		}
	}

	if( 1 ){
		bool flag = true;
		for( size_t i=0; i<matrixSize; i++ ){
			flag = true;
			for( size_t j=0; j<resultVec.size(); j++ ){
				if( find( resultVec[j].begin(), resultVec[j].end(), i ) != resultVec[j].end() ){
					flag = false;
				}
			}
			if( flag == true ){
				vector< size_t > temp;
				temp.push_back(i);
				resultVec.push_back( temp );
			}
		}
	}

	if(0){
		cout<<"A:"<<endl;
		for( size_t i=0; i<A.size(); i++ ){
			cout<<A[i][0]<<" "<<A[i][1]<<endl;
		}
	}

	t2 = clock();
	float time = ( (float)t2 - ( float )t1 ) / CLOCKS_PER_SEC ;
	cout<<"fastHierarchical_limit_alpha running time:"<<time<<endl;

	resultVec = sortCluster( resultVec );
	printCluster( resultVec );

	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * @brief 	compute the triangle area constructed by edge1, edge2, edge3
 * @param 	edge1 Edge 1
 * @param 	edge2 Edge 2
 * @param 	edge3 Edge 3
 * @return 	triangel area
 */
template<typename _Tp >
_Tp
triangleArea(_Tp edge1, _Tp edge2, _Tp edge3 ){
	if(edge1< 0 || edge2< 0 || edge3<0 ){
		cerr<<"error!!! edge length should be positive value!!! "<<endl;
	}
	_Tp area = 0;
	_Tp s = (edge1 + edge2 + edge3 )/2;
	area = sqrt(s* (s- edge1)* (s-edge2 )* (s-edge3 ) );
	return area;
}

/**
 * @brief 	compute the tetrahedron volume using Cayley-Menger determinants
 */
template<typename _Tp>
_Tp
tetrahedralVolumn2(const _Tp a, const _Tp b, const _Tp c, const _Tp p,
		const _Tp q, const _Tp r){
	if(a< 0 || b< 0 || c<0 || p <0 || q< 0 || r<0 ){
		cerr<<"error!!! edge length should be positive value!!! "<<endl;
	}
	_Tp volume = 0;
	_Tp p1 = (a*p) * (a*p) * (-a*a + b*b + c*c - p*p + q*q + r*r);
	_Tp p2 = (b*q) * (b*q) * (a*a  - b*b + c*c + p*p - q*q + r*r);
	_Tp p3 = (c*r) * (c*r) * (a*a  + b*b - c*c + p*p + q*q - r*r);
	_Tp P = (a*b*r)*(a*b*r) + (a*c*q)*(a*c*q) + (b*c*p)*(b*c*p) + (p*q*r)*(p*q*r);
	volume = sqrt(p1 + p2 + p3 - P)/12;
	return volume;
}

/**
 * @brief 	compute the tetrahedron volume using Cayley-Menger determinants
 */
template<typename _Tp>
_Tp
tetrahedralVolumn(const _Tp a, const _Tp b, const _Tp c, const _Tp d,
		const _Tp e, const _Tp f,  const bool squareDis = false){
	if(a< 0 || b< 0 || c<0 || d <0 || e< 0 || f<0 ){
		cerr<<"error!!! edge length should be positive value!!! "<<endl;
	}
	if( a == 0 || b == 0 || c == 0 || d == 0 || e == 0 || f == 0 ){
		return _Tp(0);
	}

	if( 0 ){
		cout<<"a:"<<a<<" b:"<<b<<" c:"<<c<<" d:"<<d<<" e:"<<e<<" f:"<<f<<endl;
	}
/*
	_Tp volume = 0;
	_Tp a2 = a*a;
	_Tp b2 = b*b;
	_Tp c2 = c*c;
	_Tp d2 = d*d;
	_Tp e2 = e*e;
	_Tp f2 = f*f;
	volume = ( a2 + b2 + c2 + d2 + e2 + f2 ) * ( a2*f2 + b2*e2 + c2*d2 )
			- 2 * a2 * f2 * ( a2 + f2 )
			- 2 * b2 * e2 * ( b2 + e2)
			- 2 * c2 * d2 * ( c2 + d2 )
			- ( a2 + f2 )*(b2 + e2 )*(c2 + d2 )/2
			+ ( a2 - f2 )*(b2 - e2 )*(c2 - d2 )/2;
	if( 1 ){
		cout<<"( a2 + b2 + c2 + d2 + e2 + f2 ) * ( a2*f2 + b2*e2 + c2*d2 ): "<<( a2 + b2 + c2 + d2 + e2 + f2 ) * ( a2*f2 + b2*e2 + c2*d2 )<<endl;
		cout<<"2 * a2 * f2 * ( a2 + f2 ):"<<2 * a2 * f2 * ( a2 + f2 )<<endl;
		cout<<"2 * b2 * e2 * ( b2 + e2):"<<2 * b2 * e2 * ( b2 + e2)<<endl;
		cout<<"2 * c2 * d2 * ( c2 + d2 ):"<<2 * c2 * d2 * ( c2 + d2 )<<endl;
		cout<<"( a2 + f2 )*(b2 + e2 )*(c2 + d2 )/2:"<<( a2 + f2 )*(b2 + e2 )*(c2 + d2 )/2<<endl;
		cout<<"( a2 - f2 )*(b2 - e2 )*(c2 - d2 )/2:"<< ( a2 - f2 )*(b2 - e2 )*(c2 - d2 )/2<<endl;
		cout<<"volume:"<<volume<<endl;
	}
	volume /= 144;
	if( volume < 0 ){
		volume = 0;
	}else{
		volume = sqrt(volume );
	}
	return volume;
	*/

	_Tp volume = 0;
	_Tp A, B, C, D, E, F;
	if( !squareDis ){
	A = a*a;
	B = b*b;
	C = c*c;
	D = d*d;
	E = e*e;
	F = f*f;
	}else{
		A = a;
		B = b;
		C = c;
		D = d;
		E = e;
		F = f;
	}

	if(0){
		cout<<"-A*B*C:"<<-A*B*C<<" (D-B)*(F-C)*(E-A):"<<(D-B)*(F-C)*(E-A)<<" (D-A)*(F-B)*(E-C):"<<(D-A)*(F-B)*(E-C);
		cout<<" B*(E-A)*(E-C):"<<B*(E-A)*(E-C)<<" C*(D-A)*(D-B):"<<C*(D-A)*(D-B)<<" A*(F-B)*(F-C):"<<A*(F-B)*(F-C)<<endl;
	}

	_Tp col1 = -A*B*C + (D-B)*(F-C)*(E-A) + (D-A)*(F-B)*(E-C) + B*(E-A)*(E-C) + C*(D-A)*(D-B) + A*(F-B)*(F-C);
	_Tp col2 = A*B*C + C*(D-B)*(F-C) + B*(F-B)*(E-C) + B*C*( E-C ) + C*B*(D-B) - A*(F-B)*(F-C);
	_Tp col3 = -A*B*C - A*C*(D-A) - A*C*(F-C) + B*(E-A)*(E-C) - C*(D-A)*(E-C) - A*(E-A)*(F-C);
	_Tp col4 = A*B*C + A*(D-A)*(F-B) + B*(E-A)*(D-B) - C*(D-A)*(D-B) + A*B*(F-B) + A*B*(E-A);

	double cold1 = -A*B*C + (D-B)*(F-C)*(E-A) + (D-A)*(F-B)*(E-C) + B*(E-A)*(E-C) + C*(D-A)*(D-B) + A*(F-B)*(F-C);

	if( 0 ){
		cout<<"cold1:"<< -A*B*C + (D-B)*(F-C)*(E-A) + (D-A)*(F-B)*(E-C)+ B*(E-A)*(E-C)  + A*(F-B)*(F-C)<<endl;
	}

	volume = -col1 + col2 - col3 + col4;
	if( 0 ){
		cout<<"A:"<<A<<" B:"<<B<<" C:"<<C<<" D:"<<D<<" E:"<<E<<" F:"<<F<<endl;
		cout<<"col1:"<<col1<<endl;
		cout<<"col2:"<<col2<<endl;
		cout<<"col3:"<<col3<<endl;
		cout<<"col4:"<<col4<<endl;
		cout<<"volume:"<<volume<<endl;
	}
	volume /= 288;
	if( volume < 0 ){
		volume = 0;
	}else {
		volume = sqrt(volume );
	}

	//check if the value of volume is too small compared with the edge length, if true, set it to 0
	if(1){
		vector<_Tp> edgeV;
		edgeV.push_back(a);
		edgeV.push_back(b);
		edgeV.push_back(c);
		edgeV.push_back(d);
		edgeV.push_back(e);
		edgeV.push_back(f);
		_Tp maxEdge = *max_element( edgeV.begin(), edgeV.end() );
		if( volume < ( maxEdge*maxEdge*maxEdge/10000 ) ){
			volume = 0;
		}
	}
	return volume;

}

/**
 * @brief		The iterative clustering algorithm proposed by wanglincong and xushutan
 * @param	inputDis The input distance matrix deposited in one dimension
 * @param	limiteDis The distance criteria for the algorithm
 * @param	clusterResult	The clustering result
 * @return	true if success false if failed
 * @note		The 'input' contains the pairwise distance of all structures in one dimension vector,
 * 					where $ i*(i+1 )/2 + j$th (i>j) item corresponding to the [i][j]th item in two dimension vector
 *
 *The correspondance between the one dimension 'inputDis' and two dimension distance matrix is as below
 * 	|    0                 0                  0                  0                0     |
 * 	| inputDis[0]     0                  0                  0                0     |
 * 	| inputDis[1] inputDis[2]       0                  0                0     |
 * 	| inputDis[3] inputDis[4] inputDis[5]         0                0     |
 * 	| inputDis[6] inputDis[7] inputDis[8]  inputDis[9]        0     |
 */
template<typename _Tp>
bool
iterative( const vector<_Tp> inputDis, const _Tp limiteDis,
				vector< vector<size_t> >& clusterResult  ){
	cout<<"---iterative---"<<endl;

	size_t structNum = ( sqrt( 8 * inputDis.size() + 1 ) - 1 )/2 + 1;

	if(  inputDis.size()  != ( ( structNum - 1) * ( structNum - 1 +1 ) ) / 2 ){
		cerr<<"---error in k_means, please check if  structNum is correct---"<<endl;
		return false;
	}

	vector< vector<_Tp> > inputMatrix(structNum, vector<_Tp>( structNum ));

	// initialize input distance matrix
	inputMatrix[0][0] = 0;
	for( size_t i=1; i<structNum; i++ ){
		for( size_t j=0; j<i; j++ ){
			inputMatrix[i][j] = inputMatrix[j][i] = inputDis[ i*(i-1)/2 + j  ];
		}
		inputMatrix[i][i] = 0;
	}

	clusterResult.clear();
	clusterResult = iterative( inputMatrix, limiteDis );

	return true;
}

/**
 * @brief compute according to seed points recursively
 * @param inputMatrix 	The pairwise distance matrix
 * @param limiteDis		The distance criteria
 * @param inputSeqVec	The data index vector for clustering
 *
 * @note This function cluster data according to seedpoints recursively.
 * The seedpoints is defined as the endpoints that construct the
 * largest tetrahedral
 */
template<typename _Tp>
vector<vector<size_t> >
tetrahedral( const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis,
		const vector<size_t>& inputSeqVec ){
	if(0){
		cout<<endl;
		cout<<"---tetrahedral---"<<endl;
	}
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for tetrahedral clustering!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//distance
	_Tp		dis = 0;
	_Tp 	maxDis = 0;

	//area
	_Tp		area = 0;
	_Tp		maxArea = 0;

	//volume
	_Tp		volume = 0;
	_Tp 	maxVolume = 0;

	//seed points
	size_t p0, p1, p2, p3;
	p0= p1= p2= p3= 0;

	//temp cluster corresponding to p0, p1, p2, p3
	vector<size_t > cluster0, cluster1, cluster2, cluster3;

	//
	vector<vector<size_t> > resultVec;
	resultVec.clear();

	if(inputSeqVec.size() == 1 ){
		resultVec.push_back(vector<size_t >(1, inputSeqVec[0]) );
		if(0){
			cout<<"resultVec:"<<inputSeqVec[0]<<"  done!!!"<<endl;
		}
		return resultVec;
	}

	if(0){
//		cout<<"inputMatrix.size():"<<inputMatrix.size()<<endl;
//		cout<<"inputSeqVec.size():"<<inputSeqVec.size()<<endl;
		cout<<"inputSeqVec: ";
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			cout<<inputSeqVec[i]<<" ";
		}
		cout<<endl;
	}
	//-----------------------------------------two seedpoint
	for(size_t i=1; i<inputSeqVec.size(); i++ ){
		for(size_t j=0; j<i; j++ ){
			dis = inputMatrix[inputSeqVec[i] ][ inputSeqVec[j]];
			if(dis > maxDis ){
				maxDis = dis;

				p0 = inputSeqVec[i ];
				p1 = inputSeqVec[j ];
			}
		}
	}

	if(0){
		cout<<"p0:"<<p0<<" p1:"<<p1<<" ="<<maxDis<<" limiteDis:"<<limiteDis<<endl;
	}

	// all in one cluster
	if(maxDis < limiteDis ){
		resultVec.push_back(inputSeqVec );
		if(0){
			cout<<"resultVec:"<<endl;
			for(size_t i=0; i<resultVec.size(); i++){
				cout<<i<<": ";
				for(size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<" done!!!"<<endl;
			}
		}
		return 	resultVec;
	}
	if(0){
		cout<<"1"<<endl;
	}

	//-----------------------------------------three seedpoints
	if( inputSeqVec.size() == 2 ){
		resultVec.push_back(vector<size_t>(1, inputSeqVec[0]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[1]) );
//		cout<<"resultVec:"<<inputSeqVec[0]<<"  done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[1]<<"  done!!!"<<endl;
		return resultVec;
	}
	//search for p2 that construct the largest triangle with p0,p1
	for(size_t i=0; i<inputSeqVec.size(); i++ ){
		area = triangleArea(inputMatrix[p0][p1], inputMatrix[p0][inputSeqVec[i] ],
				inputMatrix[p1][inputSeqVec[i] ]);
		if(0){
			cout<<inputMatrix[p0][p1]<<" "<<inputMatrix[p0][inputSeqVec[i] ]<<" "
					<< inputMatrix[p1][inputSeqVec[i] ]<<" area:"<<area<<endl;
		}
		if(maxArea < area ){
			maxArea = area;
			p2 = inputSeqVec[i ];
		}
	}
	if(1){
		cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" maxArea="
				<<maxArea<<" limiteDis:"<<limiteDis<<endl;
		cout<<"inputMatrix[p0 ][p2 ]:"<<inputMatrix[p0 ][p2 ]<<
				"  inputMatrix[p1 ][p2 ]:"<<inputMatrix[p1 ][p2 ]<<endl;
	}

	//seperate the data into two clusters, then cluster the two clusters recursively
	if(inputMatrix[p0 ][p2 ]  < limiteDis || inputMatrix[p1 ][p2 ] < limiteDis){
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i] ] < inputMatrix[p1 ][inputSeqVec[i] ] ){
				cluster0.push_back(inputSeqVec[i] );
			}else{
				cluster1.push_back(inputSeqVec[i] );
			}
		}

		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral(inputMatrix, limiteDis, cluster0 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );
//		resultVec.insert( resultVec.end(), cluster0);

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );
//		resultVec.insert( resultVec.end(), cluster1);

		if(1){
			cout<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
		}

		return resultVec;
	}
	if(1){
		cout<<"3 step"<<endl;
	}

	//---------------------------------------------four seedpoints
	if(inputSeqVec.size() == 3 ){
		resultVec.push_back(vector<size_t>(1, inputSeqVec[0]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[1]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[2]));
//		cout<<"resultVec:"<<inputSeqVec[0]<<" done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[1]<<" done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[2]<<" done!!!"<<endl;
		return resultVec;
	}
	//search for p3 that construct the largest tetrahedral with p0, p1, p2
	for(size_t i=0; i< inputSeqVec.size();  i++){
		volume = tetrahedralVolumn(inputMatrix[p0 ][p1 ], inputMatrix[p0][p2], inputMatrix[p1][p2],
				inputMatrix[p0][inputSeqVec[i]], inputMatrix[p1][inputSeqVec[i]],
				inputMatrix[p2][inputSeqVec[i]] );
		if(volume > maxVolume ){
			maxVolume = volume;
			p3 = inputSeqVec[i];
		}
	}
	if(0){
		cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" p3:"<<p3<<" maxVolume="
				<<maxVolume<<" limiteDis:"<<limiteDis<<endl;
		cout<<"inputMatrix[p0 ][p3 ]:"<<inputMatrix[p0 ][p3 ]<<
				"  inputMatrix[p1 ][p3 ]:"<<inputMatrix[p1 ][p3 ]<<
				"  inputMatrix[p2 ][p3 ]:"<<inputMatrix[p2 ][p3 ]<<endl;
	}

	//at least one edge is less than limiteDis
	if(inputMatrix[p0][p3] < limiteDis || inputMatrix[p1][p3] < limiteDis || inputMatrix[p2][p3]< limiteDis ){
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster0.push_back(inputSeqVec[i] );
			}else if( inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster1.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] ){
				cluster2.push_back(inputSeqVec[i] );
			}
		}

		if(0){
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
		}

//		resultVec.push_back(cluster0 );
//		resultVec.push_back(cluster1 );
//		resultVec.push_back(cluster2 );

		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral(inputMatrix, limiteDis, cluster0 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster2 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

	}else{
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster0.push_back(inputSeqVec[i] );
			}else if( inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster1.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster2.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			         inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster3.push_back(inputSeqVec[i] );
			}
		}

		if(0){
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster3:";
			for(size_t i=0; i<cluster3.size(); i++ ){
				cout<<cluster3[i]<<" ";
			}
			cout<<endl;
		}

		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral(inputMatrix, limiteDis, cluster0 );
		if(0){
			cout<<"tempVec:";
			for(size_t i=0; i<tempVec.size(); i++){
				for(size_t j=0; j<tempVec[i].size(); j++ )
					cout<<tempVec[i][j]<<" ";
				cout<<endl;
			}
			cout<<endl;
		}
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster2 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral(inputMatrix, limiteDis, cluster3 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );
	}
}

/**
 * @brief		compute the center for input group based on distance inputMatrix
 * @param inputMatrix
 * @param group
 * @return
 */
template<typename _Tp>
size_t
computeCenter( const vector<vector<_Tp> >& inputMatrix, const vector<size_t>& group ){
	//error
	if(inputMatrix.empty() || group.empty() ){
		cerr<<"---no input data for computeCenter!!!---"<<endl;
	}

	size_t center;
	_Tp minDis = 10e10;

	for( size_t i=0; i<group.size(); i++ ){
		_Tp allDis = 0;
		for( size_t j=0; j<group.size(); j++ ){
			allDis += inputMatrix [ group[ i ] ][ group[ j ] ];
		}
		if( allDis < minDis ){
			minDis = allDis;
			center = group[ i ];
		}
	}

	return center;
}

/**
 * @brief		redistribute the data in initialCluster
 * @param inputMatrix
 * @param initialCluster
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
resetCluster( const vector< vector<_Tp> >& inputMatrix, const vector<vector<size_t> >& initialCluster ){
	if(0){
		cout<<"---resetCluster---"<<endl;
	}

	size_t clusterNum = initialCluster.size();
	vector< size_t > centerVec( clusterNum, 0 );
	for( size_t i=0; i<clusterNum; i++ ){
		centerVec[i] = computeCenter( inputMatrix, initialCluster[i] );
	}

	vector< vector<size_t> > resultVec( clusterNum );

	//all data in initialCluster
	for( size_t i=0; i<clusterNum; i++ ){
		for( size_t j=0; j<initialCluster[ i ].size(); j++ ){
			vector< _Tp > disToClusterVec;

			//contain the distance of a data to all the centers
			for(size_t k=0; k<clusterNum; k++ ){
				disToClusterVec.push_back( inputMatrix[ initialCluster[i][j] ][  centerVec[k] ] );
			}

			//identify which center the data is closet to
			size_t index = size_t( min_element( disToClusterVec.begin(), disToClusterVec.end() ) - disToClusterVec.begin() );

			//contain each data to it`s closet center cluster
			resultVec[ index ].push_back( initialCluster[i][j] );
		}
	}

	return resultVec;
}

/**
 * @
 * @param inputMatrix
 * @param initialVec
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
rearrange( const vector<vector<_Tp> >& inputMatrix, const vector<vector<size_t> >& initialVec, size_t loopNum=10 ){
	if(0){
		cout<<"---rearrange---"<<endl;
	}

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for rearrange clustering!!!---"<<endl;
	}

	vector<vector<size_t> > resultVec( initialVec );
	size_t numOfCluster = initialVec.size();
	vector< size_t > centerOldVec( numOfCluster, 0 );
	vector< size_t > centerNewVec( numOfCluster, 0 );

	bool flag = true;
	while( loopNum-- && flag ){
		if(0){
			cout<<"loopNum:"<<loopNum<<endl;
		}
		flag = false;
		for( size_t i=0; i<numOfCluster; i++ ){
			//compute the center for each cluster, this likes the k-means
			centerOldVec[i] = computeCenter( inputMatrix, resultVec[i] );
			if(0){
				cout<<"numOfCluster:  "<<numOfCluster<<"  i:"<<i<<endl;
			}
		}
		if(0){
			cout<<"1"<<endl;
		}

		//reset the cluster distribution of input data in resultVec
		resultVec = resetCluster( inputMatrix, resultVec );

		if(0){
			cout<<"2"<<endl;
		}

		//compute the new center and identify whether it is the same as the centerOldVec
		for( size_t i=0; i<numOfCluster; i++ ){
			centerNewVec[i] = computeCenter( inputMatrix, resultVec[i] );
		}

		for( size_t i=0; i<numOfCluster; i++ ){
			if( centerNewVec[i] != centerOldVec[i] ){
				flag = true;
			}
			if(0){
				cout<<"numOfCluster:  "<<numOfCluster<<"  i:"<<i<<endl;
			}
		}
	}
	return resultVec;
}

/**
 * @note
 * @param inputMatrix
 * @param limiteDis
 * @param inputSeqVec
 * @param squareDis
 * @return
 */
template<typename _Tp>
vector<vector<size_t> >
tetrahedral_means( const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis,
		const vector<size_t>& inputSeqVec, const bool squareDis = false ){
	if(0){
		cout<<endl;
		cout<<"---tetrahedral_means---"<<endl;
	}
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for tetrahedral_means clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}

	//distance
	_Tp		dis = 0;
	_Tp 	maxDis = 0;

	//area
	_Tp		area = 0;
	_Tp		maxArea = 0;

	//volume
	_Tp		volume = 0;
	_Tp 	maxVolume = 0;

	//seed points
	size_t p0, p1, p2, p3;
	p0= p1= p2= p3= 0;

	//temp cluster corresponding to p0, p1, p2, p3
	vector<size_t > cluster0, cluster1, cluster2, cluster3;

	//
	vector<vector<size_t> > resultVec;
	resultVec.clear();

	if(inputSeqVec.size() == 1 ){
		resultVec.push_back(vector<size_t >(1, inputSeqVec[0]) );
		if(0){
			cout<<"resultVec:"<<inputSeqVec[0]<<"  done!!!"<<endl;
		}
		return resultVec;
	}

	if(0){
//		cout<<"inputMatrix.size():"<<inputMatrix.size()<<endl;
//		cout<<"inputSeqVec.size():"<<inputSeqVec.size()<<endl;
		cout<<"inputSeqVec: ";
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			cout<<inputSeqVec[i]<<" ";
		}
		cout<<endl;
	}

	//-----------------------------------------two seedpoint
	for(size_t i=1; i<inputSeqVec.size(); i++ ){
		for(size_t j=0; j<i; j++ ){
			dis = inputMatrix[inputSeqVec[i] ][ inputSeqVec[j]];
			if(dis > maxDis ){
				maxDis = dis;

				p0 = inputSeqVec[i ];
				p1 = inputSeqVec[j ];
			}
		}
	}

	if(0){
		cout<<"-----------------------------------------two seedpoint"<<endl;
		cout<<"p0:"<<p0<<" p1:"<<p1<<"  maxDis:"<<maxDis<<" limiteDis:"<<limiteDis<<endl;
	}

	// all in one cluster
	if(maxDis < limiteDis ){
		resultVec.push_back(inputSeqVec );
		if(0){
			cout<<"resultVec:"<<endl;
			for(size_t i=0; i<resultVec.size(); i++){
				cout<<i<<": ";
				for(size_t j=0; j<resultVec[i].size(); j++ ){
					cout<<resultVec[i][j]<<" ";
				}
				cout<<" done!!!"<<endl;
			}
		}
		return 	resultVec;
	}
	if(0){
		cout<<"1"<<endl;
	}

	//-----------------------------------------three seedpoints
	if( inputSeqVec.size() == 2 ){
		resultVec.push_back(vector<size_t>(1, inputSeqVec[0]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[1]) );
//		cout<<"resultVec:"<<inputSeqVec[0]<<"  done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[1]<<"  done!!!"<<endl;
		return resultVec;
	}
	//search for p2 that construct the largest triangle with p0,p1
	for(size_t i=0; i<inputSeqVec.size(); i++ ){
		area = triangleArea(inputMatrix[p0][p1], inputMatrix[p0][inputSeqVec[i] ],
				inputMatrix[p1][inputSeqVec[i] ]);
		if(0){
			cout<<inputMatrix[p0][p1]<<" "<<inputMatrix[p0][inputSeqVec[i] ]<<" "
					<< inputMatrix[p1][inputSeqVec[i] ]<<" area:"<<area<<endl;
		}
		if(maxArea < area ){
			maxArea = area;
			p2 = inputSeqVec[i ];
		}
	}
	if( 0 ){
		cout<<"-----------------------------------------three seedpoints"<<endl;
		cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" maxArea="
				<<maxArea<<" limiteDis:"<<limiteDis<<endl;
		cout<<"inputMatrix[p0 ][p2 ]:"<<inputMatrix[p0 ][p2 ]<<
				"  inputMatrix[p1 ][p2 ]:"<<inputMatrix[p1 ][p2 ]<<endl;
	}

	//seperate the data into two clusters, then cluster the two clusters recursively
	if(inputMatrix[p0 ][p2 ]  < limiteDis || inputMatrix[p1 ][p2 ] < limiteDis){
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i] ] < inputMatrix[p1 ][inputSeqVec[i] ] ){
				cluster0.push_back(inputSeqVec[i] );
			}else{
				cluster1.push_back(inputSeqVec[i] );
			}
		}

		if(0){
			if( cluster0.empty() || cluster1.empty() ){
				cout<<"cluster0.empty()  cluster1.empty() "<<endl;
			}
		}
		if(0){
			cout<<endl;
			cout<<"---before arrange---"<<endl;
			cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" maxArea="
					<<maxArea<<" limiteDis:"<<limiteDis<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
		}
		//rearrange according to distribution
		vector<vector< size_t > > cluster0and1;
		cluster0and1.push_back( cluster0 );
		cluster0and1.push_back( cluster1 );
		cluster0and1 = rearrange( inputMatrix, cluster0and1 );
		cluster0 = cluster0and1[0];
		cluster1 = cluster0and1[1];

		if(0){
			cout<<endl;
			cout<<"---after arrange---"<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
		}
		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster0 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		if(0){
			cout<<endl;
			cout<<"---after arrange---"<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
		}

		return resultVec;
	}
	if(0){
		cout<<"3"<<endl;
	}

	//---------------------------------------------four seedpoints
	if(inputSeqVec.size() == 3 ){
		resultVec.push_back(vector<size_t>(1, inputSeqVec[0]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[1]) );
		resultVec.push_back(vector<size_t>(1, inputSeqVec[2]));
//		cout<<"resultVec:"<<inputSeqVec[0]<<" done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[1]<<" done!!!"<<endl;
//		cout<<"resultVec:"<<inputSeqVec[2]<<" done!!!"<<endl;
		return resultVec;
	}
	//search for p3 that construct the largest tetrahedral with p0, p1, p2
	for(size_t i=0; i< inputSeqVec.size();  i++){
			volume = tetrahedralVolumn( inputMatrix[p0][p1], inputMatrix[p0][p2], inputMatrix[p0][ inputSeqVec[i] ], inputMatrix[ p1 ][p2],
					inputMatrix[ p1 ][ inputSeqVec[i] ], inputMatrix[ p2 ][ inputSeqVec[i] ], squareDis);

		if(volume > maxVolume ){
			maxVolume = volume;
			p3 = inputSeqVec[i];
		}
		if(0){
			cout<<i<<": "<<volume<<endl;
		}
	}
	if(0){
		cout<<"---------------------------------------------four seedpoints"<<endl;
		cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" p3:"<<p3<<" maxVolume="
				<<maxVolume<<" limiteDis:"<<limiteDis<<endl;
	}
	if( 0 ){
		cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" p3:"<<p3<<" maxVolume="
				<<maxVolume<<" limiteDis:"<<limiteDis<<endl;
		cout<<"inputMatrix[p0 ][p3 ]:"<<inputMatrix[p0 ][p3 ]<<
				"  inputMatrix[p1 ][p3 ]:"<<inputMatrix[p1 ][p3 ]<<
				"  inputMatrix[p2 ][p3 ]:"<<inputMatrix[p2 ][p3 ]<<endl;
	}

	//at least one edge is less than limiteDis, all there is no tetrahedron, which means all data are on a plane
	if(inputMatrix[p0][p3] < limiteDis || inputMatrix[p1][p3] < limiteDis || inputMatrix[p2][p3]< limiteDis ||
			0 == maxVolume ){
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster0.push_back(inputSeqVec[i] );
			}else if( inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster1.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] ){
				cluster2.push_back(inputSeqVec[i] );
			}
		}

		if(0){
			cout<<endl;
			cout<<"---before arrange---"<<endl;
			cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" p3:"<<p3<<" maxVolume="
					<<maxVolume<<" limiteDis:"<<limiteDis<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
		}

//		resultVec.push_back(cluster0 );
//		resultVec.push_back(cluster1 );
//		resultVec.push_back(cluster2 );

		if(0){
			if( cluster0.empty() || cluster1.empty() || cluster2.empty() ){
				cout<<"cluster0.empty()  cluster1.empty() "<<endl;
			}
		}

		//rearrange according to distribution
		vector<vector< size_t > > cluster012;
		cluster012.push_back( cluster0 );
		cluster012.push_back( cluster1 );
		cluster012.push_back( cluster2 );
		cluster012 = rearrange( inputMatrix, cluster012 );
		cluster0 = cluster012[0];
		cluster1 = cluster012[1];
		cluster2 = cluster012[2];

		if(0){
			cout<<endl;
			cout<<"---after rearrage---"<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
		}

		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster0 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster2 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

	}else{
		for(size_t i=0; i<inputSeqVec.size(); i++ ){
			if( inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]] &&
			    inputMatrix[p0][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster0.push_back(inputSeqVec[i] );
			}else if( inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]] &&
			          inputMatrix[p1][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster1.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			         inputMatrix[p2][inputSeqVec[i]]< inputMatrix[p3][inputSeqVec[i]]){
				cluster2.push_back(inputSeqVec[i] );
			}else if(inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p0][inputSeqVec[i]] &&
			         inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p1][inputSeqVec[i]] &&
			         inputMatrix[p3][inputSeqVec[i]]< inputMatrix[p2][inputSeqVec[i]]){
				cluster3.push_back(inputSeqVec[i] );
			}
		}

		if(0){
			cout<<endl;
			cout<<"---before arrange---"<<endl;
			cout<<"p0:"<<p0<<" p1:"<<p1<<" p2:"<<p2<<" p3:"<<p3<<" maxVolume="
					<<maxVolume<<" limiteDis:"<<limiteDis<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster3:";
			for(size_t i=0; i<cluster3.size(); i++ ){
				cout<<cluster3[i]<<" ";
			}
			cout<<endl;
		}

		//rearrange according to distribution
		vector<vector< size_t > > cluster0123;
		cluster0123.push_back( cluster0 );
		cluster0123.push_back( cluster1 );
		cluster0123.push_back( cluster2 );
		cluster0123.push_back( cluster3 );
		cluster0123 = rearrange( inputMatrix, cluster0123 );
		cluster0 = cluster0123[0];
		cluster1 = cluster0123[1];
		cluster2 = cluster0123[2];
		cluster3 = cluster0123[3];

		if(0){
			cout<<endl;
			cout<<"---after arrange---"<<endl;
			cout<<"cluster0:";
			for(size_t i=0; i<cluster0.size(); i++ ){
				cout<<cluster0[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster1:";
			for(size_t i=0; i<cluster1.size(); i++ ){
				cout<<cluster1[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster2:";
			for(size_t i=0; i<cluster2.size(); i++ ){
				cout<<cluster2[i]<<" ";
			}
			cout<<endl;
			cout<<"cluster3:";
			for(size_t i=0; i<cluster3.size(); i++ ){
				cout<<cluster3[i]<<" ";
			}
			cout<<endl;
		}

		//recursivley clsutering
		vector<vector<size_t> > tempVec;
		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster0 );
		if(0){
			cout<<"tempVec:";
			for(size_t i=0; i<tempVec.size(); i++){
				for(size_t j=0; j<tempVec[i].size(); j++ )
					cout<<tempVec[i][j]<<" ";
				cout<<endl;
			}
			cout<<endl;
		}
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster1 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster2 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );

		tempVec = tetrahedral_means(inputMatrix, limiteDis, cluster3 );
		resultVec.insert( resultVec.end(), tempVec.begin(), tempVec.end() );
	}
}

/**
 * @brief		our proposed algorithm
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector<vector<size_t> >
iterative(const vector<vector<_Tp > >& inputMatrix,
		const _Tp limiteDis, const bool& splitFlag = false ){
	cout<<endl<<"---iterative---"<<endl;
	clock_t t1, t2;
	t1 = clock();

	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for iterative clustering!!!---"<<endl;

	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}
	if(0){
		if(!checkMatric(inputMatrix) ){
			cerr<<"---check matric---";
		}
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t	matrixSize = inputMatrix.size();

	vector<size_t >		inputSeqVec;
	for(size_t i=0; i<inputMatrix.size(); i++ ){
		inputSeqVec.push_back(i );
	}
	if(0){
		cout<<"inputSeqVec.size()::"<<inputSeqVec.size()<<endl;
		cout<<"numer of largest distance:"<<computeNumOfLargestDis(inputMatrix )<<endl;
	}
	resultVec = tetrahedral(inputMatrix, limiteDis, inputSeqVec );
	if(1){
		cout<<"1111iterative running time"<<": "<< double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}

	if( splitFlag ){
		cout<<"cluster number before split:"<<resultVec.size()<<endl;
		resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
		cout<<"split running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	resultVec = sortCluster(resultVec);
	printCluster( resultVec );
	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

//	printCluster( resultVec );

	if(0){
		string	fileName = "iterative.txt";
		printfCluster( resultVec, fileName );
	}
//	checkClusterDis( resultVec, inputMatrix, limiteDis );
//	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );
	return resultVec;
}

/**
 * @brief    iterative clustering that have no mergeable clusters in the final result.
 * @param inputMatrix
 * @param limiteDis
 * @param squareDis
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
iterative_split(const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis, const bool squareDis = false ){
	//temp vector
	clock_t t1, t2;
	t1 = clock();

	if(1){
		cout<<endl<<"---iterative_split---"<<endl;
	}
	vector<vector<size_t> > 	resultVec;

	//same as iterative_means
	resultVec = iterative( inputMatrix,  limiteDis);

	//  split mergeable clusters
	//	resultVec = splitCluster_rude( resultVec, inputMatrix, limiteDis );
	resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
	cout<<"iterative_means_split running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;

	//get result
	resultVec = sortCluster( resultVec );
	printCluster( resultVec );
	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * @brief    iterative clustering with k-means like step in it and without mergeable clusters in final result.
 * @param inputMatrix
 * @param limiteDis
 * @param squareDis
 * @return
 */
template<typename _Tp>
vector<vector<size_t> >
iterative_means_split(const vector<vector<_Tp > >& inputMatrix,
		const _Tp limiteDis, const bool squareDis ){
	//temp vector
	clock_t t1, t2;
	t1 = clock();

	if(1){
		cout<<endl<<"---iterative_means_split---"<<endl;
	}
	vector<vector<size_t> > 	resultVec;

	printDisMatrix( inputMatrix );

	//same as iterative_means
	resultVec = iterative_means( inputMatrix,  limiteDis, squareDis );

	//  split mergeable clusters
	//	resultVec = splitCluster_rude( resultVec, inputMatrix, limiteDis );
	resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
	cout<<"iterative_means_split running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;

	resultVec = mergeCluster_rude( resultVec, inputMatrix, limiteDis );
	cout<<"resultVec.size(): "<<resultVec.size()<<endl;

	//get result
	resultVec = sortCluster( resultVec );
	printCluster( resultVec );
	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * @brief 		after each step in iterative, the clusters are redistributed as k-means do
 * @param inputMatrix
 * @param limiteDis
 * @param squareDis
 * @return
 */
template<typename _Tp>
vector<vector<size_t> >
iterative_means(const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis,
		const bool splitFlag, const bool squareDis ){

	if(1){
		cout<<endl<<"---iterative_means---"<<endl;
	}
	//error
	if(inputMatrix.empty() ){
		cerr<<"---no input data for iterative_kmeans clustering!!!---"<<endl;
	}
	if(limiteDis < 0 ){
		cerr<<"---please define neighbors correctly--- "<<endl;
	}
	if(0){
		if(!checkMatric(inputMatrix) ){
			//		cerr<<"---warnning check distance matric failed---"<<endl;;
		}
	}

	//temp vector
	vector<vector<size_t> > 	resultVec;

	//data number
	size_t	matrixSize = inputMatrix.size();

	vector<size_t >		inputSeqVec;
	for(size_t i=0; i<inputMatrix.size(); i++ ){
		inputSeqVec.push_back(i );
	}
	if(0){
		cout<<"inputSeqVec.size()::"<<inputSeqVec.size()<<endl;
		cout<<"numer of largest distance:"<<computeNumOfLargestDis(inputMatrix )<<endl;
	}

	//check running time
	clock_t t1, t2;
	t1 = clock();

	resultVec = tetrahedral_means(inputMatrix, limiteDis, inputSeqVec, squareDis );
	if(1){
		cout<<"iterative_means running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	if(1){
		cout<<"cluster number before merge:"<<resultVec.size()<<endl;
	}

	if( splitFlag ){
		cout<<"cluster number before split:"<<resultVec.size()<<endl;
		//merge cluster
		resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
		cout<<"after merge running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
		cout<<" after split: "<<endl;
	}
	resultVec = sortCluster(resultVec);

	if(0){
		string	fileName = "iterative.txt";
		printfCluster( resultVec, fileName );
	}

	printCluster( resultVec );
//	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * @brief 	check whether all data in $clusterVec$ have distances under limited
 * @param clusterVec
 * @param disMatrix
 * @param limiteDis
 * @return true if all distances are under $limiteDis$, false if else
 */
template< typename _Tp >
bool
allOne( const vector< size_t >& clusterVec, const vector< vector<_Tp> >& disMatrix,
		const _Tp& limiteDis){

	bool flag = true;
	for( size_t i=0; i<clusterVec.size() ; i++ ){
		for( size_t j=0; j<clusterVec.size(); j++ ){
			if( disMatrix[ clusterVec[i] ][ clusterVec[j] ] > limiteDis ){
				flag = false;
				break;
			}
		}
	}
	return flag;
}

/**
 *	@brief deep first search, this is used for clustering algorithm $dsf_clustering$
 * @param disMatrix
 * @param graphVec
 * @param index
 * @param limiteDis
 * @param visitedVec
 * @param oneCluster
 */
template<typename _Tp>
void
dfs_visit( const vector< vector<_Tp> >& disMatrix,
		const vector< vector< size_t > >& graphVec, const size_t& index, const _Tp limiteDis,
		vector<bool>&  visitedVec, vector<size_t>&  oneCluster ){

	//adjacent vector of node $index$
	vector< size_t > adjV;
	adjV = graphVec[ index ];

	//cluster vector
	vector< size_t > clusterVec( oneCluster );

	for( size_t i=0; i<adjV.size(); i++ ){

		if( !visitedVec[ adjV[i] ] ){

			//add new node to cluster and check whether it belong to that cluster
			clusterVec.clear();
			clusterVec = oneCluster;
			clusterVec.push_back( adjV[i] );

			if(0){
				cout<<"  clusterVec:";
				for( size_t i=0; i<clusterVec.size(); i++ ){
					cout<<" "<<clusterVec[i];
				}
				cout<<endl;
			}

			// if the new node could be merged into the cluster then deep search that node recursively
			if( allOne( clusterVec, disMatrix, limiteDis ) ){
				//mark node $adjV[i]$ as visited
				visitedVec[ adjV[i] ] = true;
				//obtain a new node
				oneCluster = clusterVec;
				dfs_visit( disMatrix, graphVec, adjV[i], limiteDis, visitedVec, oneCluster);
				if(0){
					cout<<"  index:"<<index<<" adj:"<<adjV[i]<<endl;
				}
				if(0){
					cout<<"  oneCluster: ";
					for( size_t i=0; i<oneCluster.size(); i++ ){
						cout<<" "<<oneCluster[i];
					}
					cout<<endl;
				}
			}
		}
	}
//	return clusterVec;
}

/**
 *
 * @param disMatrix
 * @param graphVec
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector< size_t > >
dfs_clustering( const vector< vector<_Tp> >& disMatrix, const vector< vector< size_t > >& graphVec,
		const _Tp limiteDis){
	//indicate whether node has been visited
	vector< bool > visitedVec( disMatrix.size(), false );
	//
	vector< size_t > oneCluster;

	//cluster result
	vector< vector<size_t> > clusterVec;

	for( size_t i=0; i<disMatrix.size(); i++ ){
		if(0){
			cout<<"i:"<<i<<endl;
		}

		if( !visitedVec[i] ){
			//new cluster
			oneCluster.clear();
			oneCluster.push_back( i );

			//node i has been visited
			visitedVec[i] = true;

			//recursively
			dfs_visit( disMatrix, graphVec, i, limiteDis, visitedVec, oneCluster);

			//get new cluster
			clusterVec.push_back( oneCluster );

			if(0){
				cout<<"  oneCluster:";
				for( size_t i=0; i<oneCluster.size(); i++ ){
					cout<<" "<<oneCluster[i];
				}
				cout<<endl;
			}
		}
	}

	return clusterVec;
}

/**
 * @brief	   a clustering algorihtm based on deep first graph search approach
 * @param inputMatrix
 * @param limiteDis
 * @return
 * @note		during deep first search, it should be better to search nearest neighbor first!!!
 */
template<typename _Tp>
vector< vector<size_t> >
graphClustering( const vector< vector<_Tp> >& inputMatrix, const _Tp limiteDis ){
	cout<<"---graph clustering---"<<endl;
	clock_t t1, t2;
	t1 = clock();

	if( inputMatrix.empty() ){
		cerr<<"---error: not input data for clustering!!!"<<endl;
	}
	if( limiteDis < 0 ){
		cerr<<"---warning: distance criteria below 0!!!"<<endl;
	}

	//neighbor graph
	vector< vector<size_t> > neighborVec;

	//node number
	size_t nodeNum = inputMatrix.size();

	//const neighbor graph
	neighborVec.clear();
	for( size_t i=0; i<nodeNum; i++ ){
		vector< size_t > tempV;
		tempV.clear();
		for( size_t j=0; j<nodeNum; j++ ){
			if( inputMatrix[i][j] < limiteDis &&  i != j ){
				if(0){
					cout<<"inputMatrix[i][j]:"<<inputMatrix[i][j]<<" limiteDis:"<<limiteDis<<endl;
				}
				tempV.push_back( j );
			}
		}
		neighborVec.push_back( tempV );
	}

	if(0){
		cout<<"adjacent matrix:"<<endl;
		for( size_t i=0; i<neighborVec.size(); i++ ){
			cout<<"i:"<<i<<endl;
			for( size_t j=0; j<neighborVec[i].size(); j++ ){
				cout<<" "<<neighborVec[i][j];
			}
			cout<<endl;
		}
	}

	vector< vector< size_t > > resultVec;
	resultVec = dfs_clustering( inputMatrix, neighborVec, limiteDis );
	if(1){
		cout<<"dsf_clustering running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	sortCluster( resultVec );
	printCluster( resultVec );
	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	//traversal the neighbor graph
	return resultVec;
}

/**
 * @brief		each time draw a circle for each data, substract the one with largest number of neighbors in that circle, and add them
 *  				into a new cluster
 * @param inputMatrix
 * @param limiteDis
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
circleClustering( const vector< vector<_Tp> >& inputMatrix, const _Tp diameter, const bool splitFlag  ){
	cout<<"---circle clustering---"<<endl;
	clock_t t1, t2;
	t1 = clock();

	_Tp limiteDis = diameter / 2;
	if( inputMatrix.empty() ){
		cerr<<"---error: not input data for clustering!!!"<<endl;
	}
	if( limiteDis < 0 ){
		cerr<<"---warning: distance criteria below 0!!!"<<endl;
	}

	//cluster result
	vector< vector<size_t> > resultVec;

	//indicate whehter data $i$ has been removed
	vector<bool> visitedVec( inputMatrix.size(), 0 );

	//neighbor number for each data
	vector< size_t > neighborCount( inputMatrix.size(), 0 );

	//neighbor
	vector< vector<size_t> > neighborVec;

	bool flag = true;

	while( flag ){
		if(1){
			size_t numOfUnvisited = 0;
			size_t numOfVisited = 0;
			for( size_t num=0; num<visitedVec.size(); num++ ){
				if( !visitedVec[ num ] ){
					numOfUnvisited ++;
				}else{
					numOfVisited ++;
				}
			}

			char ss[20] = "";
			size_t matrixSize = visitedVec.size();
			for( size_t i =0; i<size_t(20*((float) numOfVisited/matrixSize ) ); i++ ){
				ss[i] = '-';
			}
			if(size_t(20*((float) numOfVisited/matrixSize ) )  <19  ){
				ss[ size_t(20*((float) numOfVisited/matrixSize ) ) ] = '>';
			}
			for( size_t i=size_t(20*((float) numOfVisited/matrixSize ) ) + 1; i<19; i++ ){
				ss[i] = ' ';
			}
			cerr<<"\r|"<<ss<<"|"<<" complete link running time: "<< double( clock() - t1 )/ CLOCKS_PER_SEC;

		}
		flag = false;
		neighborVec.clear();
		//get neighbors for $i$, contain them in $neighborVec$
		for( size_t i=0; i<inputMatrix.size(); i++ ){
			vector<size_t> vec;
			vec.clear();
			if( !visitedVec[i] ){
				for( size_t j=0; j<inputMatrix[i].size(); j++ ){
					if( !visitedVec[j] && inputMatrix[i][j] < limiteDis ){
						vec.push_back( j );
					}
				}
				neighborVec.push_back( vec );
			}
		}

		if( neighborVec.empty() ){
//			cout<<"---break---"<<endl;
			break;
		}

		//get the largest group of $neighborVec$
		size_t max = 0;
		size_t index = 0;
		for( size_t i=0; i<neighborVec.size(); i++ ){
			if( neighborVec[i].size() > max ){
				max = neighborVec[i].size();
				index = i;
			}
		}

		if(0){
			cout<<"---neighborVec---"<<endl;
			for( size_t i=0; i<neighborVec[ index ].size(); i++ ){
				cout<<" "<<neighborVec[ index ][ i ];
			}
			cout<<endl;
		}

		//get new cluster, and label $visitedVec$
		resultVec.push_back( neighborVec[ index ] );

		//update
		for( size_t i =0; i<neighborVec[index].size(); i++ ){
			if( visitedVec[ neighborVec[index][i] ] ){
				cerr<<"---"<<neighborVec[index][i]<<" already visited???"<<endl;
			}
			visitedVec[ neighborVec[index][i] ] = true;

			flag = true;
		}
	}
	if(1){
		cout<<"circleClustering running time:"<< float( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
	}
	sortCluster( resultVec );
	if(splitFlag){
		cout<<"cluster number before split:"<<resultVec.size()<<endl;
		resultVec = splitCluster_least( resultVec, inputMatrix, diameter );
		resultVec = sortCluster( resultVec );
		cout<<"result after split:"<<endl;
	}
	printCluster( resultVec, inputMatrix );

	//note the limiteDis is half of the cluster range
	checkClusterDis( resultVec, inputMatrix, diameter );
	checkLargestInnerDis( inputMatrix, resultVec, diameter );

	return resultVec;
}

template<typename _Tp>
vector< vector<size_t> >
biDiversive(  const vector< vector<_Tp> >& inputMatrix, const _Tp diameter,
		const bool& meanFlag = true, const bool& splitFlag = true ){
	cout<<"---biDiversive clustering---"<<endl;
	clock_t t1, t2;
	t1 = clock();

	_Tp limiteDis = diameter;
	if( inputMatrix.empty() ){
		cerr<<"---error: not input data for biDiversive clustering!!!"<<endl;
	}
	if( limiteDis < 0 ){
		cerr<<"---warning: distance criteria below 0 in biDiversive!!!"<<endl;
	}

	//cluster result
	vector< vector<size_t> > resultVec;

	//temp cluster, first contains indexes os all struct
	vector<size_t> oneClu( inputMatrix.size(), 0);
	for( size_t i=0; i<inputMatrix.size(); i++ ){
		oneClu[i] = i;
	}
	resultVec.push_back( oneClu );

	for( size_t i=0; i<resultVec.size(); i++ ){
		if(0){
			cout<<"i:"<<i<<endl;
		}

		//get diameter endpoints of cluster
		_Tp maxDis = 0;
		size_t indeBegin = 0;
		size_t indeEnd = 0;
		for( size_t m=0; m<resultVec[i].size()-1; m++ ){
			for( size_t n=m; n<resultVec[i].size(); n++ ){
				if( maxDis < inputMatrix[ resultVec[i][m] ][ resultVec[i][n] ] ){
					maxDis =  inputMatrix[ resultVec[i][m] ][ resultVec[i][n] ] ;
					indeBegin =  m ;
					indeEnd =  n ;
				}
			}
		}

		if(0){
			cout<<"indeBegin:"<<indeBegin<<"  indeEnd:"<<indeEnd<<"  maxDis:"<<maxDis<<endl;
		}

		//split cluster
		if( maxDis > diameter ){
			vector<size_t> clusterFirst;
			vector<size_t> clusterSec;
			for( size_t m=0; m<resultVec[i].size(); m++ ){
				if( inputMatrix[ resultVec[i][m] ][ resultVec[i][ indeBegin] ]  <
						inputMatrix[ resultVec[i][m] ][ resultVec[i][ indeEnd ]  ]  ){
					clusterFirst.push_back( resultVec[i][ m ] );
				}else{
					clusterSec.push_back( resultVec[i][ m ] );
				}
				if(0){
					cout<<" resultVec[i][m]:"<<resultVec[i][m]<<" resultVec[i][ indeBegin]:"
							<<resultVec[i][ indeBegin]<<" resultVec[i][ indeEnd ]:"<<resultVec[i][ indeEnd ]<<endl;
				}
			}

			//reset cluster $clusterFirst$ $clusterSec$
			if( meanFlag ){
				vector<vector< size_t > > cluster0and1;
				cluster0and1.push_back( clusterFirst );
				cluster0and1.push_back( clusterSec );
				cluster0and1 = rearrange( inputMatrix, cluster0and1 );
				clusterFirst = cluster0and1[0];
				clusterSec = cluster0and1[1];
			}

			if(0){
				cout<<"clusterFirst"<<endl;
				for( size_t k=0; k<clusterFirst.size(); k++ ){
					cout<<" "<<clusterFirst[k];
				}
				cout<<endl;
				cout<<"clusterSec"<<endl;
				for( size_t k=0; k<clusterSec.size(); k++ ){
					cout<<" "<<clusterSec[k];
				}
				cout<<endl;
			}

			if( !clusterFirst.empty() ){
				resultVec.push_back( clusterFirst );
			}
			if( !clusterSec.empty() ){
				resultVec.push_back( clusterSec );
			}
			resultVec.erase( resultVec.begin() + i );
			i--;
		}
	}
	if(1){
		t2 = clock();
		float time = ( (float)t2 - ( float )t1 ) / CLOCKS_PER_SEC ;
		cout<<"completeLink running time:"<<time<<endl;
	}

	//merge cluster result
	resultVec = sortCluster( resultVec );
	if(splitFlag){
		cout<<"cluster number before split:"<<resultVec.size()<<endl;
		resultVec = splitCluster_least( resultVec, inputMatrix, limiteDis );
		cout<<"result after split:"<<endl;
	}
	printCluster( resultVec );

	checkClusterDis( resultVec, inputMatrix, limiteDis );
	checkLargestInnerDis( inputMatrix, resultVec, limiteDis );

	return resultVec;
}

/**
 * Guenoche, Efficient Algorithms for Divisive Hierarchical Clustering with the Diameter Criterion(1991)
 * @param inputMatrix
 * @return
 */
template<typename _Tp>
vector< vector<size_t> >
HubertPartition( const vector< vector<_Tp> >& inputMatrix ){
	vector<size_t>  			U;                                                        //
	vector<_Tp>				lamda( inputMatrix.size(), 0 ) ;    		//largest dissimilarity of ith
	vector<_Tp>				p( inputMatrix.size(), 0 );					//
	vector<bool>				gama( inputMatrix.size() );					//color true for R false for B

	//initialization
	//note, this starts from 1 to n-1 corresponds to 2 to n described in paper
	lamda[0] = 0;
	p[0] = 0;
	gama[0] = true;
	for( size_t j=1; j< inputMatrix.size(); j++ ){
		U.push_back( j );
		lamda[j] = inputMatrix[0][j];
		p[j] = 0;
	}

//	cout<<endl<<"inputMatrix:"<<inputMatrix.size()<<endl;

	//main step
	while( ! U.empty() ){
		size_t			k;
		_Tp				max = -1e10;
		for( size_t i=0; i<U.size(); i++ ){
			if( max < lamda[ U[i] ] ){
				max = lamda[ U[i] ];
				k= U[i];
			}
		}
//		size_t k= max_element( lamda.begin(), lamda.end() ) - lamda.begin() ;
		if( gama[ p[k] ] = true ){
			gama[ k ] = false;
		}else{
			gama[ k ] = false;
		}

		U.erase( find( U.begin(), U.end(), k ) );
		if(0){
			cout<<"U:"<<U.size()<<" k:"<<k<<endl;
			for( size_t i=0; i<U.size(); i++ ){
				cout<<" "<<U[i];
			}
			cout<<endl;
		}
		for( size_t j=0; j<U.size(); j++ ){
			if( inputMatrix[k][j] > lamda[j] ){
				lamda[j] = inputMatrix[k][j];
				p[j] = k;
			}
		}

	}
	vector<size_t>			v1, v2;
	for( size_t i=0; i<gama.size(); i++ ){
		if( gama[i] == true ){
			v1.push_back( i );
		}else{
			v2.push_back(i);
		}
	}
	vector< vector<size_t> >	all;
	all.push_back( v1 );
	all.push_back( v2 );
//	cout<<"v1: ";
//	for( size_t i=0; i<v1.size(); i++ ){
//		cout<<v1[i]<<" ";
//	}
//	cout<<endl;
//	cout<<"v2: ";
//	for( size_t i=0; i<v2.size(); i++ ){
//		cout<<v2[i]<<" ";
//	}
//	cout<<endl;
	return all;
}

template< typename _Tp >
void
HubertAlgorithmIterative(		const vector< vector<_Tp> >& matrix,
													const _Tp diameter,
													vector<size_t>& indexVec,
													vector< vector<size_t> >&  resultVec ){

	if( matrix.empty() ){
		cerr<<"---error: not input data for biDiversive clustering!!!"<<endl;
		return ;
	}
	if( diameter < 0 ){
		cerr<<"---warning: distance criteria below 0 in biDiversive!!!"<<endl;
		return ;
	}



	if( indexVec.size() == 1 ){
		resultVec.push_back( indexVec );
		return ;
	}

	if(0){
		cout<<"matrix:"<<matrix.size()<<endl;
		for( size_t i=0; i<indexVec.size(); i++ ){
			cout<<indexVec[i]<<" ";
		}
		cout<<endl;
	}

	vector< vector<_Tp> >			inputMatrix( indexVec.size(), vector<_Tp>( indexVec.size(), 0 ) );
	for( size_t i=0; i<indexVec.size(); i++ ){
		inputMatrix[i][i] = 0;
		for( size_t j=i+1; j<indexVec.size(); j++ ){
			inputMatrix[i][j] = matrix[ indexVec[i] ][ indexVec[j] ];
		}
	}

	_Tp				largest_dis = largestDis( inputMatrix );
	if( largest_dis < diameter ){
		resultVec.push_back( indexVec );
		return ;
	}

	vector< vector<size_t> > 	subCluster = HubertPartition( inputMatrix );
	vector<size_t>					sub1( subCluster[0].size() ) ;
	vector<size_t>					sub2( subCluster[1].size() ) ;
	vector<vector<_Tp> >		subMatrix1( sub1.size(), vector<_Tp>(sub1.size()) );
	vector<vector<_Tp> >		subMatrix2( sub2.size(), vector<_Tp>(sub2.size()) );
	for( size_t i=0; i<subCluster[0].size(); i++ ){
		sub1[i] = indexVec[ subCluster[0][i] ];
	}

	for( size_t i=0; i<subCluster[1].size(); i++ ){
		sub2[i] = indexVec[ subCluster[1][i] ];
	}

	if(0){
		cout<<"sub1:";
		for( size_t i=0; i<sub1.size(); i++ ){
			cout<<sub1[i]<<" ";
		}
		cout<<endl;
		cout<<"sub2:";
		for( size_t i=0; i<sub2.size(); i++ ){
			cout<<sub2[i]<<" ";
		}
		cout<<endl;
	}

	HubertAlgorithmIterative( matrix, diameter, sub1, resultVec );

	HubertAlgorithmIterative( matrix, diameter, sub2, resultVec );
}

template< typename _Tp >
vector< vector<size_t> >
Hubert(		const vector< vector<_Tp> >& matrix,
												const _Tp diameter){
	vector<size_t>					indexVec( matrix.size() );
	vector< vector<size_t> >	resultVec;
	for( size_t i=0; i<indexVec.size(); i++ ){
		indexVec[i] = i;
	}
	printDisMatrix( matrix );
	HubertAlgorithmIterative( matrix, diameter, indexVec, resultVec );
	sortCluster( resultVec );
	cout<<"diameter:"<<diameter<<"  cluster result:"<<endl;
	for( size_t i=0; i<resultVec.size(); i++ ){
		cout<<i<<": ";
		for( size_t j=0; j<resultVec[i].size(); j++ ){
			cout<<resultVec[i][j]<<" ";
		}
		cout<<endl;
	}
	return									resultVec;
}
