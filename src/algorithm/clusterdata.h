/*
 * clusterdata.h
 *
 *  Created on: Apr 11, 2013
 *      Author: stan
 */

#ifndef CLUSTERDATA_H_
#define CLUSTERDATA_H_

#include<vector>

using namespace std;

template<typename _Tp>
void
printDisMatrix( const vector< vector<_Tp> >& disMatrix );

template<typename _Tp>
vector<vector<size_t> >
singleLink(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis );

template<typename _Tp>
vector<vector<size_t> >
singleLink_limit(const vector<vector<_Tp> >& inputMatrix, const _Tp& limiteDis );

template<typename _Tp>
vector< vector< size_t > >
completeLink( const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis, const bool& splitFlag = false );

template<typename _Tp>
vector<vector<size_t> >
iterative_means(const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis,
		const bool& splitFlag = false, const bool& squareDis = false );

template<typename _Tp>
vector<vector<size_t> >
iterative_means_split(const vector<vector<_Tp > >& inputMatrix, const _Tp limiteDis, const bool squareDis = false );

template<typename _Tp>
vector< vector< size_t > >
fastHierarchical( const vector< vector< _Tp > >& inputMatrix, const _Tp& limiteDis );

template<typename _Tp, typename _Size>
vector<vector<_Size> >
k_means( const vector<vector<_Tp> >& inputMatrix, const _Size clusterNum );

template<typename _Tp>
vector< vector<size_t> >
graphClustering( const vector< vector<_Tp> >& inputMatrix, const _Tp limiteDis );

template<typename _Tp>
vector< vector<size_t> >
circleClustering( const vector< vector<_Tp> >& inputMatrix, const _Tp limiteDis  );

//template<typename _Tp>
//_Tp
//getMinPair( vector<vector<_Tp > >& distMatrix , vector<size_t> groupVec );

/**
 * note that the template could not be initialized in a .h file and implemented in another .cpp file, but could be
 * overcome in this funcky way, that is include a .cpp file in the last of the .h file.
 * @see http://stackoverflow.com/questions/3040480/c-template-function-compiles-in-header-but-not-implementation
 */
#include"clusterdata.hpp"

#endif /* CLUSTERDATA_H_ */
