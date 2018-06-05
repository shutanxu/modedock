/*
 * comparator.h
 *
 *  Created on: Sep 28, 2013
 *      Author: stan
 */

#ifndef COMPARATOR_H_
#define COMPARATOR_H_

#include<vector>
//#include<pair>
#include <utility>

using namespace std;

template<class T1, class T2 >
struct
cmpPairFirst{
	bool	operator()( const std::pair<T1, T2>& left, const std::pair<T1, T2>& right ){
		return		left.first < right.first;
	}
};

template<typename T, typename T2>
struct
cmpPairSecond{
	bool	operator()( const pair<T, T2>& left, const pair<T, T2>& right ){
		return		left.second < right.second;
	}
};

#endif /* COMPARATOR_H_ */
