/*
 * graph.h
 *
 *  Created on: Nov 4, 2013
 *      Author: stan
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include"../mol/atom.h"
#include"../mol/bond.h"

using namespace std;

//vector< vector<size_t> >
//detectRing( const Molecular& mol );

class Node{
public:
    Node(const Atom a, const int d){
        at = a;
        depth = d;
    }
    Atom at;
    int depth;
};

bool
visited( const vector<pair<size_t, bool> >& visitedVec,
		 const size_t& index );

void
setVisited( vector<pair<size_t, bool> >& visitedVec,
		    const size_t& index );

bool
allVisited( const vector<pair<size_t, bool> >& 	visitedVec );

void
dfsMol( const vector<Bond>&          bdVec,
		const size_t&                index,
		vector<pair<size_t, bool> >& visitedVec,
		vector<size_t>&              historyVec );

vector<Atom>
dfsMol(	const vector<Atom>&          atomVec,
		const vector<Bond>&          bdVec,
		const	size_t&              directionStartIndex,
		const	size_t&              directionNextIndex );

void
dfsMol( const vector<Bond>&          bdVec,
		const size_t&                index,
		const size_t&                fatherIndex,
		vector<pair<size_t, bool> >& visitedVec,
		vector<size_t>&              historyVec,
		vector<vector<size_t> >&     ringVec);

vector<Atom>
bfsMol( const vector<Atom>&          atomVec,
		const vector<Bond>&          bdVec,
		const size_t&                startAtomIndex );

vector<vector<size_t> >
connectMatrix(const vector<Atom> atVec, const vector<Bond> bd);

vector<Bond>
getBfsBonds( const vector<Bond>& bdVec, const size_t& startAtomIndex );

#endif /* GRAPH_H_ */
