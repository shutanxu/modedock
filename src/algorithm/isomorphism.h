#ifndef ISOMORPHISM_H
#define ISOMORPHISM_H

#include<iostream>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#include<queue>
//#include"../mol/molecular.h"

using namespace std;

struct BfsNode{
    BfsNode(int i, int d):index(i),depth(d){}
    int index;
    int depth;
};

/**
 * @brief bfs_mol broad first search graph, the return value is nodes in each level
 * @param nodeNameVec
 * @param edgeVec
 * @return
 */
vector<vector<int> >
bfs_mol( const vector<string> nodeNameVec,
         const vector<pair<int, int> > edgeVec,
         const int startIndex );

bool
vec_equal( vector<string> vec1, vector<string> vec2 );

vector<string>
getStringVec(  const vector<string> nodeNameVec,
               const vector<int> indexVec );

bool
bfs_equal(  const vector<string> nodeNameVec1,
            const vector<pair<int, int> > edgeVec1,
            const int startIndex1,
            const vector<string> nodeNameVec2,
            const vector<pair<int, int> > edgeVec2,
            const int startIndex2 );

vector<vector<int> >
isomorphism(  const vector<string> nodeNameVec1,
              const vector<pair<int, int> > edgeVec1,
              const vector<string> nodeNameVec2,
              const vector<pair<int, int> > edgeVec2 );

#endif // ISOMORPHISM_H

