#include"isomorphism.h"

using namespace std;

vector<vector<int> >
bfs_mol( const vector<string> nodeNameVec,
         const vector<pair<int, int> > edgeVec,
         const int startIndex ){

    int nodeNum = nodeNameVec.size();
    int edgeNum = edgeVec.size();

    vector<bool> visited(nodeNum, false);

    BfsNode n0(startIndex, 0);
    queue<BfsNode> q;
    q.push(n0);
    visited[startIndex] = true;

    vector<BfsNode> allNodes;

    while( !q.empty() ){
        BfsNode newNode = q.front();
        q.pop();
        allNodes.push_back(newNode);
        int index = newNode.index;
        int depth= newNode.depth;
//        cout<<"index:"<<index<<"  depth:"<<depth<<endl;
        //get all edges from 'index'
        for( int i=0; i<edgeNum; i++ ){
            if( edgeVec[i].first == index ){
                int newIndex = edgeVec[i].second;
                if( !visited[newIndex] ){
                    visited[newIndex]=true;
                    BfsNode n(newIndex, depth+1);
                    q.push(n);
                }
            }
            if( edgeVec[i].second == index ){
                int newIndex = edgeVec[i].first;
                if( !visited[newIndex] ){
                    visited[newIndex]=true;
                    BfsNode n(newIndex, depth+1);
                    q.push(n);
                }
            }
        }
    }

    int depth = allNodes.back().depth;
    vector<vector<int> > bfsResult(depth+1);
    for( int i=0; i<allNodes.size(); i++ ){
        bfsResult[ allNodes[i].depth ].push_back( allNodes[i].index );
    }

    if(0){
        cout<<"startIndex:"<<startIndex<<endl;
        for( size_t i=0; i<bfsResult.size(); i++ ){
            for( size_t j=0; j<bfsResult[i].size(); j++ ){
                cout<<bfsResult[i][j]<<" ";
            }
            cout<<endl;
        }
    }
    return bfsResult;
}

/**
 * @brief vec_equal return true if vec1 and vec2 have same num of same strings
 * @param vec1
 * @param vec2
 * @return
 */
bool
vec_equal( vector<string> vec1, vector<string> vec2 ){
    if( vec1.size() != vec2.size() ){
        return false;
    }else{
        vector<bool> exist(vec2.size(), false);
        for( size_t i=0; i<vec1.size(); i++ ){
            bool flag=false;
            for( size_t j=0; j<vec2.size(); j++ ){
                if( vec1[i] == vec2[j] && exist[j]==false ){
                    exist[j]=true;
                    flag=true;
                    break;
                }
            }
            if( !flag ){
                return false;
            }
        }
        return true;
    }
}

vector<string>
getStringVec(  const vector<string> nodeNameVec, const vector<int> indexVec ){
    vector<string> strVec;
    strVec.clear();
    for( size_t i=0; i<indexVec.size(); i++ ){
        if( indexVec[i] > nodeNameVec.size() ){
            return strVec;
        }
        strVec.push_back( nodeNameVec[ indexVec[i] ] );
    }
    return strVec;
}

bool
bfs_equal(  const vector<string> nodeNameVec1,
            const vector<pair<int, int> > edgeVec1,
            const int startIndex1,
            const vector<string> nodeNameVec2,
            const vector<pair<int, int> > edgeVec2,
            const int startIndex2 ){

    vector<vector<int> > bfs1=bfs_mol(nodeNameVec1, edgeVec1, startIndex1);
    vector<vector<int> > bfs2=bfs_mol(nodeNameVec2, edgeVec2, startIndex2);

    if(0){
        cout<<"nodeA:"<<startIndex1<<endl;
        for( size_t i=0; i<bfs1.size(); i++ ){
            for( size_t j=0; j<bfs1[i].size(); j++ ){
                cout<<bfs1[i][j]<<" ";
            }
            cout<<endl;
        }
        cout<<"nodeB:"<<startIndex2<<endl;
        for( size_t i=0; i<bfs2.size(); i++ ){
            for( size_t j=0; j<bfs2[i].size(); j++ ){
                cout<<bfs2[i][j]<<" ";
            }
            cout<<endl;
        }
    }

    if( bfs1.size() != bfs2.size() ){
        return false;
    }else{
        for( size_t depth=0; depth<bfs1.size(); depth++ ){
            if( bfs1[depth].size() != bfs2[depth].size() ){
                return false;
            }else{
                vector<string> lay1=getStringVec(nodeNameVec1, bfs1[depth]);
                vector<string> lay2=getStringVec(nodeNameVec2, bfs2[depth]);
                if( !vec_equal(lay1, lay2) ){
                    return false;
                }
            }
        }
        return true;
    }
}

vector<vector<int> >
isomorphism(  const vector<string> nodeNameVec1,
              const vector<pair<int, int> > edgeVec1,
              const vector<string> nodeNameVec2,
              const vector<pair<int, int> > edgeVec2 ){
    vector<vector<int> > matrix(nodeNameVec1.size(), vector<int>(nodeNameVec2.size()));
    for( size_t i=0; i<nodeNameVec1.size(); i++ ){
        for( size_t j=0; j<nodeNameVec2.size(); j++ ){
//            cout<<"equal: "<<i<<" "<<j<<endl;
            bool flag = bfs_equal( nodeNameVec1, edgeVec1, i, nodeNameVec2, edgeVec2, j);
            if( flag ){
                matrix[i][j] = 1;
            }else{
                matrix[i][j] = 0;
            }
        }
    }

    return matrix;
}
