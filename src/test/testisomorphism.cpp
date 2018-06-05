#include"../algorithm/isomorphism.h"
#include"../mol/molecular.h"
#include"../algorithm/graph.h"

using namespace std;

void
testIsomorphism2(){
    vector<string> nameVec, nameVec2;
    vector<pair<int, int> > edgeVec, edgeVec2;

    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("C"));
    nameVec.push_back(string("O"));
    nameVec.push_back(string("O"));
    nameVec.push_back(string("C"));

    edgeVec.push_back(pair<int, int>(0, 1));
    edgeVec.push_back(pair<int, int>(0, 5));
    edgeVec.push_back(pair<int, int>(1, 2));
    edgeVec.push_back(pair<int, int>(5, 6));
    edgeVec.push_back(pair<int, int>(2, 3));
    edgeVec.push_back(pair<int, int>(5, 4));
    edgeVec.push_back(pair<int, int>(4, 3));
    edgeVec.push_back(pair<int, int>(3, 7));
    edgeVec.push_back(pair<int, int>(7, 8));
    edgeVec.push_back(pair<int, int>(8, 9));

    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("C"));
    nameVec2.push_back(string("O"));
    nameVec2.push_back(string("O"));
    nameVec2.push_back(string("C"));

    edgeVec2.push_back(pair<int, int>(0, 1));
    edgeVec2.push_back(pair<int, int>(0, 5));
    edgeVec2.push_back(pair<int, int>(1, 2));
    edgeVec2.push_back(pair<int, int>(2, 3));
    edgeVec2.push_back(pair<int, int>(2, 7));
    edgeVec2.push_back(pair<int, int>(3, 4));
    edgeVec2.push_back(pair<int, int>(4, 5));
    edgeVec2.push_back(pair<int, int>(4, 6));
    edgeVec2.push_back(pair<int, int>(7, 8));
    edgeVec2.push_back(pair<int, int>(8, 9));

    vector<vector<int> > bfsR=bfs_mol(nameVec, edgeVec, 0);
    for( size_t i=0; i<bfsR.size(); i++ ){
        for( size_t j=0; j<bfsR[i].size(); j++ ){
            cout<<bfsR[i][j]<<" ";
        }
        cout<<endl;
    }

    vector<vector<int> > matrix=isomorphism(nameVec,
                                            edgeVec, nameVec2, edgeVec2);

    for( size_t i=0; i<matrix.size(); i++ ){
        for( size_t j=0; j<matrix[i].size(); j++ ){
            cout.width(5);
            cout<<left<<matrix[i][j];
        }
        cout<<endl;
    }
}

void
testIsomorphism(){
//    string fileName1 = "1f8d_ligand.mol2";
//    string fileName2 = "final_lig0.pdb";
//    Molecular lig1;
//    Molecular lig2;
//    lig1 = readMolecularFile( fileName1 ).front();
//    lig2 = readMolecularFile( fileName2 ).front();

    string fileName = "1f8d_ligand.mol2";
    string fileName0 = "final_lig0.pdb";
    string fileName1 = "final_lig1.pdb";
    string fileName2 = "final_lig2.pdb";
    string fileName3 = "final_lig3.pdb";
    string fileName4 = "final_lig4.pdb";
    string fileName5 = "final_lig5.pdb";
    string fileName6 = "final_lig6.pdb";
    string fileName7 = "final_lig7.pdb";
    string fileName8 = "final_lig8.pdb";
    string fileName9 = "final_lig9.pdb";
    string fileName10 = "final_lig10.pdb";
    string fileName11 = "final_lig11.pdb";
    string fileName12 = "final_lig12.pdb";
    string fileName13 = "final_lig13.pdb";
    string fileName14 = "final_lig14.pdb";
    string fileName15 = "final_lig15.pdb";
    string fileName16 = "final_lig16.pdb";
    string fileName17 = "final_lig17.pdb";
    string fileName18 = "final_lig18.pdb";
    string fileName19 = "final_lig19.pdb";
    string fileName20 = "final_lig20.pdb";
    string fileName21 = "final_lig21.pdb";
    string fileName22 = "final_lig22.pdb";
    string fileName23 = "final_lig23.pdb";
    string fileName24 = "final_lig24.pdb";

    Molecular lig, lig0, lig1, lig2, lig3, lig4, lig5, lig6, lig7, lig8, lig9, lig10, lig11, lig12, lig13, lig14, lig15,
            lig16, lig17, lig18, lig19, lig20, lig21, lig22, lig23, lig24  ;

    lig = readMolecularFile( fileName ).front();
    lig0 = readMolecularFile( fileName0 ).front();
    lig1 = readMolecularFile( fileName1 ).front();
    lig2 = readMolecularFile( fileName2 ).front();
    lig3 = readMolecularFile( fileName3 ).front();
    lig4 = readMolecularFile( fileName4 ).front();
    lig5 = readMolecularFile( fileName5 ).front();
    lig6 = readMolecularFile( fileName6 ).front();
    lig7 = readMolecularFile( fileName7 ).front();
    lig8 = readMolecularFile( fileName8 ).front();
    lig9 = readMolecularFile( fileName9 ).front();
    lig10 = readMolecularFile( fileName10 ).front();
    lig11 = readMolecularFile( fileName11 ).front();
    lig12 = readMolecularFile( fileName12 ).front();
    lig13 = readMolecularFile( fileName13 ).front();
    lig14 = readMolecularFile( fileName14 ).front();
    lig15 = readMolecularFile( fileName15 ).front();
    lig16 = readMolecularFile( fileName16 ).front();
    lig17 = readMolecularFile( fileName17 ).front();
    lig18 = readMolecularFile( fileName18 ).front();
    lig19 = readMolecularFile( fileName19 ).front();
    lig20 = readMolecularFile( fileName20 ).front();
    lig21 = readMolecularFile( fileName21 ).front();
    lig22 = readMolecularFile( fileName22 ).front();
    lig23 = readMolecularFile( fileName23 ).front();
    lig24 = readMolecularFile( fileName24 ).front();


    //----------------------------------------------------------------
    vector<Atom> atVec1 = lig1.get_atomVec();
    vector<Bond> bdVec1 = lig1.get_bondVec();
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

    if(0){
        for( size_t i=0; i<atVecNoH1.size(); i++ ){
            Atom a=atVecNoH1[i];
            cout<<i<<"-> ";
            a.print();
        }
        for( size_t i=0; i<connectVecNoH1.size(); i++ ){
            cout<<i<<"-> "<<connectVecNoH1[i].first<<" "<<connectVecNoH1[i].second<<endl;
        }
    }

    //--------------------------------------------------------------
    vector<Atom> atVec2 = lig2.get_atomVec();
    vector<Bond> bdVec2 = lig2.get_bondVec();
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

    if(0){
        for( size_t i=0; i<atVecNoH2.size(); i++ ){
            Atom a=atVecNoH2[i];
            cout<<i<<"-> ";
            a.print();
        }
        for( size_t i=0; i<connectVecNoH2.size(); i++ ){
            cout<<i<<"-> "<<connectVecNoH2[i].first<<" "<<connectVecNoH2[i].second<<endl;
        }
    }

    //--------------------------------------------------------------
    vector< vector<int> > matrix = isomorphism( atNameVec1, connectVecNoH1,
                                                atNameVec2, connectVecNoH2);

    //min dist; below may not be correct RMSD
    float allDis=0;
    for( size_t i=0; i<matrix.size(); i++ ){
        Atom a1=atVecNoH1[i];
        vector<float> disVec;
        for( size_t j=0; j<matrix[i].size(); j++ ){
            if( matrix[i][j] == 1 ){
                Atom a2 = atVecNoH2[j];
                float dis = atomDis(a1, a2);
                disVec.push_back(dis);
            }
        }
        float minDis;
        if( disVec.empty() ){
            minDis=10e10;
        }else{
            minDis = *min_element( disVec.begin(), disVec.end() );
        }
        allDis += minDis;
    }
    float averageDis = allDis / (float)matrix.size();
    cout<<"isomorphism dis:"<<averageDis<<endl;
    float dis = isomorphismRMSD(lig1,lig2);
    cout<< isomorphismRMSD(lig,lig0)<<endl;
    cout<< isomorphismRMSD(lig,lig1)<<endl;
    cout<< isomorphismRMSD(lig,lig2)<<endl;
    cout<< isomorphismRMSD(lig,lig3)<<endl;
    cout<< isomorphismRMSD(lig,lig4)<<endl;
    cout<< isomorphismRMSD(lig,lig5)<<endl;
    cout<< isomorphismRMSD(lig,lig6)<<endl;
    cout<< isomorphismRMSD(lig,lig7)<<endl;
    cout<< isomorphismRMSD(lig,lig8)<<endl;
    cout<< isomorphismRMSD(lig,lig9)<<endl;
    cout<< isomorphismRMSD(lig,lig10)<<endl;
    cout<< isomorphismRMSD(lig,lig11)<<endl;
    cout<< isomorphismRMSD(lig,lig12)<<endl;
    cout<< isomorphismRMSD(lig,lig13)<<endl;
    cout<< isomorphismRMSD(lig,lig14)<<endl;
    cout<< isomorphismRMSD(lig,lig15)<<endl;
    cout<< isomorphismRMSD(lig,lig16)<<endl;
    cout<< isomorphismRMSD(lig,lig17)<<endl;
    cout<< isomorphismRMSD(lig,lig18)<<endl;
    cout<< isomorphismRMSD(lig,lig19)<<endl;
    cout<< isomorphismRMSD(lig,lig20)<<endl;
    cout<< isomorphismRMSD(lig,lig21)<<endl;
    cout<< isomorphismRMSD(lig,lig22)<<endl;
    cout<< isomorphismRMSD(lig,lig23)<<endl;
    cout<< isomorphismRMSD(lig,lig24)<<endl;
}

//void
//testIsomorphism(int argc, char** argv){
//    string name1(argv[1]);
//    Molecular mol1 = readMolecularFile( argv[1] ).front();
//    Molecular mol2 = readMolecularFile( argv[2] ).front();
//    float dis = isomorphismRMSD(mol1, mol2);
//    float olDis = overlayRMSD(mol1, mol2);
//    cout<<"dis:"<<dis<<endl;
//    cout<<"olDis:"<<olDis<<endl;
//    return;
//}

void
testIsomorphism(int argc, char** argv){
    Molecular mol = readMolecularFile("corina_ligand.mol2").front();
    vector<vector<size_t> > matrix;
    matrix = connectMatrix( mol.get_atomVec(), mol.get_bondVec() );

    return ;
}
