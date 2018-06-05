/*
 * testDock.cpp
 *
 *  Created on: Oct 4, 2014
 *      Author: stan
 */


#include"testDock.h"
#include"../dock/deDock.h"
#include"../forcefield/gaff.h"

using namespace std;

int
testDock( int argc, char** argv ){
    string fileNameMol = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/1ROB_protein.mol2";
    string fileNameLig = "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1rob/ligand_corina.mol2";

    string goldParFile("ff.dat");
    Sybyl gold(goldParFile);

    Molecular mol;
    Molecular lig;
    mol = readMolecularFile( fileNameMol ).front();
    lig = readMolecularFile( fileNameLig ).front();

    vector<Molecular> modelVec;

    return 0;
}

int
test_DE_Dock(){

	string parameterFile;

//	Molecular findMol = deRigidLigandDock(parameterFile);

	return 0;
}

int
test_DE_flexLigandDock(){

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1com/");
	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1cqp/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1g9v/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1gm8/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1gpk/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1hnn/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1hpo/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1hvy/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1hwi/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1ia1/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1ig3/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1jla/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1k3u/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1ke5/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1l2s/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1kzk/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1l7f/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1lpz/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1m2z/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1meh/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1mmv/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1n1m/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1n2j/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1n46/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1nav/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1of1/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1of6/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1opk/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1owe/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1oyt/");

//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1p62/");
//	string parameterFile("/home/stan/workspace/aDock_glsl_0.0.1/dockResult/1pmn/");

//	DeDock dock(parameterFile);
//	dock.init();
//	dock.run();

	return 0;
}

void
test_MODEDock( int argc, char** argv ){
    string par(argv[1]);
    cout<<"par:"<<par<<endl;
    DeDock dock("");
    if( !dock.init()){
        cout<<"done!!!"<<endl;
        return;
    }
    if(0){
        vector<Atom> receptorAtoms = dock.getBindingSiteAtoms();
        double refHbEnergy = dock.getHBondEnergy(receptorAtoms, dock.getReferLigand());
        cout<<"refHbEnergy:"<<refHbEnergy<<endl;
        return;
    }
    if( !dock.initFlexLigandPopulation() ){
        return;
    }
    dock.deFlexLigandDock();
    cout<<"done!!!"<<endl;
}

void
test_innerVDW2(){

    string fileNameLig1 = "gold_soln_corina_ligand_m1_1.mol2";
    string fileNameLig2 = "gold_soln_corina_ligand_m1_2.mol2";
    string fileNameLig3 = "gold_soln_corina_ligand_m1_3.mol2";
    string fileNameLig4 = "gold_soln_corina_ligand_m1_4.mol2";
    string fileNameLig5 = "gold_soln_corina_ligand_m1_5.mol2";
    string fileNameLig6 = "gold_soln_corina_ligand_m1_6.mol2";
    string fileNameLig7 = "gold_soln_corina_ligand_m1_7.mol2";
    string fileNameLig8 = "gold_soln_corina_ligand_m1_8.mol2";
    string fileNameLig9 = "gold_soln_corina_ligand_m1_9.mol2";
    string fileNameLig10 = "gold_soln_corina_ligand_m1_10.mol2";

    Molecular lig1 = readMolecularFile( fileNameLig1 ).front();
    Molecular lig2 = readMolecularFile( fileNameLig2 ).front();
    Molecular lig3 = readMolecularFile( fileNameLig3 ).front();
    Molecular lig4 = readMolecularFile( fileNameLig4 ).front();
    Molecular lig5 = readMolecularFile( fileNameLig5 ).front();
    Molecular lig6 = readMolecularFile( fileNameLig6 ).front();
    Molecular lig7 = readMolecularFile( fileNameLig7 ).front();
    Molecular lig8 = readMolecularFile( fileNameLig8 ).front();
    Molecular lig9 = readMolecularFile( fileNameLig9 ).front();
    Molecular lig10 = readMolecularFile( fileNameLig10 ).front();

    string goldParFile("ff.dat");
    Sybyl gold(goldParFile);

    float vdw1 = gold.computeVDW( lig1 );
    float vdw2 = gold.computeVDW( lig2 );
    float vdw3 = gold.computeVDW( lig3 );
    float vdw4 = gold.computeVDW( lig4 );
    float vdw5 = gold.computeVDW( lig5 );
    float vdw6 = gold.computeVDW( lig6 );
    float vdw7 = gold.computeVDW( lig7 );
    float vdw8 = gold.computeVDW( lig8 );
    float vdw9 = gold.computeVDW( lig9 );
    float vdw10 = gold.computeVDW( lig10 );
    cout<<vdw1<<endl;
    cout<<vdw2<<endl;
    cout<<vdw3<<endl;
    cout<<vdw4<<endl;
    cout<<vdw5<<endl;
    cout<<vdw6<<endl;
    cout<<vdw7<<endl;
    cout<<vdw8<<endl;
    cout<<vdw9<<endl;
    cout<<vdw10<<endl;
}

void
test_innerVDW(){
    string fileNameLig = "Frog.mol2";
    vector<Molecular> ligVec = readMolecularFile( fileNameLig );
    for( size_t i=0; i<504; i++ ){
        Molecular lig = ligVec[i];

        string goldParFile("ff.dat");
        Sybyl gold(goldParFile);
        float vdw = gold.computeVDW( lig );
        cout<<" "<<vdw<<endl;
    }

    return;
}

void
test_obenergy(){
    ForceFieldGAFF gaff;
    gaff.ParseParamFile();

//    string fileName = "corina_ligand_bcc_gaff.mol2";
    string fileName1 ="gaff_1.mol2";
    string fileName2 ="gaff_2.mol2";
    string fileName3 ="gaff_3.mol2";
    string fileName4 ="gaff_4.mol2";
    string fileName5 ="gaff_5.mol2";
    string fileName6 ="gaff_6.mol2";
    Molecular lig1 = readMolecularFile( fileName1 ).front();
    Molecular lig2 = readMolecularFile( fileName2 ).front();
    Molecular lig3 = readMolecularFile( fileName3 ).front();
    Molecular lig4 = readMolecularFile( fileName4 ).front();
    Molecular lig5 = readMolecularFile( fileName5 ).front();
    Molecular lig6 = readMolecularFile( fileName6 ).front();
    gaff.calculateVDW(lig1);
    gaff.calculateVDW(lig2);
    gaff.calculateVDW(lig3);
    gaff.calculateVDW(lig4);
    gaff.calculateVDW(lig5);
    gaff.calculateVDW(lig6);
    /*
    ifstream iff;
    iff.open("junk");
    if( !iff ){
        return;
    }
    string oneLine;
    vector<string> strVec;
    float val = 0;
    int i=0;
    while( getline(iff, oneLine) ){

        strVec = tokenBySpace( oneLine );
        if( strVec.size()>5 &&
            strVec[0] == "TOTAL" &&
            strVec[1] == "VAN" &&
            strVec[2] == "DER" &&
            strVec[3] == "WAALS" &&
            strVec[4] == "ENERGY" &&
            strVec[5] == "=" ){
            val = atof( strVec[6].c_str() );
            cout<<val<<endl;
        }
    }
    */
}
