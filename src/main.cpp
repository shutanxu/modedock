/*
 * main.cpp
 *
 *  Created on: Aug 25, 2014
 *      Author: stan
 */

#include <QApplication>

#include<iostream>

#include"win/window.h"

#include"test/clusterHelix.h"
#include"test/testGold.h"
#include"test/DEtest.h"
#include"test/testHBond.h"
#include"test/testDock.h"
#include"test/testMarchingCubes.h"
#include"test/testisomorphism.h"
#include"test/testHBond.h"

using namespace std;

int
main(int argc, char** argv){

//  testIsomorphism(argc, argv);
//  test_innerVDW();
//  test_obenergy();
//	testDE();
//	testGoldVDW( argc, argv );
//  testHBond(argc, argv);
//  testDock(argc, argv);
//	test_DE_flexLigandDock();
//  test_MODEDock(argc, argv);
//    clusterFragments();
//    clusterHelix2();
//    return 0;


	if(0){
        clusterHelix();
	}

//    cout<<argc<<" "<<argv[1]<<endl;
    if( argc == 2 /*&& argv[1]=="dock.par"*/ ){
        test_MODEDock(argc, argv);
    }else if(1){
		QApplication application(argc, argv);
		MainWindow *window = new MainWindow(argc, argv);
		window->show();
		return application.exec();
	}
}
