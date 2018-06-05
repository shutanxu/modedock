/*
 * testDock.h
 *
 *  Created on: Oct 4, 2014
 *      Author: stan
 */

#ifndef TESTDOCK_H_
#define TESTDOCK_H_

#include <qapplication.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<cstring>
#include<sstream>
#include<stdlib.h>
#include<omp.h>

#include"../mol/molecular.h"
#include"../forcefield/sybyl.h"
#include"goldLigand.h"

using namespace std;

int
test_DE_Dock();

int
test_DE_flexLigandDock();

int
testDock( int argc, char** argv );

void
test_MODEDock( int argc, char** agrv );

void
test_innerVDW();

void
test_obenergy();

void
test_innerVDW2();

#endif /* TESTDOCK_H_ */
