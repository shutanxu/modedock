/*
 * testHBond.h
 *
 *  Created on: Nov 24, 2013
 *      Author: stan
 */

#ifndef TESTHBOND_H_
#define TESTHBOND_H_

#include <qapplication.h>
#include<stdlib.h>
#include<stdio.h>
#include<string>

//#include"../molecular/overlay.h"
#include"../mol/molecular.h"
#include"../forcefield/sybyl.h"
#include"goldLigand.h"

using namespace std;

int
testHBond4( int argc, char** argv );

int
testHBond2( int argc, char** argv );

int
testHBond( int argc, char** argv );

#endif /* TESTHBOND_H_ */
