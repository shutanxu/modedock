/*
 * testSS.h
 *
 *  Created on: May 10, 2014
 *      Author: stan
 */

#ifndef TESTSS_H_
#define TESTSS_H_

#include <qapplication.h>
#include<stdlib.h>
#include<stdio.h>
#include"../view/molViewer.h"

using namespace std;

int
drawSecondaryStructure( int argc, char** argv ){
	QApplication application(argc,argv);

	string						fileName = "/home/stan/Desktop/protein.mol2";
	string						dsspFile = "/home/stan/workspace/dssp-2.2.1/temp.dssp";
	Molecular					molM( fileName );
//	SStructure				ss( dsspFile, fileName );
//	SStructure				ss( molM );

	MolecularViewer					mol;
//	mol.set_SecStructure( ss );
	mol.show_lines();
//	mol.show_ballStick();
//	mol.show_labelAtomIndex();
//	mol.show_labelResidue();
//	mol.show_secondaryStructure();
//	mol.show_vdwBall();
	mol.show();
	return application.exec();
}

#endif /* TESTSS_H_ */
