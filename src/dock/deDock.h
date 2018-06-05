/*
 * de.h
 *
 *  Created on: Nov 4, 2014
 *      Author: stan
 */

#ifndef DE_H_
#define DE_H_

/*
 * PDE: A Pareto–frontier Differential Evolution Approach for
   Multi-objective Optimization Problems

   Hussein A. Abbass
 */

#include<QGLWidget>
#include<QtOpenGL/QGLFunctions>
#include<QtOpenGL/QGLShaderProgram>
#include <QtWidgets>
#include<QThread>
#include<GL/gl.h>
#include<string>
#include<fstream>
#include<iostream>

#include"../mol/molecular.h"
#include"../algorithm/graph.h"
#include"../forcefield/sybyl.h"

using namespace std;
using namespace Qt;

class DockPar{
public:
	DockPar(const string name){ readFile(name); }
	bool  readFile(const string fileName);
	string molFileName;
	string initLigName;
	string referLigName;
	string torFileName;
	string sybylFileName;
	Coord  bindCenter;
	float  bindRadius;
	int    ligCenterAt;
	int    popSize;
	int    generation;
	int    countMax;
	float  CR;           //crossover rate
};

class DeDock:/*public QWidget, */public QThread{
	Q_OBJECT
public:

	DeDock(){}
	DeDock(const string dir):parDir(dir){}

	~DeDock(){};
    bool                       init();
	bool  readParFile();
	bool  initFlexLigandPopulation();
	bool  deFlexLigandDock();
    vector<Molecular> getFinalPop(){ return finalPop; }
    vector<Atom> getBindingSiteAtoms(){ return getBindSiteAtoms(mol, bindCenter, bindRadius); }
    Molecular getReferLigand(){ return refLig; }
signals:
    void emitString(QString);
    void finishRun();
protected:
    void run();

private:
	bool                       shuffleLigand();
	bool                       sampleRotableBonds();
	void                       sampleRotableBonds( const Molecular& ligand,
			                    const size_t bdIndex,
			                    vector<float>& angleVec,
			                    vector<vector<float> >& allAngleVec,
			                    vector<Molecular>& ligRotateVec );
	Molecular                   getMolFromPosVector(  const Molecular& mol,
			                        const int centerAtomIndex,
			                        const vector<float>& posVec,
                                    const vector<Bond> rotableBondVec,
                                    bool computeRotableBonds = true );
	bool                        ligandInBindingSite( const Molecular& mol,
			                        const Coord& bindCenter,
			                        const float& bindRadius );
	vector<Atom>                getBindSiteAtoms( const Molecular& mol,
			                        const Coord& bindingSite,
			                        const float& radius );
	vector<Bond>                get_bfsRotableBonds( const Molecular& ligand,
			                        const int& centerAtomIndex );
	double                      generateGaussianNoise(const double &variance);
	bool                        getMolFromRotateBonds(const Molecular& mol,
			                        const vector<Bond>& bfsRotableBonds,
		                            const vector<float>& rotateAngleVec,
		                            Molecular& output );
	Molecular                   getTranslateRotateMol( const Molecular& mol,
			                        const int centerAtomIndex,
			                        const vector<float>& posVec );

public:
    double
    getHBondEnergy( const Atom& donor, const Atom& acceptor );
    double
	getHBondEnergy( const Atom& donor, const vector<Atom>& acceptorAtomVec );

	double
	getHBondEnergy( const vector<Atom>& donorVec, const vector<Atom>& acceptorAtomVec );

	double
	getHBondEnergy( const vector<Atom>& receptorAtomVec, const Molecular& ligand );
	bool
	dominatesHbVDW( Sybyl sybyl, const vector<Atom> receptorAtomVec, const Molecular& lig1,
			const Molecular& lig2, const bool& flag=false);

	bool
	dominatesVDW(Sybyl sybyl, const vector<Atom> receptorAtomVec, const Molecular& lig1,
			const Molecular& lig2, const bool& flag=false);

	bool
	dominatesESP(Sybyl sybyl, const vector<Atom> receptorAtomVec, const Molecular& lig1,
			const Molecular& lig2, const bool& flag=false);

private:
	string                     molFileName;
	string                     initLigName;
	string                     referLigName;
	string                     torFileName;
	string                     sybylFileName;
	Coord                      bindCenter;
	float                      bindRadius;
	int                        ligCenterAt;
	int                        popSize;
	int                        generation;
	int                        countMax;
	float                      CR;           //crossover rate

private:
	string                     parDir;
	vector<vector<float> >     posPos;  // coding for each ligand
	vector<Molecular>          population;
    vector<Molecular>          finalPop;
	Molecular                  mol;
	Molecular                  ligand;
	int                        ligCenterAtomIndex;
	Molecular                  refLig;
	vector<Bond>               rotableBonds;
	Sybyl                      sybyl;
};

#endif /* DE_H_ */
