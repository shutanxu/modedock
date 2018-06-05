/*
 * autoDock.h
 *
 *  Created on: Jan 3, 2014
 *      Author: stan
 */

#ifndef AUTODOCK_H_
#define AUTODOCK_H_

#include<iostream>
#include<stdlib.h>
#include<algorithm>
#include<vector>

#include"../mol/atom.h"
#include"../mol/molecular.h"

using namespace std;

//static AD_ATOM_PARM AD_ATOM_PARM_NULL;


struct
AD_ATOM_PARM{
	AD_ATOM_PARM(){}
	AD_ATOM_PARM( const vector<string>& strVec );
	string				atomType;
	float				Rii;
	float				espii;
	float				vol;
	float				solpar;
	float				Rij_hb;
	float				espij_hb;
	float				hbond;
	int					rec_index;
	int					map_index;
	int					bond_index;
};


struct
AD_PARM{
	void				readParamFile( const string& fileName );

	float										FE_coeff_vdw;
	float										FE_coeff_hbond;
	float										FE_coeff_estat;
	float										FE_coeff_desolv;
	float										FE_coeff_tors;
	vector<AD_ATOM_PARM>		atomParmVec;
};


/**
 * DOCKED: USER                              			x       y       z     vdW  Elec       q    Type
DOCKED: USER                                 	_______ _______ _______ _____ _____    ______ ____
DOCKED: ROOT
DOCKED: ATOM      1  C9  LIG d   1       2.829  11.960  41.473 -0.37 -0.09    -0.038 A
DOCKED: ATOM      2  C10 LIG d   1       1.489  12.151  41.191 -0.45 +0.02    +0.007 A
 */
class
AutoDock_ATOM:public Atom{
public:
	AutoDock_ATOM(){}
	AutoDock_ATOM(   const size_t& id, const string& na, const string& ty, const Coord& co,
			const float& vd, const float& el, const float& qq ):Atom( id, na, ty, co ), vdw(vd), elec(el), q(qq){ set_charge( qq ); }

	float							get_vdw()const { return vdw; }
	float							get_elec()const{ return elec ; }
private:
	float							vdw;
	float							elec;
	float							q;
};

struct
AutoDock_Result{
	AutoDock_Result( const size_t& ind, const float& bE, const float& fi, const float& vhd, const float& es,
			const float& ti, const float& tor, const float& ub, const vector<AutoDock_ATOM>& at ):ligand_index(ind),
			bindEnergy(bE), final_inter(fi), vdw_hbond_desolv(vhd), electrostatic(es), total_internal(ti), torsional(tor),
			unbound(ub), atomVec(at){}
	size_t										ligand_index ;
	float										bindEnergy;
	float										final_inter ;
	float										vdw_hbond_desolv;
	float										electrostatic;
	float										total_internal;
	float										torsional ;
	float										unbound ;
	vector<AutoDock_ATOM>		atomVec;
};

class
AutoDock{
public:
	AutoDock(){}
	AutoDock( const string& fileName ):ffpFileName( fileName ){ atom_parm.readParamFile( fileName ); }
	float					computeVDW_cut(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
	float					computeVDW_cut(   const Atom& at1, const Atom& at2  );
	float					computeVDW( const Molecular& molA, const Molecular& molB );
	float					computeVDW(  const Atom& at1, const Atom& at2 );
	float					computeVDW(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
	float					computeHBondEnergy(  const Molecular& molA, const Molecular& molB );
	float					computeHBondEnergy(  const Atom& at1, const Atom& at2 );
	float					computeHBondEnergy(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
	float					computeElectrostatic(  const Atom& at1, const Atom& at2 );
	float					computeElectrostatic(  const vector<Atom>& atVec1, const vector<Atom>& atVec2 );
private:
	AD_ATOM_PARM			getAtomParm( const Atom& at );
private:
	string							ffpFileName;
	AD_PARM						atom_parm;
};

vector<AutoDock_Result>
readFile_dlg(  const string& fileName );

#endif /* AUTODOCK_H_ */
