/*
 * goldLigand.h
 *
 *  Created on: Sep 22, 2013
 *      Author: stan
 */

#ifndef GOLDLIGAND_H_
#define GOLDLIGAND_H_

#include"../mol/molecular.h"

using namespace std;

class GoldLigand:public Molecular{

public:
	GoldLigand(){}
	GoldLigand( const string& file ):Molecular( file ){ readGoldScore( file ); }
	GoldLigand( const float& fit, const float& inVdw, const float& inTor, const float& inHb, const float& inCor,
			const float& inEner, const float& exHb, const float& exVdw, const float& rmsd): fitness( fit ), inter_vdw(inVdw),
			inter_tors( inTor ), inter_hb( inHb ), inter_corre( inCor ), inter_energy( inEner ), ext_hb( exHb ), ext_vdw( exVdw ), rmsd_ref( rmsd ){};
	~GoldLigand();
	int			readGoldScore( const string& fileName );
	float		get_fitness()const{ return fitness; }
	float		get_inter_vdw()const{ return inter_vdw; }
	float		get_inter_tors()const{ return inter_tors; }
	float		get_inter_hb()const{ return inter_hb; }
	float		get_inter_corre()const{ return inter_corre; }
	float		get_inter_energy()const{ return inter_energy; }
	float		get_ext_hb()const{ return ext_hb; }
	float		get_ext_vdw()const{ return ext_vdw; }
	float		get_rmsd_ref()const{ return rmsd_ref; }

private:
	int			readMol2( const string& fileName );
	float		fitness;
	float		inter_vdw;
	float		inter_tors;
	float		inter_hb;
	float		inter_corre;
	float		inter_energy;
	float		ext_hb;
	float		ext_vdw;
	float		rmsd_ref;
};

vector<GoldLigand>
read_goldRescoreFile( const string& fileName );

int
printGoldResult( const vector<GoldLigand>& goldLigandVec );

int
printGoldHeader();

int
print( const vector<GoldLigand>& goldLigandVec,
		const vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix );

int
print( const vector<GoldLigand>& goldLigandVec,
		const vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix,
		const	string& fileName );

int
print( const vector<GoldLigand>& goldLigandVec,
		const vector< vector<size_t> >& clusterVec,
		const	vector< vector<float> >& distMatrix,
		const	vector<float>& overlayDisVec,
		const	string& fileName );

#endif /* GOLDLIGAND_H_ */
