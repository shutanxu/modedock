#ifndef GAFF_H
#define GAFF_H

#include<iostream>
#include<stdlib.h>
#include<vector>
#include<string>

#include"../algorithm/tokenBySpace.h"
#include"../algorithm/graph.h"
#include"../mol/atom.h"
#include"../mol/bond.h"
#include"../mol/molecular.h"
#include"forcefield.h"

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

using namespace std;

class ForceFieldGAFF{
public:
    bool ParseParamFile();
    float calculateVDW( const Atom a1, const Atom a2 );
    float calculateVDW( const Molecular mol);
private:
    vector<ForceFieldParameter> ffpropparams;
    vector<ForceFieldParameter> ffbondparams;
    vector<ForceFieldParameter> ffangleparams;
    vector<ForceFieldParameter> fftorsionparams;
    vector<ForceFieldParameter> ffoopparams;
    vector<ForceFieldParameter> ffhbondparams;
    vector<ForceFieldParameter> ffvdwparams;
    vector<ForceFieldParameter> ffchargeparams;
};

#endif // GAFF_H
