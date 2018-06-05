#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include<iostream>
#include<stdlib.h>
#include<vector>

#include"../algorithm/tokenBySpace.h"
#include"../mol/atom.h"
#include"../mol/bond.h"
#include"../mol/molecular.h"

#define KCAL_TO_KJ	4.1868

using namespace std;

/*
 * derived from openbabel forcefield.h
 * @brief The ForceFieldParameter class
 */
class ForceFieldParameter{
public:
    int a,b,c,d;          //! Used to store integer atom types
    string sa,sb,sc,sd;   //! used to store string atom types
    vector<int> ipar;     //! Used to store integer type parameters (bondtypes, multiplicity, ...)
    vector<float> dpar;   //! Used to store double type parameters (force constants, bond lengths, angles, ...)

    ForceFieldParameter& operator=(const ForceFieldParameter& p){
        if(this != &p){
            a=p.a;
            b=p.b;
            c=p.c;
            d=p.d;
            sa=p.sa;
            sb=p.sb;
            sc=p.sc;
            sd=p.sd;
            ipar=p.ipar;
            dpar=p.dpar;
        }
        return *this;
    }

    void clear(){ a=b=c=d=0; sa=sb=sc=sd=""; ipar.clear(); dpar.clear(); }
};

float
computeInnerVDW( const vector<Atom> atVec, const vector<Bond> bdVec );

#endif // FORCEFIELD_H
