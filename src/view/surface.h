/*
 * surface.h
 *
 *  Created on: Mar 12, 2014
 *      Author: stan
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#include<algorithm>
#include <cmath>

#include"../mol/atom.h"

using namespace std;

const double  PI=3.14159;

class	surfaceEdge{
public:
	surfaceEdge( const Atom& a1, const Atom& a2 ):first(a1), end(a2){}
	surfaceEdge( const Atom& a1, const Atom& a2, const Coord& c ):first(a1), end(a2), probe(c){}
	Atom	get_firstAtom()const{return first;}
	Atom	get_secAtom()const{return end;}
	Coord   get_coord()const{ return probe; }
	void	print()const{ cout<<first.get_index()<<" "<<first.get_name()
						  <<" - "<<end.get_index()<<end.get_name()<<endl; }
private:
	Atom	      first;
	Atom	      end;
	Coord		  probe;
};

bool
findEdge( const vector<surfaceEdge>& edgeVec, const surfaceEdge& edge );

inline bool
operator==( const surfaceEdge& eg1, const surfaceEdge& eg2 ){

    Coord c1 = eg1.get_coord();
    Coord c2 = eg2.get_coord();

    if( c1 == c2 ){
        Atom a11 = eg1.get_firstAtom();
        Atom a12 = eg1.get_secAtom();
        Atom a21 = eg2.get_firstAtom();
        Atom a22 = eg2.get_secAtom();

        if( a11 == a21 && a12 == a22 ){
            return true;
        }else if( a11 == a22 && a12 == a21 ){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

//--------------------------------------------------------------------------------------------------------------------

class reducedSurface{
public:
	reducedSurface( const Atom& a1, const Atom& a2, const Atom& a3,
			const Coord& probe ):at1(a1), at2(a2), at3(a3){
		probeCoordVec.push_back( probe );

        index1 = a1.get_index();
        index2 = a2.get_index();
        index3 = a3.get_index();
	}
	reducedSurface( const Atom& a1, const Atom& a2, const Atom& a3,
			const vector<Coord>& probeVec ):at1(a1), at2(a2), at3(a3),
			 probeCoordVec(probeVec){
//		computeSurfaceNormal();
	}

	Atom	        get_at1()const{ return at1; }
	Atom	        get_at2()const{ return at2; }
	Atom	        get_at3()const{ return at3; }
	Coord	        get_normal()const{ return normal; }
	Coord	        get_probeCoord()const{ return probeCoordVec.front(); }
	vector<Coord>   get_probeCoordVec()const{ return probeCoordVec; }
	void	        print();

public:
    size_t          index1;    //index of at1
    size_t          index2;
    size_t          index3;

private:
	Atom	at1;
	Atom	at2;
	Atom	at3;
	Coord	normal;					// vector normal to the face
	vector<Coord> probeCoordVec;
//	Coord	probeCoord; 			// coordinates of the fixed position of the probe center
};

class sphereSurface{
public:
	sphereSurface( const Atom& a1, const Atom& a2, const Atom& a3,
			const Coord& p, const float& pR );

	void           set_arch12( const vector<Coord>& c){ arch_12 = c; }
	void           set_arch13( const vector<Coord>& c){ arch_13 = c; }
	void           set_arch23( const vector<Coord>& c){ arch_23 = c; }

	Atom           get_at1()const{ return at1; }
	Atom           get_at2()const{ return at2; }
	Atom           get_at3()const{ return at3; }
	Coord          get_probe()const{ return probe; }
	float          get_probeRadius()const{ return probeRadius; }
	vector<Coord>  get_arch12()const{ return arch_12; }
	vector<Coord>  get_arch13()const{ return arch_13; }
	vector<Coord>  get_arch23()const{ return arch_23; }
	void           print();
	bool           crash;
private:
	Atom           at1;
	Atom           at2;
	Atom           at3;
	Coord          probe;
	float          probeRadius;
	vector<Coord>  arch_12;
	vector<Coord>  arch_13;
	vector<Coord>  arch_23;

};

/**         at1
 *         /   \
 *        /  A  \
 *     at3 ----- at4
 *        \  B   /
 *         \    /
 *          at2
 *
 *  at3, at4 are the shared edge
 */
class reducedSurface_connect{
public:
	reducedSurface_connect(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4,
			const Coord& a, const Coord& b):at1(a1),at2(a2),at3(a3),at4(a4),probeCoordA(a),
			probeCoordB(b){
		probeCrash=false;
		probeRadius = getCoordDis( at3.get_coord(), probeCoordA )-at3.get_vdwRadius();

	}

	Atom			get_at1()const{ return at1; }
	Atom			get_at2()const{ return at2; }
	Atom			get_at3()const{ return at3; }
	Atom			get_at4()const{ return at4; }
	Coord			get_probeCoordA()const{ return probeCoordA; }
	Coord			get_probeCoordB()const{ return probeCoordB; }
	vector<Coord>	get_coordVec()const{ return coordVec; }
//	void            set_probeCrash(){ probeCrash = true; }
	void			set_crashCoordA( const vector<Coord> c ){ crashCoordVecA = c; }
	void			set_crashCoordB( const vector<Coord> c ){ crashCoordVecB = c; }
	vector<Coord>   get_crashCoordA()const{ return crashCoordVecA; }
	vector<Coord>   get_crashCoordB()const{ return crashCoordVecB; }
	double          get_probeRadius()const{ return probeRadius; }
//	bool            probeCrash(){ if( probeCrash==true){return true;}else{return false;}}

	bool            computeReducedSurfaceArch(const vector<Atom>& neighbors);
	vector<Coord>   getCrashProbeCoordConcave(const Coord& c1, const Coord& c2);
	vector<Coord>   getCrashProbeCoordConvex(const Coord& c1, const Coord& c2);
	void   			print()const;

	bool            probeCrash;
private:
	Atom	      	at1;
	Atom			at2;
	Atom			at3;
	Atom			at4;
	Coord			probeCoordA; 	// coordinates of the fixed position of the probe center
	Coord			probeCoordB; 	// coordinates of the fixed position of the probe center
	vector<Coord>   coordVec;       // rolling trajectory from probeCoordA to probeCoordB
	vector<Coord>   crashCoordVecA;
	vector<Coord>   crashCoordVecB;
	double          probeRadius;
};


vector<sphereSurface>
get_sphereSurface( const vector<reducedSurface>& rsVec,
		const vector<reducedSurface_connect>& rscVec );
inline bool
operator==( const reducedSurface& rs1, const reducedSurface& rs2 ){
    size_t a11 = rs1.index1;
    size_t a12 = rs1.index2;
    size_t a13 = rs1.index3;
    size_t a21 = rs2.index1;
    size_t a22 = rs2.index2;
    size_t a23 = rs2.index3;

    if( rs1.get_probeCoord() == rs2.get_probeCoord() ){
        if( a11 == a21 && a12 == a22 && a13 == a23 ){
            return true;
        }else if( a11 == a21 && a12 == a23 && a13 == a22){
            return true;
        }else if( a11 == a22 && a12 == a21 && a13 == a23 ){
            return true;
        }else if( a11 == a22 && a12 == a23 && a13 == a21 ){
            return true;
        }else if( a11 == a23 && a12 == a21 && a13 == a22 ){
            return true;
        }else if( a11 == a23 && a12 == a22 && a13 == a21 ){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }

//    Atom a11 = rs1.get_at1();
//    Atom a12 = rs1.get_at2();
//    Atom a13 = rs1.get_at3();
//    Atom a21 = rs2.get_at1();
//    Atom a22 = rs2.get_at2();
//    Atom a23 = rs2.get_at3();

//    if( a11 == a21 && a12 == a22 && a13 == a23 ){
//        return true;
//    }else if( a11 == a21 && a12 == a23 && a13 == a22){
//        return true;
//    }else if( a11 == a22 && a12 == a21 && a13 == a23 ){
//        return true;
//    }else if( a11 == a22 && a12 == a23 && a13 == a21 ){
//        return true;
//    }else if( a11 == a23 && a12 == a21 && a13 == a22 ){
//        return true;
//    }else if( a11 == a23 && a12 == a22 && a13 == a21 ){
//        return true;
//    }else{
//        return false;
//    }

}

inline bool
operator!=( const reducedSurface& rs1, const reducedSurface& rs2 ){
	return		!( rs1 == rs2 );
}

inline bool
operator==( const sphereSurface& s1, const sphereSurface& s2 ){
    if( s1.get_probe() == s2.get_probe() ){
        Atom a11 = s1.get_at1();
        Atom a12 = s1.get_at2();
        Atom a13 = s1.get_at3();
        Atom a21 = s2.get_at1();
        Atom a22 = s2.get_at2();
        Atom a23 = s2.get_at3();
        if( ( a11 == a21 || a11 == a22 || a11 == a23 ) &&
            ( a12 == a21 || a12 == a22 || a12 == a23 ) &&
            ( a13 == a21 || a13 == a22 || a13 == a23 ) ){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

inline bool
operator!=( const sphereSurface& s1, const sphereSurface& s2 ){
	return !( s1 == s2 );
}

inline bool
operator==( const  reducedSurface_connect& s1, const  reducedSurface_connect& s2 ){
    Coord c11 = s1.get_probeCoordA();
    Coord c12 = s1.get_probeCoordB();
    Coord c21 = s2.get_probeCoordA();
    Coord c22 = s2.get_probeCoordB();

    if( ( c11 == c21 || c11 == c22 ) &&
        ( c12 == c21 || c12 == c22 ) ){
        Atom a11 = s1.get_at1();
        Atom a12 = s1.get_at2();
        Atom a13 = s1.get_at3();
        Atom a14 = s1.get_at4();
        Atom a21 = s2.get_at1();
        Atom a22 = s2.get_at2();
        Atom a23 = s2.get_at3();
        Atom a24 = s2.get_at4();

        if( ( a11 == a21 || a11 == a22 ) &&
            ( a12 == a21 || a12 == a22 ) &&
            ( a13 == a23 || a13 == a24 ) &&
            ( a14 == a23 || a14 == a24 ) ){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

inline bool
operator!=( const  reducedSurface_connect& s1, const  reducedSurface_connect& s2 ){
	return !( s1 == s2 );
}

double
getAngle( const Coord& init, const Coord& p1, const Coord& p2  );

Coord
rotateAroundAxis(const Coord& input, const Coord& axis, const double& angle);

Coord
getTriangleVerticalPoint( const Coord& a, const Coord& b, const Coord& p );

bool
atomExistInRS( const Atom& at, const vector<reducedSurface>& rsVec );

bool
probeCollision( const Coord& probeCoord, const double probeRadius, const vector<Atom>& atomVec );

bool
probeMayExist( const Atom at1, const Atom at2, const Atom at3,  const double probeRadius=1.5 );

reducedSurface
getFirstReducedSurface( const vector<Atom>& atomVec , const double probeRadius=1.5 );

pair<Coord, Coord>
getFixedProbe( const Atom at1, const Atom at2, const Atom at3, const double probeRadius=1.5 );

pair<Coord, Coord>
getFixedProbe( const Coord co1, const Coord co2, const Coord co3, const double probeRadius=1.5 );

/*
Coord
getFixedProbe( const Atom at1, const Atom at2, const Atom at3,
		const vector<Atom>& atomVec, const double probeRadius=1.5  );
*/
bool
getFixedProbe( const Atom at1, const Atom at2, const Atom at3,
		const vector<Atom>& atomVec, const double probeRadius, Coord& probe  );

bool
getFixedProbe( const Atom at1, const Atom at2, const Atom at3,
		const vector<Atom>& atomVec, const double probeRadius, vector<Coord>& probeVec );

vector<reducedSurface>
fixReducedSurface( const vector<reducedSurface>& rsVec );

void
getReducedSurface(  const vector<Atom>& atomVec,
                    const double probeRadius,
		vector<reducedSurface>& rsVec,
		vector<sphereSurface>& ssVec,
		vector<reducedSurface_connect>& rscVec );

void
getReducedResSurface(  const vector<Atom>& atomVec,
                    const double probeRadius,
        vector<reducedSurface>& rsVec,
        vector<sphereSurface>& ssVec,
        vector<reducedSurface_connect>& rscVec );

bool
updateSphereSurface( vector<reducedSurface_connect>& rscVec,
		vector<sphereSurface>& ssVec );

vector<reducedSurface>
getReducedSurface( const vector<Atom>& atomVec , const double probeRadius=1.5 );

Coord
getSurfaceNormal( const Atom at1, const Atom at2, const Atom at3,  Coord& probe );

vector<Coord>
getArchPoints( const Coord& c1, const Coord& c2, const Coord& p,
		const double& probeRadius, const int& numOfOut );

/*
 * probe rotate around the axis of a1-a2 from p1 to p2
 */
vector<Coord>
getReducedSurfaceArch( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut );

vector<Coord>
getReducedSurfaceArch( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut, bool concave );

vector<Coord>
getReducedSurfaceArch( const reducedSurface_connect& rsc, const int& numOut );

void
get_SES_arch( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut,
		vector<Coord>& out1, vector<Coord>& out2);

void
get_SES_toric_reentrant_face( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut,
		const bool& concave, vector<Coord>& coordVec,
		vector<Coord>& normalVec);

void
get_SES_toric_reentrant_face( reducedSurface_connect& rsc,
		const int& numOut, vector<Coord>& coordVec, vector<Coord>& normalVec );

vector<Coord>
get_sphere_arch( const Coord& center, const double radius,
		const Coord& p1, const Coord& p2, const int numOut );

void
get_SES_spheric_reentrant_face(const sphereSurface& ss,
		vector<Coord>& coordVec, vector<Coord>& normalVec);

void
get_SES_spheric_reentrant_face(  const Atom& a1, const Atom& a2, const Atom& a3,
		const Coord& probe, const int& numOut,
		vector<Coord>& coordVec, vector<Coord>& normalVec);

vector<reducedSurface_connect>
get_reducedSurface_connect( const vector<reducedSurface>& rsVec );

vector<reducedSurface_connect>
get_reducedSurface_connect( const vector<reducedSurface>& rsVec );


#endif /* SURFACE_H_ */
