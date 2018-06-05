/*
 * atom.cpp
 *
 *  Created on: Sep 21, 2013
 *      Author: stan
 */

#include<string>
#include"atom.h"
#include"../algorithm/convert.h"

using namespace std;

//---------------------------------------------------------------------------------------------------------------------

Atom::Atom(){
	index = 0;
	name = "";
	type = "";
	Coord	co( 0,0,0);
	coord = co;
	sequ = 0;
	residueName = "";
	residueIndex = 0;
	chainName = "";
	charge = 0;
	vdwRadius = 0;
}

Atom::Atom(const vector<string>& strVec) {
	if (strVec.size() < 9) {
		cerr << "---can not initialize atom " << strVec[0] << endl;
	}

	if( strVec[0] == "ATOM" ||   strVec[0] == "HETATM" ){
		index            = atof( strVec[1].c_str() );
		name             = strVec[2];
		residueName      = strVec[3];

		residueIndex     = atoi( strVec[5].c_str() );
		coord.x          = atof( strVec[6].c_str() );
		coord.y          = atof( strVec[7].c_str() );
		coord.z          = atof( strVec[8].c_str() );
		type             = strVec.back();
	}else{
		index            = atof(strVec[0].c_str());
		name             = strVec[1].c_str() ;
		coord.x          = atof(strVec[2].c_str());
		coord.y          = atof(strVec[3].c_str());
		coord.z          = atof(strVec[4].c_str());
		type             = strVec[5].c_str() ;
		sequ             = atoi( strVec[6].c_str() );
		string           resid = strVec[7].c_str();
		residueName      = resid.substr( 0, 3 );
		if( resid.length() <= 3 ){
			residueIndex = sequ;
		}else{
			residueIndex = atoi ( resid.substr( 3  ).c_str() ) ;
		}
		charge           = atof( strVec[8].c_str() );

	}
	set_vdwRadiusDefault();
}

void
Atom::readMol2( const vector<string>& strVec ){
	index            = atof(strVec[0].c_str());
	name             = strVec[1].c_str() ;
	coord.x          = atof(strVec[2].c_str());
	coord.y          = atof(strVec[3].c_str());
	coord.z          = atof(strVec[4].c_str());
	type             = strVec[5].c_str() ;
	sequ             = atoi( strVec[6].c_str() );
	string           resid = strVec[7].c_str();
	residueName      = resid.substr( 0, 3 );
//	residueIndex     = atoi ( resid.substr( 3  ).c_str() ) ;
	if( resid.length() <= 3 ){
		residueIndex = sequ;
	}else{
		residueIndex = atoi ( resid.substr( 3  ).c_str() ) ;
	}
	charge           = atof( strVec[8].c_str() );

	set_vdwRadiusDefault();
}

Atom::Atom( const string& oneLine ){
	vector<string>	strVec = tokenBySpace( oneLine.substr(0, 6) );
	if(  strVec[0] == "ATOM" ||   strVec[0] == "HETATM" ){
		index 	 			= atof( tokenBySpace(oneLine.substr(6, 5) )[0].c_str() );
		name				= tokenBySpace( oneLine.substr(11, 5 ) )[0];
		residueName	= tokenBySpace(oneLine.substr( 17, 3 ) )[0];
//		cout<<"sub"<<oneLine.substr( 21, 1 ) <<endl;
		chainName		= tokenBySpace(oneLine.substr( 21, 1 ) )[0];
		residueIndex	= atoi( tokenBySpace( oneLine.substr( 22, 4 ) )[0].c_str() );
		coord.x				= atof( tokenBySpace( oneLine.substr( 29, 9 ) )[0].c_str() );
		coord.y				= atof( tokenBySpace( oneLine.substr( 38, 9 ) )[0].c_str() );
		coord.z				= atof( tokenBySpace( oneLine.substr( 46, 9 ) )[0].c_str() );
//		type					= strVec.back();
		if( oneLine.length() > 77 ){
			if(   tokenBySpace( oneLine.substr( 76, 3 ) )[0] == "" ){
				type				= name.substr(0,1);
			}else{
				type				=  tokenBySpace( oneLine.substr( 76, 3 ) )[0] ;
			}
		}else{
			type				= name.substr(0,1);
		}
	}
	if( chainName == "" ){
		chainName = string("N");
	}
	if( chainName == "" ){
		chainName = string("A");
	}
	set_vdwRadiusDefault();
}

bool
Atom::isResidue(){
	for( int i=0; i<20; i++ ){
		string			s( ResidueType[i] );
		if( residueName == s ){
			return true;
		}
	}
	return false;
}

bool
Atom::isDNA(){
	for( int i=0; i<8; i++ ){
		string			s( DNAtype[i] );
		if( residueName == s ){
			return true;
		}
	}
	return false;
}

bool
Atom::isRNA(){
	for( int i=0; i<8; i++ ){
		string			s( RNAtype[i] );
		if( residueName == s ){
			return true;
		}
	}
	return false;
}

bool
Atom::isLigand(){
	if( ! isResidue() && ! isDNA() &&  ! isRNA() && residueName != "HOH" ){
		return true;
	}
	return false;
}

int Atom::print()const {

	cout.width(6);
	cout<<left<<"ATOM";
	cout.width(5);
	cout <<right<< index<<"  ";
	cout.width(4);
	cout <<left<< name;
	cout.width(3);
	cout<<right<<residueName<<" ";
	cout.width(1);
	cout<<left<<chainName ;
	cout.width(4);
	cout<<right<<residueIndex;
	cout.width(3);
	cout<<"";
	cout.width(9);
	cout <<right<< float2string(coord.x).substr( 0, 7 );
	cout.width(8);
	cout <<right<< float2string(coord.y).substr( 0, 7 );
	cout.width(8);
	cout <<right<< float2string(coord.z).substr( 0, 7 );
	cout.width(24);
	cout <<right<<type.substr(0, type.find('.'))<<endl;

}

void
Atom::set_vdwRadiusDefault(){

	if( type.find( "H") != string::npos ){
		vdwRadius = 1.20;
	}else if( type.find( "Cu") != string::npos  ){
		vdwRadius =  1.4;
	}else if( type.find( "Cl") !=string::npos  ){
		vdwRadius =  1.75;
	}else if( type.find( "C") != string::npos ){
		vdwRadius =  1.70;
	}else if( type.find( "N") != string::npos ){
		vdwRadius =  1.55;
	}else if( type.find( "O") != string::npos ){
		vdwRadius =  1.52;
	}else if( type.find( "F" ) != string::npos ){
		vdwRadius =  1.47;
	}else if( type.find("P") != string::npos ){
		vdwRadius =  1.80;
	}else if( type.find( "S") != string::npos  ){
		vdwRadius =  1.80;
	}else{
		vdwRadius =  0 ;
	}
}

bool
operator==( const Atom& at1, const Atom& at2 ){

//    return		( at1.get_index() == at2.get_index() &&
//                  at1.get_name() == at2.get_name() &&
//                  at1.get_type() == at2.get_type() &&
//                  at1.get_coord() == at2.get_coord());

	return		( at1.get_index() == at2.get_index() &&
			at1.get_name() == at2.get_name() &&
			at1.get_type() == at2.get_type() );

}

Atom&
Atom::operator =(const Atom& at){
	index = at.get_index();
	name = at.get_name();
	type = at.get_type();
	coord = at.get_coord();
	sequ = at.get_sequ();
	residueName = at.get_residueName();
	residueIndex = at.get_residueIndex();
	chainName = at.get_chainName();
	charge = at.get_charge();
	vdwRadius = at.get_vdwRadius();
	return *this;
}

vector<float>
Atom::get_color(){
	vector<float> color(3);
	if( type.at(0) == 'C' ){
		color[0] = 0;
		color[1] = 1;
		color[2] = 0;
	}else	if( type.at(0) == 'N' ){
		color[0] = 0;
		color[1] = 0;
		color[2] = 1;

	}else 	if( type.at(0) == 'O' ){
		color[0] = 1;
		color[1] = 0.1;
		color[2] = 0.1;

	}else 	if( type.at(0) == 'H' ){
		color[0] = 1;
		color[1] = 1;
		color[2] = 1;

	}else 	if( type.at(0) == 'S' ){
		color[0] = 1;
		color[1] = 1;
		color[2] = 0;

	}else if(  type.at(0) == 'F'  ){
		color[0] = 1;
		color[1] = 0.7;
		color[2] = 0.7;

	}else if( type == "FE" ){
		color[0] = 1;
		color[1] = 0.498;
		color[2] = 0;

	}else{
		color[0] = 1;
		color[1] = 0.75;
		color[2] = 0.79;

	}
	return color;
}

//------------------------------------------------------------------------
Atom searchAtom(const vector<Atom>& atomVec, const size_t& index) {
	for (size_t i = 0; i < atomVec.size(); i++) {
		if (atomVec[i].get_index() == index) {
			return	atomVec[i];
		}
	}
	cerr<<"--- warnning!!! can not find atom "<<index<<endl;
}

Atom searchAtom(const vector<Atom>& atomVec,
		const size_t& index, const string& name,
		const size_t& residIndex, const string& chainName) {

	for (size_t i = 0; i < atomVec.size(); i++) {
		if (	atomVec[i].get_index() == index &&
				atomVec[i].get_name() == name &&
				atomVec[i].get_residueIndex() == residIndex &&
				atomVec[i].get_chainName() == chainName ) {
//			cout<<"atom"<<endl;
//			atomVec[i].print();
			return	atomVec[i];
		}
	}
	cerr<<" can not find atom "<<index<<endl;
}

Atom
searchAtom( const vector<Atom>& atomVec, const string& atomName,
		const size_t& atomIndex, const size_t& residIndex ){
	for (size_t i = 0; i < atomVec.size(); i++) {
		if (atomVec[i].get_residueIndex() == residIndex &&
				atomIndex == atomVec[i].get_index() &&
				atomVec[i].get_name() == atomName) {
			return	atomVec[i];
		}
	}
	cerr<<" can not find atom "<<atomName<<" in resid:"<<residIndex<<endl;
	Atom		at;
	return		at;
}

vector<Atom>
searchAtom( const vector<Atom>& atomVec, const string& atomName ){
	vector<Atom> atVec;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( atomVec[i].get_name() == atomName ){
			atVec.push_back( atomVec[i] );
		}
	}
	return atVec;
}

Coord
centerOfAtoms( const vector<Atom>& atomVec ){
	float			aveX, aveY, aveZ;
	aveX = aveY = aveZ = 0;

	for( size_t i=0; i<atomVec.size(); i++ ){
		aveX += atomVec[i].get_coord().x;
		aveY += atomVec[i].get_coord().y;
		aveZ += atomVec[i].get_coord().z;
	}

	if( !atomVec.empty() ){
		aveX /= atomVec.size();
		aveY /= atomVec.size();
		aveZ /= atomVec.size();
	}
	return Coord( aveX, aveY, aveZ );
}

//bool
//atomsCrash( const Atom& at1, const Atom& at2 ){
//	string n1=at1.get_type();
//	string n2=at2.get_type();
//	string t1=n1.substr(0,n1.find("."));
//	string t2=n2.substr(0,n2.find("."));
//	float dis=atomDis(at1, at2);
//	if( t1==string("H") )
//
//	return false;
//}

bool
atomsCrash(const Atom& at1, const Atom& at2 ){
	float dis = atomDis( at1, at2 );
    if( dis < (at1.get_vdwRadius() + at2.get_vdwRadius() ) ){
		return true;
	}
	return false;
}

/**
 * atom crash
 */
bool
atomsCrash( const vector<Atom>& atVec1, const vector<Atom>& atVec2 , bool computeHydrogen  ){
	for( size_t i=0; i<atVec1.size(); i++ ){
		for( size_t j=0; j<atVec2.size(); j++ ){
            if( computeHydrogen ){
                if(atomsCrash(atVec1[i], atVec2[j])){
                    return true;
                }
            }else if( atVec1[i].get_type() != "H" &&  atVec2[j].get_type() != "H" ){
                if(atomsCrash(atVec1[i], atVec2[j])){
                    return true;
                }
            }
		}
	}
	return false;
}

bool
atomsCrash( const Atom& at, const vector<Atom>& atVec, bool computeHydrogen ){
	for( size_t i=0; i<atVec.size(); i++ ){
        if( computeHydrogen ){
            if( atomsCrash(at, atVec[i]) )
                return true;
        }else if( atVec[i].get_type() != "H" ){
            if( atomsCrash(at, atVec[i]) ){
                return true;
            }
        }
	}
	return false;
}

float
atomDis( const Atom& at1,  const Atom& at2 ){
	float			x1 = at1.get_coord().x;
	float			y1 = at1.get_coord().y;
	float			z1 = at1.get_coord().z;

	float			x2 = at2.get_coord().x;
	float			y2 = at2.get_coord().y;
	float			z2 = at2.get_coord().z;

	float			d1 = x1-x2;
	float			d2 = y1-y2;
	float			d3 = z1-z2;

	return			sqrt( d1*d1 + d2*d2 + d3*d3 );
}

float
atomVecRMSD( const vector<Atom>& at1,  const vector<Atom>& at2 ){
	float dis = 0;
	if(at1.size() != at2.size() ){
		cout<<"error!!! can not compute rmsd. different atom vector length: "
				<<at1.size()<<" "<<at2.size()<<endl;
		for( size_t i=0; i<at1.size(); i++ ){
			at1[i].print();
		}
		for( size_t i=0; i<at2.size(); i++ ){
			at2[i].print();
		}
	}
	for( size_t i=0; i<at1.size(); i++ ){
		dis += atomDis( at1[i], at2[i] );
	}
	dis /= float(at1.size());
	return dis;
}

float
atomSquareDis( const Atom& at1,  const Atom& at2 ){
	float			x1 = at1.get_coord().x;
	float			y1 = at1.get_coord().y;
	float			z1 = at1.get_coord().z;

	float			x2 = at2.get_coord().x;
	float			y2 = at2.get_coord().y;
	float			z2 = at2.get_coord().z;

	float			d1 = x1-x2;
	float			d2 = y1-y2;
	float			d3 = z1-z2;

	return			( d1*d1 + d2*d2 + d3*d3 );
}

float
atomVecSquareRMSD( const vector<Atom>& at1,  const vector<Atom>& at2 ){
	float dis = 0;
	if(at1.size() != at2.size() ){
		cout<<"error!!! can not compute rmsd. different atom vector length: "
				<<at1.size()<<" "<<at2.size()<<endl;

		cout<<"1 "<<at1.size()<<endl;
		for( size_t i=0; i<at1.size(); i++ ){
			at1[i].print();
		}
		cout<<endl;
		cout<<"2 "<<at2.size()<<endl;
		for( size_t i=0; i<at2.size(); i++ ){
			at2[i].print();
		}
	}
	for( size_t i=0; i<at1.size(); i++ ){
		dis += atomSquareDis( at1[i], at2[i] );
	}
	dis /= float(at1.size());
	return dis;
}

float
atomAngle( const Atom& dir1Start, const Atom& dir1End, const Atom& dir2Start, const Atom& dir2End ){
	Coord dir1, dir2, crossMultiply;
	dir1.x = dir1End.get_coord().x - dir1Start.get_coord().x;
	dir1.y = dir1End.get_coord().y - dir1Start.get_coord().y;
	dir1.z = dir1End.get_coord().z - dir1Start.get_coord().z;

	dir2.x = dir2End.get_coord().x - dir2Start.get_coord().x;
	dir2.y = dir2End.get_coord().y - dir2Start.get_coord().y;
	dir2.z = dir2End.get_coord().z - dir2Start.get_coord().z;

//	crossMultiply.x = 		dir1.y * dir2.z - dir2.y * dir1.z ;
//	crossMultiply.y = - ( dir1.x * dir2.z - dir2.x * dir1.z );
//	crossMultiply.z = 		dir1.x * dir2.y - dir2.x * dir1.y;

	float			dotMultiply = dir1.x * dir2.x + dir1.y*dir2.y + dir1.z*dir2.z;

	if(dir1.getDis() == 0 || dir2.getDis() == 0 ){
		return 0;
	}

//	float			angle = asin( (float)crossMultiply.getDis() / ( dir1.getDis() * dir2.getDis() ) );
	float		angle = acos( dotMultiply / (  dir1.getDis() * dir2.getDis() ) );
	angle = angle * 180.0 / 3.14159265;
	return 		angle;
}

float
atomAngle( const Atom& end1, const Atom& vertexAt, const Atom& end2 ){
	Coord dir1, dir2, crossMultiply;
	dir1.x = end1.get_coord().x - vertexAt.get_coord().x;
	dir1.y = end1.get_coord().y - vertexAt.get_coord().y;
	dir1.z = end1.get_coord().z - vertexAt.get_coord().z;

	dir2.x = end2.get_coord().x - vertexAt.get_coord().x;
	dir2.y = end2.get_coord().y - vertexAt.get_coord().y;
	dir2.z = end2.get_coord().z - vertexAt.get_coord().z;

//	crossMultiply.x = 		dir1.y * dir2.z - dir2.y * dir1.z ;
//	crossMultiply.y = - ( dir1.x * dir2.z - dir2.x * dir1.z );
//	crossMultiply.z = 		dir1.x * dir2.y - dir2.x * dir1.y;

	float			dotMultiply = dir1.x * dir2.x + dir1.y*dir2.y + dir1.z*dir2.z;

	if(dir1.getDis() == 0 || dir2.getDis() == 0 ){
		return 0;
	}

//	float			angle = asin( (float)crossMultiply.getDis() / ( dir1.getDis() * dir2.getDis() ) );
	float		angle = acos( dotMultiply / (  dir1.getDis() * dir2.getDis() ) );
	angle = angle * 180.0 / 3.14159265;
	return 		angle;
}

Atom
findClosestAtom( const  Atom& atom1, const vector<Atom>& atomVec2 ){
	if( atomVec2.empty() ){
		cerr<<" no atom in vector to compute distance"<<endl;
		throw;
	}

	float			maxDis = 10e10;
	float			dis = 0;
	pair<Atom, Atom>			closestAtom;

	for( size_t j=0; j<atomVec2.size(); j++ ){
		dis = atomDis( atom1, atomVec2[j] );
		if( dis < maxDis ){
			maxDis = dis;
			closestAtom = make_pair( atom1, atomVec2[j] );
		}
	}
	return closestAtom.second;
}

pair<Atom, Atom>
findClosestAtom( const vector<Atom>& atomVec1, const vector<Atom>& atomVec2 ){
	if( atomVec1.empty() || atomVec2.empty() ){
		cerr<<" no atom in vector to compute distance"<<endl;
		cout<<atomVec1.size()<<" "<<atomVec2.size()<<endl;
		throw;
	}

	float			maxDis = 10e10;
	float			dis = 0;
	pair<Atom, Atom>			closestAtom;
	for( size_t i=0; i<atomVec1.size(); i++ ){
		for( size_t j=0; j<atomVec2.size(); j++ ){
			dis = atomDis( atomVec1[i], atomVec2[j] );
			if( dis < maxDis ){
				maxDis = dis;
				closestAtom = make_pair( atomVec1[i], atomVec2[j] );
			}
		}
	}
//	cout<<"!!! "<<closestAtom.first.get_index()<<" : "<<closestAtom.second.get_index()<<endl;
	return closestAtom;
}

float
atomsMaxRange( const vector<Atom>& atomVec ){
	float max = -10e10;
	float min = 10e10;

	for( size_t i=0; i<atomVec.size(); i++ ){
		max = max > atomVec[i].get_coord().x ? max : atomVec[i].get_coord().x ;
		max = max > atomVec[i].get_coord().y ? max : atomVec[i].get_coord().y ;
		max = max > atomVec[i].get_coord().z ? max : atomVec[i].get_coord().z ;

		min = min < atomVec[i].get_coord().x ? min : atomVec[i].get_coord().x ;
		min = min < atomVec[i].get_coord().y ? min : atomVec[i].get_coord().y ;
		min = min < atomVec[i].get_coord().z ? min : atomVec[i].get_coord().z ;
	}

	return max-min;
}

bool
setAttachedH( Atom& at, const vector<Atom> atVec){
    vector<Atom> hVec;
    hVec.clear();
    for( size_t i=0; i<atVec.size(); i++ ){
        Atom a=atVec[i];
        if( a.get_type().substr(0,1) == "H" ){
            if( atomDis(at, a) < 1.7 ){
                hVec.push_back(a);
            }
        }
    }
    at.set_attachedH(hVec);
    if( hVec.empty() ){
        return false;
    }else{
        return true;
    }
}

bool
setAttachedC( Atom& at, const vector<Atom> atVec){
    vector<Atom> cVec;
    cVec.clear();
    for( size_t i=0; i<atVec.size(); i++ ){
        Atom a=atVec[i];
        if( a.get_type().substr(0,1) == "C" || a.get_type().substr(0,1) == "P" ){
            if( atomDis(at, a) < 1.7 ){
                cVec.push_back(a);
            }
        }
    }
    at.set_attachedC(cVec);
    if( cVec.empty() ){
        return false;
    }else{
        return true;
    }
}

float
atomsMinX(  const vector<Atom>& atomVec ){
	float minx = 10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		minx = minx<c.x ? minx:c.x;
	}
	return minx;
}
float
atomsMaxX(  const vector<Atom>& atomVec ){
	float maxx = -10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		maxx = maxx>c.x ? maxx:c.x;
	}
	return maxx;
}
float
atomsMinY(  const vector<Atom>& atomVec ){
	float miny = 10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		miny = miny<c.y ? miny:c.y;
	}
	return miny;
}
float
atomsMaxY(  const vector<Atom>& atomVec ){
	float maxy = -10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		maxy = maxy>c.y ? maxy:c.y;
	}
	return maxy;
}
float
atomsMinZ(  const vector<Atom>& atomVec ){
	float minz = 10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		minz = minz<c.z ? minz:c.z;
	}
	return minz;
}
float
atomsMaxZ(  const vector<Atom>& atomVec ){
	float maxz = -10e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		Coord c = atomVec[i].get_coord();
		maxz = maxz>c.z ? maxz:c.z;
	}
	return maxz;
}
