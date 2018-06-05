/*
 * surface.cpp
 *
 *  Created on: Mar 12, 2014
 *      Author: stan
 */

#include<time.h>
#include"surface.h"

float saveTime = 0;

bool
findEdge( const vector<surfaceEdge>& edgeVec, const surfaceEdge& edge ){
	for( size_t i=0; i<edgeVec.size(); i++ ){
		if( ( edgeVec[i].get_firstAtom() == edge.get_firstAtom() && edgeVec[i].get_secAtom() == edge.get_secAtom() ) ||
				( edgeVec[i].get_firstAtom() == edge.get_secAtom() && edgeVec[i].get_firstAtom() == edge.get_secAtom() ) ){
			return true;
		}
		return false;
	}
}

//-------------------------------------------------------------------------------------------------------------------

void
reducedSurface::print(){
	at1.print();
	at2.print();
	at3.print();
	cout<<"probe:";
	for( size_t i=0; i<probeCoordVec.size(); i++ ){
		probeCoordVec[i].print();
	}
	cout<<endl;
}

//--------------------------------------------------------------------------------------------------------------------

sphereSurface::sphereSurface( const Atom& a1, const Atom& a2, const Atom& a3, const Coord& p,
		const float& pR ): at1(a1),at2(a2),at3(a3),probe(p),probeRadius(pR){

	crash = false;

	Coord c1 = at1.get_coord();
	Coord c2 = at2.get_coord();
	Coord c3 = at3.get_coord();

	arch_12 = getArchPoints( c1, c2, probe, probeRadius, 10 );
	arch_13 = getArchPoints( c1, c3, probe, probeRadius, 10 );
	arch_23 = getArchPoints( c2, c3, probe, probeRadius, 10 );

}


void
sphereSurface::print(){
	cout<<"a1:";
	at1.print();
	cout<<"a2:";
	at2.print();
	cout<<"a3:";
	at3.print();
	cout<<"probe:";
	probe.print();
	cout<<endl;
}

//---------------------------------------------------------------------------------------------------------------------

bool
reducedSurface_connect::computeReducedSurfaceArch( const vector<Atom>& neighbors ){
	Atom a1 = at1;
	Atom a2 = at2;
	Atom a3 = at3;
	Atom a4 = at4;

	Coord c3 = a3.get_coord();
	Coord c4 = a4.get_coord();

	Coord p1 = probeCoordA;
	Coord p2 = probeCoordB;
	Coord center =  getTriangleVerticalPoint(a3.get_coord(), a4.get_coord(), p1);
	Coord axis = cross( p1-center, p2-center );

	double dis_p1center = getCoordDis(p1, center);
	double dis_p1a1 = getCoordDis(p1, a1.get_coord());
	double dis_p2a2 = getCoordDis(p2, a2.get_coord());
	double dis_p1a3 = getCoordDis(p1, a3.get_coord());
	double dis_p1a4 = getCoordDis(p1, a4.get_coord());

	double probeRadius =  dis_p1a1 - a1.get_vdwRadius();
	double angle_p1p2 = getAngle(center, p1, p2);

	double angleRange = angle_p1p2;

//	int numOut = (angleRange*180/3.14159)/4;

	vector<Coord> coVec;

	double angleStep = (PI/180.0)*5.0;

	bool flag = true;

	coVec.push_back(p1);
	for( double angle=angleStep; angle<angleRange; angle=angle+angleStep ){
		Coord pp1 = p1 - center;
		Coord pp2 = p2 - center;
		Coord c = rotateAroundAxis( pp1, axis, angle );
		Coord coo = c+center;
		coVec.push_back( coo );

		for( size_t i=0; i<neighbors.size(); i++ ){
			double dis = getCoordDis( neighbors[i].get_coord(), coo );
			if( dis < ( probeRadius + neighbors[i].get_vdwRadius() - 0.0001 ) ){
//				cout<<"dis:"<<dis<<" +:"<<( probeRadius + neighbors[i].get_vdwRadius())<<endl;
//				neighbors[i].print();
				flag = false;
				break;
			}
		}
		if( !flag ){
			break;
		}
	}
	coVec.push_back(p2);

	if( flag ){
		coordVec = coVec;
		return true;
	}

	//--------------------------------------------------------
	flag = true;
	coVec.clear();
	angleRange = 2*PI - angle_p1p2;

	for( double ag = angleStep; ag < angleRange; ag+=angleStep ){

		double angle = 2*PI-ag;

		Coord pp1 = p1 - center;
		Coord pp2 = p2 - center;
		Coord c = rotateAroundAxis( pp1, axis, angle );
		Coord coo = c+center;
		coVec.push_back( coo );

		for( size_t i=0; i<neighbors.size(); i++ ){
			double dis = getCoordDis( neighbors[i].get_coord(), coo );
			if( dis < ( probeRadius + neighbors[i].get_vdwRadius() - 0.0001 ) ){
//				cout<<"dis:"<<dis<<" +:"<<( probeRadius + neighbors[i].get_vdwRadius())<<endl;
//				neighbors[i].print();
				flag = false;
				break;
			}
		}
		if( !flag ){
			break;
		}
	}
	coVec.push_back(p2);
	if( flag ){
		coordVec = coVec;
		return true;
	}

	return false;
}

vector<Coord>
reducedSurface_connect::getCrashProbeCoordConcave(const Coord& c1, const Coord& c2){
	vector<Coord>   crashCoordVec;
	Coord  center = getTriangleVerticalPoint(probeCoordA, probeCoordB, c1);
	double angleRange = getAngle(center, c1, c2);
	Coord  axis = cross( c1-center, c2-center );
	double angleStep = PI/180.0;
	crashCoordVec.clear();
	crashCoordVec.push_back(c1);
	for( double angle = angleStep; angle<angleRange; angle+=angleStep ){
		Coord c = ( rotateAroundAxis(c1-center, axis, angle ) + center);
		crashCoordVec.push_back(c);
	}
	crashCoordVec.push_back(c2);
	return crashCoordVec;
}

vector<Coord>
reducedSurface_connect::getCrashProbeCoordConvex(const Coord& c1, const Coord& c2){
	vector<Coord>   crashCoordVec;
	Coord  center = getTriangleVerticalPoint(probeCoordA, probeCoordB, c1);
	double angleRange = 2*PI - getAngle(center, c1, c2);

	Coord  axis = cross( c1-center, c2-center );
	double angleStep = PI/180.0;
	crashCoordVec.clear();
	crashCoordVec.push_back(c1);

	for( double ang = angleStep; ang<angleRange; ang+=angleStep ){
		double angle = 2*PI - ang;
		Coord c = ( rotateAroundAxis(c1-center, axis, angle ) + center);
		crashCoordVec.push_back(c);
	}
	crashCoordVec.push_back(c2);
	return crashCoordVec;
}

void
reducedSurface_connect::print()const{
	at1.print();
	at2.print();
	at3.print();
	at4.print();
	cout<<"probe1: ";
	probeCoordA.print();
	cout<<"probe2: ";
	probeCoordB.print();
	cout<<"probe:"<<coordVec.size()<<endl;
	cout<<endl;
}

//---------------------------------------------------------------------------------------------------------------------

vector<sphereSurface>
get_sphereSurface( const vector<reducedSurface>& rsVec,
		const vector<reducedSurface_connect>& rscVec ){
	vector<sphereSurface> sVec;
	for( size_t i=0; i<rsVec.size(); i++ ){
		vector<Coord>  pVec = rsVec[i].get_probeCoordVec();
		for( size_t j=0; j<pVec.size(); j++ ){
			Coord probe = pVec[j];
			Atom  a1 = rsVec[i].get_at1();
			Atom  a2 = rsVec[i].get_at2();
			Atom  a3 = rsVec[i].get_at3();

			float dis = getCoordDis( a1.get_coord(), probe );
			float probeRadius = dis - a1.get_vdwRadius();

			sphereSurface s( a1, a2, a3, probe, probeRadius );
			if( find(sVec.begin(), sVec.end(), s ) == sVec.end() ){
				sVec.push_back(s);
			}
		}
	}

	//--------------------------------------
	for( size_t i=0; i<sVec.size(); i++ ){
//	for( size_t i=134; i<135; i++ ){

		Coord c1 = sVec[i].get_at1().get_coord();
		Coord c2 = sVec[i].get_at2().get_coord();
		Coord c3 = sVec[i].get_at3().get_coord();
		Coord probe = sVec[i].get_probe();
		float probeRadius = sVec[i].get_probeRadius();
		vector<Coord> arch_12 = getArchPoints( c1, c2, probe, probeRadius, 10 );
		vector<Coord> arch_13 = getArchPoints( c1, c3, probe, probeRadius, 10 );
		vector<Coord> arch_23 = getArchPoints( c2, c3, probe, probeRadius, 10 );

		bool flag = false;

		double disToSame = 0.001;

		if(0){
			for( size_t j=0; j<rscVec.size(); j++ ){

				if( sVec[i].get_probe() == rscVec[j].get_probeCoordA() && rscVec[j].probeCrash ){
					vector<Coord> updateCoordVec = rscVec[j].get_crashCoordA();
					if( sVec[i].get_at1() != rscVec[j].get_at3() && sVec[i].get_at1() != rscVec[j].get_at4() ){
						cout<<"-1-"<<endl;
						if( getCoordDis(arch_23.front(), updateCoordVec.front()) <
							getCoordDis(arch_23.front(), updateCoordVec.back()) ){
							arch_23 = updateCoordVec;
						}else{
							arch_23 = inverseCoordVec( updateCoordVec );
						}
						flag = true;
					}else if( sVec[i].get_at2() != rscVec[j].get_at3() && sVec[i].get_at2() != rscVec[j].get_at4() ){
						cout<<"-2-"<<endl;
						if( getCoordDis(arch_13.front(), updateCoordVec.front())<
							getCoordDis(arch_13.front(), updateCoordVec.back()) ){
							arch_13 = updateCoordVec;
						}else{
							arch_13 = inverseCoordVec( updateCoordVec );
						}
						flag = true;
					}else if( sVec[i].get_at3() != rscVec[j].get_at3() && sVec[i].get_at3() != rscVec[j].get_at4() ){
						cout<<"-3-"<<endl;
						if( getCoordDis(arch_12.front(), updateCoordVec.front())<
							getCoordDis(arch_12.front(), updateCoordVec.back()) ){
							arch_12 = updateCoordVec;
						}else{
							arch_12 = inverseCoordVec( updateCoordVec );
						}
						flag = true;
					}
				}else if( sVec[i].get_probe() == rscVec[j].get_probeCoordB() && rscVec[j].probeCrash ){

					vector<Coord> updateCoordVec = rscVec[j].get_crashCoordB();

					if( sVec[i].get_at1() != rscVec[j].get_at3() && sVec[i].get_at1() != rscVec[j].get_at4() ){
						cout<<"-4-"<<endl;
						if( getCoordDis(arch_23.front(), updateCoordVec.front())<
							getCoordDis(arch_23.front(), updateCoordVec.back()) ){
							arch_23 = updateCoordVec;
						}else{
							arch_23 = inverseCoordVec( updateCoordVec );
						}
						flag = true;
					}else if( sVec[i].get_at2() != rscVec[j].get_at3() && sVec[i].get_at2() != rscVec[j].get_at4() ){
						cout<<"-5-"<<endl;
						if( (getCoordDis( arch_13.front(), updateCoordVec.front() ))<
							(getCoordDis( arch_13.front(), updateCoordVec.back() )) ){
							arch_13 = updateCoordVec;
						}else{
							arch_13 = inverseCoordVec( updateCoordVec );
						}
						if(0){
							cout<<"after:"<<endl;
							for( size_t k=0; k<arch_13.size(); k++ ){
								arch_13[k].print();
							}
						}
						flag = true;
					}else if( sVec[i].get_at3() != rscVec[j].get_at3() && sVec[i].get_at3() != rscVec[j].get_at4() ){
						cout<<"-6-"<<endl;
						if(1){
							cout<<"before:"<<endl;
							for( size_t k=0; k<arch_12.size(); k++ ){
								arch_12[k].print();
							}
						}
						if( getCoordDis(arch_12.front(), updateCoordVec.front())<
							getCoordDis(arch_12.front(), updateCoordVec.back()) ){
							arch_12 = updateCoordVec;
						}else{
							arch_12 = inverseCoordVec( updateCoordVec );
						}
						if(1){
							cout<<"after:"<<endl;
							for( size_t k=0; k<arch_12.size(); k++ ){
								arch_12[k].print();
							}
						}
						flag = true;
					}
				}
				if( flag ){
					cout<<"j:"<<j<<endl;
					rscVec[j].print();
					break;
				}
			}

		}
		if( flag ){
			sVec[i].crash = true;
		}

//		if(flag){
//			sVec[i].print();
//		}
//		cout<<"arch_12:"<<arch_12.size()<<endl;
//		cout<<"arch_13:"<<arch_13.size()<<endl;
//		cout<<"arch_23:"<<arch_23.size()<<endl;

		sVec[i].set_arch12( arch_12 );
		sVec[i].set_arch13( arch_13 );
		sVec[i].set_arch23( arch_23 );
	}

	return sVec;
}

double
getAngle( const Coord& init, const Coord& p1, const Coord& p2  ){
	if( p1 == p2 ){
		return 0;
	}

	double a = getCoordDis( p1, init)*getCoordDis(p2, init );
	double b = (p1-init)*(p2-init);

	double ba = b/a;
	if( abs( ba - 1 ) < 0.0001 ){
		ba = 1;
	}
	if( abs( ba + 1 ) < 0.0001 ){
		ba = -1;
	}
	double angle = acos(ba);

	return angle;
}

Coord
rotateAroundAxis(const Coord& input, const Coord& axis, const double& angle){
	double dis = sqrt(axis.x*axis.x+axis.y*axis.y+axis.z*axis.z);
	if( dis == 0 ){
		cerr<<"!!!no axis for rotating around"<<endl;
		return input;
	}
	Coord u ((double)axis.x/dis, (double)axis.y/dis, (double)axis.z/dis);

	double a11 = cos(angle)+u.x*u.x*(1-cos(angle));
	double a12 = u.x*u.y*(1-cos(angle))-u.z*sin(angle);
	double a13 = u.x*u.z*(1-cos(angle))+u.y*sin(angle);
	double a21 = u.y*u.x*(1-cos(angle))+u.z*sin(angle);
	double a22 = cos(angle)+u.y*u.y*(1-cos(angle));
	double a23 = u.y*u.z*(1-cos(angle))-u.x*sin(angle);
	double a31 = u.z*u.x*(1-cos(angle))-u.y*sin(angle);
	double a32 = u.z*u.y*(1-cos(angle))+u.x*sin(angle);
	double a33 = cos(angle)+u.z*u.z*(1-cos(angle));

	double x = a11*input.x+a12*input.y+a13*input.z;
	double y = a21*input.x+a22*input.y+a23*input.z;
	double z = a31*input.x+a32*input.y+a33*input.z;

	Coord c(x, y, z);
	return c;
}

/**
 * vertical point on edge
 */
Coord
getTriangleVerticalPoint( const Coord& a, const Coord& b, const Coord& p ){
	double disab = getCoordDis( a, b );
	double disap = getCoordDis( a, p );
	double disbp = getCoordDis( b, p );

	if(disab == 0 || disap ==0){
		throw;
	}
	double cosa = ((b-a)*(p-a))/(disab*disap);
	double disav = cosa*disap;
	double frac = disav/disab;

	double x = a.x + frac*(b.x-a.x);
	double y = a.y + frac*(b.y-a.y);
	double z = a.z + frac*(b.z-a.z);
	Coord c(x,y,z);
	return c;
}

reducedSurface
getFirstReducedSurface( const vector<Atom>& atomVec, const double probeRadius ){
	vector<reducedSurface> rsVec;
	// find the atom whose x-coordinate minus its radius is minimal
	Atom						a1;								//leftMostAtom;
	size_t						index1 = 0;
	double						xMinimum = 1e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( xMinimum > ( atomVec[i].get_coord().x - atomVec[i].get_vdwRadius() ) ){
			xMinimum =  atomVec[i].get_coord().x - atomVec[i].get_vdwRadius() ;
			index1 = i;
		}
	}
	a1 = atomVec[ index1 ] ;

	Atom						a2;
	size_t						index2 = 0;
	double						miniDis = 1e10;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( i != index1 ){
			double				dis = atomDis( atomVec[i], a1 );
			if( dis < ( a1.get_vdwRadius() + atomVec[i].get_vdwRadius() + 2 * probeRadius ) ){
				if( miniDis > dis ){
					miniDis = dis ;
					index2 = i;
				}
			}
		}
	}
	a2 = atomVec[ index2 ];

	Atom						a3 = atomVec[0] ;
	Coord						probe;
	vector<Coord>				probeVec;
	for( size_t i=0; i<atomVec.size(); i++ ){
		if( getFixedProbe( a1, a2, atomVec[i], atomVec, probeRadius, probeVec ) ){
			a3 = atomVec[i];
			break;
		}
	}
//	probe.print();
	reducedSurface		rs( a1, a2, a3, probeVec );
	return			rs;
}

bool
atomExistInRS( const Atom& at, const vector<reducedSurface>& rsVec ){
	for( size_t i=0; i<rsVec.size(); i++ ){
		size_t		index = at.get_index();
		if(		index == rsVec[i].get_at1().get_index() ||
				index == rsVec[i].get_at2().get_index() ||
				index == rsVec[i].get_at3().get_index()   ){
			return true;
		}
	}
	return false;
}

/**
 * probe have collision with atomVec, true indicates collision
 * @param probeCoord
 * @param probeRadius
 * @param atomVec
 * @return
 */
bool
probeCollision( const Coord& probeCoord, const double probeRadius, const vector<Atom>& atomVec ){
	for( size_t i=0; i<atomVec.size(); i++ ){
		double		dis = getCoordDis( probeCoord, atomVec[i].get_coord() );
		if( dis < (probeRadius + atomVec[i].get_vdwRadius() ) ){
			return		true;
		}
	}
	return false;
}

bool
probeMayExist( const Atom at1, const Atom at2, const Atom at3,  const double probeRadius ){
	double			dis12 = atomDis( at1, at2 );
	if( dis12 > ( at1.get_vdwRadius() + at2.get_vdwRadius() + probeRadius ) ){
		return		false;
	}
	double			dis13 = atomDis( at1, at3 );
	if( dis13 > ( at1.get_vdwRadius() + at3.get_vdwRadius() + probeRadius ) ){
		return		false;
	}
	double			dis23 = atomDis( at2, at3 );
	if( dis23 > ( at2.get_vdwRadius() + at3.get_vdwRadius() + probeRadius ) ){
		return		false;
	}
	return			true ;
}

/**
 * get the coordinates of fixed probe from three atoms, which are all tangential to the probe.
 * Given three atoms, there will be two fixed probe around the three atoms.
 * @param at1
 * @param at2
 * @param at3
 * @param probeRadius
 * @return
 */
pair<Coord, Coord>
getFixedProbe( const Atom at1, const Atom at2, const Atom at3, const double probeRadius ){
	//
	double		dis12 = atomDis( at1, at2 );
	double		dis13 = atomDis( at1, at3 );
	double		dis23 = atomDis( at2, at3 );
	if( 	dis12 > ( at1.get_vdwRadius() + at2.get_vdwRadius() + 2 * probeRadius ) ||
			dis13 > ( at1.get_vdwRadius() + at3.get_vdwRadius() + 2 * probeRadius ) ||
			dis23 > ( at2.get_vdwRadius() + at3.get_vdwRadius() + 2 * probeRadius ) 	){
		throw 	(int)0 ;
	}

	double		x = 0 ;			// coordinates of the probe
	double		y = 0 ;
	double		z = 0 ;

	double		r1 = ( at1.get_vdwRadius() + probeRadius) * ( at1.get_vdwRadius() + probeRadius ) ;
	double		r2 = ( at2.get_vdwRadius() + probeRadius) * ( at2.get_vdwRadius() + probeRadius ) ;
	double		r3 = ( at3.get_vdwRadius() + probeRadius) * ( at3.get_vdwRadius() + probeRadius ) ;

	double		x1 = at1.get_coord().x ;
	double		y1 = at1.get_coord().y ;
	double		z1 = at1.get_coord().z ;

	double		x2 = at2.get_coord().x ;
	double		y2 = at2.get_coord().y ;
	double		z2 = at2.get_coord().z ;

	double		x3 = at3.get_coord().x ;
	double		y3 = at3.get_coord().y ;
	double		z3 = at3.get_coord().z ;

	if( ( (2*y2-2*y1)*(2*x3-2*x2) - (2*y3-2*y2)*(2*x2-2*x1) ) == 0 ){
		throw		(int)1;
	}

	//				y = a * z + b ;
	double		a =  - ( (2*z2 - 2*z1)*( 2*x3-2*x2 ) - (2*z3 - 2*z2)*(2*x2 - 2*x1) )/( (2*y2-2*y1)*(2*x3-2*x2) - (2*y3-2*y2)*(2*x2-2*x1) );
	double		b =  ( (r1-r2 + x2*x2-x1*x1 + y2*y2-y1*y1 + z2*z2-z1*z1)*(2*x3-2*x2) - ( r2-r3 + x3*x3-x2*x2 + y3*y3-y2*y2+z3*z3-z2*z2 )*(2*x2-2*x1) ) /
							( (2*y2-2*y1)*(2*x3-2*x2) - (2*y3-2*y2)*(2*x2-2*x1) );

	if( ( (2*x2-2*x1)*(2*y3-2*y2) - (2*x3-2*x2)*(2*y2-2*y1) ) == 0 ){
		throw		(int)1;
	}

	//				x = m * z + n ;
	double		m =- ( (2*z2-2*z1)*(2*y3-2*y2) - (2*z3-2*z2)*(2*y2-2*y1) ) / ( (2*x2-2*x1)*(2*y3-2*y2) - (2*x3-2*x2)*(2*y2-2*y1) ) ;
	double		n = ( (r1-r2 + x2*x2-x1*x1 + y2*y2-y1*y1 + z2*z2 - z1*z1 )*(2*y3-2*y2) - ( r2-r3 + x3*x3 -x2*x2 + y3*y3 -y2*y2 + z3*z3 - z2*z2 )*(2*y2-2*y1) ) /
							( (2*x2-2*x1)*(2*y3-2*y2) - (2*x3-2*x2)*(2*y2-2*y1) );

	// (x - x3)^2 + (y-y3)^2 + (z-z3)^2 = r1
	// (mz+n-x3)^2 + (az+b-y3)^2 +(z-z3)^2 = r1
	double		sq =  (2*m*(n-x3) + 2*a*(b-y3) - 2*z3 )*(2*m*(n-x3) + 2*a*(b-y3) - 2*z3 ) - 4*(m*m + a*a + 1)*( (n-x3)*(n-x3) + (b-y3)*(b-y3) + z3*z3 - r3 )  ;

	if( sq < 0 ){
		throw		(int)1;
	}else{
		z = ( -( 2*m*(n-x3) + 2*a*(b-y3) - 2*z3 )   +  sqrt(sq) ) /
				( 2 * ( m* m + a*a + 1 ) ) ;
		y = a * z + b;
		x = m * z + n ;

		Coord		co1( (double)x, (double)y, (double)z );

		z = ( -( 2*m*(n-x3) + 2*a*(b-y3) - 2*z3 )   -  sqrt( sq ) ) /
				( 2 * ( m* m + a*a + 1 ) ) ;
		y = a * z + b;
		x = m * z + n ;
		Coord		co2( (double)x, (double)y, (double)z );
		pair<Coord, Coord>		coo = make_pair( co1, co2 );

		return					coo ;
	}
}

pair<Coord, Coord>
getFixedProbe( const Coord co1, const Coord co2, const Coord co3, const double probeRadius ){
	Atom		a1;
	Atom		a2;
	Atom		a3;
	a1.set_coord( co1.x, co1.y, co1.z );
	a1.set_vdwRadius( (double)1.2 );
	a2.set_coord( co2.x, co2.y, co2.z );
	a2.set_vdwRadius( (double)1.2 );
	a3.set_coord( co3.x, co3.y, co3.z );
	a3.set_vdwRadius( (double)1.2 );
	pair<Coord, Coord>		coordPair;
	coordPair = getFixedProbe( a1, a2, a3, probeRadius );
	return	coordPair;
}

/**
 * get the probe, and make sure it does not collide with any atoms in the atomVec
 * @param at1
 * @param at2
 * @param at3
 * @param atomVec
 * @param probeRadius
 * @return
 */
bool
getFixedProbe( const Atom at1, const Atom at2, const Atom at3,
		const vector<Atom>& atomVec, const double probeRadius, vector<Coord>& probeVec ){

	if(1){
        double dis12 = getCoordDis( at1.get_coord(), at2.get_coord() );
		double dis13 = getCoordDis( at1.get_coord(), at3.get_coord() );
		double dis23 = getCoordDis( at2.get_coord(), at3.get_coord() );
		if( ( at1.get_vdwRadius()+at2.get_vdwRadius()+2*probeRadius < dis12 ) ||
			( at1.get_vdwRadius()+at3.get_vdwRadius()+2*probeRadius < dis13 ) ||
			( at2.get_vdwRadius()+at3.get_vdwRadius()+2*probeRadius < dis23 )){
			return 0;
		}
	}

	pair<Coord, Coord>	probePair ;
	try{
		probePair	= getFixedProbe( at1, at2, at3, probeRadius );
	}
	catch( int e ){
		if( e == 0 ) return false;
		if( e == 1 ) return	false;
		return	false;
	}

	bool 							first = true;				// indicate if the first of probePair is valid
	bool 							second = true;

	for( size_t i=0; i<atomVec.size(); i++ ){
		if( getCoordDis( probePair.first, atomVec[i].get_coord() ) + 0.0001 < ( probeRadius + atomVec[i].get_vdwRadius() ) ){
			first = false;
		}
		if( getCoordDis( probePair.second, atomVec[i].get_coord() )+ 0.0001 < ( probeRadius + atomVec[i].get_vdwRadius() ) ){
			second = false;
		}
	}
	probeVec.clear();
	if( first ){
		probeVec.push_back(probePair.first);
	}
	if( second ){
		probeVec.push_back(probePair.second);
	}
	if( first || second ){
		return		true;
	}else{
		return		false;
	}
}

bool
getFixedProbe( const Atom at1, const Atom at2, const Atom at3,
		const vector<Atom>& atomVec, const double probeRadius, Coord& probe  ){

	if(1){
		double dis12 = getCoordDis( at1.get_coord(), at2.get_coord() );
		double dis13 = getCoordDis( at1.get_coord(), at3.get_coord() );
		double dis23 = getCoordDis( at2.get_coord(), at3.get_coord() );
		if( ( at1.get_vdwRadius()+at2.get_vdwRadius()+2*probeRadius < dis12 ) ||
			( at1.get_vdwRadius()+at3.get_vdwRadius()+2*probeRadius < dis13 ) ||
			( at2.get_vdwRadius()+at3.get_vdwRadius()+2*probeRadius < dis23 )){
			return 0;
		}
	}

	pair<Coord, Coord>	probePair ;
	try{
		probePair	= getFixedProbe( at1, at2, at3, probeRadius );
	}
	catch( int e ){
		if( e == 0 ) return false;
		if( e == 1 ) return	false;
		return	false;
	}

	bool 							first = true;				// indicate if the first of probePair is valid
	bool 							second = true;

	for( size_t i=0; i<atomVec.size(); i++ ){
		if( getCoordDis( probePair.first, atomVec[i].get_coord() ) + 0.0001 < ( probeRadius + atomVec[i].get_vdwRadius() ) ){
			first = false;
		}
		if( getCoordDis( probePair.second, atomVec[i].get_coord() )+ 0.0001 < ( probeRadius + atomVec[i].get_vdwRadius() ) ){
			second = false;
		}
	}
	if( first ){
		probe = probePair.first ;
		return		true;
	}else if( second ){
		probe = probePair.second ;
		return		true;
	}else{

		return		false;
	}
	return			false;
}

void
getReducedSurface(const vector<Atom>& atomVec,
		const double probeRadius,
		vector<reducedSurface>& rsVec,
		vector<sphereSurface>& ssVec,
		vector<reducedSurface_connect>& rscVec){

	vector<reducedSurface> saveRsVec;
	rsVec.push_back( getFirstReducedSurface( atomVec, probeRadius));
	saveRsVec.push_back( getFirstReducedSurface( atomVec, probeRadius));

	vector< vector<Atom> > neighborAtomVec;
    double neighborDis = 6.0;
	for( size_t i=0; i<atomVec.size(); i++ ){
		vector<Atom> aVec;
		aVec.push_back( atomVec[i] );
		for( size_t j=0; j<atomVec.size(); j++ ){
			double dis = getCoordDis( atomVec[i].get_coord(), atomVec[j].get_coord() );
            if( dis < neighborDis ){
				aVec.push_back( atomVec[j] );
			}
		}
		neighborAtomVec.push_back( aVec );
	}
	vector<surfaceEdge> edgeVec;

	time_t	t1;
	t1 = clock();

    float saveT = 0;

	while( !saveRsVec.empty() ){
//        cout<<"\r ->:"<<saveRsVec.size();
		reducedSurface rs = saveRsVec.back();
		saveRsVec.pop_back();

		Atom a1 = rs.get_at1();
		Atom a2 = rs.get_at2();
		Atom a3 = rs.get_at3();

		vector<Coord> probePreVec = rs.get_probeCoordVec();

		vector<Atom>    neighborAtoms;
		vector<surfaceEdge> tempEdge;

		if( probePreVec.size() == 2 ){
			neighborAtoms.clear();
			for( size_t j=0; j<neighborAtomVec.size(); j++ ){
				if( a1 == neighborAtomVec[j][0] || a2== neighborAtomVec[j][0]||
						a3== neighborAtomVec[j][0] ){
					neighborAtoms = neighborAtomVec[j];
					break;
				}
			}
			reducedSurface_connect rsc3( a3, a3, a1, a2, probePreVec[0], probePreVec[1] );
			if( rsc3.computeReducedSurfaceArch(neighborAtoms) ){
				if( find(rscVec.begin(), rscVec.end(), rsc3) == rscVec.end() ){
					rscVec.push_back(rsc3);
				}
			}
			reducedSurface_connect rsc2( a2, a2, a1, a3, probePreVec[0], probePreVec[1] );
			if( rsc2.computeReducedSurfaceArch(neighborAtoms) ){
				if( find(rscVec.begin(), rscVec.end(), rsc2) == rscVec.end() ){
					rscVec.push_back(rsc2);
				}
			}
			reducedSurface_connect rsc1( a1, a1, a2, a3, probePreVec[0], probePreVec[1] );
			if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
				if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
					rscVec.push_back(rsc1);
				}
			}
		}

		//edge a1-a2;
		for( size_t n =0; n < probePreVec.size(); n++ ){
			surfaceEdge e1( a1, a2, probePreVec[n] );
			tempEdge.clear();
			neighborAtoms.clear();
			for( size_t j=0; j<neighborAtomVec.size(); j++ ){
                if( a1 == neighborAtomVec[j][0] ||
                    a2 == neighborAtomVec[j][0] ){
					neighborAtoms = neighborAtomVec[j];
					break;
				}
			}
			for( size_t j=0; j<neighborAtoms.size(); j++ ){
				vector<Coord> probeVec;
				if( neighborAtoms[j] != a1 && neighborAtoms[j] != a2 && neighborAtoms[j] != a3 &&
						getFixedProbe( a1, a2, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
					reducedSurface surf( a1, a2, neighborAtoms[j], probeVec );

					if( find( edgeVec.begin(), edgeVec.end(), e1 ) == edgeVec.end() ){
						tempEdge.push_back(e1);
						rsVec.push_back( surf );
						saveRsVec.push_back( surf );

						double minDis = 10e10;
						size_t index, indexP;
						index = indexP = 0;
						for( size_t m=0; m<probeVec.size(); m++ ){
							for( size_t i=0; i<probePreVec.size(); i++ ){
								double dis= getCoordDis(probeVec[m], probePreVec[i]);
								if( minDis > dis ){
									minDis = dis;
									index = i;
									indexP = m;
								}
							}
						}

						reducedSurface_connect rsc( a3, neighborAtoms[j], a1, a2,
								probePreVec[index], probeVec[indexP] );

						if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
							if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
								rscVec.push_back(rsc);
							}
						}else{

						}

						if( probeVec.size()==2 && probePreVec.size()==2 ){
							if( index == 0 ){
								index = 1;
							}else if( index == 1 ){
								index = 0;
							}
							if( indexP == 0 ){
								indexP =1;
							}else if( indexP == 1 ){
								indexP = 0;
							}
							reducedSurface_connect rsc1( a3, neighborAtoms[j], a1, a2,
									probePreVec[index], probeVec[indexP] );
							if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
								if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
									rscVec.push_back(rsc1);
								}
							}else{
							}
						}
					}
				}
			}
			edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
		}

		//edge a1-a3;
		for( size_t n =0; n < probePreVec.size(); n++ ){
			surfaceEdge e2( a1, a3, probePreVec[n] );
			tempEdge.clear();
			neighborAtoms.clear();
			for( size_t j=0; j<neighborAtomVec.size(); j++ ){
				if( a1 == neighborAtomVec[j][0] || a3 == neighborAtomVec[j][0] ){
					neighborAtoms = neighborAtomVec[j];
					break;
				}
			}
			for( size_t j=0; j<neighborAtoms.size(); j++ ){
				vector<Coord> probeVec;
				if( neighborAtoms[j] != a1 && neighborAtoms[j] != a3 && neighborAtoms[j] != a2 &&
						getFixedProbe( a1, a3, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
					reducedSurface surf( a1, a3, neighborAtoms[j], probeVec );

					if( find( edgeVec.begin(), edgeVec.end(), e2 ) == edgeVec.end() ){
						tempEdge.push_back(e2);

						rsVec.push_back( surf );
						saveRsVec.push_back( surf );

						double minDis = 10e10;
						size_t index, indexP;
						index = indexP = 0;
						for( size_t m=0; m<probeVec.size(); m++ ){
							for( size_t i=0; i<probePreVec.size(); i++ ){
								double dis= getCoordDis(probeVec[m], probePreVec[i]);
								if( minDis > dis ){
									minDis = dis;
									index = i;
									indexP = m;
								}
							}
						}
						reducedSurface_connect rsc( a2, neighborAtoms[j], a1, a3,
								probePreVec[index], probeVec[indexP] );
						if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
							if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
								rscVec.push_back(rsc);
							}
//							rsc.print();
						}else{
//							cout<<"failed!!!"<<endl;
							//						rsc.print();
						}

						if( probeVec.size()==2 && probePreVec.size()==2 ){
							if( index == 0 ){
								index = 1;
							}else if( index == 1 ){
								index = 0;
							}
							if( indexP == 0 ){
								indexP =1;
							}else if( indexP == 1 ){
								indexP = 0;
							}
							reducedSurface_connect rsc1( a2, neighborAtoms[j], a1, a3,
									probePreVec[index], probeVec[indexP] );
							if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
								if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
									rscVec.push_back(rsc1);
								}
//								rsc1.print();
							}else{
//								cout<<"failed!!!"<<endl;
								//							rsc.print();
							}
						}
					}
				}
			}
			edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
		}

		//edge a2-a3;
		for( size_t n =0; n < probePreVec.size(); n++ ){
			surfaceEdge e3( a2, a3, probePreVec[n] );
			tempEdge.clear();
			neighborAtoms.clear();
			for( size_t j=0; j<neighborAtomVec.size(); j++ ){
				if( a2 == neighborAtomVec[j][0] || a3 == neighborAtomVec[j][0] ){
					neighborAtoms = neighborAtomVec[j];
					break;
				}
			}
			for( size_t j=0; j<neighborAtoms.size(); j++ ){
				vector<Coord> probeVec;
				if( neighborAtoms[j] != a2 && neighborAtoms[j] != a3 && neighborAtoms[j] != a1 &&
						getFixedProbe( a2, a3, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
					reducedSurface surf( a2, a3, neighborAtoms[j], probeVec );
					//				if( find( rsVec.begin(), rsVec.end(), surf ) == rsVec.end() ){
					if( find( edgeVec.begin(), edgeVec.end(), e3 ) == edgeVec.end() ){
						tempEdge.push_back(e3);

						rsVec.push_back( surf );
						saveRsVec.push_back( surf );
//						rs.print();
//						surf.print();

						double minDis = 10e10;
						size_t index, indexP;
						index = indexP = 0;
						for( size_t m=0; m<probeVec.size(); m++ ){
							for( size_t i=0; i<probePreVec.size(); i++ ){
								double dis= getCoordDis(probeVec[m], probePreVec[i]);
								if( minDis > dis ){
									minDis = dis;
									index = i;
									indexP = m;
								}
							}
						}
						reducedSurface_connect rsc( a1, neighborAtoms[j], a2, a3,
								probePreVec[index], probeVec[indexP] );
						if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
							if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
								rscVec.push_back(rsc);
							}
//							rsc.print();
						}else{
//							cout<<"failed!!!"<<endl;
							//						rsc.print();
						}

						if( probeVec.size()==2 && probePreVec.size()==2 ){
							if( index == 0 ){
								index = 1;
							}else if( index == 1 ){
								index = 0;
							}
							if( indexP == 0 ){
								indexP =1;
							}else if( indexP == 1 ){
								indexP = 0;
							}
							reducedSurface_connect rsc1( a1, neighborAtoms[j], a2, a3,
									probePreVec[index], probeVec[indexP] );
							if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
								if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
									rscVec.push_back(rsc1);
								}
//								rsc1.print();
							}else{
//								cout<<"failed!!!"<<endl;
								//							rsc.print();
							}
						}
					}
				}
			}
			edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
		}
	}
	vector<reducedSurface> tempVec;
    // time consumming
    for( size_t i=0; i<rsVec.size(); i++ ){
        if( find( tempVec.begin(), tempVec.end(), rsVec[i] ) == tempVec.end() ){
            tempVec.push_back( rsVec[i] );
        }
    }
    rsVec = tempVec;

	// get sphere surface
	ssVec.clear();
	for( size_t i=0; i<rsVec.size(); i++ ){
		vector<Coord>  pVec = rsVec[i].get_probeCoordVec();
		for( size_t j=0; j<pVec.size(); j++ ){
			Coord probe = pVec[j];
			Atom  a1 = rsVec[i].get_at1();
			Atom  a2 = rsVec[i].get_at2();
			Atom  a3 = rsVec[i].get_at3();

			float dis = getCoordDis( a1.get_coord(), probe );
			float probeRadius = dis - a1.get_vdwRadius();

			sphereSurface s( a1, a2, a3, probe, probeRadius );
			if( find(ssVec.begin(), ssVec.end(), s ) == ssVec.end() ){
				ssVec.push_back(s);
			}
		}
	}

	//update sphere surface archs with reducedSurface_connect,
	// edge crash!!!
	updateSphereSurface( rscVec, ssVec );

	cout<<"running time:"<<double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
}

void
getReducedResSurface(const vector<Atom>& atomVec,
        const double probeRadius,
        vector<reducedSurface>& rsVec,
        vector<sphereSurface>& ssVec,
        vector<reducedSurface_connect>& rscVec){

    vector<reducedSurface> saveRsVec;
    rsVec.push_back( getFirstReducedSurface( atomVec, probeRadius));
    saveRsVec.push_back( getFirstReducedSurface( atomVec, probeRadius));

    vector< vector<Atom> > neighborAtomVec;
    double neighborDis = 7.0;
    for( size_t i=0; i<atomVec.size(); i++ ){
        vector<Atom> aVec;
        aVec.push_back( atomVec[i] );
        for( size_t j=0; j<atomVec.size(); j++ ){
            double dis = getCoordDis( atomVec[i].get_coord(), atomVec[j].get_coord() );
            int resIndex1 = (int)atomVec[i].get_residueIndex();
            int resIndex2 = (int)atomVec[j].get_residueIndex();

            //			if( dis < neighborDis ){
            if( dis < neighborDis && abs(resIndex1-resIndex2) < 2 ){
                aVec.push_back( atomVec[j] );
            }
        }
        neighborAtomVec.push_back( aVec );
    }
    vector<surfaceEdge> edgeVec;
    vector<reducedSurface>  twoProbeSurfVec;

    time_t	t1;
    t1 = clock();

    while( !saveRsVec.empty() ){
//		cout<<"\r ->:"<<saveRsVec.size();
        reducedSurface rs = saveRsVec.back();
        saveRsVec.pop_back();

        Atom a1 = rs.get_at1();
        Atom a2 = rs.get_at2();
        Atom a3 = rs.get_at3();

        vector<Coord> probePreVec = rs.get_probeCoordVec();

        vector<Atom>    neighborAtoms;
        vector<surfaceEdge> tempEdge;


        if( probePreVec.size() == 2 ){
            neighborAtoms.clear();
            for( size_t j=0; j<neighborAtomVec.size(); j++ ){
                if( a1 == neighborAtomVec[j][0] || a2== neighborAtomVec[j][0]||
                        a3== neighborAtomVec[j][0] ){
                    neighborAtoms = neighborAtomVec[j];
                    break;
                }
            }
            reducedSurface_connect rsc3( a3, a3, a1, a2, probePreVec[0], probePreVec[1] );
            if( rsc3.computeReducedSurfaceArch(neighborAtoms) ){
                if( find(rscVec.begin(), rscVec.end(), rsc3) == rscVec.end() ){
                    rscVec.push_back(rsc3);
                }
            }
            reducedSurface_connect rsc2( a2, a2, a1, a3, probePreVec[0], probePreVec[1] );
            if( rsc2.computeReducedSurfaceArch(neighborAtoms) ){
                if( find(rscVec.begin(), rscVec.end(), rsc2) == rscVec.end() ){
                    rscVec.push_back(rsc2);
                }
            }
            reducedSurface_connect rsc1( a1, a1, a2, a3, probePreVec[0], probePreVec[1] );
            if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
                if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
                    rscVec.push_back(rsc1);
                }
            }
        }

        //edge a1-a2;
        for( size_t n =0; n < probePreVec.size(); n++ ){
            surfaceEdge e1( a1, a2, probePreVec[n] );
            tempEdge.clear();
            neighborAtoms.clear();
            for( size_t j=0; j<neighborAtomVec.size(); j++ ){
                if( e1.get_firstAtom() == neighborAtomVec[j][0] ||
                        e1.get_secAtom() == neighborAtomVec[j][0] ){
                    neighborAtoms = neighborAtomVec[j];
                    break;
                }
            }
            for( size_t j=0; j<neighborAtoms.size(); j++ ){
                vector<Coord> probeVec;
                if( neighborAtoms[j] != a1 && neighborAtoms[j] != a2 && neighborAtoms[j] != a3 &&
                        getFixedProbe( a1, a2, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
                    reducedSurface surf( a1, a2, neighborAtoms[j], probeVec );

                    if( find( edgeVec.begin(), edgeVec.end(), e1 ) == edgeVec.end() ){
                        tempEdge.push_back(e1);

                        rsVec.push_back( surf );
                        saveRsVec.push_back( surf );

                        double minDis = 10e10;
                        size_t index, indexP;
                        index = indexP = 0;
                        for( size_t m=0; m<probeVec.size(); m++ ){
                            for( size_t i=0; i<probePreVec.size(); i++ ){
                                double dis= getCoordDis(probeVec[m], probePreVec[i]);
                                if( minDis > dis ){
                                    minDis = dis;
                                    index = i;
                                    indexP = m;
                                }
                            }
                        }

                        reducedSurface_connect rsc( a3, neighborAtoms[j], a1, a2,
                                probePreVec[index], probeVec[indexP] );

                        if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
                            if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
                                rscVec.push_back(rsc);
                            }
                        }else{

                        }

                        if( probeVec.size()==2 && probePreVec.size()==2 ){
                            if( index == 0 ){
                                index = 1;
                            }else if( index == 1 ){
                                index = 0;
                            }
                            if( indexP == 0 ){
                                indexP =1;
                            }else if( indexP == 1 ){
                                indexP = 0;
                            }
                            reducedSurface_connect rsc1( a3, neighborAtoms[j], a1, a2,
                                    probePreVec[index], probeVec[indexP] );
                            if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
                                if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
                                    rscVec.push_back(rsc1);
                                }
                            }else{
                            }
                        }

                    }
                }
            }
            edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
        }

        //edge a1-a3;
        for( size_t n =0; n < probePreVec.size(); n++ ){
            surfaceEdge e2( a1, a3, probePreVec[n] );
            tempEdge.clear();
            neighborAtoms.clear();
            for( size_t j=0; j<neighborAtomVec.size(); j++ ){
                if( a1 == neighborAtomVec[j][0] || a3 == neighborAtomVec[j][0] ){
                    neighborAtoms = neighborAtomVec[j];
                    break;
                }
            }
            for( size_t j=0; j<neighborAtoms.size(); j++ ){
                vector<Coord> probeVec;
                if( neighborAtoms[j] != a1 && neighborAtoms[j] != a3 && neighborAtoms[j] != a2 &&
                        getFixedProbe( a1, a3, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
                    reducedSurface surf( a1, a3, neighborAtoms[j], probeVec );

                    if( find( edgeVec.begin(), edgeVec.end(), e2 ) == edgeVec.end() ){
                        tempEdge.push_back(e2);

                        rsVec.push_back( surf );
                        saveRsVec.push_back( surf );

                        double minDis = 10e10;
                        size_t index, indexP;
                        index = indexP = 0;
                        for( size_t m=0; m<probeVec.size(); m++ ){
                            for( size_t i=0; i<probePreVec.size(); i++ ){
                                double dis= getCoordDis(probeVec[m], probePreVec[i]);
                                if( minDis > dis ){
                                    minDis = dis;
                                    index = i;
                                    indexP = m;
                                }
                            }
                        }
                        reducedSurface_connect rsc( a2, neighborAtoms[j], a1, a3,
                                probePreVec[index], probeVec[indexP] );
                        if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
                            if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
                                rscVec.push_back(rsc);
                            }
//							rsc.print();
                        }else{
//							cout<<"failed!!!"<<endl;
                            //						rsc.print();
                        }

                        if( probeVec.size()==2 && probePreVec.size()==2 ){
                            if( index == 0 ){
                                index = 1;
                            }else if( index == 1 ){
                                index = 0;
                            }
                            if( indexP == 0 ){
                                indexP =1;
                            }else if( indexP == 1 ){
                                indexP = 0;
                            }
                            reducedSurface_connect rsc1( a2, neighborAtoms[j], a1, a3,
                                    probePreVec[index], probeVec[indexP] );
                            if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
                                if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
                                    rscVec.push_back(rsc1);
                                }
//								rsc1.print();
                            }else{
//								cout<<"failed!!!"<<endl;
                                //							rsc.print();
                            }
                        }
                    }
                }
            }
            edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
        }

        //edge a2-a3;
        for( size_t n =0; n < probePreVec.size(); n++ ){
            surfaceEdge e3( a2, a3, probePreVec[n] );
            tempEdge.clear();
            neighborAtoms.clear();
            for( size_t j=0; j<neighborAtomVec.size(); j++ ){
                if( a2 == neighborAtomVec[j][0] || a3 == neighborAtomVec[j][0] ){
                    neighborAtoms = neighborAtomVec[j];
                    break;
                }
            }
            for( size_t j=0; j<neighborAtoms.size(); j++ ){
                vector<Coord> probeVec;
                if( neighborAtoms[j] != a2 && neighborAtoms[j] != a3 && neighborAtoms[j] != a1 &&
                        getFixedProbe( a2, a3, neighborAtoms[j], neighborAtoms, probeRadius, probeVec )	){
                    reducedSurface surf( a2, a3, neighborAtoms[j], probeVec );
                    //				if( find( rsVec.begin(), rsVec.end(), surf ) == rsVec.end() ){
                    if( find( edgeVec.begin(), edgeVec.end(), e3 ) == edgeVec.end() ){
                        tempEdge.push_back(e3);

                        rsVec.push_back( surf );
                        saveRsVec.push_back( surf );
//						rs.print();
//						surf.print();

                        double minDis = 10e10;
                        size_t index, indexP;
                        index = indexP = 0;
                        for( size_t m=0; m<probeVec.size(); m++ ){
                            for( size_t i=0; i<probePreVec.size(); i++ ){
                                double dis= getCoordDis(probeVec[m], probePreVec[i]);
                                if( minDis > dis ){
                                    minDis = dis;
                                    index = i;
                                    indexP = m;
                                }
                            }
                        }
                        reducedSurface_connect rsc( a1, neighborAtoms[j], a2, a3,
                                probePreVec[index], probeVec[indexP] );
                        if( rsc.computeReducedSurfaceArch(neighborAtoms) ){
                            if( find(rscVec.begin(), rscVec.end(), rsc) == rscVec.end() ){
                                rscVec.push_back(rsc);
                            }
//							rsc.print();
                        }else{
//							cout<<"failed!!!"<<endl;
                            //						rsc.print();
                        }

                        if( probeVec.size()==2 && probePreVec.size()==2 ){
                            if( index == 0 ){
                                index = 1;
                            }else if( index == 1 ){
                                index = 0;
                            }
                            if( indexP == 0 ){
                                indexP =1;
                            }else if( indexP == 1 ){
                                indexP = 0;
                            }
                            reducedSurface_connect rsc1( a1, neighborAtoms[j], a2, a3,
                                    probePreVec[index], probeVec[indexP] );
                            if( rsc1.computeReducedSurfaceArch(neighborAtoms) ){
                                if( find(rscVec.begin(), rscVec.end(), rsc1) == rscVec.end() ){
                                    rscVec.push_back(rsc1);
                                }
//								rsc1.print();
                            }else{
//								cout<<"failed!!!"<<endl;
                                //							rsc.print();
                            }
                        }
                    }
                }
            }
            edgeVec.insert(edgeVec.end(), tempEdge.begin(), tempEdge.end());
        }
    }

    vector<reducedSurface> tempVec;
    for( size_t i=0; i<rsVec.size(); i++ ){
        if( find( tempVec.begin(), tempVec.end(), rsVec[i] ) == tempVec.end() ){
            tempVec.push_back( rsVec[i] );
        }
    }
    rsVec = tempVec;

    // get sphere surface
    ssVec.clear();
    for( size_t i=0; i<rsVec.size(); i++ ){
        vector<Coord>  pVec = rsVec[i].get_probeCoordVec();
        for( size_t j=0; j<pVec.size(); j++ ){
            Coord probe = pVec[j];
            Atom  a1 = rsVec[i].get_at1();
            Atom  a2 = rsVec[i].get_at2();
            Atom  a3 = rsVec[i].get_at3();

            float dis = getCoordDis( a1.get_coord(), probe );
            float probeRadius = dis - a1.get_vdwRadius();

            sphereSurface s( a1, a2, a3, probe, probeRadius );
            if( find(ssVec.begin(), ssVec.end(), s ) == ssVec.end() ){
                ssVec.push_back(s);
            }
        }
    }

    //update sphere surface archs with reducedSurface_connect,
    // edge crash!!!
    updateSphereSurface( rscVec, ssVec );

    cout<<"running time:"<<double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;
}

bool
updateSphereSurface( vector<reducedSurface_connect>& rscVec,
		vector<sphereSurface>& ssVec ){

	for( size_t i=0; i<rscVec.size(); i++ ){
		Atom a1 = rscVec[i].get_at1();
		Atom a2 = rscVec[i].get_at2();
		Atom a3 = rscVec[i].get_at3();
		Atom a4 = rscVec[i].get_at4();

		Coord c1 = a1.get_coord();
		Coord c2 = a2.get_coord();
		Coord c3 = a3.get_coord();
		Coord c4 = a4.get_coord();

		Coord p1 = rscVec[i].get_probeCoordA();
		Coord p2 = rscVec[i].get_probeCoordB();

		double probeRadius = rscVec[i].get_probeRadius();

		vector<Coord> arch1 = getArchPoints( c3, c4, p1, probeRadius, 10 );
		vector<Coord> arch2 = getArchPoints( c3, c4, p2, probeRadius, 10 );

		vector<Coord> aVec, bVec;

		for( size_t j=1; j<arch1.size()-1; j++ ){
			double dis0 = getCoordDis( arch1[j-1], arch2[j-1] );
			double dis1 = getCoordDis( arch1[j], arch2[j] );
			double dis2 = getCoordDis( arch1[j+1], arch2[j+1] );
			double  val = (dis1-dis0)*(dis2-dis1);

			if( val < 0 && abs( dis1-dis0 )>0.01 && abs( dis2-dis1 )>0.01 ){
				rscVec[i].probeCrash = true;
				vector<Coord> cVec;
				vector<Coord> cVec1 = rscVec[i].getCrashProbeCoordConcave(arch1[j],
						                  arch1[ arch1.size()-j ] );
				vector<Coord> cVec2 = rscVec[i].getCrashProbeCoordConvex(arch1[j],
						                  arch1[ arch1.size()-j ] );

				vector<Coord> probeArch = rscVec[i].get_coordVec();
				Coord centerProbeArch = probeArch[ probeArch.size()/2 ];

				if(  getCoordDis( cVec1[cVec1.size()/2], centerProbeArch ) >
					 getCoordDis( cVec2[cVec2.size()/2], centerProbeArch ) ){
					cVec = cVec1;
				}else{
					cVec = cVec2;
				}

//				cVec = cVec1;     //!!!!!!!!!!!!!!!!!!

				//new arch
				for( size_t k=0; k<j; k++ ){
					aVec.push_back( arch1[k] );
					bVec.push_back( arch2[k] );
				}
				aVec.insert( aVec.end(), cVec.begin(), cVec.end() );
				bVec.insert( bVec.end(), cVec.begin(), cVec.end() );
				for( size_t k=arch1.size()-j; k<arch1.size(); k++ ){
					aVec.push_back( arch1[k] );
				}
				for( size_t k=arch2.size()-j; k<arch2.size(); k++ ){
					bVec.push_back( arch2[k] );
				}

				break;
			}
		}

		if( rscVec[i].probeCrash ){
			for( size_t j=0; j<ssVec.size(); j++ ){
				Atom sa1 = ssVec[j].get_at1();
				Atom sa2 = ssVec[j].get_at2();
				Atom sa3 = ssVec[j].get_at3();
				Coord probe = ssVec[j].get_probe();
				if( sa1 == a3 && sa2 == a4 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch12( aVec );
					}else{
						ssVec[j].set_arch12( bVec );
					}
				}else if( sa1 == a4 && sa2 == a3 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch12( inverseCoordVec(aVec) );
					}else{
						ssVec[j].set_arch12( inverseCoordVec(bVec) );
					}
				}else if( sa1 == a3 && sa3 == a4 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch13( aVec );
					}else{
						ssVec[j].set_arch13( bVec );
					}
				}else if( sa1 == a4 && sa3 == a3 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch13( inverseCoordVec(aVec) );
					}else{
						ssVec[j].set_arch13( inverseCoordVec(bVec) );
					}
				}else if( sa2 == a3 && sa3 == a4 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch23( aVec );
					}else{
						ssVec[j].set_arch23( bVec );
					}
				}else if( sa2 == a4 && sa3 == a3 ){
					if( getCoordDis(probe, p1) < getCoordDis(probe, p2) ){
						ssVec[j].set_arch23( inverseCoordVec(aVec) );
					}else{
						ssVec[j].set_arch23( inverseCoordVec(bVec) );
					}
				}

			}
		}
	}
}

vector<reducedSurface>
getReducedSurface( const vector<Atom>& atomVec, const double probeRadius ){
	vector<reducedSurface>	rsVec;
	rsVec.push_back( getFirstReducedSurface( atomVec, probeRadius ) );
	vector<size_t>				surfaceAtomIndexVec;
	surfaceAtomIndexVec.push_back( rsVec[0].get_at1().get_index() );
	surfaceAtomIndexVec.push_back( rsVec[0].get_at2().get_index() );
	surfaceAtomIndexVec.push_back( rsVec[0].get_at3().get_index() );

	// edges that could be extends as edge of new surface
	vector<surfaceEdge>			edgeVec;
	surfaceEdge					se1( rsVec[0].get_at1(), rsVec[0].get_at2() );
	surfaceEdge					se2( rsVec[0].get_at1(), rsVec[0].get_at3() );
	surfaceEdge					se3( rsVec[0].get_at3(), rsVec[0].get_at2() );
	edgeVec.push_back( se1 );
	edgeVec.push_back( se2 );
	edgeVec.push_back( se3 );

	if(1){
		cout<<"first:"<< rsVec[0].get_at1().get_index()<<" "<< rsVec[0].get_at2().get_index()
				<<" "<< rsVec[0].get_at3().get_index()<<endl;
	}

	vector< vector<Atom> > neighborAtomVec;
	if(1){
		double neighborDis = 7.0;
		for( size_t i=0; i<atomVec.size(); i++ ){
			vector<Atom> aVec;
			aVec.push_back( atomVec[i] );
			for( size_t j=0; j<atomVec.size(); j++ ){
				double dis = getCoordDis( atomVec[i].get_coord(), atomVec[j].get_coord() );
				if( dis < neighborDis ){
					aVec.push_back( atomVec[j] );
				}
			}
			neighborAtomVec.push_back( aVec );
		}
	}

	time_t	t1;
	t1 = clock();
	while( !edgeVec.empty() ){
		cout<<"edgeVec:"<<edgeVec.size()<<endl;
		if(0){
			cout<<"edgeVec.size():"<<edgeVec.size()<<endl;
			for( size_t i=0; i<edgeVec.size(); i++ ){
				edgeVec[i].print();
			}
			cout<<endl;
		}
		surfaceEdge se =	edgeVec.back();
		edgeVec.pop_back();

		if(1){
			vector<Atom>    neighborAtoms;
			for( size_t j=0; j<neighborAtomVec.size(); j++ ){
				if( se.get_firstAtom() == neighborAtomVec[j][0] ){
					neighborAtoms = neighborAtomVec[j];
					break;
				}
			}
//			cout<<"neighbor:"<<neighborAtoms.size()<<endl;
			for( size_t j=0; j<neighborAtoms.size(); j++ ){
				Coord						probe;
				if( neighborAtoms[j] != se.get_secAtom() && neighborAtoms[j] != se.get_firstAtom() ){
					vector<Coord> probeVec;
					if( getFixedProbe( se.get_firstAtom(), se.get_secAtom(),
							neighborAtoms[j], neighborAtoms, probeRadius, probeVec ) ){
						for( size_t k=0; k<probeVec.size(); k++ ){
							reducedSurface	surf( se.get_secAtom(), se.get_firstAtom(), neighborAtoms[j], probeVec[k] );
							if( find( rsVec.begin(), rsVec.end(), surf ) == rsVec.end() ){
								rsVec.push_back( surf );
								//						surf.print();
								surfaceEdge				eg1( se.get_secAtom(), neighborAtoms[j] );
								surfaceEdge				eg2( se.get_firstAtom(), neighborAtoms[j] );

								if( find( edgeVec.begin(), edgeVec.end(), eg1 ) == edgeVec.end() ){
									edgeVec.push_back( eg1 );
									//							neighborAtoms[j].print();
								}
								if( find( edgeVec.begin(), edgeVec.end(), eg2 ) == edgeVec.end() ){
									edgeVec.push_back( eg2 );
									//							neighborAtoms[j].print();
								}
								break;
							}
						}
					}
				}
			}
		}
	}

	cout<<"running time:"<<double( clock() - t1 )/ CLOCKS_PER_SEC<<endl;

	if(0){
		for( size_t i=0; i<rsVec.size(); i++ ){
			cout<<"i:"<<i<<endl;
			rsVec[i].print();
		}
	}
	return		rsVec ;
}

/**
 * get the arch between c1 and c2, p is the probe with radius probeRadius
 */
vector<Coord>
getArchPoints( const Coord& c1, const Coord& c2, const Coord& p,
		const double& probeRadius, const int& numOfOut ){

	Coord axis = cross( c1-p, c2-p );
	double dis_c1p = getCoordDis(c1, p);

	double x1 = (c1.x-p.x)*probeRadius/dis_c1p + p.x;
	double y1 = (c1.y-p.y)*probeRadius/dis_c1p + p.y;
	double z1 = (c1.z-p.z)*probeRadius/dis_c1p + p.z;

	Coord p1(x1, y1, z1);
	p1 = p1-p;

	double dis_c2p = getCoordDis(c2, p);

	double x2 = (c2.x-p.x)*probeRadius/dis_c2p + p.x;
	double y2 = (c2.y-p.y)*probeRadius/dis_c2p + p.y;
	double z2 = (c2.z-p.z)*probeRadius/dis_c2p + p.z;

	Coord p2(x2, y2, z2);
	p2 = p2-p;

	double angleRange = getAngle(p, c1, c2);
	double angleStep = (PI/180.0)*5.0;

	vector<Coord> cVec;
	cVec.push_back( p1+p );
	for( double angle = angleStep; angle < angleRange; angle+=angleStep ){
		Coord c = rotateAroundAxis(p1, axis, angle);
		c = c+p;
		cVec.push_back(c);
	}
	cVec.push_back(p2+p);
	if(0){
		cout<<"arch;"<<endl;
		cout<<"c1:";
		c1.print();
		cout<<"c2:";
		c2.print();
		cout<<"p:";
		p.print();
		for( size_t i=0; i<cVec.size(); i++ ){
			cout<<i<<" ";
			cVec[i].print();
		}
	}
	return cVec;
}

/*
 * probe rotate around the axis of a1-a2 from p1 to p2,
 * this is the arch of reduced surface
 */
vector<Coord>
getReducedSurfaceArch( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut ){
//	cout<<"#########"<<endl;
	vector<Coord> coordVec;
	double stepX = (p2.x-p1.x)/(double)(numOut-1);
	double stepY = (p2.y-p1.y)/(double)(numOut-1);
	double stepZ = (p2.z-p1.z)/(double)(numOut-1);

	if( stepX >= stepY && stepX >= stepZ ){
		double d1 = getCoordDis( a1.get_coord(), p1 );
		double d2 = getCoordDis( a2.get_coord(), p1 );
		Coord c1 = a1.get_coord();
		Coord c2 = a2.get_coord();
		coordVec.push_back(p1);
		for(int i=1; i<numOut; i++){
			double x = p1.x+stepX*i;
			double m = d1*d1-(x-c1.x)*(x-c1.x);
			double n = d2*d2-(x-c2.x)*(x-c2.x);
			double a = (c1.z-c2.z)/(c2.y-c1.y);
			double b = (c2.y*c2.y - c1.y*c1.y + c2.z*c2.z -
					c1.z*c1.z + m - n )/(2*c2.y-2*c1.y);
			double A = a*a+1;
			double B = 2*a*(b-c1.y)-2*c1.z;
			double C = (b-c1.y)*(b-c1.y)+c1.z*c1.z-m;

			double z1 = (-B+sqrt(B*B-4*A*C))/(2*A);
			double z2 = (-B-sqrt(B*B-4*A*C))/(2*A);
			double y1 = a*z1+b;
			double y2 = a*z2+b;
			double y,z;

			Coord co1(x, y1, z1);
			Coord co2(x, y2, z2);
			double dis1 = getCoordDis( co1, p1 );
			double dis2 = getCoordDis( co2, p1 );
			if( dis1 < dis2 ){
				y = y1;
				z = z1;
			}else{
				y = y2;
				z = z2;
			}
			Coord co(x, y, z);
//			co1.print();
//			co2.print();
			coordVec.push_back(co);
		}
	}else if(  stepY >= stepX && stepY >= stepZ ){
		double d1 = getCoordDis( a1.get_coord(), p1 );
		double d2 = getCoordDis( a2.get_coord(), p1 );
		Coord c1 = a1.get_coord();
		Coord c2 = a2.get_coord();
		coordVec.push_back(p1);
		for(int i=1; i<numOut; i++){
			double y = p1.y+stepY*i;
			double m = d1*d1-(y-c1.y)*(y-c1.y);
			double n = d2*d2-(y-c2.y)*(y-c2.y);
			double a = (c1.z-c2.z)/(c2.x-c1.x);
			double b = (c2.x*c2.x - c1.x*c1.x + c2.z*c2.z -
					c1.z*c1.z + m - n )/(2*c2.x-2*c1.x);
			double A = a*a+1;
			double B = 2*a*(b-c1.x)-2*c1.z;
			double C = (b-c1.x)*(b-c1.x)+c1.z*c1.z-m;

			double z1 = (-B+sqrt(B*B-4*A*C))/(2*A);
			double z2 = (-B-sqrt(B*B-4*A*C))/(2*A);
			double x1 = a*z1+b;
			double x2 = a*z2+b;
			double x,z;

			Coord co1(x1, y, z1);
			Coord co2(x2, y, z2);
			double dis1 = getCoordDis( co1, p1 );
			double dis2 = getCoordDis( co2, p1 );
			if( dis1 < dis2 ){
				x = x1;
				z = z1;
			}else{
				x = x2;
				z = z2;
			}
			Coord co(x, y, z);
//			co1.print();
//			co2.print();
			coordVec.push_back(co);
		}
	}else 	if( stepZ >= stepX && stepZ >= stepY ){
		double d1 = getCoordDis( a1.get_coord(), p1 );
		double d2 = getCoordDis( a2.get_coord(), p1 );
		Coord c1 = a1.get_coord();
		Coord c2 = a2.get_coord();
		coordVec.push_back(p1);
		for(int i=1; i<numOut; i++){
			double z = p1.z+stepZ*i;
			double m = d1*d1-(z-c1.z)*(z-c1.z);
			double n = d2*d2-(z-c2.z)*(z-c2.z);
			double a = (c1.x-c2.x)/(c2.y-c1.y);
			double b = (c2.y*c2.y - c1.y*c1.y + c2.x*c2.x -
					c1.x*c1.x + m - n )/(2*c2.y-2*c1.y);
			double A = a*a+1;
			double B = 2*a*(b-c1.y)-2*c1.x;
			double C = (b-c1.y)*(b-c1.y)+c1.x*c1.x-m;

			double x1 = (-B+sqrt(B*B-4*A*C))/(2*A);
			double x2 = (-B-sqrt(B*B-4*A*C))/(2*A);
			double y1 = a*x1+b;
			double y2 = a*x2+b;
			double y,x;

			Coord co1(x1, y1, z);
			Coord co2(x2, y2, z);
			double dis1 = getCoordDis( co1, p1 );
			double dis2 = getCoordDis( co2, p1 );
			if( dis1 < dis2 ){
				y = y1;
				x = x1;
			}else{
				y = y2;
				x = x2;
			}
			Coord co(x,y,z);
//			co1.print();
//			co2.print();
			coordVec.push_back(co);
		}
	}
//	cout<<"----------"<<endl;
//	cout<<"reduced surface arch"<<endl;
//	for( size_t i=0; i<coordVec.size(); i++ ){
//		coordVec[i].print();
//	}
	return coordVec;
}

vector<Coord>
getReducedSurfaceArch( const reducedSurface_connect& rsc, const int& num ){
	Atom a1 = rsc.get_at1();
	Atom a2 = rsc.get_at2();
	Atom a3 = rsc.get_at3();
	Atom a4 = rsc.get_at4();

//	Coord c1 = a1.get_coord();
//	Coord c2 = a2.get_coord();
	Coord c3 = a3.get_coord();
	Coord c4 = a4.get_coord();

	Coord p1 = rsc.get_probeCoordA();
	Coord p2 = rsc.get_probeCoordB();
	Coord center = getTriangleVerticalPoint(a3.get_coord(), a4.get_coord(), p1);
//	cout<<"center:";
//	center.print();
	double dis_p1center = getCoordDis(p1, center);
	double dis_p1a1 = getCoordDis(p1, a1.get_coord());
	double dis_p2a2 = getCoordDis(p2, a2.get_coord());
	double dis_p1a3 = getCoordDis(p1, a3.get_coord());
	double dis_p1a4 = getCoordDis(p1, a4.get_coord());

	double angle_p1p2 = getAngle(center, p1, p2);
	double angle_a1a2 = getAngle(center, a1.get_coord(), a2.get_coord());
//	cout<<"angle_p1p2: "<<angle_p1p2<<endl;
//	cout<<"angle_a1a2: "<<angle_a1a2<<endl;
	double angleRange = 0;
	if(angle_p1p2 > angle_a1a2){
		angleRange = (2*PI) - angle_p1p2;
	}else{
		angleRange = angle_p1p2;
	}
	if(0){
		if( angle_p1p2>PI ){
			angleRange = 2*PI-angle_p1p2;
		}else{
			angleRange = angle_p1p2;
		}
	}

	int numOut = angleRange*180/3.14159;
//	int numOut = 10;

	vector<Coord> coVec;

	double angleStep = angleRange/(numOut-1);
	coVec.push_back(p1);
	double degree=PI/180.0;
	Coord preCoord = p1;
//	for( double angle=degree; angle<angleRange; angle=angle+degree ){
	for( int i=1; i<numOut-1; i++ ){
		double angle = i*angleStep;
		double cos_angle = cos(angle);
		double a = (center.y-p1.y)/(p1.x-center.x);
		double b = (center.z-p1.z)/(p1.x-center.x);
		double c = ( center.x*(p1.x-center.x)+center.y*(p1.y-center.y)+center.z*(p1.z-center.z)
				+ dis_p1center*dis_p1center*cos_angle)/( p1.x-center.x );

		double m = -(2*b*(c-c3.x)-2*b*(c-c4.x)+2*c4.z-2*c3.z)/( 2*a*(c-c3.x)-2*a*(c-c4.x)+2*c4.y-2*c3.y );
		double n = ( -(c-c3.x)*(c-c3.x)+(c-c4.x)*(c-c4.x)-c3.y*c3.y+c4.y*c4.y-c3.z*c3.z+c4.z*c4.z
				+ dis_p1a3*dis_p1a3 - dis_p1a4*dis_p1a4)/ ( 2*a*(c-c3.x)-2*a*(c-c4.x)+2*c4.y-2*c3.y );

		double A = (a*m+b)*(a*m+b)+m*m+1;
		double B = 2*(a*n+c-c3.x)*(a*m+b)+2*(n-c3.y)*m-2*c3.z;
		double Y = (a*n+c-c3.x)*(a*n+c-c3.x) + (n-c3.y)*(n-c3.y)+c3.z*c3.z-dis_p1a3*dis_p1a3;

		double z1 = (-B + sqrt(B*B-4*A*Y))/(2*A);
		double z2 = (-B - sqrt(B*B-4*A*Y))/(2*A);

		double y1 = m*z1+n;
		double y2 = m*z2+n;

		double x1 = a*y1+b*z1+c;
		double x2 = a*y2+b*z2+c;

		Coord c1(x1, y1, z1);
		Coord c2(x2, y2, z2);
		Coord coo;

		double dis1 = getCoordDis( c1, a1.get_coord() );
		double dis2 = getCoordDis( c1, a2.get_coord() );

		double anglePre = getAngle( center, c1, preCoord );
//		cout<<"anglePre:"<<anglePre<<endl;
//		cout<<"angleStep:"<<angleStep<<endl;
		if( dis1 >= dis_p1a1 && dis2 >= dis_p2a2 && abs(anglePre-angleStep)<0.001 ){
			coo = c1;
			preCoord = c1;
		}else{
			coo = c2;
			preCoord = c2;
		}

		coVec.push_back(coo);
	}
	coVec.push_back(p2);
	if(0){
		cout<<"p1:";
		p1.print();
		for( size_t i=0; i<coVec.size(); i++ ){
			cout<<i<<" ";
			coVec[i].print();
		}
		cout<<"p2:";
		p2.print();
	}
	return coVec;
}

/*
 * get the SES arch on the two contact atoms a1, a2;
 * stored in out1, out2
 */
void
get_SES_arch( const Atom& a1, const Atom& a2,
		const Coord& p1, const Coord& p2, const int& numOut,
		vector<Coord>& out1, vector<Coord>& out2){
	vector<Coord> reducedArch = getReducedSurfaceArch(a1, a2, p1, p2, numOut);
	vector<Coord> o1, o2;
	Coord c1 = a1.get_coord();
	Coord c2 = a2.get_coord();
	double d1 = getCoordDis( c1, p1 );
	double d2 = getCoordDis( c2, p1 );
	double rate1 = a1.get_vdwRadius()/d1;
	double rate2 = a2.get_vdwRadius()/d2;
	out1.clear();
	out2.clear();
	for( size_t i=0; i<reducedArch.size(); i++ ){
		Coord p = reducedArch[i];
		double x1 = (p.x-c1.x)*rate1 + c1.x;
		double y1 = (p.y-c1.y)*rate1 + c1.y;
		double z1 = (p.z-c1.z)*rate1 + c1.z;

		double x2 = (p.x-c2.x)*rate2 + c2.x;
		double y2 = (p.y-c2.y)*rate2 + c2.y;
		double z2 = (p.z-c2.z)*rate2 + c2.z;

		Coord co1((double)x1, (double)y1, (double)z1);
		Coord co2((double)x2, (double)y2, (double)z2);

		out1.push_back(co1);
		out2.push_back(co2);
	}
}

void
get_SES_toric_reentrant_face( reducedSurface_connect& rsc,
		const int& numOut, vector<Coord>& coordVec, vector<Coord>& normalVec ){
	Atom a1 = rsc.get_at3();
	Atom a2 = rsc.get_at4();

	Coord p1 = rsc.get_probeCoordA();
	Coord p2 = rsc.get_probeCoordB();
	coordVec.clear();
	normalVec.clear();

	vector<Coord> reducedArch = rsc.get_coordVec();

	Coord c1 = a1.get_coord();
	Coord c2 = a2.get_coord();
	double probeRadius = getCoordDis(c1, p1)-a1.get_vdwRadius();

	bool archCrash = false;
	int  crashIndex = 0;
	vector<Coord> arch1 = getArchPoints(c2, c1, reducedArch.front(), probeRadius, numOut );
	vector<Coord> arch2 = getArchPoints(c2, c1, reducedArch.back(), probeRadius, numOut );

	for( size_t i=1; i<arch1.size()-1; i++ ){

		double dis0 = getCoordDis( arch1[i-1], arch2[i-1] );
		double dis1 = getCoordDis( arch1[i], arch2[i] );
		double dis2 = getCoordDis( arch1[i+1], arch2[i+1] );
		double  val = (dis1-dis0)*(dis2-dis1);

		if( val < 0 && abs(dis1-dis0)>0.01 && abs(dis2-dis1)>0.01 ){
//			cout<<"crash!!!"<<endl;
//			rsc.print();
			archCrash = true;
			crashIndex = i;
			rsc.probeCrash = true;
			vector<Coord> cVec;
			vector<Coord> cVec1 = rsc.getCrashProbeCoordConcave(arch1[i], arch1[ arch1.size()-i ] );
			vector<Coord> cVec2 = rsc.getCrashProbeCoordConvex(arch1[i], arch1[ arch1.size()-i ] );
			if( getCoordDis( cVec1[cVec1.size()/2], rsc.get_at1().get_coord() ) <
				getCoordDis( cVec2[cVec1.size()/2], rsc.get_at1().get_coord() ) ){
				cVec = cVec1;
			}else{
				cVec = cVec2;
			}
			vector<Coord> aVec, bVec;

			for( size_t j=0; j<i; j++ ){
				aVec.push_back( arch1[j] );
				bVec.push_back( arch2[j] );
			}
			aVec.insert( aVec.end(), cVec.begin(), cVec.end() );
			bVec.insert( bVec.end(), cVec.begin(), cVec.end() );
			for( size_t j=arch1.size()-i; j<arch1.size(); j++ ){
				aVec.push_back( arch1[j] );
			}
			for( size_t j=arch2.size()-i; j<arch2.size(); j++ ){
				bVec.push_back( arch2[j] );
			}

//			rsc.set_crashCoordA(aVec);
//			rsc.set_crashCoordB(bVec);
			break;
		}
	}


	for( size_t i=0; i<reducedArch.size()-1; i++ ){
//	for( size_t i=reducedArch.size()-5; i<reducedArch.size()-1; i++ ){

		// arch of toric reentrant, this is the arch between atom a1, a2, not on a1 or a2
		Coord probe1 = reducedArch[i];
		vector<Coord> archP1 = getArchPoints(c1, c2, probe1, probeRadius, numOut );
		Coord probe2 = reducedArch[i+1];
		vector<Coord> archP2 = getArchPoints(c1, c2, probe2, probeRadius, numOut );
//		cout<<"archP1:"<<archP1.size()<<endl;
//		cout<<"archP2:"<<archP2.size()<<endl;

		if( !archCrash ){
			for( size_t j=0; j<archP1.size()-1; j++ ){
				Coord c11 = archP1[j];
				Coord c21 = archP2[j];
				Coord c12 = archP1[j+1];
				Coord c22 = archP2[j+1];
				Coord n11 = (probe1-archP1[j]).getNormalizeCoord();
				Coord n21 = (probe2-archP2[j]).getNormalizeCoord();
				Coord n12 = (probe1-archP1[j+1]).getNormalizeCoord();
				Coord n22 = (probe2-archP2[j+1]).getNormalizeCoord();

				coordVec.push_back(c11);
				coordVec.push_back(c21);
				coordVec.push_back(c22);
				coordVec.push_back(c11);
				coordVec.push_back(c12);
				coordVec.push_back(c22);
				normalVec.push_back(n11);
				normalVec.push_back(n21);
				normalVec.push_back(n22);
				normalVec.push_back(n11);
				normalVec.push_back(n12);
				normalVec.push_back(n22);
			}
		}else{
			for( int j=0; j<crashIndex; j++ ){
				Coord c11 = archP1[j];
				Coord c21 = archP2[j];
				Coord c12 = archP1[j+1];
				Coord c22 = archP2[j+1];
				Coord n11 = (probe1-archP1[j]).getNormalizeCoord();
				Coord n21 = (probe2-archP2[j]).getNormalizeCoord();
				Coord n12 = (probe1-archP1[j+1]).getNormalizeCoord();
				Coord n22 = (probe2-archP2[j+1]).getNormalizeCoord();

				coordVec.push_back(c11);
				coordVec.push_back(c21);
				coordVec.push_back(c22);
				coordVec.push_back(c11);
				coordVec.push_back(c12);
				coordVec.push_back(c22);
				normalVec.push_back(n11);
				normalVec.push_back(n21);
				normalVec.push_back(n22);
				normalVec.push_back(n11);
				normalVec.push_back(n12);
				normalVec.push_back(n22);
			}
			for( size_t j=archP1.size()-crashIndex-2; j<archP1.size()-1; j++ ){
				Coord c11 = archP1[j];
				Coord c21 = archP2[j];
				Coord c12 = archP1[j+1];
				Coord c22 = archP2[j+1];
				Coord n11 = (probe1-archP1[j]).getNormalizeCoord();
				Coord n21 = (probe2-archP2[j]).getNormalizeCoord();
				Coord n12 = (probe1-archP1[j+1]).getNormalizeCoord();
				Coord n22 = (probe2-archP2[j+1]).getNormalizeCoord();

				coordVec.push_back(c11);
				coordVec.push_back(c21);
				coordVec.push_back(c22);
				coordVec.push_back(c11);
				coordVec.push_back(c12);
				coordVec.push_back(c22);
				normalVec.push_back(n11);
				normalVec.push_back(n21);
				normalVec.push_back(n22);
				normalVec.push_back(n11);
				normalVec.push_back(n12);
				normalVec.push_back(n22);
			}
		}
	}
}

/**
 * p1, p2 are on the sphere with 'center' and 'raidus'
 * this function is to get the points on arch p1-p2
 */
vector<Coord>
get_sphere_arch( const Coord& center, const double radius,
		const Coord& p1, const Coord& p2, const int numOut ){

	vector<Coord> cooVec;
	double angleRange = getAngle( center, p1, p2 );
	double angleStep = angleRange/(double(numOut));

	Coord pp1 = p1-center;
	Coord pp2 = p2-center;
	Coord axis = cross(pp1, pp2);

	cooVec.push_back(p1);
//	for( double angle = angleStep; angle<=angleRange; angle+=angleStep ){
	for( int i=1; i<numOut-1; i++ ){
		double angle = i*angleStep;
		Coord p = rotateAroundAxis(pp1, axis, angle );
		Coord c = p+center;
		cooVec.push_back(c);
	}
	cooVec.push_back(p2);

	if(0){
		cout<<"@@@@@@@@@@@"<<endl;
		cout<<"numOut:"<<numOut<<endl;
		for( size_t i=0; i<cooVec.size(); i++ ){
			cout<<i<<" ";
			cooVec[i].print();
		}
	}

	return cooVec;
}

void
get_SES_spheric_reentrant_face(const sphereSurface& ss,
		vector<Coord>& coordVec, vector<Coord>& normalVec){
	Atom a1 = ss.get_at1();
	Atom a2 = ss.get_at2();
	Atom a3 = ss.get_at3();
	Coord probe = ss.get_probe();

	Coord c1 = a1.get_coord();
	Coord c2 = a2.get_coord();
	Coord c3 = a3.get_coord();
	double probeRadius = getCoordDis(c1, probe)-a1.get_vdwRadius();

	double dis12 = getCoordDis(c1,c2);
	double dis13 = getCoordDis(c1,c3);
	double dis23 = getCoordDis(c3,c2);

	vector<Coord> archP1,archP2;

	if( dis12 <= dis13 && dis12 <= dis23 ){
		archP1 = inverseCoordVec(ss.get_arch13());
		archP2 = inverseCoordVec(ss.get_arch23());
	}else if( dis13 <= dis23 && dis13 <= dis12 ){

		archP1 = inverseCoordVec(ss.get_arch12());
		archP2 = ss.get_arch23();
	}else	if( dis23 <= dis13 && dis23 <= dis12 ){
		archP1 = ss.get_arch12();
		archP2 = ss.get_arch13();
	}

	vector<Coord> line1;
	line1.push_back(archP1.front());



	for( int i=1; i<archP1.size(); i++ ){
//	for( int i=1; i<10; i++ ){

		size_t p2Index1 = (i-1)*(archP2.size()-1)/(archP1.size()-1);
		size_t p2Index2 = i*(archP2.size()-1)/(archP1.size()-1);

		line1 = get_sphere_arch( probe, probeRadius, archP1[i-1], archP2[p2Index1], i );

		vector<Coord> line2 = get_sphere_arch( probe, probeRadius,
				archP1[i], archP2[p2Index2], i+1 );

		for( size_t j=0; j<line2.size()-1; j++ ){

			coordVec.push_back( line2[j]);
			coordVec.push_back( line2[j+1]);
			coordVec.push_back( line1[j] );
			Coord n20 = (line2[j]-probe).getNormalizeCoord();
			Coord n21 = (line2[j+1]-probe).getNormalizeCoord();
			Coord n10 = (line1[j]-probe).getNormalizeCoord();
			normalVec.push_back(n20);
			normalVec.push_back(n21);
			normalVec.push_back(n10);
		}
		for( size_t j=0; j<line1.size()-1; j++ ){
			coordVec.push_back( line1[j]);
			coordVec.push_back( line1[j+1]);
			coordVec.push_back( line2[j+1] );
			Coord n10 = (line1[j]-probe).getNormalizeCoord();
			Coord n11 = (line1[j+1]-probe).getNormalizeCoord();
			Coord n21 = (line2[j+1]-probe).getNormalizeCoord();
			normalVec.push_back(n10);
			normalVec.push_back(n11);
			normalVec.push_back(n21);
		}
	}
}

void
get_SES_spheric_reentrant_face(  const Atom& a1, const Atom& a2, const Atom& a3,
		const Coord& probe, const int& numOut,
		vector<Coord>& coordVec, vector<Coord>& normalVec){
	Coord c1 = a1.get_coord();
	Coord c2 = a2.get_coord();
	Coord c3 = a3.get_coord();
	double probeRadius = getCoordDis(c1, probe)-a1.get_vdwRadius();

	double dis12 = getCoordDis(c1,c2);
	double dis13 = getCoordDis(c1,c3);
	double dis23 = getCoordDis(c3,c2);

	vector<Coord> archP1,archP2;
	if( dis12 <= dis13 && dis12 <= dis23 ){

		archP1 = getArchPoints( c3, c1, probe, probeRadius, numOut );
		archP2 = getArchPoints( c3, c2, probe, probeRadius, numOut );
	}else if( dis13 <= dis23 && dis13 <= dis12 ){

		archP1 = getArchPoints( c2, c1, probe, probeRadius, numOut );
		archP2 = getArchPoints( c2, c3, probe, probeRadius, numOut );
	}else	if( dis23 <= dis13 && dis23 <= dis12 ){

		archP1 = getArchPoints( c1, c2, probe, probeRadius, numOut );
		archP2 = getArchPoints( c1, c3, probe, probeRadius, numOut );
	}

//	vector<Coord> archP1 = getArchPoints( c1, c2, probe, probeRadius, numOut );
//	vector<Coord> archP2 = getArchPoints( c1, c3, probe, probeRadius, numOut );
	vector<Coord> line1;
	line1.push_back(archP1.front());


	for( int i=1; i<archP1.size(); i++ ){
//	for( int i=archP1.size()-1; i<archP1.size(); i++ ){

		size_t p2Index1 = (i-1)*(archP2.size()-1)/(archP1.size()-1);
		size_t p2Index2 = i*(archP2.size()-1)/(archP1.size()-1);

		line1 = get_sphere_arch( probe, probeRadius, archP1[i-1], archP2[p2Index1], i );

		vector<Coord> line2 = get_sphere_arch( probe, probeRadius, archP1[i], archP2[p2Index2], i+1 );

		for( size_t j=0; j<line2.size()-1; j++ ){

			coordVec.push_back( line2[j]);
			coordVec.push_back( line2[j+1]);
			coordVec.push_back( line1[j] );
			Coord n20 = (line2[j]-probe).getNormalizeCoord();
			Coord n21 = (line2[j+1]-probe).getNormalizeCoord();
			Coord n10 = (line1[j]-probe).getNormalizeCoord();
			normalVec.push_back(n20);
			normalVec.push_back(n21);
			normalVec.push_back(n10);
//			normalVec.push_back(n21);

//			cout<<"i:"<<i<<" j:"<<j<<endl;
//			line2[j].print();
//			line2[j+1].print();
//			line1[j].print();
//			cout<<endl;
		}
//		cout<<"##"<<endl;
		for( size_t j=0; j<line1.size()-1; j++ ){
			coordVec.push_back( line1[j]);
			coordVec.push_back( line1[j+1]);
			coordVec.push_back( line2[j+1] );
			Coord n10 = (line1[j]-probe).getNormalizeCoord();
			Coord n11 = (line1[j+1]-probe).getNormalizeCoord();
			Coord n21 = (line2[j+1]-probe).getNormalizeCoord();
			normalVec.push_back(n10);
			normalVec.push_back(n11);
			normalVec.push_back(n21);
//			normalVec.push_back(n11);

//			line1[j].print();
//			line1[j+1].print();
//			line2[j+1].print();
//			cout<<endl;
		}
//		line1 = line2;
	}
}

bool
reducedSphereShare( const reducedSurface& r1, const reducedSurface& r2,
		Atom& a, Atom& b, Atom& c, Atom&d ){
	Atom a1 = r1.get_at1();
	Atom a2 = r1.get_at2();
	Atom a3 = r1.get_at3();

	Atom aa = r2.get_at1();
	Atom ab = r2.get_at2();
	Atom ac = r2.get_at3();

	vector<Atom> sVec;

	int count = 0;
	if( a1 == aa ){
		count++;
		sVec.push_back(a1);
		sVec.push_back(aa);
	}
	if( a1 == ab ){
		count++;
		sVec.push_back(a1);
		sVec.push_back(ab);
	}
	if( a1 == ac ){
		count++;
		sVec.push_back(a1);
		sVec.push_back(ac);
	}
	if( a2 == aa ){
		count++;
		sVec.push_back(a2);
		sVec.push_back(aa);
	}
	if( a2 == ab ){
		count++;
		sVec.push_back(a2);
		sVec.push_back(ab);
	}
	if( a2 == ac ){
		count++;
		sVec.push_back(a2);
		sVec.push_back(ac);
	}
	if( a3 == aa ){
		count++;
		sVec.push_back(a3);
		sVec.push_back(aa);
	}
	if( a3 == ab ){
		count++;
		sVec.push_back(a3);
		sVec.push_back(ab);
	}
	if( a3 == ac ){
		count++;
		sVec.push_back(a3);
		sVec.push_back(ac);
	}

	if( count == 2 ){
		if( find( sVec.begin(), sVec.end(), a1 ) == sVec.end() ){
			a = a1;
		}
		if( find( sVec.begin(), sVec.end(), a2 ) == sVec.end() ){
			a = a2;
		}
		if( find( sVec.begin(), sVec.end(), a3 ) == sVec.end() ){
			a = a3;
		}
		if( find( sVec.begin(), sVec.end(), aa ) == sVec.end() ){
			b = aa;
		}
		if( find( sVec.begin(), sVec.end(), ab ) == sVec.end() ){
			b = ab;
		}
		if( find( sVec.begin(), sVec.end(), ac ) == sVec.end() ){
			b = ac;
		}

		c = sVec[0];
		d = sVec[2];
		return 1;
	}else if(count==3){
//		if()
		return 0;
	}else{
		return 0;
	}
}

vector<reducedSurface_connect>
get_reducedSurface_connect( const vector<reducedSurface>& rsVec ){

	vector<reducedSurface_connect> rcVec;


	for( size_t i=0; i<rsVec.size(); i++ ){
		for( size_t j=0; j<rsVec.size(); j++ ){
			if( i!=j ){
				Atom a, b, c, d;
				if( reducedSphereShare(rsVec[i], rsVec[j], a, b, c, d) ){
					reducedSurface_connect r(a,b,c,d,
							rsVec[i].get_probeCoord(),
							rsVec[j].get_probeCoord());
					bool flag = true;
					for( size_t k=0; k<rcVec.size(); k++ ){
						if( ( rcVec[k].get_at1() == a || rcVec[k].get_at1() == b ) &&
							( rcVec[k].get_at2() == a || rcVec[k].get_at2() == b ) &&
							( rcVec[k].get_at3() == c || rcVec[k].get_at3() == d ) &&
							( rcVec[k].get_at4() == c || rcVec[k].get_at4() == d ) &&
							( rcVec[k].get_probeCoordA() == rsVec[i].get_probeCoord() ||
							  rcVec[k].get_probeCoordA() == rsVec[j].get_probeCoord()) &&
							( rcVec[k].get_probeCoordB() == rsVec[i].get_probeCoord() ||
							  rcVec[k].get_probeCoordB() == rsVec[j].get_probeCoord() )){
							flag = false;
						}
					}
					if( flag ){
						rcVec.push_back(r);
					}
				}
			}
		}
	}

	return rcVec;
}

/**
 * the vertical point 'c' on line c1-c2 that makes c-c3 vertical to line c1-c2
 */
Coord
get_vertical_point( const Coord& c1, const Coord& c2, const Coord& c3 ){
	double a = (c2.z-c1.z)/(c2.x-c1.x);
	double b = c1.z-c1.x*(c2.z-c1.z)/(c2.x-c1.x);
	double c = (c2.y-c1.y)/(c2.x-c1.x);
	double d = c1.y-c1.x*(c2.y-c1.y)/(c2.x-c1.x);
	double A = 1+a*a+c*c;
	double B = -c1.x-c3.x+c*(2*d-c1.y-c3.y)+a*(2*b-c1.z-c3.z);
	double C = c1.x*c3.x+(d-c1.y)*(d-c3.y)+(b-c1.z)*(b-c3.z);

	double x1 = (-B+sqrt(B*B-4*A*C))/(2*A);
	double x2 = (-B-sqrt(B*B-4*A*C))/(2*A);
	double y1 = a*x1+b;
	double y2 = a*x2+b;
	double z1 = c*x1+d;
	double z2 = c*x2+d;

	double x,y,z;
	if(x1==c1.x || x1==c2.x){
		x=x2;
		y=y2;
		z=z2;
	}else{
		x=x1;
		y=y1;
		z=z1;
	}
	Coord co(x,y,z);
	return co;
}

