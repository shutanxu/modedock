/*
 * stick.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: stan
 */




#include"stick.h"

QVector3D
rotateAroundAxis( const QVector3D input, const QVector3D axis, const double angle ){
	double dis = sqrt(axis.x()*axis.x()+axis.y()*axis.y()+axis.z()*axis.z());
	if( dis == 0 ){
		cerr<<"!!!no axis for rotating around"<<endl;
		return input;
	}
	QVector3D  u((double)axis.x()/dis, (double)axis.y()/dis, (double)axis.z()/dis);

	double a11 = cos(angle)+u.x()*u.x()*(1-cos(angle));
	double a12 = u.x()*u.y()*(1-cos(angle))-u.z()*sin(angle);
	double a13 = u.x()*u.z()*(1-cos(angle))+u.y()*sin(angle);
	double a21 = u.y()*u.x()*(1-cos(angle))+u.z()*sin(angle);
	double a22 = cos(angle)+u.y()*u.y()*(1-cos(angle));
	double a23 = u.y()*u.z()*(1-cos(angle))-u.x()*sin(angle);
	double a31 = u.z()*u.x()*(1-cos(angle))-u.y()*sin(angle);
	double a32 = u.z()*u.y()*(1-cos(angle))+u.x()*sin(angle);
	double a33 = cos(angle)+u.z()*u.z()*(1-cos(angle));

	double x = a11*input.x()+a12*input.y()+a13*input.z();
	double y = a21*input.x()+a22*input.y()+a23*input.z();
	double z = a31*input.x()+a32*input.y()+a33*input.z();

	QVector3D c(x, y, z);
	return c;
}

bool
getStick(const QVector3D p1, const QVector3D p2,
		const QColor c1, const QColor c2,
		const float ratio, const float radius,
		vector<QVector3D>& trigVertex, vector<QVector3D>& trigNormal,
		vector<QColor>& trigColor ){

	float stickRadius = radius;
	trigVertex.clear();
	trigNormal.clear();
	trigColor.clear();
	QVector3D mid = p1 + (p2-p1)*ratio;
	QVector3D s1,s2, sMid;

	if( p1.z() != p2.z() && p1.y() != p2.y() ){
		s1.setX(p1.x());
		double z=sqrt(stickRadius*stickRadius*(p2.y()-p1.y())*(p2.y()-p1.y())/
				( (p2.z()-p1.z())*(p2.z()-p1.z())+(p2.y()-p1.y())*(p2.y()-p1.y()) ))+p1.z() ;
		s1.setZ(z);
		s1.setY( p1.y() - (z-p1.z())*(p2.z()-p1.z())/(p2.y()-p1.y()) );
	}else if( p1.z() != p2.z() && p2.x()!=p1.x() ){
		s1.setY(p1.y());
		double z=sqrt(stickRadius*stickRadius*(p2.x()-p1.x())*(p2.x()-p1.x())/
				( (p2.z()-p1.z())*(p2.z()-p1.z())+(p2.x()-p1.x())*(p2.x()-p1.x()) ))+p1.z() ;
		s1.setZ(z);
		s1.setX( p1.x() - (z-p1.z())*(p2.z()-p1.z())/(p2.x()-p1.x()) );
	}else if( p1.y() != p2.y() && p2.x()!=p1.x() ){
		s1.setZ(p1.z());
		double y=sqrt(stickRadius*stickRadius*(p2.x()-p1.x())*(p2.x()-p1.x())/
				( (p2.y()-p1.y())*(p2.y()-p1.y())+(p2.x()-p1.x())*(p2.x()-p1.x()) ))+p1.y() ;
		s1.setY(y);
		s1.setX( p1.x() - (y-p1.y())*(p2.y()-p1.y())/(p2.x()-p1.x()) );
	}

	s2.setX( s1.x() + p2.x()-p1.x() );
	s2.setY( s1.y() + p2.y()-p1.y() );
	s2.setZ( s1.z() + p2.z()-p1.z() );

	//stick
	sMid=s1+(s2-s1)*ratio;
	QVector3D axis=p2-p1;
	QVector3D start1=s1;
	QVector3D start2=s2;

	for( double angle=0; angle<=360.0; angle+=20.0 ){
		double ang=angle*3.14159/180.0;
		QVector3D t1=rotateAroundAxis((start1-p1), axis, ang )+p1;
		QVector3D t2=rotateAroundAxis((start2-p1), axis, ang )+p1;
		QVector3D tMid=t1+(t2-t1)*ratio;

		QVector3D nt1=t1-p1;
		QVector3D nt2=t2-p2;
		QVector3D ns1=s1-p1;
		QVector3D ns2=s2-p2;
		QVector3D ntMid=tMid-mid;
		QVector3D nsMid=sMid-mid;

		trigVertex.push_back(t1);
		trigVertex.push_back(s1);
		trigVertex.push_back(sMid);
		trigNormal.push_back(-nt1);
		trigNormal.push_back(-ns1);
		trigNormal.push_back(-nsMid);
		trigColor.push_back(c1);
		trigColor.push_back(c1);
		trigColor.push_back(c1);

		trigVertex.push_back(t1);
		trigVertex.push_back(sMid);
		trigVertex.push_back(tMid);
		trigNormal.push_back(-nt1);
		trigNormal.push_back(-nsMid);
		trigNormal.push_back(-ntMid);
		trigColor.push_back(c1);
		trigColor.push_back(c1);
		trigColor.push_back(c1);

		trigVertex.push_back(tMid);
		trigVertex.push_back(sMid);
		trigVertex.push_back(s2);
		trigNormal.push_back(-ntMid);
		trigNormal.push_back(-nsMid);
		trigNormal.push_back(-ns2);
		trigColor.push_back(c2);
		trigColor.push_back(c2);
		trigColor.push_back(c2);

		trigVertex.push_back(tMid);
		trigVertex.push_back(s2);
		trigVertex.push_back(t2);
		trigNormal.push_back(-ntMid);
		trigNormal.push_back(-ns2);
		trigNormal.push_back(-nt2);
		trigColor.push_back(c2);
		trigColor.push_back(c2);
		trigColor.push_back(c2);

		s1=t1;
		s2=t2;
		sMid=tMid;
	}

	//head1
	s1=start1;
    for( double angle=0; angle<= 180.0 ; angle+=20.0 ){
		double ang=angle*3.14159/180.0;
		QVector3D t1=rotateAroundAxis((start1-p1), axis, ang )+p1;
		QVector3D axis1=QVector3D::crossProduct(p2-p1,s1-p1);
		QVector3D axis2=QVector3D::crossProduct(p2-p1,t1-p1);

		QVector3D ss1=s1;
		QVector3D tt1=t1;

		for( double angle2=0; angle2<=180; angle2+=30 ){
			double ang2=angle2*3.14159/180.0;
			QVector3D v1=rotateAroundAxis((s1-p1),axis1,ang2)+p1;
			QVector3D v2=rotateAroundAxis((t1-p1),axis2,ang2)+p1;
			QVector3D n0=ss1-p1;
			QVector3D n1=tt1-p1;
			QVector3D n2=v1-p1;
			QVector3D n3=v2-p1;

            if( angle2 <= 90 ){
				trigVertex.push_back(ss1);
				trigVertex.push_back(tt1);
				trigVertex.push_back(v1);
				trigNormal.push_back(-n0);
				trigNormal.push_back(-n1);
				trigNormal.push_back(-n2);
				trigColor.push_back(c1);
				trigColor.push_back(c1);
				trigColor.push_back(c1);

				trigVertex.push_back(v1);
				trigVertex.push_back(tt1);
				trigVertex.push_back(v2);
				trigNormal.push_back(-n2);
				trigNormal.push_back(-n1);
				trigNormal.push_back(-n3);
				trigColor.push_back(c1);
				trigColor.push_back(c1);
				trigColor.push_back(c1);
			}else{
				trigVertex.push_back(ss1);
				trigVertex.push_back(v1);
				trigVertex.push_back(tt1);

				trigNormal.push_back(-n0);
				trigNormal.push_back(-n2);
				trigNormal.push_back(-n1);

				trigColor.push_back(c1);
				trigColor.push_back(c1);
				trigColor.push_back(c1);

				trigVertex.push_back(v1);
				trigVertex.push_back(v2);
				trigVertex.push_back(tt1);

				trigNormal.push_back(-n2);
				trigNormal.push_back(-n3);
				trigNormal.push_back(-n1);

				trigColor.push_back(c1);
				trigColor.push_back(c1);
				trigColor.push_back(c1);
			}
			ss1=v1;
			tt1=v2;
		}
        s1=t1;
	}

	//head2
	s1=start2;
	for( double angle=0; angle<=180.0; angle+=20.0 ){
		double ang=angle*3.14159/180.0;
		QVector3D t1=rotateAroundAxis((start2-p2), axis, ang )+p2;
		QVector3D axis1=QVector3D::crossProduct(p1-p2,s1-p2);
		QVector3D axis2=QVector3D::crossProduct(p1-p2,t1-p2);

		QVector3D ss1=s1;
		QVector3D tt1=t1;

		for( double angle2=0; angle2<=180; angle2+=30 ){
			double ang2=angle2*3.14159/180.0;
			QVector3D v1=rotateAroundAxis((s1-p2),axis1,ang2)+p2;
			QVector3D v2=rotateAroundAxis((t1-p2),axis2,ang2)+p2;
			QVector3D n0=ss1-p2;
			QVector3D n1=tt1-p2;
			QVector3D n2=v1-p2;
			QVector3D n3=v2-p2;

			if( angle2>90 ){
				trigVertex.push_back(ss1);
				trigVertex.push_back(tt1);
				trigVertex.push_back(v1);
				trigNormal.push_back(-n0);
				trigNormal.push_back(-n1);
				trigNormal.push_back(-n2);
				trigColor.push_back(c2);
				trigColor.push_back(c2);
				trigColor.push_back(c2);

				trigVertex.push_back(v1);
				trigVertex.push_back(tt1);
				trigVertex.push_back(v2);
				trigNormal.push_back(-n2);
				trigNormal.push_back(-n1);
				trigNormal.push_back(-n3);
				trigColor.push_back(c2);
				trigColor.push_back(c2);
				trigColor.push_back(c2);
			}else if(angle2<=90){
				trigVertex.push_back(ss1);
				trigVertex.push_back(v1);
				trigVertex.push_back(tt1);
				trigNormal.push_back(-n0);
				trigNormal.push_back(-n2);
				trigNormal.push_back(-n1);
				trigColor.push_back(c2);
				trigColor.push_back(c2);
				trigColor.push_back(c2);

				trigVertex.push_back(v1);
				trigVertex.push_back(v2);
				trigVertex.push_back(tt1);

				trigNormal.push_back(-n2);
				trigNormal.push_back(-n3);
				trigNormal.push_back(-n1);

				trigColor.push_back(c2);
				trigColor.push_back(c2);
				trigColor.push_back(c2);
			}
			ss1=v1;
			tt1=v2;
		}
		s1=t1;
	}

	return true;
}
