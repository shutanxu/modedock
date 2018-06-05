/*
 * b_spline.h
 *
 *  Created on: Sep 22, 2014
 *      Author: stan
 */

#ifndef B_SPLINE_H_
#define B_SPLINE_H_

#include<iostream>
#include<vector>
#include<algorithm>
#include<stdlib.h>
#include<math.h>
#include<QVector3D>

using namespace std;

float
bt( const float& t );

vector<QVector3D>
get_bSplineVertex(const vector<QVector3D>& P, const int& numOfVertex );


int
factorial(int n);

float
B_factor( const int& i, const int&n, const float& t );

vector<QVector3D>
get_BezierSpline( const vector<QVector3D>& P, const int& precision );

vector<QVector3D>
get_interNormal( const vector<QVector3D>& P, const int& numOfOutputNormal );

vector<QVector3D>
get_circleVertex(const QVector3D p0, const QVector3D p1,
		const QVector3D p2);

void
get_turnTrigVertex( const vector<QVector3D>& boneVec ,
		const int& precision,
		vector<QVector3D>& trigVertexVec,
		vector<QVector3D>& trigVertexNormalVec);

void
get_triangleVertex( const vector<QVector3D>& boneVec ,
		const vector<QVector3D>& normVec,
		const int& precision,
		vector<QVector3D>& trigVertexVec,
		vector<QVector3D>& trigVertexNormalVec);

#endif /* B_SPLINE_H_ */
