/*
 * stick.h
 *
 *  Created on: Dec 2, 2014
 *      Author: stan
 */

#ifndef STICK_H_
#define STICK_H_

#include<QGLWidget>
#include<QtOpenGL/QGLFunctions>
#include<QtOpenGL/QGLShaderProgram>
#include<QVector3D>
#include<GL/gl.h>
#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

QVector3D
rotateAroundAxis( const QVector3D input,
		const QVector3D axis,
		const double angle );

bool
getStick(const QVector3D p1, const QVector3D p2,
		const QColor c1, const QColor c2,
		const float ratio,const float radius,
		vector<QVector3D>& trigVertex,
		vector<QVector3D>& trigNormal,
		vector<QColor>& trigColor );


#endif /* STICK_H_ */
