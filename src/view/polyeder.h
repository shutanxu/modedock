/*
 * polyeder.h
 *
 *  Created on: Sep 3, 2014
 *      Author: stan
 */

#ifndef POLYEDER_H_
#define POLYEDER_H_

#include <QtOpenGL/QGLFunctions>
#include <QtOpenGL/QGLShaderProgram>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>

#include<iostream>
#include<vector>
#include<algorithm>
#include<math.h>

using namespace std;
using std::vector;

struct TriangleSurf{
	double v1[3];
	double v2[3];
	double v3[3];
	double normal[3];
};

class Polyeder{
public:
//	constexpr static size_t dim3 = 3;
//	constexpr static double  FOURPI  = 12.56637;
//	constexpr static double  YVERTEX = 0.8506508;
//	constexpr static double  ZVERTEX  = 0.5257311;
//	constexpr static double  EPSlocal = 0.00001;
	size_t dim3 = 3;
	double  FOURPI  = 12.56637;
	double  YVERTEX = 0.8506508;
	double  ZVERTEX  = 0.5257311;
	double  EPSlocal = 0.00001;

	Polyeder(){;}
	Polyeder(const size_t order);
	~Polyeder(){;}

	double* cross( double *v1, double *v2 );
	double dot( double *v1, double *v2 );

	void Triangle( double *const x1, double *const x2,
			double *const x3, const size_t level );
	void polyederReset();
	vector<TriangleSurf> get_surfTriangles(){ return allTriangles;}

	vector<float>    get_allVertex();
	vector<int>      get_triangleVertexIndex(){ return triangleVertexIndex; };
private:
	double distance( double* a, double* b );
	size_t order;
	vector<double> p;
	vector<double> wp;
	vector<bool> accept;
	size_t iter;
	vector<TriangleSurf> allTriangles;
	vector<int> triangleVertexIndex;
};

#endif /* POLYEDER_H_ */
