/*
 * polyeder.cpp
 *
 *  Created on: Sep 3, 2014
 *      Author: stan
 */

#include"polyeder.h"
#include "../basicMath/Math.h"

Polyeder::Polyeder(const size_t ord) :order(ord)
{  /* GENERATES ALL 12 VERTICES OF ICOSAHEDRON */
	double v[ dim3 * 12];  // 36 elements
	double v1[dim3];
	double v2[dim3];
	double v3[dim3];
	size_t i, j, k, s, level;

	double a = YVERTEX;
	double b = ZVERTEX;
	k = 0;
	for ( i = 0; i < 2; i++) {
		a = -a;
		for (j = 0; j < 2; j++, k++) {
			b = -b;
			v[dim3 * k] = 0.0;
			v[dim3 * k+1] = a;
			v[dim3 * k+2] = b;
			k++;
			v[dim3 *k] = b;
			v[dim3 *k + 1] = 0.0;
			v[dim3 *k + 2] = a;
			k++;
			v[ dim3 * k ] = a;
			v[ dim3 * k + 1] = b;
			v[ dim3 * k + 2] = 0.0;
		}
	}

	if(0){
		for( size_t i=0; i<12; i++ ){
			cout.width(4);
			cout<<left<<i;
			cout.width(10);
			cout<<left<<v[i*3];
			cout.width(10);
			cout<<v[i*3+1];
			cout.width(10);
			cout<<v[i*3+2]<<endl;
		}
	}

	level = order;
	/* GET ALL 20 FACES OF ICOSAHEDRON */
	//0
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 0 + s];
		v2[s] = v[dim3 * 2 + s];
		v3[s] = v[dim3 * 1 + s];
	}
	Triangle(v1, v2, v3, level);

	//1
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 0 + s];
		v2[s] = v[dim3 * 1 + s];
		v3[s] = v[dim3 * 4 + s];
	}
	Triangle(v1, v2, v3, level);

	//2
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 0 + s];
		v2[s] = v[dim3 * 3 + s];
		v3[s] = v[dim3 * 2 + s];
	}
	Triangle(v1, v2, v3, level);

	//3
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 0 + s];
		v2[s] = v[dim3 * 8 + s];
		v3[s] = v[dim3 * 3 + s];
	}
	Triangle(v1, v2, v3, level);

	//4
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 0 + s];
		v2[s] = v[dim3 * 4 + s];
		v3[s] = v[dim3 * 8 + s];
	}
	Triangle(v1, v2, v3, level);

	//5
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 1 + s];
		v2[s] = v[dim3 * 2 + s];
		v3[s] = v[dim3 * 5 + s];
	}
	Triangle(v1, v2, v3, level);

	//6
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 1 + s];
		v2[s] = v[dim3 * 6 + s];
		v3[s] = v[dim3 * 4 + s];
	}
	Triangle(v1, v2, v3, level);

	//7
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 1 + s];
		v2[s] = v[dim3 * 5 + s];
		v3[s] = v[dim3 * 6 + s];
	}
	Triangle(v1, v2, v3, level);

	//8
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 2 + s];
		v2[s] = v[dim3 * 3 + s];
		v3[s] = v[dim3 * 7 + s];
	}
	Triangle(v1, v2, v3, level);

	//9
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 2 + s];
		v2[s] = v[dim3 * 7 + s];
		v3[s] = v[dim3 * 5 + s];
	}
	Triangle(v1, v2, v3, level);

	//10
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 3 + s];
		v2[s] = v[dim3 * 10 + s];
		v3[s] = v[dim3 * 7 + s];
	}
	Triangle(v1, v2, v3, level);

	//11
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 3 + s];
		v2[s] = v[dim3 * 8 + s];
		v3[s] = v[dim3 * 10 + s];
	}
	Triangle(v1, v2, v3, level);

	//12
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 4 + s];
		v2[s] = v[dim3 * 6 + s];
		v3[s] = v[dim3 * 11 + s];
	}
	Triangle(v1, v2, v3, level);

	//13
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 4 + s];
		v2[s] = v[dim3 * 11 + s];
		v3[s] = v[dim3 * 8 + s];
	}
	Triangle(v1, v2, v3, level);

	//14
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 5 + s];
		v2[s] = v[dim3 * 9 + s];
		v3[s] = v[dim3 * 6 + s];
	}
	Triangle(v1, v2, v3, level);

	//15
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 5 + s];
		v2[s] = v[dim3 * 7 + s];
		v3[s] = v[dim3 * 9 + s];
	}
	Triangle(v1, v2, v3, level);

	//16
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 6 + s];
		v2[s] = v[dim3 * 9 + s];
		v3[s] = v[dim3 * 11 + s];
	}
	Triangle(v1, v2, v3, level);

	//17
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 7 + s];
		v2[s] = v[dim3 * 10 + s];
		v3[s] = v[dim3 * 9 + s];
	}
	Triangle(v1, v2, v3, level);

	//18
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 8 + s];
		v2[s] = v[dim3 * 11 + s];
		v3[s] = v[dim3 * 10 + s];
	}
	Triangle(v1, v2, v3, level);

	//19
	for( s =0; s<dim3; s++){
		v1[s] = v[dim3 * 9 + s];
		v2[s] = v[dim3 * 10 + s];
		v3[s] = v[dim3 * 11 + s];
	}
	Triangle(v1, v2, v3, level);

	//-------------------------------------------------------
	a = 0.0;
	for (i = 0; i < wp.size(); i++)
		a += wp[i];
	a = FOURPI / a;
	for (i = 0; i < wp.size(); i++)
		wp[i] *= a;
	accept = vector<bool>(wp.size(), false);
}


void Polyeder::Triangle(double *const x1, double *const x2, double *const x3,
		const size_t level ) {
	unsigned int  k, level1;
	double xnorm = 1.0;
	double x4[dim3], x5[dim3], x6[dim3];

	if (level > 0) {
		level1 = level - 1;
		for (k = 0; k < dim3 ; k++) {
			x4[k] = x1[k] + x2[k];
			x5[k] = x2[k] + x3[k];
			x6[k] = x1[k] + x3[k];
		}
		Norm(x4, &xnorm);
		Norm(x5, &xnorm);
		Norm(x6, &xnorm);
		Triangle(x1, x4, x6, level1);
		Triangle(x4, x5, x6, level1);
		Triangle(x4, x2, x5, level1);
		Triangle(x5, x3, x6, level1);
		return;
	}
	for (k = 0; k <dim3; k++)
		x6[k] = x1[k] + x2[k] + x3[k];
	Norm(x6, &xnorm);
	for (k=0; k < dim3 ;k++)
		p.push_back(x6[k]);
	Diff(x3, x1, x5);
	Diff(x2, x1, x4);
	Cross(x5, x4, x6);
	Norm(x6, &xnorm);
	wp.push_back(0.50 * xnorm );

	double v12[3];
	double v13[3];
	double c[3];
	for(int i=0; i<3; i++){
		v12[i] = (x2[i] - x1[i]);
		v13[i] = (x3[i] - x1[i]);
		c[i] = (x1[i] + x2[i] + x3[i])/3.0;
	}
	double cro[3];
	Cross(v12, v13, cro);
	Norm(cro, &xnorm);
	double d = dot( cro, c );

	TriangleSurf triangle, triangleBack;
	for( size_t i=0; i<3; i++ ){
		triangle.v1[i] = x1[i];
		triangle.v2[i] = x2[i];
		triangle.v3[i] = x3[i];
		triangle.normal[i] = cro[i];
	}

	allTriangles.push_back(triangle);
}

double*
Polyeder::cross( double *p1, double *p2 ){
	double c[3];
	c[0] = c[1] = c[2] = 0;
	c[0] = p1[1] * p2[2] - p2[1] * p1[2];
	c[1] = p2[0] * p1[2] - p2[2] * p1[0];
	c[2] = p1[0] * p2[1] - p2[0] * p1[1];
	cout<<"croo:"<<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
	return c;
}

double
Polyeder::dot( double *p1, double *p2 ){
	double v = 0;
	for(int i=0; i<3; i++){
		v += p1[i] * p2[i];
	}
	return v;
}


double
Polyeder::distance( double* a, double* b ){
	double dis = (a[0] - b[0])*(a[0] - b[0]) +
			(a[1] - b[1])*(a[1] - b[1]) +
			(a[2] - b[2])*(a[2] - b[2]);
	return sqrt( dis );
}

void Polyeder::polyederReset() {
	for(unsigned int i=0;i< wp.size();i++)
		accept[i]= false;
}

vector<float>
Polyeder::get_allVertex(){
	vector<vector<double> >  vertex;
	for( size_t i=0; i<allTriangles.size(); i++ ){
		double *v1 = allTriangles[i].v1 ;
		double *v2 = allTriangles[i].v2 ;
		double *v3 = allTriangles[i].v3 ;

		//-------------------------------------------
		bool   exist ;
		size_t index;

		exist = false;
		index = 0;
		for( size_t j=0; j<vertex.size(); j++ ){
            if(     vertex[j][0] == v1[0] &&
            		vertex[j][1] == v1[1] &&
            		vertex[j][2] == v1[2]){
            	exist = true;
            	index = j;
            	break;
            }
		}

		if( exist ){
			triangleVertexIndex.push_back( (int)index );
		}else{
			vector<double> v(3);
			v[0] = v1[0];
			v[1] = v1[1];
			v[2] = v1[2];
			vertex.push_back(v);
			triangleVertexIndex.push_back( (int)vertex.size()-1 );
		}

		//-------------------------------------------
		exist = false;
		index = 0;
		for( size_t j=0; j<vertex.size(); j++ ){
            if(     vertex[j][0] == v2[0] &&
            		vertex[j][1] == v2[1] &&
            		vertex[j][2] == v2[2]){
            	exist = true;
            	index = j;
            	break;
            }
		}

		if( exist ){
			triangleVertexIndex.push_back( (int)index);
		}else{
			vector<double> v(3);
			v[0] = v2[0];
			v[1] = v2[1];
			v[2] = v2[2];
			vertex.push_back(v);
			triangleVertexIndex.push_back( (int)vertex.size()-1 );
		}

		//-------------------------------------------
		exist = false;
		index = 0;
		for( size_t j=0; j<vertex.size(); j++ ){
            if(     vertex[j][0] == v3[0] &&
            		vertex[j][1] == v3[1] &&
            		vertex[j][2] == v3[2]){
            	exist = true;
            	index = j;
            	break;
            }
		}

		if( exist ){
			triangleVertexIndex.push_back((int)index);
		}else{
			vector<double> v(3);
			v[0] = v3[0];
			v[1] = v3[1];
			v[2] = v3[2];
			vertex.push_back(v);
			triangleVertexIndex.push_back( (int)vertex.size()-1 );
		}
	}

	vector<float> allVertex;
	for( size_t i=0; i<vertex.size(); i++ ){
		allVertex.push_back( (float)vertex[i][0] );
		allVertex.push_back( (float)vertex[i][1] );
		allVertex.push_back( (float)vertex[i][2] );
	}

	return allVertex;
}
