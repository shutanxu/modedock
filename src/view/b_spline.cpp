/*
 * b_spline.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: stan
 */

#include"b_spline.h"

float
bt( const float& t ){
	float b = 0;
	if( t>=0 && t<=1 ){
		b = t * t * t / 6;
	}else if( t >=1 && t <=2 ){
		b = ( -3*(t-1)*(t-1)*(t-1) + 3*(t-1)*(t-1) + 3*(t-1)+1 )/6;
	}else if( t>=2 && t<=3 ){
		b = ( 3*(t-2)*(t-2)*(t-2) - 6*(t-2)*(t-2) + 4 )/6;
	}else if( t>=3 && t<=4 ){
		b = ( -(t-3)*(t-3)*(t-3)+3*(t-3)*(t-3)-3*(t-3)+1 )/6;
	}
	return b;
}

vector<QVector3D>
get_bSplineVertex(const vector<QVector3D>& P, const int& numOfVertex ){
	if( P.size() <= 2 ){
		cerr<<"error: must have at least two vertex for computing B-spline"<<endl;
	}
	vector<QVector3D>  allVertex;

	for( int i=0; i<numOfVertex; i++ ){
		float t = (float)(i*(P.size()-2))/(float)numOfVertex;
		QVector3D vertex(0, 0, 0);
		for( int j=0; j<P.size(); j++ ){
			vertex += (P[i]*bt(t - j));
			cout<<t<<" "<<bt(t - j)<<endl;
		}
		allVertex.push_back(vertex);
		if(1){
			cout<<endl;
			cout<<vertex.x()<<" "<<vertex.y()<<" "<<vertex.z()<<endl;
		}
	}
	return allVertex;
}

int
factorial(int n)
{
	int i=0,fact=1;
	if(n<=1)
	{
		return(1);
	}
	else
	{
		for(i=1;i<=n;i++)
		{
			fact=fact*i;
		}
		return(fact);
	}
}

float
B_factor( const int& i, const int&n, const float& t ){

	float b = ( factorial(n) * pow(t, i) * pow( (1-t), (n-i) ) ) /
			( factorial(i) * factorial(n-i) );
	return b;
}

vector<QVector3D>
get_BezierSpline( const vector<QVector3D>& P, const int& precision ){
//	cout<<"pow(0, 0):"<<pow(0, 0)<<endl;
	int n = (int)P.size()-1;
	vector<QVector3D>   allVec;
	int numOfVertex = (P.size()-1)*precision;
	for( int m=0; m<=numOfVertex; m++ ){
		float t = float(m)/(float)numOfVertex;
		QVector3D vec;
		for( int i=0; i<=n; i++ ){
			vec += P[i]*B_factor(i, n, t);
		}
		allVec.push_back( vec );
		if(0){
			cout<<vec.x()<<" "<<vec.y()<<" "<<vec.z()<<endl;
		}
	}
	return allVec;
}

vector<QVector3D>
get_interNormal( const vector<QVector3D>& P, const int& precision ){
	QVector3D averageNorm(0, 0, 0);
	for( size_t i=0; i<P.size(); i++ ){
		averageNorm += P[i];
	}
	averageNorm.normalize();
	cout<<"ave norm:"<<averageNorm.x()<<" "<<averageNorm.y()<<" "<<averageNorm.z()<<endl;
	vector<QVector3D>  allVec;
	int numOfOutputNormal = (P.size()-1)*precision;
	for( size_t i=0; i<numOfOutputNormal; i++ ){
		allVec.push_back( averageNorm );
	}
	return allVec;
}

/*
 *   p0   p1
 *   ------O
 */
vector<QVector3D>
get_circleVertex(const QVector3D p0, const QVector3D p1,
		const QVector3D p2){

	QVector3D v12 = p2 - p1;
	QVector3D v10 = p0 - p1;
	QVector3D n0 = QVector3D::crossProduct(v10, v12);
	n0 *= 10000;
	n0.normalize();
	QVector3D n2 = QVector3D::crossProduct(v10, n0);
	n2 *= 10000;
	n2.normalize();
	QVector3D n4 = QVector3D::crossProduct(v10, n2);
	n4 *= 10000;
	n4.normalize();
	QVector3D n6 = QVector3D::crossProduct(v10, n4);
	n6 *= 10000;
	n6.normalize();

	QVector3D n1 = (n0+n2)/2;
	n1 *= 10000;
	n1.normalize();
	QVector3D n3 = (n2+n4)/2;
	n3 *= 10000;
	n3.normalize();
	QVector3D n5 = (n4+n6)/2;
	n5 *= 10000;
	n5.normalize();
	QVector3D n7 = (n6+n0)/2;
	n7 *= 10000;
	n7.normalize();

	vector<QVector3D> vec;
	vec.push_back(n0);
	vec.push_back(n1);
	vec.push_back(n2);
	vec.push_back(n3);
	vec.push_back(n4);
	vec.push_back(n5);
	vec.push_back(n6);
	vec.push_back(n7);

	return vec;
}

vector<QVector3D>
get_interpolate( const vector<QVector3D>& P, const int& precision ){
	int n = (int)P.size()-1;
	vector<QVector3D>   allVec;
	allVec.clear();
	int numOfVertex = (P.size()-1)*precision;

	for( size_t i=1; i<P.size(); i++ ){
		for( size_t j=0; j<precision; j++ ){
			QVector3D pa = P[i-1];
			QVector3D pb = P[i];
			QVector3D d = (pb-pa)/precision;
			allVec.push_back( pa + j*d );
		}
	}
	return allVec;
}

void
get_turnTrigVertex( const vector<QVector3D>& boneVec ,
		const int& precision,
		vector<QVector3D>& trigVertexVec,
		vector<QVector3D>& trigVertexNormalVec){
	trigVertexVec.clear();
	trigVertexNormalVec.clear();
	vector<QVector3D>  boneVertexVec = get_BezierSpline( boneVec, precision );
//	vector<QVector3D>  boneVertexVec = get_interpolate( boneVec, precision );

	vector<QVector3D> v0 = get_circleVertex( boneVertexVec[0], boneVertexVec[1], boneVertexVec[2] );
	for( size_t i=1; i<(boneVertexVec.size()-1); i++ ){

		vector<QVector3D> v1 = get_circleVertex( boneVertexVec[i-1], boneVertexVec[i], boneVertexVec[i+1] );
		for( size_t j=0; j<v1.size(); j++ ){
			trigVertexVec.push_back(v0[j]*0.001 + boneVertexVec[i-1]  );
			trigVertexVec.push_back(v0[(j+1)%v1.size()]*0.001 + boneVertexVec[i-1] );
			trigVertexVec.push_back(v1[j]*0.001 + boneVertexVec[i] );
			trigVertexNormalVec.push_back( v0[j] );
//			trigVertexNormalVec.push_back( v0[(j+1)%v1.size()] );
			trigVertexNormalVec.push_back( v0[j] );
			trigVertexNormalVec.push_back( v1[j] );

			trigVertexVec.push_back(v1[j]*0.001 + boneVertexVec[i]  );
			trigVertexVec.push_back(v1[(j+1)%v1.size()]*0.001 + boneVertexVec[i] );
			trigVertexVec.push_back(v0[(j+1)%v1.size()]*0.001 + boneVertexVec[i-1] );
			trigVertexNormalVec.push_back( v1[j] );
			trigVertexNormalVec.push_back( v1[j] );
//			trigVertexNormalVec.push_back( v0[(j+1)%v0.size()] );
			trigVertexNormalVec.push_back( v0[j] );
		}
		v0 = v1;
	}
}

/*
 *   2----6----10
 *   :\   :\   :\
 *   : \  : \  : \
 *   0--\-4--\-8  \
 *    \  \    \    \
 *     \  3----7----11
 *      \ :    :    :
 *       \:    :    :
 *        1----5----9
 */
void
get_triangleVertex( const vector<QVector3D>& boneVec ,
		const vector<QVector3D>& normVec,
		const int& precision,
		vector<QVector3D>& trigVertexVec,
		vector<QVector3D>& trigVertexNormalVec){

	vector<QVector3D> vertexVec, vertexNormalVec;
	vertexVec.clear();
	vertexNormalVec.clear();

	int range = precision * ( boneVec.size() - 1 );

	if( boneVec.size() != normVec.size() ){
		cerr<<"error !!! please input the same number of vertex and its normal "<<endl;
	}

	vector<QVector3D>  boneVertex = get_BezierSpline( boneVec, precision );
	vector<QVector3D>  boneNormal = get_interNormal( normVec, precision );
	vector<QVector3D>  crossNormal;

	for( size_t i=0; i<(precision * ( boneVec.size() - 2 )) ; i++ ){
		QVector3D bv;
		if( i != range-1 ){
			bv = boneVertex[i+1]-boneVertex[i];
		}else{
			bv = boneVertex[i]-boneVertex[i-1];
		}
		QVector3D hv = QVector3D::crossProduct(bv, boneNormal[i]);
		crossNormal.push_back(hv);
		hv.normalize();
		QVector3D a0 = boneVertex[i] - boneNormal[i]*0.01 + hv * 0.002;
		QVector3D a1 = boneVertex[i] + boneNormal[i]*0.01 + hv * 0.002;
		QVector3D a2 = boneVertex[i] - boneNormal[i]*0.01 - hv * 0.002;
		QVector3D a3 = boneVertex[i] + boneNormal[i]*0.01 - hv * 0.002;
		vertexVec.push_back( a0 );
		vertexVec.push_back( a1 );
		vertexVec.push_back( a2 );
		vertexVec.push_back( a3 );
	}
	// vertex in arrowhead
	for( size_t i=(precision * ( boneVec.size() - 2 ));
			i<(precision * ( boneVec.size() - 1 )) ; i++ ){
		QVector3D bv;
		if( i != range-1 ){
			bv = boneVertex[i+1]-boneVertex[i];
		}else{
			bv = boneVertex[i]-boneVertex[i-1];
		}
		QVector3D hv = QVector3D::crossProduct(bv, boneNormal[i]);
		crossNormal.push_back(hv);
		hv.normalize();
		if( i==(precision * ( boneVec.size() - 1 )) ){
			QVector3D a0 = boneVertex[i] - boneNormal[i]*0.01 + hv * 0.002;
			QVector3D a1 = boneVertex[i] + boneNormal[i]*0.01 + hv * 0.002;
			QVector3D a2 = boneVertex[i] - boneNormal[i]*0.01 - hv * 0.002;
			QVector3D a3 = boneVertex[i] + boneNormal[i]*0.01 - hv * 0.002;
			vertexVec.push_back( a0 );
 			vertexVec.push_back( a1 );
			vertexVec.push_back( a2 );
			vertexVec.push_back( a3 );
		}else{
			size_t s = i-(precision * ( boneVec.size() - 2 ));
			float scale = 2 * (float)(precision - s)/(float)precision;
			QVector3D a0 = boneVertex[i] - boneNormal[i]*0.01*scale + hv * 0.002;
			QVector3D a1 = boneVertex[i] + boneNormal[i]*0.01*scale + hv * 0.002;
			QVector3D a2 = boneVertex[i] - boneNormal[i]*0.01*scale - hv * 0.002;
			QVector3D a3 = boneVertex[i] + boneNormal[i]*0.01*scale - hv * 0.002;
			vertexVec.push_back( a0 );
			vertexVec.push_back( a1 );
			vertexVec.push_back( a2 );
			vertexVec.push_back( a3 );
		}
	}

	// surface 0123
	QVector3D n = boneVertex[0] - boneVertex[1];
	n.normalize();
	trigVertexVec.push_back( vertexVec[0] );
	trigVertexVec.push_back( vertexVec[1] );
	trigVertexVec.push_back( vertexVec[2] );
	trigVertexVec.push_back( vertexVec[1] );
	trigVertexVec.push_back( vertexVec[2] );
	trigVertexVec.push_back( vertexVec[3] );
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);

	if(0){
		cout<<"left"<<endl;
		cout<<"0:"<<vertexVec[0].x()<<" "<<vertexVec[0].y()<<" "<<vertexVec[0].z()<<endl;
		cout<<"2:"<<vertexVec[2].x()<<" "<<vertexVec[2].y()<<" "<<vertexVec[2].z()<<endl;
		cout<<"1:"<<vertexVec[1].x()<<" "<<vertexVec[1].y()<<" "<<vertexVec[1].z()<<endl;
		cout<<"3:"<<vertexVec[3].x()<<" "<<vertexVec[3].y()<<" "<<vertexVec[3].z()<<endl;
		cout<<"n:"<<n.x()<<" "<<n.y()<<" "<<n.z()<<endl;
	}

	// surface right 8-9-10-11
	n = boneVertex[range-1] - boneVertex[range-2];
	n.normalize();
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 0] );
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 1] );
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 2] );
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 1] );
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 2] );
	trigVertexVec.push_back( vertexVec[ (range-1)*4 + 3] );
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);
	trigVertexNormalVec.push_back(n);

	if(0){
		cout<<"right"<<endl;
		cout<<"0:"<<vertexVec[ (range-1)*4 + 0].x()<<" "<<vertexVec[ (range-1)*4 + 0].y()<<" "<<vertexVec[ (range-1)*4 + 0].z()<<endl;
		cout<<"2:"<<vertexVec[ (range-1)*4 + 2].x()<<" "<<vertexVec[ (range-1)*4 + 2].y()<<" "<<vertexVec[ (range-1)*4 + 2].z()<<endl;
		cout<<"1:"<<vertexVec[ (range-1)*4 + 1].x()<<" "<<vertexVec[ (range-1)*4 + 1].y()<<" "<<vertexVec[ (range-1)*4 + 1].z()<<endl;
		cout<<"3:"<<vertexVec[ (range-1)*4 + 3].x()<<" "<<vertexVec[ (range-1)*4 + 3].y()<<" "<<vertexVec[ (range-1)*4 + 3].z()<<endl;
		cout<<"n:"<<n.x()<<" "<<n.y()<<" "<<n.z()<<endl;
	}

	//
	for( size_t i=0; i<range-1; i++ ){
		//top 2-3-6-7-10-11
		trigVertexVec.push_back( vertexVec[ i*4 + 2] );
		trigVertexNormalVec.push_back( -crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 3] );
		trigVertexNormalVec.push_back( -crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 6] );
		trigVertexNormalVec.push_back( -crossNormal[i+1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 3] );
		trigVertexNormalVec.push_back( -crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 6] );
		trigVertexNormalVec.push_back( -crossNormal[i+1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 7] );
		trigVertexNormalVec.push_back( -crossNormal[i+1] );

		if(0){
			cout<<"top"<<endl;
			cout<<"2:"<<vertexVec[ i*4 + 2].x()<<" "<<vertexVec[ i*4 + 2].y()<<" "<<vertexVec[ i*4 + 2].z()<<endl;
			cout<<"6:"<<vertexVec[ i*4 + 6].x()<<" "<<vertexVec[ i*4 + 6].y()<<" "<<vertexVec[ i*4 + 6].z()<<endl;
			cout<<"3:"<<vertexVec[ i*4 + 3].x()<<" "<<vertexVec[ i*4 + 3].y()<<" "<<vertexVec[ i*4 + 3].z()<<endl;
			cout<<"7:"<<vertexVec[ i*4 + 7].x()<<" "<<vertexVec[ i*4 + 7].y()<<" "<<vertexVec[ i*4 + 7].z()<<endl;
			cout<<"n:"<<n.x()<<" "<<n.y()<<" "<<n.z()<<endl;
		}

		//bottom 0-1-4-5-8-9
		trigVertexVec.push_back( vertexVec[ i*4 + 0] );
		trigVertexNormalVec.push_back( crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 1] );
		trigVertexNormalVec.push_back( crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 4] );
		trigVertexNormalVec.push_back( crossNormal[i+1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 1] );
		trigVertexNormalVec.push_back( crossNormal[i] );
		trigVertexVec.push_back( vertexVec[ i*4 + 4] );
		trigVertexNormalVec.push_back( crossNormal[i+1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 5] );
		trigVertexNormalVec.push_back( crossNormal[i+1] );

		if(0){
			cout<<"top"<<endl;
			cout<<"2:"<<vertexVec[ i*4 + 2].x()<<" "<<vertexVec[ i*4 + 2].y()<<" "<<vertexVec[ i*4 + 2].z()<<endl;
			cout<<"6:"<<vertexVec[ i*4 + 6].x()<<" "<<vertexVec[ i*4 + 6].y()<<" "<<vertexVec[ i*4 + 6].z()<<endl;
			cout<<"3:"<<vertexVec[ i*4 + 3].x()<<" "<<vertexVec[ i*4 + 3].y()<<" "<<vertexVec[ i*4 + 3].z()<<endl;
			cout<<"7:"<<vertexVec[ i*4 + 7].x()<<" "<<vertexVec[ i*4 + 7].y()<<" "<<vertexVec[ i*4 + 7].z()<<endl;
			cout<<"n:"<<n.x()<<" "<<n.y()<<" "<<n.z()<<endl;
		}

		//near 1-3-5-7-9-11
		trigVertexVec.push_back( vertexVec[ i*4 + 1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 3] );
		trigVertexVec.push_back( vertexVec[ i*4 + 7] );
		trigVertexVec.push_back( vertexVec[ i*4 + 1] );
		trigVertexVec.push_back( vertexVec[ i*4 + 7] );
		trigVertexVec.push_back( vertexVec[ i*4 + 5] );

		if( i==(precision * ( boneVec.size() - 2 ) - 1 ) ){
			// normal of arrow head bottom
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
		}else{
			trigVertexNormalVec.push_back( boneNormal[i] );
			trigVertexNormalVec.push_back( boneNormal[i] );
			trigVertexNormalVec.push_back( boneNormal[i+1] );
			trigVertexNormalVec.push_back( boneNormal[i] );
			trigVertexNormalVec.push_back( boneNormal[i+1] );
			trigVertexNormalVec.push_back( boneNormal[i+1] );
		}

		//far 0-2-6-4-8-10
		trigVertexVec.push_back( vertexVec[ i*4 + 0] );
		trigVertexVec.push_back( vertexVec[ i*4 + 2] );
		trigVertexVec.push_back( vertexVec[ i*4 + 6] );
		trigVertexVec.push_back( vertexVec[ i*4 + 0] );
		trigVertexVec.push_back( vertexVec[ i*4 + 6] );
		trigVertexVec.push_back( vertexVec[ i*4 + 4] );

		if( i==(precision * ( boneVec.size() - 2 ) - 1) ){
			// normal of arrow head bottom
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
			trigVertexNormalVec.push_back( boneVertex[i-1]-boneVertex[i]  );
		}else{
			trigVertexNormalVec.push_back( -boneNormal[i] );
			trigVertexNormalVec.push_back( -boneNormal[i] );
			trigVertexNormalVec.push_back( -boneNormal[i+1] );
			trigVertexNormalVec.push_back( -boneNormal[i] );
			trigVertexNormalVec.push_back( -boneNormal[i+1] );
			trigVertexNormalVec.push_back( -boneNormal[i+1] );
		}
	}
}
