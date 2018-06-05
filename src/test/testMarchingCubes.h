/*
 * testMarchingCubes.h
 *
 *  Created on: Jul 31, 2014
 *      Author: stan
 */

#ifndef TESTMARCHINGCUBES_H_
#define TESTMARCHINGCUBES_H_

#include"../mol/molecular.h"
#include"../marchingCube/MarchingCubes.h"
//#include"../mol/vdwSurface.h"

using namespace std;
/*
// Compute data
void compute_ball( MarchingCubes &mc )
//-----------------------------------------------------------------------------
{
  float x,y,z      ;
  float sx,sy,sz   ;
  float tx,ty,tz   ;

  float r,R ;
  r = 1.85f ;
  R = 4 ;

  sx     = (float) mc.size_x() / 16 ; // 4
  sy     = (float) mc.size_y() / 16 ; // 4
  sz     = (float) mc.size_z() / 16 ; // 4
  tx     = (float) mc.size_x() / (2*sx) ; // 8
  ty     = (float) mc.size_y() / (2*sy) ; // 8
  tz     = (float) mc.size_z() / (2*sz) ; // 8

  for( int k = 0 ; k < mc.size_z() ; k++ )
  {
    z = ( (float) k ) / sz  - tz ;

    for( int j = 0 ; j < mc.size_y() ; j++ )
    {
      y = ( (float) j ) / sy  - ty ;

      for( int i = 0 ; i < mc.size_x() ; i++ )
      {
        x = ( (float) i ) / sx - tx ;
        mc.set_data( (float) (
          // cushin
//                      z*z*x*x - z*z*z*z - 2*z*x*x + 2*z*z*z + x*x - z*z - (x*x - z)*(x*x - z) - y*y*y*y - 2*x*x*y*y - y*y*z*z + 2*y*y*z + y*y ) ,
          // sphere
//                      ( (x-2)*(x-2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x+2)*(x+2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x-2)*(x-2) + (y+2)*(y+2) + (z-2)*(z-2) - 1 )) ,
          //  plane
//                      x+y+z -3) ,
          // cassini
//                      (x*x + y*y + z*z + 0.45f*0.45f)*(x*x + y*y + z*z + 0.45f*0.45f) - 16*0.45f*0.45f*(x*x + z*z) - 0.5f*0.5f ),
          // blooby
          //           x*x*x*x - 5*x*x+ y*y*y*y - 5*y*y + z*z*z*z - 5*z*z + 11.8 ),
          //  chair
          //            x*x+y*y+z*z-0.95f*25)*(x*x+y*y+z*z-0.95f*25)-0.8f*((z-5)*(z-5)-2*x*x)*((z+5)*(z+5)-2*y*y ,
          // cyclide
          //            ( x*x + y*y + z*z + b*b - d*d ) * ( x*x + y*y + z*z + b*b - d*d ) - 4 * ( ( a*x - c*d ) * ( a*x - c*d ) + b*b * y*y ),
          // 2 torus
//                      ( ( x*x + y*y + z*z + R*R - r*r ) * ( x*x + y*y + z*z + R*R - r*r ) - 4 * R*R * ( x*x + y*y ) ) *
//                      ( ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) * ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) - 4 * R*R * ( (y+R)*(y+R) + z*z ) ) ) ,
          // mc case
          // - 26.5298*(1-x)*(1-y)*(1-z) + 81.9199*x*(1-y)*(1-z) - 100.68*x*y*(1-z) + 3.5498*(1-x)*y*(1-z)
          // + 24.1201*(1-x)*(1-y)*  z   - 74.4702*x*(1-y)*  z   + 91.5298*x*y*  z  - 3.22998*(1-x)*y*  z  ),
          // Drip
                    x*x + y*y - 0.5*( 0.995*z*z + 0.005 - z*z*z ) +0.0025 ),  // -0.0754+0.01, -0.0025 + 0.01, grid 40^3, [-1.5,1.5]

          i,j,k ) ;
      }
    }
  }
}

int testMC2 (int argc, char **argv){
	QApplication application(argc,argv);
	string							fileNameA ( "/mnt/2t/xushutan/project/gold_dock/test1024/normal/1com/ligand_corina.mol2" ) ;
	Molecular						molA( fileNameA );
	molA = readMolecularFile( fileNameA ).front();

	MarchingCubes mc ;
	//	mc.set_resolution( 60,60,60 ) ;

	//	mc.init_all() ;
	//	compute_VdwMarchingCube( molA , mc ) ;
	mc.set_resolution( 60, 60, 60 );
	compute_ball( mc );
	mc.run() ;
	mc.clean_temps() ;

	cout<<"trig:"<<mc.ntrigs()<<endl;
	cout<<"ver:"<<mc.nverts()<<endl;
	cout<<"-----------------"<<endl;

	vector<TriangleSurf>			trigVec;
	float		maxX, maxY, maxZ;
	float		minX, minY, minZ;
	float		aveX, aveY, aveZ;
	maxX = maxY = maxZ = -10e10;
	minX = minY = minZ = 10e10;
	aveX = aveY = aveZ = 0;
	if(1){
		cout<<"trig:"<<mc.ntrigs()<<endl;
		cout<<"ver:"<<mc.nverts()<<endl;
		for( int i=0; i<mc.ntrigs(); i++ ){
			//			cout<<"0:"<<i<<endl;
			//			cout<<"1:"<<i<<endl;
			const Triangle	*trig = mc.trig(i);
			const Vertex *v1 = mc.vert( trig->v1 );
			const Vertex *v2 = mc.vert( trig->v2 );
			const Vertex *v3 = mc.vert( trig->v3 );

			Coord	c1( v1->x, v1->y, v1->z );
			Coord	c2( v2->x, v2->y, v2->z );
			Coord	c3( v3->x, v3->y, v3->z );

			Coord	n1( v1->nx, v1->ny, v1->nz );
			Coord	n2( v2->nx, v2->ny, v2->nz );
			Coord	n3( v3->nx, v3->ny, v3->nz );

			maxX = maxX > v1->x ? maxX:v1->x ;
			maxY = maxY > v1->y ? maxY:v1->y ;
			maxZ = maxZ > v1->z ? maxZ:v1->z ;

			minX = minX < v1->x ? minX:v1->x ;
			minY = minY < v1->y ? minY:v1->y ;
			minZ = minZ < v1->z ? minZ:v1->z ;

			aveX += v1->x ;
			aveY += v1->y ;
			aveZ += v1->z ;

//			TriangleSurf	triangle( c1, n1, c2, n2, c3, n3 );
//			trigVec.push_back( triangle );
		}
		aveX = aveX /(float)mc.ntrigs();
		aveY = aveY /(float)mc.ntrigs();
		aveZ = aveZ /(float)mc.ntrigs();
		cout<<"ave: "<<aveX<<" "<<aveY<<" "<<aveZ<<endl;
		cout<<"max: "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
		cout<<"min: "<<minX<<" "<<minY<<" "<<minZ<<endl;
	}

	mc.clean_all() ;
	return application.exec()  ;
	return 0;
}
*/
/*

struct TriangleSurf{
	TriangleSurf( const Coord& va, const Coord& na, const Coord& vb, const Coord& nb, const Coord& vc, const Coord& nc ):v1(va), n1(na), v2(vb), n2(nb), v3(vc), n3( nc ){}
	Coord		v1;
	Coord		n1;
	Coord		v2;
	Coord		n2;
	Coord		v3;
	Coord		n3;
};

class plotMC : public QGLViewer {
public:
	plotMC(){}
	plotMC( const vector<TriangleSurf>& tri ):surf(tri){}
	plotMC( MarchingCubes &m ){ mc = m ; }
    plotMC( const vector<pair<float, float> >& data ):inputCoordVec( data ){};
    plotMC( const vector<pair<float, float> >& data,
    		const string& xlab, const string& ylab):inputCoordVec( data ), xlabel(xlab), ylabel(ylab){};

protected :
  virtual void draw();
  virtual void init();
  virtual QString helpString() const;

private:
  void 			dataInit();
  void 			drawData();
  void 			drawAxis();
  void			testDraw();

  vector< pair<float, float> > 	inputCoordVec;
  vector< pair<float, float> > 	viewerCoordVec;
  string											xlabel;
  string											ylabel;
//  float											minX;
//  float											maxX;
//  float											minY;
//  float											maxY;
  GLuint 										dataList;			//    GLuint modelList;
  float											marginRatio;
  vector<	qglviewer::Vec> 			proj;
  MarchingCubes 							mc ;
  vector<TriangleSurf>				surf;

	float	aveX ;
	float	aveY ;
	float	aveZ ;
	float	maxX;
	float	maxY ;
	float	maxZ;
	float	minX ;
	float	minY ;
	float	minZ ;
	float	rangeX;
	float	rangeY ;
	float	rangeZ ;
};


void
plotMC::dataInit(){

	if(inputCoordVec.empty()){
		  return;
	  }

	  minX = minY = -1e10;
	  maxX = maxY = 1e10;
	  for( size_t i=0; i<inputCoordVec.size(); i++ ){
		  minX = minX < inputCoordVec[i].first ? minX:inputCoordVec[i].first;
		  maxX = maxX > inputCoordVec[i].first ? maxX:inputCoordVec[i].first;
		  minY = minY < inputCoordVec[i].second ? minY:inputCoordVec[i].second;
		  maxY = maxY > inputCoordVec[i].second ? minY:inputCoordVec[i].second;
	  }

	  marginRatio = 0.1;
	//  float				xval = *std::max_element( inputCoordVec.begin(), inputCoordVec.end(), cmpPairFirst<float, float>() );
//	  std::sort(inputCoordVec.begin(), inputCoordVec.end(), cmpPairFirst<float, float>());
	  float				xRange = inputCoordVec.back().first - inputCoordVec.front().first;
	  float				xMiddle = (  inputCoordVec.back().first + inputCoordVec.front().first ) /2;
//	  std::sort(inputCoordVec.begin(), inputCoordVec.end(), cmpPairSecond<float, float>() );
	  float				yRange =inputCoordVec.back().second - inputCoordVec.front().second;
	  float				yMiddle = (  inputCoordVec.back().second + inputCoordVec.front().second ) /2;

	  viewerCoordVec.clear();
	  for( size_t i=0; i<inputCoordVec.size(); i++ ){
		  float			x = ( inputCoordVec[i].first - xMiddle ) / (xRange/2) * ( 1-marginRatio );
		  float			y = ( inputCoordVec[i].second - yMiddle ) / (yRange/2) * ( 1-marginRatio );
		  viewerCoordVec.push_back( make_pair(x, y) );
	  }
}

void plotMC::init()
{
	// Restore previous viewer state.
	//  restoreStateFromFile();

	dataInit();

	glutInitWindowPosition(200, 200);
	glutInitWindowSize(700,700);

	// Light default parameters
	const GLfloat light_ambient[4]  = {1.0, 1.0, 1.0, 1.0};
	const GLfloat light_specular[4] = {1.0, 1.0, 1.0, 1.0};
	const GLfloat light_diffuse[4]  = {1.0, 1.0, 1.0, 1.0};

	glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glClearDepth(1.0);
	glShadeModel (GL_SMOOTH);
	glShadeModel (GL_FLAT);
	glDisable (GL_DEPTH_TEST);
	glMatrixMode(GL_MODELVIEW);

	camera()->lookAt( sceneCenter() );
	camera()->setPosition( Vec(0, 0, 10) );
	camera()->setType( Camera::ORTHOGRAPHIC );
	camera()->setType( Camera::PERSPECTIVE );
	camera()->showEntireScene();
	camera()->setRevolveAroundPoint( sceneCenter() );
	setMouseTracking( false );

	setAxisIsDrawn( false );
}

void
plotMC::drawData(){

	//draw data
	glColor3f(1.0, 0.0, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if(1){
		//		return 0;//		return 0;
		aveX = 20.2176 ;
		aveY = 13.5567 ;
		aveZ = 5.38177 ;
		maxY = 27 ;
		maxX =20 ;
		maxZ =9 ;
		minX = 12 ;
		minY = 7 ;
		minZ = 2 ;
		rangeX = 2*( maxX - minX );
		rangeY = 2*(maxY - minY );
		rangeZ =2*( maxZ - minZ );

		for( size_t i=0; i<surf.size(); i++ ){
			glBegin( GL_TRIANGLES );
			float				normalX =  surf[i].n1.x + surf[i].n2.x + surf[i].n3.x ;
			float				normalY =  surf[i].n1.y + surf[i].n2.y + surf[i].n3.y ;
			float				normalZ =  surf[i].n1.z + surf[i].n2.z + surf[i].n3.z ;
//			glNormal3f( surf[i].n1.x, surf[i].n1.y,  surf[i].n1.z );
			glNormal3f( normalX, normalY, normalZ );
			glVertex3f( (surf[i].v1.x - aveX)/rangeX, ( surf[i].v1.y - aveY)/rangeY, (surf[i].v1.z-aveZ)/rangeZ );
			glVertex3f( (surf[i].v2.x - aveX)/rangeX, ( surf[i].v2.y - aveY)/rangeY, (surf[i].v2.z-aveZ)/rangeZ );
			glVertex3f( (surf[i].v3.x - aveX)/rangeX, ( surf[i].v3.y - aveY)/rangeY, (surf[i].v3.z-aveZ)/rangeZ );
			glEnd();

			cout<<"i:"<<i<<endl;
			cout<<"coord:"<< (surf[i].v1.x - aveX)/rangeX<<" "<< ( surf[i].v1.y - aveY)/rangeY<<" "<< (surf[i].v1.z-aveZ)/rangeZ<<endl;
			cout<<"coord:"<< (surf[i].v2.x - aveX)/rangeX<<" "<< ( surf[i].v2.y - aveY)/rangeY<<" "<< (surf[i].v2.z-aveZ)/rangeZ<<endl;
			cout<<"coord:"<< (surf[i].v3.x - aveX)/rangeX<<" "<< ( surf[i].v3.y - aveY)/rangeY<<" "<< (surf[i].v3.z-aveZ)/rangeZ<<endl<<endl;
		}
	}
	if(0){
		cout<<"trig:"<<mc.ntrigs()<<endl;
		cout<<"ver:"<<mc.nverts()<<endl;
	}
}

void
plotMC::drawAxis(){
	//draw axis
	glBegin(GL_LINES);
	glColor3f( 0, 0, 0 );
	glVertex2f( -1, -1 );
	glVertex2f( -1, 1 );
	glVertex2f( -1, -1 );
	glVertex2f( 1, -1 );
	glVertex2f( -1, -1 );
	glVertex2f( 1, 1 );
	glEnd();
}

void
plotMC::draw(){
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glClear (GL_COLOR_BUFFER_BIT );
	glColor3f( 0.0, 0.0, 0.0 );

	glLineWidth( 1.0f );
	drawData();
	drawAxis();
}

void
plotMC::testDraw(){
	for( size_t i=0; i<100; i++ ){
		glBegin( GL_LINES );

		glEnd();
	}
}

QString plotMC::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
//  return text;
}


//_____________________________________________________________________________
// Compute data
void compute_data( MarchingCubes &mc )
//-----------------------------------------------------------------------------
{
  float x,y,z      ;
  float sx,sy,sz   ;
  float tx,ty,tz   ;

  float r,R ;
  r = 1.85f ;
  R = 4 ;

  sx     = (float) mc.size_x() / 16 ;
  sy     = (float) mc.size_y() / 16 ;
  sz     = (float) mc.size_z() / 16 ;
  tx     = (float) mc.size_x() / (2*sx) ;
  ty     = (float) mc.size_y() / (2*sy) + 1.5f ;
  tz     = (float) mc.size_z() / (2*sz) ;

  for( int k = 0 ; k < mc.size_z() ; k++ )
  {
    z = ( (float) k ) / sz  - tz ;

    for( int j = 0 ; j < mc.size_y() ; j++ )
    {
      y = ( (float) j ) / sy  - ty ;

      for( int i = 0 ; i < mc.size_x() ; i++ )
      {
        x = ( (float) i ) / sx - tx ;
        mc.set_data( (float) (
          // cushin
          //            z*z*x*x - z*z*z*z - 2*z*x*x + 2*z*z*z + x*x - z*z - (x*x - z)*(x*x - z) - y*y*y*y - 2*x*x*y*y - y*y*z*z + 2*y*y*z + y*y ,
          // sphere
          //            ( (x-2)*(x-2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x+2)*(x+2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x-2)*(x-2) + (y+2)*(y+2) + (z-2)*(z-2) - 1 )) ,
          //  plane
                      x+y+z -3) ,
          // cassini
          //            (x*x + y*y + z*z + 0.45f*0.45f)*(x*x + y*y + z*z + 0.45f*0.45f) - 16*0.45f*0.45f*(x*x + z*z) - 0.5f*0.5f ,
          // blooby
//                     x*x*x*x - 5*x*x+ y*y*y*y - 5*y*y + z*z*z*z - 5*z*z + 11.8 ),
          //  chair
          //            x*x+y*y+z*z-0.95f*25)*(x*x+y*y+z*z-0.95f*25)-0.8f*((z-5)*(z-5)-2*x*x)*((z+5)*(z+5)-2*y*y ,
          // cyclide
//                      ( x*x + y*y + z*z + b*b - d*d ) * ( x*x + y*y + z*z + b*b - d*d ) - 4 * ( ( a*x - c*d ) * ( a*x - c*d ) + b*b * y*y )  ,
          // 2 torus
//                      ( ( x*x + y*y + z*z + R*R - r*r ) * ( x*x + y*y + z*z + R*R - r*r ) - 4 * R*R * ( x*x + y*y ) ) *
//                      ( ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) * ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) - 4 * R*R * ( (y+R)*(y+R) + z*z ) ) ) ,
          // mc case
          // - 26.5298*(1-x)*(1-y)*(1-z) + 81.9199*x*(1-y)*(1-z) - 100.68*x*y*(1-z) + 3.5498*(1-x)*y*(1-z)
          // + 24.1201*(1-x)*(1-y)*  z   - 74.4702*x*(1-y)*  z   + 91.5298*x*y*  z  - 3.22998*(1-x)*y*  z  ),
          // Drip
          //          x*x + y*y - 0.5*( 0.995*z*z + 0.005 - z*z*z ) +0.0025 ),  // -0.0754+0.01, -0.0025 + 0.01, grid 40^3, [-1.5,1.5]

          i,j,k ) ;
      }
    }
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________
// Compute data
void compute_ball( MarchingCubes &mc )
//-----------------------------------------------------------------------------
{
  float x,y,z      ;
  float sx,sy,sz   ;
  float tx,ty,tz   ;

  float r,R ;
  r = 1.85f ;
  R = 4 ;

  sx     = (float) mc.size_x() / 16 ; // 4
  sy     = (float) mc.size_y() / 16 ; // 4
  sz     = (float) mc.size_z() / 16 ; // 4
  tx     = (float) mc.size_x() / (2*sx) ; // 8
  ty     = (float) mc.size_y() / (2*sy) ; // 8
  tz     = (float) mc.size_z() / (2*sz) ; // 8

  for( int k = 0 ; k < mc.size_z() ; k++ )
  {
    z = ( (float) k ) / sz  - tz ;

    for( int j = 0 ; j < mc.size_y() ; j++ )
    {
      y = ( (float) j ) / sy  - ty ;

      for( int i = 0 ; i < mc.size_x() ; i++ )
      {
        x = ( (float) i ) / sx - tx ;
        mc.set_data( (float) (
          // cushin
//                      z*z*x*x - z*z*z*z - 2*z*x*x + 2*z*z*z + x*x - z*z - (x*x - z)*(x*x - z) - y*y*y*y - 2*x*x*y*y - y*y*z*z + 2*y*y*z + y*y ) ,
          // sphere
//                      ( (x-2)*(x-2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x+2)*(x+2) + (y-2)*(y-2) + (z-2)*(z-2) - 1 ) * ( (x-2)*(x-2) + (y+2)*(y+2) + (z-2)*(z-2) - 1 )) ,
          //  plane
//                      x+y+z -3) ,
          // cassini
//                      (x*x + y*y + z*z + 0.45f*0.45f)*(x*x + y*y + z*z + 0.45f*0.45f) - 16*0.45f*0.45f*(x*x + z*z) - 0.5f*0.5f ),
          // blooby
          //           x*x*x*x - 5*x*x+ y*y*y*y - 5*y*y + z*z*z*z - 5*z*z + 11.8 ),
          //  chair
          //            x*x+y*y+z*z-0.95f*25)*(x*x+y*y+z*z-0.95f*25)-0.8f*((z-5)*(z-5)-2*x*x)*((z+5)*(z+5)-2*y*y ,
          // cyclide
          //            ( x*x + y*y + z*z + b*b - d*d ) * ( x*x + y*y + z*z + b*b - d*d ) - 4 * ( ( a*x - c*d ) * ( a*x - c*d ) + b*b * y*y ),
          // 2 torus
//                      ( ( x*x + y*y + z*z + R*R - r*r ) * ( x*x + y*y + z*z + R*R - r*r ) - 4 * R*R * ( x*x + y*y ) ) *
//                      ( ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) * ( x*x + (y+R)*(y+R) + z*z + R*R - r*r ) - 4 * R*R * ( (y+R)*(y+R) + z*z ) ) ) ,
          // mc case
          // - 26.5298*(1-x)*(1-y)*(1-z) + 81.9199*x*(1-y)*(1-z) - 100.68*x*y*(1-z) + 3.5498*(1-x)*y*(1-z)
          // + 24.1201*(1-x)*(1-y)*  z   - 74.4702*x*(1-y)*  z   + 91.5298*x*y*  z  - 3.22998*(1-x)*y*  z  ),
          // Drip
                    x*x + y*y - 0.5*( 0.995*z*z + 0.005 - z*z*z ) +0.0025 ),  // -0.0754+0.01, -0.0025 + 0.01, grid 40^3, [-1.5,1.5]

          i,j,k ) ;
      }
    }
  }
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// test function
int testMC3 (int argc, char **argv)
//-----------------------------------------------------------------------------
{
	QApplication application(argc,argv);
	MarchingCubes mc ;
	mc.set_resolution( 64,64,64 ) ;

	mc.init_all() ;
//	compute_data( mc ) ;
	compute_ball( mc ) ;
	mc.run() ;
	mc.clean_temps() ;

//	  mc.writePLY("test.ply") ;
	vector<TriangleSurf>			trigVec;
	float		maxX, maxY, maxZ;
	float		minX, minY, minZ;
	float		aveX, aveY, aveZ;
	maxX = maxY = maxZ = -10e10;
	minX = minY = minZ = 10e10;
	aveX = aveY = aveZ = 0;
	if(1){
		cout<<"trig:"<<mc.ntrigs()<<endl;
		cout<<"ver:"<<mc.nverts()<<endl;
		for( int i=0; i<mc.ntrigs(); i++ ){
//			cout<<"0:"<<i<<endl;
//			cout<<"1:"<<i<<endl;
			const Triangle	*trig = mc.trig(i);
			const Vertex *v1 = mc.vert( trig->v1 );
			const Vertex *v2 = mc.vert( trig->v2 );
			const Vertex *v3 = mc.vert( trig->v3 );

			Coord	c1( v1->x, v1->y, v1->z );
			Coord	c2( v2->x, v2->y, v2->z );
			Coord	c3( v3->x, v3->y, v3->z );

			Coord	n1( v1->nx, v1->ny, v1->nz );
			Coord	n2( v2->nx, v2->ny, v2->nz );
			Coord	n3( v3->nx, v3->ny, v3->nz );

			maxX = maxX > v1->x ? maxX:v1->x ;
			maxY = maxY > v1->y ? maxY:v1->y ;
			maxZ = maxZ > v1->z ? maxZ:v1->z ;

			minX = minX < v1->x ? minX:v1->x ;
			minY = minY < v1->y ? minY:v1->y ;
			minZ = minZ < v1->z ? minZ:v1->z ;

			aveX += v1->x ;
			aveY += v1->y ;
			aveZ += v1->z ;

			TriangleSurf	triangle( c1, n1, c2, n2, c3, n3 );
			trigVec.push_back( triangle );
		}
		aveX = aveX /(float)mc.ntrigs();
		aveY = aveY /(float)mc.ntrigs();
		aveZ = aveZ /(float)mc.ntrigs();
		cout<<"ave: "<<aveX<<" "<<aveY<<" "<<aveZ<<endl;
		cout<<"max: "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
		cout<<"min: "<<minX<<" "<<minY<<" "<<minZ<<endl;
	}


	plotMC			pl( trigVec );
	pl.show();
	mc.clean_all() ;
	return application.exec()  ;
}
//_____________________________________________________________________________

int testMC (int argc, char **argv){
	QApplication application(argc,argv);
	string							fileNameA ( "/mnt/2t/xushutan/project/gold_dock/1a9u_diverse500/ligand_corina.mol2" ) ;
	Molecular						molA( fileNameA );
	molA = readMolecularFile( fileNameA ).front();

	MarchingCubes mc ;
//	mc.set_resolution( 60,60,60 ) ;

//	mc.init_all() ;
	compute_VdwMarchingCube( molA.get_atomVec() , mc ) ;
	mc.run() ;
	mc.clean_temps() ;

	cout<<"trig:"<<mc.ntrigs()<<endl;
	cout<<"ver:"<<mc.nverts()<<endl;
	cout<<"-----------------"<<endl;

	vector<TriangleSurf>			trigVec;
	float		maxX, maxY, maxZ;
	float		minX, minY, minZ;
	float		aveX, aveY, aveZ;
	maxX = maxY = maxZ = -10e10;
	minX = minY = minZ = 10e10;
	aveX = aveY = aveZ = 0;
	if(1){
		cout<<"trig:"<<mc.ntrigs()<<endl;
		cout<<"ver:"<<mc.nverts()<<endl;
		for( int i=0; i<mc.ntrigs(); i++ ){
//			cout<<"0:"<<i<<endl;
//			cout<<"1:"<<i<<endl;
			const Triangle	*trig = mc.trig(i);
			const Vertex *v1 = mc.vert( trig->v1 );
			const Vertex *v2 = mc.vert( trig->v2 );
			const Vertex *v3 = mc.vert( trig->v3 );

			Coord	c1( v1->x, v1->y, v1->z );
			Coord	c2( v2->x, v2->y, v2->z );
			Coord	c3( v3->x, v3->y, v3->z );

			Coord	n1( v1->nx, v1->ny, v1->nz );
			Coord	n2( v2->nx, v2->ny, v2->nz );
			Coord	n3( v3->nx, v3->ny, v3->nz );

			maxX = maxX > v1->x ? maxX:v1->x ;
			maxY = maxY > v1->y ? maxY:v1->y ;
			maxZ = maxZ > v1->z ? maxZ:v1->z ;

			minX = minX < v1->x ? minX:v1->x ;
			minY = minY < v1->y ? minY:v1->y ;
			minZ = minZ < v1->z ? minZ:v1->z ;

			aveX += v1->x ;
			aveY += v1->y ;
			aveZ += v1->z ;

			TriangleSurf	triangle( c1, n1, c2, n2, c3, n3 );
			trigVec.push_back( triangle );
		}
		aveX = aveX /(float)mc.ntrigs();
		aveY = aveY /(float)mc.ntrigs();
		aveZ = aveZ /(float)mc.ntrigs();
		cout<<"ave: "<<aveX<<" "<<aveY<<" "<<aveZ<<endl;
		cout<<"max: "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
		cout<<"min: "<<minX<<" "<<minY<<" "<<minZ<<endl;
	}

	plotMC			pl( trigVec );
	pl.show();
	mc.clean_all() ;
	return application.exec()  ;
	return 0;
}
*/

#endif /* TESTMARCHINGCUBES_H_ */
