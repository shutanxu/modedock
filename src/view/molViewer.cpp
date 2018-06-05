/*
 * molViewer.cpp
 *
 *  Created on: Aug 25, 2014
 *      Author: stan
 */

#include <QMouseEvent>
#include <QTimer>
#include <QVector2D>
#include <QVector3D>
#include <QOpenGLShaderProgram>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QObject>
#include <QFile>
#include <QDateTime>
#include <QFileSystemWatcher>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions>
#include<math.h>
#include <qmath.h>
#include <QFileInfo>
#include <QTime>
#include <QStaticText>

#include"molViewer.h"

MolecularViewer::MolecularViewer(QWidget *parent):QOpenGLWidget(parent){
//    setFormat( QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer) );

    m_trackBalls = new TrackBall(0.05f, QVector3D(0, 1, 0), TrackBall::Sphere);
    translation = QVector2D();
    m_distExp= 2000;
    m_trackMover= new trackMover();
    showAxis=false;

    showLine=showStick=showSphere=showSurface=showResSurface=false;

    QPointF p;
    m_trackBalls->push( pixelPosToViewPos(p), QQuaternion());
    rotation = m_trackBalls->rotation();
    update();
}

MolecularViewer::~MolecularViewer(){
    delete m_trackBalls;
}

void
MolecularViewer::initializeGL(){

    initializeOpenGLFunctions();

    glClearColor(0, 0, 0, 1);
    glEnable(GL_DEPTH_TEST);

    glEnable( GL_POINT_SMOOTH );
    glEnable(GL_LINE_SMOOTH);

    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glFrontFace(GL_CW);
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    initializeAxis();
    initializeLine();
    initializeStick();
    initializeSphere();
    initializeSurface();
    initializeResSurface();
}

void MolecularViewer::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if( showAxis ){
        drawAxis();
    }

    if(  !molVec.empty() ){

    	if( true ){
    		drawLine();
    	}
        if( !dotsData.empty() ){
            drawBindingSite();
        }
        if( !stickBondVec.empty() && showStick ){
            drawStick();
        }
        if( !sphereModeAtomsVec.empty() && showSphere ){
            drawSphere();
        }
        if( !surfaceAtomVec.empty() && showSurface ){
            drawSphere();
            drawSurface();
        }
        if( !resSurfaceAtomVec.empty() && showResSurface ){
            drawSphere();
            drawResSurface();
        }
    }

}

void
MolecularViewer::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    orthographic.setToIdentity();
    double hScale = static_cast<double>(h) / static_cast<double>(500);
    double wScale = static_cast<double>(w) / static_cast<double>(500);
    orthographic.ortho(-wScale, wScale, -hScale, hScale, -100, 100);

    window_normalised_matrix.setToIdentity();
    window_normalised_matrix.translate(w/2.0, h/2.0);

    window_painter_matrix.setToIdentity();
    window_painter_matrix.translate(w/2.0, h/2.0);

    projection.setToIdentity();
    projection.perspective(45.f, qreal(w)/qreal(h), 0.1f, 100.f);
}

void
MolecularViewer::setMolVec( const vector<Molecular> mol ){
    molVec=mol;
    vector<Atom> atomVec;
    vector<Bond> bdVec;

    for( size_t i=0; i<molVec.size(); i++ ){

        vector<Atom> atVec=molVec[i].get_atomVec();

        for( size_t j=0; j<atVec.size(); j++ ){
            atomVec.push_back( atVec[j] );
        }
        vector<Bond> bd=molVec[i].get_bondVec();

        for( size_t j=0; j<bd.size(); j++ ){
            bdVec.push_back( bd[j] );
        }
    }

    noOfAtom = atomVec.size();

    float range = atomsMaxRange( atomVec );
    maxCoordRange = range;
    atomsCenter = centerOfAtoms( atomVec );
}

void
MolecularViewer::initializeAxis(){
    string vshaderFile = string("./shader/axis_vs.glsl");
    string fshaderFile = string("./shader/axis_fs.glsl");

    if( !axisProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !axisProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !axisProgram.link() ){
        cerr<<"can not link shader#!!!"<<endl;
    }
    if( !axisProgram.bind() ){
        cerr<<"can not bind shader##!!!"<<endl;
    }
}

void
MolecularViewer::initializeLine(){
    string vshaderFile = string("./shader/line_vs.glsl");
    string fshaderFile = string("./shader/line_fs.glsl");

    if( !lineProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !lineProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !lineProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !lineProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }
}

void
MolecularViewer::updateLine( const vector<Molecular> mol ){

    vector<Bond> bdVec;

    for( size_t i=0; i<mol.size(); i++ ){
        vector<Bond> bd=mol[i].get_bondVec();
        for( size_t j=0; j<bd.size(); j++ ){
            bdVec.push_back( bd[j] );
        }
    }
    vector<float> coordVec, colorVec;

    for( size_t i=0; i<bdVec.size(); i++ ){
        Atom a1 = bdVec[i].get_firstAtom();
        Coord c1 = a1.get_coord();
        Atom a2 = bdVec[i].get_secAtom();
        Coord c2 = a2.get_coord();

        if( a1.get_type() == a2.get_type() ){
            coordVec.push_back((c1.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c1.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c1.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a1.get_color()[0] );
            colorVec.push_back( a1.get_color()[1] );
            colorVec.push_back( a1.get_color()[2] );

            coordVec.push_back((c2.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c2.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c2.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a2.get_color()[0] );
            colorVec.push_back( a2.get_color()[1] );
            colorVec.push_back( a2.get_color()[2] );
        }else{
            Coord med;
            med.x = ( a1.get_coord().x + a2.get_coord().x )/2;
            med.y = ( a1.get_coord().y + a2.get_coord().y )/2;
            med.z = ( a1.get_coord().z + a2.get_coord().z )/2;

            coordVec.push_back((c1.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c1.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c1.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a1.get_color()[0] );
            colorVec.push_back( a1.get_color()[1] );
            colorVec.push_back( a1.get_color()[2] );

            coordVec.push_back((med.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((med.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((med.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a1.get_color()[0] );
            colorVec.push_back( a1.get_color()[1] );
            colorVec.push_back( a1.get_color()[2] );

            coordVec.push_back((c2.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c2.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c2.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a2.get_color()[0] );
            colorVec.push_back( a2.get_color()[1] );
            colorVec.push_back( a2.get_color()[2] );

            coordVec.push_back((med.x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((med.y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((med.z-atomsCenter.z)/maxCoordRange);
            colorVec.push_back( a2.get_color()[0] );
            colorVec.push_back( a2.get_color()[1] );
            colorVec.push_back( a2.get_color()[2] );
        }
    }

     lineVertexBuf.create();
     lineVertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
     lineVertexBuf.bind();
     lineVertexBuf.allocate( &coordVec[0], coordVec.size()*sizeof(float) );

     lineColorBuf.create();
     lineColorBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
     lineColorBuf.bind();
     lineColorBuf.allocate( &colorVec[0], colorVec.size()*sizeof(float) );

     lineNum = coordVec.size()/3;
}

void
MolecularViewer::initializeStick(){
    string vshaderFile = string("./shader/stick_vs.glsl");
    string fshaderFile = string("./shader/stick_fs.glsl");

    if( !stickProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !stickProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !stickProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !stickProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }
}

void
MolecularViewer::updateStick( const vector<Molecular> mol ){
    if(mol.empty()){
        stickBondVec.clear();
        return;
    }

    vector<Bond> bdVec;
    for( size_t i=0; i<mol.size(); i++ ){
      vector<Bond> b= mol[i].get_bondVec();
      for( size_t j=0; j<b.size(); j++ ){
          bdVec.push_back( b[j] );
      }
    }
    stickBondVec=bdVec;

    vector<QVector3D>    oneVertexVec,oneNormalVec;
    vector<QColor>       oneColorVec;
    vector<QColor>       stick_colorVec;

    stick_colorVec.clear();
    stick_vertexVec.clear();
    stick_normalVec.clear();

    for( size_t i=0; i<bdVec.size(); i++ ){

        oneVertexVec.clear();
        oneNormalVec.clear();
        oneColorVec.clear();
        Atom a1 = bdVec[i].get_firstAtom();
        Coord c1 = a1.get_coord();
        QVector3D p1( ( c1.x-atomsCenter.x)/maxCoordRange,
                               ( c1.y-atomsCenter.y)/maxCoordRange,
                               ( c1.z-atomsCenter.z)/maxCoordRange );

        QColor color1( a1.get_color()[0], a1.get_color()[1], a1.get_color()[2] );
        Atom a2 = bdVec[i].get_secAtom();
        Coord c2 = a2.get_coord();
        QVector3D p2( (c2.x-atomsCenter.x)/maxCoordRange,
                                (c2.y-atomsCenter.y)/maxCoordRange,
                                (c2.z-atomsCenter.z)/maxCoordRange );

        QColor color2( a2.get_color()[0], a2.get_color()[1], a2.get_color()[2] );

        float radius = 0.3/maxCoordRange;
        getStick(p1,p2,color1,color2,0.5, radius, oneVertexVec,oneNormalVec, oneColorVec);

        stick_vertexVec.insert(stick_vertexVec.end(),oneVertexVec.begin(),oneVertexVec.end());
        stick_normalVec.insert(stick_normalVec.end(),oneNormalVec.begin(),oneNormalVec.end());
        stick_colorVec.insert(stick_colorVec.end(),oneColorVec.begin(),oneColorVec.end());
    }

    vector<float> colorVec;
    for( size_t i=0; i<stick_colorVec.size(); i++ ){
        colorVec.push_back(stick_colorVec[i].red());
        colorVec.push_back(stick_colorVec[i].green());
        colorVec.push_back(stick_colorVec[i].blue());
    }

    stickVertexBuf.create();
    stickVertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    stickVertexBuf.bind();
    stickVertexBuf.allocate( &stick_vertexVec[0], stick_vertexVec.size()*sizeof(QVector3D) );

    stickNormalBuf.create();
    stickNormalBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    stickNormalBuf.bind();
    stickNormalBuf.allocate( &stick_normalVec[0], stick_normalVec.size()*sizeof(QVector3D) );

    stickColorBuf.create();
    stickColorBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    stickColorBuf.bind();
    stickColorBuf.allocate( &colorVec[0], colorVec.size()*sizeof(float) );
}

void
MolecularViewer::initializeSphere(){

    string vshaderFile = string("./shader/sphere_vs.glsl");
    string fshaderFile = string("./shader/sphere_fs.glsl");

    if( !sphereProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !sphereProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !sphereProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !sphereProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }
}

void
MolecularViewer::updateSphere( const vector<Molecular> mol ){

    if( mol.empty() ){
        sphereModeAtomsVec.clear();
        surfaceAtomVec.clear();
        return;
    }

    vector<Atom> atomVec;

    for( size_t i=0; i<mol.size(); i++ ){
      vector<Atom> at= mol[i].get_atomVec();
      for( size_t j=0; j<at.size(); j++ ){
          atomVec.push_back( at[j] );
      }
    }
    sphereModeAtomsVec = atomVec;

    Polyeder ply(3);
    surfs = ply.get_surfTriangles();

    float           vertex[surfs.size()*3][3];
    for( size_t i=0; i<surfs.size(); i++ ){
        vertex[i*3][0] = surfs[i].v1[0];
        vertex[i*3][1] = surfs[i].v1[1];
        vertex[i*3][2] = surfs[i].v1[2];

        vertex[i*3+1][0] = surfs[i].v2[0];
        vertex[i*3+1][1] = surfs[i].v2[1];
        vertex[i*3+1][2] = surfs[i].v2[2];

        vertex[i*3+2][0] = surfs[i].v3[0];
        vertex[i*3+2][1] = surfs[i].v3[1];
        vertex[i*3+2][2] = surfs[i].v3[2];
    }

    vector<float> allVertex = ply.get_allVertex();
    vector<int> triangleVertexIndex = ply.get_triangleVertexIndex();

    sphereVertexBuf.create();
    sphereVertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    sphereVertexBuf.bind();
    sphereVertexBuf.allocate( &allVertex[0], allVertex.size()*sizeof(float) );

    QOpenGLBuffer indbuf( QOpenGLBuffer::IndexBuffer);
    sphereIndexBuf = indbuf;
    sphereIndexBuf.create();
    sphereIndexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw );
    sphereIndexBuf.bind();
    sphereIndexBuf.allocate( &triangleVertexIndex[0], triangleVertexIndex.size()*sizeof(int) );
}

void
MolecularViewer::initializeSurface(){
    glDisable(GL_CULL_FACE);
    string vshaderFile = string("./shader/surface_vs.glsl");
    string fshaderFile = string("./shader/surface_fs.glsl");

    if( !surfaceProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !surfaceProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !surfaceProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !surfaceProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }
}

void
MolecularViewer::updateSurface(const vector<Molecular> mol){

    if( mol.empty() ){
        surfaceAtomVec.clear();
        sphereModeAtomsVec.clear();
        return;
    }

    updateSphere(mol);

    vector<Atom> atomVec;
    for( size_t i=0; i<mol.size(); i++ ){
      vector<Atom> at= mol[i].get_atomVec();
      for( size_t j=0; j<at.size(); j++ ){
          atomVec.push_back( at[j] );
      }
    }
    surfaceAtomVec = atomVec;

    vector<reducedSurface>  rsVec;
    vector<sphereSurface>   ssVec;
    float		probeRadius = 1.4;

    vector<float> coord;
    vector<float> normal;

    surfaceTrigNum = coord.size() / 9;

    //----------------------arch surface--------------------------------
    coord.clear();
    normal.clear();

    vector<reducedSurface_connect> rscVec;
    getReducedSurface(atomVec, probeRadius, rsVec, ssVec, rscVec);

    for( size_t i=0; i<rscVec.size(); i++ ){
        Atom a1 = rscVec[i].get_at1();
        Atom a2 = rscVec[i].get_at2();
        Atom a3 = rscVec[i].get_at3();
        Atom a4 = rscVec[i].get_at4();

        vector<Coord> coo, nor;
        get_SES_toric_reentrant_face( rscVec[i], 10, coo, nor);

        for( size_t i=0; i<coo.size(); i++ ){
            coord.push_back( (coo[i].x- atomsCenter.x)/maxCoordRange );
            coord.push_back( (coo[i].y- atomsCenter.y)/maxCoordRange );
            coord.push_back( (coo[i].z- atomsCenter.z)/maxCoordRange );
            normal.push_back( -nor[i].x );
            normal.push_back( -nor[i].y );
            normal.push_back( -nor[i].z );
        }
    }
    surfaceTrigNum = coord.size() / 9;
    for( size_t i=0; i<ssVec.size(); i++ ){
        vector<Coord> coo, nor;
        get_SES_spheric_reentrant_face( ssVec[i], coo, nor );
        for( size_t k=0; k<coo.size(); k++ ){
            coord.push_back( (coo[k].x- atomsCenter.x)/maxCoordRange );
            coord.push_back( (coo[k].y- atomsCenter.y)/maxCoordRange );
            coord.push_back( (coo[k].z- atomsCenter.z)/maxCoordRange );
            normal.push_back( nor[k].x );
            normal.push_back( nor[k].y );
            normal.push_back( nor[k].z );
        }
    }
    surfaceTrigNum = coord.size() / 9;

    surfaceVertexBuf.create();
    surfaceVertexBuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    surfaceVertexBuf.bind();
    surfaceVertexBuf.allocate( &coord[0], coord.size()*sizeof(float) );

    surfaceNormalBuf.create();
    surfaceNormalBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    surfaceNormalBuf.bind();
    surfaceNormalBuf.allocate(&normal[0], normal.size()*sizeof(float));
}


void
MolecularViewer::initializeResSurface(){
    glDisable(GL_CULL_FACE);
    string vshaderFile = string("./shader/resSurface_vs.glsl");
    string fshaderFile = string("./shader/resSurface_fs.glsl");

    if( !resSurfaceProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !resSurfaceProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !resSurfaceProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !resSurfaceProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }
}

void
MolecularViewer::updateResSurface(const vector<Molecular> mol){

    if( mol.empty() ){
        resSurfaceAtomVec.clear();
        sphereModeAtomsVec.clear();
        return;
    }

    updateSphere(mol);

    vector<Atom> atomVec;
    for( size_t i=0; i<mol.size(); i++ ){
      vector<Atom> at= mol[i].get_atomVec();
      for( size_t j=0; j<at.size(); j++ ){
          atomVec.push_back( at[j] );
      }
    }
    resSurfaceAtomVec = atomVec;

    vector<reducedSurface>  rsVec;
    vector<sphereSurface>   ssVec;
    float		probeRadius = 1.4;

    vector<float> coord;
    vector<float> normal;
    vector<float> color;
    coord.clear();
    normal.clear();
    color.clear();

    resSurfaceTrigNum = coord.size() / 9;

    //----------------------arch surface--------------------------------
    vector<reducedSurface_connect> rscVec;
    getReducedResSurface(atomVec, probeRadius, rsVec, ssVec, rscVec);

    for( size_t i=0; i<rscVec.size(); i++ ){

        vector<Coord> coo, nor;
        get_SES_toric_reentrant_face( rscVec[i], 10, coo, nor);
        vector<float> tempColor=getColorFromResidue(rscVec[i].get_at3().get_residueIndex());

        for( size_t i=0; i<coo.size(); i++ ){
            coord.push_back( (coo[i].x- atomsCenter.x)/maxCoordRange );
            coord.push_back( (coo[i].y- atomsCenter.y)/maxCoordRange );
            coord.push_back( (coo[i].z- atomsCenter.z)/maxCoordRange );
            normal.push_back( -nor[i].x );
            normal.push_back( -nor[i].y );
            normal.push_back( -nor[i].z );
            color.push_back(tempColor[0]);
            color.push_back(tempColor[1]);
            color.push_back(tempColor[2]);
        }
    }
    resSurfaceTrigNum = coord.size() / 9;
    for( size_t i=0; i<ssVec.size(); i++ ){
        vector<Coord> coo, nor;
        get_SES_spheric_reentrant_face( ssVec[i], coo, nor );
        vector<float> tempColor=getColorFromResidue(ssVec[i].get_at1().get_residueIndex());
        for( size_t k=0; k<coo.size(); k++ ){
            coord.push_back( (coo[k].x- atomsCenter.x)/maxCoordRange );
            coord.push_back( (coo[k].y- atomsCenter.y)/maxCoordRange );
            coord.push_back( (coo[k].z- atomsCenter.z)/maxCoordRange );
            normal.push_back( nor[k].x );
            normal.push_back( nor[k].y );
            normal.push_back( nor[k].z );
            color.push_back(tempColor[0]);
            color.push_back(tempColor[1]);
            color.push_back(tempColor[2]);

        }
    }
    resSurfaceTrigNum = coord.size() / 9;

    resSurfaceVertexBuf.create();
    resSurfaceVertexBuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    resSurfaceVertexBuf.bind();
    resSurfaceVertexBuf.allocate(&coord[0], coord.size()*sizeof(float) );

    resSurfaceNormalBuf.create();
    resSurfaceNormalBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    resSurfaceNormalBuf.bind();
    resSurfaceNormalBuf.allocate(&normal[0], normal.size()*sizeof(float));

    resSurfaceColorBuf.create();
    resSurfaceColorBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    resSurfaceColorBuf.bind();
    resSurfaceColorBuf.allocate(&color[0], color.size()*sizeof(float));
}

void
MolecularViewer::initializeHelix(){
    string vshaderFile = string("./shader/helix.vs");
    string fshaderFile = string("./shader/helix.fs");

    if( !helixProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !helixProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !helixProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !helixProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    //light
    helixProgram.setUniformValue("lightDirect", QVector3D(0.0, 0.0, 0.5));
    helixProgram.setUniformValue("ambientIntensity", float(0.1) );
    helixProgram.setUniformValue("diffuseIntensity", float(1.0) );
    helixProgram.setUniformValue("specularIntensity", float(0.5) );

    //----------------------------------------------------------------

    //-- atom
    string fileName("1A9U.pdb");

    vector<Molecular> molVec = readMolecularFile(fileName);
    vector<Atom> atomVec = molVec.front().get_atomVec();
    noOfAtom = atomVec.size();
    sphereModeAtomsVec = atomVec;

    float maxRange = atomsMaxRange( atomVec );

    atomsCenter = centerOfAtoms( atomVec );

    //----------------------------------------------------------------
    QVector3D catom0;
    QVector3D catom1;
    QVector3D catom2;
    QVector3D catom3;
    QVector3D catom4;
    QVector3D catom5;
    QVector3D catom6;

    QVector3D cat1, cat2, cat3, cat4;
    QVector3D catt0, catt1, catt2, catt3, catt4;

    QVector3D T1, T2;

    QVector3D CN1, CN2, N1, N2, N3, N4;

    QVector3D A, A1, A2, A3;
    QVector3D B, B1, B2, B3;
    QVector3D C, C1, C2, C3;
    QVector3D D, D1, D2, D3;

    vector<vector<Atom> > helixAtomVec ;
    vector<Atom>          atVec;
    for( size_t i=0; i<atomVec.size(); i++ ){
        if( atomVec[i].get_residueIndex() > 61 && atomVec[i].get_residueIndex() < 76 ){
            if( atomVec[i].get_name() == "CA" ){
                atVec.push_back( atomVec[i] );
            }
        }
    }
    helixAtomVec.push_back( atVec );
    Coord center = centerOfAtoms( atVec );
    quadsCoord.clear();
    quadsNormal.clear();

    for( size_t i=0; i<helixAtomVec.size(); i++ ){
        for( size_t j=0; j<helixAtomVec[i].size() - 3; j++ ){
            catom1.setX( helixAtomVec[i][j].get_coord().x - center.x );
            catom1.setY( helixAtomVec[i][j].get_coord().y - center.y );
            catom1.setZ( helixAtomVec[i][j].get_coord().z - center.z );

            catom2.setX( helixAtomVec[i][j+1].get_coord().x - center.x);
            catom2.setY( helixAtomVec[i][j+1].get_coord().y - center.y);
            catom2.setZ( helixAtomVec[i][j+1].get_coord().z - center.z);

            catom3.setX( helixAtomVec[i][j+2].get_coord().x - center.x);
            catom3.setY( helixAtomVec[i][j+2].get_coord().y - center.y);
            catom3.setZ( helixAtomVec[i][j+2].get_coord().z - center.z);

            catom4.setX( helixAtomVec[i][j+3].get_coord().x - center.x);
            catom4.setY( helixAtomVec[i][j+3].get_coord().y - center.y);
            catom4.setZ( helixAtomVec[i][j+3].get_coord().z - center.z);

            CN1 = catom2 - catom1;
            CN2 = catom2 - catom3;
            N2 = QVector3D::crossProduct(CN2, CN1) ;

            N1 = N2;

            CN1 = catom3 - catom2;
            CN2 = catom3 - catom4;
            N3 = QVector3D::crossProduct(CN2, CN1);

            N4 = N3;

            if( fabs(N1.length() ) >1.0E-10  ){
                N1.normalize();
            }
            if( fabs(N2.length() ) >1.0E-10 ){
                N2.normalize();
            }
            if( fabs(N3.length() ) >1.0E-10  ){
                N3.normalize();
            }
            if( fabs(N4.length() ) >1.0E-10  ){
                N4.normalize();
            }

            cat1 = catom1 -  N1;
            cat2 = catom2 -  N2;
            cat3 = catom3 -  N3;
            cat4 = catom4 -  N4;
            T1 = cat3 - cat1;
            T2 = cat4 - cat2;
            D = cat2;
            C = T1;
            B = 3*cat3 - 3*cat2 - 2*T1 - T2;
            A = T2 + T1 - 2*cat3 + 2*cat2;

            catt1 = catom1 +   N1;
            catt2 = catom2 +   N2;
            catt3 = catom3 +   N3;
            catt4 = catom4 +   N4;
            T1 = catt3 - catt1;
            T2 = catt4 - catt2;
            D1 = catt2;
            C1 = T1;
            B1 = 3*catt3 - 3*catt2 - 2*T1 - T2;
            A1 = T2 + T1 - 2*catt3 + 2*catt2;

            for( float t = 0.0; t<1.0; t = t+0.01 ){
                QVector3D p, p1, p2, p3;
                QVector3D q, q1, q2, q3;
                float	t1 = t + 0.01;

                // outer surface of helix
                p = A * t * t * t + B * t * t + C * t + D;
                p1 = A1 * t * t * t + B1 * t * t + C1 * t + D1;
                p2 = A1 * t1 * t1 * t1 + B1 * t1 * t1 + C1 * t1 + D1;
                p3 = A * t1 * t1 * t1 + B * t1 * t1 + C * t1 + D;
                QVector3D	normal = QVector3D::crossProduct( p - p2, p - p1 );
                if( fabs(normal.length() ) >1.0E-10 ){
                    normal.normalize();
                }
                glBegin( GL_QUADS );
                glNormal3f( normal.x(), normal.y(), normal.z() );
                glVertex3f( p.x() / maxRange, p.y()/maxRange, p.z()/maxRange );
                glVertex3f( p1.x() / maxRange, p1.y()/maxRange, p1.z()/maxRange );
                glVertex3f( p2.x() / maxRange, p2.y()/maxRange, p2.z()/maxRange );
                glVertex3f( p3.x() / maxRange, p3.y()/maxRange, p3.z()/maxRange );
                glEnd();

                //---------------------
                quadsCoord.push_back( p.x() / maxRange );
                quadsCoord.push_back( p.y() / maxRange );
                quadsCoord.push_back( p.z() / maxRange );
                quadsCoord.push_back( p1.x() / maxRange );
                quadsCoord.push_back( p1.y() / maxRange );
                quadsCoord.push_back( p1.z() / maxRange );
                quadsCoord.push_back( p3.x() / maxRange );
                quadsCoord.push_back( p3.y() / maxRange );
                quadsCoord.push_back( p3.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );

                quadsCoord.push_back( p2.x() / maxRange );
                quadsCoord.push_back( p2.y() / maxRange );
                quadsCoord.push_back( p2.z() / maxRange );
                quadsCoord.push_back( p3.x() / maxRange );
                quadsCoord.push_back( p3.y() / maxRange );
                quadsCoord.push_back( p3.z() / maxRange );
                quadsCoord.push_back( p1.x() / maxRange );
                quadsCoord.push_back( p1.y() / maxRange );
                quadsCoord.push_back( p1.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );

                //---------------------

                // inner surface of helix
                q = p - 0.15 * normal;
                q1 = p1 - 0.15 * normal;
                q2 = p2 - 0.15 * normal;
                q3 = p3 - 0.15 * normal;
                normal = QVector3D::crossProduct( q - q1, q-q2 );
                if( fabs(normal.length() ) >1.0E-10 ){
                    normal.normalize();
                }
                glBegin( GL_QUADS );
                glNormal3f( normal.x(), normal.y(), normal.z() );
                glVertex3f( q.x() / maxRange, q.y()/maxRange, q.z()/maxRange );
                glVertex3f( q1.x() / maxRange, q1.y()/maxRange, q1.z()/maxRange );
                glVertex3f( q2.x() / maxRange, q2.y()/maxRange, q2.z()/maxRange );
                glVertex3f( q3.x() / maxRange, q3.y()/maxRange, q3.z()/maxRange );
                glEnd();

                //---------------------
                quadsCoord.push_back( q.x() / maxRange );
                quadsCoord.push_back( q.y() / maxRange );
                quadsCoord.push_back( q.z() / maxRange );
                quadsCoord.push_back( q1.x() / maxRange );
                quadsCoord.push_back( q1.y() / maxRange );
                quadsCoord.push_back( q1.z() / maxRange );
                quadsCoord.push_back( q3.x() / maxRange );
                quadsCoord.push_back( q3.y() / maxRange );
                quadsCoord.push_back( q3.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );

                quadsCoord.push_back( q2.x() / maxRange );
                quadsCoord.push_back( q2.y() / maxRange );
                quadsCoord.push_back( q2.z() / maxRange );
                quadsCoord.push_back( q3.x() / maxRange );
                quadsCoord.push_back( q3.y() / maxRange );
                quadsCoord.push_back( q3.z() / maxRange );
                quadsCoord.push_back( q1.x() / maxRange );
                quadsCoord.push_back( q1.y() / maxRange );
                quadsCoord.push_back( q1.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                //---------------------

                // lower edge
                normal = QVector3D::crossProduct( q - q3, q - p );
                if( fabs(normal.length() ) >1.0E-10 ){
                    normal.normalize();
                }
                glBegin( GL_QUADS );
                glNormal3f( normal.x(), normal.y(), normal.z() );
                glVertex3f( q.x() / maxRange, q.y()/maxRange, q.z()/maxRange );
                glVertex3f( p.x() / maxRange, p.y()/maxRange, p.z()/maxRange );
                glVertex3f( p3.x() / maxRange, p3.y()/maxRange, p3.z()/maxRange );
                glVertex3f( q3.x() / maxRange, q3.y()/maxRange, q3.z()/maxRange );
                glEnd();

                //---------------------
                quadsCoord.push_back( q.x() / maxRange );
                quadsCoord.push_back( q.y() / maxRange );
                quadsCoord.push_back( q.z() / maxRange );
                quadsCoord.push_back( p.x() / maxRange );
                quadsCoord.push_back( p.y() / maxRange );
                quadsCoord.push_back( p.z() / maxRange );
                quadsCoord.push_back( q3.x() / maxRange );
                quadsCoord.push_back( q3.y() / maxRange );
                quadsCoord.push_back( q3.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );

                quadsCoord.push_back( p3.x() / maxRange );
                quadsCoord.push_back( p3.y() / maxRange );
                quadsCoord.push_back( p3.z() / maxRange );
                quadsCoord.push_back( q3.x() / maxRange );
                quadsCoord.push_back( q3.y() / maxRange );
                quadsCoord.push_back( q3.z() / maxRange );
                quadsCoord.push_back( p.x() / maxRange );
                quadsCoord.push_back( p.y() / maxRange );
                quadsCoord.push_back( p.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                //---------------------

                // upper edge
                normal = QVector3D::crossProduct( q1 - p1,  q1 - q2 );
                if( fabs(normal.length() ) >1.0E-10 ){
                    normal.normalize();
                }
                glBegin( GL_QUADS );
                glNormal3f( normal.x(), normal.y(), normal.z() );
                glVertex3f( q1.x() / maxRange, q1.y()/maxRange, q1.z()/maxRange );
                glVertex3f( p1.x() / maxRange, p1.y()/maxRange, p1.z()/maxRange );
                glVertex3f( p2.x() / maxRange, p2.y()/maxRange, p2.z()/maxRange );
                glVertex3f( q2.x() / maxRange, q2.y()/maxRange, q2.z()/maxRange );
                glEnd();

                //---------------------
                quadsCoord.push_back( q1.x() / maxRange );
                quadsCoord.push_back( q1.y() / maxRange );
                quadsCoord.push_back( q1.z() / maxRange );
                quadsCoord.push_back( p1.x() / maxRange );
                quadsCoord.push_back( p1.y() / maxRange );
                quadsCoord.push_back( p1.z() / maxRange );
                quadsCoord.push_back( q2.x() / maxRange );
                quadsCoord.push_back( q2.y() / maxRange );
                quadsCoord.push_back( q2.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );

                quadsCoord.push_back( p2.x() / maxRange );
                quadsCoord.push_back( p2.y() / maxRange );
                quadsCoord.push_back( p2.z() / maxRange );
                quadsCoord.push_back( q2.x() / maxRange );
                quadsCoord.push_back( q2.y() / maxRange );
                quadsCoord.push_back( q2.z() / maxRange );
                quadsCoord.push_back( p1.x() / maxRange );
                quadsCoord.push_back( p1.y() / maxRange );
                quadsCoord.push_back( p1.z() / maxRange );

                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                quadsNormal.push_back( normal.x() );
                quadsNormal.push_back( normal.y() );
                quadsNormal.push_back( normal.z() );
                //---------------------
            }
        }
    }

    glGenBuffers(1, &helixVBO);
    glBindBuffer(GL_ARRAY_BUFFER, helixVBO);

    glBufferData( GL_ARRAY_BUFFER, (quadsCoord.size()*sizeof(float)+quadsNormal.size()*sizeof(float)),
            NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, quadsCoord.size()*sizeof(float), &quadsCoord[0]);
    glBufferSubData(GL_ARRAY_BUFFER, quadsCoord.size()*sizeof(float),
            quadsNormal.size()*sizeof(float), &quadsNormal[0]);

}

void
MolecularViewer::initializeSheet(){
    string vshaderFile = string("./shader/sheet.vs");
    string fshaderFile = string("./shader/sheet.fs");

    if( !sheetProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !sheetProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !sheetProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !sheetProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    //light
    sheetProgram.setUniformValue("lightDirect", QVector3D(0.0, 0.0, 0.5));
    sheetProgram.setUniformValue("ambientIntensity", float(0.1) );
    sheetProgram.setUniformValue("diffuseIntensity", float(1.0) );
    sheetProgram.setUniformValue("specularIntensity", float(0.5) );

    //----------------------------------------------------------------

    //-- atom
//	string fileName("ligand.mol2");
    string fileName("1A9U.pdb");
//	string fileName("4CYX.pdb");

    vector<Molecular> molVec = readMolecularFile(fileName);
    vector<Atom> atomVec = molVec.front().get_atomVec();
    noOfAtom = atomVec.size();
//    cout<<"noOfAtom: "<<noOfAtom<<endl;
    sphereModeAtomsVec = atomVec;

    float maxRange = atomsMaxRange( atomVec );

    atomsCenter = centerOfAtoms( atomVec );

    vector<Residue>  residueVec = molVec.front().get_residueVec();

    //----------------------------------------------------------------

    vector<QVector3D>   normalVec;
    for( size_t i=0; i<residueVec.size(); i++ ){
        if( residueVec[i].get_index() >= 36 && residueVec[i].get_index() <= 42 ){
            Coord n;
            if( residueVec[i].get_index() %2 == 0 ){
                n = residueVec[i].get_O().get_coord() - residueVec[i].get_C().get_coord() ;
            }else{
                n = residueVec[i].get_C().get_coord() - residueVec[i].get_O().get_coord() ;
            }
            QVector3D normal( n.x, n.y, n.z );
            normalVec.push_back( normal );
        }
    }

    vector<vector<Atom> > sheetAtomVec ;
    vector<Atom>          atVec;
    for( size_t i=0; i<atomVec.size(); i++ ){
        if( atomVec[i].get_residueIndex() >= 36 && atomVec[i].get_residueIndex() <= 42 ){
            if( atomVec[i].get_name() == "CA" ){
                atVec.push_back( atomVec[i] );
            }
        }
    }
    sheetAtomVec.push_back( atVec );
    Coord center = centerOfAtoms( atVec );
    sheetCoords.clear();
    sheetNormals.clear();

    vector<QVector3D> boneVec;
    for( size_t i=0; i<atVec.size(); i++ ){
        Coord c = atVec[i].get_coord();
        QVector3D q((c.x-atomsCenter.x)/ maxRange,
                (c.y-atomsCenter.y)/ maxRange,
                (c.z-atomsCenter.z)/ maxRange);

        boneVec.push_back(q);
    }

    get_triangleVertex( boneVec, normalVec, 20, sheetVertexVec, sheetVertexNormalVec );

    glGenBuffers(1, &sheetVBO);
    glBindBuffer(GL_ARRAY_BUFFER, sheetVBO);
    glBufferData( GL_ARRAY_BUFFER,
            (sheetVertexVec.size()+sheetVertexNormalVec.size())*3*sizeof(float),
            NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sheetVertexVec.size()*3*sizeof(float),
            &sheetVertexVec[0]);
    glBufferSubData(GL_ARRAY_BUFFER, sheetVertexVec.size()*3*sizeof(float),
            sheetVertexNormalVec.size()*3*sizeof(float), &sheetVertexNormalVec[0]);

}

void
MolecularViewer::initializeTurn(){
    string vshaderFile = string("./shader/turn.vs");
    string fshaderFile = string("./shader/turn.fs");

    if( !turnProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !turnProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !turnProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !turnProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    //light
    turnProgram.setUniformValue("lightDirect", QVector3D(0.0, 0.0, 0.5));
    turnProgram.setUniformValue("ambientIntensity", float(0.1) );
    turnProgram.setUniformValue("diffuseIntensity", float(1.0) );
    turnProgram.setUniformValue("specularIntensity", float(0.5) );

    //----------------------------------------------------------------

    //-- atom
//	string fileName("ligand.mol2");
    string fileName("1A9U.pdb");
//	string fileName("4CYX.pdb");

    vector<Molecular> molVec = readMolecularFile(fileName);
    vector<Atom> atomVec = molVec.front().get_atomVec();
    noOfAtom = atomVec.size();
//    cout<<"noOfAtom: "<<noOfAtom<<endl;
    sphereModeAtomsVec = atomVec;

    float maxRange = atomsMaxRange( atomVec );

    atomsCenter = centerOfAtoms( atomVec );

    vector<Residue>  residueVec = molVec.front().get_residueVec();

    //----------------------------------------------------------------

    vector<Atom>          atVec;
    for( size_t i=0; i<atomVec.size(); i++ ){
        if( atomVec[i].get_residueIndex() >= 166 && atomVec[i].get_residueIndex() <= 179 ){
            if( atomVec[i].get_name() == "CA" ){
                atVec.push_back( atomVec[i] );
            }
        }
    }
    Coord center = centerOfAtoms( atVec );

    vector<QVector3D> boneVec;
    for( size_t i=0; i<atVec.size(); i++ ){
        Coord c = atVec[i].get_coord();
        QVector3D q((c.x-atomsCenter.x)/ maxRange,
                (c.y-atomsCenter.y)/ maxRange,
                (c.z-atomsCenter.z)/ maxRange);

        boneVec.push_back(q);
    }

    get_turnTrigVertex( boneVec, 20, turnVertexVec, turnVertexNormalVec );

    glGenBuffers(1, &turnVBO);
    glBindBuffer(GL_ARRAY_BUFFER, turnVBO);
    glBufferData( GL_ARRAY_BUFFER,
            (turnVertexVec.size()+turnVertexNormalVec.size())*3*sizeof(float),
            NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, turnVertexVec.size()*3*sizeof(float),
            &turnVertexVec[0]);
    glBufferSubData(GL_ARRAY_BUFFER, turnVertexVec.size()*3*sizeof(float),
            turnVertexNormalVec.size()*3*sizeof(float), &turnVertexNormalVec[0]);

}

void
MolecularViewer::initializeMesh(){
    string vshaderFile = string("./shader/mesh.vs");
    string fshaderFile = string("./shader/mesh.fs");

    if( !meshProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !meshProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !meshProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !meshProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    //light
    meshProgram.setUniformValue("lightDirect", QVector3D(0.0, 0.0, -0.5));
    meshProgram.setUniformValue("ambientIntensity", float(0.1) );
    meshProgram.setUniformValue("diffuseIntensity", float(1.0) );
    meshProgram.setUniformValue("specularIntensity", float(0.5) );

    //------------------------------------------------------

    //-- atom
    string fileName("ligand.mol2");

    vector<Molecular> molVec = readMolecularFile(fileName);
    vector<Atom> atomVec = molVec.front().get_atomVec();
    maxCoordRange = atomsMaxRange( atomVec );
    atomsCenter = centerOfAtoms( atomVec );
    vector<reducedSurface>  rsVec;
    float		probeRadius = 1.4;

    rsVec =     getReducedSurface( atomVec, probeRadius );

    vector<QVector3D> meshVec;
    meshNum = rsVec.size()*3*2;
    vector<float> coordVec;
    for( size_t i=0; i<rsVec.size(); i++ ){
        coordVec.push_back( (rsVec[i].get_at1().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at1().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at1().get_coord().z - atomsCenter.z)/maxCoordRange );

        coordVec.push_back( (rsVec[i].get_at2().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at2().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at2().get_coord().z - atomsCenter.z)/maxCoordRange );

        coordVec.push_back( (rsVec[i].get_at1().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at1().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at1().get_coord().z - atomsCenter.z)/maxCoordRange );

        coordVec.push_back( (rsVec[i].get_at3().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at3().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at3().get_coord().z - atomsCenter.z)/maxCoordRange );

        coordVec.push_back( (rsVec[i].get_at2().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at2().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at2().get_coord().z - atomsCenter.z)/maxCoordRange );

        coordVec.push_back( (rsVec[i].get_at3().get_coord().x - atomsCenter.x)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at3().get_coord().y - atomsCenter.y)/maxCoordRange );
        coordVec.push_back( (rsVec[i].get_at3().get_coord().z - atomsCenter.z)/maxCoordRange );
    }
    if(0){
        coordVec.clear();
        Atom a1 =  searchAtom( atomVec, 33 );
        Atom a2 =  searchAtom( atomVec, 34 );
        Atom a3 =  searchAtom( atomVec, 35 );
        for( size_t k=0; k<rsVec.size(); k++ ){
            vector<Coord> archVec = getArchPoints( rsVec[k].get_at1().get_coord(),
                    rsVec[k].get_at2().get_coord(), rsVec[k].get_probeCoord(),
                    probeRadius, 10 );
            for( size_t i=0; i<archVec.size()-1; i++ ){
                coordVec.push_back( (archVec[i].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i].z - atomsCenter.z)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].z - atomsCenter.z)/maxCoordRange );
            }

            archVec = getArchPoints( rsVec[k].get_at1().get_coord(),
                    rsVec[k].get_at3().get_coord(), rsVec[k].get_probeCoord(),
                    probeRadius, 10 );
            for( size_t i=0; i<archVec.size()-1; i++ ){
                coordVec.push_back( (archVec[i].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i].z - atomsCenter.z)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].z - atomsCenter.z)/maxCoordRange );
            }

            archVec = getArchPoints( rsVec[k].get_at2().get_coord(),
                    rsVec[k].get_at3().get_coord(), rsVec[k].get_probeCoord(),
                    probeRadius, 10 );
            for( size_t i=0; i<archVec.size()-1; i++ ){
                coordVec.push_back( (archVec[i].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i].z - atomsCenter.z)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].x - atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].y - atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (archVec[i+1].z - atomsCenter.z)/maxCoordRange );
            }
        }
        meshNum = coordVec.size();
    }
    if(0){
        //reduced surface arch
        coordVec.clear();
        Atom a1 =  searchAtom( atomVec, 52 );
        Atom a2 =  searchAtom( atomVec, 62 );
        Atom a3 =  searchAtom( atomVec, 16 );
        Atom a4 =  searchAtom( atomVec, 57 );

        Coord p1;
        getFixedProbe( a1, a2, a3, atomVec, 1.5, p1 );
        Coord p2;
        getFixedProbe( a1, a2, a4, atomVec, 1.5, p2 );

        vector<Coord> c = getReducedSurfaceArch(a1, a2, p1, p2, 10);
        for( size_t i=0; i<c.size()-1; i++ ){
            coordVec.push_back((c[i].x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c[i].y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c[i].z-atomsCenter.z)/maxCoordRange);
            coordVec.push_back((c[i+1].x-atomsCenter.x)/maxCoordRange);
            coordVec.push_back((c[i+1].y-atomsCenter.y)/maxCoordRange);
            coordVec.push_back((c[i+1].z-atomsCenter.z)/maxCoordRange);
        }
        meshNum = coordVec.size();
    }
    if(0){
        // SES arch
        coordVec.clear();
        Atom a1 =  searchAtom( atomVec, 29 );
        Atom a2 =  searchAtom( atomVec, 81 );
        Atom a3 =  searchAtom( atomVec, 75 );
        Atom a4 =  searchAtom( atomVec, 79 );

        Coord p1;
        getFixedProbe( a1, a2, a3, atomVec, probeRadius, p1 );
        Coord p2;
        getFixedProbe( a1, a2, a4, atomVec, probeRadius, p2 );

        vector<Coord> c1, c2;
        get_SES_arch(a1, a2, p1, p2, 10, c1, c2);
        for( size_t i=0; i<c1.size()-1; i++ ){
            coordVec.push_back( (c1[i].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (c1[i].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (c1[i].z-atomsCenter.z)/maxCoordRange );
            coordVec.push_back( (c1[i+1].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (c1[i+1].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (c1[i+1].z-atomsCenter.z)/maxCoordRange );
        }
        for( size_t i=0; i<c2.size()-1; i++ ){
            coordVec.push_back( (c2[i].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (c2[i].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (c2[i].z-atomsCenter.z)/maxCoordRange );
            coordVec.push_back( (c2[i+1].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (c2[i+1].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (c2[i+1].z-atomsCenter.z)/maxCoordRange );
        }
        meshNum = coordVec.size();

    }
    if(0){
        cout<<"---mesh---"<<endl;
        // SES arch
        coordVec.clear();
        Atom a1 =  searchAtom( atomVec, 29 );
        Atom a2 =  searchAtom( atomVec, 75 );
        Atom a3 =  searchAtom( atomVec, 81 );
        Coord c1 = a1.get_coord();
        Coord c2 = a2.get_coord();
        Coord c3 = a3.get_coord();

        Coord p1;
        getFixedProbe( a1, a2, a3, atomVec, probeRadius, p1 );
        int numOut = 10;
        vector<Coord> archP1 = getArchPoints( c1, c2, p1, probeRadius, numOut );
        vector<Coord> archP2 = getArchPoints( c1, c3, p1, probeRadius, numOut );
        vector<Coord> archP3 = getArchPoints( c2, c3, p1, probeRadius, numOut );

        if(0){

            c1.print();
            c2.print();
            c3.print();
            p1.print();
            cout<<"sphere"<<endl;
            cout<<"archP1:"<<endl;
            for( size_t i=0; i<archP1.size(); i++ ){
                archP1[i].print();
            }
            cout<<"archP2:"<<endl;
            for( size_t i=0; i<archP2.size(); i++ ){
                archP2[i].print();
            }
            cout<<"archP3:"<<endl;
            for( size_t i=0; i<archP3.size(); i++ ){
                archP3[i].print();
            }
            cout<<"dis12:"<<getCoordDis( archP1.back(), archP1.front() )<<endl;
            cout<<"dis13:"<<getCoordDis( archP2.back(), archP2.front() )<<endl;
            cout<<"dis23:"<<getCoordDis( archP3.back(), archP3.front() )<<endl;

        }

        for( size_t i=0; i<archP1.size()-1; i++ ){
            coordVec.push_back( (archP1[i].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP1[i].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP1[i].z-atomsCenter.z)/maxCoordRange );
            coordVec.push_back( (archP1[i+1].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP1[i+1].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP1[i+1].z-atomsCenter.z)/maxCoordRange );
        }
        for( size_t i=0; i<archP2.size()-1; i++ ){
            coordVec.push_back( (archP2[i].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP2[i].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP2[i].z-atomsCenter.z)/maxCoordRange );
            coordVec.push_back( (archP2[i+1].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP2[i+1].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP2[i+1].z-atomsCenter.z)/maxCoordRange );
        }
        for( size_t i=0; i<archP3.size()-1; i++ ){
            coordVec.push_back( (archP3[i].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP3[i].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP3[i].z-atomsCenter.z)/maxCoordRange );
            coordVec.push_back( (archP3[i+1].x-atomsCenter.x)/maxCoordRange );
            coordVec.push_back( (archP3[i+1].y-atomsCenter.y)/maxCoordRange );
            coordVec.push_back( (archP3[i+1].z-atomsCenter.z)/maxCoordRange );
        }

        if(1){
            size_t k=8;
            vector<Coord> line1 = get_sphere_arch( p1, probeRadius, archP1[k-1], archP2[k-1], k );
            vector<Coord> line2 = get_sphere_arch( p1, probeRadius, archP1[k], archP2[k], k+1 );
            cout<<"line1"<<endl;
            for( size_t i=0; i<line1.size()-1; i++ ){
                coordVec.push_back( (line1[i].x-atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (line1[i].y-atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (line1[i].z-atomsCenter.z)/maxCoordRange );
                coordVec.push_back( (line1[i+1].x-atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (line1[i+1].y-atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (line1[i+1].z-atomsCenter.z)/maxCoordRange );

            }
            cout<<"line2"<<endl;

            for( size_t i=0; i<line2.size()-1; i++ ){
                coordVec.push_back( (line2[i].x-atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (line2[i].y-atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (line2[i].z-atomsCenter.z)/maxCoordRange );
                coordVec.push_back( (line2[i+1].x-atomsCenter.x)/maxCoordRange );
                coordVec.push_back( (line2[i+1].y-atomsCenter.y)/maxCoordRange );
                coordVec.push_back( (line2[i+1].z-atomsCenter.z)/maxCoordRange );
            }
        }
        meshNum = coordVec.size();
    }
    glGenBuffers( 1, &meshVBO );
    glBindBuffer( GL_ARRAY_BUFFER, meshVBO );
    glBufferData( GL_ARRAY_BUFFER, coordVec.size()*sizeof(float),
            &coordVec[0], GL_STATIC_DRAW );
}

void
MolecularViewer::initializeHBond(){
    string vshaderFile = string("./shader/line.vs");
    string fshaderFile = string("./shader/dashline.fs");

    if( !dashlineProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            vshaderFile.c_str() ) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !dashlineProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !dashlineProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !dashlineProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    //light
    dashlineProgram.setUniformValue("lightDirect", QVector3D(0.0, 0.0, -0.5));
    dashlineProgram.setUniformValue("ambientIntensity", float(0.1) );
    dashlineProgram.setUniformValue("diffuseIntensity", float(1.0) );
    dashlineProgram.setUniformValue("specularIntensity", float(0.5) );

    //---atom
    string fileName("ligand.mol2");

    vector<Molecular> molVec = readMolecularFile(fileName);
    vector<Atom> atomVec = molVec.front().get_atomVec();
    vector<Bond> bdVec = molVec.front().get_bondVec();
    noOfAtom = atomVec.size();
    cout<<"noOfAtom: "<<noOfAtom<<endl;
    sphereModeAtomsVec = atomVec;

    float range = atomsMaxRange( atomVec );
    maxCoordRange = range;
    atomsCenter = centerOfAtoms( atomVec );

    vector<float> coordVec, colorVec;
    cout<<"bond vector:"<<bdVec.size()<<endl;

    string           goldParFile = "ff.dat";
    Sybyl            sybyl( goldParFile );

    vector<Atom> donorVec = sybyl.getHydrogenBondDonorAtoms(atomVec);
    vector<Atom> acceptorVec = sybyl.getHydrogenBondAcceptorAtoms(atomVec);

    float maxHbondDis = 3.0;

    for( size_t i=0; i<donorVec.size(); i++ ){
        for( size_t j=0; j<acceptorVec.size(); j++ ){
            float dis = getCoordDis( donorVec[i].get_coord(), acceptorVec[j].get_coord() );
            if( dis < maxHbondDis ){
                donorVec[i].print();
                acceptorVec[j].print();
                cout<<endl;
            }
        }
    }
}

//#####################################################################
void
MolecularViewer::drawAxis(){
    glEnable( GL_BLEND );
    glEnable(GL_LINE_SMOOTH);

    vector<QVector3D> axis;
    vector<QVector3D> axisColor;
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0, 0, 0.5) );
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0, 0.5, 0) );
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0.5, 0, 0) );

    axisColor.push_back( QVector3D(1, 0, 0) );
    axisColor.push_back( QVector3D(1, 0, 0) );
    axisColor.push_back( QVector3D(0, 1, 0) );
    axisColor.push_back( QVector3D(0, 1, 0) );
    axisColor.push_back( QVector3D(0, 0, 1) );
    axisColor.push_back( QVector3D(0, 0, 1) );

    QOpenGLBuffer vertexBuf;//( QOpenGLBuffer::VertexBuffer );
    vertexBuf.create();
    vertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vertexBuf.bind();
    vertexBuf.allocate( &axis[0], 6*sizeof(QVector3D) );

    QOpenGLBuffer colorBuf( QOpenGLBuffer::VertexBuffer );
    colorBuf.create();
    colorBuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    colorBuf.bind();
    colorBuf.allocate(&axisColor[0], 6*sizeof(QVector3D));

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));
    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    axisProgram.bind();
    vertexBuf.bind();
    axisProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    axisProgram.enableAttributeArray("position");
    axisProgram.setAttributeBuffer( "position", GL_FLOAT, 0, 3,  sizeof(QVector3D) );
    colorBuf.bind();
    axisProgram.enableAttributeArray("vcolor");
    axisProgram.setAttributeBuffer("vcolor", GL_FLOAT, 0, 3, sizeof(QVector3D) );
    glDrawArrays(GL_LINES, 0, axis.size() );

    axisProgram.release();
    vertexBuf.destroy();
    colorBuf.destroy();

}

void
MolecularViewer::drawTest( const QSize &windowSize ){
//    glEnableClientState( GL_VERTEX_ARRAY );
    string vshaderFile = string("./shader/test_vs.glsl");
    string fshaderFile = string("./shader/test_fs.glsl");

    QOpenGLShaderProgram testProgram;

    if( !testProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !testProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !testProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !testProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    vector<QVector3D> axis;
    vector<QVector3D> axisColor;
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0, 0, 0.5) );
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0, 0.5, 0) );
    axis.push_back( QVector3D(0, 0, 0) );
    axis.push_back( QVector3D(0.5, 0, 0) );

    axisColor.push_back( QVector3D(1, 0, 0) );
    axisColor.push_back( QVector3D(1, 0, 0) );
    axisColor.push_back( QVector3D(0, 1, 0) );
    axisColor.push_back( QVector3D(0, 1, 0) );
    axisColor.push_back( QVector3D(0, 0, 1) );
    axisColor.push_back( QVector3D(0, 0, 1) );

    QOpenGLBuffer vertexBuf;
    QOpenGLBuffer colorBuf;

    vertexBuf.create();
    colorBuf.create();

    colorBuf.bind();
    colorBuf.allocate(&axisColor[0], 6*sizeof(QVector3D));
    vertexBuf.bind();
    vertexBuf.allocate( &axis[0], 6*sizeof(QVector3D) );

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    glPointSize(20);

    testProgram.bind();
    testProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );

    int texcoordLocation = testProgram.attributeLocation("a_position");
    testProgram.enableAttributeArray(texcoordLocation);
    testProgram.setAttributeBuffer( texcoordLocation, GL_FLOAT, 0, 3, sizeof(QVector3D) );

    glDrawArrays( GL_POINTS, 0, axis.size() );

    testProgram.release();
    vertexBuf.destroy();
    colorBuf.destroy();
}

void
MolecularViewer::drawBindingSite( ){
    if( dotsData.empty() ){
        return;
    }

    glEnable( GL_BLEND );
    glDisable( GL_CULL_FACE );

    float x= dotsData[0];
    float y= dotsData[1];
    float z= dotsData[2];
    float r= dotsData[3];

    QOpenGLShaderProgram dotsProgram;
    string vshaderFile = string("./shader/bs_vs.glsl");
    string fshaderFile = string("./shader/bs_fs.glsl");

    if( !dotsProgram.addShaderFromSourceFile(QOpenGLShader::Vertex,
            QString::fromStdString(vshaderFile)) ){
        cerr<<"can not load vertex shader!!!"<<endl;
    }
    if( !dotsProgram.addShaderFromSourceFile(QOpenGLShader::Fragment,
            QString::fromStdString(fshaderFile)) ){
        cerr<<"can not load fragment shader!!!"<<endl;
    }
    if( !dotsProgram.link() ){
        cerr<<"can not link shader!!!"<<endl;
    }
    if( !dotsProgram.bind() ){
        cerr<<"can not bind shader!!!"<<endl;
    }

    Polyeder ply(3);
    vector<TriangleSurf> sphere = ply.get_surfTriangles();

    float           vertex[sphere.size()*3][3];
    for( size_t i=0; i<sphere.size(); i++ ){
        vertex[i*3][0] = sphere[i].v1[0];
        vertex[i*3][1] = sphere[i].v1[1];
        vertex[i*3][2] = sphere[i].v1[2];

        vertex[i*3+1][0] = sphere[i].v2[0];
        vertex[i*3+1][1] = sphere[i].v2[1];
        vertex[i*3+1][2] = sphere[i].v2[2];

        vertex[i*3+2][0] = sphere[i].v3[0];
        vertex[i*3+2][1] = sphere[i].v3[1];
        vertex[i*3+2][2] = sphere[i].v3[2];
    }

    vector<float> allVertex = ply.get_allVertex();
    vector<int> triangleVertexIndex = ply.get_triangleVertexIndex();

    QOpenGLBuffer vertexBuf, indexBuf;
    vertexBuf.create();
    vertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vertexBuf.bind();
    vertexBuf.allocate( &allVertex[0], allVertex.size()*sizeof(float) );

    QOpenGLBuffer indbuf( QOpenGLBuffer::IndexBuffer);
    indexBuf = indbuf;
    indexBuf.create();
    indexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw );
    indexBuf.bind();
    indexBuf.allocate( &triangleVertexIndex[0], triangleVertexIndex.size()*sizeof(int) );

    QMatrix4x4 matrix, matrix2;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    dotsProgram.bind();
    dotsProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    dotsProgram.setUniformValue( "mv_matrix",  matrixR );

    vertexBuf.bind();
    dotsProgram.enableAttributeArray( "position" );
    dotsProgram.setAttributeBuffer( "position", GL_FLOAT, 0, 3, 3*sizeof( float )  );
    indexBuf.bind();
    QVector3D center( (x-atomsCenter.x)/maxCoordRange,
                ( y-atomsCenter.y)/maxCoordRange,
                ( z-atomsCenter.z)/maxCoordRange );
    QVector3D color(1, 1, 0);

    dotsProgram.setUniformValue( "center", center );
    dotsProgram.setUniformValue( "color", color );
    dotsProgram.setUniformValue( "radius", r/maxCoordRange );

    glDrawElements( GL_TRIANGLES, sphere.size() * 3, GL_UNSIGNED_INT, 0 );
    dotsProgram.release();
}

void
MolecularViewer::drawLine(){
    glEnable( GL_BLEND );
    glEnable(GL_LINE_SMOOTH);

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    glLineWidth(2.0);
    lineProgram.bind();
    lineProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    lineVertexBuf.bind();
    lineProgram.enableAttributeArray("position");
    lineProgram.setAttributeBuffer( "position", GL_FLOAT, 0, 3, 3*sizeof(float) );
    lineColorBuf.bind();
    lineProgram.enableAttributeArray("vcolor");
    lineProgram.setAttributeBuffer("vcolor", GL_FLOAT, 0, 3, 3*sizeof(float) );
    glDrawArrays(GL_LINES, 0, lineNum );

    lineProgram.release();
}

void
MolecularViewer::drawStick(){
    glDisable(GL_BLEND);

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    stickProgram.bind();
    stickProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    stickProgram.setUniformValue( "mv_matrix",  matrixR );

    stickVertexBuf.bind();
    stickProgram.enableAttributeArray("position");
    stickProgram.setAttributeBuffer("position", GL_FLOAT, 0, 3, sizeof(QVector3D) );

    stickNormalBuf.bind();
    stickProgram.enableAttributeArray("normal");
    stickProgram.setAttributeBuffer( "normal", GL_FLOAT, 0, 3, 3*sizeof(float) );

    stickColorBuf.bind();
    stickProgram.enableAttributeArray("color");
    stickProgram.setAttributeBuffer( "color", GL_FLOAT, 0, 3, 3*sizeof(float) );

    glDrawArrays( GL_TRIANGLES, 0, stick_vertexVec.size() );
    stickProgram.release();
}

void
MolecularViewer::drawSphere(){
    glDisable(GL_DIFFUSE);
    glDisable( GL_BLEND );
    glDisable(GL_CULL_FACE);
    QMatrix4x4 matrix, matrix2;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix2.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));
    matrix2.scale( exp(m_distExp / 1200.0f));
    sphereProgram.bind();
    sphereProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    sphereProgram.setUniformValue( "mv_matrix",  matrix2 );

    sphereVertexBuf.bind();
    sphereProgram.enableAttributeArray( "position" );
    sphereProgram.setAttributeBuffer( "position", GL_FLOAT, 0, 3, 3*sizeof( float )  );
    sphereIndexBuf.bind();
    for( size_t i=0; i<sphereModeAtomsVec.size() ; i++ ){
//        for( size_t i=0; i<1 ; i++ ){
        Coord c =sphereModeAtomsVec[i].get_coord();
        QVector3D center( (c.x-atomsCenter.x)/maxCoordRange,
                (c.y-atomsCenter.y)/maxCoordRange,
                (c.z-atomsCenter.z)/maxCoordRange );
        QVector3D color;

        //color gradually changed according to the residue index
        if( !resSurfaceAtomVec.empty() ){
            vector<float> tempColor=getColorFromResidue(sphereModeAtomsVec[i].get_residueIndex());
            color.setX(tempColor[0]);
            color.setY(tempColor[1]);
            color.setZ(tempColor[2]);
        }else if( !surfaceAtomVec.empty() ){
            //surface mode all yellow
            color.setX(1);
            color.setY(1);
            color.setZ(0);
        }else{
            color.setX( sphereModeAtomsVec[i].get_color()[0] );
            color.setY( sphereModeAtomsVec[i].get_color()[1] );
            color.setZ( sphereModeAtomsVec[i].get_color()[2]);
        }

        sphereProgram.setUniformValue( "center", center );
        sphereProgram.setUniformValue( "color", color );
        sphereProgram.setUniformValue( "radius",
                sphereModeAtomsVec[i].get_vdwRadius() /( maxCoordRange) );

        glDrawElements( GL_TRIANGLES, surfs.size() * 3, GL_UNSIGNED_INT, 0 );
    }
    sphereProgram.release();
}

void
MolecularViewer::drawSheet(){
    glBindBuffer(GL_ARRAY_BUFFER, sheetVBO);
    QMatrix4x4 matrix;

    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    sheetProgram.bind();
    sheetProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    sheetProgram.setUniformValue( "mv_matrix",  matrixR );
    sheetProgram.setUniformValue("color", QVector3D(0.0f, 1.0f, 0.0f));

    sheetProgram.enableAttributeArray(0);
    sheetProgram.enableAttributeArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0 );
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0,
            (const void *)(sheetVertexVec.size()*3*sizeof(float)) );

    glDrawArrays( GL_TRIANGLES, 0, sheetVertexVec.size()*2 );

    sheetProgram.disableAttributeArray(0);
    sheetProgram.disableAttributeArray(1);
    sheetProgram.release();
}

void
MolecularViewer::drawTurn(){
    glBindBuffer(GL_ARRAY_BUFFER, turnVBO);
    QMatrix4x4 matrix;

    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    turnProgram.bind();
    turnProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    turnProgram.setUniformValue( "mv_matrix",  matrixR );
    turnProgram.setUniformValue("color", QVector3D(1.0f, 0.0f, 1.0f));

    turnProgram.enableAttributeArray(0);
    turnProgram.enableAttributeArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0 );
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0,
            (const void *)(turnVertexVec.size()*3*sizeof(float)) );

    glDrawArrays( GL_TRIANGLES, 0, turnVertexVec.size()*2 );

    turnProgram.disableAttributeArray(0);
    turnProgram.disableAttributeArray(1);
    turnProgram.release();
}

void
MolecularViewer::drawHelix(){
    glBindBuffer( GL_ARRAY_BUFFER, helixVBO );

    QMatrix4x4 matrix;

    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    helixProgram.bind();
    helixProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    helixProgram.setUniformValue( "mv_matrix",  matrixR );
    helixProgram.setUniformValue("color", QVector3D(1.0f, 0.0f, 1.0f));

    helixProgram.enableAttributeArray(0);
    helixProgram.enableAttributeArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0 );
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0,
            (const void *)(quadsCoord.size()*sizeof(float)) );

    glDrawArrays( GL_TRIANGLES, 0, 26664 );
    helixProgram.disableAttributeArray(0);
    helixProgram.disableAttributeArray(1);
    helixProgram.release();
}

void
MolecularViewer::drawMesh(){
    glBindBuffer( GL_ARRAY_BUFFER, meshVBO );

    QMatrix4x4 matrix;

    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    meshProgram.bind();
    meshProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    meshProgram.setUniformValue( "mv_matrix",  matrixR );
    meshProgram.enableAttributeArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0 );

    glDrawArrays(GL_LINES, 0, meshNum);
    meshProgram.disableAttributeArray(0);
    meshProgram.release();
}

void
MolecularViewer::drawSurface(){
    glDisable(GL_BLEND);

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));
    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    surfaceProgram.bind();
    surfaceProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    surfaceProgram.setUniformValue( "mv_matrix",  matrixR );

    surfaceVertexBuf.bind();
    surfaceProgram.enableAttributeArray("position");
    surfaceProgram.setAttributeBuffer("position",GL_FLOAT,0,3,3*sizeof(float));

    surfaceNormalBuf.bind();
    surfaceProgram.enableAttributeArray("normal");
    surfaceProgram.setAttributeBuffer("normal",GL_FLOAT,0,3,3*sizeof(float));

    glDrawArrays( GL_TRIANGLES, 0, surfaceTrigNum*3 );
    surfaceProgram.release();
}

void
MolecularViewer::drawResSurface(){
    glDisable(GL_BLEND);

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));
    QMatrix4x4 matrixR;
    matrixR.rotate(rotation);
    matrixR.scale( exp(m_distExp / 1200.0f));

    resSurfaceProgram.bind();
    resSurfaceProgram.setUniformValue( "mvp_matrix",  orthographic * matrix );
    resSurfaceProgram.setUniformValue( "mv_matrix",  matrixR );

    resSurfaceVertexBuf.bind();
    resSurfaceProgram.enableAttributeArray("position");
    resSurfaceProgram.setAttributeBuffer("position",GL_FLOAT,0,3,3*sizeof(float));

    resSurfaceNormalBuf.bind();
    resSurfaceProgram.enableAttributeArray("normal");
    resSurfaceProgram.setAttributeBuffer("normal",GL_FLOAT,0,3,3*sizeof(float));

    resSurfaceColorBuf.bind();
    resSurfaceProgram.enableAttributeArray("color");
    resSurfaceProgram.setAttributeBuffer("color",GL_FLOAT,0,3,3*sizeof(float));

    glDrawArrays( GL_TRIANGLES, 0, resSurfaceTrigNum*3 );
    resSurfaceProgram.release();
}

void
MolecularViewer::labelAtomIndex(){
    if(molVec.empty()){
        return;
    }
    glDisable(GL_BLEND);

    vector<Atom> atomVec;
    vector<Bond> bdVec;

    for( size_t i=0; i<molVec.size(); i++ ){
        vector<Atom> atVec=molVec[i].get_atomVec();
        for( size_t j=0; j<atVec.size(); j++ ){
            atomVec.push_back( atVec[j] );
        }
        vector<Bond> bd=molVec[i].get_bondVec();
        for( size_t j=0; j<bd.size(); j++ ){
            bdVec.push_back( bd[j] );
        }
    }

    QMatrix4x4 matrix;
    matrix.translate(translation.x(),
                     -translation.y(), 0.0);

    matrix.scale( exp(m_distExp / 3200.0f));
    cout<<"dist:"<<m_distExp<<endl;

    QMatrix4x4 matrixT;
    matrixT.translate(translation.x(),
                      -translation.y(), 0.0);
    matrixT.rotate(rotation);
    matrixT.scale( exp(m_distExp / 3200.0f));

    QMatrix4x4 view;
    view.lookAt(QVector3D(0,1,1), QVector3D(0,0,0), QVector3D(0,1,0));

    QPainter p(this);
    p.setWorldTransform((window_normalised_matrix).toTransform());
    p.setTransform( matrix.toTransform(), true);
    p.setPen(QPen(Qt::white));

    for( size_t i=0; i<atomVec.size(); i+=3 ){
        QStaticText text(atomVec[i].get_name().c_str());
        Coord c=atomVec[i].get_coord();
        float x= ((c.x-atomsCenter.x)/maxCoordRange);
        float y=-((c.y-atomsCenter.y)/maxCoordRange);
        float z=((c.z-atomsCenter.z)/maxCoordRange);
        QVector4D vec(x,y,z,1);
        QVector4D vecNew = orthographic*matrixT*vec;
        p.drawStaticText(vecNew.x()*500, vecNew.y()*500, text);
//        p.drawStaticText(x*500, y*500, text)
        cout<<vecNew.x()*500<<" "<<vecNew.y()*500<<" "<<text.text().toStdString()<<endl;
        cout<<x*500<<" "<<y*500<<endl<<endl;
    }

    if(1){
        QStaticText test("hello");
        p.drawStaticText(100,100,test);
    }
    //-------------------------------------------------------------
}

QPointF MolecularViewer::pixelPosToViewPos(const QPointF& p) {
    return QPointF(2.0 * float(p.x()) / width() - 1.0,
                   1.0 - 2.0 * float(p.y()) / height());
}

void MolecularViewer::mousePressEvent(QMouseEvent *event)
{

    if (event->buttons() & Qt::LeftButton) {
        m_trackBalls->push(pixelPosToViewPos(event->localPos()), QQuaternion() );
        rotation = m_trackBalls->rotation();
    }else if (event->buttons() & Qt::RightButton) {
        m_trackMover->push(pixelPosToViewPos(event->localPos()) );
    }
}

void MolecularViewer::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton) {
        m_trackBalls->move(pixelPosToViewPos(event->localPos() ), QQuaternion() );
        rotation = m_trackBalls->rotation();
        update();
    }else if (event->buttons() & Qt::RightButton) {

        m_trackMover->move(pixelPosToViewPos(event->localPos()) );
        translation =m_trackMover->translation();
        update();
    }else{
        m_trackBalls->release(pixelPosToViewPos(event->localPos()), QQuaternion() );
    }
}

void MolecularViewer::mouseReleaseEvent( QMouseEvent *event ){}

void MolecularViewer::mouseDoubleClickEvent(QMouseEvent *event){
    float clickPosX=(event->x()/float(width())-0.5)*2;
    float clickPosY=-(event->y()/float(height())-0.5)*2;

    vector<Atom>	atomVec;
    for( size_t i=0; i<molVec.size(); i++ ){
        vector<Atom> atVec=molVec[i].get_atomVec();
        for( size_t j=0; j<atVec.size(); j++ ){
            atomVec.push_back( atVec[j] );
        }
    }

    vector<QVector3D>   vertexVec;

    QMatrix4x4 matrix;
    matrix.translate(translation.x(), translation.y(), 0.0);
    matrix.rotate(rotation);
    matrix.scale( exp(m_distExp / 1200.0f));

    size_t atSize=atomVec.size();
    for( size_t i=0; i<atSize; i++ ){
    	Coord	c = atomVec[i].get_coord();
        float   x=(c.x-atomsCenter.x)/maxCoordRange;
        float   y=(c.y-atomsCenter.y)/maxCoordRange;
        float   z=(c.z-atomsCenter.z)/maxCoordRange;

        QVector3D v(x,y,z);
        vertexVec.push_back(v);
        QVector3D proj=orthographic*matrix*v;

        float clickToNodeDist=(clickPosX-proj.x())*(clickPosX-proj.x())+
                (clickPosY-proj.y())*(clickPosY-proj.y());
        if(clickToNodeDist<0.01){
        	atomsCenter=c;
        	break;
        }

        QPointF currPos = QPointF(0,0);
        trackMover *tm=new trackMover();
        m_trackMover=tm;
        translation = QVector2D(0,0);
    }
    updateLine( molVec );
    updateStick( molVec );
    updateSphere( molVec );
    updateSurface( molVec );
    updateResSurface( molVec );

    update();
}

void MolecularViewer::wheelEvent(QWheelEvent *event){
    m_distExp -= event->delta();
    update();
    event->accept();
}

vector<float>
MolecularViewer::getColorFromResidue(int residueIndex){
    vector<float> color(3,0);
    color[0] = (sin(float(residueIndex*4)*3.14/180.0))/2.0+0.5;
    color[1] = (cos(float(residueIndex*5)*3.14/180.0))/2.0+0.5;
    color[2] = (sin(float(residueIndex*4)*3.14/180.0))/2.0+0.5;
    return color;
}
