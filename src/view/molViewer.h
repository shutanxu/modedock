/*
 * molViewer.h
 *
 *  Created on: Aug 25, 2014
 *      Author: stan
 */

#ifndef MOLVIEWER_H_
#define MOLVIEWER_H_

#define GL_GLEXT_PROTOTYPES

#include<QGLWidget>
#include<QOpenGLWidget>
#include<QOpenGLFunctions>
#include<QOpenGLShaderProgram>
#include<GL/gl.h>
#include<string>
#include<fstream>
#include<iostream>

#include"../mol/molecular.h"
#include"trackball.h"
#include"polyeder.h"
#include"b_spline.h"
#include"trackMover.h"
#include"surface.h"
#include"../forcefield/sybyl.h"
#include"stick.h"

using namespace std;

class MolecularViewer:public QOpenGLWidget, protected QOpenGLFunctions{
    Q_OBJECT
public:
    MolecularViewer( QWidget *parent = 0 );
    ~MolecularViewer();
    void    setMolVec( const vector<Molecular> mol );
    void    clear(){ molVec.clear();  }
    void    updateAxis(const bool b){ if(b){showAxis=true;}else{showAxis=false;}}
    void    updateLine( const vector<Molecular> mol );
    void    updateStick( const vector<Molecular> mol );
    void    updateSphere( const vector<Molecular> mol );
    void    updateSurface( const vector<Molecular> mol );
    void    updateResSurface( const vector<Molecular> mol );
    void	viewLine( bool f){ showLine = f; }
    void	viewStick( bool f){ showStick = f; }
    void	viewSphere( bool f){ showSphere = f; }
    void	viewSurface( bool f){ showSurface = f; }
    void	viewResSurface( bool f){ showResSurface = f; }
protected:
    void initializeGL();
    void resizeGL( int width, int height );
    void paintGL();
    void mousePressEvent( QMouseEvent *event );
    void mouseMoveEvent( QMouseEvent *event );
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseDoubleClickEvent( QMouseEvent *event );
    void wheelEvent( QWheelEvent *e );

private slots:
    void sDrawDots(vector<float> vec){ dotsData=vec;drawBindingSite(); /*update();*/ };
    void sClearDots(){ dotsData.clear(); }
private:
    QPointF pixelPosToViewPos(const QPointF& p) ;
    void    initializeAxis();
//    void    initializeDots();
    void    initializeLine();
    void    initializeStick();
    void    initializeSphere();

    void    initializeHelix();
    void    initializeSheet();
    void    initializeTurn();
    void    initializeMesh();
    void    initializeSurface();
    void    initializeResSurface();
    void    initializeHBond();
    void    drawTest( const QSize &windowSize );
    void    drawAxis();
    void    drawBindingSite();
    void    drawLine();
    void    drawStick();
    void    drawSphere();
    void    drawHelix();
    void    drawSheet();
    void    drawTurn();
    void    drawMesh();
    void    drawSurface();
    void    drawResSurface();
    void    labelAtomIndex();

    vector<float>    getColorFromResidue(int residueIndex);

    vector<Molecular>    molVec;
    bool                 showAxis;
    vector<Bond>         stickBondVec;
    vector<Atom>         sphereModeAtomsVec;
    vector<Atom>         surfaceAtomVec;
    vector<Atom>         resSurfaceAtomVec;
    trackMover           *m_trackMover;
    QOpenGLShaderProgram     program;
    QOpenGLShaderProgram
                         axisProgram,
                         lineProgram,
                         stickProgram,
                         sphereProgram,
                         helixProgram,
                         sheetProgram,
                         turnProgram,
                         meshProgram,
                         surfaceProgram,
                         resSurfaceProgram,
                         dashlineProgram;  
    QMatrix4x4           orthographic;

    QMatrix4x4           window_normalised_matrix;
    QMatrix4x4           window_painter_matrix;
    QMatrix4x4           projection;

    QVector2D            mousePressPosition;
    QQuaternion          rotation;
    QVector2D            translation;
    int                  m_distExp;
    TrackBall*           m_trackBalls;

    bool				 showLine, showStick, showSphere, showSurface, showResSurface;

    QOpenGLBuffer  lineVertexBuf, lineColorBuf;
    QOpenGLBuffer  stickVertexBuf, stickNormalBuf, stickColorBuf ;
    vector<QVector3D>    stick_vertexVec, stick_normalVec;
    QOpenGLBuffer sphereVertexBuf,sphereIndexBuf;
    QOpenGLBuffer surfaceVertexBuf,surfaceNormalBuf,surfaceColorBuf ;
    QOpenGLBuffer resSurfaceVertexBuf,resSurfaceNormalBuf,resSurfaceColorBuf ;

    GLuint               lineVBO,
                         stickVBO,
                         sphereVBO,
                         sphereIBO,
                         helixVBO,
                         sheetVBO,
                         turnVBO,
                         meshVBO,
                         surfaceVBO,
                         resSurfaceVBO;
    vector<TriangleSurf> surfs;
    int                  noOfAtom;

    float                maxCoordRange;
    Coord                atomsCenter;

    vector<float>        quadsCoord;
    vector<float>        quadsNormal;

    vector<float>        sheetCoords;
    vector<float>        sheetNormals;
    vector<QVector3D>    sheetVertexVec, sheetVertexNormalVec;
    vector<QVector3D>    turnVertexVec, turnVertexNormalVec;
    int                  lineNum;
    int					 meshNum;
    int                  surfaceTrigNum;
    int                  resSurfaceTrigNum;

    vector<float>        dotsData;
};


#endif /* MOLVIEWER_H_ */
