/*
 * molTreeWidget.h
 *
 *  Created on: May 31, 2014
 *      Author: stan
 */

#ifndef MOLTREEWIDGET_H_
#define MOLTREEWIDGET_H_

#include<QGLWidget>
#include<QtOpenGL/QGLFunctions>
#include<QtOpenGL/QGLShaderProgram>
#include<GL/gl.h>
#include<string>
#include<fstream>
#include<iostream>
#include<qtreewidget.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<cstring>
#include<sstream>
#include<stdlib.h>
#include<algorithm>
#include<cmath>
#include<cstdlib>
//#include<glut.h>
#include <QDialog>
#include <QtGui>
#include<time.h>

#include"../mol/coord.h"
#include"../mol/molecular.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
class MyTreeWidgetItem: public QTreeWidgetItem{
//	Q_OBJECT;

public:
    explicit MyTreeWidgetItem(int type = Type):QTreeWidgetItem(type){ Atom at; atom = at; }
    explicit MyTreeWidgetItem(QTreeWidget *view, int type = Type):QTreeWidgetItem( view, type ){}

	void     set_atom( const Atom& at ){ atom = at; }
	Atom     get_atom()const{ return atom; }
private:
	Atom     atom;
};

class MyTreeWidget : public QTreeWidget
{
	Q_OBJECT;

public:
	MyTreeWidget()
	{
		setContextMenuPolicy(Qt::CustomContextMenu);

		connect(this,
				SIGNAL(customContextMenuRequested(const QPoint&)),
				SLOT(onCustomContextMenuRequested(const QPoint&)));
	}

private slots:
void        inverse_lines();
void        inverse_sticks();
void        inverse_vdwBall();
void        inverse_stickBall();
void        onCustomContextMenuRequested(const QPoint& pos);
void        showContextMenu( QTreeWidgetItem* item, const QPoint& globalPos) ;

private:
QAction*          showLineAction;
QAction*          showStickAction;
QAction*          showStickBallAction;
QAction*          showVdwBallAction;
QAction*          labelAtom;
QTreeWidgetItem*  item;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MolTreeViewer: public QGLWidget{
	Q_OBJECT
public:
	MolTreeViewer(){};
	~MolTreeViewer(){ delete treeWidget; };

	MolTreeViewer( const vector<Molecular>& molVec, QWidget *parent=0 );
    void               setMolVec( const vector<Molecular>& atVec );
	void               addMol( );
	vector<Molecular>  getMolVec()const{ return molVec; }
	QTreeWidget*       get_treeWidget(){ return treeWidget; }
	void               updateChildItemState(  MyTreeWidgetItem* item );
	vector< vector< pair<Atom, Atom> > >  getVisibleLineModeAtomPairVec();			// first level is mol index
	vector< vector< pair<Atom, Atom> > >  getVisibleStickModeAtomPairVec();		    // first level is mol index
	vector< vector< Atom > >              getVisibleVDWballModeAtomPairVec();		// first level is mol index
	vector< vector< pair<Atom, Atom> > >  getVisibleStickBallModeAtomPairVec();		// first level is mol index
private slots:
//	void												showContextMenu(const QPoint &pos );
private:
	void         addChain( MyTreeWidgetItem *parent,
                           const size_t& moleIndex,
                           const Chain& ch );
	void         addSubChain( MyTreeWidgetItem *parent,
                              const size_t& moleIndex,
                              const Chain& ch,
                              const SubChain& sub );
	void         addResidue(  MyTreeWidgetItem *parent,
                              const size_t& moleIndex,
                              const Chain& ch,
                              const SubChain& sub ,
                              const Residue& res);
	void         addAtom(  MyTreeWidgetItem *parent,
                           const size_t& moleIndex,
                           const Chain& ch,
                           const SubChain& sub ,
                           const Residue& res,
                           const Atom& at );
	vector<Molecular>          molVec;
	MyTreeWidget*              treeWidget;
	vector<MyTreeWidgetItem*>  allWidgetItemVec;
};
*/


//----------------------------------------------------------------------------
class MolTreeViewer: public QTreeView {
    Q_OBJECT
public:
    MolTreeViewer();
	void setMolVec( const vector<Molecular> mol );
    void clear(){ molVec.clear(); standardModel->clear(); }
    vector<Molecular> getVisibleMols();
    QStandardItemModel *standardModel;
public slots:
    void molChanged(QStandardItem* it);
signals:
    void emitMolChanged( );
private:
    void createActions();
	vector<Molecular> molVec;
    vector<int> visibleMolIndex;
};
#endif /* MOLTREEWIDGET_H_ */
