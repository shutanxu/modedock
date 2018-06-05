/*
 * window.h
 *
 *  Created on: Aug 25, 2014
 *      Author: stan
 */

#ifndef WINDOW_H_
#define WINDOW_H_

#include <QMainWindow>
#include<QGridLayout>
#include <QtWidgets>

#include"../view/molViewer.h"
#include"rmsdDialog.h"
#include"molTreeViewer.h"
#include"dockConfDialog.h"
#include"../dock/deDock.h"

using namespace std;

class MainWindow:public QMainWindow{
	Q_OBJECT
public:
	MainWindow(int argc, char** argv);
	~MainWindow();

private slots:
    void          updateViewMols();
	void          openFile();
	void          openFile( const QString str );
	void          saveFile();
    void          clear();
	void          viewAxis();
	void          viewStick();
	void          viewSphere();
	void          viewSurface();
    void          viewResSurface();
    void          computeLigandRMSD();
    void          computeProteinRMSD();
	void          dockConfig();
	void          dockRun();
    void          updateDockResultLigs( );
    void		  updateMol(vector<Molecular> molVec);
signals:
    void  	emitRun();

private:
    bool  loadFile(const QString &fileName);
    bool  okToContinue();
    void  createActions();
    void  createMenus();
    void  updateWindow();

    LigandRMSDDialog*		ligandRmsd;
    ProteinRMSDDialog*		proteinRmsd;
    DeDock*					deDock;

	string             	molFileName;
	vector<Molecular>  	molVec;
	MolecularViewer*	molViewer;
	MolTreeViewer*		molTreeViewer;
	QTextEdit*			termEdit;

    QMenu           *fileMenu;
    QAction         *openAction;
    QAction         *saveAction;
    QAction         *clearAction;

    QMenu           *editMenu;

    QMenu           *viewMenu;
    QAction         *viewAxisAction;
    QAction         *viewStickAction;
    QAction         *viewSphereAction;
    QAction         *viewSurfaceAction;
    QAction         *viewResSurfaceAction;

    QMenu           *toolsMenu;
    QAction         *ligandRmsdAction;
    QAction         *proteinRmsdAction;

    QMenu           *dockMenu;
    QAction         *dockAction;
    ConfigDialog    *confDlg;

    QMenu *			helpMenu;

    QAction*		showLineAction;
    QAction*		showStickAction;
    QAction*		showBallStickAction;
    QAction*		showVdwBallAction;
    QAction*		showSecStrucAction;
    QAction*		addHydrogenAction;
    QAction*		computeVdwAction;
};

#endif /* WINDOW_H_ */
