/*
 * window.cpp
 *
 *  Created on: Aug 25, 2014
 *      Author: stan
 */

#define GLEW_STATIC
#define GL_GLEXT_PROTOTYPES

#include <QGridLayout>
#include <QSplitter>
#include<QGLFunctions>
#include"window.h"

using namespace std;

MainWindow::MainWindow(int argc, char** argv){

    GLuint *axisVBO;
//    glGenBuffers(1, axisVBO);

	setWindowTitle(tr("mol"));
    resize(700, 700);
	QWidget *centralWidget = new QWidget;

	molViewer = new MolecularViewer();
	QGridLayout *layout = new QGridLayout;
	layout->addWidget( molViewer, 0, 0);
	centralWidget->setLayout(layout);
	setCentralWidget(centralWidget);

	molTreeViewer	= new MolTreeViewer();
    molVec.clear();

    termEdit = new QTextEdit;
	termEdit->setReadOnly(true);
//	termEdit->append(tr("xushutan, Jilin University, xushutan@gmail.com"));

	createActions();
	createMenus();

	if(1){
        QSplitter	*mainSplitter = new QSplitter( Qt::Vertical );
        QSplitter       *upSplitter =  new QSplitter( Qt::Horizontal );
        upSplitter->addWidget( molTreeViewer );
        upSplitter->addWidget( molViewer );
        upSplitter->setStretchFactor(0, 8);
        upSplitter->setStretchFactor(1, 20);
        mainSplitter->addWidget(upSplitter);
        mainSplitter->addWidget(termEdit);
        mainSplitter->setStretchFactor(0, 20);
        mainSplitter->setStretchFactor(1, 2);

        setCentralWidget( mainSplitter );
	}

    QDesktopWidget *desktop = QApplication::desktop();
    move( (desktop->width() - this->width())/2, (desktop->height()-this->height())/2 );
}

MainWindow::~MainWindow(){

}

void
MainWindow::updateWindow(){

//    QSplitter		*mainSplitter = new QSplitter( Qt::Vertical );
//	QSplitter       *upSplitter =  new QSplitter( Qt::Horizontal );
//	upSplitter->addWidget( molTreeViewer );
//    upSplitter->addWidget( molViewer );
//	upSplitter->setStretchFactor(0, 8);
//	upSplitter->setStretchFactor(1, 20);
//	mainSplitter->addWidget(upSplitter);
//	mainSplitter->addWidget(termEdit);
//	mainSplitter->setStretchFactor(0, 20);
//	mainSplitter->setStretchFactor(1, 2);

//	setCentralWidget( mainSplitter );
}

void
MainWindow::createActions(){
	openAction = new QAction( tr("&Open..."), this );
    openAction->setStatusTip(tr("Open an existing molecular file"));
	connect( openAction, SIGNAL( triggered() ), this, SLOT( openFile()  ) );

    clearAction = new QAction( tr("&Clear"), this );
    clearAction->setStatusTip(tr("Clear all the moleculars"));
    connect( clearAction, SIGNAL( triggered() ), this, SLOT( clear() ) );

	viewAxisAction=new QAction(tr("&Axis"), this);
	viewAxisAction->setStatusTip(tr("Show Axis"));
	viewAxisAction->setCheckable(true);
	connect(viewAxisAction, SIGNAL(triggered()), this, SLOT(viewAxis()));

	viewStickAction=new QAction(tr("&Stick"), this);
	viewStickAction->setStatusTip( tr("Stick Model") );
	viewStickAction->setCheckable(true);
    viewStickAction->setChecked(false);
	connect(viewStickAction, SIGNAL(triggered()), this, SLOT(viewStick()));

	viewSphereAction=new QAction(tr("&Sphere"),this);
	viewSphereAction->setStatusTip(tr("Sphere Model"));
	viewSphereAction->setCheckable(true);
	connect(viewSphereAction, SIGNAL(triggered()), this, SLOT(viewSphere()));

	viewSurfaceAction=new QAction(tr("&Surface"),this);
	viewSurfaceAction->setStatusTip(tr("Surface Model"));
	viewSurfaceAction->setCheckable(true);
	connect(viewSurfaceAction, SIGNAL(triggered()), this, SLOT(viewSurface()));

    viewResSurfaceAction=new QAction(tr("&ResSurface"),this);
    viewResSurfaceAction->setStatusTip(tr("ResSurface Model"));
    viewResSurfaceAction->setCheckable(true);
    connect(viewResSurfaceAction, SIGNAL(triggered()), this, SLOT(viewResSurface()));


    ligandRmsdAction = new QAction( tr("&Ligand RMSD"),this);
    ligandRmsdAction->setStatusTip(tr("RMSD"));
    connect(ligandRmsdAction, SIGNAL(triggered()), this, SLOT(computeLigandRMSD()));

    proteinRmsdAction = new QAction( tr("&Protein RMSD"),this);
    proteinRmsdAction->setStatusTip(tr("RMSD"));
    connect(proteinRmsdAction, SIGNAL(triggered()), this, SLOT(computeProteinRMSD()));
	//---------------------------------------------------------------------------

	dockAction = new QAction( tr("&Dock"), this );
	dockAction->setStatusTip(tr("Config docking parameters."));
	connect( dockAction, SIGNAL( triggered() ), this, SLOT(dockConfig() ) );

	//----------------------- config dock --------------------------------------------------------
	confDlg = new ConfigDialog;
	connect( confDlg, SIGNAL( loadP(QString) ), this, SLOT(openFile(QString)) );
	connect( confDlg, SIGNAL( loadB(vector<float>) ), molViewer, SLOT( sDrawDots(vector<float>) ) );
    connect( confDlg, SIGNAL( loadL(QString) ), this, SLOT(openFile(QString)) );
    connect( confDlg, SIGNAL( dockRun() ), this, SLOT(dockRun()) );
    connect( confDlg, SIGNAL( dockRun() ), molViewer, SLOT( sClearDots() ) );
    //-- textEdit
//    connect( deDock, SIGNAL(emitString(QString)), termEdit, SLOT(append(QString)) );
//    connect( this, SIGNAL(emitRun()), deDock, SLOT(run()) );
    connect( molTreeViewer,  SIGNAL(emitMolChanged( )),
             this, SLOT(updateViewMols( ) ) );
}

void
MainWindow::createMenus(){

	fileMenu =  menuBar()->addMenu( tr("&File") );
	fileMenu->addAction( openAction );
    fileMenu->addAction( clearAction );

	editMenu =  menuBar()->addMenu( tr("&Edit") );

	viewMenu = menuBar()->addMenu( tr("&View") );
	viewMenu->addAction( viewAxisAction );
	viewMenu->addAction( viewStickAction );
	viewMenu->addAction( viewSphereAction );
	viewMenu->addAction( viewSurfaceAction );
    viewMenu->addAction( viewResSurfaceAction );

	toolsMenu =  menuBar()->addMenu( tr("&Tools") );
	toolsMenu->addAction(ligandRmsdAction);
	toolsMenu->addAction(proteinRmsdAction);

	dockMenu = menuBar()->addMenu( tr("&Dock") );
	dockMenu->addAction( dockAction );
}

void
MainWindow::updateViewMols(){

    molViewer->updateLine( molTreeViewer->getVisibleMols() );
    molViewer->update();

    if( viewStickAction->isChecked() ){
        vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
        molViewer->updateStick(mol);
        molViewer->update();
    }
    if( viewSphereAction->isChecked() ){
        vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
        molViewer->updateSphere(mol);
        molViewer->update();
    }
    if( viewSurfaceAction->isChecked() ){
        vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
        molViewer->updateSurface(mol);
        molViewer->update();
    }
    if( viewResSurfaceAction->isChecked() ){
        vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
        molViewer->updateResSurface(mol);
        molViewer->update();
    }
}

void
MainWindow::updateMol(vector<Molecular> molVec){
    if( ! molVec.empty() ){

        molTreeViewer->setMolVec( molVec );
        molTreeViewer->update();

        molViewer->setMolVec( molVec );
        molViewer->updateLine(molVec );
        molViewer->update();

        updateWindow();
	}
}

void
MainWindow::updateDockResultLigs( ){
    vector<Molecular> mol = deDock->getFinalPop();
    molVec.insert( molVec.end(), mol.begin(), mol.end() );

    molViewer->setMolVec( molVec );
    molViewer->updateLine(molVec );
    molViewer->update();

    molTreeViewer->setMolVec( molVec );
    molTreeViewer->update();

    updateWindow();
}

bool MainWindow::loadFile(const QString &fileName )
{

	cout<<"file:"<<fileName.toStdString()<<endl;
	vector<Molecular>				mol;
    mol = readMolecularFile( fileName.toStdString().c_str() );
	termEdit->append(tr((string("read molecular: ")).c_str())
			+fileName);

    if( ! mol.empty() ){

		for( size_t i=0; i<mol.size(); i++ ){
			molVec.push_back( mol[i] );
		}

        molTreeViewer->setMolVec( molVec );
        molTreeViewer->update();

        molViewer->setMolVec( molVec );

        molViewer->updateLine(molVec );

        molViewer->update();

        updateWindow();
	}
    return true;
}

void
MainWindow::openFile(){

	QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open Molecular File"), ".",
			tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
	if (!fileName.isEmpty()){
		loadFile(fileName);
	}
}

void
MainWindow::openFile( const QString str ){

	QString fileName = str;
	if (!fileName.isEmpty()){
		loadFile(fileName);
	}
}

void
MainWindow::clear(){
    viewAxisAction->setChecked(false);
    viewStickAction->setChecked(false);
    viewSphereAction->setChecked(false);
    viewSurfaceAction->setChecked(false);
    viewResSurfaceAction->setChecked(false);

    molVec.clear();

    vector<Molecular> mol;
    molViewer->updateStick( mol );
    molViewer->updateSphere( mol );
    molViewer->updateSurface( mol );
    molViewer->updateResSurface( mol );

    molViewer->clear();
    molViewer->update();
    molTreeViewer->clear();
    molTreeViewer->update();
}

void
MainWindow::saveFile(){

}

void
MainWindow::viewAxis(){
		if( !viewAxisAction->isChecked() ){
			viewAxisAction->setChecked(false);
			molViewer->updateAxis(false);
			molViewer->update();
		}else{
			viewAxisAction->setChecked(true);
			molViewer->updateAxis(true);
			molViewer->update();
		}
}

void
MainWindow::viewStick(){
	if( !molVec.empty() ){
        if( viewStickAction->isChecked() ){
            viewStickAction->setChecked(true);
            vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
            molViewer->viewStick(true);
            molViewer->updateStick(mol);
			molViewer->update();
		}else{
            viewStickAction->setChecked(false);
			vector<Molecular> mol;
			molViewer->viewStick(false);
			molViewer->updateStick( mol );
			molViewer->update();
		}
	}
}

void
MainWindow::viewSphere(){
	if( !molVec.empty() ){
        if( viewSphereAction->isChecked() ){
			viewSphereAction->setChecked(true);
            vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
            molViewer->viewSphere(true);
            molViewer->updateSphere(mol);
//			molViewer->updateSphere( molVec );
			molViewer->update();
		}else{
			viewSphereAction->setChecked(false);
			vector<Molecular> mol;
			molViewer->viewSphere(false);
			molViewer->updateSphere( mol );
			molViewer->update();
		}
	}
}

void
MainWindow::viewSurface(){
	if( !molVec.empty() ){
        if( viewSurfaceAction->isChecked() ){
			viewSurfaceAction->setChecked(true);
            vector<Molecular> mol = molTreeViewer->getVisibleMols() ;

            molViewer->viewSurface(true);
            molViewer->updateSurface(mol);
			molViewer->update();

		}else{
			viewSurfaceAction->setChecked(false);
			vector<Molecular> mol;

			molViewer->viewSurface(false);
			molViewer->updateSurface( mol );
			molViewer->update();

		}
	}
}

void
MainWindow::viewResSurface(){
    if( !molVec.empty() ){
        if( viewResSurfaceAction->isChecked() ){
            viewResSurfaceAction->setChecked(true);
            vector<Molecular> mol = molTreeViewer->getVisibleMols() ;
            molViewer->viewResSurface(true);
            molViewer->updateResSurface(mol);
            molViewer->update();
        }else{
            viewResSurfaceAction->setChecked(false);
            vector<Molecular> mol;
            molViewer->viewResSurface(false);
            molViewer->updateResSurface( mol );
            molViewer->update();
        }
    }
}

void
MainWindow::computeLigandRMSD(){
    ligandRmsd = new LigandRMSDDialog;
    ligandRmsd->move(200, 200);
    ligandRmsd->show();
}

void
MainWindow::computeProteinRMSD(){
    proteinRmsd = new ProteinRMSDDialog;
    proteinRmsd->move(200, 200);
    proteinRmsd->show();

    connect( proteinRmsd, SIGNAL(signalOlyProtein(vector<Molecular>)),
    		this, SLOT(updateMol(vector<Molecular>)) );
}

void
MainWindow::dockConfig(){
	confDlg->open(); // non-block call

}

void
MainWindow::dockRun(){
	cout<<"---docking---"<<endl;
	deDock = new DeDock(string(""));
    connect( deDock, SIGNAL(emitString(QString)), termEdit, SLOT(append(QString)) );
    connect( deDock, SIGNAL(finishRun()), this, SLOT(updateDockResultLigs()) );
	deDock->start();
}
