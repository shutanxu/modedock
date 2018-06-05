/*
 * dockConfDialog.cpp
 *
 *  Created on: Jan 18, 2015
 *      Author: stan
 */

#include"dockConfDialog.h"
#include<fstream>

ProteinPage::ProteinPage(QWidget *parent)
    : QWidget(parent)
{
    QGroupBox *configGroup = new QGroupBox(tr("Protein configuration"));

    QLabel    *proteinLabel = new QLabel(tr("Protein:"));
    proteinEdit = new QLineEdit;
    QPushButton *proteinLoad = new QPushButton(tr("Load"));

    connect( proteinLoad, SIGNAL(clicked()), this, SLOT(loadFile()) );

    QHBoxLayout *serverLayout = new QHBoxLayout;
    serverLayout->addWidget(proteinLabel);
    serverLayout->addWidget(proteinEdit);
    serverLayout->addWidget(proteinLoad);

    QVBoxLayout *configLayout = new QVBoxLayout;
    configLayout->addLayout(serverLayout);
    configGroup->setLayout(configLayout);

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));
    QPushButton *nextButton = new QPushButton(tr("Next"));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(sNext()));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
//    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(configGroup);
    mainLayout->addStretch(120);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);
}

void
ProteinPage::loadFile(){
	QString file = QFileDialog::getOpenFileName(this,
			tr("Open Molecular File"), ".",
			tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
	if (!file.isEmpty()){
		proteinEdit->insert( file );
		fileName = file.toStdString();
		emit loadProtein( file );
	}
}

//-------------------------------------------------------------

BindingSitePage::BindingSitePage(QWidget *parent)
    : QWidget(parent)
{

    QGroupBox *packageGroup = new QGroupBox(tr("Binding Site"));

    QLabel *xLabel=new QLabel( tr("center x:") );
    xEdit=new QLineEdit;

    QLabel *yLabel=new QLabel( tr("y:") );
    yEdit=new QLineEdit;

    QLabel *zLabel=new QLabel( tr("z:") );
    zEdit=new QLineEdit;

    QHBoxLayout *centerLayout = new QHBoxLayout;
    centerLayout->addWidget(xLabel);
    centerLayout->addWidget(xEdit);
    centerLayout->addWidget(yLabel);
    centerLayout->addWidget(yEdit);
    centerLayout->addWidget(zLabel);
    centerLayout->addWidget(zEdit);

    QLabel *rLabel=new QLabel(tr("radius:"));
    rEdit=new QLineEdit;
    QHBoxLayout *rLayout = new QHBoxLayout;
    rLayout->addWidget(rLabel);
    rLayout->addWidget(rEdit);

    QPushButton *updateButton = new QPushButton(tr("view binding site"));
    connect(updateButton, SIGNAL(clicked()), this, SLOT(getData()) );

    QVBoxLayout *packageLayout = new QVBoxLayout;
    packageLayout->addLayout(centerLayout);
    packageLayout->addLayout(rLayout);
    packageGroup->setLayout(packageLayout);

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));
    QPushButton *nextButton = new QPushButton(tr("Next"));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(sNext()));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(packageGroup);
    mainLayout->addWidget(updateButton);
    mainLayout->addStretch(12);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);
}

void
BindingSitePage::getData(){
	x=xEdit->text().toFloat();
	y=yEdit->text().toFloat();
	z=zEdit->text().toFloat();
	r=rEdit->text().toFloat();

	vector<float> d;
	d.push_back(x);
	d.push_back(y);
	d.push_back(z);
	d.push_back(r);

	emit bindingSite(d);
}

//-------------------------------------------------------------------------------

LigandPage::LigandPage(QWidget *parent)
    : QWidget(parent)
{
    QGroupBox *configGroup = new QGroupBox(tr("Ligand configuration"));

    QLabel    *ligandLabel = new QLabel(tr("Init Ligand:"));
    ligandEdit = new QLineEdit;
    QPushButton *ligandLoad = new QPushButton(tr("Load"));
    connect( ligandLoad, SIGNAL(clicked()), this, SLOT(loadFile()) );

    QLabel *ligandLabelRef=new QLabel( tr("Reference Ligand:") );
    refLigandEdit=new QLineEdit;
    QPushButton *refLigandLoad = new QPushButton(tr("Load"));
    connect( refLigandLoad, SIGNAL(clicked()), this, SLOT(loadRefFile()) );

    QHBoxLayout *serverLayout = new QHBoxLayout;
    serverLayout->addWidget(ligandLabel);
    serverLayout->addWidget(ligandEdit);
    serverLayout->addWidget(ligandLoad);

    QHBoxLayout *refLayout = new QHBoxLayout;
    refLayout->addWidget(ligandLabelRef);
    refLayout->addWidget(refLigandEdit);
    refLayout->addWidget(refLigandLoad);

    QVBoxLayout *configLayout = new QVBoxLayout;
    configLayout->addLayout(serverLayout);
    configLayout->addLayout(refLayout);
    configGroup->setLayout(configLayout);

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));
    QPushButton *nextButton = new QPushButton(tr("Next"));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(sNext()));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(configGroup);
    mainLayout->addStretch(12);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);
}

void
LigandPage::loadFile(){
	QString file = QFileDialog::getOpenFileName(this,
			tr("Open Molecular File"), ".",
			tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
	if (!file.isEmpty()){
		ligandEdit->insert( file );
		fileName = file.toStdString();
        emit loadLigand( file );
	}
}

void
LigandPage::loadRefFile(){
    QString file = QFileDialog::getOpenFileName(this,
            tr("Open Molecular File"), ".",
            tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
    if (!file.isEmpty()){
        refLigandEdit->insert( file );
        refFileName = file.toStdString();
        emit loadRefLigand( file );
    }
}

//--------------------------------------------------------------------------------

ScorePage::ScorePage(QWidget *parent)
    : QWidget(parent)
{
    QGroupBox *packagesGroup = new QGroupBox(tr("Scoring function"));

    QLabel *ffLabel=new QLabel(tr("Force field:"));
    QComboBox *ffCombo=new QComboBox;
//    ffCombo->addItem( tr("Amber") );
//    ffCombo->addItem( tr("Charmm") );
    ffCombo->addItem( tr("Sybyl") );

    QHBoxLayout *ffLayout=new QHBoxLayout;
    ffLayout->addWidget(ffLabel);
    ffLayout->addWidget(ffCombo);

//    QPushButton *startQueryButton = new QPushButton(tr("Set"));

    QGridLayout *packagesLayout = new QGridLayout;
//    packagesLayout->addLayout(ffLayout);
    packagesGroup->setLayout(packagesLayout);

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));
    QPushButton *nextButton = new QPushButton(tr("Next"));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(sNext()));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(packagesGroup);
    mainLayout->addLayout(ffLayout);
    mainLayout->addSpacing(12);
//    mainLayout->addWidget(startQueryButton);
    mainLayout->addStretch(12);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);
}

AlgorithmPage::AlgorithmPage(QWidget *parent)
    : QWidget(parent)
{
    QGroupBox *configGroup = new QGroupBox(tr("Algorithm configuration"));

    QLabel    *popLabel = new QLabel(tr("Population size:"));
    popEdit = new QLineEdit(tr("25"));

    QLabel    *genLabel = new QLabel( tr("Maximum generation:") );
    genEdit = new QLineEdit(tr("100"));

    QLabel    *croLabel = new QLabel(tr("Crossover rate:"));
    croEdit = new QLineEdit(tr("0.5"));

    QGridLayout *packagesLayout = new QGridLayout;
    packagesLayout->addWidget(popLabel,0,0);
    packagesLayout->addWidget(popEdit,0,1);
    packagesLayout->addWidget(genLabel,1,0);
    packagesLayout->addWidget(genEdit,1,1);
    packagesLayout->addWidget(croLabel,2,0);
    packagesLayout->addWidget(croEdit,2,1);

    configGroup->setLayout(packagesLayout);

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));
    QPushButton *nextButton = new QPushButton(tr("Next"));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(sNext()));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(configGroup);
    mainLayout->addStretch(12);
    mainLayout->addLayout(buttonsLayout);

    setLayout(mainLayout);
}

void
AlgorithmPage::sPrevious(){
	popSize = popEdit->text().toInt();
	genSize = genEdit->text().toInt();
	crossRate = croEdit->text().toFloat();
	emit previous();
}

void
AlgorithmPage::sNext(){
	popSize = popEdit->text().toInt();
	genSize = genEdit->text().toInt();
	crossRate = croEdit->text().toFloat();
	emit next();
}

//---------------------------------------------------------------------------

RunPage::RunPage(QWidget *parent)
    : QWidget(parent)
{
    QGroupBox *packagesGroup = new QGroupBox(tr("Finished configuration"));
    QGridLayout *packagesLayout = new QGridLayout;
    packagesGroup->setLayout(packagesLayout);

    QPushButton *runButton = new QPushButton(tr("Start docking"));
    connect( runButton, SIGNAL(clicked()), this, SLOT(dock()));
    QVBoxLayout *mainLayout = new QVBoxLayout;

    QPushButton *previousButton = new QPushButton(tr("Previous"));
    connect(previousButton, SIGNAL(clicked()), this, SLOT(sPrevious()));

    QLabel    *popLabel = new QLabel(tr("              "));

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
//    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(popLabel);

    mainLayout->addSpacing(12);
    mainLayout->addWidget(packagesGroup);
    mainLayout->addWidget(runButton);
    mainLayout->addStretch(1);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);
}

//---------------------------------------------------------------------

ConfigDialog::ConfigDialog()
{
    contentsWidget = new QListWidget;
    contentsWidget->setViewMode(QListView::IconMode);
    contentsWidget->setIconSize(QSize(96, 84));
    contentsWidget->setMovement(QListView::Static);
    contentsWidget->setMaximumWidth(128);
    contentsWidget->setSpacing(12);

    pagesWidget = new QStackedWidget;
    proteinPage = new ProteinPage;
    pagesWidget->addWidget( proteinPage);
    connect( proteinPage, SIGNAL(loadProtein(QString)), this, SLOT(getProteinFile(QString)) );
    connect( proteinPage, SIGNAL(previous()), this, SLOT(previous()) );
    connect( proteinPage, SIGNAL(next()), this, SLOT(next()) );

    bindingSitePage = new BindingSitePage;
    pagesWidget->addWidget( bindingSitePage);
    connect( bindingSitePage, SIGNAL(bindingSite(vector<float>)), this, SLOT(getBindingSite(vector<float>)));
    connect( bindingSitePage, SIGNAL(previous()), this, SLOT(previous()) );
    connect( bindingSitePage, SIGNAL(next()), this, SLOT(next()) );


    ligandPage = new LigandPage;
    pagesWidget->addWidget( ligandPage );
    connect( ligandPage, SIGNAL(loadLigand(QString)), this, SLOT(getLigandFile(QString)) );
    connect( ligandPage, SIGNAL(previous()), this, SLOT(previous()) );
    connect( ligandPage, SIGNAL(next()), this, SLOT(next()) );


    scorePage = new ScorePage;
    pagesWidget->addWidget(scorePage);
    connect( scorePage, SIGNAL(previous()), this, SLOT(previous()) );
    connect( scorePage, SIGNAL(next()), this, SLOT(next()) );

    algorithmPage = new AlgorithmPage;
    pagesWidget->addWidget( algorithmPage);
    connect( algorithmPage, SIGNAL(previous()), this, SLOT(previous()) );
    connect( algorithmPage, SIGNAL(next()), this, SLOT(next()) );

    runPage = new RunPage;
    pagesWidget->addWidget(runPage);
    connect( runPage, SIGNAL(run()), this, SLOT(dock()));
    connect( runPage, SIGNAL(previous()), this, SLOT(previous()) );

    createIcons();
    contentsWidget->setCurrentRow(0);

    QPushButton *previousButton =  new QPushButton(tr("Previous"));
    connect( previousButton, SIGNAL(clicked()), this, SLOT(previous()) );

    QPushButton *nextButton =  new QPushButton(tr("Next"));
    connect( nextButton, SIGNAL(clicked()), this, SLOT( next() ) );
    QPushButton *runButton =  new QPushButton(tr("Run"));
    connect( runButton, SIGNAL(clicked()), this, SLOT( run()) );

    QPushButton *closeButton = new QPushButton(tr("Close"));
    connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));

    QHBoxLayout *horizontalLayout = new QHBoxLayout;
    horizontalLayout->addWidget(contentsWidget);
    horizontalLayout->addWidget(pagesWidget, 1);

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(previousButton);
    buttonsLayout->addWidget(nextButton);
    buttonsLayout->addWidget(runButton);
    buttonsLayout->addWidget(closeButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(horizontalLayout);
    mainLayout->addStretch(1);
    mainLayout->addSpacing(12);
//    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);

    setWindowTitle(tr("Config Dialog"));
}

void ConfigDialog::createIcons()
{
    QListWidgetItem *proteinButton = new QListWidgetItem(contentsWidget);
    proteinButton->setIcon(QIcon("images/protein.png1"));
    proteinButton->setText(tr("protein"));
    proteinButton->setTextAlignment(Qt::AlignHCenter);
    proteinButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    QListWidgetItem *bindingSiteButton = new QListWidgetItem(contentsWidget);
    bindingSiteButton->setIcon(QIcon("images/bindingSite.png1"));
    bindingSiteButton->setText(tr("bindingSite"));
    bindingSiteButton->setTextAlignment(Qt::AlignHCenter);
    bindingSiteButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    QListWidgetItem *ligandButton = new QListWidgetItem(contentsWidget);
    ligandButton->setIcon(QIcon("images/ligand.png1"));
    ligandButton->setText(tr("lignd"));
    ligandButton->setTextAlignment(Qt::AlignHCenter);
    ligandButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    QListWidgetItem *scoreButton = new QListWidgetItem(contentsWidget);
    scoreButton->setIcon(QIcon("images/score.png1"));
    scoreButton->setText(tr("score"));
    scoreButton->setTextAlignment(Qt::AlignHCenter);
    scoreButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    QListWidgetItem *algorithmButton = new QListWidgetItem(contentsWidget);
    algorithmButton->setIcon(QIcon("images/algorithm.png1"));
    algorithmButton->setText(tr("algorithm"));
    algorithmButton->setTextAlignment(Qt::AlignHCenter);
    algorithmButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    QListWidgetItem *runButton = new QListWidgetItem(contentsWidget);
    runButton->setIcon(QIcon("images/run.png1"));
    runButton->setText(tr("run"));
    runButton->setTextAlignment(Qt::AlignHCenter);
    runButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    connect(contentsWidget,
            SIGNAL(currentItemChanged(QListWidgetItem*,QListWidgetItem*)),
            this, SLOT(changePage(QListWidgetItem*,QListWidgetItem*)));
}

void ConfigDialog::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
    if (!current)
        current = previous;

    pagesWidget->setCurrentIndex(contentsWidget->row(current));
}

void ConfigDialog::previous( ){
	int ind = contentsWidget->currentRow() - 1;
	if( ind <=0 ){
		ind = 0;
	}
	contentsWidget->setCurrentRow( ind );
}

void ConfigDialog::next(){
	int ind = contentsWidget->currentRow() + 1;
	if( ind >5 ){
		ind = 5;
	}
	contentsWidget->setCurrentRow( ind );
}

void
ConfigDialog::run(){

	emit dockRun();
}

void
ConfigDialog::dock(){

	ofstream off("dock.par");

	off.width(15);
	off<<left<<"MolFile"<<(proteinPage->getFileName())<<endl;

	off.width(15);
    off<<left<<"InitLig"<<ligandPage->getInitLigandName()<<endl;

    off.width(15);
    off<<left<<"ReferLig"<<ligandPage->getRefLigandName()<<endl;

    off.width(15);
    off<<left<<"ForceField"<<"ff.dat"<<endl;;

    off.width(15);
    off<<left<<"BindCenter"<<bindingSitePage->getX()<<" "<<bindingSitePage->getY()<<" "<<bindingSitePage->getZ()<<endl;

    off.width(15);
    off<<left<<"BindRadius"<<bindingSitePage->getR()<<endl;

    off.width(15);
    off<<left<<"PopSize"<<algorithmPage->getPop()<<endl;

    off.width(15);
    off<<left<<"Generation"<<algorithmPage->getGen()<<endl;

    off.width(15);
    off<<left<<"CountMax"<<20<<endl;

    off.width(15);
    off<<left<<"CR"<<algorithmPage->getCro()<<endl;

    off.close();
    close();
	emit dockRun();

}
