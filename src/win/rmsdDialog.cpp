#include"rmsdDialog.h"

using namespace std;


LigandRMSDDialog::LigandRMSDDialog(QWidget *parent):QWidget(parent){

    QLabel *labelMol1=new QLabel(tr("lig 1:"));
    editMol1 = new QLineEdit;
    QPushButton *load1 = new QPushButton(tr("Load"));
    connect(load1, SIGNAL(clicked()), this, SLOT(loadFile1()) );

    QHBoxLayout *layout1= new QHBoxLayout;
    layout1->addWidget(labelMol1);
    layout1->addWidget(editMol1);
    layout1->addWidget(load1);

    QLabel *labelMol2=new QLabel(tr("lig 2:"));
    editMol2 = new QLineEdit;
    QPushButton *load2 = new QPushButton(tr("Load"));
    connect(load2, SIGNAL(clicked()), this, SLOT(loadFile2()) );

    QHBoxLayout *layout2= new QHBoxLayout;
    layout2->addWidget(labelMol2);
    layout2->addWidget(editMol2);
    layout2->addWidget(load2);

    QLabel *label1= new QLabel(tr("RMSD:"));
    labelRMSD = new QLabel;
    QHBoxLayout *layout3=new QHBoxLayout;
    layout3->addWidget(label1);
    layout3->addWidget(labelRMSD);

    QLabel *label2= new QLabel(tr("Overlay RMSD:"));
    labelOverlayRMSD = new QLabel;
    QHBoxLayout *layout4=new QHBoxLayout;
    layout4->addWidget(label2);
    layout4->addWidget(labelOverlayRMSD);

    QPushButton *compute=new QPushButton(tr("Compute"));
    connect(compute, SIGNAL(clicked()), this, SLOT(computeRMSD()));
    QHBoxLayout *layout5= new QHBoxLayout;
    layout5->addWidget(compute);

    QVBoxLayout *layoutAll= new QVBoxLayout;
    layoutAll->addLayout(layout1);
    layoutAll->addLayout(layout2);
    layoutAll->addLayout(layout3);
    layoutAll->addLayout(layout4);
    layoutAll->addWidget(compute);
    setLayout(layoutAll);

    QDesktopWidget *desktop = QApplication::desktop();
    move( (desktop->width() - this->width())/2, (desktop->height()-this->height())/2 );
}

void
LigandRMSDDialog::loadFile1(){
    QString file=QFileDialog::getOpenFileName(this,
                                              tr("Open Molecular File"), ".",
                                              tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
    if (!file.isEmpty()){
        editMol1->insert( file );
        fileName1 = file.toStdString();
        mol1 = readMolecularFile(fileName1).front();
    }
}

void
LigandRMSDDialog::loadFile2(){
    QString file=QFileDialog::getOpenFileName(this,
                                              tr("Open Molecular File"), ".",
                                              tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
    if (!file.isEmpty()){
        editMol2->insert( file );
        fileName2 = file.toStdString();
        mol2 = readMolecularFile(fileName2).front();
    }
}

void
LigandRMSDDialog::computeRMSD(){
    float dis = isomorphismRMSD(mol1, mol2);
    float olDis = overlayRMSD(mol1, mol2);
    labelRMSD->setText(QString("%1").arg(dis));
    labelOverlayRMSD->setText(QString("%1").arg(olDis));
    update();
    cout<<"Dis:"<<dis<<endl;
}

//-------------------------------------------------------------------------------------


ProteinRMSDDialog::ProteinRMSDDialog(QWidget *parent):QWidget(parent){

    QLabel *labelMol1=new QLabel(tr("protein 1:"));
    editMol1 = new QLineEdit;
    QPushButton *load1 = new QPushButton(tr("Load"));
    connect(load1, SIGNAL(clicked()), this, SLOT(loadFile1()) );

    QHBoxLayout *layout1= new QHBoxLayout;
    layout1->addWidget(labelMol1);
    layout1->addWidget(editMol1);
    layout1->addWidget(load1);

    QLabel *labelMol2=new QLabel(tr("protein 2:"));
    editMol2 = new QLineEdit;
    QPushButton *load2 = new QPushButton(tr("Load"));
    connect(load2, SIGNAL(clicked()), this, SLOT(loadFile2()) );

    QHBoxLayout *layout2= new QHBoxLayout;
    layout2->addWidget(labelMol2);
    layout2->addWidget(editMol2);
    layout2->addWidget(load2);

    QLabel *label1= new QLabel(tr("Pairwise CA RMSD:"));
    labelRMSD = new QLabel;
    QHBoxLayout *layout3=new QHBoxLayout;
    layout3->addWidget(label1);
    layout3->addWidget(labelRMSD);

    QLabel *label2= new QLabel(tr("Overlay CA RMSD:"));
    labelOverlayRMSD = new QLabel;
    QHBoxLayout *layout4=new QHBoxLayout;
    layout4->addWidget(label2);
    layout4->addWidget(labelOverlayRMSD);

    QPushButton *compute=new QPushButton(tr("Compute"));
    connect(compute, SIGNAL(clicked()), this, SLOT(computeRMSD()));
    QHBoxLayout *layout5= new QHBoxLayout;
    layout5->addWidget(compute);

    QVBoxLayout *layoutAll= new QVBoxLayout;
    layoutAll->addLayout(layout1);
    layoutAll->addLayout(layout2);
    layoutAll->addLayout(layout3);
    layoutAll->addLayout(layout4);
    layoutAll->addWidget(compute);
    setLayout(layoutAll);

    QDesktopWidget *desktop = QApplication::desktop();
    move( (desktop->width() - this->width())/2, (desktop->height()-this->height())/2 );
}

void
ProteinRMSDDialog::loadFile1(){
    QString file=QFileDialog::getOpenFileName(this,
                                              tr("Open Molecular File"), ".",
                                              tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
    if (!file.isEmpty()){
        editMol1->insert( file );
        fileName1 = file.toStdString();
        mol1 = readMolecularFile(fileName1).front();
    }
}

void
ProteinRMSDDialog::loadFile2(){
    QString file=QFileDialog::getOpenFileName(this,
                                              tr("Open Molecular File"), ".",
                                              tr("molecular files (*.pdb *.mol2 *.pdbqt)"));
    if (!file.isEmpty()){
        editMol2->insert( file );
        fileName2 = file.toStdString();
        mol2 = readMolecularFile(fileName2).front();
    }
}

void
ProteinRMSDDialog::computeRMSD(){
    float dis = pairRMSD_CA(mol1, mol2);
    float olDis = overlayRMSD_CA(mol1, mol2);
    string	overlayAtomName("CA");

    Molecular olyMol = overlay( mol1, mol2, overlayAtomName, true );
    vector<Molecular> molVec;
    molVec.push_back(mol1);
    molVec.push_back(olyMol);

    vector<Atom>						atomStd = mol1.get_atomVec();
    vector<Atom>						atomNew = olyMol.get_atomVec();

    emit signalOlyProtein(molVec);

    labelRMSD->setText(QString("%1").arg(dis));
    labelOverlayRMSD->setText(QString("%1").arg(olDis));
    update();
    cout<<"Dis:"<<dis<<endl;
}


