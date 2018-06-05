#ifndef RMSDDIALOG_H
#define RMSDDIALOG_H

#include <QWidget>
#include <QtWidgets>
#include<string>
#include<vector>
#include"../mol/molecular.h"

using namespace std;

class LigandRMSDDialog:public QWidget{
    Q_OBJECT

public:
    LigandRMSDDialog(QWidget *parent=0);

public slots:
    void loadFile1();
    void loadFile2();
    void computeRMSD();
signals:

private:
    QLineEdit *editMol1;
    QLineEdit *editMol2;
    string fileName1;
    string fileName2;
    Molecular mol1;
    Molecular mol2;
    QLabel *labelRMSD;
    QLabel *labelOverlayRMSD;
};

class ProteinRMSDDialog:public QWidget{
    Q_OBJECT

public:
    ProteinRMSDDialog(QWidget *parent=0);

public slots:
    void loadFile1();
    void loadFile2();
    void computeRMSD();
signals:
	void	signalOlyProtein(vector<Molecular> molVec);
private:
    QLineEdit *editMol1;
    QLineEdit *editMol2;
    string fileName1;
    string fileName2;
    Molecular mol1;
    Molecular mol2;
    QLabel *labelRMSD;
    QLabel *labelOverlayRMSD;
};

#endif // RMSDDIALOG_H
