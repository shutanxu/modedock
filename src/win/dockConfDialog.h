/*
 * dockConfDialog.h
 *
 *  Created on: Jan 18, 2015
 *      Author: stan
 */

#ifndef DOCKCONFDIALOG_H_
#define DOCKCONFDIALOG_H_

#include <QWidget>
#include <QtWidgets>
#include<string>
#include<vector>

using namespace std;

class ProteinPage : public QWidget
{
	Q_OBJECT
public:
    ProteinPage(QWidget *parent = 0);
    string getFileName()const{ return fileName; }
//    void accept();
private slots:
    void loadFile();
    void sPrevious(){ emit previous(); }
    void sNext(){ emit next(); }
signals:
    void loadProtein( QString file );
    void previous();
    void next();
private:
    string     fileName;
    QLineEdit *proteinEdit;
};

class BindingSitePage : public QWidget
{
	Q_OBJECT
public:
    BindingSitePage(QWidget *parent = 0);
//    vector<float> getBindingSite();
    float getX(){ return x; };
    float getY(){ return y; };
    float getZ(){ return z; };
    float getR(){ return r; };

private slots:
    void getData();
    void sPrevious(){ emit previous(); }
    void sNext(){ emit next(); }

signals:
    void bindingSite( vector<float> vec );
    void previous();
    void next();
private:
    QLineEdit  *xEdit;
    QLineEdit  *yEdit;
    QLineEdit  *zEdit;
    QLineEdit  *rEdit;

    float x,y,z,r;
};

class LigandPage : public QWidget
{
	Q_OBJECT
public:
    LigandPage(QWidget *parent = 0);
    string getInitLigandName(){ return fileName; }
    string getRefLigandName(){ return refFileName; }
private slots:
    void loadFile();
    void loadRefFile();
    void sPrevious(){ emit previous(); }
    void sNext(){ emit next(); }

signals:
    void loadLigand( QString file );
    void loadRefLigand( QString file );
    void previous();
    void next();
private:
    string fileName;
    string refFileName;
    QLineEdit *ligandEdit;
    QLineEdit *refLigandEdit;
};

class ScorePage : public QWidget
{
	Q_OBJECT
public:
    ScorePage(QWidget *parent = 0);
private slots:
    void sPrevious(){ emit previous(); }
    void sNext(){ emit next(); }
signals:
    void previous();
    void next();
};

class AlgorithmPage : public QWidget
{
	Q_OBJECT
public:
    AlgorithmPage(QWidget *parent = 0);
    int getPop(){ return popSize; }
    int getGen(){ return genSize; }
    float getCro(){ return crossRate; }
private slots:
    void sPrevious();//{ emit previous(); }
    void sNext();//{ emit next(); }
signals:
    void previous();
    void next();

private:
    QLineEdit *popEdit;
    QLineEdit *genEdit;
    QLineEdit *croEdit;

    int popSize;
    int genSize;
    float crossRate;
};

class RunPage : public QWidget
{
	Q_OBJECT
public:
    RunPage(QWidget *parent = 0);
private slots:
    void sPrevious(){ emit previous(); }
    void sNext(){ emit next(); }
    void dock(){ emit run(); }
signals:
    void previous();
    void next();
    void run();
};

//-------------------------------------------------------------------

class ConfigDialog : public QDialog
{
    Q_OBJECT

public:
    ConfigDialog();

public slots:
    void changePage(QListWidgetItem *current, QListWidgetItem *previous);
    void previous();
    void next();
    void run();
    void getProteinFile( QString s ){ emit loadP(s); }
    void getBindingSite( vector<float> d ){ emit loadB(d); }
    void getLigandFile( QString s ){ emit loadL(s); }
    void dock();
signals:
    void loadP( QString file );
    void loadB( vector<float> );
    void loadL( QString file );
    void dockRun();
private:
    void createIcons();

    QListWidget *contentsWidget;
    QStackedWidget *pagesWidget;

    ProteinPage     *proteinPage;
    BindingSitePage *bindingSitePage;
    LigandPage      *ligandPage;
    ScorePage       *scorePage;
    AlgorithmPage   *algorithmPage;
    RunPage         *runPage;
};


#endif /* DOCKCONFDIALOG_H_ */
