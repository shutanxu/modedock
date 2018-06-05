/*
 * de.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: stan
 */

#include"deDock.h"
#include <ctime>        // std::time

#define PI 3.14159
#define TWO_PI 6.2831853071795864769252866

//--------------------------------------------------------------------------------------------

//DeDock::DeDock( QWidget *parent ):QWidget(parent){}
//DeDock::DeDock(const string dir, QWidget *parent ):QWidget(parent){}

bool
DeDock::readParFile(){
	string  parName;
	parName=parDir+string("dock.par");

	ifstream iff;
	iff.open(parName.c_str());
	string oneLine;
	vector<string> strVec;
	int count=0;
	while( getline(iff, oneLine) ){
		strVec=tokenBySpace(oneLine);
//		cout<<"str:"<<strVec[0]<<endl;
		if( strVec[0]==string("MolFile") ){
			molFileName=strVec[1];
			count++;
		}else if( strVec[0]==string("InitLig") ){
			initLigName=strVec[1];
			count++;
		}else if( strVec[0]==string("ReferLig") ){
			referLigName=strVec[1];
			count++;
		}else if( strVec[0]==string("ForceField") ){
			sybylFileName=strVec[1];
			count++;
		}else if( strVec[0]==string("BindCenter") ){
			float x=atof(strVec[1].c_str());
			float y=atof(strVec[2].c_str());
			float z=atof(strVec[3].c_str());
			Coord c(x,y,z);
			bindCenter=c;
			count++;
		}else if( strVec[0]==string("BindRadius") ){
			bindRadius=atof(strVec[1].c_str());
			count++;
		}else if( strVec[0]==string("PopSize") ){
			popSize=atoi(strVec[1].c_str());
			count++;
		}else if( strVec[0]==string("Generation") ){
			generation=atoi(strVec[1].c_str());
			count++;
		}else if( strVec[0]==string("CountMax") ){
			countMax=atoi(strVec[1].c_str());
			count++;
		}else if( strVec[0]==string("CR") ){
			CR=atof(strVec[1].c_str());
			count++;
		}
	}
	if( count==10 ){
		return true;
	}else{
		cout<<"error!!! can not read parameter file"<<endl;
		cout<<"count:"<<count<<endl;
		return false;
	}
}
/*
 * lig1 dominates lig2 in HBond energy and VDW energy
 */
bool
DeDock::dominatesHbVDW( Sybyl sybyl, const vector<Atom> receptorAtomVec,
		const Molecular& lig1, const Molecular& lig2, const bool& flag ){

	double hb1 = getHBondEnergy( receptorAtomVec, lig1 );
	double hb2 = getHBondEnergy( receptorAtomVec, lig2 );

//	float vdw1 = sybyl.computeGoldVDW( receptorAtomVec, lig1.get_atomVec() );
//	float vdw2 = sybyl.computeGoldVDW( receptorAtomVec, lig2.get_atomVec() );

	float vdw1 = sybyl.computeCutVDW( receptorAtomVec, lig1.get_atomVec() );
	float vdw2 = sybyl.computeCutVDW( receptorAtomVec, lig2.get_atomVec() );


	if(flag){
		cout<<"hb1:";
		cout.width(10);
		cout<<left<<hb1;
		cout<<" hb2:";
		cout.width(10);
		cout<<left<<hb2;
		cout<<"vdw1:";
		cout.width(10);
		cout<<left<<vdw1<<" vdw2:"<<vdw2;
	}

	if( ( hb1 >= hb2 && vdw1 < vdw2 ) || (hb1 > hb2 && vdw1 <= vdw2) ){
		if(flag){
			cout<<"   @"<<endl;
		}
		return true;
	}else{
		if(flag){
			cout<<endl;
		}
		return false;
	}
}

bool
DeDock::dominatesVDW( Sybyl sybyl, const vector<Atom> receptorAtomVec, const Molecular& lig1,
        const Molecular& lig2, const bool& computeHydrogen ){

    float vdw1 = sybyl.computeCutVDW( receptorAtomVec, lig1.get_atomVec(), computeHydrogen );
    float vdw2 = sybyl.computeCutVDW( receptorAtomVec, lig2.get_atomVec(), computeHydrogen );

    if(0){
		cout<<"vdw1:";
		cout.width(12);
		cout<<left<<vdw1<<" vdw2:";
		cout.width(12);
		cout<<vdw2;
	}

	if(   vdw1 < vdw2 ){
        if(0){
			cout<<"   @";
		}
		return true;
	}else{
        if(0){
//			cout<<endl;
		}
		return false;
	}
}

bool
DeDock::dominatesESP( Sybyl sybyl, const vector<Atom> receptorAtomVec, const Molecular& lig1,
		const Molecular& lig2, const bool& flag ){

	float val1 = sybyl.computeESP( receptorAtomVec, lig1.get_atomVec() );
	float val2 = sybyl.computeESP( receptorAtomVec, lig2.get_atomVec() );

	if(flag){
		cout<<"esp1:";
		cout.width(12);
		cout<<left<<val1<<" esp2:";
		cout.width(12);
		cout<<val2;
	}

	if(  val1 < val2 ){
		if(flag){
			cout<<"   @";
		}
		return true;
	}else{
		if(flag){
//			cout<<endl;
		}
		return false;
	}
}

double
DeDock::getHBondEnergy(const Atom &donor, const vector<Atom> &acceptorAtomVec){

}

/*
double
DeDock::getHBondEnergy( const Atom& donor, const vector<Atom>& acceptorAtomVec ){
	vector<double> valVec;
	for( size_t i=0; i<acceptorAtomVec.size(); i++ ){
		Atom at = acceptorAtomVec[i];
		double val = 0;
		double dis = getCoordDis( donor.get_coord(), at.get_coord() );
		if( dis == 0 ){
			val = 0;
		}else if( dis < ( donor.get_vdwRadius()+at.get_vdwRadius() ) ){
			val = 0;
		}else{
			val = 1.0/dis;
		}

		valVec.push_back( val );
	}
	double maxVal = *max_element( valVec.begin(), valVec.end() );
	return maxVal;
}

double
DeDock::getHBondEnergy( const vector<Atom>& donorVec, const vector<Atom>& acceptorAtomVec ){
    if( donorVec.empty() || acceptorAtomVec.empty() ){
        return 0;
    }
    double valSum = 0;
	for( size_t i=0; i<donorVec.size(); i++ ){
		valSum += getHBondEnergy( donorVec[i], acceptorAtomVec );
	}
	return valSum;
}

double
DeDock::getHBondEnergy( const vector<Atom>& receptorAtoms, const Molecular& ligand ){

    if(0){
        for( size_t i=0; i<receptorAtoms.size(); i++ ){
            cout<<i<<" ";
            receptorAtoms[i].print();
        }
    }

	string							goldParFile = "ff.dat";
	Sybyl							goldPar( goldParFile );

	vector<Atom> hDonorAtomVec_mol = goldPar.getHydrogenBondDonorAtoms( receptorAtoms );
	vector<Atom> hDonorAtomVec_lig = goldPar.getHydrogenBondDonorAtoms( ligand.get_atomVec() );

	vector<Atom> hAcceptorAtomVec_mol = goldPar.getHydrogenBondAcceptorAtoms( receptorAtoms );
	vector<Atom> hAcceptorAtomVec_lig = goldPar.getHydrogenBondAcceptorAtoms( ligand.get_atomVec() );
	double hbEnergy = 0;
	hbEnergy += getHBondEnergy( hDonorAtomVec_mol, hAcceptorAtomVec_lig );
	hbEnergy += getHBondEnergy( hDonorAtomVec_lig, hAcceptorAtomVec_mol );

	return hbEnergy;
}
*/

/*
double
DeDock::getHBondEnergy( const Atom& donor, const Atom& acceptor ){

    if( donor.get_attachedH().empty() ){
        return 0;
    }

    float dis = atomDis( donor, acceptor );
    if( dis == 0 || dis > 4.0 ){
        return 0;
    }
    float eps_1=0;
    float eps_2=0;
    float r_eqm1=0;
    float r_eqm2=0;
    float r1=0;
    float r2=0;

    string t1=donor.get_type().substr(0,1);
    string t2=acceptor.get_type().substr(0,1);

    if( t1=="N" || t1=="O" ){
        eps_1=5.0;
        r_eqm1=1.9;
    }else if( t1=="S" ){
        eps_1=1.0;
        r_eqm1=2.5;
    }else{
        return 0;
    }

    if( t2=="N" || t2=="O" ){
        eps_2=5.0;
        r_eqm2=1.9;
    }else if( t2=="S" ){
        eps_2=1.0;
        r_eqm2=2.5;
    }else{
        return 0;
    }

    float r12 = pow( dis, 12 );
    float r10 = pow( dis, 10 );

    float epsi_12 = sqrt( eps_1 * eps_2 );
    float r_eqmXY = ( acceptor.get_vdwRadius() + donor.get_vdwRadius() ) / 2.0 ;

    float C_ij = 5 * epsi_12 * pow( r_eqmXY, 12 );
    float D_ij = 6 * epsi_12 * pow( r_eqmXY, 10 );

    float energy = ( C_ij / r12 ) - ( D_ij / r10 ) ;
    if(1){
        cout<<"h: ";
        vector<Atom> hVec=donor.get_attachedH();
        hVec[0].print();
        cout<<"d: ";
        donor.print();
        cout<<"a: ";
        acceptor.print();
        vector<Atom> cVec=acceptor.get_attachedC();
        cout<<"c: ";
        cVec[0].print();
        cout<<"dis:"<<dis<<"    energy: "<<energy<<endl;
        cout<<endl;
    }
    return energy ;
}
*/

double
DeDock::getHBondEnergy( const Atom& donor, const Atom& acceptor ){

    if( donor.get_attachedH().empty() ){
        return 0;
    }

    float dis = atomDis( donor, acceptor );
    if( dis == 0 || dis > 3.5 ){
        return 0;
    }

    vector<Atom> hVec=donor.get_attachedH();
    if( hVec.empty() ){
        return 0;
    }
    Atom atH;
    float disH=10e10;
    for( size_t i=0; i<hVec.size(); i++ ){
        float dis = atomDis(hVec[i], acceptor);
        if( disH > dis ){
            disH = dis;
            atH = hVec[i];
        }
    }
    float weightDis = 0;
    if( disH < 2.5){
        weightDis = 1;
    }else{
        weightDis=(3.5-disH);
        weightDis = weightDis < 0 ? 0 : weightDis;
    }

    // X-H--Y-R
    float angleXHY=atomAngle( donor, atH, acceptor );
    float weightAngle1=0;
    if( angleXHY > 160 ){
        weightAngle1 = 1;
    }else{
        weightAngle1 = ( angleXHY - 120)/40.0;
        weightAngle1 = weightAngle1 < 0 ? 0:weightAngle1;
    }

    // X-H--Y-R
    Atom atC=acceptor.get_attachedC().front();
    float angleHYR=atomAngle( atH, acceptor, atC );
    float weightAngle2=0;
    if( angleHYR > 100 ){
        weightAngle2=1;
    }else{
        weightAngle2=(angleHYR-90)/10.0;
        weightAngle2=weightAngle2 < 0 ? 0 : weightAngle2;
    }

    float weight=weightDis*weightAngle1*weightAngle2;
    if(0){
        cout<<"h: ";
        atH.print();
        cout<<"d: ";
        donor.print();
        cout<<"a: ";
        acceptor.print();
        vector<Atom> cVec=acceptor.get_attachedC();
        cout<<"c: ";
        cVec[0].print();
        cout<<"dis:"<<dis<<"    weight: "<<weight<<endl;
        cout<<endl;
    }
    return weight*-6.0;
}

double
DeDock::getHBondEnergy( const vector<Atom>& donorVec, const vector<Atom>& acceptorAtomVec ){
    if( donorVec.empty() || acceptorAtomVec.empty() ){
        return 0;
    }
    double valSum = 0;
    for( size_t i=0; i<donorVec.size(); i++ ){
        for( size_t j=0; j<acceptorAtomVec.size(); j++ ){
            valSum += getHBondEnergy( donorVec[i], acceptorAtomVec[j] );
        }
    }
    return valSum;
}

double
DeDock::getHBondEnergy(const vector<Atom> &receptorAtoms, const Molecular &ligand){
    string							goldParFile = "ff.dat";
    Sybyl							goldPar( goldParFile );

    vector<Atom> hDonorAtomVec_mol = goldPar.getHydrogenBondDonorAtoms( receptorAtoms );
    vector<Atom> hDonorAtomVec_lig = goldPar.getHydrogenBondDonorAtoms( ligand.get_atomVec() );

    vector<Atom> hAcceptorAtomVec_mol = goldPar.getHydrogenBondAcceptorAtoms( receptorAtoms );
    vector<Atom> hAcceptorAtomVec_lig = goldPar.getHydrogenBondAcceptorAtoms( ligand.get_atomVec() );
    if(0){
        for( size_t i=0; i<hAcceptorAtomVec_lig.size(); i++ ){
            hAcceptorAtomVec_lig[i].print();
        }
    }
    double hbEnergy = 0;
    hbEnergy += getHBondEnergy( hDonorAtomVec_mol, hAcceptorAtomVec_lig );
    hbEnergy += getHBondEnergy( hDonorAtomVec_lig, hAcceptorAtomVec_mol );

    return hbEnergy;
}

double
DeDock::generateGaussianNoise(const double &variance)
{
	static bool haveSpare = false;
	static double rand1, rand2;

	if(haveSpare)
	{
		haveSpare = false;
		return sqrt(variance * rand1) * sin(rand2);
	}

	haveSpare = true;

	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;

	return sqrt(variance * rand1) * cos(rand2);
}

Molecular
DeDock::getTranslateRotateMol( const Molecular& mol, const int centerAtomIndex,
		const vector<float>& posVec ){

	float transX = posVec[0];
	float transY = posVec[1];
	float transZ = posVec[2];
	float transDis = posVec[3];

	float rotX = posVec[4];
	float rotY = posVec[5];
	float rotZ = posVec[6];

	Coord dir(transX, transY, transZ);
	Molecular newLig = translateMol( mol, dir, transDis );

	Molecular newLig2 = rotateAroundAtom( newLig, centerAtomIndex, rotX, rotY, rotZ );

	return newLig2;
}

vector<Bond>
DeDock::get_bfsRotableBonds( const Molecular& ligand, const int& centerAtomIndex ){
	vector<Bond>      bfsBondVec = getBfsBonds(ligand.get_bondVec(), centerAtomIndex );
	vector<Bond>      rotableBonds = getRotableBonds( ligand );

	vector<Bond>      bfsRotableBonds;
	for( size_t i=0; i<bfsBondVec.size(); i++ ){
		if( std::find(rotableBonds.begin(), rotableBonds.end(), bfsBondVec[i] )
				!= rotableBonds.end() ){
			bfsRotableBonds.push_back(bfsBondVec[i]);
		}
	}

	return bfsRotableBonds;
}

vector<Atom>
DeDock::getBindSiteAtoms( const Molecular& mol, const Coord& bindingSite, const float& radius ){
	vector<Atom> atVec = mol.get_atomVec();
	vector<Atom> bindAtVec;
	for( size_t i=0; i<atVec.size(); i++ ){
		Coord c = atVec[i].get_coord();
		if( getCoordDis( bindingSite, c ) < radius ){
			bindAtVec.push_back(atVec[i]);
		}
	}

	return bindAtVec;
}

bool
DeDock::ligandInBindingSite( const Molecular& mol, const Coord& bindCenter, const float& bindRadius ){
	vector<Atom>  atVec = mol.get_atomVec();
	for( size_t i=0; i<atVec.size(); i++ ){
		if( getCoordDis( atVec[i].get_coord(), bindCenter ) > ( bindRadius-atVec[i].get_vdwRadius() ) ){
			return false;
		}
	}
	return true;
}

Molecular
DeDock::getMolFromPosVector(  const Molecular& mol, const int centerAtomIndex,
        const vector<float>& posVec, const vector<Bond> rotableBondVec, bool computeRotableBonds ){

	if( rotableBondVec.size() != posVec.size()-7 ){
		cout<<"ligand pos decomposition error"<<endl;
	}
	Molecular lig = mol ;
    if( computeRotableBonds ){
        Molecular ligTemp;
        for( size_t i=0; i<rotableBondVec.size(); i++ ){
            rotateAroundBond( lig, rotableBondVec[i].get_firstAtomIndex(),
                              rotableBondVec[i].get_secAtomIndex(), posVec[i], ligTemp );
            lig=ligTemp;
        }
    }

	size_t startIndex = rotableBondVec.size();
	float transX = posVec[startIndex];
	float transY = posVec[startIndex+1];
	float transZ = posVec[startIndex+2];
	float transDis = posVec[startIndex+3];

	Coord dir(transX, transY, transZ);
	Molecular newLig = translateMol( lig, dir, transDis );

	float rotX = posVec[startIndex+4];
	float rotY = posVec[startIndex+5];
	float rotZ = posVec[startIndex+6];

	Molecular newLig2 = rotateAroundAtom( newLig, centerAtomIndex, rotX, rotY, rotZ );

	return newLig2;
}

bool
DeDock::initFlexLigandPopulation( ){

	posPos.clear();
	population.clear();
    int countMax = 20;

	vector<Atom>      receptorAtoms = getBindSiteAtoms( mol, bindCenter, bindRadius );

    Atom              centerAtom = ligand.get_centerAtom();

	vector<Atom>      atomVec = ligand.get_atomVec();
	int count = 0;

	time_t t1,tStart;
	tStart=t1 = clock();
	int countBD=0;
	int countTr=0;
	int countRO=0;

    emit emitString( QString("initializing...") );
//    cout<<1<<endl;
	while( posPos.size() < popSize ){

		/*********************************************************************************/
		// rotate bonds
		vector<float> bdDegreeVec;
		Molecular lig_bd_rot;
		lig_bd_rot = ligand;
		Molecular lig_new_rot;
		Bond rotateBond;

		for( size_t i=0; i<rotableBonds.size(); i++ ){
			float rDegree=0;
            int countR = 0;
			do{
//				rDegree=(rand()/(float)RAND_MAX)*360.0;
				rDegree=(rand()/(float)RAND_MAX)*20.0-10;
                countR ++ ;
            }while( !rotateAroundBond(lig_bd_rot,rotableBonds[i],rDegree, lig_new_rot) && countR < 20);

            if( countR < 20 ){
                lig_bd_rot=lig_new_rot;
                bdDegreeVec.push_back(rDegree);
            }else{
                bdDegreeVec.push_back(0.0);
            }
		}

		if(0){
			stringstream ss;
			countBD++;
			ss<<countBD;
			string name= string("dockResult/init_rotBD_")+ ss.str();
			name += string(".pdb");
			lig_bd_rot.print(name);
		}

		/*********************************************************************************/
		//translate ligand
		int tempCount=posPos.size();
		int tempCount2=0;        //try times
		int count1 = 0;          //translate times
		int count2 = 0;          //rotate times
		while( tempCount==posPos.size() && tempCount2<50 ){
			tempCount2++;

			Molecular ligand_rot, ligand_trans ;
			Coord newCenter;
			Coord dir;
			float dis;
			count1 = 0;
			count2 = 0;

			float roundX, roundY, roundZ;
			do{
				float ligandRadius = lig_bd_rot.get_centerAtomRadius();
				centerAtom = lig_bd_rot.get_centerAtom();

				double ranX = (rand()/(double)RAND_MAX-0.5)*(double)(bindRadius-ligandRadius);
				double ranY = (rand()/(double)RAND_MAX-0.5)*(double)(bindRadius-ligandRadius);
				double ranZ = (rand()/(double)RAND_MAX-0.5)*(double)(bindRadius-ligandRadius);
				Coord newCenter( ranX+bindCenter.x, ranY+bindCenter.y, ranZ+bindCenter.z );
				dir = newCenter-centerAtom.get_coord();
				dis = getCoordDis(newCenter, centerAtom.get_coord() );
				ligand_trans=translateMol(lig_bd_rot, dir, dis);
				count1++;
            }while( ( atomsCrash(ligand_trans.get_centerAtom(), receptorAtoms, false ) ||
					!ligandInBindingSite(ligand_trans, bindCenter, bindRadius) ) &&
					count1 < countMax);
			if(0 && count1<countMax){
				stringstream ss;
				countTr++;
				ss<<countTr;
				string name= string("dockResult/init_trans_")+ ss.str();
				name += string(".pdb");
				ligand_trans.print(name);
			}

			//rotate ligand
			if( count1 < countMax ){
				do{
					roundX = ( rand()/(double)RAND_MAX )*360;
					roundY = ( rand()/(double)RAND_MAX )*360;
					roundZ = ( rand()/(double)RAND_MAX )*360;
					ligand_rot = rotateAroundAtom(ligand_trans, ligCenterAtomIndex, roundX, roundY, roundZ);
					count2++;
                }while( ( ( atomsCrash( ligand_rot.get_atomVec(), receptorAtoms, false ) ) ||
						!ligandInBindingSite(ligand_rot, bindCenter, bindRadius) ) &&
						count2 < countMax);
			}
			if(0 && count2<400){
				stringstream ss;
				countRO++;
				ss<<countRO;
				string name= string("dockResult/init_")+ ss.str();
				name += string(".pdb");
				ligand_rot.print(name);
			}
			/*********************************************************************************/

			if( count1 < countMax && count2 < countMax ){
				vector<float> ligandPos;
				ligandPos.insert( ligandPos.end(), bdDegreeVec.begin(), bdDegreeVec.end() );
				ligandPos.push_back( dir.x );
				ligandPos.push_back( dir.y );
				ligandPos.push_back( dir.z );
				ligandPos.push_back( dis );
				ligandPos.push_back( roundX );
				ligandPos.push_back( roundY );
				ligandPos.push_back( roundZ );

				Molecular testlig = getMolFromPosVector(ligand, ligCenterAtomIndex,
						ligandPos, rotableBonds);

				if( ligandInBindingSite(testlig, bindCenter, bindRadius) &&
                        !atomsCrash(testlig.get_atomVec(), receptorAtoms, false  ) ){

					//incase the 'ligandPos' already exist in 'posPos'
					bool flag = true;
					for( size_t i=0; i<posPos.size(); i++ ){
						int matchCount=0;
						for( size_t j=0; j<posPos[i].size(); j++ ){
							if( posPos[i][j] == ligandPos[j] ){
								matchCount++;
							}
						}
						if( matchCount == ligandPos.size() ){
							flag = false;
							break;
						}
					}
					if( flag ){
						posPos.push_back( ligandPos );
						population.push_back( ligand_rot );
						if(0){
							stringstream ss;
							ss<<countRO;
							string name= string("dockResult/init_")+ ss.str();
							string name2= string("dockResult/init_")+ ss.str()+string("test");
							name += string(".pdb");
							name2 += string(".pdb");
							ligand_rot.print(name);
							testlig.print(name2);
							countRO++;
						}
					}
				}else{
                    if(0){
						stringstream ss;
						ss<<countRO;
						string name= string("dockResult/init_rotate")+ ss.str();
						string name2= string("dockResult/init_testlig")+ ss.str();
						name += string(".pdb");
						name2 += string(".pdb");
						ligand_rot.print(name);
						testlig.print(name2);
						countRO++;
					}
				}
			}
		}
		if(1){
			stringstream ss;
			ss<<"posPos:";
			ss.width(5);
			ss<<left<<posPos.size();
			ss<<"popSize:";
			ss.width(5);
			ss<<popSize;
			ss<<"count1:";
			ss.width(5);
			ss<<count1;
			ss<<"count2:";
			ss.width(5);
			ss<<count2;
			ss<<"tempCount2:";
			ss.width(5);
			ss<<tempCount2;
			string line = ss.str();
			emitString(QString(line.c_str()));
            cout<<line<<endl;
            if( posPos.size() == 0 ){
                return false;
            }
		}
	}
	stringstream ss;
    ss<<"init running time:"<<double( clock() - tStart )/ (CLOCKS_PER_SEC * 60 )<<" minutes"<<endl;
	emitString(QString(ss.str().c_str()));
    cout<<"init running time:"<<double( clock() - tStart )/ (CLOCKS_PER_SEC * 60 )<<" minutes"<<endl;
	return true;
}

bool
DeDock::deFlexLigandDock( ){

	string logFile;
	logFile=parDir+string("dock.log");
	ofstream off;
	off.open(logFile.c_str());

	off<<"ligCenterAtomIndex:"<<ligCenterAtomIndex<<endl;

    if(0){
		for( size_t i=0; i<posPos.size(); i++ ){
			stringstream ss;
			ss<<i;
			string name= parDir+string("init_")+ ss.str();
			name += string(".pdb");
			Molecular mmol = getMolFromPosVector(ligand, ligCenterAtomIndex, posPos[i], rotableBonds );
			mmol.print(name);
		}
	}
	vector<Atom>  receptorAtoms = getBindSiteAtoms(mol, bindCenter, bindRadius);
	vector<Bond>  bfsRotableBonds = get_bfsRotableBonds( ligand, ligCenterAtomIndex );

	double refHBondEnergy=getHBondEnergy(receptorAtoms, refLig );
    double refVdwEnergy=sybyl.computeCutVDW(receptorAtoms, refLig.get_atomVec(), false);
    double refESP=sybyl.computeESP(receptorAtoms, refLig.get_atomVec());
    cout<<"dist between init and refer:"<<isomorphismRMSD( ligand, refLig, true );
    vector<Molecular> population(popSize);
    time_t tStart=clock();
	int N = posPos[0].size();     //dimension
	int k=0;
    bool earlyTerminate = false;
    while( k < generation && !earlyTerminate ){
		k++;

        if(1){
			cout<<endl;
			cout<<molFileName<<endl;
			cout<<"generation:"<<k<<" "<<generation<<endl;
			cout<<"refer ligand HB:"<<refHBondEnergy<<" VDW:"<<refVdwEnergy<<" ESP:"<<refESP<<endl;
			cout<<"crash:"<<atomsCrash( ligand.get_atomVec(), receptorAtoms )<<endl;
			cout<<"inner:"<<ligandInBindingSite(ligand, bindCenter, bindRadius)<<endl;
		}

		off<<endl;
		off<<"generation:"<<k<<" "<<generation<<endl;
		off<<"refer ligand HB:"<<refHBondEnergy<<" VDW:"<<refVdwEnergy<<endl;

		if(1){
			emit emitString(QString(molFileName.c_str()));
			emit emitString(QString("generation %2").arg(k)+QString(" %2").arg(generation));
//			emit emitString(QString("refer ligand HB:%2").arg(refHBondEnergy)+ QString(" VDW:%2").arg(refVdwEnergy));

            stringstream ss;
            ss<<"running time:"<<double( clock() - tStart )/ (CLOCKS_PER_SEC*60)<<" minutes"<<endl;
            emitString(QString(ss.str().c_str()));
            off<<"running time:"<<double( clock() - tStart )/ (CLOCKS_PER_SEC*60)<<" minutes"<<endl;
            cout<<"running time:"<<double( clock() - tStart )/ (CLOCKS_PER_SEC*60)<<" minutes"<<endl;
        }

		for( size_t ligIndex = 0; ligIndex<popSize; ligIndex++ ){
			srand (time(NULL));
			size_t r1 = rand()%popSize;
			size_t r2 = rand()%popSize;
			size_t r3 = rand()%popSize;
			size_t i_rand = rand()%N;
			vector<float> newPos(N);
			newPos=posPos[ligIndex];

            Molecular preLig = getMolFromPosVector( ligand, ligCenterAtomIndex,
                    posPos[ligIndex], rotableBonds, false );


			if(0){
				for( size_t i=0; i<posPos[ligIndex].size(); i++ ){
					cout.width(12);
					cout<<posPos[ligIndex][i];
				}
				cout<<endl;
			}

            Molecular newLig = preLig ;
			int count = 0;
			//----------------------------------------------------------------------------
			//rotate around bonds
			Molecular ligand_bd;
			ligand_bd = ligand;
			Molecular lig_new_rot;
			int countRotBd=0;
            bool considerRotableBonds=false; // time consumming to compute rotable bonds
            if( considerRotableBonds ){
                for( size_t i=0; i<rotableBonds.size(); i++ ){
                    float ra=rand()/(float)RAND_MAX;
                    if( ra>0.5 ){
                        float rDegree=posPos[ligIndex][i]+ ( ( rand()/(float)RAND_MAX )*20 -10 );
                        newPos[i]=rDegree;
                    }else{
                        newPos[i]=posPos[ligIndex][i];
                    }
                    rotateAroundBond(ligand_bd, rotableBonds[i], newPos[i], lig_new_rot );
                    ligand_bd=lig_new_rot;
                }
            }

			//----------------------------------------------------------------------------
			//translate and rotate ligand
			bool updateFlag = false;
			int co=0;
            countMax = 10;
            while( !updateFlag && co<10 ){
				co++;

				count = 0;
				do{
					count++;
					vector<float> molRotateTransVec(7);
					if( (N-rotableBonds.size()) != 7 ){
						cout<<"!!! error, not appropriate coding"<<endl;
						off<<"!!! error, not appropriate coding"<<endl;
					}
					for( size_t i=rotableBonds.size(); i<N; i++ ){
						if( rand()/((double)RAND_MAX)<CR /*&& i < i_rand*/ ){
							double gn = generateGaussianNoise(1);
							float x = posPos[r3][i] + gn*(posPos[r1][i]-posPos[r2][i]) ;
							newPos[i] = x;
						}else{
							newPos[i] = posPos[ligIndex][i];
						}
						molRotateTransVec[i-rotableBonds.size()] = newPos[ i ];
					}
					newLig = getTranslateRotateMol( ligand_bd, ligCenterAtomIndex, molRotateTransVec );
				}while( count < countMax &&
                        ( ( atomsCrash( newLig.get_atomVec(), receptorAtoms, false ) ) ||
								( !ligandInBindingSite(newLig, bindCenter, bindRadius) ) ) );

				if( count < countMax ){

					if( dominatesVDW( sybyl, receptorAtoms, newLig, preLig, false) ){
						posPos[ligIndex] = newPos;
                        population[ligIndex]=newLig;
						updateFlag=true;
					}
				}else{
//					dominatesVDW( sybyl, receptorAtoms, newLig, preLig, false);
				}
			}

			if(1){

				stringstream lineStream;

                lineStream<<"index";
                lineStream.width(4);
                lineStream<<left<<ligIndex<<" ";
				lineStream<<"count:";
				lineStream.width(4);
				lineStream<<count<<" ";
				lineStream.width(4);
				lineStream<<co<<" ";


				off<<"index:";
				off.width(4);
				off<<left<<ligIndex<<" ";
				off<<"count:";
				off.width(4);
				off<<count<<" ";
				off.width(4);
				off<<co<<" ";

                float dis=isomorphismRMSD( refLig, preLig );
//				float dis1=overlayRMSD( refLig, preLig );
				float dis1=0;

				lineStream<<"RMSD:";
				lineStream.width(8);
				lineStream<<left<<dis<<" ";
				lineStream<<"olRMSD:";
				lineStream.width(12);
				lineStream<<left<<dis1<<" ";

				off<<"RMSD:";
				off.width(8);
				off<<left<<dis<<" ";
				off<<"olRMSD:";
				off.width(12);
				off<<left<<dis1<<" ";

                float hbEnergy=getHBondEnergy( receptorAtoms, preLig );
//                float hbEnergy= 0;
				lineStream<<"HB:";
				lineStream.width(12);
				lineStream<<hbEnergy<<" ";

				off<<"HB:";
				off.width(12);
				off<<hbEnergy<<" ";

//				float esp=sybyl.computeESP(receptorAtoms, preLig.get_atomVec());
                float esp=0;
				lineStream<<"ESP:";
				lineStream.width(12);
				lineStream<<esp<<" ";

				off<<"ESP:";
				off.width(12);
				off<<esp<<" ";

                float vdw=sybyl.computeCutVDW(receptorAtoms, preLig.get_atomVec(), false );
//                float vdw= 0;
				lineStream<<"vdw:";
				lineStream.width(12);
				lineStream<<vdw;

				off<<"vdw:";
				off.width(12);
				off<<vdw;
                if( 0/*updateFlag*/ ){
					lineStream<<" update";
					off<<" update";
					if( countRotBd >0 ){
						lineStream<<" rb:"<<countRotBd;
						off<<" rb:"<<countRotBd<<endl;
					}else{
						off	<<endl;
					}
				}else{
					off<<endl;
				}

				string line = lineStream.str();
				emitString(QString(line.c_str()));
                cout<<line<<endl;
			}
		}
		vector<Molecular> allMol;
        if(0){
			for( size_t i=0; i<posPos.size(); i++ ){
				stringstream ss, s2;
				ss<<i;
				s2<<k;
				string name= parDir+string("gen")+ s2.str() + string("_pos") + ss.str();
				name += string(".pdb");
				Molecular mmol = getMolFromPosVector(ligand, ligCenterAtomIndex, posPos[i], rotableBonds);
				mmol.print(name);
				allMol.push_back( mmol );
			}
		}
	}

	if(1){
        finalPop.clear();
		for( size_t i=0; i<posPos.size(); i++ ){
			stringstream ss;
			ss<<i;
            string name= parDir + string("final_lig")+ ss.str();
			name += string(".pdb");
			Molecular mmol = getMolFromPosVector(ligand, ligCenterAtomIndex, posPos[i], rotableBonds);
            mmol.set_name(name);
            finalPop.push_back(mmol);
			mmol.print(name);
		}
	}
	off.close();
    emit finishRun();
    cout<<"####################"<<endl;
	return true;
}

void
DeDock::sampleRotableBonds( const Molecular& ligand,
		const size_t bdIndex,
		vector<float>& angleVec,
		vector<vector<float> >& allAngleVec,
		vector<Molecular>& ligRotateVec ){

	if( bdIndex > rotableBonds.size()-1 ){
		float tolDis = -2.0;
		if( !ligand.innerCrash(tolDis) ){
			allAngleVec.push_back(angleVec);
			ligRotateVec.push_back(ligand);
		}
		return;
	}

	for( float angle=0; angle<360; angle+=90 ){
		Bond rotBd = rotableBonds[bdIndex];
		Molecular newLig;
		if(rotateAroundBond(ligand, rotBd, angle, newLig) ){
			angleVec.push_back(angle);
		}else{
			newLig = ligand;
			angleVec.push_back(angle);
		}
		sampleRotableBonds( newLig, bdIndex+1, angleVec,
				allAngleVec, ligRotateVec );
	}
	return;
}

bool
DeDock::sampleRotableBonds(){
	//sample all the possible rotable bonds
	vector<Molecular> rotateLigVec;
	rotateLigVec.clear();
	vector<float> angleVec;
	vector< vector<float> > allAngleVec;

	float sc=0;
	vector<float> scoreVec;
	sampleRotableBonds( ligand, 0, angleVec,
			allAngleVec, rotateLigVec );

	vector<float> allOverlayDis;
	vector< pair<float, float> > disScoreVec;
	vector< pair<float, float> > scoreDisVec;
	if(1){
		for( size_t i=0; i<rotateLigVec.size(); i++ ){
			stringstream ss;
			ss<<i;
			string name=parDir + string("init_rotBD_")+ss.str();
			name+=string(".pdb");
			Molecular olMol;//=overlay(refLig, rotateLigVec[i]);
			olMol.print(name);
			float overDis = 0; //overlayRMSD(refLig, olMol);
			allOverlayDis.push_back(overDis);
//			cout<<i<<"     overDis:"<<overDis<<endl;
			float vdw=sybyl.computeVDW(rotateLigVec[i]);
			cout<<"ligand vdw:"<<vdw<<endl;
//			pair<float, float> p(overDis, scoreVec[i]);
//			pair<float, float> p2( scoreVec[i], overDis);
			pair<float, float> p(overDis, vdw);
			pair<float, float> p2( vdw, overDis);

			disScoreVec.push_back(p);
			scoreDisVec.push_back(p2);
		}
	}

	sort( disScoreVec.begin(), disScoreVec.end() );
	sort( scoreDisVec.begin(), scoreDisVec.end() );
	cout<<"dis:"<<endl;
	for( size_t i=0; i<disScoreVec.size(); i++ ){
		cout<<disScoreVec[i].first<<endl;
	}
	cout<<"score:"<<endl;
	for( size_t i=0; i<disScoreVec.size(); i++ ){
		cout<<disScoreVec[i].second<<endl;
	}

	cout<<"!!!!!!!"<<endl;
	return true;
}

bool
DeDock::shuffleLigand(){
	Molecular tempMol=ligand;
	Molecular newMol;

	for( size_t i=0; i<rotableBonds.size(); i++){
		float degree= (rand()/(float)RAND_MAX)*360.0;
		if( rotateAroundBond(tempMol, rotableBonds[i], degree, newMol) ){
			tempMol=newMol;
		}
	}

	float dirX=(rand()/(float)RAND_MAX);
	float dirY=(rand()/(float)RAND_MAX);
	float dirZ=(rand()/(float)RAND_MAX);
	Coord dir(dirX, dirY, dirZ);
	float dis=(rand()/(float)RAND_MAX)*50.0;
	tempMol=translateMol(tempMol, dir, dis);

	double roundX = ( rand()/(double)RAND_MAX )*360;
	double roundY = ( rand()/(double)RAND_MAX )*360;
	double roundZ = ( rand()/(double)RAND_MAX )*360;
	tempMol = rotateAroundAtom(tempMol, ligCenterAtomIndex, roundX, roundY, roundZ);


	string newName("shuffLigand.pdb");
	newName = parDir+newName;
	tempMol.print(newName);
	ligand=tempMol;
	return true;
}

bool
DeDock::init(){
	if( !readParFile() ){
		return false;
	}
    emit emitString( QString("start docking...") );
    mol = readMolecularFile( molFileName ).front();
	refLig = readMolecularFile( referLigName ).front();
    ligand = readMolecularFile( initLigName ).front();

    if(0){
        float dis = overlayRMSD( refLig, ligand );
        cout<<"dis:"<<dis<<endl;
        return false;
    }

//    vector<Atom> atVec = ligand.get_atomVec();
//    float r= atomsMaxRange( atVec )/2.0;
//    if( r > bindRadius ){
//        emitString(QString( "binding site radius is too small for initial ligand" ));
//        cout<<"ligand too large, exceeds the binding site radius"<<endl;
//        return false;
//    }

	ligCenterAtomIndex = ligand.get_centerAtom().get_index();
	rotableBonds = ligand.get_rotableBondVec();
	Sybyl			  sy( sybylFileName );
	sybyl = sy;

    if(0){
		cout<<"rotable bonds"<<endl;
		cout<<rotableBonds.size()<<endl;
		for( size_t i=0; i<rotableBonds.size(); i++ ){
			cout<<i<<endl;
			rotableBonds[i].print();
		}
	}
	if(0){
		sampleRotableBonds();
	}
	return true;
}

void
DeDock::run( ){
    if( !init() ){
        return;
    }
    if( !initFlexLigandPopulation() ){
        return;
    }
	deFlexLigandDock();

}
