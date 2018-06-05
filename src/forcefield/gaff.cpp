#include<string.h>
#include"gaff.h"


using namespace std;

bool
ForceFieldGAFF::ParseParamFile(){
    vector<string> vs;
//    int BUFF_SIZE=50;
    char buffer[BUFF_SIZE];
    ForceFieldParameter parameter;

    ifstream ifs;
    ifs.open("gaff.dat");
    if(!ifs){
        cout<<"can not open gaff.dat!!!"<<endl;
        return false;
    }
    ifs.getline(buffer, BUFF_SIZE);
    ifs.getline(buffer, BUFF_SIZE);
    cout<<"!!!!"<<endl;
    while( !ifs.fail() && strlen(buffer) !=0 ){
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa = vs[0]; //KNDSYM
        parameter.dpar.push_back(atof(vs[1].c_str())); // AMASS
        parameter.dpar.push_back(atof(vs[2].c_str())); // ATPOL [A^3]
        ffpropparams.push_back(parameter);
        ifs.getline(buffer,BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE);
    ifs.getline(buffer, BUFF_SIZE);
    while( !ifs.fail() && strlen(buffer) != 0 ){ //read block 4 (bonds)
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0]; //IBT
        parameter.sb=vs[1]; //JBT
        parameter.dpar.push_back(atof(vs[2].c_str())); // RK [kcal/mol/(A^2)]
        parameter.dpar.push_back(atof(vs[3].c_str())); // REQ [A]
        ffbondparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE);
    while( !ifs.fail() && strlen(buffer) != 0 ){  //read block 5 (angles)
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0]; //ITT
        parameter.sb=vs[1]; //JTT
        parameter.sc=vs[2]; //KTT
        parameter.dpar.push_back(atof(vs[3].c_str())); // TK [kcal/mol/(rad**2)]
        parameter.dpar.push_back(atof(vs[4].c_str())); // TEQ [degrees]
        ffangleparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE); //read block 6 (dihedrals)
    while( !ifs.fail() && strlen(buffer) != 0 ){
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0];  //IPT
        parameter.sb=vs[1];  //JPT
        parameter.sc=vs[2];  //KPT
        parameter.sd=vs[3];  //LPT
        parameter.dpar.push_back(atof(vs[4].c_str())); // IDIVF
        parameter.dpar.push_back(atof(vs[5].c_str())); // PK
        parameter.dpar.push_back(atof(vs[6].c_str())); // GAMMA [degrees]
        parameter.dpar.push_back(atof(vs[7].c_str())); // PN
        fftorsionparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE);
    while( !ifs.fail() && strlen(buffer) !=0 ){//read block 7(impropers)
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0]; //IPT
        parameter.sb=vs[1]; //JPT
        parameter.sc=vs[2]; //KPT
        parameter.sd=vs[3]; //LPT
        parameter.dpar.push_back(atof(vs[4].c_str())); // PK
        parameter.dpar.push_back(atof(vs[5].c_str())); // GAMMA
        parameter.dpar.push_back(atof(vs[6].c_str())); // PN
        ffoopparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while( !ifs.fail() && strlen(buffer) !=0 ){
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0]; //KT1
        parameter.sb=vs[1]; //KT2
        parameter.dpar.push_back(atof(vs[2].c_str())); // A
        parameter.dpar.push_back(atof(vs[3].c_str())); // B
        ffhbondparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE); //next line
    while(!ifs.fail() && strlen(buffer)!=0) //read block 9 (equivalent atoms for vdW potential)
    {
        // not used. We assume the types are all listed individually in block 10
        ifs.getline(buffer, BUFF_SIZE);
    }
    ifs.getline(buffer, BUFF_SIZE); //next line
    ifs.getline(buffer, BUFF_SIZE); //'RE' means vdW radius and well-depth are read
    while( !ifs.fail() && strlen(buffer) !=0 ){ //read block 10 (vdWaals)
        // van der Waals potental EDEP * ( (R/r)^12 - 2*(R/r)^6 )
        tokenize(vs, buffer," -\n\t");
        parameter.clear();
        parameter.sa=vs[0]; //IBT
        parameter.dpar.push_back(atof(vs[1].c_str())); //R
        parameter.dpar.push_back(atof(vs[2].c_str())); //EDEP(kcal/mol)
        ffvdwparams.push_back(parameter);
        ifs.getline(buffer, BUFF_SIZE);
    }
    if(ifs){
        ifs.close();
    }

    if(0){
        for( size_t i=0; i<ffpropparams.size(); i++ ){
            cout<<ffpropparams[i].sa<<" "<<ffpropparams[i].dpar[0]<<" "<<ffpropparams[i].dpar[1]<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<ffbondparams.size(); i++ ){
            cout<<ffbondparams[i].sa<<" ";
            cout<<ffbondparams[i].sb<<" ";
            cout<<ffbondparams[i].dpar[0]<<" ";
            cout<<ffbondparams[i].dpar[1]<<" ";
            cout<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<ffangleparams.size(); i++ ){
            cout<<ffangleparams[i].sa<<" ";
            cout<<ffangleparams[i].sb<<" ";
            cout<<ffangleparams[i].sc<<" ";
            cout<<ffangleparams[i].dpar[0]<<" ";
            cout<<ffangleparams[i].dpar[1]<<" ";
            cout<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<fftorsionparams.size(); i++ ){
            cout<<fftorsionparams[i].sa<<" ";
            cout<<fftorsionparams[i].sb<<" ";
            cout<<fftorsionparams[i].sc<<" ";
            cout<<fftorsionparams[i].sd<<" ";
            cout<<fftorsionparams[i].dpar[0]<<" ";
            cout<<fftorsionparams[i].dpar[1]<<" ";
            cout<<fftorsionparams[i].dpar[2]<<" ";
            cout<<fftorsionparams[i].dpar[3]<<" ";
            cout<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<ffoopparams.size(); i++ ){
            cout<<ffoopparams[i].sa<<" ";
            cout<<ffoopparams[i].sb<<" ";
            cout<<ffoopparams[i].sc<<" ";
            cout<<ffoopparams[i].sd<<" ";
            cout<<ffoopparams[i].dpar[0]<<" ";
            cout<<ffoopparams[i].dpar[1]<<" ";
            cout<<ffoopparams[i].dpar[2]<<" ";
            cout<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<ffhbondparams.size(); i++ ){
            cout<<ffhbondparams[i].sa<<" ";
            cout<<ffhbondparams[i].sb<<" ";
            cout<<ffhbondparams[i].dpar[0]<<" ";
            cout<<ffhbondparams[i].dpar[1]<<" ";
            cout<<endl;
        }
        cout<<endl;
        for( size_t i=0; i<ffvdwparams.size(); i++ ){
            cout<<ffvdwparams[i].sa<<" ";
            cout<<ffvdwparams[i].dpar[0]<<" ";
            cout<<ffvdwparams[i].dpar[1]<<" ";
            cout<<endl;
        }
        cout<<endl;
    }

    return true;
}

float
ForceFieldGAFF::calculateVDW(const Atom a1, const Atom a2){
    float dis = atomDis(a1, a2);
    float r1,r2,e1,e2;
    r1=r2=e1=e2=0;
    float rVDWab,eab;
    rVDWab=eab=0;
    for( size_t i=0; i<ffvdwparams.size(); i++ ){
        if( a1.get_type() == ffvdwparams[i].sa ){
            r1=ffvdwparams[i].dpar[0];
            e1=ffvdwparams[i].dpar[1];
        }
        if( a2.get_type() == ffvdwparams[i].sa ){
            r2=ffvdwparams[i].dpar[0];
            e2=ffvdwparams[i].dpar[1];
        }
    }
    rVDWab=(r1+r2);
    eab=sqrt(e1*e2);

    double term=rVDWab/dis;
    double term6=term*term*term;
    term6=term6*term6;
    double term12=term6*term6;

    float energy=eab*(term12-2*term6) * KCAL_TO_KJ;
    return energy;
}

float
ForceFieldGAFF::calculateVDW(const Molecular mol){
    vector<Atom> atVec = mol.get_atomVec();
    vector<Bond> bdVec = mol.get_bondVec();
    float sum=0;

    vector<vector<size_t> > matrix = connectMatrix(atVec, bdVec);

    for( size_t i=0; i<atVec.size()-1; i++ ){
        for( size_t j=i+1; j<atVec.size(); j++ ){
            Atom a1=atVec[i];
            Atom a2=atVec[j];
            if( matrix[i][j]>3 ){
                float val=calculateVDW(a1, a2);
                sum+=val;
            }else if( matrix[i][j] == 3 ){
                float val=calculateVDW(a1, a2);
                sum+=val/2.0 ;
            }
        }
    }
    cout<<"vdw energy:"<<sum<<endl;
    return sum;
}
