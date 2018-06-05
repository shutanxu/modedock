/*
 * overlay.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: stan
 */


#include"overlay.h"
#include"../basicMath/geomBasic.h"

/**
 * http://en.wikipedia.org/wiki/Wahba%27s_problem
 * http://en.wikipedia.org/wiki/Kabsch_algorithm
 * @param structCoordA			reference coordinate, array type
 * @param structCoordB
 * @param debug
 * @return
 */

vector< vector<double> >
overlayMatrix(const vector<float>&structCoordA,
         const vector<float>&structCoordB, const bool& debug ){

    vector<vector<double> >		stdMatrix( 3, vector<double>(3, 0) );
    stdMatrix[0][0] = stdMatrix[1][1] = stdMatrix[2][2] = 1;

    vector<double>				arrA( structCoordA.size(), 0 );
    vector<double>				arrB( structCoordB.size(), 0 );
    vector<double>				centerA( 3, 0 );
    vector<double> 				centerB( 3, 0 );
    size_t						numOfAtom = structCoordA.size() / 3;

    for ( unsigned int i = 0; i < numOfAtom ; i++){
        for ( unsigned int j = 0; j< 3 ; j++) {
            centerA[ j ] += structCoordA[ 3 * i + j];
            centerB[ j ] += structCoordB[ 3 * i + j];
        }
    }
    for ( unsigned  int j = 0; j<3; j++) {
            centerA[j]  /= numOfAtom;
            centerB[j]  /= numOfAtom;
    }

    for (unsigned int i=0; i<numOfAtom; i++){
        for (unsigned int j = 0; j<3; j++){
            arrA[ 3*i+j ] = structCoordA[ 3 * i + j] - centerA[j] ;
            arrB[ 3*i+j ] = structCoordB[ 3 * i + j] - centerB[j] ;
        }
    }
    double  arrU[ 3 * 3 ];
    for (unsigned int i = 0; i< 3 ; i++ ){
        for  (unsigned int j = 0; j< 3; j++) {
            arrU[ i * 3 + j ] = 0.0;
            for (unsigned int k = 0; k< numOfAtom; k++)
                arrU[ i * 3 + j ]  += arrB [k * 3 + j ] * arrA[ k * 3 + i ];
        }
    }

/*
 * A = U S V^T
*/
    gsl_matrix_view 		U				= gsl_matrix_view_array (arrU, 3, 3);
    gsl_matrix *            V 				= gsl_matrix_alloc (3, 3);
    gsl_vector *			sValues 		= gsl_vector_alloc ( 3);
    gsl_vector *			work 			= gsl_vector_alloc ( 3 );
    gsl_matrix *			rotMatrix 		= gsl_matrix_alloc (3, 3);
    gsl_matrix *			deterMatrix 	= gsl_matrix_alloc (3, 3);
    gsl_matrix * 			unitMatrix 	= gsl_matrix_alloc( 3, 3 );
    gsl_matrix * 			tempMatrix 	= gsl_matrix_alloc( 3, 3 );

    if( debug ){
        cout<<"centerA: "<<centerA[0]<<" "<<centerA[1]<<" "<<centerA[2]<<endl;
        cout<<"centerB: "<<centerB[0]<<" "<<centerB[1]<<" "<<centerB[2]<<endl;
        cout<<"correlate matrix:"<<endl;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                cout.width(15);
                cout<<gsl_matrix_get( &U.matrix, i, j );
            }
            cout<<endl;
        }
    }

    gsl_linalg_SV_decomp (&U.matrix, V, sValues, work);
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, &U.matrix, V,  0, deterMatrix );
    double					deterM = gsl_linalg_LU_det( deterMatrix, 3 );
    double					deterU = gsl_linalg_LU_det( &U.matrix, 3 );
    double					deterV = gsl_linalg_LU_det( V, 3 );

    for( size_t i=0; i<3; i++ ){
        for( size_t j=0; j<3; j++ ){
            gsl_matrix_set( unitMatrix, i, j, 0 );
        }
    }
    double					d = 0;
    if( deterU*deterV > 0 ){
        d = 1;
    }else if( deterU*deterV < 0 ){
        d = -1;
    }
    gsl_matrix_set( unitMatrix, 0, 0, 1 );
    gsl_matrix_set( unitMatrix, 1, 1, 1 );
//	gsl_matrix_set( unitMatrix, 2, 2, deterU*deterV );
    gsl_matrix_set( unitMatrix, 2, 2, d );

    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, &U.matrix, unitMatrix, 0, tempMatrix );
    gsl_matrix_transpose( V );
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, tempMatrix, V,  0, rotMatrix );

    if( debug ){
        cout<<"V:"<<endl;
        cout<<"deterU:"<<deterU<<endl;
        cout<<"deterV:"<<deterV<<endl;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                cout.width(15);
                cout<<gsl_matrix_get( V, i, j );
            }
            cout<<endl;
        }
        cout<<"U.matrix:"<<endl;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                cout.width(15);
                cout<<gsl_matrix_get( &U.matrix, i, j );
            }
            cout<<endl;
        }
        cout<<"unitMatrix:"<<endl;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                cout.width(15);
                cout<<gsl_matrix_get( unitMatrix, i, j );
            }
            cout<<endl;
        }
    }

    vector< vector<double> >		rotateMatrix( 3, vector<double>(3, 0) );
    for( size_t i=0; i<3; i++ ){
        for( size_t j=0; j<3; j++ ){
            rotateMatrix[i][j] = gsl_matrix_get( rotMatrix, i, j );
        }
    }

    //-----------------------------------------------------------------------------------------------
    vector< vector<double> >		rotateMatrix2( 3, vector<double>(3, 0) );
    //-----------------------------------------------------------------------------------------------

    gsl_matrix_free (V);
    gsl_matrix_free (rotMatrix);
    gsl_matrix_free (deterMatrix);
    gsl_matrix_free (unitMatrix);
    gsl_matrix_free (tempMatrix);
    gsl_vector_free (work);
    gsl_vector_free (sValues);

    if(debug){
        cout<<"rotate matrix:"<<endl;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                cout.width( 15 );
                cout<<rotateMatrix[i][j];
            }
            cout<<endl;
        }
    }

    if(1){
        matrix rotM;
        double rmsd;
        bool isz = false;
        vector<float> vecA, vecB;
        for( size_t i=0; i<arrA.size(); i++ ){
            vecA.push_back( arrA[i] );
            vecB.push_back( arrB[i] );
        }
        rotM = rotMatrixByRMSD4C( vecA, vecB, isz, rmsd);
        size_t c = rotM.getColumnDimension();
        size_t r = rotM.getRowDimension();
        vector<double> array = rotM.getArray();
        vector< vector<double> >		rotateMatrixNew( 3, vector<double>(3, 0) );
        float sum=0;
        for( size_t i=0; i<3; i++ ){
            for( size_t j=0; j<3; j++ ){
                rotateMatrixNew[i][j] = array[i*3+j];
                sum+=array[i*3+j];
            }
        }
        if(sum!=0){
            rotateMatrix=rotateMatrixNew;
        }

    }

    return		rotateMatrix;
}

vector<float>
overlayCoord( const vector<float>&structCoordA,
              const vector<float>&structCoordB,
              const bool& debug ){

    if( structCoordA.size() != structCoordB.size() ){
        cout<<"error"<<endl;
//        cout<<structCoordA.size() <<" "<<structCoordB.size()<<endl;
//        return 0;
    }
    vector< vector<double> > rotMatrix=overlayMatrix(structCoordA, structCoordB);
    size_t numOfAtom = structCoordA.size()/3;

    vector<float> centerA(3,0);
    vector<float> centerB(3,0);

    for( size_t i=0; i<numOfAtom; i++ ){
        for( size_t j=0; j<3; j++ ){
            centerA[j]+=structCoordA[3*i+j];
            centerB[j]+=structCoordB[3*i+j];
        }
    }
    for( size_t i=0; i<3; i++ ){
        centerA[i]/=numOfAtom;
        centerB[i]/=numOfAtom;
    }

    vector<float> coordNewVec;
    for(size_t i=0; i<numOfAtom; i++){
        vector<float> coordNew(3,0);
        vector<float> coordOld(3,0);

        coordOld[0]=structCoordB[3*i]-centerB[0];
        coordOld[1]=structCoordB[3*i+1]-centerB[1];
        coordOld[2]=structCoordB[3*i+2]-centerB[2];

        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNew[m] += rotMatrix[m][n] * coordOld[n];
            }
        }
        coordNew[0]+=centerA[0];
        coordNew[1]+=centerA[1];
        coordNew[2]+=centerA[2];

        coordNewVec.push_back(coordNew[0]);
        coordNewVec.push_back(coordNew[1]);
        coordNewVec.push_back(coordNew[2]);
    }
    return coordNewVec;
}

float
overlayRMSD( const vector<float>&structCoordA,
const vector<float>&structCoordB, const bool& debug ){
    if( structCoordA.size() != structCoordB.size() ){
        cout<<"error"<<endl;
//        cout<<structCoordA.size() <<" "<<structCoordB.size()<<endl;
        return 0;
    }
    vector< vector<double> > rotMatrix=overlayMatrix(structCoordA, structCoordB);
    size_t numOfAtom = structCoordA.size()/3;

    vector<float> centerA(3,0);
    vector<float> centerB(3,0);

    for( size_t i=0; i<numOfAtom; i++ ){
        for( size_t j=0; j<3; j++ ){
            centerA[j]+=structCoordA[3*i+j];
            centerB[j]+=structCoordB[3*i+j];
        }
    }
    for( size_t i=0; i<3; i++ ){
        centerA[i]/=numOfAtom;
        centerB[i]/=numOfAtom;
    }

    vector<float> coordNewVec;
    for(size_t i=0; i<numOfAtom; i++){
        vector<float> coordNew(3,0);
        vector<float> coordOld(3,0);

        coordOld[0]=structCoordB[3*i]-centerB[0];
        coordOld[1]=structCoordB[3*i+1]-centerB[1];
        coordOld[2]=structCoordB[3*i+2]-centerB[2];

        for( size_t m=0; m<3; m++ ){
            for( size_t n=0; n<3; n++ ){
                coordNew[m] += rotMatrix[m][n] * coordOld[n];
            }
        }
        coordNew[0]+=centerA[0];
        coordNew[1]+=centerA[1];
        coordNew[2]+=centerA[2];

        coordNewVec.push_back(coordNew[0]);
        coordNewVec.push_back(coordNew[1]);
        coordNewVec.push_back(coordNew[2]);
    }

    float val=0;
    for( size_t i=0; i<numOfAtom; i++ ){
        val += (coordNewVec[i*3]-structCoordA[i*3])*(coordNewVec[i*3]-structCoordA[i*3])+
               (coordNewVec[i*3+1]-structCoordA[i*3+1])*(coordNewVec[i*3+1]-structCoordA[i*3+1])+
               (coordNewVec[i*3+2]-structCoordA[i*3+2])*(coordNewVec[i*3+2]-structCoordA[i*3+2]);
    }

    val/=(3*numOfAtom);
    val = sqrt(val);

    return val;
}
