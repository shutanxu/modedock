/**
 * A class for converting NOE intensity into distance range
 */
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "normalPDF.h"
using namespace std;

normalPDF::normalPDF() {
	// TODO Auto-generated constructor stub
}

void normalPDF::computeParam( bool printPar ) {
	if ( N == 0 ){
		cerr<<"No data points \n";
		exit(1);
	}
	mu = 0.0;
	double sumSquare = 0.0;
	for ( unsigned int i =0; i<N; i++){
		mu += data[i];
		sumSquare += data[i] * data[i];
	}
	mu /= N;
	sumSquare /= N;
	sigma =sumSquare - mu * mu;
	if (printPar)
		printf("mu=%15.6f%, sigma=%15.6f\n", mu, sigma);
}


void normalPDF::print(){
	cout<<"------------------- Data --------------------------\n";
	for (int i = 0; i < N; i++)
		printf("%15.6f\n", data[i]);
	cout<<endl;
	cout<<"Normal parameters: ";
	printf("mu=%15.6f%, sigma=%15.6f\n", mu, sigma);
}

normalPDF::~normalPDF() {
	delete [] data;
	// TODO Auto-generated destructor stub
}
