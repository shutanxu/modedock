/*
 * normalPDF.h
 *
 *  Created on: Oct 14, 2008
 *      Author: lwang4
 */

#ifndef NORMALPDF_H_
#define NORMALPDF_H_

class normalPDF {
public:
	normalPDF();
	normalPDF( double *const dat, size_t no): data(dat), N(no) {;}

	void computeParam(bool printPar=false);

	virtual ~normalPDF();
    virtual void print();

public:
	double *data;
	size_t N;  //the number of points
	// the output
	double mu;
	double sigma;
//	double sigmaSquare;
//	double deviation;
};

#endif /* NORMALPDF_H_ */
