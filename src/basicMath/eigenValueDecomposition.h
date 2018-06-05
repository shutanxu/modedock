#ifndef EIGENVALUEDECOMPOSITION_H_
#define EIGENVALUEDECOMPOSITION_H_

#include "matrix.h"

class eigenValueDecomposition{
public:
	eigenValueDecomposition();
	eigenValueDecomposition (matrix *Arg, bool issymmetric);
	eigenValueDecomposition (const size_t dim, bool issym, double ** const data):
		n(dim), issymmetric(issym), V(data) { ; };

	virtual ~eigenValueDecomposition();
	
	bool decomposeSymmetric();

	matrix* getV() const {return new matrix(V, n, n);}
	double ** getData() const {return V;};
	double *getRealEigenValues() const {return d;};
	double *getImagEigenValues() const {return e;};
	matrix* getD() const ;
	
private:
	void tred2();
	void tql2();
	void orthes();
	void cdiv(double xr, double xi, double yr, double yi) ;
	void hqr2 () ;
		
public:
	// input
	int n;
	bool issymmetric;
    double **V;
    // intermediate variables or output
	double *d;
	double *e;
	double **H;
	double *ort;
};


/**
 * must be initialized first by the constructor eigenValueDecomposition (const size_t dim, bool issym, double **data)
 */
inline bool eigenValueDecomposition::decomposeSymmetric(){
	  d = new double[n];
	  e = new double[n];
	  // Tridiagonalize.
	  tred2();
	  // Diagonalize.
	  tql2();
	  return true;
}

#endif /*EIGENVALUEDECOMPOSITION_H_*/
