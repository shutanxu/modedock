/*
 * linreg.h
 *
 *  Created on: Sep 22, 2013
 *      Author: lincong
 */

#ifndef LINREG_H_
#define LINREG_H_


template <class T> class linearRegression {
	friend ostream& operator<<(ostream& out, linearRegression& lr) {
		if (lr.haveData())
			out << "f(x) = " << lr.getA()
			<< " + ( " << lr.getB()
			<< " * x )";
		return out;
	}

public:

	// This is also the default constructor

	linearRegression( pair<T, T> *p = 0, long size = 0) {
		long i;
		a = b = sumX = sumY = sumXsquared = sumYsquared = sumXY = 0.0;
		n = 0L;
		if (size > 0L) // if size greater than zero there are data arrays
			for (n = 0, i = 0L; i < size; i++)
				addPoint(p[i]);
	}

	// Constructor using arrays of x values and y values
	linearRegression(T *x, T *y, long size = 0) {
		long i;
		a = b = sumX = sumY = sumXsquared = sumYsquared = sumXY = 0.0;
		n = 0L;
		if (size > 0L) // if size greater than zero there are data arrays
			for (n = 0, i = 0L; i < size; i++)
				addXY(x[i], y[i]);
	}
	virtual ~linearRegression() {;}

	void addXY(const T& x, const T& y);
	void addPoint(const pair<T, T> & p)  { addXY(p.first, p.second); }

	// Must have at least 3 points to calculate
	// standard error of estimate.  Do we have enough data?
	int haveData() const { return (n > 2 ? 1 : 0); }
	long items() const { return n; }

	T getA() const { return a; }
	T getB() const { return b; }

	T getCoefDeterm() const  { return coefD; }
	T getCoefCorrel() const { return coefC; }
	T getStdErrorEst() const { return stdError; }
	virtual T estimateY(T x) const { return (a + b * x); }

protected:
	long n;             // number of data points input so far
	T sumX, sumY;  // sums of x and y
	T sumXsquared; // sum of x squares
	T sumYsquared; // sum y squares
	T sumXY;       // sum of x*y

	T a, b;        // coefficients of f(x) = a + b*x
	T coefD;       // coefficient of determination
	T coefC;       // coefficient of correlation
	T stdError;    // standard error of estimate

	void Calculate();   // calculate coefficients
};

template<class T> void linearRegression::addXY(const T& x, const T& y) {
	n++;
	sumX += x;
	sumY += y;
	sumXsquared += x * x;
	sumYsquared += y * y;
	sumXY += x * y;
	Calculate();
}

template<class T> void linearRegression::Calculate() {
	if (haveData())     {
		if (fabs( T(n) * sumXsquared - sumX * sumX) > DBL_EPSILON)         {
			b = ( T(n) * sumXY - sumY * sumX) /
			( T(n) * sumXsquared - sumX * sumX);
			a = (sumY - b * sumX) / T(n);

			T sx = b * ( sumXY - sumX * sumY / T(n) );
			T sy2 = sumYsquared - sumY * sumY / T(n);
			T sy = sy2 - sx;

			coefD = sx / sy2;
			coefC = sqrt(coefD);
			stdError = sqrt(sy / T(n - 2));
		}        else        {
			a = b = coefD = coefC = stdError = 0.0;
		}
	}
}

#endif /* LINREG_H_ */
