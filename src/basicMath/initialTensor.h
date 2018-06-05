/**
 * preCond: three dimensional
 * There is an implementation in CGAL of an algorithm for computing the minimum enclosing ellipsoid
 * we should use it, instead of the initial tensor here.
 */

#ifndef INITIALTENSOR_H_
#define INITIALTENSOR_H_

#include <vector>

#include "eigenValueDecomposition.h"
//#include "../molecule/bio/proteinStruct/Atom.h"
using namespace std;

class initialTensor{
public:
	const static size_t dim = 3;

	initialTensor();
	initialTensor( const vector<float>& coord): coordVec(coord) {;}
//	initialTensor( const vector<Atom*>& atomVec);
	initialTensor(const vector<float>& t, const vector<float>& c,const size_t i, const vector<float> &fp):
		tensor(t), center(c), fpIndex(i), farthestPoint(fp){};

	bool initialization();

	matrix principalMoment( vector< float > &eigenValue, vector<double>& eigenVector );
	virtual ~initialTensor();

	vector<float> getTensor() const {return tensor;}
	vector<float> getCenter()const {return center;}
	size_t getFPindex() const {return fpIndex;};
	vector<float> getFP()const { return farthestPoint;}
	float getFarthestDistance () const { return farthestDistance; }

public:
	// input
	vector<float> coordVec;
	//output
	vector<float> tensor; // a 3 x 3 array
	vector<float> center;
	size_t fpIndex;
	vector<float> farthestPoint;
	float farthestDistance;
};

#endif /*INITIALTENSOR_H_*/
