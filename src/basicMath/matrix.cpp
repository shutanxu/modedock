#include "matrix.h"
using namespace std;

float** similarityTransformZ( float **a, float sinT, float cosT);

matrix::matrix(){}
matrix::~matrix(){
//	for(int i =0; i < m; i++)
//		delete [] A[i];
//	delete [] A;
	A.clear();
}
 /** Construct an m-by-n matrix of zeros.
 @param m    Number of rows.
 @param n    Number of colums.
 */
matrix::matrix (const size_t m, const size_t n) {
	this->m = m;
	this->n = n;
	A = vector<double>(n * m);
}

 /** Construct an m-by-n constant matrix.
 @param m    Number of rows.
 @param n    Number of colums.
 @param s    Fill the matrix with this scalar value.
 */
 matrix::matrix (const size_t m, const size_t n, double s) {
    this->m = m;
    this->n = n;
    A = vector<double>(n * m , s);
 }
 /** Construct a matrix from a 2-D array.
 @param A    Two-dimensional array of doubles.
 @exception  IllegalArgumentException All rows must have the same length
 */
 matrix::matrix(double**  a) {
	 A = vector<double>(n * m, 0.0);
	 for (unsigned int i = 0; i < m; i++)
		 for (unsigned int j = 0; j < n; j++)
			 A [ n*i+j ]= a[i][j];
 }

 /** Construct a matrix quickly without checking arguments.
 @param A    Two-dimensional array of doubles.
 @param m    Number of rows.
 @param n    Number of colums.
 */
 matrix::matrix( float **a, const size_t m, const size_t n) {
    this->m = m;
    this->n= n;
//    cout<<"m ="<<m<<endl;
    A = vector<double>(n * m , 0.0);
//    cout<<A [ 0)<<"  "<<A [ n * m - 1)<<endl;
    for (unsigned int i = 0; i < m; i++)
    	for (unsigned int j = 0; j < n; j++)
    		A [ n*i+j ]= a[i][j];
 }

 matrix::matrix( double **a, const size_t mm, const size_t nn) {
     m = mm;
     n= nn;
     A = vector<double>(n * m , 0.0);
//     cout<<A [ 0)<<endl;
     for (unsigned int i = 0; i < m; i++)
    	 for (unsigned int j = 0; j < n; j++)
    		 A [ n*i+j ]= a[i][j];
  }

 /** Construct a matrix from a one-dimensional packed array
 @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
 @param m    Number of rows.
 @exception  IllegalArgumentException Array length must be a multiple of m.
 */
 matrix::matrix(double* vals, const size_t m, const size_t size) {
	 this->m = m;
	 n = (m != 0 ?  size / m : 0);
	 if (m*n != size) {
		 cerr<<"Array length must be a multiple of m."<<endl;
		 exit(0);
	 }
	 A = vector<double>(size);
	 for (unsigned int i = 0; i<size; i++)
		 A [ i ] = vals[i];
 }

 matrix::matrix( const matrix& a) {
     m = a.m;
     n= a.n;
     A = a.A;
  }

 matrix& matrix::operator=( const matrix& a) {
      m = a.m;
      n= a.n;
      A = a.A;
      return *this;
   }

//
// /**
//     Make a deep copy of a matrix
// */
// matrix matrix::copy () {
//    matrix X = new matrix(m,n);
//    double** C = X.getArray();
//    for (unsigned int i = 0; i < m; i++) {
//       for (unsigned int j = 0; j < n; j++) {
//          C[i][j] = A[i][j];
//       }
//    }
//    return X;
// }


  /** Copy the internal two-dimensional array.
 @return     Two-dimensional array copy of matrix elements.
   */
  double** matrix::getArrayCopy () const {
	  double** mat = new double*[m];
	  for(unsigned int i = 0; i<m; i++)
		  mat[i] = new double[n];
	  for (unsigned int i = 0; i < m; i++)
		  for (unsigned int j = 0; j < n; j++)
			  mat[i][j] = A [ n * i + j ];
	  return mat;
  }

  /** -a of every elements
   */
 matrix matrix::uminus () const {
	 vector<double> C ( m * n , 0.0);
	 for (unsigned int i = 0; i < m*n ; i++)
			 C[ i ] = -A [ i ];
	 return matrix(C, m, n);
 }

 /** C = A + B
 @param B    another matrix
 @return     A + B
 */
 matrix matrix::plus (const matrix& B)  const{
	 if( B.m != m || B.n!= n){
		 cerr<<"Dimension must match when adds two matrices"<<endl;
		 exit(0);
	 }
	 vector<double> C ( m * n , 0.0);
	 for (unsigned int i = 0; i < m*n ; i++)
		 C[ i ] =  A [ i ] + B.A[i];
	 return matrix(C, m, n);
 }

 /** C = A - B
 @param B    another matrix
 @return     A - B
 */

 matrix matrix::minus (const matrix& B) const {
	 if( B.m != m || B.n!= n){
		 cerr<<"Dimension must match when adds two matrices"<<endl;
		 exit(0);
	 }
	 vector<double> C ( m * n , 0.0);
	 for (unsigned int i = 0; i < m*n ; i++)
		 C[ i ] =  A [ i ] - B.A[i];
	 return matrix(C, m, n);
 }

 /** Multiply a matrix by a scalar, C = s*A
 @param s    scalar
 @return     s*A
 */
 matrix matrix::times (double s) const{
	 vector<double> C ( m * n, 0.0);
	 for (unsigned int i = 0; i < m*n ; i++)
		 C[ i ] =  s*A [ i ];
	 return matrix(C, m, n);
 }

// /** Linear algebraic matrix multiplication, A * B
// @param B    another matrix
// @return     matrix product, A * B
// @exception  IllegalArgumentException matrix inner dimensions must agree.
// */
//  matrix matrix::times (const matrix& B) const {
//	 if (B.m != n) {
//		 cerr<<"matrix inner dimensions must agree."<<endl;
//		 exit(0);
//	 }
////	 int mb = B.m;
//	 int nb = B.n;
////	 vector<double> arrayB = B.getArray();
//	 vector<double>C(m * nb , 0.0) ;
////	 double Bcolj[n];
////	 double* Arowi = new double[n];
//	 for (unsigned int i = 0; i < m; i++) {
//		 for (unsigned int j = 0; j < nb; j++) {
//			 for (int k = 0; k < n; k++)
//				 C[ i * nb + j ] += A [ i * n + k ] * B.A[ k * nb + j ];  //C_{ij} = A_{ik} * B_{kj}
//		 }
//	 }
//	 return matrix(C, m, nb);
// }

  matrix matrix::changeHandness() const{
 	 vector<double> A2 (m * n );
 	 for( unsigned int i = 0; i < m; i++){
 		A2 [ n*i ] = -A [ n * i ];
 		A2 [ n*i + 1 ] = A [ n*i + 1 ];
 		A2 [ n*i  + 2 ] = A [ n*i + 2] ;
// 		 A2[i][0] = -A[i][0];
// 		 A2[i][1] =  A[i][1];
// 		 A2[i][2] =  A[i][2];
 	 }
 	 return matrix(A2, m, n);
  }


//  vector<double> matrix::times (const vector<double>& B, int len) const {
//  	 if ( len != n) {
//  		 cerr<<"matrix inner dimensions must agree."<<endl;
//  		 exit(0);
//  	 }
//  	 vector<double> vec(m, 0.0);
//  	 for (int j = 0; j < m; j++) {
//  		 for (int k = 0; k < n; k++) {
//  			 vec [ j ]+= A [ j * n + k ] * B [ k ]; //A[j][k] * B[k];
//  		 }
//  	 }
//  	 return vec;
//   }



/** One norm
 @return    maximum column sum.
 */
 double  matrix::norm1 () {
	 double f = 0;
	 for (unsigned int j = 0; j < n; j++) {
		 double s = 0;
		 for (unsigned int i = 0; i < m; i++) {
			 s += fabs(A [ i * n + j ] ) ; //A[i][j]);
		 }
		 f = max(f, s);
	 }
	 return f;
 }

 double matrix::trace () {
	 double t = 0;
	 for (unsigned int i = 0; i < min(m,n); i++) {
		 t += A.at (i * n + i);
	 }
	 return t;
 }

 double matrix::det (const matrix& V) const{
	 int m = V.m;
	 int n = V.n;
	 if (m != 3 || n != 3){
		 cerr<<" Det Only for 3 x 3 matrix "<<endl;
		 exit(1);
	 }
	 //double **arrV = V.getArray();
	 return (V.A [ 0 ] * V.A [ 4 ] * V.A [ 8 ] + V.A [ 3] * V.A [ 7] * V.A [ 2]  + V.A [ 6] * V.A [ 5] * V.A [ 1]
			 - V.A [ 2]  * V.A [ 4] *  V.A [ 6] -V.A [ 5]  * V.A [ 7] *  V.A [ 0] -V.A [ 8]  * V.A [ 3] *  V.A [ 1 ]  );
 }

 /** Generate identity matrix
 @param m    Number of rows.
 @param n    Number of colums.
 @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
 */

 void matrix::print() const{
	 for(unsigned int i = 0; i<m; i++){
		 for(unsigned int j = 0; j<n; j++)
			 printf("%18.10f", A [ i * n + j ]);
		 cout<<endl;
	 }
	 cout<<endl;
 }

/**
 * This one must be corrected for stopping the memory leaking !!!
 */
 float** similarityTransformZ(float **a,  float sinT, float cosT){
	 // $ B = R_z(\theta} A $
	 const unsigned int m = 3;
	 const unsigned int n = 3;
	 float **b = new float*[m];
	 float **c = new float*[m];
	 for(unsigned int j = 0; j < m; j++){
		 b[j] = new float[n];
		 c[j] = new float[n];
	 }
	 b[0][0] = a[0][0] * cosT + a[1][0] * sinT;
	 b[0][1] = a[0][1] * cosT + a[1][1] * sinT;
	 b[0][2] = a[0][2] * cosT + a[1][2] * sinT;
	 b[1][0] = - a[0][0] * sinT + a[1][0] * cosT;
	 b[1][1] = - a[0][1] * sinT + a[1][1] * cosT;
	 b[1][2] = - a[0][2] * sinT + a[1][2] * cosT;
	 b[2][0] = a[2][0];
	 b[2][1] = a[2][1];
	 b[2][2] = a[2][2];
	 // the matrix $C = A^{-1} B$
	 c[0][0] = a[0][0] * b[0][0] + a[1][0] * b[1][0] +  a[2][0] * b[2][0];
	 c[0][1] = a[0][0] * b[0][1] + a[1][0] * b[1][1] +  a[2][0] * b[2][1];
	 c[0][2] = a[0][0] * b[0][2] + a[1][0] * b[1][2] +  a[2][0] * b[2][2];
	 c[1][0] = a[0][1] * b[0][0] + a[1][1] * b[1][0] +  a[2][1] * b[2][0];
	 c[1][1] = a[0][1] * b[0][1] + a[1][1] * b[1][1] +  a[2][1] * b[2][1];
	 c[1][2] = a[0][1] * b[0][2] + a[1][1] * b[1][2] +  a[2][1] * b[2][2];
	 c[2][0] = a[0][2] * b[0][0] + a[1][2] * b[1][0] +  a[2][2] * b[2][0];
	 c[2][1] = a[0][2] * b[0][1] + a[1][2] * b[1][1] +  a[2][2] * b[2][1];
	 c[2][2] = a[0][2] * b[0][2] + a[1][2] * b[1][2] +  a[2][2] * b[2][2];
//	 for(int j = 0; j < m; j++)
//		 delete[]  b[j];
//	 delete [] b;
	 return c;
 }
//   /** Print the matrix to stdout.   Line the elements up in columns
//   * with a Fortran-like 'Fw.d' style format.
// @param w    Column width.
// @param d    Number of digits after the decimal.
// */
//  public void print (int w, int d) {
//    print(new PrintWriter(System.out,true),w,d); }
//
// /** Print the matrix to the output stream.   Line the elements up in
//   * columns with a Fortran-like 'Fw.d' style format.
// @param output Output stream.
// @param w      Column width.
// @param d      Number of digits after the decimal.
// */
//
// public void print (PrintWriter output, int w, int d) {
//    DecimalFormat format = new DecimalFormat();
//    format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
//    format.setMinimumIntegerDigits(1);
//    format.setMaximumFractionDigits(d);
//    format.setMinimumFractionDigits(d);
//    format.setGroupingUsed(false);
//    print(output,format,w+2);
// }
//
//  public void print (NumberFormat format, int width) {
//	print(new PrintWriter(System.out,true),format,width); }

 // DecimalFormat is a little disappointing coming from Fortran or C's printf.
 // Since it doesn't pad on the left, the elements will come out different
 // widths.  Consequently, we'll pass the desired column width in as an
 // argument and do the extra padding ourselves.

 /** Print the matrix to the output stream.  Line the elements up in columns.
   * Use the format object, and right justify within columns of width
   * characters.
   * Note that is the matrix is to be read back in, you probably will want
   * to use a NumberFormat that is set to US Locale.
 @param output the output stream.
 @param format A formatting object to format the matrix elements
 @param width  Column width.
 @see java.text.DecimalFormat#setDecimalFormatSymbols
 */
//
// public void print (PrintWriter output, NumberFormat format, int width) {
//    output.println();  // start on new line.
//    for (unsigned int i = 0; i < m; i++) {
//       for (unsigned int j = 0; j < n; j++) {
//          String s = format.format(A[i][j]); // format the number
//          int padding = max(1,width-s.length()); // At _least_ 1 space
//          for (int k = 0; k < padding; k++)
//             output.print(' ');
//          output.print(s);
//       }
//       output.println();
//    }
//    output.println();   // end with blank line.
// }
  /** Singular Value Decomposition
 @return     SingularValueDecomposition
 @see SingularValueDecomposition
 */
//  public SingularValueDecomposition svd () {
//    return new SingularValueDecomposition(this);
// }

 matrix matrix::rotationMat(double theta, string axis) {
 	 vector<double> A(9, 0.0); // 3*3 - 1
 	 if (axis=="+x"){
 		 A [ 0] =  1.0;
 		 A [ 4 ]= cos(theta);
 		 A [ 5] =  sin(theta);
 		 A [ 7] =  -sin(theta);
 		 A [ 8] =  cos(theta);
 	 } else if (axis=="-x"){
 		 A [ 0] =  1.0;
 		 A [ 4] =  cos(theta);
 		 A [ 5] =  -sin(theta);
 		 A [ 7] =  sin(theta);
 		 A [ 8] =  cos(theta);
 	 }else if (axis=="+y"){
 		 A [ 0] =  cos(theta);
 		 A [ 2] =  -sin(theta);
 		 A [ 4] =  1.0;
 		 A [ 6] =  sin(theta);
 		 A [ 8] =  cos(theta);
 	 }else if (axis=="-y"){
 		 A [ 0] =  cos(theta);
 		 A [ 2] =  sin(theta);
 		 A [ 4] =  1.0;
 		 A [ 6] =  -sin(theta);
 		 A [ 8] =  cos(theta);
 	 }else if (axis=="+z"){
 		 A [ 0] =  cos(theta);
 		 A [ 1] =  sin(theta);
 		 A [ 3] =  -sin(theta);
 		 A [ 4] =  cos(theta);
 		 A [ 8] =  1.0;
 	 }else if (axis=="-z"){
 		 A [ 0] =  cos(theta);
 		 A [ 1] =  -sin(theta);
 		 A [ 3] =  sin(theta);
 		 A [ 4] =  cos(theta);
 		 A [ 8] =  1.0;
 	 }
 	 return matrix(A, 3, 3);
  }

 matrix matrix::rotationMatCS(double cosTheta, string axis) {
  	 vector<double> A(9, 0.0); // 3*3 - 1
  	 if (axis=="+x"){
  		 A [ 0] =  1.0;
  		 A [ 4 ]= cosTheta;
  		 A [ 5] =  sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 7] =  -A[5];
  		 A [ 8] =  cosTheta;
  	 } else if (axis=="-x"){
  		 A [ 0] =  1.0;
  		 A [ 4] =  cosTheta;
  		 A [ 5] =  - sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 7] =  -A[5];
  		 A [ 8] =  cosTheta;
  	 }else if (axis=="+y"){
  		 A [ 0] =  cosTheta;
  		 A [ 2] =  -sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 4] =  1.0;
  		 A [ 6] =  - A [ 2];
  		 A [ 8] =  cosTheta;
  	 }else if (axis=="-y"){
  		 A [ 0] =  cosTheta;
  		 A [ 2] =  sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 4] =  1.0;
  		 A [ 6] =  - A [ 2];
  		 A [ 8] =  cosTheta;
  	 }else if (axis=="+z"){
  		 A [ 0] =  cosTheta;
  		 A [ 1] =  sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 3] =  -A[1];
  		 A [ 4] =  cosTheta;
  		 A [ 8] =  1.0;
  	 }else if (axis=="-z"){
  		 A [ 0] =  cosTheta;
  		 A [ 1] =  -sqrt(1.0 - cosTheta* cosTheta );
  		 A [ 3] =  -A[1];
  		 A [ 4] =  cosTheta;
  		 A [ 8] =  1.0;
  	 }
  	 return matrix(A, 3, 3);
   }


