#include<functional>
#include<iostream>
#include<vector>
#include <cstdlib>
#include <cstdio>

//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_multifit_nlin.h>

#include "../commons/point3D.h"

struct data1d {
	size_t n;
	double *t;
	double *y;
	double *sigma;
};

struct data2d {
	size_t n;
	double *t1;
	double *t2;
	double *z;
	double *sigma;
};

int gaussb_f (const gsl_vector * x, void *data,	gsl_vector * f) {
	size_t n = ((struct data1d *)data)->n;
	double *t = ( (struct data1d *)dat )->t;
	double *y = ((struct data1d *)data)->y;
	double *sigma4Fit = ((struct data1d *) data)->sigma;

	double A = gsl_vector_get (x, 0);
	double mu = gsl_vector_get (x, 1);
	double sigma  = gsl_vector_get (x, 2);
	double sigmaSquare = sigma * sigma;
	double b = gsl_vector_get (x, 3);

	size_t i;
	double Yi;

	/* Model $Yi = b +  A \exp( - \frac{ ( t_i - \mu )^2 } { 4\pi  \sigma^2} ) */
	for (i = 0; i < n; i++) 	{
		Yi = A * exp ( - ( t[i] - mu ) * ( t[i] - mu )  / sigmaSquare ) + b;
		gsl_vector_set (f, i, (Yi - y[i] ) / sigma4Fit[i] );
	}

	return GSL_SUCCESS;
}

int gaussb_df (const gsl_vector * x, void *data,	gsl_matrix * J) {
	size_t n = ((struct data1d *)data)->n;
	double *t = ( (struct data1d *)dat )->t;
	double *sigma = ((struct data1d *) data)->sigma;

	double A = gsl_vector_get (x, 0);
	double mu = gsl_vector_get (x, 1);
	double sigma  = gsl_vector_get (x, 2);
	double sigmaS = sigma * sigma;
	double b = gsl_vector_get (x, 3);

	size_t i;

	/**
	 * Ignoring 4 \pi
	 * y^'_A = \exp { - \frac{ ( t_i - \mu)^ 2 }{ sigma^2 }
	 * y^'_{\mu} = A \frac{ 2 (t_i - \mu)}{ sigma^2}  \exp { - \frac{ ( t_i - \mu)^ 2 }{ sigma^2 }
	 * y^'_{\sigma^2} = A \frac{ (t_i - \mu)^2}{ sigma^4}  \exp { - \frac{ ( t_i - \mu)^ 2 }{ sigma^2 }
	 * y^'_b = 1.0
	 * let
	 * Let c = \frac{t_i - \mu}{sigma^2},  e = \exp { -c (t_i - \mu) }, then
	 *  y^'_A = e;
	 *  y^'_{\mu}  = 2 A ce
	 *   y^'_{\sigma^2} = A c^2 e
	 * The real function for taking derivatives is
	 * where fi = (y - yi)/sigma[i],
	 */
	double s = 0.0, e = 0.0, c = 0.0;
	for (i = 0; i < n; i++) 	{
		s = sigma[i];
		c = ( t[i] - mu )  / sigmaS;
		e =  exp( - c * ( t[i] - mu ) ) ;
		gsl_matrix_set (J, i, 0, e / s);
		gsl_matrix_set (J, i, 1, 2.0 * A * c * e / s);
		gsl_matrix_set (J, i, 2, A * c * c * e / s);
		gsl_matrix_set (J, i, 3, 1/s);

	}
	return GSL_SUCCESS;
}

int gaussb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	gaussb_f (x, data, f);
	gaussb_df (x, data, J);

	return GSL_SUCCESS;
}

int expb_f (const gsl_vector * x, void *data,  gsl_vector * f) {
  size_t n = ((struct data1d *)data)->n;
  double *y = ((struct data1d *)data)->y;
  double *sigma = ((struct data1d *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++) {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set ( f, i, ( Yi - y[i] ) / sigma[i]);
    }

  return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J ) {
	size_t n = ((struct data1d *)data)->n;
	double *sigma = ((struct data1d *) data)->sigma;

	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);

	size_t i;

	for (i = 0; i < n; i++)  {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double t = i;
		double s = sigma[i];
		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, e/s);
		gsl_matrix_set (J, i, 1, -t * A * e/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
  expb_f (x, data, f);
  expb_df (x, data, J);
  return GSL_SUCCESS;
}



/**
 *
 * #include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "expfit.c"

#define N 40

void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = N;
  const size_t p = 3;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double y[N], sigma[N];
  struct data d = { n, y, sigma};
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.0, 0.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  /  This is the data to be fitted

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t)
                 + gsl_ran_gaussian (r, 0.1);
      sigma[i] = 0.1;
      printf ("data: %u %g %g\n", i, y[i], sigma[i]);
    };

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  {
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  }

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
  return 0;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
}
 *
 */


inline int lorentz_f ( const gsl_vector *x, void *dat, gsl_vector *f ) {

	size_t n = ( (struct data1d *)dat )->n;
	double *t = ( (struct data1d *)dat )->t;
	double *y = ( (struct data1d *)dat )->y;
	double *sigma = ( (struct data1d *) dat )->sigma;

	double A = gsl_vector_get (x, 0);
	double mu = gsl_vector_get (x, 1);
	double gamma = gsl_vector_get (x, 2);

	size_t i;

	for (i = 0; i < n; i++) 	{
		/* Model $Yi = \frac{ A \gamma^2} { ( t_i - \mu )^2 + \gamma^2 */
//		double t = points[i][0];
		double Yi = ( A * gamma * gamma ) / (  (t[i] - mu ) * (t[i] - mu)  +  gamma * gamma);
		gsl_vector_set (f, i, (Yi - y[i] ) / sigma[i]);
	}

	return GSL_SUCCESS;
}

inline int lorentz_df (const gsl_vector *x, void *dat,	 gsl_matrix *J) {
	size_t n = ( (struct data1d *)dat )->n;
	double *t = ( (struct data1d *)dat )->t;
	double *sigma = ( (struct data1d *) dat )->sigma;

	double A = gsl_vector_get (x, 0);
	double mu = gsl_vector_get (x, 1);
	double gamma = gsl_vector_get (x, 2);

	size_t i;

	for ( i = 0; i < n; i++) 	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double s = sigma[i];
		double a = ( t[i] - mu ) * ( t[i] - mu )  +  gamma * gamma;
		//		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, gamma*gamma / a / s);
		gsl_matrix_set (J, i, 1, 2.0 * A * gamma * gamma * ( t[i] - mu ) / a / a / s);
		gsl_matrix_set (J, i, 2, 2.0 * A * gamma *  ( t[i] - mu ) *  ( t[i] - mu )   / a / a / s);
	}
	return GSL_SUCCESS;
}

inline int lorentz_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	lorentz_f (x, data, f);
	lorentz_df (x, data, J);

	return GSL_SUCCESS;
}


//void print_state (size_t iter, gsl_multifit_fdfsolver * s);
inline void print_state (size_t N, size_t iter, gsl_multifit_fdfsolver * s) {
	cout <<"iter "<< iter<<": "<<endl;
	for (int i = 0; i < N; i++){
		printf("%15.8f%15.8f%15.8f\n",
				gsl_vector_get (s->x, 3*i ), gsl_vector_get (s->x, i*3+1), gsl_vector_get (s->x, i*3+2) );
	}
	cout<<"|f(x)| = "<<gsl_blas_dnrm2 (s->f)<<endl;
//	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
//			"|f(x)| = %g\n",
//			iter,
//			gsl_vector_get (s->x, 0),
//			gsl_vector_get (s->x, 1),
//			gsl_vector_get (s->x, 2),
//			gsl_blas_dnrm2 (s->f));
}

inline void print_state (size_t iter, gsl_multifit_fdfsolver * s) {
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f " "|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0),
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2),
			gsl_blas_dnrm2 (s->f)  );
}


inline int lorentzFit (const size_t p, const vector<pair<float, float>*> & onePkVec,
		double *x_init) {
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	const size_t n = onePkVec.size();
	unsigned int i, iter = 0;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double t[n], y[n], sigma[n];

	for ( i = 0; i < n; i++){
		t[i] = onePkVec.at( i )->first;
		y[i] = onePkVec.at( i )->second;
		sigma[i] = 1.0;
		printf ("data: %g %g %g\n", t[i], y[i], sigma[i]);
	}
	struct data1d d = { n, t, y, sigma};
	gsl_multifit_function_fdf f;
//	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &lorentz_f;
	f.df = &lorentz_df;
	f.fdf = &lorentz_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set ( s, &f, &x.vector );

	print_state ( iter, s );

	do	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

//		printf ("status = %s\n", gsl_strerror (status));
//
		print_state (iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-6, 1e-6);
	} while (status == GSL_CONTINUE && iter < 5000);

	gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

		printf ("Intensity      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("mu = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("gamma      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	}

	printf ("status = %s\n", gsl_strerror (status));
	x_init[0] = gsl_vector_get (s->x, 0);
	x_init[1] = gsl_vector_get (s->x, 1);
	x_init[2] = gsl_vector_get (s->x, 2);
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return 0;
}
//

inline int lorentz_2f (const gsl_vector *x, void *data, gsl_vector *f ) {

	size_t n = ((struct data1d *)data)->n;
	double *t = ((struct data1d *)data)->t;
	double *y = ((struct data1d *)data)->y;
	double *sigma = ((struct data1d *) data)->sigma;

	double A1 = gsl_vector_get (x, 0);
	double A2 = gsl_vector_get (x, 1);
	double mu1 = gsl_vector_get (x, 2);
	double mu2 = gsl_vector_get (x, 3);
	double gamma1 = gsl_vector_get (x, 4);
	double gamma2 = gsl_vector_get (x, 5);

	size_t i;

	for (i = 0; i < n; i++) 	{
		/* Model $Yi = \frac{ A \gamma^2} { ( t_i - \mu )^2 + \gamma^2 */
//		double t = points[i][0];
		double Yi =  A1 * gamma1 * gamma1  / (  ( t[i] - mu1 ) * ( t[i] - mu1 )  +  gamma1 * gamma1 )
							   +  A2 * gamma2 * gamma2  / (  ( t[i] - mu2 ) * ( t[i] - mu2 )  +  gamma2 * gamma2 );
		gsl_vector_set (f, i, (Yi - y[i] ) / sigma[i]);
	}

	return GSL_SUCCESS;
}

inline int lorentz_2df (const gsl_vector * x, void *data,	gsl_matrix * J) {
	size_t n = ((struct data1d *)data)->n;
	double *t = ((struct data1d *)data)->t;
	double *sigma = ((struct data1d *) data)->sigma;

	double x0 = gsl_vector_get (x, 0); //A1
	double x1 = gsl_vector_get (x, 1);  //A2
	double x2 = gsl_vector_get (x, 2);  //mu1
	double x3 = gsl_vector_get (x, 3);  //mu2
	double x4 = gsl_vector_get (x, 4);  //gamma1
	double x5 = gsl_vector_get (x, 5);  //gamma2

	size_t i;

	for ( i = 0; i < n; i++) 	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double s = sigma[i];
		double b1 = (t[i] - x2);
		double a1 = b1 * b1  +  x4 * x4;
		double b2 = (t[i] - x3);
		double a2 = b2 * b2  +  x5 * x5;
		//		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, x4 * x4  / a1 / s);
		gsl_matrix_set (J, i, 1, x5 * x5  / a2 / s);
		gsl_matrix_set (J, i, 2, 2.0 * x0 * x4 * x4 * b1  / a1 / a1 / s);
		gsl_matrix_set (J, i, 3, 2.0 * x1 * x5 * x5 * b2 / a2  / a2 / s);
		gsl_matrix_set (J, i, 4, 2.0 * x0 * x4 *  b1 * b1  / a1 / a1 / s);
		gsl_matrix_set (J, i, 5, 2.0 * x1 * x5 *  b2 * b2  / a2 / a2 / s);
	}
	return GSL_SUCCESS;
}

inline int lorentz_2fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	lorentz_2f (x, data, f);
	lorentz_2df (x, data, J);

	return GSL_SUCCESS;
}


inline int lorentzFit2Peaks (const size_t p, const vector<pair<float, float>*>& onePkVec, double *x_init) {

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	const size_t n = onePkVec.size();
	unsigned int i, iter = 0;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double t[n], y[n], sigma[n];

	for ( i = 0; i < n; i++){
		t[i] = onePkVec.at( i )->first;
		y[i] = onePkVec.at( i )->second;
		sigma[i] = 1.0;
		printf ("data: %g %g %g\n", t[i], y[i], sigma[i]);
	}
	struct data1d d = { n, t, y, sigma};
	gsl_multifit_function_fdf f;
//	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &lorentz_2f;
	f.df = &lorentz_2df;
	f.fdf = &lorentz_2fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	print_state (iter, s);

	do	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		printf ("status = %s\n", gsl_strerror (status));
		print_state (iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-5, 1e-5);
	} while (status == GSL_CONTINUE && iter < 2000);

	gsl_multifit_covar ( s->J, 0.0, covar );

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

		printf ("Intensity1      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("mu1 = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
		printf ("gamma1     = %.5f +/- %.5f\n", FIT(4), c*ERR(4));

		printf ("Intensity2      = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("mu2 = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
		printf ("gamma2      = %.5f +/- %.5f\n", FIT(5), c*ERR(5));
	}

	printf ("status = %s\n", gsl_strerror (status));
	x_init[0] = gsl_vector_get (s->x, 0);
	x_init[1] = gsl_vector_get (s->x, 1);
	x_init[2] = gsl_vector_get (s->x, 2);
	x_init[3] = gsl_vector_get (s->x, 3);
	x_init[4] = gsl_vector_get (s->x, 4);
	x_init[5] = gsl_vector_get (s->x, 5);
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return 0;
}


inline int lorentz_Nf ( const gsl_vector *x, void *data, gsl_vector *f ) {

	size_t n = ((struct data1d *)data)->n;
	double *t = ((struct data1d *)data)->t;
	double *y = ((struct data1d *)data)->y;
	double *sigma = ((struct data1d *) data)->sigma;

	const unsigned int N = (x->size) / 3;
//	cout<<"N = "<<N<<endl;
	double A[N];
	double mu[N];
	double gamma[N];
	for (unsigned int i =0; i < N; i++){
		A[i] = gsl_vector_get (x, 3 * i );
		mu[i] = gsl_vector_get (x, 3* i  + 1 );
		gamma[i] = gsl_vector_get (x, 3* i  + 2 );
		cout<<A[i]<<"  "<<mu[i]<<"  "<<gamma[i]<<endl;
//		double A1 = gsl_vector_get (x, 0);
//		double A2 = gsl_vector_get (x, 1);
//		double mu1 = gsl_vector_get (x, 2);
//		double mu2 = gsl_vector_get (x, 3);
//		double gamma1 = gsl_vector_get (x, 4);
//		double gamma2 = gsl_vector_get (x, 5);
	}

	size_t i;
	double Yi = 0.0;
	for ( i = 0; i < n; i++ ) 	{
		/* Model $Yi = \frac{ A \gamma^2} { ( t_i - \mu )^2 + \gamma^2 */
//		double t = points[i][0];
		Yi = 0.0;
		for( unsigned j = 0; j < N; j++) {
			Yi += A[j] * gamma[j] * gamma[j]  / ( ( t[i] - mu[j] ) * ( t[i] - mu[j] ) +gamma[j] * gamma[j] ) ;
//							   +  A2 * gamma2 * gamma2  / (  ( t[i] - mu2 ) * ( t[i] - mu2 )  +  gamma2 * gamma2 );
		}
		gsl_vector_set ( f, i, (Yi - y[i] ) / sigma[i] );
	}

	return GSL_SUCCESS;
}

inline int lorentz_Ndf (const gsl_vector *x, void *data,	gsl_matrix * J) {
	size_t n = ((struct data1d *)data)->n;
	double *t = ((struct data1d *)data)->t;
	double *sigma = ((struct data1d *) data)->sigma;

	const unsigned int N =  ( x->size ) / 3;  //there are three parameters
	double A[N];
	double mu[N];
	double gamma[N];
	for (unsigned int j =0; j < N; j++){
		A[j] = gsl_vector_get (x, 3 * j );
		mu[j] = gsl_vector_get (x, 3*j + 1 );
		gamma[j] = gsl_vector_get (x, 3* j  + 2 );

//		A[j] = gsl_vector_get (x, j * N );
//		mu[j] = gsl_vector_get (x, j * N + 1);
//		gamma[j] = gsl_vector_get (x, j * N + 2);

		//		double A1 = gsl_vector_get (x, 0);
		//		double A2 = gsl_vector_get (x, 1);
		//		double mu1 = gsl_vector_get (x, 2);
		//		double mu2 = gsl_vector_get (x, 3);
		//		double gamma1 = gsl_vector_get (x, 4);
		//		double gamma2 = gsl_vector_get (x, 5);
	}

//	double x0 = gsl_vector_get (x, 0); //A1
//	double x1 = gsl_vector_get (x, 1);  //A2
//	double x2 = gsl_vector_get (x, 2);  //mu1
//	double x3 = gsl_vector_get (x, 3);  //mu2
//	double x4 = gsl_vector_get (x, 4);  //gamma1
//	double x5 = gsl_vector_get (x, 5);  //gamma2

	size_t i;

	double a[N];
	double b[N];
	for ( i = 0; i < n; i++) 	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double s = sigma[i];
		for ( unsigned int j = 0; j < N; j++) {
			b[j] = t[i] - mu[ j ];
			a[j] = b[j] * b[j] + gamma[ j ] * gamma[ j ];
			gsl_matrix_set (J, i, j * 3 + 0, gamma[j] * gamma[j]  / a[ j ] / s);
			gsl_matrix_set (J, i, j * 3 + 1,  2.0 * A[j] * gamma[j] * gamma[j] * b[j]  / a[j] / a[j] / s);
			gsl_matrix_set (J, i, j * 3 + 2,  2.0 * A[j] * gamma[j] * b[j] * b[j]  / a[j] / a[j] / s);
//		   double b1 = ( t[i] - x2);
//			double a1 = b1 * b1  +  x4 * x4;
//			double b2 = ( t[i] - x3);
//			double a2 = b2 * b2  +  x5 * x5;
//			//		double e = exp(-lambda * t);
//			gsl_matrix_set (J, i, 0, x4 * x4  / a1 / s);
//			gsl_matrix_set (J, i, 1, x5 * x5  / a2 / s);
//			gsl_matrix_set (J, i, 2, 2.0 * x0 * x4 * x4 * b1  / a1 / a1 / s);
//			gsl_matrix_set (J, i, 3, 2.0 * x1 * x5 * x5 * b2 / a2  / a2 / s);
//			gsl_matrix_set (J, i, 4, 2.0 * x0 * x4 *  b1 * b1  / a1 / a1 / s);
//			gsl_matrix_set (J, i, 5, 2.0 * x1 * x5 *  b2 * b2  / a2 / a2 / s);
		}
	}
	return GSL_SUCCESS;
}

inline int lorentz_Nfdf (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {
	lorentz_Nf ( x, data, f );
	lorentz_Ndf ( x, data, J );

	return GSL_SUCCESS;
}


inline int lorentzFit2NPeaks (const size_t p, const vector<pair<float, float>*>& onePkVec, double *x_init) {

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	const size_t n = onePkVec.size();
	unsigned int i, iter = 0;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double t[n], y[n], sigma[n];

	for ( i = 0; i < n; i++){
		t[i] = onePkVec.at( i )->first;
		y[i] = onePkVec.at( i )->second;
		sigma[i] = 1.0;
		printf ("data: %g %g %g\n", t[i], y[i], sigma[i] );
	}
	struct data1d d = { n, t, y, sigma};
	gsl_multifit_function_fdf f;
//	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &lorentz_Nf;
	f.df = &lorentz_Ndf;
	f.fdf = &lorentz_Nfdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;  //gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	const int N = p / 3;
	print_state (N, iter, s);

	do	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		printf ("status = %s\n", gsl_strerror (status));
		print_state ( N, iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x,1e-7, 1e-7);
	} while (status == GSL_CONTINUE && iter < 5000);

	gsl_multifit_covar ( s->J, 0.0, covar );

	cout<<endl;
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		for ( int i = 0; i < N; i++){
			cout<<"Peak: "<<i<<endl;
			printf ("Intensity      = %.5f +/- %.5f\n", FIT( i * 3 ), c*ERR( i * 3 ));
			printf ("mu = %.5f +/- %.5f\n", FIT( i * 3 + 1), c*ERR( i * 3 + 1));
			printf ("gamma     = %.5f +/- %.5f\n\n", FIT( i * 3 + 2), c*ERR( i * 3 + 2));
		}
	}

	printf ("status = %s\n", gsl_strerror (status));
	for ( unsigned int i =0; i<p; i++)
		x_init[i] = gsl_vector_get (s->x, i);

//	x_init[0] = gsl_vector_get (s->x, 0);
//	x_init[1] = gsl_vector_get (s->x, 1);
//	x_init[2] = gsl_vector_get (s->x, 2);
//	x_init[3] = gsl_vector_get (s->x, 3);
//	x_init[4] = gsl_vector_get (s->x, 4);
//	x_init[5] = gsl_vector_get (s->x, 5);

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return 0;
}


/**
 * Fit to 2 dimensional Lorentz-Cauchy line shape
 */
inline int twoDimLorentz_f (const gsl_vector *x, void *data, gsl_vector *f ) {

	size_t n = ((struct data2d *)data)->n;
	double *t1 = ((struct data2d *)data)->t1;
	double *t2 = ((struct data2d *)data)->t2;
	double *z = ((struct data2d *)data)->z;
	double *sigma = ((struct data2d *) data)->sigma;

	double A = gsl_vector_get (x, 0);
	double mu1 = gsl_vector_get (x, 1);
	double mu2 = gsl_vector_get (x, 2);
	double gamma1 = gsl_vector_get (x, 3);
	double gamma2 = gsl_vector_get (x, 4);

	size_t i;

	for (i = 0; i < n; i++) 	{
		/* Model $Yi = \frac{ A \gamma^2} { ( t_i - \mu )^2 + \gamma^2 */
//		double t = points[i][0];
		double Zi =  A * gamma1 * gamma1 * gamma2* gamma2 / (  ( t1[i] - mu1 ) * ( t1[i] - mu1 )  +  gamma1 * gamma1 )
							   /   (  ( t2[i] - mu2 ) * ( t2[i] - mu2 )  +  gamma2 * gamma2 );
		gsl_vector_set (f, i, (Zi - z[i] ) / sigma[i]);
	}

	return GSL_SUCCESS;
}

inline int twoDimLorentz_df (const gsl_vector *x, void *data,	gsl_matrix * J) {
	size_t n = ((struct data2d *)data)->n;
	double *t1 = ((struct data2d *)data)->t1;
	double *t2 = ((struct data2d *)data)->t2;
	double *sigma = ((struct data2d *) data)->sigma;

	double x0 = gsl_vector_get (x, 0); //A
	double x1 = gsl_vector_get (x, 1);  //mu1
	double x2 = gsl_vector_get (x, 2);  //mu2
	double x3 = gsl_vector_get (x, 3);  //gamma1
	double x4 = gsl_vector_get (x, 4);  //gamma2

	size_t i;

	for ( i = 0; i < n; i++) 	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double s = sigma[i];
		double xx3 = x3 * x3;
		double xx4 = x4 * x4;
		double b1 =  t1[i] - x1 ;
		double a1 = b1 * b1  +  xx3;
		double b2 =  t2[i] - x2 ;
		double a2 = b2 * b2  +  xx4;
		double c = x0 * xx3 * xx4;
		//		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, xx3 * xx4  / a1 / a2 / s);
		gsl_matrix_set (J, i, 1, 2.0 * b1 * c  / a2 / a1 / a1 / s );
		gsl_matrix_set (J, i, 2, 2.0 * b2 * c  / a1 / a2 / a2 / s);
		gsl_matrix_set (J, i, 3, 2.0 * x0 * x3 * xx4 * b1 * b1  / a1 / a1 / a2 / s);
		gsl_matrix_set (J, i, 4, 2.0 * x0 * x4 * xx3 * b2 * b2  / a1 / a2 / a2 / s);
	}
	return GSL_SUCCESS;
}

inline int twoDimLorentz_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	twoDimLorentz_f (x, data, f);
	twoDimLorentz_df (x, data, J);

	return GSL_SUCCESS;
}
//
/**
 * Which kind of data we want to pass?
 * a point3D data by extracting all the points above a threshold, adding the contour at the threshold
 */
inline int twoDimLorentzFit (const size_t p, const vector<point3D<double, double, double>* > &onePkVec, double *x_init) {

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	const size_t n = onePkVec.size();
	unsigned int i, iter = 0;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double t1[n], t2[n], y[n], sigma[n];

	for ( i = 0; i < n; i++){
		t1[i] = onePkVec.at( i )->getX();
		t2[i] = onePkVec.at( i )->getY();
		y[i] =  onePkVec.at( i )->getZ();
		sigma[i] = 1.0;
		printf ("data: %8g%8g %8g %g\n", t1[i], t2[i], y[i], sigma[i]);
	}
	struct data2d d = { n, t1, t2, y, sigma};

	gsl_multifit_function_fdf f;
	//	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &twoDimLorentz_f;
	f.df = & twoDimLorentz_df;
	f.fdf = & twoDimLorentz_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	print_state (iter, s);

	do	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		printf ("status = %s\n", gsl_strerror (status));
		print_state (iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-7, 1e-7);
	} while (status == GSL_CONTINUE && iter < 5000);

	gsl_multifit_covar ( s->J, 0.0, covar );

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

		printf ("Intensity      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("mu1 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("mu2 = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
		printf ("gamma1     = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
		printf ("gamma2      = %.5f +/- %.5f\n", FIT(4), c*ERR(4));
	}

	printf ("status = %s\n", gsl_strerror (status));
	x_init[0] = gsl_vector_get (s->x, 0);
	x_init[1] = gsl_vector_get (s->x, 1);
	x_init[2] = gsl_vector_get (s->x, 2);
	x_init[3] = gsl_vector_get (s->x, 3);
	x_init[4] = gsl_vector_get (s->x, 4);
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return status;
//	cout<<status<<endl;
//	if (gsl_strerror (status) =="success")
//		return 0;
//	else
//		return -1;
}

/**
 * Fit to 2 dimensional Lorentz-Cauchy line shape
 */
inline int twoDimLorentz2_f (const gsl_vector *x, void *data, gsl_vector *f ) {

	size_t n = ((struct data2d *)data)->n;
	double *t1 = ((struct data2d *)data)->t1;
	double *t2 = ((struct data2d *)data)->t2;
	double *z = ((struct data2d *)data)->z;
	double *sigma = ((struct data2d *) data)->sigma;

	double A1 = gsl_vector_get (x, 0);
	double mu1 = gsl_vector_get (x, 1);
	double mu2 = gsl_vector_get (x, 2);
	double gamma1 = gsl_vector_get (x, 3);
	double gamma2 = gsl_vector_get (x, 4);

	double A2    = gsl_vector_get (x, 5);
	double mu3 = gsl_vector_get (x, 6);
	double mu4 = gsl_vector_get (x, 7);
	double gamma3 = gsl_vector_get (x, 8);
	double gamma4 = gsl_vector_get (x, 9);

	size_t i;

	for (i = 0; i < n; i++) 	{
		/* Model $Yi = \frac{ A \gamma^2} { ( t_i - \mu )^2 + \gamma^2 */
//		double t = points[i][0];
		double Zi =  A1 * gamma1 * gamma1 * gamma2* gamma2 / (  ( t1[i] - mu1 ) * ( t1[i] - mu1 )  +  gamma1 * gamma1 )
							        /   (  ( t2[i] - mu2 ) * ( t2[i] - mu2 )  +  gamma2 * gamma2 )
							   + A2 * gamma3 * gamma3 * gamma4* gamma4 / (  ( t1[i] - mu3) * ( t1[i] - mu3 )  +  gamma3 * gamma3 )
							  							   /   (  ( t2[i] - mu4 ) * ( t2[i] - mu4 )  +  gamma4 * gamma4)
							   ;
		gsl_vector_set (f, i, (Zi - z[i] ) / sigma[i]);
	}

	return GSL_SUCCESS;
}

inline int twoDimLorentz2_df (const gsl_vector *x, void *data,	gsl_matrix * J) {
	size_t n = ((struct data2d *)data)->n;
	double *t1 = ((struct data2d *)data)->t1;
	double *t2 = ((struct data2d *)data)->t2;
	double *sigma = ((struct data2d *) data)->sigma;

	double x0 = gsl_vector_get (x, 0); //A1
	double x1 = gsl_vector_get (x, 1);  //mu1
	double x2 = gsl_vector_get (x, 2);  //mu2
	double x3 = gsl_vector_get (x, 3);  //gamma1
	double x4 = gsl_vector_get (x, 4);  //gamma2
	double x5 = gsl_vector_get (x, 5); //A2
	double x6 = gsl_vector_get (x, 6);  //mu3
	double x7 = gsl_vector_get (x, 7);  //mu4
	double x8 = gsl_vector_get (x, 8);  //gamma3
	double x9 = gsl_vector_get (x, 9);  //gamma4

	size_t i;

	for ( i = 0; i < n; i++) 	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double s = sigma[i];
		// first dimension
		double xx3 = x3 * x3;
		double xx4 = x4 * x4;
		double b1 =  t1[i] - x1 ;
		double a1 = b1 * b1  +  xx3;
		double b2 =  t2[i] - x2 ;
		double a2 = b2 * b2  +  xx4;
		double c1 = x0 * xx3 * xx4;
		gsl_matrix_set (J, i, 0, xx3 * xx4  / a1 / a2 / s);
		gsl_matrix_set (J, i, 1, 2.0 * b1 * c1  / a2 / a1 / a1 / s );
		gsl_matrix_set (J, i, 2, 2.0 * b2 * c1  / a1 / a2 / a2 / s);
		gsl_matrix_set (J, i, 3, 2.0 * x0 * x3 * xx4 * b1 * b1  / a1 / a1 / a2 / s);
		gsl_matrix_set (J, i, 4, 2.0 * x0 * x4 * xx3 * b2 * b2  / a1 / a2 / a2 / s);
		// second dimension
		double xx8 = x8 * x8;
		double xx9 = x9 * x9;
		double b3 =  t1[i] - x6 ;
		double a3 = b3 * b3  +  xx8;
		double b4 =  t2[i] - x7 ;
		double a4 = b4 * b4  +  xx9;
		double c2 = x5 * xx8 * xx9;

		gsl_matrix_set (J, i, 5, xx8 * xx9  / a3 / a4 / s);
		gsl_matrix_set (J, i, 6, 2.0 * b3 * c2  / a4 / a3 / a3 / s );
		gsl_matrix_set (J, i, 7, 2.0 * b4 * c2  / a3 / a4 / a4 / s);
		gsl_matrix_set (J, i, 8, 2.0 * x5 * x8 * xx9 * b3 * b3  / a3 / a3 / a4 / s);
		gsl_matrix_set (J, i, 9, 2.0 * x5 * x9 * xx8 * b4 * b4  / a3 / a4 / a4 / s);

	}
	return GSL_SUCCESS;
}

inline int twoDimLorentz2_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	twoDimLorentz2_f (x, data, f);
	twoDimLorentz2_df (x, data, J);

	return GSL_SUCCESS;
}
//
/**
 * Which kind of data we want to pass?
 * a point3D data by extracting all the points above a threshold, adding the contour at the threshold
 */
inline int twoDimLorentzFit2 (const size_t p, const vector<point3D<double, double, double>*>& onePkVec,
		double *x_init) {

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	const size_t n = onePkVec.size();
	unsigned int i, iter = 0;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double t1[n], t2[n], y[n], sigma[n];

	for ( i = 0; i < n; i++){
		t1[i] = onePkVec.at( i )->getX();
		t2[i] = onePkVec.at( i )->getY();
		y[i] =  onePkVec.at( i )->getZ();
		sigma[i] = 1.0;
		printf ("data: %8g%8g %8g %g\n", t1[i], t2[i], y[i], sigma[i]);
	}
	struct data2d d = { n, t1, t2, y, sigma};

	gsl_multifit_function_fdf f;
	//	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

	f.f = &twoDimLorentz2_f;
	f.df = & twoDimLorentz2_df;
	f.fdf = & twoDimLorentz2_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	print_state (iter, s);

	do	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		printf ("status = %s\n", gsl_strerror (status));
		print_state (iter, s);

		if (status)
			break;

		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-7, 1e-7);
	} while (status == GSL_CONTINUE && iter < 5000);

	gsl_multifit_covar ( s->J, 0.0, covar );

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

		printf ("Intensity1      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("mu1 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("mu2 = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
		printf ("gamma1     = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
		printf ("gamma2      = %.5f +/- %.5f\n", FIT(4), c*ERR(4));

		printf ("Intensity2      = %.5f +/- %.5f\n", FIT(5), c*ERR(5));
		printf ("mu3 = %.5f +/- %.5f\n", FIT(6), c*ERR(6));
		printf ("mu4 = %.5f +/- %.5f\n", FIT(7), c*ERR(7));
		printf ("gamma3     = %.5f +/- %.5f\n", FIT(8), c*ERR(8));
		printf ("gamma4      = %.5f +/- %.5f\n", FIT(9), c*ERR(9));
	}

	printf ("status = %s\n", gsl_strerror (status));
	for ( unsigned int i = 0; i<p ; i++)
		x_init[i] = gsl_vector_get (s->x, i);
//	x_init[0] = gsl_vector_get (s->x, 0);
//	x_init[1] = gsl_vector_get (s->x, 1);
//	x_init[2] = gsl_vector_get (s->x, 2);
//	x_init[3] = gsl_vector_get (s->x, 3);
//	x_init[4] = gsl_vector_get (s->x, 4);
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return 0;
}



















