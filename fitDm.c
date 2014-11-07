// algorithm to fit phase shift and DM  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "ptimeT.h"

/* Paraboloid centered on (p[0],p[1]), with
 * scale factors (p[2],p[3]) and minimum p[4] */
double my_f (const gsl_vector *v, void *params)
{
	double x, y;
	double *p = (double *)params;
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);
	return p[2] * (x - p[0]) * (x - p[0]) +
		p[3] * (y - p[1]) * (y - p[1]) + p[4];
}

double chiSquare (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;
	//freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	double phaseNchn;

	double s, c00, nu0;
	double c001, c002, c003;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	freqRef = sqrt(c00/(2.0*nu0));

	chi2 = 0.0;
	double P, PS, S;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

int miniseNelderMead (params *param, double ini_guess, double *phase, double *dmFit)
{
	double psrFreq = ((params *)param)->psrFreq;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;

	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	double guess, dmGuess;
	//guess = 0.88;
	guess = ini_guess;
	dmGuess = param->dm;
	//dmGuess = 0.005;
	printf ("phase guess: %.10lf(%.10lf); DM guess: %.5lf\n", ((guess/3.1415926)/(psrFreq*2.0))*1.0e+6, guess, dmGuess);	

	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, guess);
	gsl_vector_set (x, 1, dmGuess);

	// test function
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 0.0001);

	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = chiSquare;
	minex_func.params = param;

	// test function
	//minex_func.f = my_f;
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	//minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 2);
	printf ("Fit for phase shift and DM.\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-6);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		printf ("%5d %10.3e %10.3e f() = %7.3f size = %.6f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval, size);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = gsl_vector_get (s->x, 1);
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

int covariance (void *param, double phase, double dm)
{
	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;

	freqRef = ((params *)param)->freqRef;

	int i,j;

	double s;
	double c00,c11,c01,nu0;
	double c001,c002,c003;
	double c111,c112;
	double c121;
	double phaseNchn;

	c00 = 0.0;
	nu0 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}
	freqRef = sqrt(c00/(2.0*nu0));

	c00 = 0.0;
	c11 = 0.0;
	c01 = 0.0;
	nu0 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		c111 = 0.0;
		c112 = 0.0;
		c121 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c111 += ((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c112 += pow(((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef))),2.0)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c121 += ((j+1)*(j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		c11 += 2.0*((-c111*c111+c001*c112)/s)/(rms[i]*rms[i]);
		c01 += 2.0*((-c111*c002+c001*c121)/s)/(rms[i]*rms[i]);
		//nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	printf ("phase error: %lf; DM error: %lf\n",  ((sqrt(2.0/fabs(c00))/3.1415926)/(psrFreq*2.0))*1.0e+6, sqrt(2.0/fabs(c11)));
	printf ("c00: %lf; c11: %lf; c12: %lf\n", c00, c11, c01);
	printf ("freqRef: %lf\n", freqRef);
	
	return 0;
}

double chiSquareTest (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = 15.9898;

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	int psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	chi2 = 0.0;
	double P, PS, S;
	double phaseNchn;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//phaseNchn = phase;
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(nfreq[nchn/2]*nfreq[nchn/2]));
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

int miniseNelderMeadTest (params *param, double ini_guess, double *phase, double *dmFit)
{
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;

	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	double guess;
	guess = 0.490874;
	//guess = ini_guess;
	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, guess);
	//printf ("initial guess: %lf\n", ((guess/3.1415926)/(param->psrFreq*2.0))*1.0e+6);
	//gsl_vector_set (x, 1, param->dm);

	// test function
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	
	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, 0.01);

	/* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = chiSquareTest;
	minex_func.params = param;

	// test function
	//minex_func.f = my_f;
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	//minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 1);
	printf ("Fit for phase shift and DM (Test).\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-6);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		printf ("%5d %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (s->x, 0), s->fval, size);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = 15.9898;
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

double chiSquare2 (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;
	freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	double phaseNchn;

	chi2 = 0.0;
	double P, PS, S;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

void dfChiSquare2 (const gsl_vector *x, void *param, gsl_vector *df)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;

	freqRef = ((params *)param)->freqRef;

	int i,j;

	double s;
	double c0,c1;
	double c01,c02;
	double c11;
	double phaseNchn;

	c0 = 0.0;
	c1 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c01 = 0.0;
		c02 = 0.0;
		c11 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c01 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c02 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c11 += ((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		}
		c0 += 2.0*((c01*c02)/s)/(rms[i]*rms[i]);
		c1 += 2.0*((c01*c11)/s)/(rms[i]*rms[i]);
		//nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	//printf ("phase error: %lf; DM error: %lf\n",  ((sqrt(2.0/fabs(c00))/3.1415926)/(psrFreq*2.0))*1.0e+6, sqrt(2.0/fabs(c11)));
	//printf ("c00: %lf; c11: %lf; c12: %lf\n", c00, c11, c01);
	//printf ("freqRef: %lf\n", freqRef);
	
	gsl_vector_set(df, 0, c0);
	gsl_vector_set(df, 1, c1);
}

void fdfChiSquare2 (const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = chiSquare2(x, params);
	dfChiSquare2(x, params, df);
}

int miniseD (params *param, double ini_guess, double *phase, double *dmFit)
{
	double psrFreq = ((params *)param)->psrFreq;

	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;

	gsl_multimin_function_fdf my_func;

	/* Initialize method and iterate */
	my_func.n = 2;
	my_func.f = chiSquare2;
	my_func.df = dfChiSquare2;
	my_func.fdf = fdfChiSquare2;
	my_func.params = param;

	/* Starting point */
	double guess, dmGuess;
	//guess = 0.497010;
	guess = ini_guess;
	dmGuess = param->dm;
	//dmGuess = 0.0;
	printf ("phase guess: %.10lf(%.10lf); DM guess: %.5lf\n", ((guess/3.1415926)/(psrFreq*2.0))*1.0e+6, guess, dmGuess);	

	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, guess);
	gsl_vector_set (x, 1, dmGuess);

	
	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc (T, 2);
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-6);

	printf ("Fit for phase shift and DM.\n");

	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		if (status)
			break;
		status = gsl_multimin_test_gradient (s->gradient, 1e-3);
		if (status == GSL_SUCCESS)
			printf ("Minimum found at:\n");
		printf ("%5d %.5f %.5f %10.5f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->f);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = gsl_vector_get (s->x, 1);
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free(x);

	return status;
}

