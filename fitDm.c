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
	int psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;

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
		phaseNchn = phase + (K*dm)/(nfreq[i]*nfreq[i]*psrFreq);
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += p_s[i][j]*p_s[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

int miniseNelderMead (params *param, double guess, double *phase, double *dmFit)
{
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;

	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, guess);
	gsl_vector_set (x, 1, param->dm);

	// test function
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 1.0);

	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = chiSquare;
	minex_func.params = param;

	// test function
	//minex_func.f = my_f;
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	//minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 2);
	printf ("here\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		//printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval, size);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = gsl_vector_get (s->x, 1);
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

