#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"
#include "simulatePseudoBB.h"

typedef struct cmplx {
  float real;
  float imag;
} cmplx;

cmplx createCmplxAmpPhs(double amp,double phs);
double cmod2(cmplx val);
cmplx cconj(cmplx val);
double calcCrossTerm4(cmplx a,cmplx b,cmplx c,cmplx d);
double calcCrossTerm3(cmplx a,cmplx b,cmplx c);
cmplx cmult(cmplx a,cmplx b);
cmplx cadd(cmplx a,cmplx b);
cmplx cmult3(cmplx a,cmplx b,cmplx c);
cmplx cmult4(cmplx a,cmplx b,cmplx c,cmplx d);
cmplx add9(cmplx a,cmplx b,cmplx c,cmplx d,cmplx e,cmplx f,cmplx g,cmplx h,cmplx i);

void simulatePseudoBB(int nbin,float *ret_aa,float *ret_bb,float *ret_cr,float *ret_ci)
{
  double bw=256; // MHz;
  long npts = 10000;
  long i,j;
  long seed = TKsetSeed();
  double magEx = 40;
  //  double magEy = 5;
  double magSx = 60;
  double magSy = 60;
  double magNx = 45;
  double magNy = 5;
  cmplx AA[npts],BB[npts];
  double AAav,BBav;
  double rABav,iABav;
  cmplx A,B;

  cmplx t1,t2,t3,t4;
  cmplx tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9;

  //
  // Feed parameters
  //
  cmplx ax;
  cmplx ay;
  cmplx c;
  cmplx d;
  cmplx AsB;

  cmplx *Ex,*Ey,*Sx,*Sy,*Nx,*Ny;

  //  ax = createCmplxAmpPhs(1.01,0.09817);
  // ay = createCmplxAmpPhs(1,0.0);
  //c  = createCmplxAmpPhs(0.05432,1.0981);
  //d  = createCmplxAmpPhs(0.0143,1.0981);

  ax = createCmplxAmpPhs(0.5,0.2);
  ay = createCmplxAmpPhs(0.6,0.2);
  c  = createCmplxAmpPhs(0.02,0.652);
  d  = createCmplxAmpPhs(0.03,0.250);
  
  Ex = (cmplx *)malloc(sizeof(cmplx)*npts);
  Ey = (cmplx *)malloc(sizeof(cmplx)*npts);
  Sx = (cmplx *)malloc(sizeof(cmplx)*npts);
  Sy = (cmplx *)malloc(sizeof(cmplx)*npts);
  Nx = (cmplx *)malloc(sizeof(cmplx)*npts);
  Ny = (cmplx *)malloc(sizeof(cmplx)*npts);

  for (j=0;j<nbin;j++)
    {
      AAav = 0;
      BBav = 0;
      rABav = 0;
      iABav = 0;

      for (i=0;i<npts;i++)
	{
	  // Ex and Ey should be set from the known Stokes
	  if (j > 50 && j <200)
	    {
	      //	      Ex[i] = createCmplxAmpPhs(magEx*TKgaussDev(&seed),2*M_PI*TKranDev(&seed));
	      //	      Ey[i] = createCmplxAmpPhs(magEy*TKgaussDev(&seed),2*M_PI*TKranDev(&seed));
	      Ex[i].real =  sqrt(magEx)*TKgaussDev(&seed);
	      Ex[i].imag =  sqrt(magEx)*TKgaussDev(&seed);
	      Ey[i].real =  Ex[i].real;
	      Ey[i].imag =  Ex[i].imag;

	    }
	  else
	    {
	      Ex[i] = createCmplxAmpPhs(0,0);
	      Ey[i] = createCmplxAmpPhs(0,0);
	    }
	  // Assume random and uncorrelated sky noise
	  Sx[i].real =  sqrt(magSx)*TKgaussDev(&seed);
	  Sx[i].imag =  sqrt(magSx)*TKgaussDev(&seed);
	  Sy[i].real =  sqrt(magSy)*TKgaussDev(&seed);
	  Sy[i].imag =  sqrt(magSy)*TKgaussDev(&seed);

	  //	  Sy[i] = createCmplxAmpPhs(magSy*TKgaussDev(&seed),2*M_PI*TKranDev(&seed));
	  
	  // Assume independent noise from the LNAs
	  //	  Nx[i] = createCmplxAmpPhs(magNx*TKgaussDev(&seed),2*M_PI*TKranDev(&seed));
	  //	  Ny[i] = createCmplxAmpPhs(magNy*TKgaussDev(&seed),2*M_PI*TKranDev(&seed));
	  Nx[i].real =  sqrt(magNx)*TKgaussDev(&seed);
	  Nx[i].imag =  sqrt(magNx)*TKgaussDev(&seed);
	  Ny[i].real =  sqrt(magNy)*TKgaussDev(&seed);
	  Ny[i].imag =  sqrt(magNy)*TKgaussDev(&seed);

	  
	  t1 = cadd(Ex[i],Sx[i]);
	  t2 = cadd(Ey[i],Sy[i]);
	  tt1 = cmult(ax,t1);
	  tt2 = cmult(c,t2);
	  t3 = cadd(tt1,tt2);
	  A = cadd(t3,Nx[i]);

	  t1 = cadd(Ey[i],Sy[i]);
	  t2 = cadd(Ex[i],Sx[i]);
	  tt1 = cmult(ay,t1);
	  tt2 = cmult(d,t2);
	  t3 = cadd(tt1,tt2);
	  B = cadd(t3,Ny[i]);

	  AsB   = cmult(cconj(A),B);

	  AAav  += cmod2(A);
	  BBav  += cmod2(B);
	  rABav += AsB.real;
	  iABav += AsB.imag;

	}
      ret_aa[j] = AAav/(double)npts;
      ret_bb[j] = BBav/(double)npts;
      ret_cr[j] = rABav/(double)npts;
      ret_ci[j] = iABav/(double)npts;
      printf("Returning %g %g %g %g\n",ret_aa[j],ret_bb[j],ret_cr[j],ret_ci[j]);
    }
  free(Ex); free(Ey);
  free(Sx); free(Sy);
  free(Nx); free(Ny);

  //  free(aa); free(bb); free(cr); free(ci);
  
  
  
}


cmplx add9(cmplx a,cmplx b,cmplx c,cmplx d,cmplx e,cmplx f,cmplx g,cmplx h,cmplx i)
{
  cmplx ret;
  ret.real = a.real+b.real+c.real+d.real+e.real+f.real+g.real+h.real+i.real;
  ret.imag = a.imag+b.imag+c.imag+d.imag+e.imag+f.imag+g.imag+h.imag+i.imag;
  return ret;
}
cmplx cmult4(cmplx a,cmplx b,cmplx c,cmplx d)
{
  cmplx ret;
  ret = cmult(a,b);
  ret = cmult(ret,c);
  ret = cmult(ret,d);
  return ret;
}

cmplx cmult3(cmplx a,cmplx b,cmplx c)
{
  cmplx ret;
  ret = cmult(a,b);
  ret = cmult(ret,c);
  return ret;
}

cmplx cadd(cmplx a,cmplx b)
{
  cmplx c;
  c.real = a.real+b.real;
  c.imag = a.imag+b.imag;
  return c;
}

// Calculate 2*Re(a*b*c*d)
double calcCrossTerm4(cmplx a,cmplx b,cmplx c,cmplx d)
{
  cmplx tt,tt2,tt3,tt4;
  double ret;

  tt = cmult(a,b);
  tt2 = cmult(tt,c);
  tt3 = cmult(tt2,d);
  ret = 2*tt3.real;
  return ret;
}

// Calculate 2*Re(a*b*c)
double calcCrossTerm3(cmplx a,cmplx b,cmplx c)
{
  cmplx tt,tt2,tt3,tt4;
  double ret;

  tt = cmult(a,b);
  tt2 = cmult(tt,c);
  ret = 2*tt2.real;
  return ret;
}


cmplx cmult(cmplx a,cmplx b)
{
  cmplx c;

  c.real = a.real*b.real-a.imag*b.imag;
  c.imag = a.imag*b.real+b.imag*a.real;
  return c;
}

cmplx cconj(cmplx val)
{
  cmplx ret;
  ret.real = val.real;
  ret.imag = -val.imag;
  return ret;
    
}
double cmod2(cmplx val)
{
  double ret;
  ret = pow(val.real,2)+pow(val.imag,2);
  return ret;
}


cmplx createCmplxAmpPhs(double amp,double phs)
{
  cmplx ret;
  ret.real = amp*cos(phs);
  ret.imag = amp*sin(phs);
  return ret;
}
