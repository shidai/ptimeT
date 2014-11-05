#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"

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

int main()
{
  double bw=256; // MHz;
  long npts = 10000;
  long i,j;
  long seed = TKsetSeed();
  double magEx = 5;
  double magEy = 0;
  double magSx = 1;
  double magSy = 1;
  double magNx = 1;
  double magNy = 1;
  cmplx AA[npts],BB[npts];
  double AAav,BBav;
  double rABav,iABav;
  double mean1,mean2,mean3,mean4;
  double mean1_on,mean2_on,mean3_on,mean4_on;
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
  int nc=0,nc_on=0;

  //  ax = createCmplxAmpPhs(0.1,0.2);
  //ay = createCmplxAmpPhs(0.2,0.2);
  //c  = createCmplxAmpPhs(0.05,0.123);
  //d  = createCmplxAmpPhs(0.05,0.123);

  ax = createCmplxAmpPhs(1000,0.0);
  ay = createCmplxAmpPhs(1000,0.0);
  c  = createCmplxAmpPhs(0.0,0.0);
  d  = createCmplxAmpPhs(0.0,0.0);
  
  Ex = (cmplx *)malloc(sizeof(cmplx)*npts);
  Ey = (cmplx *)malloc(sizeof(cmplx)*npts);
  Sx = (cmplx *)malloc(sizeof(cmplx)*npts);
  Sy = (cmplx *)malloc(sizeof(cmplx)*npts);
  Nx = (cmplx *)malloc(sizeof(cmplx)*npts);
  Ny = (cmplx *)malloc(sizeof(cmplx)*npts);

  mean1=mean2=mean3=mean4=0.0;
  mean1_on=mean2_on=mean3_on=mean4_on=0.0;
  for (j=0;j<256;j++)
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

	  // AA
	  // Auto terms 
	  /*	  AA[i] = cmod2(ax)*cmod2(Ex[i]) + 
	    cmod2(ax)*cmod2(Sx[i]) +
	    cmod2(c)*cmod2(Ey[i]) +
	    cmod2(c)*cmod2(Sy[i]) +
	    cmod2(Nx[i]);
	  // Cross terms
	  AA[i] += calcCrossTerm4(ax,cconj(c),Ex[i],cconj(Ey[i]))
	    + calcCrossTerm4(ax,cconj(c),Sx[i],cconj(Sy[i]))
	    + calcCrossTerm4(ax,cconj(c),Ex[i],cconj(Sy[i]))
	    + calcCrossTerm4(ax,cconj(c),Sx[i],cconj(Ey[i]))
	    + calcCrossTerm4(c,cconj(c),Ey[i],cconj(Sy[i]))
	    + calcCrossTerm4(ax,cconj(ax),Ex[i],cconj(Sx[i]))
	    + calcCrossTerm3(ax,Ex[i],cconj(Nx[i]))
	    + calcCrossTerm3(ax,Sx[i],cconj(Nx[i]))
	    + calcCrossTerm3(c,Ey[i],cconj(Nx[i]))
	    + calcCrossTerm3(c,Sy[i],cconj(Nx[i]));

	  // BB
	  // Auto terms 
	  BB[i] = cmod2(ay)*cmod2(Ey[i]) + 
	    cmod2(ay)*cmod2(Sy[i]) +
	    cmod2(d)*cmod2(Ex[i]) +
	    cmod2(d)*cmod2(Sx[i]) +
	    cmod2(Ny[i]);
	  // Cross terms
	  BB[i] += calcCrossTerm4(ay,cconj(d),Ey[i],cconj(Ex[i]))
	    + calcCrossTerm4(ay,cconj(d),Sy[i],cconj(Sx[i]))
	    + calcCrossTerm4(ay,cconj(d),Ey[i],cconj(Sx[i]))
	    + calcCrossTerm4(ay,cconj(d),Sy[i],cconj(Ex[i]))
	    + calcCrossTerm4(d,cconj(d),Ex[i],cconj(Sx[i]))
	    + calcCrossTerm4(ay,cconj(ay),Ey[i],cconj(Sy[i]))
	    + calcCrossTerm3(ay,Ey[i],cconj(Ny[i]))
	    + calcCrossTerm3(ay,Sy[i],cconj(Ny[i]))
	    + calcCrossTerm3(d,Ex[i],cconj(Ny[i]))
	    + calcCrossTerm3(d,Sx[i],cconj(Ny[i]));
	  
	  AAav+= AA[i];
	  BBav+= BB[i];

	  t1 = cadd(cconj(Ex[i]),cconj(Sx[i]));
	  t2 = cadd(Ey[i],Sy[i]);
	  t3 = cadd(Ex[i],Sx[i]);
	  t4 = cadd(cconj(Ey[i]),cconj(Sy[i]));

	  tt1 = cmult4(cconj(ax),ay,t1,t2);
	  tt2 = cmult4(cconj(ax),d,t1,t3);
	  tt3 = cmult3(cconj(ax),t1,Ny[i]);
	  tt4 = cmult4(cconj(c),ay,t4,t2);
	  tt5 = cmult4(cconj(c),d,t4,t3);
	  tt6 = cmult3(cconj(c),t4,Ny[i]);
	  tt7 = cmult3(cconj(Nx[i]),ay,t2);
	  tt8 = cmult3(cconj(Nx[i]),d,t3);
	  tt9 = cmult(cconj(Nx[i]),Ny[i]);
	  //	  printf("Have: %g %g %g %g\n",tt9.real,t2.real,t3.real,t4.real);
	  //	  printf("Have2: %g %g %g %g\n",tt9.imag,t2.imag,t3.imag,t4.imag);

	  AsB = add9(tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9);
	  //	  printf("ASB %g %g\n",AsB.real,AsB.imag);
	  rABav += AsB.real;
	  iABav += AsB.imag;*/
	  
	  t1 = cadd(Ex[i],Sx[i]);
	  t2 = cadd(Ey[i],Sy[i]);
	  tt1 = cmult(ax,t1);
	  tt2 = cmult(c,t2);
	  t3 = cadd(tt1,tt2);
	  /*	  if (j==100)
	    {
	      printf("Have expect %g %g\n",t1.real*ax.real,t1.imag*ax.real);
	      printf("Have %g %g %g %g\n",t1.real,t2.real,tt1.real,tt1.imag);
	      printf("Have2 %g %g %g %g\n",ax.real,ax.imag,Ex[i].real,Ex[i].imag);
	      exit(1);
	      } */
	  //	  printf("Adding %g %g %g %g\n",t3.real,Nx[i].real,t3.imag,Nx[i].imag);
	  A = cadd(t3,Nx[i]);

	  t1 = cadd(Ey[i],Sy[i]);
	  t2 = cadd(Ex[i],Sx[i]);
	  tt1 = cmult(ay,t1);
	  tt2 = cmult(d,t2);
	  t3 = cadd(tt1,tt2);
	  B = cadd(t3,Ny[i]);
	  //	  printf("%g %g %g %g\n",A.real,A.imag,B.real,B.imag);
	  //	  AA[i] = cmult(A,cconj(A));
	  //	  BB[i] = cmult(B,cconj(B));
	  AsB   = cmult(cconj(A),B);
	  //	  printf("%g %g\n",AsB.real,AsB.imag);
	  AAav  += cmod2(A);
	  BBav  += cmod2(B);
	  rABav += AsB.real;
	  iABav += AsB.imag;

	  //      printf("%g\n",AA[i]);
	}
      if (j <= 50 || j >= 200)
	{
	  mean1 += AAav/npts;
	  mean2 += BBav/npts;
	  mean3 += rABav/npts;
	  mean4 += iABav/npts;	  
	  nc++;
	}
      else 
	{
	  mean1_on += AAav/npts;
	  mean2_on += BBav/npts;
	  mean3_on += rABav/npts;
	  mean4_on += iABav/npts;	  
	  nc_on++;
	}
      //      exit(1);
            printf("%g %g %g %g\n",AAav/(double)npts,BBav/(double)npts,rABav/(double)npts,iABav/(double)npts);
    }
  printf("Means: %g %g %g %g\n",mean1/(double)nc,mean2/(double)nc,mean3/(double)nc,mean4/(double)nc);
  {
    double AA,BB,CR,CI;
    AA = (mean1_on/(double)nc_on-mean1/(double)nc);
    BB = (mean2_on/(double)nc_on-mean2/(double)nc);
    CR = (mean3_on/(double)nc_on-mean3/(double)nc);
    CI = (mean4_on/(double)nc_on-mean4/(double)nc);
    printf("Stokes: %g %g %g %g %g\n",AA+BB,AA-BB,2*CR,2*CI,
	   sqrt(pow(AA-BB,2)+pow(2*CR,2)+pow(2*CI,2)));
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
  c.imag = a.imag*b.real+a.real*b.imag;
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
