/*
 *
 * Currently working on:
 * - scintillation
 *
 * Program to simulate data sets obtained for the PPTA project using the DFB systems.
 *
 * Compilation: ./compile
 *
 * 
 * Inputs:
 *  - pulsar timing model for the pulsar to simulate (with exact parameters)
 *  - pulsar timing model used in the simulated observations 
 *  - nbin   number of phase bins
 *  - nchan  number of channels 
 *  - nsub   number of subintegrations (currently only simulate 1)
 *  - list of primary header parameters
 * 
 * Requirements
 *  - requires up-to-date psrheader.fits file in current directory

 * Output data sets:
 *
 * FITS files containing simulated folded calibration signals for each observation
 * FITS files containing simulated folded pulsar data for each observation for a given pulsar
 * Simulated flux calibration set of files using Hydra A
 *
 *
 * Real issues
 *
 * Polarisation calibration (DONE)
 * Faraday rotation
 * DM variations (DONE)
 * Pulse jitter
 * Pulse shape changes
 * Frequency-dependence of pulse shape
 * Scintillation
 * Incorrect timing model
 * RFI
 * Digitization
 * Machine precision issues 
 *
 * --
 * Want to be able to fix the random number seed so the results are reproducible
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "simulatePseudoBB.h"

#define MAX_NBIN 4096
#define MAX_NCHAN 4096
#define MAX_RFI_FILES 10

long seed;

typedef struct dataStruct {
  char fname[128];
  float aa[MAX_NBIN*MAX_NCHAN],bb[MAX_NBIN*MAX_NCHAN];
  float cr[MAX_NBIN*MAX_NCHAN],ci[MAX_NBIN*MAX_NCHAN]; // Received data
  char  primaryHeaderParams[1024]; // Filename of file containing primary header parameters

  // Source parameters
  long double period;     // Period of source in seconds

  // Observing setup parameters
  int   obsMode; // 1 = levcal, 2 = psr, 3 = fluxcal
  double cFreq;  // Centre frequency for observation (MHz)
  char   srcName[128]; // Name of source
  double scanLength; // Length of observation in seconds
  int    stt_imjd;   // Start date integer (MJD)
  double stt_smjd;   // Start time (s)
  double stt_offs;   // Start time offset (s)
  double tsubRequested; // Time for each subintegration requested (s)
  double tsub;       // Time for each subintegration [modified to give an integral number of pulses per subint] (s)
  double pa;         // Parallactic angle (deg)

  // Backend details
  double obsBW;  // Bandwidth in MHz
  int    nchan;  // Number of frequency channels
  int    nbin;   // Number of bins in profile
  int    nsub;   // Number of subintegrations
  int    npol;   // Number of polarisations

  double scaleAA; // Scaling factor for AA
  double scaleBB; // Scaling factor for BB
  double scaleCR; // Scaling factor for real part of AB
  double scaleCI; // Scaling factor for imaginary part of AB

  //  double whiteNoiseA;
  //  double whiteNoiseB;
  double receiverNoise;
  double whiteNoiseI;

  // Cal details
  double cal_freq; // Frequency of cal (Hz)
  double cal_phs;  // Cal phase (wrt start time)
  double cal_dcyc; // Cal duty cycle
  double cal_amp;   // Amplitude of cal 
  // Pulsar parameters
  char   template_I[128]; // Filename of Von Mises template for total intensity
  char   template_II[128]; // Filename of Von Mises template for unpolarised power
  char   compTemplateI[MAX_NCHAN][128];
  int    nCompTemplate;
  char   template_Q[128]; // Filename of Von Mises template for Stokes Q
  char   template_U[128]; // Filename of Von Mises template for Stokes U
  char   template_V[128]; // Filename of Von Mises template for Stokes V
  char   exact_ephemeris[128]; // Ephemeris used for the simulation of arrival times
  double phaseOffset[MAX_NCHAN]; // Phase offset from start of observation to pulse arrival time
  double dm;              // DM value extracted from the parameter file
  double dmset;           // DM value given by user
  // Scintillation parameters
  double scint_freqbw;
  double scint_ts;
  
  // Feed parameters
  char   feedParamFile[1024];
  double absGain[MAX_NCHAN];
  double diffGain[MAX_NCHAN];
  double diffPhase[MAX_NCHAN];
  double cc_eps1[MAX_NCHAN],cc_eps2[MAX_NCHAN];
  double cc_phi1[MAX_NCHAN],cc_phi2[MAX_NCHAN];
  double para;

  // RFI
  int nRFI;
  char  rfiFile[MAX_RFI_FILES][128];
  float rfiVal[MAX_RFI_FILES][MAX_NCHAN];
  float rfiBin[MAX_RFI_FILES][MAX_NCHAN];
  int   rfiNbinNum;
  
} dataStruct;

void addRFI(dataStruct *data);
void createFitsFile(fitsfile *fptr,dataStruct *data);
void writeHeaderParameters(fitsfile *fptr,char *fname,dataStruct *data);
void removeTables(fitsfile *fptr,dataStruct *data);
void writeNewSub(fitsfile *fptr,dataStruct *data,int subint,long double timeFromStart);
void simulateData(dataStruct *data,float **scint,int subint);
void loadData(dataStruct *data);
void simulateCal(dataStruct *data);
void simulatePsr(dataStruct *data,float **scint,int subint);
void convertStokesAABBCRCI(double *stokes_act,double *obsact);
void createM_amp(double gain,double phase,double m_amp[4][4]);
void createM_pa(double pa,double m_pa[4][4]);
void createM_feed(double gamma,double m_feed[4][4]);
void multVectMat(double m[4][4],double *vIn,double *vOut);
void mult4x4(double a[4][4],double b[4][4],double c[4][4]);
int readTemplate(char *fname,double *concentration,double *centre,double *height);
int readTemplateDf(char *fname,double *concentration,double *centre,double *height,
		   double *concentration_df,double *centre_df,double *height_df,double *fiducial,int nbin);
void writeEphemeris(fitsfile *fptr,dataStruct *data);
void writePredictor(fitsfile *fptr,char *fname);
void runTempo2(dataStruct *data);
void simulateFeed(dataStruct *data,double *stokes_act,double *stokes_obs,int fchan);
long double calcLocalSiderealTime(long double t,dataStruct *data);
double calculatePara(dataStruct *data);
void readFeedParameters(dataStruct *data);
void initialiseData(dataStruct *data);
void simulateOffPulseNoise(double *simWhiteAA,double *simWhiteBB,double *simWhiteCR,double *simWhiteCI,dataStruct *data,int nChan,long *seed);
void readData(char *fname,dataStruct *data);
void readScintillationFile(float **scint);
void convertAABBCRCI_stokes(double aa,double bb,double cr,double ci,double *stokes);

/*
int simulatePPTA(void)
{
  dataStruct *data;
  T2Predictor pred;
  int i;
  int havePrimaryHeaderParams=0;
  int status=0;
  fitsfile *fptr;
  FILE *fin;
  char file[128];
  int  subint=1;
  char obsTable[128]="obsTable";
  char signalFile[128]="NULL";
  char line[128];
  char temp[128];
  long double timeFromStart,tsubLength;
  long double phase0,mjd0,freq,f0,toff;

  // Scintillation simulation parameters
  float **scint;
  int nx=16384;
  //int ny=4024;
  int ny=1024;

  // Set random number seed
  seed  = TKsetSeed();

  scint = (float **)malloc(sizeof(float *)*ny);
  for (i=0;i<ny;i++)
    scint[i] = (float *)malloc(sizeof(float)*nx);

  readScintillationFile(scint);


  data = (dataStruct *)malloc(sizeof(dataStruct));

  initialiseData(data);

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-phead")==0)
	{
	  havePrimaryHeaderParams=1;
	  strcpy(data->primaryHeaderParams,argv[++i]);
	}
      else if (strcmp(argv[i],"-signal")==0)
	strcpy(signalFile,argv[++i]);
      else if (strcmp(argv[i],"-obsTable")==0)
	strcpy(obsTable,argv[++i]);
      else
	{
	  printf("Unknown command line argument: %s\n",argv[i]);
	}
    }
  if (!(fin = fopen(obsTable,"r")))
    {
      printf("Unable to open file >%s<\n",obsTable);
      printf("Must provide an observation file using -obsTable\n");
      exit(1);
    }
  while (!feof(fin))
    {
      if (fscanf(fin,"%s",line)==1)
	{
	  if (strcasecmp(line,"SRC")==0)
	    fscanf(fin,"%s",data->srcName);
	  else if (strcasecmp(line,"PHEAD")==0)
	    {
	      havePrimaryHeaderParams=1;
	      fscanf(fin,"%s",data->primaryHeaderParams);
	    }
	  else if (strcasecmp(line,"FILE")==0)
	    fscanf(fin,"%s",data->fname);
	  else if (strcasecmp(line,"TYPE")==0)
	    {
	      fscanf(fin,"%s",temp);
	      if (strcasecmp(temp,"PSR")==0)
		data->obsMode = 2;
	      else if (strcasecmp(temp,"PSR_LOAD")==0)
		data->obsMode = 3;
	      else if (strcasecmp(temp,"CAL")==0)
		data->obsMode = 1;
	      else
		{
		  printf("Error: Unknown source type: %s\n",temp);
		  exit(1);
		}
	    }
	  else if (strcasecmp(line,"DMSET")==0)
	    fscanf(fin,"%lf",&(data->dmset));
	  else if (strcasecmp(line,"SCINT_TS")==0)
	    fscanf(fin,"%lf",&(data->scint_ts));
	  else if (strcasecmp(line,"SCINT_FREQBW")==0)
	    fscanf(fin,"%lf",&(data->scint_freqbw));
	  else if (strcasecmp(line,"STT_IMJD")==0)
	    fscanf(fin,"%d",&(data->stt_imjd));
	  else if (strcasecmp(line,"STT_SMJD")==0)
	    fscanf(fin,"%lf",&(data->stt_smjd));
	  else if (strcasecmp(line,"STT_OFFS")==0)
	    fscanf(fin,"%lf",&(data->stt_offs));
	  else if (strcasecmp(line,"FEED_PARAM")==0)
	    fscanf(fin,"%s",data->feedParamFile);
	  else if (strcasecmp(line,"TSUB")==0)
	    fscanf(fin,"%lf",&(data->tsubRequested));
	  else if (strcasecmp(line,"NSUB")==0)
	    fscanf(fin,"%d",&(data->nsub));
	  else if (strcasecmp(line,"NCHAN")==0)
	    fscanf(fin,"%d",&(data->nchan));
	  else if (strcasecmp(line,"NBIN")==0)
	    fscanf(fin,"%d",&(data->nbin));
	  else if (strcasecmp(line,"EXACT_EPHEMERIS")==0)
	    fscanf(fin,"%s",data->exact_ephemeris);
	  else if (strcasecmp(line,"WHITE_NOISE_I")==0)
	    fscanf(fin,"%lf",&(data->whiteNoiseI));
	  else if (strcasecmp(line,"RECEIVER_NOISE")==0)
	    fscanf(fin,"%lf",&(data->receiverNoise));
	  else if (strcasecmp(line,"TEMPLATE_II")==0)
	    fscanf(fin,"%s",data->template_II);
	  else if (strcasecmp(line,"TEMPLATE_I")==0)
	    fscanf(fin,"%s",data->template_I);
	  else if (strcasecmp(line,"RFI")==0)
	    fscanf(fin,"%s",data->rfiFile[(data->nRFI)++]);
	  else if (strcasecmp(line,"COMP_TEMPLATE_I")==0)
	    {
	      fscanf(fin,"%s",data->compTemplateI[data->nCompTemplate]);
	      (data->nCompTemplate)++;
	    }
	  else if (strcasecmp(line,"TEMPLATE_Q")==0)
	    fscanf(fin,"%s",data->template_Q);
	  else if (strcasecmp(line,"TEMPLATE_U")==0)
	    fscanf(fin,"%s",data->template_U);
	  else if (strcasecmp(line,"TEMPLATE_V")==0)
	    fscanf(fin,"%s",data->template_V);
	  else if (strcasecmp(line,"CAL_FREQ")==0)
	    fscanf(fin,"%lf",&(data->cal_freq));	    
	  else if (strcasecmp(line,"CAL_PHS")==0)
	    fscanf(fin,"%lf",&(data->cal_phs));	    
	  else if (strcasecmp(line,"CAL_DCYC")==0)
	    fscanf(fin,"%lf",&(data->cal_dcyc));	    
	  else if (strcasecmp(line,"CAL_AMP")==0)
	    fscanf(fin,"%lf",&(data->cal_amp));	    
	  else if (strcasecmp(line,"PERIOD")==0)
	    fscanf(fin,"%Lf",&(data->period));	    
	  else if (strcasecmp(line,"END_OBS")==0)
	    {
	      if (havePrimaryHeaderParams==0)
		{
		  printf("Error: must use PHEAD to provide the filename of a file containing the primary header parameters\n");
		  exit(1);
		}

	      if (data->nCompTemplate > 0 && (data->nCompTemplate!=data->nchan))
		{
		  printf("WARNING: nchan != nCompTemplate --- changing nchan to %d\n",data->nCompTemplate);
		  data->nchan  = data->nCompTemplate;
		}
	      
	      sprintf(file,"!%s(psrheader.fits)",data->fname);
	      fits_create_file(&fptr,file,&status);
	      fits_report_error(stdout,status);
	      
	      // Write all the primary header information
	      createFitsFile(fptr,data);
	      
	      // For a pulsar observation: write the ephemeris into the file
	      // and extract the DM
	      if (data->obsMode == 2 || data->obsMode == 3)
		writeEphemeris(fptr,data);
	      // Read feed parameteres
	      readFeedParameters(data);

	      // For a pulsar observation: write the predictor file
	      if (data->obsMode == 2 || data->obsMode == 3)
		{
		  int ret;


		  T2Predictor_Init(&pred);
		  runTempo2(data);
		  writePredictor(fptr,"t2pred.dat");
		  // Calculate phase offset for the pulse arrival time
		  if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
		    {
		      printf("Error: unable to read predictor\n");
		      exit(1);
		    }


		}
	      mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + data->tsubRequested*0.5/86400.0L;
	      data->period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
	      data->tsub = ((int)(data->tsubRequested/data->period+0.5))*data->period;
	      printf("At this point %.10f\n",(double)data->tsub);
	      // Iterate once more	      
	      mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + data->tsub*0.5/86400.0L;
	      data->period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
	      data->tsub = ((int)(data->tsubRequested/data->period+0.5))*data->period;
	      printf("At this point2 %.10f\n",(double)data->tsub);

	      timeFromStart = 0;
		// GOT HERE
	      for (subint=1;subint<=data->nsub;subint++)
		{
		  if (data->obsMode == 2 || data->obsMode == 3)
		    {
		      toff = (int)(data->tsub/2.0/data->period+0.5)*data->period;
		      mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + (timeFromStart + toff)/86400.0L; // Centre of subint
		      data->period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
		      data->tsub = ((int)(data->tsubRequested/data->period+0.5))*data->period;
		      toff = (int)(data->tsub/2.0/data->period+0.5)*data->period;
		      printf("At this point3 %.10f %.15Lf\n",(double)data->tsub,data->period);
		      mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + (timeFromStart + toff)/86400.0L; // Centre of subint
		      // Iterate again
		      data->period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
		      data->tsub = ((int)(data->tsubRequested/data->period+0.5))*data->period;
		      printf("At this point4 %.10f\n",(double)data->tsub);
		      toff = (int)(data->tsub/2.0/data->period+0.5)*data->period;
		      mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + (timeFromStart + toff)/86400.0L; // Centre of subint

		      for (i=0;i<data->nchan;i++)
			{
			  f0 = data->cFreq + fabs(data->obsBW)/2.0; // Highest frequency
			  freq = f0-fabs(data->obsBW/(double)data->nchan)*i;
			  phase0 = T2Predictor_GetPhase(&pred,mjd0,freq);
			  data->phaseOffset[i] = (phase0-floor(phase0));
			}
		      //		  printf("Phase %.16f %.16Lf\n",data->phaseOffset,phase0);
		      data->period      = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
		      printf("PHASE OFFSET: %.15Lf %.15Lf %.15Lf %g\n",mjd0,phase0,data->period,data->phaseOffset[data->nchan-1]);
		    }

		  // New subintegration
		  if (strcmp(signalFile,"NULL")==0)
		    simulateData(data,scint,subint);
		  else
		    readData(signalFile,data);
		  writeNewSub(fptr,data,subint,timeFromStart+toff);
		  timeFromStart += data->tsub;

		}
	      printf("Closing file\n");
	      fits_close_file(fptr,&status);
	      initialiseData(data);
	      if (data->obsMode == 2 || data->obsMode == 3)
		T2Predictor_Destroy(&pred);
	    }
	}
    }
  fclose(fin);
  printf("Goodbye\n");
  
  for (i=0;i<ny;i++)
    free(scint[i]);
  free(scint);

  free(data);
}
*/

void simulateData(dataStruct *data,float **scint,int subint)
{
  if (data->obsMode==1)
    simulateCal(data);
  else if (data->obsMode==2)
    simulatePsr(data,scint,subint);
  else if (data->obsMode==3)
    loadData(data);

  // Add in RFI
  addRFI(data);
}

void addRFI(dataStruct *data)
{
  //  int nRFI;
  //  char  rfiFile[MAX_RFI_FILES][128];
  //  float rfiVal[MAX_RFI_FILES][MAX_NCHAN];
  //  float rfiBin[MAX_RFI_FILES][MAX_NCHAN];
  //  int   rfiNbinNum;
  float failureTimes[5000];
  float rfiAmp[5000];
  int nf=0;
  int i,j,k;
  FILE *fin;
  int nPhase;
  double tf=0;

  //  long seed = TKsetSeed();
  int ibin;

  for (i=0;i<data->nRFI;i++)
    {
      if (!(fin = fopen(data->rfiFile[i],"r")))
	{
	  printf("Unable to open RFI file: %s\n",data->rfiFile[i]);
	  exit(1);
	}
      for (k=0;k<data->nchan;k++)
	fscanf(fin,"%f",&(data->rfiVal[i][k])); // Only one bin
    }

  // Process this subintegration
  nPhase = data->tsub/data->period;
  for (i=0;i<data->nRFI;i++)
    {
      do {
	tf += 0.3+TKgaussDev(&seed)*0.01;
	failureTimes[nf] = tf;
	rfiAmp[nf] = 1+TKgaussDev(&seed)*0.1;
	nf++;
      } while (tf < data->tsub);
      printf("Simulating %d RFI spikes\n",nf);

      for (j=0;j<nf;j++)
	{
	  ibin = ((failureTimes[j]/data->period)-(int)(failureTimes[j]/data->period))*data->nbin;
	  printf("ibun = %d\n",ibin);
	  for (k=0;k<data->nchan;k++)
	    {
	      data->aa[k*data->nbin+ibin] = sqrt(pow(data->aa[k*data->nbin+ibin],2)+pow(rfiAmp[j]*data->rfiVal[i][k]/2.0,2)); // SHOULD DO THIS PROPERLY
	      data->bb[k*data->nbin+ibin] = sqrt(pow(data->bb[k*data->nbin+ibin],2)+pow(rfiAmp[j]*data->rfiVal[i][k]/2.0,2)); // SHOULD DO THIS PROPERLY
	      //	  data->cr[k*data->nbin+10] += data->rfiVal[i][k];
	      //	  data->ci[k*data->nbin+10] += data->rfiVal[i][k];
	    }
	}
    }

}

void initialiseData(dataStruct *data)
{
  int i;
  data->obsMode = -1;
  strcpy(data->fname,"NOT SET");
  strcpy(data->primaryHeaderParams,"NOT SET");
  strcpy(data->srcName,"NOT SET");
  data->period = 0.0L;
  data->nRFI = 0;
  data->scanLength=0.0;
  data->stt_imjd = 0;
  data->stt_smjd = 0;
  data->stt_offs = 0;
  data->tsub     = 0;
  data->pa       = 0;
  data->obsBW    = 0;
  data->nchan    = 0;
  data->nbin     = 0;
  data->nsub     = 0;
  data->npol     = 0;
  data->scaleAA  = 0.0;
  data->scaleBB  = 0.0;
  data->scaleCR  = 0.0;
  data->scaleCI  = 0.0;
  data->receiverNoise = 0.0;
  data->whiteNoiseI = 0.0;
  data->cal_freq = 0.0;
  data->cal_phs  = 0.0;
  data->cal_dcyc = 0.0;
  data->cal_amp = 1.0;
  data->scint_ts  = 0.0;
  data->scint_freqbw = 0.0;
  strcpy(data->template_I,"NOT SET");
  strcpy(data->template_II,"NOT SET");
  strcpy(data->template_Q,"NOT SET");
  strcpy(data->template_U,"NOT SET");
  strcpy(data->template_V,"NOT SET");
  strcpy(data->exact_ephemeris,"NOT SET");
  for (i=0;i<MAX_NCHAN;i++)
    data->phaseOffset[i] = 0.0;
  data->dm = 0.0;
  data->dmset = 0.0;
  strcpy(data->feedParamFile,"NOT SET");
  data->para=0.0;
  data->nCompTemplate = 0;

  // Now set some defaults -- SHOULD FIX
  data->cFreq = 1369.0;
  data->obsBW = -256;
    //  data->cFreq = 1700;
    //  data->obsBW = -2500;
  //  data->cFreq = 2550;
  //  data->obsBW = -4096;
  data->scanLength = 60;
  data->npol     = 4;
  data->scaleAA  = 1;
  data->scaleBB  = 1;
  data->scaleCR  = 1;
  data->scaleCI  = 1;
  data->whiteNoiseI = 1;
}

void loadData(dataStruct *data)
{
  int i,k,j;
  FILE *fin,*finFile;
  double f0,fc;
  int fileOpen;
  double f1[4096],f2[4096];
  char fname[4096][128];
  int nFiles=0;
  float dummy;

  finFile = fopen("1022_wideband.dat","r");
  while (!feof(finFile))
    {
      if (fscanf(finFile,"%lf %lf %s",&f1[nFiles],&f2[nFiles],fname[nFiles])==3)
	nFiles++;
    }
  fclose(finFile);
  printf("nFiles = %d\n",nFiles);
  f0 = data->cFreq + fabs(data->obsBW)/2.0; // Highest frequency
  for (k=0;k<data->nchan;k++)  
    {
      fc = f0 - fabs(data->obsBW/(double)data->nchan)*k;
      fileOpen=0;
      for (j=0;j<nFiles;j++)
	{
	  if (fc >= f1[j] && fc < f2[j])
	    {
	      fin = fopen(fname[j],"r");
	      fileOpen=1;
	    }
	}
      for (i=0;i<data->nbin;i++)
	{
	  if (fileOpen==1)
	    {
	      fscanf(fin,"%f %f %f %f",&dummy,&dummy,&dummy,&(data->aa[k*data->nbin+i]));
	    }
	  else
	    data->aa[k*data->nbin+i] = 0;
	}
      if (fileOpen==1)
	fclose(fin);
    }
}

void simulatePsr(dataStruct *data,float **scint,int subint)
{
  int nCompII,nCompQ,nCompU,nCompV,nCompI;
  double concentrationII[128],centreII[128],heightII[128];
  double concentrationI[128],centreI[128],heightI[128];
  double concentrationI_df[128],centreI_df[128],heightI_df[128];
  double height,centre,concentration;
  double fiducial_f[128];
  double concentrationQ[128],centreQ[128],heightQ[128];
  double concentrationU[128],centreU[128],heightU[128];
  double concentrationV[128],centreV[128],heightV[128];
  double stokes_act[4],obsact[4],stokes_obs[4];
  double II,I,Q,U,V;
  int i,j,k;
  long double x0minO = data->phaseOffset[0]*(double)data->nbin,x0min;
  double dm,tdm,bdm;
  int totIn=0;
  double df;
  int newSession =1 ;
  double scintScale=1;
  double fc;
  double boxX,boxY0,dt,ratio,ratio2,rf,sdiff,f2;
  double f0,f1,dstep,scint_ts_f,scint_freqbw_f;
  double tint = data->tsub;
  int nx=16384;
  int k2,l,nc;
  static int xs=0;
  double simFlux,scint_ts_f0,fa,fb;

  if (strcmp(data->template_I,"NOT SET")==0 && data->nCompTemplate == 0)// Have full Stokes
    {
      printf("Reading polarisation components\n");
      nCompII = readTemplate(data->template_II,concentrationII,centreII,heightII);
      nCompQ = readTemplate(data->template_Q,concentrationQ,centreQ,heightQ);
      nCompU = readTemplate(data->template_U,concentrationU,centreU,heightU);
      nCompV = readTemplate(data->template_V,concentrationV,centreV,heightV);
    }
  else
    {
      totIn=1;
      data->npol=1;

      if (data->nCompTemplate == 0)
	{
	  //      nCompI = readTemplate(data->template_I,concentrationI,centreI,heightI);
	  nCompI = readTemplateDf(data->template_I,concentrationI,centreI,heightI,concentrationI_df,centreI_df,heightI_df,fiducial_f,data->nbin);
	}
    }
  if (data->dmset!=0)
    dm = data->dmset;
  else
    dm = data->dm;
  f0 = data->cFreq + fabs(data->obsBW)/2.0; // Highest frequency

  //  printf("DM value = %g\n",dm);
  // Note currently DM smear being put in through predictor

  data->pa = calculatePara(data);
  printf("Simulating pulsar\n");
  for (k=0;k<data->nchan;k++)
    {      
      if (data->nCompTemplate > 0)
	{
	  nCompI = readTemplateDf(data->compTemplateI[k],concentrationI,centreI,heightI,concentrationI_df,centreI_df,heightI_df,fiducial_f,data->nbin);
	}
      x0minO = data->phaseOffset[k]*(double)data->nbin;
      x0min = x0minO;

      // Scale for scintillation
      if (data->scint_ts > 0)
	{
	  fc = f0 - fabs(data->obsBW/(double)data->nchan)*k;
	  // Simulation parameters
	  rf = 20;
	  sdiff = 4.3;
	  
	  // Deal with scintillation
	  scint_ts_f     = data->scint_ts*pow(fc/1440,1.2);
	  scint_freqbw_f = data->scint_freqbw*pow(fc/1440,4.4);
	  ratio = pow(rf/sdiff,2);
	  ratio2 = fc*1e6/data->scint_freqbw;
	  f2 = pow(10,log10(ratio2/ratio)/(-3.4));
	  
	  scint_ts_f0 = data->scint_ts*pow(1440/f2/1440,1.2);
	  dt = scint_ts_f0/sdiff;
	  boxX = tint/dt;
	  //	  printf("boxX = %g\n",boxX);
	  //		  printf("dt %g %g %g %g %g %g %g\n",dt,scint_ts_f0,sdiff,f2,ratio2,ratio,data->scint_freqbw); 
	  // Now find correct section in the simulation for the particular receiver band
	  fa = f2*(fc-fabs(data->obsBW/(double)data->nchan))/1440; // Measurements made at 1.4GHz
	  fb = f2*(fc)/1440;
	  if (fa > 1.6) fa = 1.6; // Limit of simulation
	  if (fa < 0.4) 
	    {
	      fb = 0.4+(fb-fa);
	      fa = 0.4;
	    }
	  
	  if (fb > 1.6) fb = 1.6; // Limit of simulation
	  dstep = (1.6-0.4)/1024.0;
	  simFlux=0.0;
	  nc=0;
	  if (newSession==1) // This is only to ensure that the frequencies are taken from the same time section
	    {
	      if (subint==1) // This checks whether we are processing adjacent subints
		xs = (int)(TKranDev(&seed)*(nx-3*86400.0/dt)); // Must fix
	      xs += (int)((subint-1)*data->tsub/dt+0.5);
	      newSession=0;
	    }
	  //	  printf("xs = %d dxs = %d %g %g\n",xs,(int)((subint-1)*data->tsub/dt+0.5),fa,fb);
	  //	  printf("#Have %d %g %g %g %g %g %d %d %d\n",k,f0,f1,(double)boxX,dt,simFlux,nc,(int)((f0-0.4)/dstep+0.5),(int)((f1-0.4)/dstep+0.5));
	  for (k2=xs;k2<xs+boxX;k2++)
	    {
		   	      printf("Processing %d %d %d\n",l,k2,xs);
	      for (l=(int)((fa-0.4)/dstep);l<(int)((fb-0.4)/dstep+1.0);l++)
		{
		   	      printf("Processing %d %d %d\n",l,k2,xs);
		   	      printf("Processing\n");
		   	      printf("Processing %g\n",scint[l][k2]);
		  simFlux+=scint[l][k2]; // Should think carefully about this averaging
		  		  printf("Done process: %d\n",nc);
		  nc++;
		}
	    }
	  simFlux/=(double)nc;
	  //	  printf("simFlux %d %g %g %g %d %g\n",k,fa,fb,simFlux,nc,dstep);
	  scintScale = simFlux; //*4;
	  //	  printf("#Have %d %g %g %g %g %g %g %d %d %d\n",k,fc,fa,fb,(double)boxX,dt,simFlux,nc,(int)((fa-0.4)/dstep+0.5),(int)((fb-0.4)/dstep+0.5)-1);
	}
      

      for (i=0;i<data->nbin;i++)
	{
	  if (totIn==0)
	    {
	      double chanbw = fabs(data->obsBW/(double)data->nchan);
	      double noise = data->receiverNoise/sqrt(data->tsub*chanbw);
	      data->aa[k*data->nbin+i] = TKgaussDev(&seed)*data->receiverNoise;
	      data->bb[k*data->nbin+i] = TKgaussDev(&seed)*data->receiverNoise;
	      data->cr[k*data->nbin+i] = TKgaussDev(&seed)*data->receiverNoise/sqrt(2); // MUST FIX
	      data->ci[k*data->nbin+i] = TKgaussDev(&seed)*data->receiverNoise/sqrt(2); 
	    
	      II = I = Q = U = V = 0.0;
	      for (j=0;j<nCompII;j++)
		II += heightII[j]*exp(concentrationII[j]*(cos((i+x0min-centreII[j])/(double)data->nbin*2*M_PI)-1));
	      for (j=0;j<nCompQ;j++)
		Q += heightQ[j]*exp(concentrationQ[j]*(cos((i+x0min-centreQ[j])/(double)data->nbin*2*M_PI)-1));
	      for (j=0;j<nCompU;j++)
		U += heightU[j]*exp(concentrationU[j]*(cos((i+x0min-centreU[j])/(double)data->nbin*2*M_PI)-1));
	      for (j=0;j<nCompV;j++)
		V += heightV[j]*exp(concentrationV[j]*(cos((i+x0min-centreV[j])/(double)data->nbin*2*M_PI)-1));
	      I = sqrt(II*II + Q*Q + U*U + V*V);
	      
	      stokes_act[0] = I;
	      stokes_act[1] = Q;
	      stokes_act[2] = U;
	      stokes_act[3] = V;
	      
	      simulateFeed(data,stokes_act,stokes_obs,k);
	      
	      convertStokesAABBCRCI(stokes_obs,obsact);


	      data->aa[k*data->nbin+i] += obsact[0];
	      data->bb[k*data->nbin+i] += obsact[1];	  
	      data->cr[k*data->nbin+i] += obsact[2];
	      data->ci[k*data->nbin+i] += obsact[3];

	      data->aa[k*data->nbin+i]*=data->scaleAA*data->absGain[k];
	      data->bb[k*data->nbin+i]*=data->scaleBB*data->absGain[k];
	      data->cr[k*data->nbin+i]*=data->scaleCR*data->absGain[k];
	      data->ci[k*data->nbin+i]*=data->scaleCI*data->absGain[k];
	      //	      printf("result: %g %g %g %g\n",data->aa[k*data->nbin+i],data->bb[k*data->nbin+i],
	      //		     data->cr[k*data->nbin+i],data->ci[k*data->nbin+i]);
	    }
	  else
	    {
	      I = 0.0;
	      for (j=0;j<nCompI;j++)
		{
		  df = f0-fabs(data->obsBW/(double)data->nchan)*k-fiducial_f[j];
		  height = heightI[j] + heightI_df[j]*df;
		  concentration = concentrationI[j]  + concentrationI_df[j]*df;
		  centre = centreI[j] + centreI_df[j]*df;
		  I += height*exp(concentration*(cos((i+x0min-centre)/(double)data->nbin*2*M_PI)-1));
		}
	      //	      if (k > 100 && k <200)
	      //		data->aa[k*data->nbin+i] = 0;
	      //	      else
		data->aa[k*data->nbin+i] = I*scintScale+TKgaussDev(&seed)*data->whiteNoiseI;
	    }
       }
    }
  printf("Complete simulating pulsar\n");
}

void readData(char *fname,dataStruct *data)
{
  int j;
  printf("In here\n");
  data->period = 1.0/11.123L;

  for (j=0;j<data->nchan;j++)
    {
      simulatePseudoBB(data->nbin,&(data->aa[j*data->nbin]),
		      &(data->bb[j*data->nbin]),
		      &(data->cr[j*data->nbin]),
		      &(data->ci[j*data->nbin]));
    }
  for (j=0;j<data->nbin;j++)
    printf("DOING %g %g\n",data->aa[j],data->bb[j]);
}

void simulateCal(dataStruct *data)
{
  int i,j;
  int nbin,nchan;
  int iOn,iOff;
  double stokes_act[4],obsact[4],stokes_obs[4];
  double tvect[4];
  double scale=10;
  //  long seed = TKsetSeed();
  double simWhiteAA,simWhiteBB,simWhiteCR,simWhiteCI;
  double chanbw,noise;
  double I,Q,U,V;
  double calAmp=10;
  double zeroNoise=5;
  // The cal is affected by differential gain and phase and cross coupling
  
  data->period   = 1.0/11.123L;  // MUST FIX THIS

  iOn = (int)(data->nbin*data->cal_phs+0.5);
  iOff = (int)(data->nbin*data->cal_phs+data->cal_dcyc*data->nbin+0.5);
  printf("Ion = %d, ioff = %d\n",iOn,iOff);
  nbin = data->nbin;
  nchan = data->nchan;

  chanbw = fabs(data->obsBW/(double)data->nchan);
  noise = data->receiverNoise/sqrt(data->tsub*chanbw);
  for (j=0;j<nchan;j++)
    {


      for (i=0;i<nbin;i++)
	{
	  data->aa[j*data->nbin+i] = TKgaussDev(&seed)*noise;
	  data->bb[j*data->nbin+i] = TKgaussDev(&seed)*noise;
	  data->cr[j*data->nbin+i] = TKgaussDev(&seed)*noise/sqrt(2); 
	  data->ci[j*data->nbin+i] = TKgaussDev(&seed)*noise/sqrt(2); 

	  //	  simulateOffPulseNoise(&simWhiteAA,&simWhiteBB,&simWhiteCR,&simWhiteCI,data,j,&seed);
	  //      stokes_act[0] = data->cal_amp;
	  //      stokes_act[1] = 0;
	  //      stokes_act[2] = data->cal_amp;
	  //      stokes_act[3] = 0;

	  //	  simulateFeed(data,stokes_act,stokes_obs,j);      
	  //	  convertStokesAABBCRCI(stokes_obs,obsact);

	  //	  data->aa[j*nbin+i] = simWhiteAA;
	  //	  data->bb[j*nbin+i] = simWhiteBB;
	  //	  data->cr[j*nbin+i] = simWhiteCR;
	  //	  data->ci[j*nbin+i] = simWhiteCI;

	  I = Q = U = V = 0;
	  
	  // Work out if the cal if on or not in this bin
	  /*	  if (i > iOn && i < iOff)
		  {
		  I = U = 1;
		  Q = V = 0;
		  }
		  stokes_act[0] = I;
		  stokes_act[1] = Q;
		  stokes_act[2] = U;
		  stokes_act[3] = V; 

		  simulateFeed(data,stokes_act,stokes_obs,j); */

	  if (i >= iOn && i < iOff)
	    {
	      double amp = calAmp;
	      // Overshoot
	      if (i == iOn) amp=calAmp*1.2;
	      data->aa[j*data->nbin+i] += amp/2.0;
	      data->bb[j*data->nbin+i] += amp/2.0;
	      data->cr[j*data->nbin+i] += amp/2.0;
	      data->ci[j*data->nbin+i] += 0.0;
	    }
	  convertAABBCRCI_stokes(data->aa[j*data->nbin+i],data->bb[j*data->nbin+i],
				 data->cr[j*data->nbin+i],data->ci[j*data->nbin+i],
				 stokes_act);

	  simulateFeed(data,stokes_act,stokes_obs,j);      
	  convertStokesAABBCRCI(stokes_obs,obsact);


	  data->aa[j*data->nbin+i] = 50 + data->absGain[j] * obsact[0];
	  data->bb[j*data->nbin+i] = 50 + data->absGain[j] * obsact[1];	  
	  data->cr[j*data->nbin+i] = data->absGain[j]*obsact[2];
	  data->ci[j*data->nbin+i] = data->absGain[j]*obsact[3];
		       
	  //	      data->aa[j*data->nbin+i]*=data->scaleAA*data->absGain[j];
	  //	      data->bb[j*data->nbin+i]*=data->scaleBB*data->absGain[j];
	  //	      data->cr[j*data->nbin+i]*=data->scaleCR*data->absGain[j];
	  //	      data->ci[j*data->nbin+i]*=data->scaleCI*data->absGain[j];

	}
    }

}

void writeNewSub(fitsfile *fptr,dataStruct *data,int subint,long double timeFromStart)
{
  int status=0;
  int dval_0=0;
  int dval_1=1;
  double chan_bw = (data->obsBW/data->nchan); 
  double tbin = (double)(data->period/data->nbin);
  int colnum;
  int totIn;

  if (strcmp(data->template_I,"NOT SET")==0 && data->nCompTemplate == 0)// Have full Stokes
    totIn=0;
  else
    totIn=1;
  printf("totIn = %d\n",totIn);
  fits_movnam_hdu(fptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  if (status) {fits_report_error(stdout,status); exit(1);}
  fits_update_key(fptr, TINT, (char *)"NAXIS2", &subint, NULL, &status );
  if (status) {fits_report_error(stdout,status); exit(1);}
  if (subint==1)
    {
      fits_update_key(fptr, TSTRING, (char *)"INT_TYPE", (char *)"TIME", NULL, &status );
      fits_update_key(fptr, TSTRING, (char *)"INT_UNIT", (char *)"SEC", NULL, &status );
      fits_update_key(fptr, TSTRING, (char *)"SCALE", (char *)"FluxDen", NULL, &status );
      if (totIn==0)
	fits_update_key(fptr, TSTRING, (char *)"POL_TYPE", (char *)"AABBCRCI", NULL, &status );
      else
	fits_update_key(fptr, TSTRING, (char *)"POL_TYPE", (char *)"INTEN", NULL, &status );
      fits_update_key(fptr,TINT, (char *)"NPOL",&(data->npol),NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NBIN",&(data->nbin),NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NBIN_PRD",&(data->nbin),NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"PHS_OFFS",&dval_0,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NBITS",&dval_1,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"ZERO_OFF",&dval_0,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NSUBOFFS",&dval_0,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NCHAN",&(data->nchan),NULL,&status);  
      fits_update_key(fptr,TDOUBLE, (char *)"CHAN_BW",&chan_bw,NULL,&status);  
      fits_update_key(fptr,TDOUBLE, (char *)"TBIN",&tbin,NULL,&status);  
      fits_update_key(fptr,TDOUBLE, (char *)"DM",&(data->dm),NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"RM",&dval_0,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NCHNOFFS",&dval_0,NULL,&status);  
      fits_update_key(fptr,TINT, (char *)"NSBLK",&dval_1,NULL,&status);  
      // We need to understand what this parameter really is
      //      fits_update_key(fptr,TSTRING, (char *)"EPOCHS",(char *)"VALID",NULL,&status);  
    }
  // Now write the information for this new subint
  {
    int indxval = 0;
    int nchan = data->nchan;
    int npol  = data->npol;
    int nbin  = data->nbin;
    long naxes[4];
    int naxis=3;
    float dat_freq[nchan],dat_wts[nchan],dat_offs[nchan*npol],dat_scl[nchan*npol];
    double rajd,decjd;
    double lst_sub,para;
    long double mjd;
    double f0;
    double tsub = data->tsub;
    double offsSub = (double)timeFromStart;
    int i,j;

    // MUST FIX
    //    rajd = 69.3158537729834; // Hardcoded to 0437
    //    decjd = -47.2523961499288; // Hardcoded to 0437
        rajd = 155.74168; // Hardcoded to 1022
        decjd = 10.03132; // Hardcoded to 1022

    // Now calculate the Local Sidereal Time (LST)
    mjd = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L + subint*tsub/86400.0;
    lst_sub = (double)calcLocalSiderealTime(mjd,data)*60.0*60.0;
    fits_get_colnum(fptr,CASEINSEN,"LST_SUB",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&lst_sub,&status);

    fits_get_colnum(fptr,CASEINSEN,"TSUBINT",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&tsub,&status);

    fits_get_colnum(fptr,CASEINSEN,"OFFS_SUB",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&offsSub,&status);

    fits_get_colnum(fptr,CASEINSEN,"RA_SUB",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&rajd,&status);

    fits_get_colnum(fptr,CASEINSEN,"DEC_SUB",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&decjd,&status);


    // NOTE THAT THIS SHOULD BE THE PAR_ANG AT THE SUBINT CENTRE -- CURRENTLY AT THE START
    //
    para = calculatePara(data);
    fits_get_colnum(fptr,CASEINSEN,"PAR_ANG",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&para,&status);

    // MUST CHECK IF THE POS ANGLE SHOULD BE INDENTICAL TO THE PAR_ANG
    fits_get_colnum(fptr,CASEINSEN,"POS_ANG",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&para,&status);

  f0 = data->cFreq + fabs(data->obsBW)/2.0; // Highest frequency
    
    for (i=0;i<nchan;i++)
      {
	dat_freq[i] = f0-fabs(data->obsBW/(double)data->nchan)*i; // Must fix
	dat_wts[i] = 1;
      }
    for (i=0;i<nchan*npol;i++)
      {
	dat_offs[i] = 0;
	dat_scl[i] = 1;
      }
    
    fits_get_colnum(fptr,CASEINSEN,"INDEXVAL",&colnum,&status);
    fits_report_error(stdout,status);
    fits_write_col(fptr,TINT,colnum,subint,1,1,&indxval,&status);
    fits_report_error(stdout,status); 

    // Write the data
    fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);  
    fits_modify_vector_len (fptr, colnum, (nchan*npol*nbin), &status); 
    if (status) {fits_report_error(stdout,status); exit(1);}

    naxes[0] = nbin;
    naxes[1] = nchan;
    naxes[2] = npol;
    fits_delete_key(fptr, "TDIM18", &status); // THIS SHOULD NOT BE HARDCODED
    fits_write_tdim(fptr, colnum, naxis, naxes, &status);
    
    fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
    // Calculate scaling parameters
    {
      double aa_mean,aa_min,aa_max;
      double bb_mean,bb_min,bb_max;
      double cr_mean,cr_min,cr_max;
      double ci_mean,ci_min,ci_max;
      double scaleAA,scaleBB,scaleCR,scaleCI;
      double offsAA,offsBB,offsCR,offsCI;
      for (j=0;j<nchan;j++)
	{
	  aa_mean=bb_mean=cr_mean=ci_mean=0.0;
	  for (i=0;i<nbin;i++)
	    {
	      if (i==0)
		{
		  aa_min = aa_max = data->aa[j*nbin+i];
		  bb_min = bb_max = data->bb[j*nbin+i];
		  cr_min = cr_max = data->cr[j*nbin+i];
		  ci_min = ci_max = data->ci[j*nbin+i];
		}
	      else
		{
		  if (aa_min > data->aa[j*nbin+i]) aa_min = data->aa[j*nbin+i];
		  if (aa_max < data->aa[j*nbin+i]) aa_max = data->aa[j*nbin+i];
		  
		  if (bb_min > data->bb[j*nbin+i]) bb_min = data->bb[j*nbin+i];
		  if (bb_max < data->bb[j*nbin+i]) bb_max = data->bb[j*nbin+i];
		  
		  if (cr_min > data->cr[j*nbin+i]) cr_min = data->cr[j*nbin+i];
		  if (cr_max < data->cr[j*nbin+i]) cr_max = data->cr[j*nbin+i];
		  
		  if (ci_min > data->ci[j*nbin+i]) ci_min = data->ci[j*nbin+i];
		  if (ci_max < data->ci[j*nbin+i]) ci_max = data->ci[j*nbin+i];
		  
		}
	      aa_mean += data->aa[j*nbin+i];
	      bb_mean += data->bb[j*nbin+i];
	      cr_mean += data->cr[j*nbin+i];
	      ci_mean += data->ci[j*nbin+i];
	    }
	  scaleAA = (aa_max-aa_min)/2.0/16384.0;
	  scaleBB = (bb_max-bb_min)/2.0/16384.0;
	  scaleCR = (cr_max-cr_min)/2.0/16384.0;
	  scaleCI = (ci_max-ci_min)/2.0/16384.0;
	  offsAA = aa_max-16384.0*scaleAA;	  
	  if (totIn==0)
	    {
	      offsBB = bb_max-16384.0*scaleBB;
	      offsCR = cr_max-16384.0*scaleCR;
	      offsCI = ci_max-16384.0*scaleCI;
	    }
	  for (i=0;i<nbin;i++)
	    {
	      //	      printf("Have %g %g %g %g %g %g %g %g (%g %g %g %g)\n",scaleAA,scaleBB,scaleCR,scaleCI,offsAA,offsBB,offsCR,offsCI,data->aa[j*nbin+i],data->bb[j*nbin+i],data->cr[j*nbin+i],data->ci[j*nbin+i]);
	      /*	      data->aa[j*nbin+i] = (data->aa[j*nbin+i]-aa_mean/(double)nbin)*scaleAA;
	      data->bb[j*nbin+i] = (data->bb[j*nbin+i]-bb_mean/(double)nbin)*scaleBB;
	      data->cr[j*nbin+i] = (data->cr[j*nbin+i]-cr_mean/(double)nbin)*scaleCR;
	      data->ci[j*nbin+i] = (data->ci[j*nbin+i]-ci_mean/(double)nbin)*scaleCI;*/
	      //
	      data->aa[j*nbin+i] = (data->aa[j*nbin+i]-offsAA)/scaleAA;
	      if (totIn==0)
		{
		  data->bb[j*nbin+i] = (data->bb[j*nbin+i]-offsBB)/scaleBB;
		  data->cr[j*nbin+i] = (data->cr[j*nbin+i]-offsCR)/scaleCR;
		  data->ci[j*nbin+i] = (data->ci[j*nbin+i]-offsCI)/scaleCI;
		}
	      //	      printf("Have %g %g %g %g %g\n",data->aa[j*nbin+i],data->bb[j*nbin+i],scaleAA,scaleBB,bb_mean);
	    }	
	  if (totIn==0)
	    {
	      // Note strange ordering of frequency and pol information
	      dat_scl[j]   = scaleAA; dat_offs[j]   = offsAA; //scaleAA*aa_mean/(double)nbin;
	      dat_scl[j+nchan] = scaleBB; dat_offs[j+nchan] = offsBB; //scaleBB*bb_mean/(double)nbin;
	      dat_scl[j+2*nchan] = scaleCR; dat_offs[j+2*nchan] = offsCR; //scaleCR*cr_mean/(double)nbin;
	      dat_scl[j+3*nchan] = scaleCI; dat_offs[j+3*nchan] = offsCI; //scaleCI*ci_mean/(double)nbin;
	    }
	  else
	    {
	      dat_scl[j]   = scaleAA; dat_offs[j]   = scaleAA*aa_mean/(double)nbin;
	    }
	}
    }
      //      printf("Writing %g\n",data->aa[i]);

    if (totIn==0){
      float vals[nbin*nchan*4];
      for (j=0;j<nchan;j++)
	{
	  for (i=0;i<nbin;i++)
	    {
	      vals[j*nbin*4+i] = data->aa[j*nbin+i];
	      vals[j*nbin*4+nbin+i] = data->bb[j*nbin+i];
	      vals[j*nbin*4+2*nbin+i] = data->cr[j*nbin+i];
	      vals[j*nbin*4+3*nbin+i] = data->ci[j*nbin+i];
	    }
	}
      //     fits_write_col(fptr,TFLOAT,colnum,subint,1,nbin*nchan*4,vals,&status);
            fits_write_col(fptr,TFLOAT,colnum,subint,1,nbin*nchan,data->aa,&status);
	    fits_write_col(fptr,TFLOAT,colnum,subint,nbin*nchan+1,nbin*nchan,data->bb,&status);
	    fits_write_col(fptr,TFLOAT,colnum,subint,2*nbin*nchan+1,nbin*nchan,data->cr,&status);
	    fits_write_col(fptr,TFLOAT,colnum,subint,3*nbin*nchan+1,nbin*nchan,data->ci,&status);
    }
    else
      fits_write_col(fptr,TFLOAT,colnum,subint,1,nbin*nchan,data->aa,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan, &status); 
    fits_write_col(fptr,TFLOAT,colnum,subint,1,nchan,dat_freq,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan, &status); 
    fits_write_col(fptr,TFLOAT,colnum,subint,1,nchan,dat_wts,&status);

    fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan*npol, &status); 
    fits_write_col(fptr,TFLOAT,colnum,subint,1,nchan*npol,dat_offs,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan*npol, &status); 
    fits_write_col(fptr,TFLOAT,colnum,subint,1,nchan*npol,dat_scl,&status);


  }
}

// Routine to create the output FITS file
void createFitsFile(fitsfile *fptr,dataStruct *data)
{
  int status=0;

  // Remove unnecessary tables
  removeTables(fptr,data);
  writeHeaderParameters(fptr,data->primaryHeaderParams,data);
}

void removeTables(fitsfile *fptr,dataStruct *data)
{
  int newHDUtype;
  int status=0;

  if (data->obsMode != 3)
    {fits_movnam_hdu(fptr, BINARY_TBL, "FLUX_CAL", 0, &status);  fits_delete_hdu(fptr, &newHDUtype, &status);}
  fits_movnam_hdu(fptr, BINARY_TBL, "COHDDISP", 0, &status);  fits_delete_hdu(fptr, &newHDUtype, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, "POLYCO", 0, &status);  fits_delete_hdu(fptr, &newHDUtype, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, "CAL_POLN", 0, &status);  fits_delete_hdu(fptr, &newHDUtype, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, "FEEDPAR", 0, &status);  fits_delete_hdu(fptr, &newHDUtype, &status);
  if (status)
    {
      fits_report_error(stdout,status);
      exit(1);
    }
}

void writeHeaderParameters(fitsfile *fptr,char *fname,dataStruct *data)
{
  FILE *fin;
  char keyword[128],setVal[128];
  double lst;
  long double mjd;
  int status=0;
  double obsBW;

  fits_movabs_hdu( fptr, 1, NULL, &status );
  fits_write_date(fptr, &status);
  
  if (!(fin = fopen(fname,"r")))
    {
      printf("Error: Unable to open file >%s<\n",fname);
      exit(1);
    }

  while (!feof(fin))
    {
      if (fscanf(fin,"%s %s",keyword,setVal)==2)
	{
	  //	  printf("Read: %s %s\n",keyword,setVal);
	  fits_update_key(fptr, TSTRING, keyword, setVal, NULL, &status );
	  fits_report_error(stdout,status);
	}
    }
  if (data->obsMode==1)
    {
      fits_update_key(fptr, TSTRING, (char *)"OBS_MODE", (char *)"LEVCAL", NULL, &status );
      fits_update_key(fptr, TSTRING, (char *)"CAL_MODE", (char *)"SYNC", NULL, &status );
      fits_update_key(fptr,TDOUBLE, (char *)"CAL_FREQ",&(data->cal_freq),NULL,&status);
      fits_update_key(fptr,TDOUBLE, (char *)"CAL_DCYC",&(data->cal_dcyc),NULL,&status);
      fits_update_key(fptr,TDOUBLE, (char *)"CAL_PHS",&(data->cal_phs),NULL,&status);
    }
  else
    {
      fits_update_key(fptr, TSTRING, (char *)"OBS_MODE", (char *)"PSR", NULL, &status );
      fits_update_key(fptr, TSTRING, (char *)"CAL_MODE", (char *)"OFF", NULL, &status );
    }
  fits_report_error(stdout,status);

  fits_update_key(fptr,TSTRING, (char *)"SRC_NAME",&(data->srcName),NULL,&status);
  fits_update_key(fptr,TDOUBLE, (char *)"OBSFREQ",&(data->cFreq),NULL,&status);
  obsBW = fabs(data->obsBW);
  fits_update_key(fptr,TDOUBLE, (char *)"OBSBW",&(obsBW),NULL,&status);
  fits_update_key(fptr,TINT, (char *)"OBSNCHAN",&(data->nchan),NULL,&status);
  fits_update_key(fptr,TDOUBLE, (char *)"SCANLEN",&(data->scanLength),NULL,&status);

  fits_update_key(fptr,TINT, (char *)"STT_IMJD",&(data->stt_imjd),NULL,&status);
  fits_update_key(fptr,TDOUBLE, (char *)"STT_SMJD",&(data->stt_smjd),NULL,&status);
  fits_update_key(fptr,TDOUBLE, (char *)"STT_OFFS",&(data->stt_offs),NULL,&status);
  fclose(fin);

  // Now calculate the Local Sidereal Time (LST)
  mjd = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L;
  lst = (double)calcLocalSiderealTime(mjd,data)*60.0*60.0;
  fits_update_key(fptr,TDOUBLE, (char *)"STT_LST",&lst,NULL,&status);

  // SHOULD FIX THIS
  // HARDCODE FOR 0437
  {
    //    char raj[128] = "04:37:15.810";
    //    char decj[128] = "-47:15:08.600";
    char raj[128] = "10:22:15.810";
    char decj[128] = "+10:01:08.600";
    fits_update_key(fptr, TSTRING, (char *)"RA", &raj, NULL, &status );
    fits_update_key(fptr, TSTRING, (char *)"DEC", &decj, NULL, &status );
  }

}
void writePredictor(fitsfile *fptr,char *fname)
{
  FILE *fin;
  char line[128];
  int colnum;
  int arr=1;
  int status=0;
  int linenum=41;
  char *temp = &(line[0]);
  fits_movnam_hdu(fptr,BINARY_TBL,(char *)"T2PREDICT",0,&status);
  if (status) {fits_report_error(stdout,status); exit(1);}
  linenum=1;
  fits_get_colnum(fptr,CASEINSEN,"PREDICT",&colnum,&status);
  fits_report_error(stdout,status);

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open file >%s<\n",fname);
      exit(1);      
    }
  while (!feof(fin))
    {
      if (fgets(line,128,fin))
	{
	  line[strlen(line)-1]='\0';
	  fits_write_col(fptr,TSTRING,colnum,linenum++,1,1,&temp,&status);
	  if (status) {fits_report_error(stdout,status);  exit(1);}
	}
    }
  fclose(fin);
}

void writeEphemeris(fitsfile *fptr,dataStruct *data)
{
  FILE *fin;
  char line[128];
  int colnum;
  int arr=1;
  int status=0;
  int linenum=41;
  char word1[128],word2[128];
  char *temp = &(line[0]);
  fits_movnam_hdu(fptr,BINARY_TBL,(char *)"PSRPARAM",0,&status);
  if (status) {fits_report_error(stdout,status); exit(1);}
  linenum=1;
  fits_get_colnum(fptr,CASEINSEN,"PARAM",&colnum,&status);
  fits_report_error(stdout,status);

  if (!(fin = fopen(data->exact_ephemeris,"r")))
    {
      printf("Unable to open file >%s<\n",data->exact_ephemeris);
      exit(1);      
    }
  while (!feof(fin))
    {
      if (fgets(line,128,fin))
	{
	  line[strlen(line)-1]='\0';
	  fits_write_col(fptr,TSTRING,colnum,linenum++,1,1,&temp,&status);
	  if (status) {fits_report_error(stdout,status);  exit(1);}
	  if (sscanf(line,"%s %s",word1,word2)==2)
	    {
	      if (strcasecmp(word1,"DM")==0)
		sscanf(word2,"%lf",&(data->dm));
	    }
	}
    }
  fclose(fin);
}

void convertStokesAABBCRCI(double *stokes_act,double *obsact)
{
  double I,Q,U,V,A,B,C,D;

  I = stokes_act[0];
  Q = stokes_act[1];
  U = stokes_act[2];
  V = stokes_act[3];

  A = 0.5*(I+Q);
  B = 0.5*(I-Q);
  C = U/2.0;
  D = V/2.0;

  obsact[0] = A;
  obsact[1] = B;
  obsact[2] = C;
  obsact[3] = D;

}

void createM_amp(double gain,double phase,double m_amp[4][4])
{
  m_amp[0][0] = 1.0;  m_amp[0][1] = gain/2.0; m_amp[0][2] = 0.0;        m_amp[0][3] = 0.0;
  m_amp[1][0] = gain/2.0; m_amp[1][1] = 1;    m_amp[1][2] = 0.0;        m_amp[1][3] = 0.0;
  m_amp[2][0] = 0.0;  m_amp[2][1] = 0.0;  m_amp[2][2] = cos(phase); m_amp[2][3] = -sin(phase);
  m_amp[3][0] = 0.0;  m_amp[3][1] = 0.0;  m_amp[3][2] = sin(phase); m_amp[3][3] = cos(phase);

}

void createM_feed(double gamma,double m_feed[4][4])
{
  gamma*=M_PI/180.0;
  m_feed[0][0] = 1.0;  m_feed[0][1] = 0.0;  m_feed[0][2] = 0.0;        m_feed[0][3] = 0.0;
  m_feed[1][0] = 0.0;  m_feed[1][1] = cos(2*gamma);  m_feed[1][2] = 0.0;        m_feed[1][3] = sin(2*gamma);
  m_feed[2][0] = 0.0;  m_feed[2][1] = 0.0;  m_feed[2][2] = 1.0;        m_feed[2][3] = 0.0;
  m_feed[3][0] = 0.0;  m_feed[3][1] = -sin(2*gamma);  m_feed[3][2] = 0.0;        m_feed[3][3] = cos(2*gamma);
}

void createM_pa(double pa,double m_pa[4][4])
{
  pa = pa*M_PI/180.0;
  m_pa[0][0] = 1.0; m_pa[0][1] = 0; m_pa[0][2] = 0; m_pa[0][3] = 0;
  m_pa[1][0] = 0.0; m_pa[1][1] = cos(2*pa); m_pa[1][2] = sin(2*pa); m_pa[1][3] = 0;
  m_pa[2][0] = 0.0; m_pa[2][1] = -sin(2*pa); m_pa[2][2] = cos(2*pa); m_pa[2][3] = 0;
  m_pa[3][0] = 0.0; m_pa[3][1] = 0; m_pa[3][2] = 0; m_pa[3][3] = 1;
}

void createM_cc(double cc_eps1,double cc_eps2,double cc_phi1,double cc_phi2,double m_cc[4][4])
{
  double A,B,C,D;

  A = cc_eps1*cos(cc_phi1) + cc_eps2*cos(cc_phi2);
  B = cc_eps1*sin(cc_phi1) + cc_eps2*sin(cc_phi2);
  C = cc_eps1*cos(cc_phi1) - cc_eps2*cos(cc_phi2);
  D = cc_eps1*sin(cc_phi1) - cc_eps2*sin(cc_phi2);

  m_cc[0][0] = 1.0;  m_cc[0][1] = 0.0; m_cc[0][2] = A; m_cc[0][3] = B;
  m_cc[1][0] = 0.0;  m_cc[1][1] = 1.0; m_cc[1][2] = C; m_cc[1][3] = D;
  m_cc[2][0] = A;    m_cc[2][1] = -C;  m_cc[2][2] = 1; m_cc[2][3] = 0;
  m_cc[3][0] = B;    m_cc[3][1] = -D;  m_cc[3][2] = 0; m_cc[3][3] = 1;
}

void multVectMat(double m[4][4],double *vIn,double *vOut)
{
  int i,j;
  for (i=0;i<4;i++)
    {
      vOut[i]=0.0;
      for (j=0;j<4;j++)
	vOut[i] += m[i][j]*vIn[j];
    }
}

int readTemplateDf(char *fname,double *concentration,double *centre,double *height,
		   double *concentration_df,double *centre_df,double *height_df,double *fiducial_f,int nbin)
{
  FILE *fin;
  int n=0;
  float tmp;
  int nread;
  char line[1024];
  char w1[1024],w2[1024],w3[1024];
  int format=0;

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open template: %s\n",fname);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fgets(line,1024,fin)!=NULL)
	{
	  if (line[0] == '#')
	    {
	      sscanf(line,"%s %s %s\n",w1,w2,w3);
	      if (strcmp(w2,"FORMAT")==0 && strcmp(w3,"1")==0)
		format=1;
	    }
	  else
	    {
	      nread = sscanf(line,"%lf %lf %lf %lf %lf %lf %lf",&centre[n],&concentration[n],&height[n],
			     &centre_df[n],&concentration_df[n],&height_df[n],&fiducial_f[n]);
	      if (format==0)
		centre[n] *= nbin;
	      
	      if (nread==7)
		n++;
	      else if (nread==3)
		{
		  //	      printf("Read in %g %g %g\n",centre[n],concentration[n],height[n]);
		  centre_df[n] = 0;
		  concentration_df[n] = 0;
		  height_df[n] = 0;
		  fiducial_f[n] = 0;
		  n++;
		}
	    }
	}
    }
  fclose(fin);
  //  printf("Returning n = %d\n",n);
  return n;

}

int readTemplate(char *fname,double *concentration,double *centre,double *height)
{
  FILE *fin;
  int n=0;
  float tmp;

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open template: %s\n",fname);
      exit(1);
    }
  while (!feof(fin))
    {
      //      if (fscanf(fin,"%lf %lf %lf",&centre[n],&concentration[n],&height[n])==3)
      if (fscanf(fin,"%lf %lf %lf",&height[n],&concentration[n],&centre[n])==3)
	n++;
    }
  fclose(fin);
  return n;
}

void runTempo2(dataStruct *data)
{
  char execString[1024];
  double seg_length = 48000; // MUST FIX
  int nfreqcoeff = 2; // MUST FIX
  //  int ntimecoeff = 12;  
  int ntimecoeff = 24;   // MUST FIX

  double freq1 = data->cFreq - fabs(data->obsBW); 
  double freq2 = data->cFreq + fabs(data->obsBW);
  long double mjd1 = data->stt_imjd + data->stt_smjd/86400.0L - 60*60/86400.0L; // MUST FIX
  long double mjd2 = data->stt_imjd + data->stt_smjd/86400.0L + (data->nsub)*data->tsub/86400.0 + 60*60/86400.0L; // MUST FIX

  sprintf(execString,"tempo2 -pred \"PKS %Lf %Lf %g %g %d %d %g\" -f %s",mjd1,mjd2,freq1,freq2,ntimecoeff,nfreqcoeff,seg_length,data->exact_ephemeris);
  printf("Running tempo2 to get predictor\n");
  system(execString);
  printf("Complete running tempo2\n");
}

void simulateFeed(dataStruct *data,double *stokes_act,double *stokes_obs,int fchan)
{
  double m_amp[4][4]; 
  double m_cc[4][4]; // Cross coupling Mueller matrix
  double mul[4][4];  // Mueller matrix
  double m_feed[4][4];
  double temp[4][4],temp2[4][4];
  double m_pa[4][4]; // Parallactic angle Mueller matrix
  double pa;
  double tvect[4],tvect2[4],tvect3[4];

  createM_amp(data->diffGain[fchan],data->diffPhase[fchan],m_amp);
  createM_cc(data->cc_eps1[fchan],data->cc_eps2[fchan],data->cc_phi1[fchan],data->cc_phi2[fchan],m_cc);
  createM_feed(0,m_feed);

  if (data->obsMode==2)
    {
      createM_pa(data->pa,m_pa);    
      //      mult4x4(m_pa,m_cc,temp);
      //      mult4x4(temp,m_amp,mul);


      mult4x4(m_amp,m_cc,temp);
      mult4x4(temp,m_feed,temp2);
      mult4x4(temp2,m_pa,mul);

      /*      multVectMat(m_pa,stokes_act,tvect);
	      multVectMat(m_feed,tvect,tvect2);
	      multVectMat(m_cc,tvect2,tvect3);
	      multVectMat(m_amp,tvect3,stokes_obs); */
    }
  else if (data->obsMode==1)
    {
      /*      multVectMat(m_feed,stokes_act,tvect);
      multVectMat(m_cc,tvect,tvect2);      
      multVectMat(m_amp,tvect2,stokes_obs);*/

            mult4x4(m_amp,m_cc,temp);
            mult4x4(temp,m_feed,mul);
      //    multVectMat(m_amp,stokes_act,stokes_obs);
    }


      multVectMat(mul,stokes_act,stokes_obs);
}

void mult4x4(double a[4][4],double b[4][4],double c[4][4])
{
  int i,j,k;
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	{
	  c[i][j]=0.0;
	  for (k=0;k<4;k++)
	      c[i][j]+=a[i][k]*b[k][j];
	}
    }
}

// Calculate parallactic angle in degrees
double calculatePara(dataStruct *data)
{
  long double mjd,last,ha,latitude,dec,ra;
  double pa;

  // MUST FIX - HARDCODED TO PARKES
  latitude = -32.0-59.0/60.0-54.263/60.0/60.0; // Degrees 
  // MUST FIX - HARDCODED TO J0437-4715
  //  ra  = 4+37.0/60.0+15.883250/60.0/60.0; // Hours
  //  dec = -47-15/60.0-9.031863/60.0/60.0; // Degrees

  ra  = 10+37.0/60.0+15.883250/60.0/60.0; // Hours
  dec = 10-15/60.0-9.031863/60.0/60.0; // Degrees
 
  mjd = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L;
  printf("mjd = %.15f\n",(double)mjd);
  last = calcLocalSiderealTime(mjd,data);
  printf("last = %.15f\n",(double)last);
  ha = last - ra;
  printf("Hour angle = %g\n",(double)ha);
  ha = ha*180.0/12;
  pa = atan2(sin(ha*M_PI/180.0)*cos(latitude*M_PI/180.0),sin(latitude*M_PI/180.0)*cos(dec*M_PI/180.0)-cos(latitude*M_PI/180.0)*sin(dec*M_PI/180.0)*cos(ha*M_PI/180.0));
  printf("Para angle = %g rad = %g deg\n",pa,pa*180.0/M_PI);

  return pa*180.0/M_PI;
}

long double calcLocalSiderealTime(long double mjd,dataStruct *data)
{
  long double JD,D,JD0,D0,gmst,H,T,gast,eps,L,Omega,Psi,E,last;
  long double longitude;

  printf("MJD = %.15f\n",(double)mjd);
  JD = mjd + 2400000.5;
  JD0 = floor(mjd) + 2400000.5;


  D = JD - 2451545.0L;
  D0 = JD0 - 2451545.0L;
  T = D/36525.0L;
  H = (JD-JD0)*24.0L;

  printf("H = %g\n",(double)H);

  gmst = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T;
  printf("gmst = %g\n",(double)gmst);
  // MUST REDUCE TO THE RANGE of 0h to 24h
  while (gmst < 0 || gmst > 24)
    {
      if (gmst < 0) gmst += 24;
      else if (gmst > 24) gmst -= 24;
    }
  printf("gmst now = %g\n",(double)gmst);
  eps = 23.4393 - 0.0000004*D;
  L = 280.47+0.98565*D;
  Omega = 125.04 - 0.052954*D;
  Psi = -0.000319*sin(Omega*M_PI/180.0)-0.000024*sin(2*L*M_PI/180.0);
  E = Psi*cos(eps*M_PI/180.0);
  gast = gmst + E;
  printf("gast = %.5f\n",(double)gast);

  // TO FIX --- specific for Parkes
  longitude = 148.0+15/60.0+48.636/60.0/60.0; // Degrees East of Greenwich
  last = gast + longitude/15.0;
  // MUST REDUCE TO THE RANGE of 0h to 24h
  /*  while (last < 0 || last > 24)
    {
      if (last < 0) last += 24;
      else if (last > 24) last -= 24;
      }*/

  printf("Local sidereal time = %g\n",(double)last);

  return last;
}

void readFeedParameters(dataStruct *data)
{
  FILE *fin;
  int i;

  if (!(fin = fopen(data->feedParamFile,"r")))
    {
      printf("Unable to open file >%s< for the feed parameters\n",data->feedParamFile);
      for (i=0;i<data->nchan;i++)
	{
	  data->absGain[i] = 1;
	  data->diffGain[i] = 0;
	  data->diffPhase[i] = 0;
	  data->cc_eps1[i] = 0;
	  data->cc_phi1[i] = 0;
	  data->cc_eps2[i] = 0;
	  data->cc_phi2[i] = 0;
	}
    }
  else
    {
      for (i=0;i<data->nchan;i++)
	{
	  fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf",
		 &(data->absGain[i]),&(data->diffGain[i]),&(data->diffPhase[i]),
		 &(data->cc_eps1[i]),&(data->cc_phi1[i]),&(data->cc_eps2[i]),
		 &(data->cc_phi2[i]));
	}
      fclose(fin);
    }
}

void simulateOffPulseNoise(double *simWhiteAA,double *simWhiteBB,double *simWhiteCR,double *simWhiteCI,dataStruct *data,int chan,long *seed)
{
  double simI = TKgaussDev(seed);
  double stokes[4],stokes_obs[4];
  //  double scaleFactor = 0.1;
  //  double scaleFactor = 0.01;
  double obsact[4];
  double chanbw = fabs(data->obsBW/(double)data->nchan);
  double noise = data->receiverNoise/sqrt(data->tsub*chanbw);

  //  stokes[0] = 5+(simI*scaleFactor);
  // WHAT SHOULD THIS ACTUALLY BE!  WORK ON THIS
  stokes[0] = noise*simI; 
  stokes[1] = noise*simI;
  stokes[2] = noise*simI;
  stokes[3] = noise*simI;

  simulateFeed(data,stokes,stokes_obs,chan);      
  convertStokesAABBCRCI(stokes_obs,obsact);
  *simWhiteAA = obsact[0];
  *simWhiteBB = obsact[1];
  *simWhiteCR = obsact[2];
  *simWhiteCI = obsact[3];
  //  printf("Res = %g %g %g %g\n",*simWhiteAA,*simWhiteBB,*simWhiteCR,*simWhiteCI);
}


void readScintillationFile(float **scint)
{
  FILE *fin;
  int nx=16384;
  int ny=1024;
  int i,j;
  float *flt;
  flt = (float *)malloc(sizeof(float)*nx*2); // Read in real and imaginary parts of the e-field 

  if (!(fin = fopen("strong10w.spe","rb")))
    {
      printf("Unable to read scintillation file: strong10w.spe\n");
	exit(1);
    }
  // Read header information
  fread(flt,sizeof(float),nx*2,fin);  
  for (j=0;j<ny;j++)
    {  
      fread(flt,sizeof(float),nx*2,fin);  
      for (i=0;i<nx;i++)
	{
	  scint[j][i] = pow(flt[2*i],2)+pow(flt[2*i+1],2);
	}
    }
  fclose(fin);
  free(flt);
}

void convertAABBCRCI_stokes(double aa,double bb,double cr,double ci,double *stokes)
{
  stokes[0] = aa+bb;
  stokes[1] = aa-bb;
  stokes[2] = 2*cr;
  stokes[3] = 2*ci;
}
