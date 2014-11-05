/* Compilation: ./compile2

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "simulatePTA.h"

int main(int argc,char *argv[])
{
  dataStruct *data;
  fitsfile *fptr;
  int i;
  int havePrimaryHeaderParams=0;
  char obsTable[128]="obsTable";
  FILE *fin;
  char line[128];
  char temp[128];
  char file[128];
  int status=0;
  int subint;
  int processFast=0;

  data = (dataStruct *)malloc(sizeof(dataStruct));
  initialiseData(data);

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-obsTable")==0)
	strcpy(obsTable,argv[++i]);
      else if (strcmp(argv[i],"-fast")==0)
	processFast=1;
      else
	{
	  printf("Unknown command line argument: %s\n",argv[i]);
	}
    }
  // Check command line arguments
  if (!(fin = fopen(obsTable,"r")))
    {
      printf("Unable to open file >%s<\n",obsTable);
      printf("Must provide an observation file using -obsTable\n");
      exit(1);
    }
  // Start the simulation
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
	    fscanf(fin,"%lf",&(data->tsub));
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
	    {
	      fscanf(fin,"%lf",&(data->cal_freq));	    
	      data->period = 1.0/data->cal_freq;
	    }
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
	      // Now open the fits files
	      // Overwrite if already existing
	      sprintf(file,"!%s(psrheader.fits)",data->fname);
	      fits_create_file(&fptr,file,&status);
	      fits_report_error(stdout,status);
	      // Write all the primary header information
	      createFitsFile(fptr,data);
	      // For a pulsar observation: write the ephemeris into the file
	      if (data->obsMode == 2)
		writeEphemeris(fptr,data);

	      subint=1;
	      // Read feed parameteres
	      readFeedParameters(data); 
	      // For a pulsar observation: write the predictor file
	      if (data->obsMode == 2)
		{
		  T2Predictor pred;
		  int ret;
		  long double phase0,mjd0,freq,f0;

		  T2Predictor_Init(&pred);
		  runTempo2(data);
		  writePredictor(fptr,"t2pred.dat");
		  // Calculate phase offset for the pulse arrival time
		  if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
		    {
		      printf("Error: unable to read predictor\n");
		      exit(1);
		    }
		  mjd0 = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L;

		  for (i=0;i<data->nchan;i++)
		    {
		      f0 = data->cFreq + fabs(data->obsBW)/2.0; // Highest frequency
		      freq = f0-fabs(data->obsBW/(double)data->nchan)*i;
		      phase0 = T2Predictor_GetPhase(&pred,mjd0,freq);
		      data->phaseOffset[i] = (phase0-floor(phase0));
		    }
		  //		  printf("Phase %.16f %.16Lf\n",data->phaseOffset,phase0);
		  data->period      = 1.0/T2Predictor_GetFrequency(&pred,mjd0,freq);
		  T2Predictor_Destroy(&pred);
		  }

	      // Create and write a new subintegration
	      simulateData(data,processFast);
	      writeNewSub(fptr,data,subint);

	      // Close this file and get ready for the next simulation
	      fits_close_file(fptr,&status);
	      initialiseData(data);
	    }
	}
    }
  fclose(fin);
  printf("Goodbye\n");
  free(data);
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
  data->scanLength = 60;
  data->npol     = 4;
  data->scaleAA  = 1;
  data->scaleBB  = 1;
  data->scaleCR  = 1;
  data->scaleCI  = 1;
  data->whiteNoiseI = 1;
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


  {
    char raj[128]; 
    char decj[128];
    strcpy(raj,data->rajStr);
    strcpy(decj,data->decjStr);
    fits_update_key(fptr, TSTRING, (char *)"RA", &raj, NULL, &status );
    fits_update_key(fptr, TSTRING, (char *)"DEC", &decj, NULL, &status );
  }
}


void simulateData(dataStruct *data,int processFast)
{
  double tobs = 1;
  double period = data->period;
  int nbin = data->nbin;
  double tbin = period/nbin;
  long count,n;
  int nsamp = 1000; // This should be how much we're averaging over
  cmplx ex_noise[nsamp];
  cmplx ey_noise[ nsamp];
  long seed = TKsetSeed();
  double ax = 1;
  double ay = 1;
  float aa[MAX_NBIN],bb[MAX_NBIN],cr[MAX_NBIN],ci[MAX_NBIN];
  float si[MAX_NBIN],sq[MAX_NBIN],su[MAX_NBIN],sv[MAX_NBIN];
  cmplx psrx,psry;
  int i,j,k,l;
  double calval;
  cmplx cal;
  int calOn = (int)(data->nbin*data->cal_phs+0.5);
  int calOff = (int)(data->nbin*data->cal_phs+data->cal_dcyc*data->nbin+0.5);;
  double ii[MAX_NBIN],q[MAX_NBIN],u[MAX_NBIN],v[MAX_NBIN],psr_i[MAX_NBIN];
  double psr_stokes[MAX_NBIN][4];
  double sig_stokes[MAX_NBIN][4];
  double sig_i[MAX_NBIN],sig_q[MAX_NBIN],sig_u[MAX_NBIN],sig_v[MAX_NBIN];
  double aabbcrci[4];
  double ranVal;
  long double x0min;
  double m_amp[4][4]; 
  double m_cc[4][4]; // Cross coupling Mueller matrix
  double m_feed[4][4];
  double tvect[4],tvect2[4],tvect3[4];
  double pval;

  printf("nsamp = %d\n",nsamp);

  data->pa = calculatePara(data);
  printf("Parallactic angle = %g\n",data->pa);
  for (k=0;k<data->nchan;k++)
    {
      x0min = data->phaseOffset[k]*(double)data->nbin;
      printf("Simulating channel #%d\n",k+1);
      for (l=0;l<100;l++)
      //      l=0;
	{
	  printf("k = %d l= %d\n",k,l);
	  if (data->obsMode==2)
	    {
	      int nCompII,nCompQ,nCompU,nCompV;
	      double concentrationQ[128],centreQ[128],heightQ[128];
	      double concentrationU[128],centreU[128],heightU[128];
	      double concentrationV[128],centreV[128],heightV[128];
	      double concentrationII[128],centreII[128],heightII[128];
	      double m_pa[4][4];
	      
	      nCompII = readTemplate(data->template_II,concentrationII,centreII,heightII);
	      nCompQ = readTemplate(data->template_Q,concentrationQ,centreQ,heightQ);
	      nCompU = readTemplate(data->template_U,concentrationU,centreU,heightU);
	      nCompV = readTemplate(data->template_V,concentrationV,centreV,heightV);
	      
	      
	      for (i=0;i<nbin;i++)
		{
		  ii[i] = q[i] =u[i] = v[i] = 0.0;
		  for (j=0;j<nCompII;j++)
		    ii[i] += heightII[j]*exp(concentrationII[j]*(cos((i+x0min-centreII[j])/(double)data->nbin*2*M_PI)-1));
		  for (j=0;j<nCompQ;j++)
		    q[i] += heightQ[j]*exp(concentrationQ[j]*(cos((i+x0min-centreQ[j])/(double)data->nbin*2*M_PI)-1));
		  for (j=0;j<nCompU;j++)
		    u[i] += heightU[j]*exp(concentrationU[j]*(cos((i+x0min-centreU[j])/(double)data->nbin*2*M_PI)-1));
		  for (j=0;j<nCompV;j++)
		    v[i] += heightV[j]*exp(concentrationV[j]*(cos((i+x0min-centreV[j])/(double)data->nbin*2*M_PI)-1));
		  psr_stokes[i][0] = sqrt(pow(ii[i],2)+pow(q[i],2)+pow(u[i],2)+pow(v[i],2));
		  psr_stokes[i][1] = q[i];
		  psr_stokes[i][2] = u[i];
		  psr_stokes[i][3] = v[i];
		  
		  //
		  // Now account for parallactic angle 
		  //
		  createM_pa(data->pa,m_pa);    
		  multVectMat(m_pa,psr_stokes[i],tvect);
		  psr_stokes[i][0] = tvect[0];
		  psr_stokes[i][1] = tvect[1];
		  psr_stokes[i][2] = tvect[2];
		  psr_stokes[i][3] = tvect[3];
		  
		}
	    }
	  else
	    {
	      for (i=0;i<nbin;i++)
		{
		  psr_stokes[i][0] = 0.0;
		  psr_stokes[i][1] = 0.0;
		  psr_stokes[i][2] = 0.0;
		  psr_stokes[i][3] = 0.0;
		}
	    }
	  
	  
	  for (i=0;i<nbin;i++)
	    {
	      if (data->obsMode==1)
		{
		  if (i>= calOn && i<calOff)
		    {
		      if (i==calOn)
			calval=0.2; // 1.50; // overshoot  
		      else
			calval=0.1; //1.42;
		    }
		  else
		    calval=0;
		}
	      else
		calval=0;
	      
	      if (processFast==1)
		{
		  sig_stokes[i][0] = TKgaussDev(&seed)*(data->receiverNoise) + psr_stokes[i][0] +calval*data->receiverNoise;
		  sig_stokes[i][1] = TKgaussDev(&seed)*data->receiverNoise + psr_stokes[i][1];
		  sig_stokes[i][2] = TKgaussDev(&seed)*(data->receiverNoise) + psr_stokes[i][2] +calval*data->receiverNoise;
		  sig_stokes[i][3] = TKgaussDev(&seed)*data->receiverNoise + psr_stokes[i][3];
		}
	      else
		{
		  // Simulate Ex and Ey for perfect feed
		  for (count=0;count<nsamp;count++)
		    {
		      if (data->obsMode==1)
			{
			  cal.r = 1.0/sqrt(2)*calval*TKgaussDev(&seed);
			  cal.i = 0;
			}
		      else
			{
			  cal.r = cal.i = 0;
			}
		      ex_noise[count].r = (ax)*TKgaussDev(&seed) + cal.r;
		      ex_noise[count].i = (ax)*TKgaussDev(&seed) + cal.i;
		      
		      ey_noise[count].r = (ay)*TKgaussDev(&seed) + cal.r;
		      ey_noise[count].i = (ay)*TKgaussDev(&seed) + cal.i;
		    }
		  
		  // Now detect the signal
		  aa[i] = 0.0;
		  bb[i] = 0.0;
		  cr[i] = 0.0;
		  ci[i] = 0.0;
		  for (count=0;count<nsamp;count++)
		    {
		      aa[i]+=pow(ex_noise[count].r,2) + pow(ex_noise[count].i,2);
		      bb[i]+=pow(ey_noise[count].r,2) + pow(ey_noise[count].i,2);
		      cr[i]+=(ex_noise[count].r*ey_noise[count].r + ex_noise[count].i*ey_noise[count].r);
		      ci[i]+=(ex_noise[count].i*ey_noise[count].r - ex_noise[count].r*ey_noise[count].i);
		    }
		  aa[i]/=(double)nsamp;
		  bb[i]/=(double)nsamp;
		  cr[i]/=(double)nsamp;
		  ci[i]/=(double)nsamp;
		  //	  printf("Using: %g %g %g %g\n",aa[i],bb[i],cr[i],ci[i]);
		  sig_stokes[i][0] = (aa[i]+bb[i])+psr_stokes[i][0];
		  sig_stokes[i][1] = (aa[i]-bb[i])+psr_stokes[i][1];
		  sig_stokes[i][2] = 2*cr[i] + psr_stokes[i][2];
		  sig_stokes[i][3] = 2*ci[i] + psr_stokes[i][3];
		  pval = sqrt(pow(sig_stokes[i][1],2)+pow(sig_stokes[i][2],2)+pow(sig_stokes[i][3],2));
		  //	      printf("overpol = %g %g %g %g %g %g\n",sig_stokes[i][0],sig_stokes[i][1],sig_stokes[i][2],sig_stokes[i][3],pval,sig_stokes[i][0]/pval);
		}
	      // Account for cross-coupling and differential gain and phase
	      /*	  createM_amp(data->diffGain[k],data->diffPhase[k],m_amp);
			  createM_cc(data->cc_eps1[k],data->cc_eps2[k],data->cc_phi1[k],data->cc_phi2[k],m_cc);
			  createM_feed(0,m_feed);
			  
			  multVectMat(m_feed,sig_stokes[i],tvect);
			  multVectMat(m_cc,tvect,tvect2);
			  multVectMat(m_amp,tvect2,sig_stokes[i]);
	      */
	      convertStokesAABBCRCI(sig_stokes[i],aabbcrci);
	      
	      if (l==0){
		data->aa[k*nbin+i] = aabbcrci[0];
		data->bb[k*nbin+i] = aabbcrci[1];
		data->cr[k*nbin+i] = aabbcrci[2];
		data->ci[k*nbin+i] = aabbcrci[3];
	      }
	      else
		{
		data->aa[k*nbin+i] += aabbcrci[0];
		data->bb[k*nbin+i] += aabbcrci[1];
		data->cr[k*nbin+i] += aabbcrci[2];
		data->ci[k*nbin+i] += aabbcrci[3];
		}
	    }  // Repeat for next bin

	}
      //	  exit(1);

    } // Repeat for next frequency channel

  //  printf("nsamp = %d\n",nsamp);

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

// Calculate parallactic angle in degrees
double calculatePara(dataStruct *data)
{
  long double mjd,last,ha,latitude,dec,ra;
  double pa;

  // MUST FIX - HARDCODED TO PARKES
  latitude = -32.0-59.0/60.0-54.263/60.0/60.0; // Degrees 

  ra = data->raj/M_PI*12; // Hours
  dec = data->decj/M_PI*180; // Degrees
 
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
	      else if (strcasecmp(word1,"RAJ")==0)
		{
		  char rajStr[128];
		  double h,m,s;
		  strcpy(rajStr,word2);
		  sscanf(rajStr,"%lf:%lf:%lf",&h,&m,&s);
		  strcpy(data->rajStr,rajStr);
		  data->raj = (h+m/60.0+s/60.0/60.0)*M_PI/12.0;
		}
	      else if (strcasecmp(word1,"DECJ")==0)
		{
		  char decStr[128];
		  double d,m,s;
		  int sign=1;

		  strcpy(decStr,word2);
		  strcpy(data->decjStr,decStr);
		  if (word2[0]=='-') sign=-1;
		  sscanf(decStr,"%lf:%lf:%lf",&d,&m,&s);
		  d = fabs(d);
		  data->decj = sign*(d+m/60.0+s/60.0/60.0)*M_PI/180.0;
		}
	    }
	}
    }
  fclose(fin);
}

void writeNewSub(fitsfile *fptr,dataStruct *data,int subint)
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
    int i,j;

    rajd = data->raj*180.0/M_PI;
    decjd = data->decj*180.0/M_PI;

    // Now calculate the Local Sidereal Time (LST)
    mjd = data->stt_imjd + (data->stt_smjd + data->stt_offs)/86400.0L;
    lst_sub = (double)calcLocalSiderealTime(mjd,data)*60.0*60.0;
    fits_get_colnum(fptr,CASEINSEN,"LST_SUB",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&lst_sub,&status);

    fits_get_colnum(fptr,CASEINSEN,"TSUBINT",&colnum,&status);
    fits_write_col(fptr,TDOUBLE,colnum,subint,1,1,&tsub,&status);


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
	  //	  printf("Have %g %g %g %g %g %g %g %g\n",scaleAA,scaleBB,scaleCR,scaleCI,offsAA,offsBB,offsCR,offsCI);
	  for (i=0;i<nbin;i++)
	    {
	      /*	      data->aa[j*nbin+i] = (data->aa[j*nbin+i]-aa_mean/(double)nbin)*scaleAA;
	      data->bb[j*nbin+i] = (data->bb[j*nbin+i]-bb_mean/(double)nbin)*scaleBB;
	      data->cr[j*nbin+i] = (data->cr[j*nbin+i]-cr_mean/(double)nbin)*scaleCR;
	      data->ci[j*nbin+i] = (data->ci[j*nbin+i]-ci_mean/(double)nbin)*scaleCI;*/

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
	      dat_scl[j]   = scaleAA; dat_offs[j]   = scaleAA*aa_mean/(double)nbin;
	      dat_scl[j+nchan] = scaleBB; dat_offs[j+nchan] = scaleBB*bb_mean/(double)nbin;
	      dat_scl[j+2*nchan] = scaleCR; dat_offs[j+2*nchan] = scaleCR*cr_mean/(double)nbin;
	      dat_scl[j+3*nchan] = scaleCI; dat_offs[j+3*nchan] = scaleCI*ci_mean/(double)nbin;
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
      //     fits_write_col(fptr,TFLOAT,colnum,1,1,nbin*nchan*4,vals,&status);
            fits_write_col(fptr,TFLOAT,colnum,1,1,nbin*nchan,data->aa,&status);
	    fits_write_col(fptr,TFLOAT,colnum,1,nbin*nchan+1,nbin*nchan,data->bb,&status);
	    fits_write_col(fptr,TFLOAT,colnum,1,2*nbin*nchan+1,nbin*nchan,data->cr,&status);
	    fits_write_col(fptr,TFLOAT,colnum,1,3*nbin*nchan+1,nbin*nchan,data->ci,&status);
    }
    else
      fits_write_col(fptr,TFLOAT,colnum,1,1,nbin*nchan,data->aa,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan, &status); 
    fits_write_col(fptr,TFLOAT,colnum,1,1,nchan,dat_freq,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
  fits_modify_vector_len (fptr, colnum, nchan, &status); 
    fits_write_col(fptr,TFLOAT,colnum,1,1,nchan,dat_wts,&status);

    fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan*npol, &status); 
    fits_write_col(fptr,TFLOAT,colnum,1,1,nchan*npol,dat_offs,&status);
    fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
    fits_modify_vector_len (fptr, colnum, nchan*npol, &status); 
    fits_write_col(fptr,TFLOAT,colnum,1,1,nchan*npol,dat_scl,&status);


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

void runTempo2(dataStruct *data)
{
  char execString[1024];
  double seg_length = 48000;
  int nfreqcoeff = 2;
  int ntimecoeff = 12;  
  double freq1 = data->cFreq - fabs(data->obsBW); 
  double freq2 = data->cFreq + fabs(data->obsBW);
  long double mjd1 = data->stt_imjd + data->stt_smjd/86400.0L - 60*60/86400.0L; // MUST FIX
  long double mjd2 = data->stt_imjd + data->stt_smjd/86400.0L + 60*60/86400.0L; // MUST FIX

  sprintf(execString,"tempo2-dev -pred \"PKS %Lf %Lf %g %g %d %d %g\" -f %s",mjd1,mjd2,freq1,freq2,ntimecoeff,nfreqcoeff,seg_length,data->exact_ephemeris);
  printf("Running tempo2 to get predictor\n");
  system(execString);
  printf("Complete running tempo2\n");
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

void createM_pa(double pa,double m_pa[4][4])
{
  pa = pa*M_PI/180.0;
  m_pa[0][0] = 1.0; m_pa[0][1] = 0; m_pa[0][2] = 0; m_pa[0][3] = 0;
  m_pa[1][0] = 0.0; m_pa[1][1] = cos(2*pa); m_pa[1][2] = sin(2*pa); m_pa[1][3] = 0;
  m_pa[2][0] = 0.0; m_pa[2][1] = -sin(2*pa); m_pa[2][2] = cos(2*pa); m_pa[2][3] = 0;
  m_pa[3][0] = 0.0; m_pa[3][1] = 0; m_pa[3][2] = 0; m_pa[3][3] = 1;
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

void createM_amp(double gain,double phase,double m_amp[4][4])
{
  m_amp[0][0] = 1.0;  m_amp[0][1] = gain; m_amp[0][2] = 0.0;        m_amp[0][3] = 0.0;
  m_amp[1][0] = gain; m_amp[1][1] = 1;    m_amp[1][2] = 0.0;        m_amp[1][3] = 0.0;
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
