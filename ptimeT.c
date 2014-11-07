// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptimeT.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "simulatePseudoBB.h"

int main (int argc, char *argv[])
{
	int h,i,j,k;

	/*
	if (argc != 8)
	{
		printf ("Usage: ptime_time -f fname -std tname (-pt tname) -o oname -single (-multi)\n"
	            "Derive the TOAs\n"
	            "fname: data file; tname: templates; oname: output .tim; -std: standard template format; -pt: ptime template;\n"
				"-single: do freq-dependent matching and get one TOA; -multi: do freq-dependent matching and get TOAs for each channel.\n");
	    exit (0);
	}
	*/


	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	char tname[128];   // name of template
	char oname[128];   // name of output .tim
	int mode; // to distinguish different type of templates
	int tmode; // to distinguish different type of algorithm

	int index, n;
	for (i=0;i<argc;i++)
    {
		if (strcmp(argv[i],"-f") == 0)
		{
            index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-single") != 0)
			{
				n++;
		    }
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 0; // standard template format
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
		else if (strcmp(argv[i],"-single")==0)
		{
			tmode = 0; // do freq-dependent matching, and get one TOA
		}
		else if (strcmp(argv[i],"-multi")==0)
		{
			tmode = 1; // do freq-dependent matching, and get TOA for each channel
		}
    }

	// name of different extension of data files
	char name_data[50]; 
	char name_predict[50]; 
	char name_psrparam[50]; 

	char data[] = "[SUBINT]";
	char predict[] = "[T2PREDICT]";
	char psrparam[] = "[PSRPARAM]";

	// read a std
	char std[50];
	strcpy(std,tname);
	if ( mode == 0)
	{
		strcat(std, data);
	}

	/////////////////////////////////////////////////////////////////////////////////
	// open file to write toa 
	FILE *fp;
	if ((fp = fopen(oname, "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	
	fprintf (fp, "FORMAT 1\n");

	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		strcpy(fname,argv[k]);
		printf ("%s\n", fname);

		// name of different extension
		strcpy(name_data,fname);
		strcpy(name_predict,fname);
		strcpy(name_psrparam,fname);

		strcat(name_data, data);
		strcat(name_predict, predict);
		strcat(name_psrparam, psrparam);

		////////////////////////////////////////////////////
	
		double psrfreq;
		psrfreq = read_psrfreq(name_psrparam);
		printf ("psrfreq: %.15lf\n", psrfreq);
	
		double dm;
		dm = readDm(name_data);
		printf ("DM: %.4lf\n", dm);
	
		////////////////////////////////////////////////
		long int imjd, smjd;
		double offs;
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		imjd = stt_imjd(fname);
		smjd = stt_smjd(fname);
		offs = stt_offs(fname);

		nchn = get_nchan(name_data);	
		npol = get_npol(name_data);	
		nsub = get_subint(name_data);	
		nphase = get_nphase(name_data);	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		// read a std
		double s_multi[nphase*nchn*npol];
		double s_temp[nphase];

		//readfile(argv[1],&n,tt,s);
		//read_prof(std,1,s_multi,nphase);
	
		// check the channel number of template
		check_std(std,1,mode,nchn,nphase);

		read_std(std,1,s_multi,nphase,mode,nchn);

		////////////////////////////////////////////////////////////////////////////////

		double p_multi[nchn*npol*nphase];
		double p_temp[nphase];

		double rms[nchn];  // rms for each profile
		double b[nchn];  // b for each profile
		double phase, e_phase;
		long double e_dt;  
		long double t;     // TOA
		double frequency;

		// start to derive toa from different subint
		for (h = 1; h <= nsub; h++)
		{
			// simulate data

			//SNR = 500.0 + 200.0*i;
			//simulate(n,SNR,s,p_temp);

			// read profiles from data file
			read_prof(name_data,h,p_multi,nphase);
			//readfile(argv[2],&n,tt,p_multi);

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				for (j = 0; j < nphase; j++)
				{
					//printf ("%lf %lf\n", p_multi[j], s[j]);
					//s_multi[i*nphase + j] = s[j];
					p_temp[j] = p_multi[i*nphase + j];
					s_temp[j] = s_multi[i*nphase + j];
				}

				// calculate toa, rms for each channel
				get_toa (s_temp, p_temp, &phase, &e_phase, psrfreq, nphase, &rms[i], &b[i]);

				// if tmode == 1, get TOA for each channel, and transform phase shifts to MJD TOAs
				if ( tmode == 1)
				{
					form_toa (name_data, name_predict, h, i, nchn, imjd, smjd, offs, phase, e_phase, &t, &e_dt, &frequency);
					fprintf (fp, "%s  %lf  %.15Lf  %Lf  7 -f c%d\n", fname, frequency, t, e_dt*1e+6, i+1);
				}
			}

			// if tmode == 0, do freq-dependent template matching, get one phase shift
			if ( tmode == 0)
			{
				get_toa_multi (name_data, name_predict, h, s_multi, p_multi, rms, nchn, &phase, &e_phase, psrfreq, nphase, &frequency);
				//get_toa_multi(s_multi, p_multi, rms, b, nchn, &phase, &e_phase, psrfreq, nphase);
				getToaMultiDM (name_data, name_predict, h, s_multi, p_multi, rms, nchn, &phase, &e_phase, psrfreq, nphase, dm, &frequency);

				// transform phase shifts to MJD TOAs
				form_toa_multi (name_data, name_predict, h, nchn, imjd, smjd, offs, phase, e_phase, &t, &e_dt, frequency);

				fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);
			}
		}
	}

    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");

	return 0;
}
