/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  10/09/08
**    
**    File:  mr3d_dct.cc
**
*******************************************************************************
**
**    DESCRIPTION  3D Fast Curvelet program
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "writefits3d.h"
#include "fftw3.h"

void cosft2(float y[], int n, int isign);

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

/********************************/
//		Input variables			//
/********************************/

Bool Verbose=False;
bool Compute_forward=false;
bool Compute_recons=false;
float SigmaNoise=-1;
float NSigma=3;
int filter=1;
/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
	manline();

	fprintf(OUTMAN, "         [-F] Compute forward transform\n");
	fprintf(OUTMAN, "                input file must be a FITS file\n");
	manline();

	fprintf(OUTMAN, "         [-R] Compute reconstruction\n");
	fprintf(OUTMAN, "                input file must be a MR file\n");
	manline();

	fprintf(OUTMAN, "         [-g sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
	manline();
	
	fprintf(OUTMAN, "         [-s NSigma]\n");
	fprintf(OUTMAN, "             Thresholding at NSigma * noise_lvl\n");
	manline();
	
	fprintf(OUTMAN, "         [-f Filtering type]\n");
	fprintf(OUTMAN, "            (1): Hard thresolding\n");
	fprintf(OUTMAN, "             2 : Soft thresolding\n");
	manline();

	exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"vFRSB:Og:a:s:f:")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose = True;
				break;
			case 'F': 
				Compute_forward = true;
				break;
			case 'R': 
				Compute_recons = true;
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 's': 
				if (sscanf(OptArg,"%f",&NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'f': 
				if (sscanf(OptArg,"%d",&filter) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case '?': 
				usage(argv);
				break;
			default : 
				usage(argv);
				break;
		}
	} 


	/* get optional input file names from trailing 
	parameters and open files */
	if (OptInd < argc)
		strcpy(Name_Imag_In, argv[OptInd++]);
	else
		usage(argv);

	if (OptInd < argc)
		strcpy(Name_Imag_Out, argv[OptInd++]);
	else
		usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}


/***************/

void sig_dct ( fltarray& Sig, fltarray& Trans, Bool Reverse) {

    int Nx = Sig.nx();
    float Norm = sqrt(2./float(Nx));
    int isign = (Reverse == False) ? 1: -1;
    float *Ptr = Trans.buffer() - 1;
    Trans = Sig;
    cosft2 (Ptr, Nx, isign);
    
    for (int i=0; i<Nx; i++) Trans(i) *= Norm;     
}

void dct(fltarray &In, fltarray &Out, bool backward)
{
//	if(Verbose) cerr<<"dctw2d..."<<endl;
	double norm;
	int nx = In.nx();
	double* in = new double[nx];
	double* out = new double[nx];
	Out.alloc(nx);
	
	fftw_plan plan;

// dct on x axis
	norm = (double)sqrt(double(nx)*2.);
	
	for (int i = 0; i < nx; i++ )
		in[i] = In(i);
	
// Transform
	if(backward)
		plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT01, FFTW_ESTIMATE);
	else
		plan = fftw_plan_r2r_1d(nx, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
	
	fftw_execute ( plan );

	// Normalization
	for (int i = 0; i < nx; i++ )
		Out(i) = out[i]/norm;
}

int main(int argc, char *argv[])
{
	double start,stop;
	start=clock();
	
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	map<string,string>::iterator it;

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }

// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !Compute_forward) Compute_forward=true;
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   

		// Filtering
		if(SigmaNoise>0)
		{
			cout << "Filtering type : " << filter<<endl;
			cout << " Noise level = " <<  SigmaNoise    << endl;
			cout << " NSigma Threshold = " <<  NSigma    << endl;
		}
		cout << endl;  
	}

	char filename[64];
	fltarray Data;
	fltarray TabBand;
	fltarray Recons;
	int Nx;;

	if(Compute_forward)
	{
	// Input data
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);
		Nx = Data.nx();
		if(Nx != pow(2,log2(Nx)))
		{
			cerr<<"This works only for 2^n signal size."<<endl;
			cerr<<Nx<<","<<pow(2,log2(Nx))<<","<<log2(Nx)<<","<<endl;
			exit(0);
		}
		TabBand.alloc(Nx);
		dct(Data, TabBand, false);
	}
	else
	{
		strcat(Name_Imag_In,".mr");
		fits_read_fltarr(Name_Imag_In,TabBand);
		Nx = TabBand.nx();
	}
	
	if( SigmaNoise>=0)
	{
		if(filter==1)
			for(int i=0;i<Nx;i++)
				TabBand(i) *= abs(TabBand(i)) > SigmaNoise*NSigma;
		else //soft
			for(int i=0;i<Nx;i++)
				TabBand(i) = ::soft_threshold(TabBand(i),SigmaNoise*NSigma);
	}
	
	if(Compute_recons) 
	{
		Recons.alloc(Nx);
		dct(TabBand, Recons, true);
		sprintf(filename,"%s_recons.fits",Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	else
	{
		sprintf(filename,"%s.mr",Name_Imag_Out);
		writefltarr(filename, TabBand);
	}

// Free all arrays
	Data.~fltarray();
	Recons.~fltarray();
	
	exit(0);
}





