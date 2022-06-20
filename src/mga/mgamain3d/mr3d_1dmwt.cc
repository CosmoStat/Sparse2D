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
**    Date:  25/08/08
**    
**    File:  mr3d_rcur.cc
**
*******************************************************************************
**
**    DESCRIPTION  Meyer 1d wavelets program
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "MeyerWT1D.h"
#include "GetLongOptions.h"

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
bool Input_mr=false;
bool Output_mr=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
bool Noise_calib=false;
bool Normalize=false;
bool NormalizeInv=false;
float SigmaNoise=-1;
float NSigma=3;
bool force3sigma=false;
int Niter=1;
int NbrScale = 4;
filter_type FilterType = FT_HARD;

/***************************************/
 
/***************************************/
	
static void usage(char *argv[])
{
	// int i;
	fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-n NbrScale]\n");
	fprintf(OUTMAN, "         [-N NbrScale](alias)\n");
	manline();

	fprintf(OUTMAN, "         [-v] Verbose.\n");
	manline();

	fprintf(OUTMAN, "         [-F] Compute forward transform\n");
	fprintf(OUTMAN, "                input file must be a FITS file\n");
	manline();

	fprintf(OUTMAN, "         [-R] Compute reconstruction\n");
	fprintf(OUTMAN, "                input file must be a MR file\n");
	manline();

	fprintf(OUTMAN, "         [--MRin] Input transformed file\n");
	fprintf(OUTMAN, "                load the transform from a .mr file");
	manline();

	fprintf(OUTMAN, "         [--MRout] Output transformed file\n");
	fprintf(OUTMAN, "                enable to save the transform into a .mr file");
	manline();

	fprintf(OUTMAN, "         [--noMR] To prevent creating MR output\n");
	fprintf(OUTMAN, "                when there is no reconstruction either.");
	manline();

	fprintf(OUTMAN, "         [--stat] Statistic extraction\n");
	fprintf(OUTMAN, "                in file 'NameOut_stat.dat'");
	manline();

	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
	manline();

	fprintf(OUTMAN, "         [-g sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
	manline();

	fprintf(OUTMAN, "         [-s NSigma]\n");
	fprintf(OUTMAN, "             Thresholding at NSigma * noise_lvl\n");
	manline();

	fprintf(OUTMAN, "         [-f Filtering type]\n");
	fprintf(OUTMAN, "            (1): Hard thresolding or use --hard\n");
	fprintf(OUTMAN, "             2 : Soft thresolding or use --soft\n");
	fprintf(OUTMAN, "             3 : Wiener with local blocks or use --wiener\n");
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholding\n");
	manline();

	fprintf(OUTMAN, "         [--force3sigma] (not yet) Force thresholding at 3 sigma\n");
	fprintf(OUTMAN, "             and not more at finest scale\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;  

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:N:vFRScg:s:f:i:")) != -1) 
	{
		switch (c) 
		{
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{
					fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);
					exit(-1);
				}
				break;
			case 'N': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{
					fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);
					exit(-1);
				}
			case 'v': 
				Verbose = True;
				break;
			case 'F': 
				Compute_forward = true;
				break;
			case 'R': 
				Compute_recons = true;
				break;
			case 'S': 
				Noise_calib = true;
				break;
			case 'c': 
				Normalize = true;
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 's': 
				if (sscanf(OptArg,"%f",&NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'f': int plop;
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);} else FilterType = (filter_type) plop;
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
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

int main(int argc, char *argv[])
{
	double start,stop;
	start=clock();
	
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	map<string,string>::iterator it;

	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Normalize; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); cerr<<"plop"<<endl; ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--MRin");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Input_mr; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
	it = opts.find("--noise_calib");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Noise_calib; }
	it = opts.find("--Niter");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Niter; }
	it = opts.find("--hard");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_HARD; }
	it = opts.find("--soft");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_SOFT; }
	it = opts.find("--wiener");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_WIENER; }
	it = opts.find("--force3sigma");
	if(it!=opts.end()){ istringstream ss(it->second); force3sigma=true; }

// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	if(!Input_mr) Compute_forward=true;
	if(Input_mr && Compute_forward){cerr<<"Cannot input a MR file AND make the forward transform."<<endl;exit(1);}
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Input_mr) cout << "Input MR file"<<endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;
		
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale << endl;  
		if(Normalize) 		cout << " Normalizing the transform" << endl;
		
		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;
		if(SigmaNoise>0)	cout << " Noise level = " <<  SigmaNoise    << endl;
		if(SigmaNoise>0)	cout << " NSigma Threshold = " <<  NSigma    << endl;
		if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
		if(force3sigma)		cout << " Force '3sigma lvl' at all scales" << endl;
		cout << endl;  
	}

	fltarray Data;
	fltarray *TabBand;
	fltarray Recons;
	MEYER_WT1D *DataW = new MEYER_WT1D;
	char filename[64];
	
// Forward transform
	if(Compute_forward)
	{
	// Input data
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);

	// RCurvelet initialisation
		DataW->init(NbrScale, Data.nx());
		
	// Forward Transform
		DataW->transform(Data, TabBand);
	}
	else // We load the transformed cube
	{
		strcat(Name_Imag_In,".mr");
		DataW->read(Name_Imag_In,TabBand,&NormalizeInv);
		
		// Inverse normalisation if normalized
		if(NormalizeInv) DataW->normalize_self(TabBand,1);
	}
	
// Tests
	DataW->test(TabBand);

// Statistic tools
	if(Noise_calib)
	{
		// variance estimation
		DataW->noise_calibration(TabBand, Name_Imag_Out);
	}
	if(Extract_stat) DataW->extract_stat(TabBand, Name_Imag_Out);
	
// Thresholding
	if( SigmaNoise>=0 )
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
		{
			DataW->threshold(TabBand, SigmaNoise, NSigma, FilterType, force3sigma);
			for(int i=1;i<Niter;i++)
			{
				DataW->recons(TabBand, Data);
				DataW->transform(Data, TabBand);
				DataW->threshold(TabBand, SigmaNoise, NSigma, FilterType, force3sigma);
			}
		}
		else if(FilterType==FT_WIENER) DataW->wiener(TabBand, SigmaNoise, 3);
	}
	
// Save the transform
	if(Output_mr)
	{
	// Normalize the transform if not already
		if(Normalize) DataW->normalize_self(TabBand, 0);
	// Save the transform
		sprintf(filename,"%s.mr",Name_Imag_Out);
		DataW->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) DataW->normalize_self(TabBand, 1);
	}
	
	//DataW->values_at(TabBand,(char*)"coeflist.dat",Name_Imag_Out);
	
// Reconstruction and save output
	if(Compute_recons) 
	{
		//	for(int s=0;s<DataW->NbrScale3D;s++){sprintf(filename,"%s_B%d.fits",Name_Imag_Out,s);writefltarr(filename, TabBand[s]);}
		DataW->recons(TabBand,Recons);
		sprintf(filename,"%s_recons.fits",Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	
// Free memory	
	delete [] TabBand;
	delete DataW;
	
	for(int i=0;i<short_argc;i++)
		delete [] short_argv[i];
	delete [] short_argv;
	
	Data.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





