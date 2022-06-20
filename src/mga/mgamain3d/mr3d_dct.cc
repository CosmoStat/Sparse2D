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
#include "MGA_Inc.h"

#include "IM3D_DCT.h"

char* Name_Imag_In; /* input file image */
char* Name_Imag_Out; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

/********************************/
//		Input variables			//
/********************************/

Bool Verbose=False;
bool Output_mr=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
bool Noise_calib=false;
float SigmaNoise=-1;
float Alpha=0.01;
float NSigma=3;
bool force4sigma=false;
bool Normalize=false;
bool NormalizeInv=false;
int BlockSize=8;
Bool BlockOverlap=False;
bool lapped=false;
int skip_order=-1;
filter_type FilterType = FT_HARD;

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

	fprintf(OUTMAN, "         [--MRout] Output transformed file\n");
	fprintf(OUTMAN, "                enable to save the transform into a .mr file");
	manline();
	
	fprintf(OUTMAN, "         [--noMR] To prevent creating MR output\n");
	fprintf(OUTMAN, "                when there is no reconstruction either.");
	manline();
	
	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
	manline();

	fprintf(OUTMAN, "         [--overlap][-O] BlockOverlap\n");
	manline();

	fprintf(OUTMAN, "         [--lapped] Lapped transform\n");
	manline();

	fprintf(OUTMAN, "         [--BS Blocksize][-B BlockSize]\n");
	fprintf(OUTMAN, "             must be a power of 2\n");
	fprintf(OUTMAN, "             default = 8\n");
	manline();

	fprintf(OUTMAN, "         [-g sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
	manline();
	
	fprintf(OUTMAN, "         [-a alpha]\n");
	fprintf(OUTMAN, "             fdr parameter\n");
	manline();
	
	fprintf(OUTMAN, "         [-s NSigma]\n");
	fprintf(OUTMAN, "             Thresholding at NSigma * noise_lvl\n");
	manline();
	
	fprintf(OUTMAN, "         [--skip_order Order]\n");
	fprintf(OUTMAN, "             Set the coefficients of order <= order to 0\n");
	manline();
	
	fprintf(OUTMAN, "         [-f Filtering type]\n");
	fprintf(OUTMAN, "            (1): Hard thresolding or use --hard\n");
	fprintf(OUTMAN, "             2 : Soft thresolding or use --soft\n");
	fprintf(OUTMAN, "             3 : Wiener with local blocks or use --wiener\n");
	fprintf(OUTMAN, "             4 : FDR or use --fdr\n");
	fprintf(OUTMAN, "             5 : Stein Block Thresholding or use --stein\n");
	manline();

	exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

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
			case 'S': 
				Noise_calib = true;
				break;
			case 'B': 
				if (sscanf(OptArg,"%d",&BlockSize) != 1) 
				{fprintf(OUTMAN, "Error: bad Block size: %s\n", OptArg);exit(-1);}
				break;
			case 'O': 
				BlockOverlap = True;
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'a': 
				if (sscanf(OptArg,"%f",&Alpha) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 's': 
				if (sscanf(OptArg,"%f",&NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'f': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);} else FilterType = (filter_type) plop;
				break;
			case 'c': 
				Normalize = true;
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
		Name_Imag_In = strdup(argv[OptInd++]);
	else
		usage(argv);

	if (OptInd < argc)
		Name_Imag_Out = strdup(argv[OptInd++]);
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

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
	it = opts.find("--BS");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>BlockSize; }
	it = opts.find("--overlap");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	BlockOverlap=(Bool)plop; }
	it = opts.find("--lapped");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>lapped; }
	it = opts.find("--noise_calib");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Noise_calib; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--skip_order");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>skip_order; }
	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); Normalize=true; }

// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   

		cout << "BlockSize = " <<  BlockSize    << endl;
		if(Normalize) cout << "Normalizing the transform" << endl;
		if(BlockOverlap) cout << "BlockOverlaping activated" << endl;
		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(FilterType) <<endl;
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
			if(skip_order>-1)	cout << " skip_order = " <<  skip_order    << endl;
		}
		cout << endl;  
	}

	char* filename;
	fltarray Data;
	fltarray TabBand;
	fltarray Recons;
	IM3D_DCT *dct = new IM3D_DCT;
	
	if(Compute_forward)
	{
	// Input data
		Name_Imag_In = add_fits(Name_Imag_In);
		fits_read_fltarr(Name_Imag_In, Data);
		
//		dct->set_lapped(lapped);
//		dct->set_skip_order(skip_order);
		dct->init(Data.nx(),Data.ny(),Data.nz(),BlockSize,BlockOverlap);
		dct->transform(Data,TabBand);
	}
	else
	{
		Name_Imag_In = add_mr(Name_Imag_In);
		dct->read(Name_Imag_In,TabBand,&NormalizeInv);
	
	// Inverse normalisation if normalized
		if(NormalizeInv) dct->normalize_self(TabBand,1);
	}
	
// Statistic tools
	if(Noise_calib)
	{
		// variance estimation
		//dct->noise_calibration(TabBand, Name_Imag_Out);
	}
	//if(Extract_stat) dct->extract_stat(TabBand, Name_Imag_Out);
	
// Thresholding
	if( SigmaNoise>=0)
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
			dct->threshold(TabBand, SigmaNoise, NSigma, FilterType);
		else if(FilterType==FT_WIENER) dct->wiener(TabBand, SigmaNoise, 3);
		else if(FilterType==FT_FDR) dct->fdr(TabBand, Alpha, SigmaNoise);
		else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";
	}
	
// Save the transform
	if(Output_mr)
	{
	// Normalize the transform if needed
		if(Normalize) dct->normalize_self(TabBand, 0);
	// Save the transform
		filename = add_mr(Name_Imag_Out);
		dct->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) dct->normalize_self(TabBand, 1);
		
	}

// Reconstruction and save output
	if(Compute_recons) 
	{
		dct->recons(TabBand,Recons);
		filename = add_fits(Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	
// Free all arrays
	Data.~fltarray();
	Recons.~fltarray();
	delete dct;
	
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





