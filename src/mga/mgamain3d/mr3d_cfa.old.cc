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
**    File:  mr3d_cfa.cc
**
*******************************************************************************
**
**    DESCRIPTION  3D combined filtering
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "bfct.h"
#include "RCurvelet3D.h"
#include "BCurvelet3D.h"
#include "uowt.h"
#include "IM3D_DCT.h"
#include "dowt.h"

char Name_Imag_In[256]; /* input file image */
char Name_Mask[256]; /* input mask name */
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
float SigmaNoise=-1;
bool Use3sig=false;
bool force4sigma=false;
bool UseCubeSigma=false;

bool use_max=false;
bool use_min=false;
type_extraction extract_type = DEF_TYPE_EXTRACTION;
filter_type FilterType = FT_HARD;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
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
	
	fprintf(OUTMAN, "         [-f Filtering type]\n");
	fprintf(OUTMAN, "            (1): Hard thresolding or use --hard\n");
	fprintf(OUTMAN, "             2 : Soft thresolding or use --soft\n");
	fprintf(OUTMAN, "             3 : Wiener with local blocks or use --wiener\n");
	fprintf(OUTMAN, "             4 : FDR or use --fdr\n");
	fprintf(OUTMAN, "             5 : Stein Block Thresholding or use --stein\n");
	manline();

	fprintf(OUTMAN, "         [--force4sigma] Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at the finest scale\n");
	manline();

	fprintf(OUTMAN, "         [--use3sig] Use the cumulative at lvl 1-2.7e-3\n");
	fprintf(OUTMAN, "             to estimate the normalisation coefficients\n");
	fprintf(OUTMAN, "             instead of the variance\n");
	manline();

	fprintf(OUTMAN, "         [-P][--positive] Positivity constraint\n");
	manline();

	fprintf(OUTMAN, "         [--max-value] Maximum value of output data\n");
	fprintf(OUTMAN, "             255 for 8bit data.\n");
	manline();

	fprintf(OUTMAN, "         [--min-value] Minimum value of output data\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:d:vg:a:s:f:i:COB:P")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose = True;
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
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

void set_border( fltarray & Data, float min_value, float max_value)
{
	if(use_min)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) < min_value ) Data(i,j,k) = min_value;

	if(use_max)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) > max_value ) Data(i,j,k) = max_value;
}

int main(int argc, char *argv[])
{

// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	map<string,string>::iterator it;

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	
// Get short options
	filtinit(short_argc, short_argv);

// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << endl;  
	}

// **********************************************
// **********************************************
// 				  CFA Parameters
// **********************************************
// **********************************************

	bool GVerbose = true;

// BeamCurvelets				
	bool Use_BC			=true;	
	int Nscale_BC		=3;		
	int BlockSize_BC	=17;	
	Bool BlockOver_BC	=False;	
	bool BC_no_coarse	=false;	
	
// Levels
	float Nsig_orig=3.;
	
// **********************************************
// **********************************************
// 				CFA initialization
// **********************************************
// **********************************************
	
// Input and output data
	fltarray Data, Orig, FullOrig;
	fltarray Recons;

// Input data and mask
	char filename[64];
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	
// Filter type for BCurvelets and Wavelets
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
//	SB1D->Border = I_CONT;
	
// 3D transforms declaration
	BCurvelet3D *DataBC = new BCurvelet3D;
	fltarray *TabBand_BC,	*TabBand_BC_orig;

// BC Initialization
	type_linelet3d_WTtrans CurTrans = DEF_LIN3D_TRANS;
	type_wavelet3D TypeW3D = DEF_TYPE_W3D;
	DataBC->CurTrans=CurTrans;
	DataBC->BlockOverlap=false;
	DataBC->set_stat(false);
	DataBC->set_3sigma(Use3sig);
	DataBC->set_no_coarse(false);
	DataBC->init(Data.nx(), Data.ny(), Data.nz(), 3, 17, TypeW3D, SB1D);

// Forward Transform
	DataBC->transform(Data, TabBand_BC_orig, true);

// Filtering
	DataBC->threshold(TabBand_BC_orig, SigmaNoise, Nsig_orig, FT_HARD, force4sigma, UseCubeSigma);

// Reconstruction
	DataBC->recons(TabBand_BC_orig, Recons);

// Save output
	sprintf(filename,"%s_bcur_orig_recons.fits", Name_Imag_Out);
	writefltarr(filename, Recons);

	exit(0);
}





