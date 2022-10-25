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
**    File:  mr3d_fct.cc
**
*******************************************************************************
**
**    DESCRIPTION  3D Fast Curvelet program
**    ----------- 
**                 
******************************************************************************/
//get_sigma_mad inside defmath.h
#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "IM_Obj.h"
#include "Usage.h"
#include "PCur.h"
#include "FCur.h"
#include "writefits3d.h"
#include "MGA_Inc.h"
#include "FCur_Frame.h"


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
float Alpha=0.01;
float NSigma=3;
bool force3sigma=false;
//	bool UseCubeSigma=false;
int Niter=1;
int NbrScale = 3;
int NbrDir2d = 16;
Bool Isotropic=False;
bool Complex_transform=false;
filter_type FilterType = FT_HARD;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n NbrScale]\n");
    manline();
    
	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
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
	
	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
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
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholding\n");
	manline();
	
	fprintf(OUTMAN, "         [--force3sigma] Force thresholding at 3 sigma\n");
	fprintf(OUTMAN, "             and not more at finest scale\n");
	manline();

	fprintf(OUTMAN, "         [-I]\n");
	fprintf(OUTMAN, "             Isotropic WT (not yet).\n");
	manline();

	fprintf(OUTMAN, "         [-C]\n");
	fprintf(OUTMAN, "             To use the Complex transform.\n");
	fprintf(OUTMAN, "             (only transform and reconstruction available yet).\n");
	manline();

	fprintf(OUTMAN, "         [-d Ndir2d]\n");
	fprintf(OUTMAN, "             number of directions (in 2d), min=12, def=16.\n");
	manline();
	
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"In:d:vFRScg:a:s:f:i:Cx:")) != -1) 
	{
		switch (c) 
		{
			case 'I': 
				Isotropic = (Isotropic == True) ? False: True;
				break;
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);exit(-1); }
				break;
			case 'd': 
				if (sscanf(OptArg,"%d",&NbrDir2d) != 1) 
				{ fprintf(OUTMAN, "Error: bad number of directions: %s\n", OptArg); exit(-1); }
				break;
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
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'C': 
				Complex_transform = true;
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
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	map<string,string>::iterator it;

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--MRin");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Input_mr; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--hard");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_HARD; }
	it = opts.find("--soft");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_SOFT; }
	it = opts.find("--wiener");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_WIENER; }
	it = opts.find("--Niter");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Niter; }
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
		if(Complex_transform) cout << "Transform in Complex space" << endl;   
		cout << "NbrScale = " <<  NbrScale    << endl;  
		cout << "NbrDir2d = " <<  NbrDir2d    << endl;  

		if(Extract_stat)	cout << " Statistics output" << endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
			if(force3sigma)		cout << " Force '3sigma lvl' at all scales" << endl;
		}
		cout << endl;  
	}

	fltarray Data;
	char filename[64];
	
	FCur_Frame CurF;

	if(Compute_forward)
	{
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);
		
		CurF.alloc_from_coarse(NbrScale, Data.nx(), Data.ny(), Data.nz(), NbrDir2d, False, False, (Bool)!Complex_transform);
		if(SigmaNoise>=0) CurF.get_norm_coeff(NSigma);
		CurF.transform(Data, true);
	}
	else // We load the transformed cube
	{
		//strcat(Name_Imag_In,".mr");
		//Cur.read(Name_Imag_In);
	}
	
// Thresholding
	if( SigmaNoise>=0 )
		CurF.threshold(SigmaNoise, NSigma, FilterType);
	
// Save the transform
	if(Output_mr & !Complex_transform)
	{
		//strcat(Name_Imag_Out,".fits");
		//Cur.write(Name_Imag_Out);
	}

// Reconstruction
	if (Compute_recons)
	{
		CurF.recons(Data);
		
		sprintf(filename,"%s_recons.fits",Name_Imag_Out);
		writefltarr(filename, Data);
	}


// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
	Data.~fltarray();
	
	exit(0);
}





