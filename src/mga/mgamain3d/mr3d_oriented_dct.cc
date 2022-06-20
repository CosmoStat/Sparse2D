
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
**    Date:  
**    
**    File:  
**
*******************************************************************************
**
**    DESCRIPTION  see Oriented_DCT3D.h
**    ----------- 
**                 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "Oriented_DCT3D.h"

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
bool MRidl=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
bool Noise_calib=false;
bool Normalize=false;
bool NormalizeInv=false;
float SigmaNoise=-1;
float NSigma=3;
bool force4sigma=false;
int skip_order=-1;
bool threshold_coarse=false;
bool no_coarse=false;
bool no_fine=false;
int BlockSize=0;
Bool BlockOverlap=False;

filter_type FilterType = FT_HARD;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n NbrScale]\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_3DCUR_NBR_SCALE);
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
	
	fprintf(OUTMAN, "         [--MRidl] Output formated for idl\n");
	fprintf(OUTMAN, "                no header, all the bands in one vector");
	manline();
	
	fprintf(OUTMAN, "         [--noMR] To prevent creating MR output\n");
	fprintf(OUTMAN, "                when there is no reconstruction either.");
	manline();
	
	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] Set the coarse scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [--no-fine] Set the fine scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [-P][--positive] Positivity constraint\n");
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

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholding\n");
	manline();
	
	fprintf(OUTMAN, "         [--force4sigma] Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at the finest scale\n");
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
	
	fprintf(OUTMAN, "         [-x wedge extraction type]\n");
	fprintf(OUTMAN, "             0 : Back : when puting back wedges to scales\n");
	fprintf(OUTMAN, "             1 : Forward : when extracting wedges from scales\n");
	fprintf(OUTMAN, "            (2): Forward & Back : half forward and back\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"In:d:vFRScg:a:s:f:i:Cx:OB:P")) != -1) 
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
			case 'O': 
				BlockOverlap = True;
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
			case 'f': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);} else FilterType = (filter_type) plop;
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

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--BS");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>BlockSize; }
	it = opts.find("--MRin");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Input_mr; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--MRidl");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>MRidl; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
	it = opts.find("--noise_calib");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Noise_calib; }
	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Normalize; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--hard");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_HARD; }
	it = opts.find("--soft");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_WIENER; }
	it = opts.find("--force4sigma");
	if(it!=opts.end()){ istringstream ss(it->second); force4sigma=true; }
	it = opts.find("--no-fine");
	if(it!=opts.end()){ istringstream ss(it->second); no_fine=true; }
	it = opts.find("--no-coarse");
	if(it!=opts.end()){ istringstream ss(it->second); no_coarse=true; }
	it = opts.find("--threshold-coarse");
	if(it!=opts.end()){ istringstream ss(it->second); threshold_coarse=true; }
	it = opts.find("--skip-order");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>skip_order; }
	
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

		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(FilterType) <<endl;
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
		}
		cout << endl;  
	}

	fltarray Data;
	fltarray Trans;
	fltarray Recons;
	Oriented_DCT3D *DataC = new Oriented_DCT3D;
	char filename[64];
	
	bool allocTB=true;
		
	if(Compute_forward)
	{
	// Input data
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);

	// FCurvelet initialisation
		DataC->init(Data.nx(),Data.ny(),Data.nz(),BlockSize,BlockOverlap);

	// Forward Transform
		DataC->transform(Data, Trans);
		allocTB=false;
	}
/*	else // We load the transformed cube
	{
	// read data
		strcat(Name_Imag_In,".mr");
		if(MRidl) DataC->read_nohead(Name_Imag_In,TabBand,&NormalizeInv);
		else DataC->read(Name_Imag_In,TabBand,&NormalizeInv);
		
	// Noise calibration if not already done, and if useful here
		if(!DataC->isset_tabsigma())
			if( SigmaNoise>=0 || Extract_stat || Normalize || NormalizeInv)
			{
				DataC->get_norm_coeff(3, TabBand, Name_Imag_Out, allocTB);
				allocTB=false;
			}
			
	// Inverse normalisation if normalized
		if(NormalizeInv) DataC->normalize_self(TabBand,1);
	}
	
*/
//DataC->temp(Trans);
/*
// Statistic tools
	if(Noise_calib & !Complex_transform)
	{
		// variance estimation
		DataC->noise_calibration(TabBand, Name_Imag_Out);
	}
	if(Extract_stat) DataC->extract_stat(TabBand, Name_Imag_Out);
*/		
// Thresholding
	if( SigmaNoise>=0)
	{
		//if(no_coarse) DataC->set_no_coarse(true);
		//if(no_fine) DataC->set_no_fine(true);
		DataC->set_skip_order(skip_order);
		DataC->set_threshold_coarse(threshold_coarse);
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
			DataC->threshold(Trans, SigmaNoise, NSigma, FilterType);
		//else if(FilterType==FT_WIENER) DataC->wiener(TabBand, SigmaNoise, 3);
		//else if(FilterType==FT_FDR) DataC->fdr(TabBand, Alpha, SigmaNoise);
		//else if(FilterType==FT_SBT) DataC->stein_block_threshold(TabBand, SigmaNoise);
		else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";
	}

// Save the transform
	if(Output_mr)
	{
	// Normalize the transform if not already
		if(Normalize) DataC->normalize_self(Trans, 0);
	// Save the transform
		sprintf(filename,"%s.mr",Name_Imag_Out);
		DataC->write(filename, Trans, Normalize);
	// Inverse normalisation if normalized and the program is going further
		if(Normalize & Compute_recons) DataC->normalize_self(Trans, 1);
	}


// Reconstruction and save output
	if(Compute_recons) 
	{
		//if(positivity) DataC->set_positivity(true);
		DataC->recons(Trans, Recons);
		sprintf(filename,"%s_recons.fits",Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	
// Free all arrays
	delete DataC;
	Trans.~fltarray();
	Data.~fltarray();
	Recons.~fltarray();
	
	exit(0);
}





