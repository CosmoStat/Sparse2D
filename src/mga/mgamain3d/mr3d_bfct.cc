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
#include "bfct.h"

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

// Curvelet 3D Parameters
	BFCurvelet3D_params P;

Bool Verbose=False;
bool Output_mr=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Noise_calib=false;
int Niter=1;
bool Extract_stat=false;

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
	fprintf(OUTMAN, "             5 : Stein Block Thresholding or use --stein\n");
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholding\n");
	manline();
	
	fprintf(OUTMAN, "         [--force4sigma] Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at the finest scale\n");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] Set the coarse scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [--no-fine] Set the fine scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [-P][--positive] Positivity constraint\n");
	manline();

	fprintf(OUTMAN, "         [--max-value] Maximum value of output data\n");
	fprintf(OUTMAN, "             255 for 8bit data.\n");
	manline();

	fprintf(OUTMAN, "         [--min-value] Minimum value of output data\n");
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

	fprintf(OUTMAN, "         [--overlap][-O] BlockOverlap\n");
	manline();

	fprintf(OUTMAN, "         [--lapped][-L] Lapped transform, without extra redundancy\n");
	manline();

	fprintf(OUTMAN, "         [--lowmem] Low memory usage\n");
	fprintf(OUTMAN, "             Doesn't evaluate the full transform,\n");
	fprintf(OUTMAN, "             use for simple filtering (block or orientation independent)\n");
	manline();

	fprintf(OUTMAN, "         [--reuse_input] Erase the input data (from ram)\n");
	fprintf(OUTMAN, "             to put the filtered version instead\n");
	manline();

	fprintf(OUTMAN, "         [--BS Blocksize][-B BlockSize]\n");
	fprintf(OUTMAN, "             17 or 33\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_3DCUR_BLOCK_SIZE);
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:d:vFRScg:a:s:f:i:Cx:OLB:P")) != -1) 
	{
		switch (c) 
		{
			case 'n': 
				if (sscanf(OptArg,"%d",&P.NbrScale) != 1) 
				{fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);exit(-1); }
				break;
			case 'd': 
				if (sscanf(OptArg,"%d",&P.NbrDir2d) != 1) 
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
			case 'g': 
				if (sscanf(OptArg,"%f",&P.SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'a': 
				if (sscanf(OptArg,"%f",&P.Alpha) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 's': 
				if (sscanf(OptArg,"%f",&P.NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'f': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);} else P.FilterType = (filter_type) plop;
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'C': 
				P.RealData = False;
				break;
			case 'x': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: extraction type: %s\n", OptArg);exit(-1);} else P.extract_type = (type_extraction) plop;
				break;
			case 'O': 
				P.BlockOverlap = True;
				break;
			case 'L': 
				P.lapped = true;
				break;
			case 'B': 
				if (sscanf(OptArg,"%d",&P.BlockSize) != 1) 
				{fprintf(OUTMAN, "Error: bad Block size: %s\n", OptArg);exit(-1);}
				break;
			case 'P': 
				P.Positivity = true;
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
	bool tmp;
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"MRout",Output_mr);
	GetLongOption(opts,"noMR",No_mr);
	GetLongOption(opts,"noise_calib",Noise_calib);
	GetLongOption(opts,"stat",Extract_stat);
	GetLongOption(opts,"BS",P.BlockSize);
	GetLongOption(opts,"overlap",P.BlockOverlap);
	GetLongOption(opts,"lapped",P.lapped);
	GetLongOption(opts,"SigmaNoise",P.SigmaNoise);
	GetLongOption(opts,"NSigma",P.NSigma);
	GetLongOption(opts,"Niter",Niter);
	GetLongOption(opts,"hard",tmp)	&& (P.FilterType = FT_HARD);
	GetLongOption(opts,"soft",tmp)	&& (P.FilterType = FT_SOFT);
	GetLongOption(opts,"wiener",tmp)&& (P.FilterType = FT_WIENER);
	GetLongOption(opts,"fdr",tmp) 	&& (P.FilterType = FT_FDR);
	GetLongOption(opts,"stein",tmp)	&& (P.FilterType = FT_SBT);
	GetLongOption(opts,"force4sigma",P.force4sigma);
	GetLongOption(opts,"no-fine",P.no_fine);
	GetLongOption(opts,"no-coarse",P.no_coarse);
	GetLongOption(opts,"threshold-coarse",P.threshold_coarse);
	GetLongOption(opts,"positive",P.Positivity);
	GetLongOption(opts,"min-value",P.min_value) && (P.use_min=true);
	GetLongOption(opts,"max-value",P.max_value) && (P.use_max=true);
	GetLongOption(opts,"lowmem",P.lowmem);
	GetLongOption(opts,"reuse_input",P.reuse_input_data);
	GetLongOption(opts,"no-norm",P.no_norm);

// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	if(P.reuse_input_data) P.lowmem=true;
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;
		if(P.lowmem) cout << "Low memory filtering method" << endl;
		if(P.reuse_input_data) cout << "Reusing input data as output" << endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		if(P.RealData==False) cout << "Transform in Complex space" << endl;   
		cout << "NbrScale = " <<  P.NbrScale    << endl;  
		cout << "NbrDir2d = " <<  P.NbrDir2d    << endl;  
		cout << "BlockSize = " <<  P.BlockSize    << endl;
		if(P.BlockOverlap) cout << "BlockOverlaping activated" << endl;
		if(P.lapped) cout << "Lapped transform" << endl;

		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;

		// Filtering
		if(P.SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(P.FilterType) <<endl;
								cout << " Noise level = " <<  P.SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  P.NSigma    << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
			if(P.force4sigma)	cout << " Force '4sigma lvl' at finest scales" << endl;
			if(P.no_coarse)		cout << " Set the coarse scale to zero" << endl;
			if(P.no_fine)		cout << " Set the fine scale to zero : " <<  Niter    << endl;
			if(P.Positivity)	cout << " The reconstruction is positive" << endl;
			if(P.use_min)		cout << " The reconstruction is > " <<P.min_value<< endl;
			if(P.use_max)		cout << " The reconstruction is < " <<P.max_value<< endl;
		}
		cout << endl;  
	}

	fltarray Data;
	fltarray ***TabBand;
	fltarray Recons;
	BFCurvelet3D *DataC = new BFCurvelet3D;
	char* filename;
	
	bool allocTB=true;

	if(Noise_calib & P.RealData)
	{
	// Input data
		Name_Imag_In = add_fits(Name_Imag_In);
		fits_read_fltarr(Name_Imag_In, Data);

	// FCurvelet initialisation
		DataC->alloc_from_coarse(Data,P);
		
	// Normalizing coefficients measurement
		DataC->estim_normalization(Name_Imag_Out);
	}
	else if(Compute_forward && Compute_recons && !Output_mr)
	{
	// Input data
		Name_Imag_In = add_fits(Name_Imag_In);
		fits_read_fltarr(Name_Imag_In, Data);

	// FCurvelet initialisation
		DataC->set_stat(Extract_stat);
		DataC->alloc_from_coarse(Data,P);
	
	// Transform, filtering and reconstruction
		DataC->filter(Data,Recons, P, Name_Imag_Out);

	// Save output
		filename = add_fits(Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	else 
	{
		if(Compute_forward)
		{
		// Input data
			Name_Imag_In = add_fits(Name_Imag_In);
			fits_read_fltarr(Name_Imag_In, Data);

		// FCurvelet initialisation
			DataC->set_stat(Extract_stat);
			DataC->alloc_from_coarse(Data,P);

		// Forward Transform
			DataC->transform(Data, TabBand, allocTB);
			allocTB=false;
		}
		else if(P.RealData)// We load the transformed cube
		{
			Name_Imag_In = add_mr(Name_Imag_In);
			DataC->read(Name_Imag_In,TabBand, P);
		}

//		DataC->temp(TabBand);

	// Statistic tools
		if(Extract_stat) DataC->extract_stat(TabBand, Name_Imag_Out);

	// Thresholding
		if( (P.SigmaNoise>=0) & P.RealData)
		{
			if(P.FilterType==FT_HARD || P.FilterType==FT_SOFT) DataC->threshold(TabBand, P);
			else if(P.FilterType==FT_WIENER) DataC->wiener(TabBand, P);
			else if(P.FilterType==FT_FDR) DataC->fdr(TabBand, P);
			else if(P.FilterType==FT_SBT) DataC->stein_block_threshold(TabBand, P);
			else if(P.FilterType==FT_CONTRAST) DataC->enhance(TabBand, P);
			else cerr<<"Filtering method '"<<string_filter_type(P.FilterType)<<"'not implemented yet";
		}

	// Save the transform
		if(Output_mr & P.RealData)
		{
		// Save the transform
			filename = add_mr(Name_Imag_Out);
			DataC->write(filename, TabBand);
		}

	// Reconstruction and save output
		if(Compute_recons) 
		{
			DataC->recons(TabBand, Recons);
			filename = add_fits(Name_Imag_Out);
			writefltarr(filename, Recons);
		}

	// Outputs the redundancy factor of the transform
	//	DataC->redundancy();
		
		for(int i=0;i<DataC->nbr_block();i++)
		{
			for(int s=0;s<P.NbrScale;s++)
				delete [] TabBand[i][s];
			delete [] TabBand[i];
		}
		delete [] TabBand;
	}

// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
	delete DataC;
	
	Data.~fltarray();
	Recons.~fltarray();
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





