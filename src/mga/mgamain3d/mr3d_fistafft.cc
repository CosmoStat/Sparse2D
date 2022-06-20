/******************************************************************************
**                   Copyright (C) 2009 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  15/12/09
**    
**    File:  mr3d_fistafft.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**    Applies the Fista algorithm to solve 
**    
**    min_alpha 1/2 || y - A D alpha ||_2 + ||alpha||_1
**    
**    where - y is the observed signal
**    		- D is the Dictionary, denoted here by "Domain" (here 2d wavelets)
**    		- A is the degradation operator.
**    			In this file, "degrade_fftmask2d" stands for A in the case of
**    			missing data in Fourier domain, using the array "Mask".
**    		
**    calling example : 
**    		program_name --mu 1 --tolvar 1e-6 -g 0.01 -i 200 -n 7 noisy_image mask out
**    		program_name --mu 1 --tolvar 1e-6 -g 0.01 -i 200 -n 7 fourier_data out
**    		
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "bfct.h"
#include "FCur_float.h"
#include "SB_Filter_float.h"
#include "Fista.h"

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
int Niter=1;

float SigmaNoise=1;
int NbrScale=3;
int NbrDir2d=16;
bool Positivity;

bool fast=false;
float mu = -1;
float TolVar=-1;
bool no_coarse = false;
bool use_init=false; 
bool decrease=false;

fltarray Mask;
cfarray cMask;
FFTN_3D *DataFT3D = new FFTN_3D;
FFTN_2D *DataFT2D = new FFTN_2D;

/***************************************/

static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options data mask result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
	manline();

	fprintf(OUTMAN, "         [-n NbrScale]\n");
	fprintf(OUTMAN, "             Number of scales of the transform\n");
	manline();

	fprintf(OUTMAN, "         [-g sigma][--sigmanoise sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholdings\n");
	manline();

	fprintf(OUTMAN, "         [-P][--positive] Positivity constraint\n");
	manline();

	fprintf(OUTMAN, "         [--mu mu] \n");
	fprintf(OUTMAN, "             The value of the step. default=1\n");
	manline();

	fprintf(OUTMAN, "         [--tolvar tv] \n");
	fprintf(OUTMAN, "             The stop criterion.\n");
	fprintf(OUTMAN, "             Stop when |x(n+1)-x(n)|_infty < tv\n");
	fprintf(OUTMAN, "             Default = 1e-6\n");
	manline();

	fprintf(OUTMAN, "         [--fast] Set to apply Fista\n");
	manline();

	fprintf(OUTMAN, "         [--decrease] Set to use a decreasing threshold (MCA)\n");
	manline();

	fprintf(OUTMAN, "         [--use-init] Set to use the input as an initial point.\n");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] To substract the coarse scale \n");
	fprintf(OUTMAN, "             before applying the algorithm (not avail. w all transf.)\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"vn:g:i:P")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose = True;
				break;
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);exit(-1); }
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'P': 
				Positivity = true;
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
	{
		cerr<<"Not enough parameters : missing Name_In Name_mask Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
		strcpy(Name_Mask, argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing Name_mask Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
		strcpy(Name_Imag_Out, argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing Name_Out"<<endl;
		usage(argv);
	}

	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}


/***************/

void degrade_mask(fltarray &x, bool degrade, bool reverse)
{
	if(degrade) x *= Mask;
}
void degrade_fftmask2d(cfarray &x, bool degrade, bool reverse)
{
	if(!reverse)
	{
		Icomplex_f xFrame;
		xFrame.alloc(x.buffer(), x.ny(),x.nx());
		DataFT2D->fftn2d(xFrame, False, true);
		if(degrade) x *= cMask;
	}
	else
	{
		if(degrade) x *= cMask;
		Icomplex_f xFrame;
		xFrame.alloc(x.buffer(), x.ny(),x.nx());
		DataFT2D->fftn2d(xFrame, True, true);
	}
}
	
int main(int argc, char *argv[])
{
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"SigmaNoise",SigmaNoise);
	GetLongOption(opts,"Niter",Niter);
	GetLongOption(opts,"positive",Positivity);
	GetLongOption(opts,"fast",fast);
	GetLongOption(opts,"decrease",decrease);
	GetLongOption(opts,"mu",mu);
	GetLongOption(opts,"tolvar",TolVar);
	GetLongOption(opts,"use-init",use_init);
	GetLongOption(opts,"no-coarse",no_coarse);

// Get short options
	filtinit(short_argc, short_argv);

// Print parameters
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(use_init) cout << "Using input as initial point" << endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale    << endl;  

		// Filtering
		if(SigmaNoise>0)
		{
								cout << " Noise level = " <<  SigmaNoise    << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
			if(Positivity)	cout << " The reconstruction is positive" << endl;
//			if(use_min)		cout << " The reconstruction is > " <<min_value<< endl;
//			if(use_max)		cout << " The reconstruction is < " <<max_value<< endl;
		}
		cout << endl;  
	}

// Input data
	fltarray Data;
	fltarray Recons;
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	strcat(Name_Mask,".fits");
	fits_read_fltarr(Name_Mask, Mask);
	cMask.alloc(Mask.nx(),Mask.ny(),0);
	for(int j=0;j<Mask.ny();j++)
	for(int i=0;i<Mask.nx();i++)
		cMask(j,i) = complex_f(Mask(j,i),0.);
/*
// UWT initialisation
	UOWT_float *Domain;
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Domain = new UOWT_float(SB1D);
	Domain->alloc(Data,NbrScale);

// FCurvelet initialisation
	BFCurvelet3D *Domain = new BFCurvelet3D;
	BFCurvelet3D_params P;
	P.Positivity = Positivity;
	P.NbrScale = NbrScale;
	P.SigmaNoise = SigmaNoise;
	Domain->alloc_from_coarse(Data,P);

// 2D FCurvelet initialisation
	FCUR_float *Domain = new FCUR_float();
	Domain->alloc(Data, NbrScale, NbrDir2d);
	Domain->get_norm_coeff(3);

// 2D Decimated Orthogonal Wavelet transform 
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Ortho_2D_WT_float *Domain = new Ortho_2D_WT_float(*SB1D);
	Domain->alloc(Data, NbrScale);
*/
// 2D Undecimated Orthogonal Wavelet transform 
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	PAVE_2D_WT_float *Domain = new PAVE_2D_WT_float(*SB1D);
	Domain->alloc(Data, NbrScale);
/*
// 2D Isotropic Atrou Wavelet transform 
	ATROUS_2D_WT_float *Domain = new ATROUS_2D_WT_float();
	Domain->alloc(Data, NbrScale);
	Domain->AdjointRec = True;
	//Domain->set_positive_coef(true);
*/
// Fista init
	Fista F(Domain);
	F.Decreasing = decrease;
	F.Fast = fast;
	F.Threshold = SigmaNoise;
	F.MaxNiter = Niter;
	F.No_coarse = no_coarse;
	F.Positivity = Positivity;
	if(mu>0) F.Mu = mu;
	if(TolVar>0) F.TolVar = TolVar;
	strcpy(F.NameOut,Name_Imag_Out);

// Fista algorithm
	cfarray cData(Data.nx(),Data.ny(),Data.nz()), cRecons(Data.nx(),Data.ny(),Data.nz());
	for(int j=0;j<Data.ny();j++) for(int i=0;i<Data.nx();i++) cData(j,i) = complex_f(Data(j,i),0.);
	degrade_fftmask2d(cData, !use_init, false);// degrade only if not using init
	F.run(cData, cRecons, &degrade_fftmask2d);
	degrade_fftmask2d(cRecons, false, true);
	
	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	delete Domain;
	Data.~fltarray();
	Recons.~fltarray();
	exit(0);
}





