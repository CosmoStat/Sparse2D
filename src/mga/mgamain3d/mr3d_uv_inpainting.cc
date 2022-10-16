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
**    		
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "FFTN_2D.h"
#include "GetLongOptions.h"
#include "Nesterov.h"
#include "Fista.h"
#include "FCur_float.h"
#include "SB_Filter_float.h"

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

// Transforms parameters
int NbrScale=3;
int NbrDir2d = 16;

// Algorithms
enum {ALG_ISTA, ALG_FISTA, ALG_MCA, ALG_NESTEROV} typedef type_algo;
type_algo Type_algo = ALG_FISTA;

// Domains
enum {DOM_DWT, DOM_UWT, DOM_IWT, DOM_FCT, DOM_TV} typedef type_domain;
type_domain Type_domain = DOM_IWT;

// Fista/ista/mca parameters
bool fast=false;
bool decrease=false;

// Nesterov parameters
float cont=1;

// Common parameters
float SigmaNoise=1;
float mu = -1;
float TolVar=-1;
bool no_coarse = false;
bool use_init=false; 
bool Positivity = false;

bool fourier_in = true;
fltarray Mask;
cfarray cMask;
FFTN_2D *DataFT2D = new FFTN_2D;
int dimension = 2; // dimension of the data

/***************************************/

static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options data (mask) result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
	manline();

	fprintf(OUTMAN, "         [-a algorithm] Choose the algorithm\n");
	fprintf(OUTMAN, "             1 : ISTA (Iterative Soft Thresholding)\n");
	fprintf(OUTMAN, "            [2]: FISTA (Fast ISTA)\n");
	fprintf(OUTMAN, "             3 : MCA (Morphological Component Analysis)\n");
	fprintf(OUTMAN, "             4 : NESTEROV\n");
	manline();

	fprintf(OUTMAN, "         [-t transform] Choose the representation\n");
	fprintf(OUTMAN, "             1 : Decimated Orthogonal Wavelets\n");
	fprintf(OUTMAN, "             2 : Undecimated Orthogonal Wavelets\n");
	fprintf(OUTMAN, "            [3]: Isotropic Undecimated Wavelets\n");
	fprintf(OUTMAN, "             4 : Fast Curvelets\n");
	fprintf(OUTMAN, "             5 : Total Variation\n");
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

	fprintf(OUTMAN, "         [--use-init] Set to use the input as an initial point.\n");
	manline();

	fprintf(OUTMAN, "         // Nesterov specific parameter\n");
	fprintf(OUTMAN, "         [--cont c] Specify the number of continuation steps.\n");
	fprintf(OUTMAN, "             Default = 1 (no continuation)\n");
	manline();

	fprintf(OUTMAN, "         // FISTA specific parameter\n");
	fprintf(OUTMAN, "         [--no-coarse] To substract the coarse scale \n");
	fprintf(OUTMAN, "             before applying the algorithm.\n");
	fprintf(OUTMAN, "             (not availaible with all transforms)\n");
	manline();

	exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"vn:g:i:Pa:t:")) != -1) 
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
			case 'a': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad number of algorithm: %s\n", OptArg);exit(-1);}
				Type_algo = (type_algo) (plop-1);
				break;
			case 't': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad number of transform: %s\n", OptArg);exit(-1);}
				Type_domain = (type_domain) (plop-1);
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
		cerr<<"Not enough parameters : missing Name_mask Name_Out"<<endl;
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
	GetLongOption(opts,"mu",mu);
	GetLongOption(opts,"tolvar",TolVar);
	GetLongOption(opts,"use-init",use_init);
	GetLongOption(opts,"cont",cont);
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
//			if(P.use_min)		cout << " The reconstruction is > " <<P.min_value<< endl;
//			if(P.use_max)		cout << " The reconstruction is < " <<P.max_value<< endl;
		}
		cout << endl;  
	}

// Input data
	fltarray Data; cfarray cData;
	fltarray Recons; cfarray cRecons;
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	
	if(Data.nz()) dimension = 3;
cerr<<"data nx ny nl nc "<<Data.nx()<<","<<Data.ny()<<","<<Data.nl()<<","<<Data.nc()<<endl;
	cData.alloc(Data.nx(),Data.ny(),Data.nz());
	cRecons.alloc(Data.nx(),Data.ny(),Data.nz());
cerr<<"crecons nx ny nl nc "<<cRecons.nx()<<","<<cRecons.ny()<<","<<cRecons.nl()<<","<<cRecons.nc()<<endl;
	if(dimension==2) for(int j=0;j<Data.ny();j++) for(int i=0;i<Data.nx();i++) cData(j,i) = complex_f(Data(j,i),0.);
	else for(int k=0;k<Data.nz();k++) for(int j=0;j<Data.ny();j++) for(int i=0;i<Data.nx();i++) cData(i,j,k) = complex_f(Data(i,j,k),0.);
	strcat(Name_Mask,".fits");
	fits_read_fltarr(Name_Mask, Mask);
	cMask.alloc(Mask.nx(),Mask.ny(),Mask.nz());
	if(dimension==2) for(int j=0;j<Mask.ny();j++) for(int i=0;i<Mask.nx();i++) cMask(j,i) = complex_f(Mask(j,i),0.);
	else for(int k=0;k<Data.nz();k++) for(int j=0;j<Data.ny();j++) for(int i=0;i<Data.nx();i++) cMask(i,j,k) = complex_f(Mask(i,j,k),0.);

// Generic transform
	FloatTrans *Domain;

// Specific transforms
	Ortho_2D_WT_float *Domain_dwt;
	PAVE_2D_WT_float *Domain_uwt;
	ATROUS_2D_WT_float *Domain_iwt;
	FCUR_float *Domain_fct;
	PAVE_2D_WT_float *Domain_tv;

// Filters fir the wavelet transforms and TV
	FilterAnaSynt SelectFilter;
	if(Type_domain == DOM_TV)
		SelectFilter.alloc(F_HAAR);
	else
		SelectFilter.alloc(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);

	switch(Type_domain)
	{
		case DOM_DWT:
			Domain_dwt = new Ortho_2D_WT_float(*SB1D);
			Domain_dwt->alloc(Data, NbrScale);
			Domain = Domain_dwt;
			break;
		case DOM_UWT:
			Domain_uwt = new PAVE_2D_WT_float(*SB1D);
			Domain_uwt->alloc(Data, NbrScale);
			Domain = Domain_uwt;
			break;
		case DOM_IWT:
			Domain_iwt = new ATROUS_2D_WT_float();
			Domain_iwt->alloc(Data, NbrScale);
			Domain_iwt->AdjointRec = True;
			//Domain_iwt->set_positive_coef(true);
			Domain = Domain_iwt;
			break;
		case DOM_FCT:
			Domain_fct = new FCUR_float();
			Domain_fct->alloc(Data, NbrScale, NbrDir2d);
			Domain_fct->get_norm_coeff(3);
			Domain = Domain_fct;
			break;
		case DOM_TV:
			Domain_tv = new PAVE_2D_WT_float(*SB1D);
			Domain_tv->alloc(Data, 2);
			Domain = Domain_tv;
			break;
	}

	if(Type_algo == ALG_ISTA || Type_algo == ALG_FISTA || Type_algo == ALG_MCA)
	{
	// Fista init
		Fista F(Domain);
		F.Decreasing = (Type_algo == ALG_MCA);
		F.Fast = (Type_algo == ALG_FISTA);
		F.Threshold = SigmaNoise;
		F.MaxNiter = Niter;
		F.No_coarse = no_coarse;
		F.Positivity = Positivity;
		if(mu>0) F.Mu = mu;
		if(TolVar>0) F.TolVar = TolVar;
		strcpy(F.NameOut,Name_Imag_Out);

	// Fista algorithm
		degrade_fftmask2d(cData, !use_init, false);// degrade only if not using init
		F.run(cData, cRecons, &degrade_fftmask2d);
		degrade_fftmask2d(cRecons, false, true);
	}
	else if(Type_algo == ALG_NESTEROV)
	{
	// Nesterov init
		Nesterov N(Domain);
		N.MaxNiter = Niter;
		N.Continuation = cont;
		N.Threshold = SigmaNoise;
		N.Positivity = Positivity;
		if(mu>0) N.Mu = mu;
		if(TolVar>0) N.TolVar = TolVar;
		strcpy(N.NameOut,Name_Imag_Out);

	// Nesterov algorithm
		degrade_fftmask2d(cData, !use_init, false);// degrade only if not using init
		N.run(cData, cRecons, &degrade_fftmask2d);
	}
	else 
	{
		cerr<<"Unknown method"<<endl;
	}

// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	delete Domain;
	Data.~fltarray();
	Recons.~fltarray();
	exit(0);
}





