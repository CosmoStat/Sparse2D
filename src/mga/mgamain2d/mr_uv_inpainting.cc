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
**    File:  mr_uv_inpainting.cc
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
#include "FFTN_2D.h"
#include "GetLongOptions.h"
#include "writefits3d.h"
#include "Nesterov.h"
#include "Fista.h"
#include "MCA.h"
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

float SigmaNoise=1;

// Transforms parameters
int NbrScale=3;
int NbrDir2d = 16;
type_sb_filter SB_Filter = F_MALLAT_7_9;

// Algorithms
enum {ALG_ISTA, ALG_FISTA, ALG_MCA, ALG_NESTEROV} typedef type_algo;
type_algo Type_algo = ALG_FISTA;

// Domains
enum {DOM_DWT, DOM_UWT, DOM_IWT, DOM_FCT, DOM_TV} typedef type_domain;
type_domain Type_domain = DOM_IWT;

// Fista/ista/mca parameters
bool fast=false;
bool decrease=false;
bool UseMadNorm=false;

// Nesterov parameters
float cont=1;

// Common parameters
float mu = -1;
float TolVar=-1;
bool no_coarse = false;
bool use_init=false; 
bool Positivity = true;

bool fourier_in = true;
fltarray Mask;
cfarray cMask;
FFTN_2D *DataFT2D = new FFTN_2D;
float ksigma=3.;

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

	fprintf(OUTMAN, "         [-p][--positive] Remove the Positivity constraint\n");
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
	while ((c = GetOpt(argc,argv,(char*)"vn:g:i:pa:t:s:")) != -1) 
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
			case 'p': 
				Positivity = (Positivity == true) ? false: true;
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
			case 's': 
				if (sscanf(OptArg,"%f",&ksigma) != 1) 
				{fprintf(OUTMAN, "Error: bad ksigma noise level: %s\n", OptArg);exit(-1);}
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
		cerr<<"Not enough parameters : missing Name_In (Name_mask) Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
		strcpy(Name_Mask, argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing (Name_mask) Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
	{
		strcpy(Name_Imag_Out, argv[OptInd++]);
		fourier_in = false;
	}
	else
	{
		// Fourier data
		strcpy(Name_Imag_Out, Name_Mask);
		fourier_in = true;
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
   bool Normalize=true;
   
	if(reverse == false)
	{
        // Forward transform
 		DataFT2D->fftn2d(x, False, Normalize);
		if(degrade) x *= cMask;
 	}
	else
	{
        // Inverse transform
		if(degrade) x *= cMask;
 		DataFT2D->fftn2d(x, True, Normalize);
	}
}
    
/********  *MU CALCULATION ********/ 

// Helper functions. Many of these are overloaded to deal with either complex or real-valued arrays.
// fltarray is an array of floats, cfarray is an array of complex.

static double inner_product(fltarray data1, fltarray data2) {
  double sum = 0;

  for (int i=0; i<data1.nx(); ++i)        	
    for (int j=0; j<data1.ny(); ++j) {
     sum += data1(i,j)*data2(i,j);         

    }
  return sum;
}

// 2-norm
static double norm2(cfarray data) {
  complex_f sum = 0;
  for (int i=0; i<data.nx(); ++i)          
    for (int j=0; j<data.ny(); ++j) {
      sum += (data(i,j)*conj((data(i,j))));
    }

  return sqrt(sum.real());
}

static double norm2(fltarray data) {
  double sum = 0;
  for (int i=0; i<data.nx(); ++i)          
    for (int j=0; j<data.ny(); ++j) {
      sum += data(i,j)*data(i,j);
    }

  return sqrt(sum);
}

static void divide_array(cfarray& data, double num) {
   for (int i=0; i<data.nx(); ++i)          
    for (int j=0; j<data.ny(); ++j) {
      data(i,j) /= num;
    }
}

static void divide_array(fltarray& data, double num) {
   for (int i=0; i<data.nx(); ++i)          
    for (int j=0; j<data.ny(); ++j) {
      data(i,j) /= num;
    }
}

// Normalize so norm2 is 1 afterwards.
// Returns original norm.
static double normalize(cfarray& data) {
  double n2 = norm2(data);
  divide_array(data,n2);
  return n2;
}

static double normalize(fltarray& data) {
  double n2 = norm2(data);
  divide_array(data,n2);
  return n2;
}

static void copy(cfarray data1, fltarray& data2) {
  for (int i=0; i<data1.nx(); ++i)          
    for (int j=0; j<data1.ny(); ++j) data2(i,j) = data1(i,j).real();
}

static void copy(fltarray data1, cfarray& data2) {
  for (int i=0; i<data1.nx(); ++i)          
    for (int j=0; j<data1.ny(); ++j) data2(i,j) = complex_f(data1(i,j),0);
}

static void copy(fltarray data1, fltarray& data2) {
  for (int i=0; i<data1.nx(); ++i)          
    for (int j=0; j<data1.ny(); ++j) data2(i,j) = data1(i,j);
}

// This applies the matrix/operator A, in this case a masked FFT
static void image2fft(fltarray im, cfarray& fft) {
  copy(im,fft);
  degrade_fftmask2d(fft,true,false);
}

//This applies the transpose of A
static void fft2image(cfarray fft, fltarray& im) {
  degrade_fftmask2d(fft,true,true);
  copy(fft,im);
}
 
 
 // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.34.9868&rep=rep1&type=pdf
// The entry point. The mask is used to mask the FFT, and hence defines the operators, along 
// with the Fourier transform.
double mu_estimate(int nit_max) {
  fltarray init_vec, vec2;
  cfarray vec1;
  double vec_1n;
  int k;
  
  // Prepare and store. Must have mask in cMask.
  init_vec.alloc(cMask.nx(),cMask.nx(),0); vec1.alloc(cMask.nx(),cMask.nx(),0); vec2.alloc(cMask.nx(),cMask.nx(),0); 
  //copy(mask,c_mask);

  
  // Create a test image
  for (int i=0; i<init_vec.nx(); ++i)
    for (int j=0; j<init_vec.ny(); ++j) {
      init_vec(i,j) = rand();
  }

  // Run the loop
  for (k=0; k<nit_max; ++k) {
    normalize(init_vec);

    image2fft(init_vec, vec1);	

    vec_1n = normalize(vec1);
    
    fft2image(vec1, vec2);		
    
    if ( norm2(vec2) <= inner_product(vec2,init_vec) ) break;
    
    copy(vec2, init_vec);
  }
  if ( k == nit_max ) cout << "Max iterations reached in relaxation_cycle\n";
  return 2/(1.05*vec_1n*vec_1n);
}
/**********  END of mu estimate **************/

int main(int argc, char *argv[])
{
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"SigmaNoise",SigmaNoise);
	GetLongOption(opts,"Niter",Niter);
    bool nopos=false;
	GetLongOption(opts,"positive",nopos);
    if (nopos == true) Positivity = false;
 
	GetLongOption(opts,"mu",mu);
	GetLongOption(opts,"tolvar",TolVar);
	GetLongOption(opts,"use-init",use_init);
	GetLongOption(opts,"cont",cont);
	GetLongOption(opts,"no-coarse",no_coarse);
	GetLongOption(opts,"MAD",UseMadNorm);

// Get short options
	filtinit(short_argc, short_argv);

// Print parameters
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(use_init) cout << "Using input as initial point" << endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale    << endl;  
        if (UseMadNorm == true) cout << "Use MAD" << endl;
        if (Positivity== true) cout << "Use Positivity constraint" << endl;
        else cout << "NO Positivity constraint" << endl;
        
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
	cData.alloc(Data.nx(),Data.ny());
//cerr<<"data nx ny nz "<<Data.nx()<<","<<Data.ny()<<","<<Data.nz()<<endl;
	cRecons.alloc(Data.nx(),Data.ny());
	if(fourier_in)
	{
        int Symx, Symy;
        fltarray T1,T2;
        fits_read_cfarr2d(Name_Imag_In, cData);
        int nx = cData.nx();
        int ny = cData.ny();        
		float norm = 1.F/sqrt(float(nx*ny));
        
        // Renormalize the Fourier components
        for(int j=0;j<ny;j++) 
        for(int i=0;i<nx;i++) cData(i,j) *= norm;
            

        Data.alloc(nx,ny);
               
		// The fourier transform is assumed to be in IDL format :
		// not normalized, and with opposite sign in the exponentials with respect to here.
		cMask.alloc(nx,ny);

  		for(int j=0;j<ny;j++) 
        for(int i=0;i<nx;i++)
		{
            if ( abs( cData(i,j)) != 0) cMask(i,j) = complex_f(1.,0.);
       /*     Symx = nx - i;
            Symy = ny - j;
            if ((Symy >= 0) && (Symx >= 0) && (Symx < nx ) && (Symy < ny))
            {
               cData(Symx,Symy) = complex_f(cData(i,j).real(), - cData(i,j).imag());
               cMask(Symx,Symy) = cMask(i,j);
            } */
            
		//	cerr<<","<<cMask(i,j);//(Data(i,j,0) != 0);
		}
 
       // fits_write_cfarr2d ("xxv.fits", cData);
       // fits_write_cfarr2d ("xxt.fits", cMask);

		// put the data back into direct space
		degrade_fftmask2d(cData, false, true); 
       // fits_write_cfarr2d ("xxvd.fits", cData);

	}
	else
	{
    	fits_read_fltarr(Name_Imag_In, Data);
        cData.alloc(Data.nx(),Data.ny());
        
//cerr<<"cdata nx ny "<<cData.nx()<<","<<cData.ny()<<","<<cData.nz()<<endl;
		for(int j=0;j<Data.ny();j++) for(int i=0;i<Data.nx();i++) cData(i,j) = complex_f(Data(i,j),0.);
		strcat(Name_Mask,".fits");
		fits_read_fltarr(Name_Mask, Mask);
		cMask.alloc(Mask.nx(),Mask.ny(),0);
		for(int j=0;j<Mask.ny();j++) for(int i=0;i<Mask.nx();i++) cMask(i,j) = complex_f(Mask(i,j),0.);
	}
    if (mu == 0) mu = mu_estimate(100) ;
    if (Verbose) cout <<  "Estimated mu = " << mu << endl;

// Generic transform
	FloatTrans *Domain;

// Specific transforms
	Ortho_2D_WT_float *Domain_dwt;
	PAVE_2D_WT_float *Domain_uwt;
	ATROUS_2D_WT_float *Domain_iwt;
	FCUR_float *Domain_fct;
	PAVE_2D_WT_float *Domain_tv;

// Filters for the wavelet transforms and TV
	FilterAnaSynt SelectFilter;
	if(Type_domain == DOM_TV)
		SelectFilter.alloc(F_HAAR);
	else
		SelectFilter.alloc(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);

	switch(Type_domain)
	{
		case DOM_DWT:
            if (Verbose) cout <<  "Use DWT " << endl;
			Domain_dwt = new Ortho_2D_WT_float(*SB1D);
			Domain_dwt->alloc(Data, NbrScale);
			Domain = Domain_dwt;
            Domain_iwt->UseMad = UseMadNorm;
			break;
		case DOM_UWT:
            if (Verbose) cout <<  "Use UWT " << endl;
			Domain_uwt = new PAVE_2D_WT_float(*SB1D);
			Domain_uwt->alloc(Data, NbrScale);
            Domain_uwt->UseMad = UseMadNorm;
			Domain = Domain_uwt;
			break;
		case DOM_IWT:
            if (Verbose) cout <<  "Use Starlet " << endl;
			Domain_iwt = new ATROUS_2D_WT_float();
			Domain_iwt->alloc(Data, NbrScale);
			// Domain_iwt->AdjointRec = True;
            Domain_iwt->UseMad = UseMadNorm;
			//Domain_iwt->set_positive_coef(true); 
			Domain = Domain_iwt;
			break;
		case DOM_FCT:
            if (Verbose) cout <<  "Use Curvelet " << endl;
			Domain_fct = new FCUR_float();
			Domain_fct->alloc(Data, NbrScale, NbrDir2d);
			Domain_fct->get_norm_coeff(3);
			Domain = Domain_fct;
			break;
		case DOM_TV:
            if (Verbose) cout <<  "TV " << endl;
			Domain_tv = new PAVE_2D_WT_float(*SB1D);
			Domain_tv->alloc(Data, 2);
			Domain = Domain_tv;
			break;
	}

	if(Type_algo == ALG_ISTA || Type_algo == ALG_FISTA)
	{
	// Fista init
		Fista F(Domain);
		F.Decreasing = false;		// don't need this now
		F.Fast = (Type_algo == ALG_FISTA);
		F.SigmaNoise = SigmaNoise;
		F.MaxNiter = Niter;
		F.No_coarse = no_coarse;
		F.Positivity = Positivity;
        F.ksigma = ksigma;
		if(mu>0) F.Mu = mu;
		if(TolVar>0) F.TolVar = TolVar;
		strcpy(F.NameOut,Name_Imag_Out);

	// Fista algorithm
//fltarray datareal(cData.nx(),cData.ny(),cData.nz());
//datareal.init(-10);
//for(int j=0;j<cData.ny();j++) for(int i=0;i<cData.nx();i++) { datareal(i,j) = cData(i,j).real(); }//cerr<<cData(i,j)<<",";}
//writefltarr("out_data.fits", datareal);
		degrade_fftmask2d(cData, !use_init, false);// degrade only if not using init
		F.run(cData, cRecons, &degrade_fftmask2d);
		degrade_fftmask2d(cRecons, false, true);
	} else if (Type_algo == ALG_MCA)
	{
        MCA M(Domain); 		
		M.Threshold = SigmaNoise;
		M.SigmaNoise = SigmaNoise;
		M.MaxNiter = Niter;
		M.No_coarse = no_coarse;
		M.Positivity = Positivity;
        M.ksigma = ksigma;
		if(mu>0) M.Mu = mu;
		if(TolVar>0) M.TolVar = TolVar;
		strcpy(M.NameOut,Name_Imag_Out);

		degrade_fftmask2d(cData, !use_init, false);// degrade only if not using init
		M.run(cData, cRecons, &degrade_fftmask2d); 
		degrade_fftmask2d(cRecons, false, true);
	}else if(Type_algo == ALG_NESTEROV)
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





