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
**    File:  mr3d_mca.cc
**
*******************************************************************************
**
**    DESCRIPTION  3D Morphological component analysis : morphological separation
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
#include "MeyerWT1D.h"
#include "FCur_TFrame.h"

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

// Curvelet 3D Parameters
	BFCurvelet3D_params P;

Bool Verbose=False;
Bool BlockOverlap=False;
int BlockSize = 0;
float SigmaNoise=-1;
float Alpha=0.01;
float NSigma=3;
bool Use3sig=false;
bool force4sigma=false;
bool UseCubeSigma=false;

int Niter=1;
int NbrScale = 3;
int NbrDir2d = 16;

int WienerBS=3;

bool no_coarse=false;
bool no_fine=false;
bool positivity=false;
bool use_max=false; float max_value=0;
bool use_min=false; float min_value=0;

type_extraction extract_type = DEF_TYPE_EXTRACTION;
filter_type FilterType = FT_HARD;

bool ChangeTabDir=false;

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

	fprintf(OUTMAN, "         [--use3sig] Use the cumulative at lvl 1-2.7e-3\n");
	fprintf(OUTMAN, "             to estimate the normalisation coefficients\n");
	fprintf(OUTMAN, "             instead of the variance\n");
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

	fprintf(OUTMAN, "         [-I]\n");
	fprintf(OUTMAN, "             Isotropic WT (not yet).\n");
	manline();

	fprintf(OUTMAN, "         [-d Ndir2d]\n");
	fprintf(OUTMAN, "             number of directions (in 2d), min=12, def=16.\n");
	manline();
	
	fprintf(OUTMAN, "         [--overlap][-O] BlockOverlap\n");
	manline();

	fprintf(OUTMAN, "         [--lowmem] Low memory usage\n");
	fprintf(OUTMAN, "             Doesn't evaluate the full transform,\n");
	fprintf(OUTMAN, "             use for simple filtering (block independent)\n");
	manline();

	fprintf(OUTMAN, "         [--ulowmem] Ultra low memory usage\n");
	fprintf(OUTMAN, "             Doesn't evaluate the full transform,\n");
	fprintf(OUTMAN, "             use for simple filtering (band independent)\n");
	fprintf(OUTMAN, "             even when there is only one big block.\n");
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
	while ((c = GetOpt(argc,argv,(char*)"n:d:vg:a:s:f:i:COB:P")) != -1) 
	{
		switch (c) 
		{
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
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'O': 
				BlockOverlap = True;
				break;
			case 'B': 
				if (sscanf(OptArg,"%d",&BlockSize) != 1) 
				{fprintf(OUTMAN, "Error: bad Block size: %s\n", OptArg);exit(-1);}
				break;
			case 'P': 
				positivity = true;
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
	it = opts.find("--overlap");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	BlockOverlap=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--Niter");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Niter; }
	it = opts.find("--hard");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_HARD; }
	it = opts.find("--soft");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_SOFT; }
	it = opts.find("--wiener");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_WIENER; }
	it = opts.find("--fdr");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_FDR; }
	it = opts.find("--stein");
	if(it!=opts.end()){ istringstream ss(it->second); FilterType = FT_SBT; }
	it = opts.find("--use3sig");
	if(it!=opts.end()){ istringstream ss(it->second); Use3sig=true; }
	it = opts.find("--force4sigma");
	if(it!=opts.end()){ istringstream ss(it->second); force4sigma=true; }
	it = opts.find("--no-fine");
	if(it!=opts.end()){ istringstream ss(it->second); no_fine=true; }
	it = opts.find("--no-coarse");
	if(it!=opts.end()){ istringstream ss(it->second); no_coarse=true; }
	it = opts.find("--positive");
	if(it!=opts.end()){ istringstream ss(it->second); positivity=true; }
	it = opts.find("--min-value");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>min_value; use_min=true; }
	it = opts.find("--max-value");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>max_value; use_max=true; }
	
// Get short options
	filtinit(short_argc, short_argv);

// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale    << endl;  
		cout << "NbrDir2d = " <<  NbrDir2d    << endl;  
		cout << "BlockSize = " <<  BlockSize    << endl;
		if(BlockOverlap) cout << "BlockOverlaping activated" << endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(FilterType) <<endl;
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
			if(force4sigma)		cout << " Force '4sigma lvl' at finest scales" << endl;
			if(Use3sig)			cout << " Use 3sigma cumulative estimation" << endl;
			if(no_coarse)		cout << " Set the coarse scale to zero" << endl;
			if(no_fine)			cout << " Set the fine scale to zero : " <<  Niter    << endl;
			if(positivity)		cout << " The reconstruction is positive" << endl;
			if(use_min)			cout << " The reconstruction is > " <<min_value<< endl;
			if(use_max)			cout << " The reconstruction is < " <<max_value<< endl;
		}
		cout << endl;  
	}


// ********************
//        START
// ********************
	
// to print the current loop
	bool GVerbose = true;

// MCA parameters	
	int Niter = 5;
	float Lmax= 10;
	float Lmin = 3;
	
// Transform parameters
	bool Use_FC			=false;	bool Use_BC			=true;	bool Use_RC			=true;	bool Use_DC			=false;
	int Nscale_FC		=4;		int Nscale_BC		=4;		int Nscale_RC		=4;		int skip_order		=0;
	int BlockSize_FC	=0; 	int BlockSize_BC	=17;	int BlockSize_RC	=17;	int BlockSize_DC	=16;
	Bool BlockOver_FC	=False;	Bool BlockOver_BC	=True;	Bool BlockOver_RC	=True;	Bool BlockOver_DC	=False;
	bool FC_no_coarse	=true;	bool BC_no_coarse	=true;	bool RC_no_coarse	=true;	
	
	bool Use_UW			=false;	bool Use_MW			=false;	bool Use_AW			=true;	
	int Nscale_UW		=5;		int Nscale_MW		=5; 	int Nscale_AW		=4; 	
															bool AW_no_fine		=false;
															bool AW_del_single	=true;
															
	bool Use_TC			=false;
	int Nscale_TC		=4;
	int TC_Nscale1D		=3;
	int TC_NDir2d		=16;
	int TC_kill_time_f	=1;
	
// FC param init
	BFCurvelet3D_params P;
	P.SigmaNoise = SigmaNoise;
	P.NbrScale = Nscale_FC;
	P.BlockSize = BlockSize_FC;
	P.BlockOverlap = BlockOver_FC;
	P.no_coarse = FC_no_coarse;

// Print parameters
	//if (Verbose)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout<<" Iterations parameters : "<<Niter<<" in ["<<Lmax<<","<<Lmin<<"]"<<endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  P.NSigma    << endl;
		}
		cout << endl;  
	}

// Filter type for BCurvelets and Wavelets
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	
// Input data
	fltarray Data, Orig;
	char filename[64];
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	
// Transformed data
	fltarray Recons, Recons_tmp;

// 3D transforms declaration
	BFCurvelet3D *DataFC = new BFCurvelet3D;
	FCur_TFrame *DataTC = new FCur_TFrame;		fltarray ***TabBand_TC;
	BCurvelet3D *DataBC = new BCurvelet3D;
	RCurvelet3D *DataRC = new RCurvelet3D;		fltarray *TabBand_RC;
	UOWT *DataUW = new UOWT(SB1D);				fltarray *TabBand_UW;
	MEYER_WT3D *DataMW = new MEYER_WT3D;		fltarray *TabBand_MW;
	ATROUS_3D_WT *DataAW = new ATROUS_3D_WT;	fltarray *TabBand_AW;
	IM3D_DCT *DataDC = new IM3D_DCT;			fltarray TabBand_DC;
	
// Separated data
	Orig = Data;
	Recons.alloc(Data.nx(),Data.ny(),Data.nz()); Recons.init(0.);
	Recons_tmp.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_tmp.init(0.);
	
// Print the threshold levels
	float lambda;
	double maxcoef;
	if(GVerbose)
	{
		for(int iter=0;iter<Niter;iter++)
		{
		// Threshold level
			lambda = 1-iter*1.0/float(Niter-1);
			lambda = (pow((double)lambda,3.)+lambda/10.)/1.1;
			lambda = Lmin+ lambda*(Lmax-Lmin);
			cerr<<lambda<<"\t";
		}
		cerr<<endl;
	}
	
// Main separation loop
	for(int iter=0;iter<Niter;iter++)
	{
	// Threshold level
		lambda = 1-iter*1.0/float(Niter-1);
		lambda = (pow((double)lambda,3.)+lambda/10.)/1.1;
		lambda = Lmin+ lambda*(Lmax-Lmin);
	
	// Transforms
		{
		// Meyer wavelets
			if(Use_MW)
			{
				if(GVerbose) cerr<<"UWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;
			
			// MWT initialisation
				if(iter==0)
					DataMW->init(Nscale_MW, Data.nx(), Data.ny(), Data.nz(), False, False);

			// Forward Transform
				DataMW->transform(Data, TabBand_MW);
				
			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT)
					DataMW->threshold(TabBand_MW, SigmaNoise, lambda, FilterType, !force4sigma);
				else if(FilterType==FT_WIENER) DataMW->wiener(TabBand_MW, SigmaNoise, 3);
				else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";

			// Reconstruction
				DataMW->recons(TabBand_MW,Recons_tmp,False);
				
			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_mwt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_mwt_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
			
		// A Trou wavelets
			if(Use_AW)
			{
				if(GVerbose) cerr<<"ATrouWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;
			
			// AWT initialisation
				if(iter==0)
				{
					DataAW->set_no_fine(AW_no_fine);
					DataAW->alloc(TabBand_AW, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
				}

			// Forward Transform
				DataAW->transform(Data, TabBand_AW, Nscale_AW);
				
			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT)
					DataAW->threshold(TabBand_AW, SigmaNoise*lambda, FilterType==FT_SOFT);
				else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";
				if(AW_del_single) DataAW->clean_single(TabBand_AW, SigmaNoise*lambda);

			// Reconstruction
				DataAW->recons(TabBand_AW, Recons_tmp, Nscale_AW);

			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_awt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_awt_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
			
		// Undecimated wavelets
			if(Use_UW)
			{
				if(GVerbose) cerr<<"UWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;
			
			// UOWT initialisation
				if(iter==0)
					DataUW->init(TabBand_UW, Data.nx(), Data.ny(), Data.nz(),Nscale_UW);

			// Forward Transform
				DataUW->transform(Data,TabBand_UW,Nscale_UW);
				
			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT)
					DataUW->threshold(TabBand_UW, SigmaNoise, lambda, FilterType, force4sigma);
				else if(FilterType==FT_WIENER) DataUW->wiener(TabBand_UW, SigmaNoise, 3);
				else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";

			// Reconstruction
				DataUW->recons(TabBand_UW,Recons_tmp,Nscale_UW);

			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_uwt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_uwt_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
			
		// Time Curvelets : 2D Curvelets + Time Wavelets
			if(Use_TC)
			{
				Data = Orig - Recons;

				if(GVerbose) cerr<<"TCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
			// Initialization
				if(iter==0)
				{
					DataTC->alloc_from_coarse(Nscale_TC, TC_Nscale1D, Data.nx(), Data.ny(), Data.nz(), TC_NDir2d, False, False, True);
					DataTC->set_kill_time_fine(TC_kill_time_f);
					if(SigmaNoise>=0) DataTC->get_norm_coeff(3);
				}

			// Forward Transform
				DataTC->transform(Data, TabBand_TC, (iter==0));

			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT) DataTC->threshold(TabBand_TC, SigmaNoise, lambda, FilterType);

			// Reconstruction
				DataTC->recons(TabBand_TC, Recons_tmp);
				
			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_tcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_tcur_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
			
		// Beam Curvelets
			if(Use_BC)
			{
				if(GVerbose) cerr<<"BCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;

			// Initialization
				if(iter==0)
				{
					type_linelet3d_WTtrans CurTrans = DEF_LIN3D_TRANS;
					type_wavelet3D TypeW3D = DEF_TYPE_W3D;
					DataBC->CurTrans=CurTrans;
					DataBC->BlockOverlap=BlockOver_BC;
					DataBC->set_stat(false);
					DataBC->set_no_recons(false);
					DataBC->set_3sigma(Use3sig);
					DataBC->set_no_coarse(BC_no_coarse);
					DataBC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_BC, BlockSize_BC, TypeW3D, SB1D);
				}

			// Transform, filtering and reconstruction
				DataBC->filter(Data, Recons_tmp, SigmaNoise, lambda, FilterType, Alpha, WienerBS, force4sigma, UseCubeSigma, Name_Imag_Out);

			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_bcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_bcur_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}


		// Ridge Curvelets
			if(Use_RC)
			{
				Data = Orig - Recons;

				if(GVerbose) cerr<<"RCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
			// Initialization
				if(iter==0)
				{
					type_wavelet3D TypeW3D = DEF_TYPE_W3D;//2
					type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1
					DataRC->CurTrans=CurTrans;
					DataRC->BlockOverlap=BlockOver_RC;
					DataRC->set_no_coarse(RC_no_coarse);
					DataRC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, TypeW3D, Use3sig);
				}

			// Forward Transform
				DataRC->transform(Data, TabBand_RC, (iter==0));

			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT) DataRC->threshold(TabBand_RC, SigmaNoise, lambda, FilterType, force4sigma, UseCubeSigma);
				else if(FilterType==FT_WIENER) DataRC->wiener(TabBand_RC, SigmaNoise, 3, UseCubeSigma);
				else if(FilterType==FT_FDR) DataRC->fdr(TabBand_RC, Alpha, SigmaNoise);

			// Reconstruction
				DataRC->recons(TabBand_RC, Recons_tmp);
				
			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_rcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_rcur_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
			
		// Fast Curvelets
			if(Use_FC)
			{
				if(GVerbose) cerr<<"BFCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;
			
			// Initialization
				if(iter==0)
				{
					//DataFC->set_lowmem(false);
					//DataFC->set_positivity(positivity);
					//if(use_min) DataFC->set_min(min_value);
					//if(use_max) DataFC->set_max(max_value);
					//DataFC->set_no_coarse(FC_no_coarse);
					//DataFC->set_stat(true);
					//DataFC->alloc_from_coarse(Nscale_FC,Data.nx(), Data.ny(), Data.nz(), NbrDir2d, BlockSize_FC, BlockOver_FC, True, extract_type);
					DataFC->alloc_from_coarse(Data,P);

				// Noise calibration if not already done, and if useful here
					bool AllocTB = true;//(iter ==0);
					if(!DataFC->isset_tabsigma())
					{
						fltarray ***TabBand_FC;
					// Normalizing coefficients measurement
						DataFC->estim_normalization(NULL);
						DataFC->alloc_from_coarse(P);
					// Free
						DataFC->dealloc(TabBand_FC);
					}
				}
								
			// Transform, filtering and reconstruction
				//DataFC->filter(Data,Recons_tmp,SigmaNoise, lambda, FilterType, Alpha, WienerBS, force4sigma, Name_Imag_Out);
				P.NSigma = lambda;
				DataFC->filter(Data,Recons_tmp,P,Name_Imag_Out);
			maxcoef = DataFC->get_max_coef();

			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_bfct_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_bfct_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
				if(GVerbose)
				{
					double c;
					fltarray err = (Orig-Recons);
					c = err.energy() + maxcoef;
					cerr<<"  Criterion : "<<c<<" ("<<err.energy()<<"+"<<maxcoef<<")"<<endl;
				}
			}

		// Discrete cosine transform
			if(Use_DC)
			{
				if(GVerbose) cerr<<"DCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				Data = Orig - Recons;
			
			// DCT initialisation
				if(iter==0)
				{
					DataDC->set_skip_order(skip_order);
					DataDC->init(Data.nx(),Data.ny(),Data.nz(),BlockSize_DC,BlockOver_DC);
				}

			// Forward Transform
				DataDC->transform(Data,TabBand_DC);
				
			// Filtering
				if(FilterType==FT_HARD || FilterType==FT_SOFT)
					DataDC->threshold(TabBand_DC, SigmaNoise, lambda, FilterType);
				else if(FilterType==FT_WIENER) DataDC->wiener(TabBand_DC, SigmaNoise, 3);
				else if(FilterType==FT_FDR) DataDC->fdr(TabBand_DC, Alpha, SigmaNoise);
				else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";
				
			// Reconstruction
				DataDC->recons(TabBand_DC,Recons_tmp);
				
			// Update Reconstruction
				Recons += Recons_tmp;

			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_dct_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_tmp);
					sprintf(filename,"%s_dct_recons_%03d_total.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons);
				}
			}
		}
	}
	
	
// Save output
	sprintf(filename,"%s_recons.fits",Name_Imag_Out);
	writefltarr(filename, Recons);
	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
	delete DataFC;
	if(Use_TC) { for(int i=0;i<Nscale_TC;i++) { for(int j=0;j<DataTC->nbr_bands(i);j++) delete [] TabBand_TC[i][j]; delete [] TabBand_TC[i]; } delete [] TabBand_TC; }
	delete DataTC; 
	delete DataBC; 
	delete DataRC; if(Use_RC) delete [] TabBand_RC;
	delete DataAW; if(Use_AW) delete [] TabBand_AW;
	delete DataUW; if(Use_UW) delete [] TabBand_UW;
	delete DataMW; if(Use_MW) delete [] TabBand_MW;
	delete DataDC; TabBand_DC.~fltarray();
	delete SB1D;
	
	Data.~fltarray();
	Orig.~fltarray();
	Recons.~fltarray();
	Recons_tmp.~fltarray();
	
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





