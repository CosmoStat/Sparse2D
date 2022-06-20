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
	double start,stop;
	start=clock();

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

// Iteration parameters	
	int Niter	= 10;
	float Lmax	= 1.0;
	float Lmin	= 0.0;
	
// Transform parameters
	bool regul_TV = false;
	bool regul_TV_last =false;	float Nsig_last_TV	= 0.5;
	
	// Fast Curvelets			BeamCurvelets				RidCurvelets				Discrete cosine
	bool Use_FC			=false;	bool Use_BC			=false;	bool Use_RC			=false;	bool Use_DC 		=false;
	int Nscale_FC		=2;		int Nscale_BC		=3;		int Nscale_RC		=4;		int skip_order		=0;
	int BlockSize_FC	=0; 	int BlockSize_BC	=17;	int BlockSize_RC	=17;	int BlockSize_DC	=0;
	Bool BlockOver_FC	=False;	Bool BlockOver_BC	=False;	Bool BlockOver_RC	=True;	Bool BlockOver_DC	=False;
	bool FC_no_coarse	=false;	bool BC_no_coarse	=false;	bool RC_no_coarse	=false;	
	int FC_NbrDir2d		=16;
	
	// Undecim ortho wavelet	Atrou wavelets				
	bool Use_UW			=false;	bool Use_AW			=true;	
	int Nscale_UW		=4;		int Nscale_AW		=4; 	
	
// Levels
	float Nsig_orig=3.;
	float th_far = SigmaNoise/2.;
	use_max=false;  float min_value = -00.;
	use_min=false;  float max_value = 255.;
	
// **********************************************
// **********************************************
// 				CFA initialization
// **********************************************
// **********************************************
	
// FC param init
	BFCurvelet3D_params P;
	P.SigmaNoise = SigmaNoise;
	P.NbrScale = Nscale_FC;
	P.BlockSize = BlockSize_FC;
	P.BlockOverlap = BlockOver_FC;
	P.no_coarse = FC_no_coarse;
	P.NbrDir2d = FC_NbrDir2d;

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

// Input and output data
	fltarray Data, Orig, FullOrig, Mask, InvMask;
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
	BFCurvelet3D *DataFC = new BFCurvelet3D;
	BCurvelet3D *DataBC = new BCurvelet3D;
	RCurvelet3D *DataRC = new RCurvelet3D;
	UOWT *DataUW = new UOWT(SB1D);
	ATROUS_3D_WT *DataAW = new ATROUS_3D_WT;
	IM3D_DCT *DataDC = new IM3D_DCT;
	
	fltarray ***TabBand_FC,	***TabBand_FC_orig;
	fltarray *TabBand_BC,	*TabBand_BC_orig;
	fltarray *TabBand_RC,	*TabBand_RC_orig;
	fltarray *TabBand_UW,	*TabBand_UW_orig;
	fltarray *TabBand_AW,	*TabBand_AW_orig;
	fltarray TabBand_DC,	TabBand_DC_orig;
	
// TV regularisation operator init
	FilterAnaSynt SelectFilter_b(F_HAAR);
	SubBandFilter *SB1D_h = new SubBandFilter(SelectFilter_b, NORM_L2);
	SB1D_h->Border = I_CONT;
	UOWT *TV = new UOWT(SB1D_h);
	fltarray *TabBand_TV;
	if(regul_TV || regul_TV_last)
		TV->init(TabBand_TV, Data.nx(), Data.ny(), Data.nz(),2);
	
// **********************************************
// **********************************************
// 				CFA main program
// **********************************************
// **********************************************
	
// Original transforms
	if(Use_AW)
	{
	// AWT initialisation
		DataAW->alloc(TabBand_AW, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
		DataAW->alloc(TabBand_AW_orig, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
		
	// Forward Transform
		DataAW->transform(Data, TabBand_AW_orig, Nscale_AW);

	// Filtering
		DataAW->threshold(TabBand_AW_orig, SigmaNoise*Nsig_orig, FT_HARD);
		DataAW->normalize_self(TabBand_AW_orig);

		if(GVerbose)
		{
		// Reconstruction
			DataAW->normalize_self(TabBand_AW_orig, true);
			DataAW->recons(TabBand_AW_orig, Recons, Nscale_AW);
			DataAW->normalize_self(TabBand_AW_orig);

		// Save output
			sprintf(filename,"%s_awt_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
	}

	if(Use_UW)
	{
	// UOWT initialisation
		DataUW->init(TabBand_UW, Data.nx(), Data.ny(), Data.nz(),Nscale_UW);
		DataUW->init(TabBand_UW_orig, Data.nx(), Data.ny(), Data.nz(),Nscale_UW);
		
	// Forward Transform
		DataUW->transform(Data,TabBand_UW_orig,Nscale_UW);

	// Filtering
		DataUW->threshold(TabBand_UW_orig, SigmaNoise, Nsig_orig, FT_HARD, force4sigma);

		if(GVerbose)
		{
		// Reconstruction
			DataUW->recons(TabBand_UW_orig,Recons,Nscale_UW);

		// Save output
			sprintf(filename,"%s_uwt_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
	}

	if(Use_FC)
	{
	// FC Initialization
		DataFC->alloc_from_coarse(Data,P);

	// Forward Transform
		bool Alloc=true;
		DataFC->transform(Data, TabBand_FC_orig, Alloc);
		
	// Filtering
		P.SigmaNoise = SigmaNoise;
		P.NSigma = Nsig_orig;
		DataFC->threshold(TabBand_FC_orig, P);
		
		if(GVerbose)
		{
		// Reconstruction
			DataFC->recons(TabBand_FC_orig, Recons);

		// Save output
			sprintf(filename,"%s_bfct_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
	}

	if(Use_BC)
	{
	// BC Initialization
		type_linelet3d_WTtrans CurTrans = DEF_LIN3D_TRANS;
		type_wavelet3D TypeW3D = DEF_TYPE_W3D;
		DataBC->CurTrans=CurTrans;
		DataBC->BlockOverlap=BlockOver_BC;
		DataBC->set_stat(false);
		DataBC->set_no_recons(false);
		DataBC->set_3sigma(Use3sig);
		DataBC->set_no_coarse(BC_no_coarse);
		DataBC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_BC, BlockSize_BC, TypeW3D, SB1D);

	// Forward Transform
		DataBC->transform(Data, TabBand_BC_orig, true);
		
	// Filtering
		DataBC->threshold(TabBand_BC_orig, SigmaNoise, Nsig_orig, FT_HARD, force4sigma, UseCubeSigma);
		
		if(GVerbose)
		{
		// Reconstruction
			DataBC->recons(TabBand_BC_orig, Recons);

		// Save output
			sprintf(filename,"%s_bcur_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
		DataBC->normalize_self(TabBand_BC_orig, false, UseCubeSigma);
	}

	if(Use_RC)
	{
	// RC Initialization
		type_wavelet3D TypeW3D = DEF_TYPE_W3D;//2
		type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1
		DataRC->CurTrans=CurTrans;
		DataRC->BlockOverlap=BlockOver_RC;
		DataRC->set_no_coarse(RC_no_coarse);
		DataRC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, TypeW3D, Use3sig);
		
	// Forward Transform
		DataRC->transform(Data, TabBand_RC_orig, true);

	// Filtering
		DataRC->threshold(TabBand_RC_orig, SigmaNoise, Nsig_orig, FT_HARD, force4sigma, UseCubeSigma);
		DataRC->normalize_self(TabBand_RC_orig, false, UseCubeSigma);

		if(GVerbose)
		{
		// Reconstruction
			DataRC->normalize_self(TabBand_RC_orig, true, UseCubeSigma);
			DataRC->recons(TabBand_RC_orig, Recons);
			DataRC->normalize_self(TabBand_RC_orig, false, UseCubeSigma);

		// Save output
			sprintf(filename,"%s_rcur_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
	}

	if(Use_DC)
	{
	// DCT initialisation
		DataDC->set_skip_order(skip_order);
		DataDC->init(Data.nx(),Data.ny(),Data.nz(),BlockSize_DC,BlockOver_DC);
		
	// Forward Transform
		DataDC->transform(Data,TabBand_DC_orig);
		
	// Filtering
		DataDC->threshold(TabBand_DC_orig, SigmaNoise, Nsig_orig, FT_HARD);
		DataDC->normalize_self(TabBand_DC_orig);

		if(GVerbose)
		{
		// Reconstruction
			DataDC->normalize_self(TabBand_DC_orig, true);
			DataDC->recons(TabBand_DC_orig,Recons);
			DataDC->normalize_self(TabBand_DC_orig);

		// Save output
			sprintf(filename,"%s_dct_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
	}
	
	
	Recons=fltarray(Data.nx(),Data.ny(),Data.nz());
	Recons.init(0.0);
	
// Main denoising loop
	float lambda = Lmax;
	float delta = (Lmax-Lmin)/float(Niter-1);
	float pos = Nsig_orig*SigmaNoise/2.;		// positive detection in original
	for(int kk=0;kk<2;kk++) // to print the lambda
	for(int iter=0;iter<Niter;iter++)
	{
	// Threshold level
		lambda=(Niter-1-iter)*delta+Lmin;
		
	// Print the lambda, then main loop
		if(kk==0) cerr<<lambda<<" ";
		else 
		{
		// Transforms
			{
			// A Trou wavelets
				if(Use_AW)
				{
					if(GVerbose) cerr<<"ATrouWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				// Forward Transform
					DataAW->transform(Recons, TabBand_AW, Nscale_AW);

				// Filtering
					if(!regul_TV)
						DataAW->threshold(TabBand_AW, SigmaNoise*lambda, true); // true for FT_SOFT
					DataAW->normalize_self(TabBand_AW);

				// Update the coefficients
float cnt=0;
float cntt=0;
					for(int s=0;s<Nscale_AW-1;s++)
					for(int i=0;i<Data.nx();i++)
					for(int j=0;j<Data.ny();j++)
					for(int k=0;k<Data.nz();k++)
						if( abs( (TabBand_AW_orig[s])(i,j,k) ) > pos ) // !=0 but for float data
						{
							if( abs( (TabBand_AW[s])(i,j,k)-(TabBand_AW_orig[s])(i,j,k) ) > th_far )
							{
								cnt++;
								(TabBand_AW[s])(i,j,k) = (TabBand_AW_orig[s])(i,j,k);
							}
							cntt++;
						}
cerr<<"AW maj : "<<cnt<<"/"<<cntt<<endl;
					TabBand_AW[Nscale_AW-1] = TabBand_AW_orig[Nscale_AW-1];
					
				// Reconstruction
					DataAW->normalize_self(TabBand_AW, true);
					DataAW->recons(TabBand_AW, Recons, Nscale_AW);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_awt_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// Undecimated wavelets
				if(Use_UW)
				{
					if(GVerbose) cerr<<"UWT... "<<lambda<<endl;

				// Forward Transform
					DataUW->transform(Recons,TabBand_UW,Nscale_UW);

				// Filtering
					if(!regul_TV)
						DataUW->threshold(TabBand_UW, SigmaNoise, lambda, FT_SOFT, force4sigma);

				// Update the coefficients
					for(int s=0;s<7*Nscale_UW-7;s++)
					for(int i=0;i<Data.nx();i++)
					for(int j=0;j<Data.ny();j++)
					for(int k=0;k<Data.nz();k++)
						if( abs( (TabBand_UW_orig[s])(i,j,k) ) > pos ) // !=0 but for float data
						if( abs( (TabBand_UW[s])(i,j,k)-(TabBand_UW_orig[s])(i,j,k) ) > th_far )
							(TabBand_UW[s])(i,j,k) = (TabBand_UW_orig[s])(i,j,k);
					TabBand_UW[Nscale_UW-1] = TabBand_UW_orig[Nscale_UW-1];
					
				// Reconstruction
					DataUW->recons(TabBand_UW,Recons,Nscale_UW);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_uwt_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// Fast Curvelets
				if(Use_FC)
				{
					if(GVerbose) cerr<<"BFCT... "<<lambda<<endl;

				// Forward Transform
					DataFC->transform(Recons, TabBand_FC, iter==0);

				// Filtering
					if(!regul_TV)
					{
						P.NSigma = lambda;
						DataFC->threshold(TabBand_FC, P);
					}
		
					float cnt=0, cntt=0;
					
				// Update the coefficients
					for(int B=0;B<DataFC->nbr_block();B++)
					{
						for(int s=0;s<Nscale_FC-1;s++)
						for(int b=0;b<DataFC->nbr_band(s);b++)
						for(int i=0;i<TabBand_FC[B][s][b].nx();i++)
						for(int j=0;j<TabBand_FC[B][s][b].ny();j++)
						for(int k=0;k<TabBand_FC[B][s][b].nz();k++)
							if( abs( (TabBand_FC_orig[B][s][b])(i,j,k) ) > pos ) // !=0 but for float data
							{
								cntt++;
								if( abs( (TabBand_FC[B][s][b])(i,j,k)-(TabBand_FC_orig[B][s][b])(i,j,k) ) > th_far )
								{
									cnt++;
									(TabBand_FC[B][s][b])(i,j,k) = (TabBand_FC_orig[B][s][b])(i,j,k);
								}
							}
	cerr<<"FC maj : "<<cnt<<" "<<"/"<<cntt<<endl;
						TabBand_FC[B][Nscale_AW-1][0] = TabBand_FC_orig[B][Nscale_AW-1][0];
					}
					
				// Reconstruction
					DataFC->recons(TabBand_FC, Recons);
					
				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_bfct_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// BeamCurvelets
				if(Use_BC)
				{
					if(GVerbose) cerr<<"BCur..."<<lambda<<endl;

				// Forward Transform
					DataBC->transform(Recons, TabBand_BC, iter==0);

				// Filtering
					if(!regul_TV)
						DataBC->threshold(TabBand_BC, SigmaNoise, lambda, FT_SOFT, force4sigma, UseCubeSigma);
					DataBC->normalize_self(TabBand_BC, 0, UseCubeSigma);

				// Update the coefficients
					for(int s=0;s<Nscale_BC-1;s++)
					for(int i=0;i<TabBand_BC[s].nx();i++)
					for(int j=0;j<TabBand_BC[s].ny();j++)
					for(int k=0;k<TabBand_BC[s].nz();k++)
						if( abs( (TabBand_BC_orig[s])(i,j,k) ) > pos ) // !=0 but for float data
						if( abs( (TabBand_BC[s])(i,j,k)-(TabBand_BC_orig[s])(i,j,k) ) > th_far )
							(TabBand_BC[s])(i,j,k) = (TabBand_BC_orig[s])(i,j,k);
					TabBand_BC[Nscale_BC-1] = TabBand_BC_orig[Nscale_BC-1];
					
				// Reconstruction
					DataBC->normalize_self(TabBand_BC, 1, UseCubeSigma);
					DataBC->recons(TabBand_BC, Recons);
					
				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_bcur_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}


			// RidCurvelets
				if(Use_RC)
				{
					if(GVerbose) cerr<<"RCur..."<<lambda<<endl;

				// Forward Transform
					DataRC->transform(Recons, TabBand_RC, (iter==0));

				// Filtering
					if(!regul_TV)
						DataRC->threshold(TabBand_RC, SigmaNoise, lambda, FT_SOFT, force4sigma, UseCubeSigma);
					DataRC->normalize_self(TabBand_RC, false, UseCubeSigma);

				// Update the coefficients
					for(int s=0;s<Nscale_RC-1;s++)
					for(int i=0;i<TabBand_RC[s].nx();i++)
					for(int j=0;j<TabBand_RC[s].ny();j++)
					for(int k=0;k<TabBand_RC[s].nz();k++)
						if( abs( (TabBand_RC_orig[s])(i,j,k) ) > pos ) // !=0 but for float data
						if( abs( (TabBand_RC[s])(i,j,k)-(TabBand_RC_orig[s])(i,j,k) ) > th_far )
							(TabBand_RC[s])(i,j,k) = (TabBand_RC_orig[s])(i,j,k);
					TabBand_RC[Nscale_RC-1] = TabBand_RC_orig[Nscale_RC-1];
					
				// Reconstruction
					DataRC->normalize_self(TabBand_RC, true, UseCubeSigma);
					DataRC->recons(TabBand_RC, Recons);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_rcur_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// Discrete cosine transform
				if(Use_DC)
				{
					if(GVerbose) cerr<<"DCT... "<<lambda<<endl;

				// Forward Transform
					DataDC->transform(Recons,TabBand_DC);

				// Filtering
					if(!regul_TV)
						DataDC->threshold(TabBand_DC, SigmaNoise, lambda, FT_SOFT);
					DataDC->normalize_self(TabBand_DC);

				// Update the coefficients
					for(int i=0;i<TabBand_DC.nx();i++)
					for(int j=0;j<TabBand_DC.ny();j++)
					for(int k=0;k<TabBand_DC.nz();k++)
						if( abs( (TabBand_DC_orig)(i,j,k) ) > pos ) // !=0 but for float data
						if( abs( (TabBand_DC)(i,j,k)-(TabBand_DC_orig)(i,j,k) ) > th_far )
							(TabBand_DC)(i,j,k) = (TabBand_DC_orig)(i,j,k);
					
				// Reconstruction
					DataDC->normalize_self(TabBand_DC, true);
					DataDC->recons(TabBand_DC,Recons);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_dct_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}
			}
		}
	}
	
// Last TV regularisation
	if(regul_TV_last)
	{
		TV->transform(Recons,TabBand_TV,2);// Nscale=2
		TV->threshold(TabBand_TV, SigmaNoise, Nsig_last_TV, FT_HARD, force4sigma);
		TV->recons(TabBand_TV,Recons,2);
	}

// Save output
	sprintf(filename,"%s_recons.fits", Name_Imag_Out);
	writefltarr(filename, Recons);

	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
	if(Use_FC)
	{
		for(int B=0;B<DataFC->nbr_block();B++)
		{
			for(int s=0;s<Nscale_FC;s++)
			{
				delete [] TabBand_FC[B][s];
				delete [] TabBand_FC_orig[B][s];
			}
			delete [] TabBand_FC[B];
			delete [] TabBand_FC_orig[B];
		}
		delete [] TabBand_FC;
		delete [] TabBand_FC_orig;
	}
	delete DataFC;
	delete DataBC; 
	if(Use_BC)
	{
		delete [] TabBand_BC;
		delete [] TabBand_BC_orig;
	}
	delete DataRC;
	if(Use_RC) 
	{
		delete [] TabBand_RC;
		delete [] TabBand_RC_orig;
	}
	delete DataUW;
	if(Use_UW) 
	{
		delete [] TabBand_UW;
		delete [] TabBand_UW_orig;
	}
	delete DataAW;
	if(Use_AW) 
	{
		delete [] TabBand_AW;
		delete [] TabBand_AW_orig;
	}
	delete DataDC;
	TabBand_DC.~fltarray();
	TabBand_DC_orig.~fltarray();
	
	delete TV;
	if(regul_TV) 
		delete [] TabBand_TV;
	
	delete SB1D;
	delete SB1D_h;
	
	Data.~fltarray();
	Orig.~fltarray();
	FullOrig.~fltarray();
	InvMask.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





