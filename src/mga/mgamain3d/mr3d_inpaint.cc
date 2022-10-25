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
**    File:  mr3d_inpaint.cc
**
*******************************************************************************
**
**    DESCRIPTION  3D Morphological inpainting
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

// Curvelet 3D Parameters
	BFCurvelet3D_params P;

Bool Verbose=False;
Bool BlockOverlap=False;
int BlockSize = 0;
float SigmaNoise=1.;
float Alpha=0.01;
float NSigma=3;
bool Use3sig=false;
bool force4sigma=false;
bool UseCubeSigma=false;

int Niter=1;
int NbrScale = 3;
int NbrDir2d = 16;
int WienerBS=3;
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

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;//,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"vg:")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose=True;
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
		strcpy(Name_Mask, argv[OptInd++]);
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
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"SigmaNoise",SigmaNoise);
	
// Get short options
	filtinit(short_argc, short_argv);
	
	bool GVerbose = true;

// Iteration parameters	
	int Niter = 30;
	float Lmax=1000;
	float Lmin = 0;
	
// Transform parameters
	bool Use_init=false;
	bool Use_FC			=true;	bool Use_BC			=false;	bool Use_RC			=false;	bool Use_DC 		=false;
	int Nscale_FC		=3;		int Nscale_BC		=4;		int Nscale_RC		=4;		int skip_order		=0;
	int BlockSize_FC	=0; 	int BlockSize_BC	=17;	int BlockSize_RC	=17;	int BlockSize_DC	=0;
	Bool BlockOver_FC	=False;	Bool BlockOver_BC	=True;	Bool BlockOver_RC	=True;	Bool BlockOver_DC	=False;
	bool FC_no_coarse	=false;	bool BC_no_coarse	=false;	bool RC_no_coarse	=false;	
	int FC_imin			=0; 	int BC_imin			=0; 	int RC_imin			=0; 	int DC_imin			=0; 	
	int FC_imax			=1000;	int BC_imax			=1000;	int RC_imax			=1000;	int DC_imax			=1000;	
	bool FC_lowmem		=true;								
	
	bool Use_UW			=false;	bool Use_AW			=false;	
	int Nscale_UW		=4;		int Nscale_AW		=3; 	
	int UW_imin			=0; 	int AW_imin			=0; 	
	int UW_imax			=1000;	int AW_imax			=1000;	
								int AW_no_fine		=true;
								bool AW_del_single	=true;

	bool FC_decor_mask	=false; 							// weight with the FC transform of the mask
	bool TV_regul		=false;	float TV_lvl		= 10.0;	// TV penalty
	bool TV_mask_data	=false;								// Mask and data gradient l1 penalisation ||VM VX||_1
	
	use_min=true;	float min_value = 0.;
	use_max=true;	float max_value = 255;
	
// FC param init
	BFCurvelet3D_params P;
	P.SigmaNoise = SigmaNoise;
	P.NbrScale = Nscale_FC;
	P.BlockSize = BlockSize_FC;
	P.BlockOverlap = BlockOver_FC;
	P.no_coarse = FC_no_coarse;
	P.lowmem = FC_lowmem;
	P.NbrDir2d = NbrDir2d;
	
// Print parameters
	//if (Verbose)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout<<" Iterations parameters : "<<Niter<<" in ["<<Lmax<<","<<Lmin<<"]"<<endl;
		if(TV_regul) cout << "Using TV regul at lvl "<<TV_lvl<<endl;
		if(TV_mask_data) cout << "Using TV_mask_data correlation at lvl "<<TV_lvl<<endl;

		// Filtering
		if(SigmaNoise>0)
		{
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
		}
		cout << endl;  
	}

	fltarray Data, Orig, FullOrig, Mask, InvMask;
	fltarray Recons;

// Input data and mask
	char filename[64];
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	strcat(Name_Mask,".fits");
	fits_read_fltarr(Name_Mask, Mask);
	InvMask.alloc(Data.nx(),Data.ny(),Data.nz()); InvMask.init(1.);
	InvMask -= Mask;
	
// Filter type for BCurvelets and Wavelets
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	
// 3D transforms declaration
	BFCurvelet3D *DataFC = new BFCurvelet3D;
		fltarray *** TabBand_mask_FC;
	BCurvelet3D *DataBC = new BCurvelet3D;
	RCurvelet3D *DataRC = new RCurvelet3D;
		fltarray *TabBand_RC;
	
	UOWT *DataUW = new UOWT(SB1D);
		fltarray *TabBand_UW;
	ATROUS_3D_WT *DataAW = new ATROUS_3D_WT;
		fltarray *TabBand_AW;
		
	IM3D_DCT *DataDC = new IM3D_DCT;
		fltarray TabBand_DC;
	
// TV regularisation operator init
	FilterAnaSynt SelectFilter_b(F_HAAR);
	SubBandFilter *SB1D_h = new SubBandFilter(SelectFilter_b, NORM_L2);
	SB1D_h->Border = I_MIRROR;
	UOWT *TV = new UOWT(SB1D_h);
	fltarray *TabBand_TV,*TabBand_TV_mask,*TabBand_TV_temp;
	fltarray Recons_TV;
	TV->init(TabBand_TV, Data.nx(), Data.ny(), Data.nz(),2);
	Recons_TV.alloc(Data.nx(),Data.ny(),Data.nz());
	if(TV_mask_data)
	{
		UOWT *TV_mask = new UOWT(SB1D_h);
		TV_mask->init(TabBand_TV_mask, Data.nx(), Data.ny(), Data.nz(),2);
	// Estimate the gradient of the mask
		TV->transform(Mask,TabBand_TV_mask,2);// Nscale=2
		TabBand_TV_temp=new fltarray[8];
		for(int s=0;s<8;s++)
			TabBand_TV_temp[s]=TabBand_TV_mask[s];
	}

// Separated data
	Orig = Data * Mask;
	if(Use_init) FullOrig = Data;
	Mask.~fltarray();
	Recons.alloc(Data.nx(),Data.ny(),Data.nz()); Recons.init(0.);
	
// Main separation loop
	bool first=true;
	float lambda;
	for(int kk=0;kk<2;kk++)
	for(int iter=0;iter<Niter;iter++)
	{
	// Threshold level
		lambda = 1-iter*1.0/float(Niter-1);
	// fast, slow
		lambda = (pow((double)lambda,3)+lambda/25.)/1.04;
	// slow, fast, slow
		//lambda = (1.-cos(lambda*3.1416)+lambda/10.)/1.1;
		lambda = Lmin+ lambda*(Lmax-Lmin);
		if(kk==0) cerr<<lambda<<" ";
		else 
		{
		// Transforms
			{
			// A Trou wavelets
				if(Use_AW && iter>=AW_imin && iter<=AW_imax)
				{
					if(GVerbose) cerr<<"ATrouWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

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
					DataAW->recons(TabBand_AW, Recons, Nscale_AW);

				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_awt_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// Undecimated wavelets
				if(Use_UW && iter>=UW_imin && iter<=UW_imax)
				{
					if(GVerbose) cerr<<"UWT... "<<lambda<<endl;
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

				// UOWT initialisation
					if(iter==UW_imin)
						DataUW->init(TabBand_UW, Data.nx(), Data.ny(), Data.nz(),Nscale_UW);

				// Forward Transform
					DataUW->transform(Data,TabBand_UW,Nscale_UW);

				// Filtering
					if(FilterType==FT_HARD || FilterType==FT_SOFT)
						DataUW->threshold(TabBand_UW, SigmaNoise, lambda, FilterType, force4sigma);
					else if(FilterType==FT_WIENER) DataUW->wiener(TabBand_UW, SigmaNoise, 3);
					else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";

				// Reconstruction
					DataUW->recons(TabBand_UW,Recons,Nscale_UW);

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
				if(Use_FC && iter>=FC_imin && iter<=FC_imax)
				{
					if(GVerbose) cerr<<"BFCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

				// Initialization
					if(iter==FC_imin)
					{
						if(FC_decor_mask)
						{
							// May not work anymore
							if(GVerbose) cerr<<"Transforming the mask...";
							BFCurvelet3D *MaskFC = new BFCurvelet3D;
							MaskFC->alloc_from_coarse(Data,P);
							MaskFC->transform(InvMask, TabBand_mask_FC, true);
							delete MaskFC;
							
							DataFC->set_extern_norm(TabBand_mask_FC);
						}
						DataFC->alloc_from_coarse(Data,P);
					}

				// Transform, filtering and reconstruction
					P.NSigma = lambda;
					DataFC->filter(Data,Recons,P,Name_Imag_Out);
				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_bfct_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}

			// Beam Curvelets
				if(Use_BC && iter>=BC_imin && iter<=BC_imax)
				{
					if(GVerbose) cerr<<"BCur..."<<lambda<<endl;
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

				// Initialization
					if(iter==BC_imin)
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
					DataBC->filter(Data, Recons, SigmaNoise, lambda, FilterType, Alpha, WienerBS, force4sigma, UseCubeSigma, Name_Imag_Out);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_bcur_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}


			// Ridge Curvelets
				if(Use_RC && iter>=RC_imin && iter<=RC_imax)
				{
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

					if(GVerbose) cerr<<"RCur..."<<lambda<<endl;

				// Initialization
					if(iter==RC_imin)
					{
						type_wavelet3D TypeW3D = DEF_TYPE_W3D;//2
						type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1
						DataRC->CurTrans=CurTrans;
						DataRC->BlockOverlap=BlockOver_RC;
						DataRC->set_no_coarse(RC_no_coarse);
						DataRC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, TypeW3D, Use3sig);
					}

				// Forward Transform
					DataRC->transform(Data, TabBand_RC, (iter==RC_imin));

				// Filtering
					if(FilterType==FT_HARD || FilterType==FT_SOFT) DataRC->threshold(TabBand_RC, SigmaNoise, lambda, FilterType, force4sigma, UseCubeSigma);
					else if(FilterType==FT_WIENER) DataRC->wiener(TabBand_RC, SigmaNoise, 3, UseCubeSigma);
					else if(FilterType==FT_FDR) DataRC->fdr(TabBand_RC, Alpha, SigmaNoise);

				// Reconstruction
					DataRC->recons(TabBand_RC, Recons);

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
				if(Use_DC && iter>=DC_imin && iter<=DC_imax)
				{
					if(GVerbose) cerr<<"DCT... "<<lambda<<endl;
					Data = Recons*InvMask + Orig ;
					if(first && Use_init) Data += InvMask*FullOrig;
					first=false;

				// DCT initialisation
					if(iter==DC_imin)
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
					DataDC->recons(TabBand_DC,Recons);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_dct_recons_%03d.fits", Name_Imag_Out, iter);
						writefltarr(filename, Recons);
					}
				}
				
			// TV regularisation
				if(TV_regul)
				{
					TV->transform(Recons,TabBand_TV,2);// Nscale=2
					TV->threshold(TabBand_TV, SigmaNoise, TV_lvl, FT_SOFT, force4sigma);
					TV->recons(TabBand_TV,Recons,2);
				}

			// L1 norm on the cross gradient || grad(M) grad(X) ||1
				if(TV_mask_data)
				{
					float mu = 10.;
					float e = 0.9; // percentage kept when grad=mu/2
					float alpha = mu*mu*(1.-e)/(4.*e);
					TV->transform(Recons,TabBand_TV,2);// Nscale=2
					float val;
					for(int s=0;s<7;s++)
						for(int x=0;x<TabBand_TV_mask[s].nx();x++)
						for(int y=0;y<TabBand_TV_mask[s].ny();y++)
						for(int z=0;z<TabBand_TV_mask[s].nz();z++)
						{
//							val = (TabBand_TV[s])(x,y,z)*(TabBand_TV_mask[s])(x,y,z);
//							val = ((val >0) - (val <0)) * val*val/(val*val+alpha); // smoothed L1 norm derivative
//							(TabBand_TV_temp[s])(x,y,z) = mu * (TabBand_TV_mask[s])(x,y,z) * val;
							val = (TabBand_TV[s])(x,y,z)*(TabBand_TV_mask[s])(x,y,z);
							val -= ::soft_threshold(val,TV_lvl);
							(TabBand_TV_temp[s])(x,y,z) = val;
						}
					TV->recons(TabBand_TV_temp,Recons_TV,2);
					Recons = Recons - Recons_TV;

					//TV->transform(Recons,TabBand_TV,2);// Nscale=2
					//TV->threshold_on_mask(TabBand_TV, TabBand_TV_mask, SigmaNoise, lambda, FT_SOFT, force4sigma);
					//TV->recons(TabBand_TV,Recons,2);
					}

			}// end transforms
		}
	}
	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
	delete DataFC;
	delete DataBC; 
	delete DataRC; if(Use_RC) delete [] TabBand_RC;
	delete DataUW; if(Use_UW) delete [] TabBand_UW;
	delete DataDC; TabBand_DC.~fltarray();
	delete SB1D;
	
	Data.~fltarray();
	Orig.~fltarray();
	FullOrig.~fltarray();
	InvMask.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





