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
#include "Oriented_DCT3D.h"
#include "DCT_Frame.h"
#include "MeyerWT1D.h"
#include "FCur_TFrame.h"
#include "FCur_SFrame.h"
#include "MR2D1D.h"

char Name_Imag_In[256]={}; /* input file image */
char Name_Imag_Out[256]={}; /* output file name */
char Name_Mask[256]={}; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

/********************************/
//		Input variables			//
/********************************/

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
bool use_max=false;
bool use_min=false;
type_extraction extract_type = DEF_TYPE_EXTRACTION;
filter_type FilterType = FT_HARD;
enum type_mca {	MCA_RESIDUAL,	// threshold the residual : Y-sum(Xi)
				MCA_TOTAL		// threshold the residual + previous estimate of current kth transform : 
				};				//		Y-sum(Xi)+Xk = Y-sum(Xi)_{i!=k}

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
	
	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;

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

	if (OptInd < argc-1)
	{
		strcpy(Name_Mask, argv[OptInd++]);
		strcpy(Name_Imag_Out, argv[OptInd++]);
	}
	else if (OptInd < argc)
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

// ********************
//        START
// ********************
	
// to print the current loop
	bool GVerbose = true;
	filter_type WT2D1D_FilterType, FC_FilterType, BC_FilterType, RC_FilterType, DC_FilterType, OT_FilterType, DT_FilterType, TC_FilterType, SC_FilterType, UW_FilterType, MW_FilterType, AW_FilterType;

// MCA parameters	
	int Niter = 30;
	float Lmax= 30;
	float Lmin = 3;
	float lambda;
	type_mca MCA_type = MCA_TOTAL;// MCA_TOTAL MCA_RESIDUAL
	
// Transform parameters
//  FAST Curvelet 3D               Beam Curvelet               RidCurvelet                  DCT                      Oriented Texture             DCT2D per frame
	bool Use_FC			=false;	bool Use_BC			=false;	bool Use_RC			=false;	bool Use_DC			=false;	bool Use_OT			=false;	bool Use_DT			=false;	
	FC_FilterType	=FT_HARD;	BC_FilterType	=FT_HARD;	RC_FilterType	=FT_HARD;	DC_FilterType	=FT_HARD;	OT_FilterType	=FT_HARD;	DT_FilterType	=FT_HARD;	
	int Nscale_FC		=3;		int Nscale_BC		=4;		int Nscale_RC		=4;		int skip_order_DC	=3; 	int skip_order_OT	=3;		int skip_order_DT	=5;		
	int BlockSize_FC	=0; 	int BlockSize_BC	=17;	int BlockSize_RC	=17;	int BlockSize_DC	=16;	int BlockSize_OT	=17;	int BlockSize_DT	=32;	
	Bool BlockOver_FC	=False;	Bool BlockOver_BC	=True;	Bool BlockOver_RC	=True;	Bool BlockOver_DC	=True;	Bool BlockOver_OT	=True;	Bool BlockOver_DT	=True;	
	bool FC_no_coarse	=true;	bool BC_no_coarse	=true;	bool RC_no_coarse	=true;	
	bool FC_th_coarse	=false;	bool BC_keep_energy =true;
	bool FC_lowmem		=true;								
	bool FC_regul_TV	=true;
	
//  Curvelet2D+WT1D                  2D Curvelets + average
	bool Use_TC			=true;	bool Use_SC			=false;	
	TC_FilterType	=FT_HARD;	SC_FilterType	=FT_HARD;	
	int Nscale_TC		=4; 	int Nscale_SC		=4; 	
	int TC_Nscale1D		=1;  // ==> no WT1D transform								
	int TC_NDir2d		=16; 	int SC_NDir2d		=16; 	
	int TC_kill_time_f	=1; 	bool SC_th_coarse	=true;
	
//  Undecimated WT              3D Meyer WT                 Starlet 3D     
	bool Use_UW			=false;	bool Use_MW			=false;	bool Use_AW			=false;	
	UW_FilterType	=FT_HARD;	MW_FilterType	=FT_HARD;	AW_FilterType	=FT_HARD;	
	int Nscale_UW		=3;		int Nscale_MW		=3; 	int Nscale_AW		=3;		
															bool AW_no_fine		=false;	
	use_min=false;	float min_value = 0.;
	use_max=false;	float max_value = 255.;
	

    // Curvelet2D+WT1D            
    bool Use_WT2D1D			=true;
    WT2D1D_FilterType	=FT_HARD;
    int WT2D1D_NScale1D		=2;
    int WT2D1D_NScale2d		=16; 
    type_transform WT2D1D_Trans2D = TO_PAVE_BSPLINE;
    
// FC param init
	BFCurvelet3D_params P;
	P.FilterType = FC_FilterType;
	P.SigmaNoise = SigmaNoise;
	P.NbrScale = Nscale_FC;
	P.BlockSize = BlockSize_FC;
	P.BlockOverlap = BlockOver_FC;
	P.no_coarse = FC_no_coarse;
	P.threshold_coarse = FC_th_coarse;
	P.lowmem = FC_lowmem;
	P.NbrDir2d = NbrDir2d;
	
// Filter type for BCurvelets and Wavelets
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	
// Print parameters
	if (GVerbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		if(Name_Mask[0]!='\0') cout << "File Name mask = " << Name_Mask << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "MCA_TYPE = " << (int(MCA_type) == int(MCA_RESIDUAL) ? "MCA_RESIDUAL" : "MCA_TOTAL")  << endl;   

		// Filtering
		cout << " Noise level = " <<  SigmaNoise    << endl;
		cout << " NSigma Threshold = " <<  NSigma    << endl;
		cout << " Number of iterations : " <<  Niter    << endl;
		cout << endl;  
	}

// Input data
	fltarray Data, Orig, Mask, InvMask;
	char filename[64];
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	bool Use_mask = Name_Mask[0]!='\0';
	if(Use_mask)
	{
		strcat(Name_Mask,".fits");
		fits_read_fltarr(Name_Mask, Mask);
		InvMask.alloc(Data.nx(),Data.ny(),Data.nz()); InvMask.init(1.);
		InvMask -= Mask;
		Orig = Data*Mask;
	}
	else 
		Orig = Data;
	
// Transformed data
	fltarray Recons, Recons_WT2D1D, Recons_FC, Recons_TC, Recons_SC, Recons_BC, Recons_RC, Recons_UW, Recons_MW, Recons_DC, Recons_OT, Recons_DT, Recons_AW;

// 3D transforms declaration
	BFCurvelet3D *DataFC = new BFCurvelet3D;
	FCur_TFrame *DataTC = new FCur_TFrame;		fltarray ***TabBand_TC;
	FCur_SFrame *DataSC = new FCur_SFrame;		Ifloat **TabBand_SC;
	BCurvelet3D *DataBC = new BCurvelet3D;
	RCurvelet3D *DataRC = new RCurvelet3D;		fltarray *TabBand_RC;
	UOWT *DataUW = new UOWT(SB1D);				fltarray *TabBand_UW;
	MEYER_WT3D *DataMW = new MEYER_WT3D;		fltarray *TabBand_MW;
	ATROUS_3D_WT *DataAW = new ATROUS_3D_WT;	fltarray *TabBand_AW;
	IM3D_DCT *DataDC = new IM3D_DCT;			fltarray TabBand_DC;
	Oriented_DCT3D *DataOT = new Oriented_DCT3D;fltarray TabBand_OT; // Oriented Texture
	DCT_Frame *DataDT = new DCT_Frame;			fltarray TabBand_DT; // 2D DCT by frame
    // Curvelet2D+WT1D            
    MR2D1D *DataWT2D1D = new MR2D1D;            fltarray *TabBand_WT2D1D;
	
// Separated data
	Recons.alloc(Data.nx(),Data.ny(),Data.nz()); Recons.init(0.);
	Recons_MW.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_MW.init(0.);
	Recons_AW.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_AW.init(0.);
	Recons_UW.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_UW.init(0.);
	Recons_TC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_TC.init(0.);
	Recons_SC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_SC.init(0.);
	Recons_BC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_BC.init(0.);
	Recons_RC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_RC.init(0.);
	Recons_FC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_FC.init(0.);
	Recons_DC.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_DC.init(0.);
	Recons_OT.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_OT.init(0.);
	Recons_DT.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_DT.init(0.);
    Recons_WT2D1D.alloc(Data.nx(),Data.ny(),Data.nz()); Recons_WT2D1D.init(0.);

// TV regularisation operator init
	FilterAnaSynt SelectFilter_b(F_HAAR);
	SubBandFilter *SB1D_h = new SubBandFilter(SelectFilter_b, NORM_L2);
	//SB1D_h->Border = I_CONT;
	UOWT *TV = new UOWT(SB1D_h);
	fltarray *TabBand_TV;
	if(FC_regul_TV)
		TV->init(TabBand_TV, Data.nx(), Data.ny(), Data.nz(),2);

// Print the threshold levels
	if(GVerbose)
	{
		float f = Lmin*1./(Lmax-Lmin);
		for(int iter=0;iter<Niter;iter++)
		{
		// Threshold level
			lambda = 1-iter*1.0/float(Niter-1);
			lambda = (pow((double)lambda,4.)+lambda*f)/(1.+f);
			lambda = Lmin+ lambda*(Lmax-Lmin);
			cerr<<lambda<<"\t";
		}
		cerr<<endl;
	}
	
// Main separation loop
	for(int iter=0;iter<Niter;iter++)
	{
	// Threshold level
		float f = Lmin*1./(Lmax-Lmin);
		lambda = 1-iter*1.0/float(Niter-1);
		lambda = (pow((double)lambda,4.)+lambda*f)/(1.+f);
		lambda = Lmin+ lambda*(Lmax-Lmin);
	
	// Transforms
		{
        // Wavelets W2D1D
          if(Use_WT2D1D)
          {
				if(GVerbose) cerr<<"WT2D1D... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
                
				Recons = Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT + Recons_DT;
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons + Recons_WT2D1D;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons - Recons_WT2D1D) + InvMask*Recons_WT2D1D;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons;
					else // Inpainting included in MCA
						Data = Orig - Mask*Recons;
				}
                
                // W2D1D initialisation
				if(iter==0)
					DataWT2D1D->alloc(Data.nx(), Data.ny(), Data.nz(), WT2D1D_Trans2D,  WT2D1D_NScale2d);
                
                // Forward Transform
				DataWT2D1D->transform(Data);
                TabBand_WT2D1D = DataWT2D1D->TabBand;
				
                // Filtering to be developped
				// if(WT2D1D_FilterType==FT_HARD || WT2D1D_FilterType==FT_SOFT) DataUW->threshold(TabBand_WT2D1D, SigmaNoise, lambda, WT2D1D_FilterType, force4sigma);
				// else if(UW_FilterType==FT_WIENER) DataUW->wiener(TabBand_UW, SigmaNoise, 3);
				// else cerr<<"Filtering method '"<<string_filter_type(UW_FilterType)<<"' not implemented yet";
                
                // Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataWT2D1D->recons(Recons_WT2D1D);
				else // MCA_type = MCA_RESIDUAL
				{
 					DataWT2D1D->recons(Recons);
					Recons_WT2D1D+=Recons;
				}
                
                // Borders
				if(use_min | use_max) set_border(Recons_WT2D1D,min_value, max_value);
                
                // Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_wt2d1d_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_UW);
				}
			}            

            
		// Meyer wavelets
			if(Use_MW)
			{
				if(GVerbose) cerr<<"MWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT) + InvMask*Recons_MW;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC);
				}
			// MWT initialisation
				if(iter==0)
					DataMW->init(Nscale_MW, Data.nx(), Data.ny(), Data.nz(), False, False);

			// Forward Transform
				DataMW->transform(Data, TabBand_MW);
				
			// Filtering
				if(MW_FilterType==FT_HARD || MW_FilterType==FT_SOFT)
					DataMW->threshold(TabBand_MW, SigmaNoise, lambda, MW_FilterType, !force4sigma);
				else if(MW_FilterType==FT_WIENER) DataMW->wiener(TabBand_MW, SigmaNoise, 3);
				else cerr<<"Filtering method '"<<string_filter_type(MW_FilterType)<<"' not implemented yet";

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataMW->recons(TabBand_MW,Recons_MW,False);
				else // MCA_type = MCA_RESIDUAL
				{
					DataMW->recons(TabBand_MW,Recons,False);
					Recons_MW+=Recons;
				}

			// Borders
				if(use_min | use_max) set_border(Recons_MW,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_mwt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_MW);
				}
			}
			
		// Undecimated wavelets
			if(Use_UW)
			{
				if(GVerbose) cerr<<"UWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				Recons = Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT + Recons_DT;
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons + Recons_UW;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons - Recons_UW) + InvMask*Recons_UW;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons;
					else // Inpainting included in MCA
						Data = Orig - Mask*Recons;
				}

			// UOWT initialisation
				if(iter==0)
					DataUW->init(TabBand_UW, Data.nx(), Data.ny(), Data.nz(), Nscale_UW);

			// Forward Transform
				DataUW->transform(Data,TabBand_UW,Nscale_UW);
				
			// Filtering
				if(UW_FilterType==FT_HARD || UW_FilterType==FT_SOFT)
					DataUW->threshold(TabBand_UW, SigmaNoise, lambda, UW_FilterType, force4sigma);
				else if(UW_FilterType==FT_WIENER) DataUW->wiener(TabBand_UW, SigmaNoise, 3);
				else cerr<<"Filtering method '"<<string_filter_type(UW_FilterType)<<"' not implemented yet";

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataUW->recons(TabBand_UW,Recons_UW,Nscale_UW);
				else // MCA_type = MCA_RESIDUAL
				{
					DataUW->recons(TabBand_UW,Recons,Nscale_UW);
					Recons_UW+=Recons;
				}

			// Borders
				if(use_min | use_max) set_border(Recons_UW,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_uwt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_UW);
				}
			}
			
		// Time Curvelets : 2D Curvelets + Time Wavelets
			if(Use_TC)
			{
				if(GVerbose) cerr<<"TCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT) + InvMask*Recons_TC;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}

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
				if(TC_FilterType==FT_HARD || TC_FilterType==FT_SOFT) DataTC->threshold(TabBand_TC, SigmaNoise, lambda, TC_FilterType);

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataTC->recons(TabBand_TC, Recons_TC);
				else // MCA_type = MCA_RESIDUAL
				{
					DataTC->recons(TabBand_TC, Recons);
					Recons_TC+=Recons;
				}
				
			// Borders
				if(use_min | use_max) set_border(Recons_TC,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_tcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_TC);
				}
			}
			
		// Stacked Curvelets : 2D Curvelets + Time constant
			if(Use_SC)
			{
				if(GVerbose) cerr<<"SCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT) + InvMask*Recons_SC;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}

			// Initialization
				if(iter==0)
				{
					DataSC->alloc_from_coarse(Nscale_SC, Data.nx(), Data.ny(), Data.nz(), SC_NDir2d, False, False, True);
					DataSC->set_threshold_coarse(SC_th_coarse);
				}

			// Forward Transform
				DataSC->transform(Data, TabBand_SC, (iter==0));

			// Filtering
				if(SC_FilterType==FT_HARD || SC_FilterType==FT_SOFT) DataSC->threshold(TabBand_SC, SigmaNoise, lambda, SC_FilterType);

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataSC->recons(TabBand_SC, Recons_SC);
				else // MCA_type = MCA_RESIDUAL
				{
					DataSC->recons(TabBand_SC, Recons);
					Recons_SC+=Recons;
				}
				
			// Borders
				if(use_min | use_max) set_border(Recons_SC,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_scur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_SC);
				}
			}
			
		// Beam Curvelets
			if(Use_BC)
			{
				if(GVerbose) cerr<<"BCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;

				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_RC + Recons_FC + Recons_DC + Recons_OT) + InvMask*Recons_BC;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}

			// Initialization
				if(iter==0)
				{
					type_linelet3d_WTtrans CurTrans = DEF_LIN3D_TRANS;
					type_wavelet3D TypeW3D = DEF_TYPE_W3D;
					DataBC->CurTrans=CurTrans;
					DataBC->BlockOverlap=BlockOver_BC;
					DataBC->set_stat(false);
					DataBC->set_no_recons(false);
					//DataBC->set_3sigma(Use3sig);
					DataBC->set_no_coarse(BC_no_coarse);
					DataBC->set_keep_energy(BC_keep_energy);
					DataBC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_BC, BlockSize_BC, TypeW3D, SB1D);
				}

			// Transform, filtering and reconstruction
				if(MCA_type == MCA_TOTAL)
					DataBC->filter(Data, Recons_BC, SigmaNoise, lambda, BC_FilterType, Alpha, WienerBS, force4sigma, UseCubeSigma, Name_Imag_Out);
				else // MCA_type = MCA_RESIDUAL
				{
					DataBC->filter(Data, Recons, SigmaNoise, lambda, BC_FilterType, Alpha, WienerBS, force4sigma, UseCubeSigma, Name_Imag_Out);
					Recons_BC+=Recons;
				}

			// Borders
				if(use_min | use_max) set_border(Recons_BC,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_bcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_BC);
				}
			}

		// A Trou wavelets
			if(Use_AW)
			{
				if(GVerbose) cerr<<"ATrouWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT) + InvMask*Recons_AW;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}
			
			// AWT initialisation
				if(iter==0)
				{
					DataAW->set_no_fine(AW_no_fine);
					DataAW->alloc(TabBand_AW, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
				}

			// Forward Transform
				DataAW->transform(Data, TabBand_AW, Nscale_AW);
				
			// Filtering
				if(AW_FilterType==FT_HARD || AW_FilterType==FT_SOFT)
					DataAW->threshold(TabBand_AW, SigmaNoise*lambda, AW_FilterType==FT_SOFT);
				else cerr<<"Filtering method '"<<string_filter_type(AW_FilterType)<<"' not implemented yet";

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataAW->recons(TabBand_AW, Recons_AW, Nscale_AW);
				else // MCA_type = MCA_RESIDUAL
				{
					DataAW->recons(TabBand_AW, Recons, Nscale_AW);
					Recons_AW+=Recons;
				}
				

			// Borders
				if(use_min | use_max) set_border(Recons_AW,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_awt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_AW);
				}
			}
			
		// Ridge Curvelets
			if(Use_RC)
			{
				if(GVerbose) cerr<<"RCur..."<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
			// Data update
				Recons = Recons_WT2D1D+ Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT + Recons_DT;
				if(MCA_type == MCA_TOTAL) // else MCA_RESIDUAL
					Recons -= Recons_RC;
				if(!Use_mask) // No inpainting, classical MCA
					Data = Orig - Recons;
				else // Inpainting included in MCA
				{
					Data = Orig - Mask*Recons;
					if(MCA_type == MCA_TOTAL) Data += InvMask*Recons_RC;
				}

				//sprintf(filename,"%s_rcur_data_%03d.fits", Name_Imag_Out, iter);
				//writefltarr(filename, Data);
			// Initialization
				if(iter==0)
				{
					type_wavelet3D TypeW3D = DEF_TYPE_W3D;//2
					type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1
					DataRC->CurTrans=CurTrans;
					DataRC->BlockOverlap=BlockOver_RC;
					DataRC->set_no_coarse(RC_no_coarse);
					DataRC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, TypeW3D, false);// Use3sig
				}

			// Forward Transform
				DataRC->transform(Data, TabBand_RC, (iter==0));

			// Filtering
				if(RC_FilterType==FT_HARD || RC_FilterType==FT_SOFT) DataRC->threshold(TabBand_RC, SigmaNoise, lambda, RC_FilterType, force4sigma, UseCubeSigma);
				else if(RC_FilterType==FT_WIENER) DataRC->wiener(TabBand_RC, SigmaNoise, 3, UseCubeSigma);
				else if(RC_FilterType==FT_FDR) DataRC->fdr(TabBand_RC, Alpha, SigmaNoise);

			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataRC->recons(TabBand_RC, Recons_RC);
				else // MCA_type = MCA_RESIDUAL
				{
					DataRC->recons(TabBand_RC, Recons);
					Recons_RC+=Recons;
				}
				
			// Borders
				if(use_min | use_max) set_border(Recons_RC,min_value, max_value);
					
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_rcur_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_RC);
				}
			}
			
		// Fast Curvelets
			if(Use_FC)
			{
				if(GVerbose) cerr<<"BFCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;
				
			// Data update
				Recons = Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT + Recons_DT;
				if(MCA_type == MCA_TOTAL) // else MCA_RESIDUAL
					Recons -= Recons_FC;
				if(!Use_mask) // No inpainting, classical MCA
					Data = Orig - Recons;
				else // Inpainting included in MCA
					Data = Orig - Mask*Recons + InvMask*Recons_FC;

			// Initialization
				if(iter==0)
				{
					DataFC->alloc_from_coarse(Data,P);
					
				// Noise calibration if not already done, and if useful here
					bool AllocTB = true;//(iter ==0);
					//if(!DataFC->isset_tabsigma())
					//	DataFC->get_norm_coeff(3, NULL);
				}
								
			// Transform, filtering and reconstruction
				if(MCA_type == MCA_TOTAL)
				{
					P.NSigma = lambda;
					DataFC->filter(Data,Recons_FC,P,Name_Imag_Out);
				}
				else // MCA_type = MCA_RESIDUAL
				{
					P.NSigma = lambda;
					DataFC->filter(Data,Recons,P,Name_Imag_Out);
					Recons_FC+=Recons;
				}

			// TV regularisation
				if(FC_regul_TV)
				{
					TV->transform(Recons_FC,TabBand_TV,2);// Nscale=2
					TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, force4sigma);
					TV->recons(TabBand_TV,Recons_FC,2);
				}

			// Borders
				if(use_min | use_max) set_border(Recons_FC,min_value, max_value);
				
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_bfct_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_FC);
				}
				
			}

		// Discrete cosine transform
			if(Use_DC)
			{
				if(GVerbose) cerr<<"DCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_OT) + InvMask*Recons_DC;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}
			
			// DCT initialisation
				if(iter==0)
				{
					DataDC->set_skip_order(skip_order_DC);
					DataDC->init(Data.nx(),Data.ny(),Data.nz(),BlockSize_DC,BlockOver_DC);
				}

			// Forward Transform
				DataDC->transform(Data,TabBand_DC);
				
			// Filtering
				if(DC_FilterType==FT_HARD || DC_FilterType==FT_SOFT)
					DataDC->threshold(TabBand_DC, SigmaNoise, lambda, DC_FilterType);
				else if(DC_FilterType==FT_WIENER) DataDC->wiener(TabBand_DC, SigmaNoise, 3);
				else if(DC_FilterType==FT_FDR) DataDC->fdr(TabBand_DC, Alpha, SigmaNoise);
				else cerr<<"Filtering method '"<<string_filter_type(DC_FilterType)<<"' not implemented yet";
				
			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataDC->recons(TabBand_DC,Recons_DC);
				else // MCA_type = MCA_RESIDUAL
				{
					DataDC->recons(TabBand_DC,Recons);
					Recons_DC+=Recons;
				}
				
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_dct_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_DC);
				}
			}

		// OT : Oriented Texture, Oriented DCT
			if(Use_OT)
			{
				if(GVerbose) cerr<<"Oriented DCT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC) + InvMask*Recons_OT;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons_WT2D1D - Recons_MW - Recons_AW - Recons_UW - Recons_TC - Recons_SC - Recons_BC - Recons_RC - Recons_FC - Recons_DC - Recons_OT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT);
				}
			
			// Oriented DCT initialisation
				if(iter==0)
					DataOT->init(Data.nx(),Data.ny(),Data.nz(),BlockSize_OT,BlockOver_OT);

			// Forward Transform
				DataOT->transform(Data,TabBand_OT);
				
			// Filtering
				DataOT->set_skip_order(skip_order_OT);
				if(OT_FilterType==FT_HARD || OT_FilterType==FT_SOFT)
					DataOT->threshold(TabBand_OT, SigmaNoise, lambda, OT_FilterType);
				else cerr<<"Filtering method '"<<string_filter_type(OT_FilterType)<<"' not implemented yet";
				
			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataOT->recons(TabBand_OT,Recons_OT);
				else // MCA_type = MCA_RESIDUAL
				{
					DataOT->recons(TabBand_OT,Recons);
					Recons_OT+=Recons;
				}
				
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_ot_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_OT);
				}
			}
			
		// DT : DCT 2D by frame
			if(Use_DT)
			{
				if(GVerbose) cerr<<"Frame DCT2D... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				Recons = Recons_WT2D1D + Recons_MW + Recons_AW + Recons_UW + Recons_TC + Recons_SC + Recons_BC + Recons_RC + Recons_FC + Recons_DC + Recons_OT + Recons_DT;
				if(MCA_type == MCA_TOTAL)
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons + Recons_DT;
					else // Inpainting included in MCA
						Data = Orig - Mask*(Recons - Recons_DT) + InvMask*Recons_DT;
				}
				else // MCA_type = MCA_RESIDUAL
				{
					if(!Use_mask) // No inpainting, classical MCA
						Data = Orig - Recons;
					else // Inpainting included in MCA
						Data = Orig - Mask*Recons;
				}
			
			// Frame DCT initialisation
				if(iter==0)
				{
					DataDT->alloc(Data.nx(), Data.ny(), Data.nz(), BlockSize_DT, BlockOver_DT);
					DataDT->set_skip_order(skip_order_DT);
				}

			// Forward Transform
				DataDT->transform(Data,TabBand_DT, true);
				
			// Filtering
//				DataDT->set_skip_order(skip_order_DT);
				if(DT_FilterType==FT_HARD || DT_FilterType==FT_SOFT)
					DataDT->threshold(TabBand_DT, SigmaNoise, lambda, DT_FilterType);
				else cerr<<"Filtering method '"<<string_filter_type(DT_FilterType)<<"' not implemented yet";
				
			// Reconstruction
				if(MCA_type == MCA_TOTAL)
					DataDT->recons(TabBand_DT,Recons_DT);
				else // MCA_type = MCA_RESIDUAL
				{
					DataDT->recons(TabBand_DT,Recons);
					Recons_DT+=Recons;
				}
				
			// Save output
				if(GVerbose)
				{
					sprintf(filename,"%s_dt_recons_%03d.fits", Name_Imag_Out, iter);
					writefltarr(filename, Recons_DT);
				}
			}
		}
	}
	
// Save output
	Recons = Recons_WT2D1D + Recons_FC + Recons_BC + Recons_RC + Recons_UW + Recons_DC + Recons_OT;
	sprintf(filename,"%s_recons.fits",Name_Imag_Out);
	writefltarr(filename, Recons);
	
// Free all arrays
	for (int i=0;i<short_argc;i++) delete [] short_argv[i];
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
	delete DataOT; TabBand_OT.~fltarray();
	delete SB1D;
	
	Data.~fltarray();
	Orig.~fltarray();
	Recons.~fltarray();
	Recons_FC.~fltarray();
	Recons_TC.~fltarray();
	Recons_BC.~fltarray();
	Recons_RC.~fltarray();
	Recons_AW.~fltarray();
	Recons_UW.~fltarray();
	Recons_MW.~fltarray();
	Recons_DC.~fltarray();
	Recons_OT.~fltarray();
	
	stop=clock();
//	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





