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
**    File:  mr3d_cfa_poison.cc
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
#include "uowt.h"
#include "RCurvelet3D.h"
#include "BCurvelet3D.h"
#include "MR3D_Obj.h"
#include "PoissonAtrou3D.h"

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
float SigmaNoise=1;

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

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"v?")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose = True;
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
	
	// Atrou wavelets				// Meyer wavelets  		
	bool Use_AW			= false; 	bool Use_MW		 =false; 
	int Nscale_AW		= 4; 		int Nscale_MW  	 =4;	
	
	// BeamCurvelets				// RidCurvelets				
	bool Use_BC			= true; 	bool Use_RC			=false;	
	int Nscale_BC		= 3;		int Nscale_RC		=4;		
	int BlockSize_BC	= 17;		int BlockSize_RC	=17;	
	float BlockOver_BC	= 0; 		Bool BlockOver_RC	=True;	
	bool BC_no_coarse	= false; 	bool RC_no_coarse	=false;	
	type_pmeyer3D Type_PMW3D = PMW3D_SUM; //PMW3D_COHERENT; //
	
// Levels
	float Nsig_orig=3.;
	float th_far = SigmaNoise/5.;
	use_max=false;	float max_value = 255.;
	use_min=true;	float min_value = 0.;
	
// **********************************************
// **********************************************
// 				CFA initialization
// **********************************************
// **********************************************
	
	float pos = Nsig_orig*SigmaNoise/2.;		// positive detection in original

// Print parameters
	//if (Verbose)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << " Iterations parameters : "<<Niter<<" in ["<<Lmax<<","<<Lmin<<"]"<<endl;
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
	BCurvelet3D *DataBC = new BCurvelet3D;
	RCurvelet3D *DataRC = new RCurvelet3D;
	ATROUS_3D_WT *DataAW = new ATROUS_3D_WT;
	MEYER_WT3D *DataMW = new MEYER_WT3D;
// 3D Poisson VST transforms declaration
	BCurvelet3D *DataBCP = new BCurvelet3D;
	RCurvelet3D *DataRCP = new RCurvelet3D;
	POISSONATROUS_3D_WT *DataAWP = new POISSONATROUS_3D_WT;
	POISSON_MWT3D *DataMWP = new POISSON_MWT3D;
	
	fltarray *TabBand_BC,	*TabBand_BC_orig,  *TabBand_BCP_orig;
	fltarray *TabBand_RC,	*TabBand_RC_orig,  *TabBand_RCP_orig;
	fltarray *TabBand_AW,	*TabBand_AW_orig,  *TabBand_AWP_orig;
	fltarray *TabBand_MW,	*TabBand_MW_orig,  *TabBand_MWP_orig;
	
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
	float posP = Nsig_orig/2.;
	type_wavelet3D TypeW3D = W3D_MEYER;
	
	if(Use_AW)
	{
		if(Verbose) cerr<<"Init AWT"<<endl;
	// Detection part : Poisson noise
		// initialisation
		DataAWP->alloc(TabBand_AWP_orig, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
		// Forward Transform
		DataAWP->transform(Data, TabBand_AWP_orig, Nscale_AW);
		// Detection
		DataAWP->threshold(TabBand_AWP_orig, Nsig_orig, FT_HARD);
		DataAWP->normalize_self(TabBand_AWP_orig, false);
		
	// Approximation
		// initialisation
		DataAW->alloc(TabBand_AW, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
		DataAW->alloc(TabBand_AW_orig, Data.nx(), Data.ny(), Data.nz(),Nscale_AW);
		// Forward Transform
		DataAW->transform(Data, TabBand_AW_orig, Nscale_AW);
		// Filtering
		DataAW->normalize_self(TabBand_AW_orig, false);
		for(int s=0;s<Nscale_AW-1;s++)
		for(int i=0;i<Data.nx();i++)
		for(int j=0;j<Data.ny();j++)
		for(int k=0;k<Data.nz();k++)
			if( abs( (TabBand_AWP_orig[s])(i,j,k) ) < posP ) // !=0 but for float data
				(TabBand_AW_orig[s])(i,j,k) = 0;

		if(GVerbose)
		{
		// Reconstruction
			DataAW->normalize_self(TabBand_AW_orig, true);
			DataAW->recons(TabBand_AW_orig, Recons, Nscale_AW);
			DataAW->normalize_self(TabBand_AW_orig, false);

		// Save output
			sprintf(filename,"%s_awt_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
		delete [] TabBand_AWP_orig;
		delete DataAWP;
	}

	if(Use_MW)
	{
		if(Verbose) cerr<<"Init MWT"<<endl;
	// Detection part : Poisson noise
		// initialisation
		DataMWP->init(Nscale_MW, Data.nx(), Data.ny(), Data.nz(), False, False);
		// Forward Transform
		DataMWP->transform(Data, TabBand_MWP_orig);
		// Detection
		DataMWP->threshold(TabBand_MWP_orig, Nsig_orig, FT_HARD);
		DataMWP->Type_PMW3D=Type_PMW3D;
		DataMWP->normalize_self(TabBand_MWP_orig, false);
		
	// Approximation
		// initialisation
		DataMW->init(Nscale_MW, Data.nx(), Data.ny(), Data.nz(), False, False);
		// Forward Transform
		DataMW->transform(Data, TabBand_MW_orig);
		// Filtering
		DataMW->normalize_self(TabBand_MW_orig, false);
		for(int s=0;s<Nscale_MW-1;s++)
		for(int i=0;i<TabBand_MWP_orig[s].nx();i++)
		for(int j=0;j<TabBand_MWP_orig[s].ny();j++)
		for(int k=0;k<TabBand_MWP_orig[s].nz();k++)
			if( abs( (TabBand_MWP_orig[s])(i,j,k) ) < posP ) // !=0 but for float data
				(TabBand_MW_orig[s])(i,j,k) = 0;

		if(GVerbose)
		{
		// Reconstruction
			DataMWP->normalize_self(TabBand_MWP_orig, true);
			DataMWP->recons(TabBand_MWP_orig, Recons);
		// Save output
			sprintf(filename,"%s_mwtP_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
			
		// Reconstruction
			DataMW->normalize_self(TabBand_MW_orig, true);
			DataMW->recons(TabBand_MW_orig, Recons);
			DataMW->normalize_self(TabBand_MW_orig, false);

		// Save output
			sprintf(filename,"%s_mwt_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
		delete [] TabBand_MWP_orig;
		delete DataMWP;
	}

	if(Use_BC)
	{
		if(Verbose) cerr<<"Init BC"<<endl;
		type_linelet3d_WTtrans CurTrans = DEF_LIN3D_TRANS;
	
	// Detection part : Poisson noise
		// Initialization
		DataBCP->CurTrans=CurTrans;
		DataBCP->Type_PMW3D=Type_PMW3D;
		DataBCP->BlockOverlap=BlockOver_BC;
		DataBCP->set_stat(false);
		DataBCP->set_no_recons(false);
		DataBCP->set_no_coarse(BC_no_coarse);
		DataBCP->init(Data.nx(), Data.ny(), Data.nz(), Nscale_BC, BlockSize_BC, W3D_PMEYER, SB1D);
		// Forward Transform
		DataBCP->transform(Data, TabBand_BCP_orig, true);
		// Filtering
		DataBCP->threshold(TabBand_BCP_orig, 1., Nsig_orig, FT_HARD, false, false);
		DataBCP->normalize_self(TabBand_BCP_orig, false, false);
		
	// Approximation
		// Initialization
		DataBC->CurTrans=CurTrans;
		DataBC->BlockOverlap=BlockOver_BC;
		DataBC->set_stat(false);
		DataBC->set_no_recons(false);
		DataBC->set_3sigma(false);
		DataBC->set_no_coarse(BC_no_coarse);
		DataBC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_BC, BlockSize_BC, TypeW3D, SB1D);
		// Forward Transform
		DataBC->transform(Data, TabBand_BC_orig, true);
		// Filtering
		DataBC->normalize_self(TabBand_BC_orig, false, false);
		for(int s=0;s<Nscale_BC-1;s++)
		for(int i=0;i<TabBand_BCP_orig[s].nx();i++)
		for(int j=0;j<TabBand_BCP_orig[s].ny();j++)
		for(int k=0;k<TabBand_BCP_orig[s].nz();k++)
			if( abs( (TabBand_BCP_orig[s])(i,j,k) ) < posP ) // !=0 but for float data
				(TabBand_BC_orig[s])(i,j,k) = 0.;
		
		if(GVerbose)
		{
		// Reconstruction
			DataBCP->normalize_self(TabBand_BCP_orig, true, false);
			DataBCP->recons(TabBand_BCP_orig, Recons);
		// Save output
			sprintf(filename,"%s_bcurP_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
			
		// Reconstruction
			DataBC->normalize_self(TabBand_BC_orig, true, false);
			DataBC->recons(TabBand_BC_orig, Recons);
			DataBC->normalize_self(TabBand_BC_orig, false, false);
		// Save output
			sprintf(filename,"%s_bcur_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
		delete [] TabBand_BCP_orig;
		delete DataBCP;
	}

	if(Use_RC)
	{
		if(Verbose) cerr<<"Init RC"<<endl;
	// Detection part : Poisson noise
		// Initialization
		type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1
		DataRCP->CurTrans=CurTrans;
		DataRCP->BlockOverlap=BlockOver_RC;
		DataRCP->set_no_coarse(RC_no_coarse);
		DataRCP->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, W3D_PMEYER, false);
		// Forward Transform
		DataRCP->transform(Data, TabBand_RCP_orig, true);
		// Filtering
		DataRCP->threshold(TabBand_RC_orig, 1., Nsig_orig, FT_HARD, false, false);
		DataRCP->normalize_self(TabBand_RCP_orig, false, false);
	
	// Approximation
		// initialisation
		DataRC->CurTrans=CurTrans;
		DataRC->BlockOverlap=BlockOver_RC;
		DataRC->set_no_coarse(RC_no_coarse);
		DataRC->init(Data.nx(), Data.ny(), Data.nz(), Nscale_RC, 0, BlockSize_RC, TypeW3D, false);
		// Forward Transform
		DataRC->transform(Data, TabBand_RC_orig, true);
		// Filtering
		DataRC->normalize_self(TabBand_RC_orig, false, false);
		for(int s=0;s<Nscale_RC-1;s++)
		for(int i=0;i<TabBand_RC_orig[s].nx();i++)
		for(int j=0;j<TabBand_RC_orig[s].ny();j++)
		for(int k=0;k<TabBand_RC_orig[s].nz();k++)
			if( abs( (TabBand_RCP_orig[s])(i,j,k) ) < posP ) // !=0 but for float data
				(TabBand_RC_orig[s])(i,j,k) = 0.;

		if(GVerbose)
		{
		// Reconstruction
			DataRC->normalize_self(TabBand_RC_orig, true, false);
			DataRC->recons(TabBand_RC_orig, Recons);
			DataRC->normalize_self(TabBand_RC_orig, false, false);

		// Save output
			sprintf(filename,"%s_rcur_orig_recons.fits", Name_Imag_Out);
			writefltarr(filename, Recons);
		}
		delete [] TabBand_RCP_orig;
		delete DataRCP;
	}

	
	Recons=fltarray(Data.nx(),Data.ny(),Data.nz());
	Recons.init(0.0);
	
	if(Verbose) cerr<<"End Init"<<endl;
	
// Main denoising loop
	float lambda = Lmax;
	float delta = (Lmax-Lmin)/float(Niter-1);
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
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, false);
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

				if(Use_MW)
				{
					if(GVerbose) cerr<<"ATrouWT... "<<"iter "<<iter<<", lambda="<<lambda<<endl;

				// Forward Transform
					DataMW->transform(Recons, TabBand_MW);

				// Filtering
					if(!regul_TV)
						DataMW->threshold(TabBand_MW, SigmaNoise*lambda, true); // true for FT_SOFT
					DataMW->normalize_self(TabBand_MW);

				// Update the coefficients
float cnt=0;
float cntt=0;
					for(int s=0;s<Nscale_MW-1;s++)
					for(int i=0;i<TabBand_MW_orig[s].nx();i++)
					for(int j=0;j<TabBand_MW_orig[s].ny();j++)
					for(int k=0;k<TabBand_MW_orig[s].nz();k++)
						if( abs( (TabBand_MW_orig[s])(i,j,k) ) > pos ) // !=0 but for float data
						{
							if( abs( (TabBand_MW[s])(i,j,k)-(TabBand_MW_orig[s])(i,j,k) ) > th_far )
							{
								cnt++;
								(TabBand_MW[s])(i,j,k) = (TabBand_MW_orig[s])(i,j,k);
							}
							cntt++;
						}
//cerr<<"MW maj : "<<cnt<<"/"<<cntt<<endl;
					TabBand_MW[Nscale_MW-1] = TabBand_MW_orig[Nscale_MW-1];
					
				// Reconstruction
					DataMW->normalize_self(TabBand_MW, true);
					DataMW->recons(TabBand_MW, Recons);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, false);
						TV->recons(TabBand_TV,Recons,2);
					}

				// Borders
					set_border(Recons,min_value, max_value);
					
				// Save output
					if(GVerbose)
					{
						sprintf(filename,"%s_mwt_recons_%03d.fits", Name_Imag_Out, iter);
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
						DataBC->threshold(TabBand_BC, SigmaNoise, lambda, FT_SOFT, false, false);
					DataBC->normalize_self(TabBand_BC, 0, false);

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
					DataBC->normalize_self(TabBand_BC, 1, false);
					DataBC->recons(TabBand_BC, Recons);
					
				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, false);
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
						DataRC->threshold(TabBand_RC, SigmaNoise, lambda, FT_SOFT, false, false);
					DataRC->normalize_self(TabBand_RC, false, false);

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
					DataRC->normalize_self(TabBand_RC, true, false);
					DataRC->recons(TabBand_RC, Recons);

				// Borders
					set_border(Recons,min_value, max_value);
					
				// TV regularisation
					if(regul_TV)
					{
						TV->transform(Recons,TabBand_TV,2);// Nscale=2
						TV->threshold(TabBand_TV, SigmaNoise, lambda, FT_SOFT, false);
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
			}
		}
	}
	
// Last TV regularisation
	if(regul_TV_last)
	{
		TV->transform(Recons,TabBand_TV,2);// Nscale=2
		TV->threshold(TabBand_TV, SigmaNoise, Nsig_last_TV, FT_HARD, false);
		TV->recons(TabBand_TV,Recons,2);
	}

// Save output
	sprintf(filename,"%s_recons.fits", Name_Imag_Out);
	writefltarr(filename, Recons);

	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	
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
	delete DataAW;
	if(Use_AW) 
	{
		delete [] TabBand_AW;
		delete [] TabBand_AW_orig;
	}
	delete DataMW;
	if(Use_MW) 
	{
		delete [] TabBand_MW;
		delete [] TabBand_MW_orig;
	}
	
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





