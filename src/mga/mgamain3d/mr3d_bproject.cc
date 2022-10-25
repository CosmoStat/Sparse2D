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
**    Date:  20/07/08
**    
**    File:  mr3d_bcur.cc
**
*******************************************************************************
**
**    DESCRIPTION  BeamCurvelet program
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "BProject_3d2d.h"
#include "GetLongOptions.h"

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
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
bool Noise_calib=false;
bool Normalize=false;
bool NormalizeInv=false;
Bool BlockOverlap=False;
float SigmaNoise=-1;
float Alpha=0.01;
float NSigma=3;
bool Use3sig=false;
bool force4sigma=false;
bool UseCubeSigma=false;
bool no_coarse=false;
int Niter=1;
int NbrScale3D = 0;
int BlockSize = DEF_3DCUR_BLOCK_SIZE;
filter_type FilterType = FT_HARD;
type_wavelet3D TypeW3D = DEF_TYPE_W3D;
type_BProject_3d2d_WTtrans CurTrans = DEF_LIN3D_TRANS;
bool lowmem=false;
bool keep_energy=false;
int WienerBS=3;
int NbrScale = 3;
int NbrDir2d = 16;

/***************************************/

static void usage(char *argv[])
{
//					 12345678901234567890123456789012345678901234567890123456789012345678901234567890

    fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
	fprintf(OUTMAN, "         [-n NbrScale3D]\n");
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
	
	fprintf(OUTMAN, "         [--noMR] To prevent creating MR output\n");
	fprintf(OUTMAN, "                when there is no reconstruction either.");
	manline();
	
	fprintf(OUTMAN, "         [--stat] Statistic extraction\n");
	fprintf(OUTMAN, "                in file 'NameOut_stat.dat'");
	manline();

	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
	manline();

	fprintf(OUTMAN, "         [--normalize][-c] Normalize the transform\n");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] Set the coarse scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [--overlap][-O] BlockOverlap\n");
	manline();

	fprintf(OUTMAN, "         [--BS Blocksize][-B BlockSize]\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_3DCUR_BLOCK_SIZE);
	manline();

    fprintf(OUTMAN, "         [-b 3D_Beamlet_Type]\n");
    fprintf(OUTMAN, "            (1): Orthonormal BT,\n");
    fprintf(OUTMAN, "             2 : Pyramidal BT.(not implemented)\n");
    manline();

    fprintf(OUTMAN, "         [-w 3D_Wavelet_Type]\n");
    fprintf(OUTMAN, "             1 : A Trou,\n");
    fprintf(OUTMAN, "            (2): Meyer.\n");
	
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
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholding\n");
	manline();

	fprintf(OUTMAN, "         [--force4sigma] Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at finest scale\n");
	manline();

	fprintf(OUTMAN, "         [--use3sig] Use the cumulative at lvl 1-2.7e-3\n");
	fprintf(OUTMAN, "             to estimate the normalisation coefficients\n");
	fprintf(OUTMAN, "             instead of the variance\n");
	manline();

	fprintf(OUTMAN, "         [--usecubesigma] Use a table of normalisation \n");
	fprintf(OUTMAN, "             coefficients estimated in every point of each scale\n");
	manline();

	fprintf(OUTMAN, "         [--lowmem] Low memory usage\n");
	fprintf(OUTMAN, "             Doesn't evaluate the full transform,\n");
	fprintf(OUTMAN, "             use for simple filtering (block independent)\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;  
	int temp;
	
	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:vFRScOB:b:w:g:a:s:f:i:")) != -1) 
	{
		switch (c) 
		{
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale3D) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D scales parameter: %s\n", OptArg);exit(-1);}
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
			case 'c': 
				Normalize = true;
				break;
			case 'O': 
				BlockOverlap = True;
				break;
			case 'B': 
				if (sscanf(OptArg,"%d",&BlockSize) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 1D scales parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'b': 
				if (sscanf(OptArg,"%d",&temp) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D Beamlet Type: %s\n", OptArg);exit(-1);}
				CurTrans=(type_BProject_3d2d_WTtrans)temp;
				break;
			case 'w': 
				if (sscanf(OptArg,"%d",&temp) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D Wavelet Type: %s\n", OptArg);exit(-1);}
				TypeW3D=(type_wavelet3D)temp;
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
			case 'f': int plop;
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);} else FilterType = (filter_type) plop;
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
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

	it = opts.find("--BS");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>BlockSize; }
	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); Normalize=true; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--overlap");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	BlockOverlap=(Bool)plop; }
	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--MRin");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Input_mr; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
	it = opts.find("--noise_calib");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Noise_calib; }
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
	it = opts.find("--force4sigma");
	if(it!=opts.end()){ istringstream ss(it->second); force4sigma=true; }
	it = opts.find("--use3sig");
	if(it!=opts.end()){ istringstream ss(it->second); Use3sig=true; }
	it = opts.find("--usecubesigma");
	if(it!=opts.end()){ istringstream ss(it->second); UseCubeSigma=true; }
	it = opts.find("--no-coarse");
	if(it!=opts.end()){ istringstream ss(it->second); no_coarse=true; }
	it = opts.find("--lowmem");
	if(it!=opts.end()){ istringstream ss(it->second); lowmem=true; }
	it = opts.find("--energy");
	if(it!=opts.end()){ istringstream ss(it->second); keep_energy=true; }
	
// Get short options
	filtinit(short_argc, short_argv);

	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Input_mr) cout << "Input MR file"<<endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;
		
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale3D = " <<  NbrScale3D    << endl;  
		cout << "Wavelet_Type = " <<  (int)TypeW3D    << endl;  
		cout << "Beamlet_Type = " <<  (int)CurTrans    << endl;  
		cout << "BlockSize = " <<  BlockSize    << endl;
		if(BlockOverlap) cout << "BlockOverlaping activated" << endl;
		if(Normalize) cout << "Normalizing the transform" << endl;
		if(lowmem) cout << "Low memory filtering method" << endl;
		
		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;
		
		// Filtering
		if(SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(FilterType) <<endl;
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
			if(UseCubeSigma)	cout << " Using Local Sigma in bands " << endl;
			if(force4sigma)		cout << " Force '4sigma lvl' at finest scale" << endl;
			if(Use3sig)			cout << " Use 3sigma cumulative estimation" << endl;
			if(no_coarse)		cout << " Set the coarse scale to zero" << endl;
		}
		cout << endl;  
	}
	
	fltarray Data;
	Ifloat *** TabBand;
	fltarray Recons;
	fltarray Band;
	BProject_3d2d *DataC = new BProject_3d2d();
	char filename[64];
	
// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	if(!Input_mr) Compute_forward=true;
	if(Input_mr && Compute_forward){cerr<<"Cannot input a MR file AND make the forward transform."<<endl;exit(1);}

// Forward transform
	if(Compute_forward)
	{
	// Input data
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);

	// BCurvelet initialisation	
		DataC->set_BlockOverlap(BlockOverlap);
		DataC->set_LinTransf(CurTrans);
//cerr<<Data.nx()<<","<<Ny<<","<< Data.ny()<<","<< Data.nz()<<endl;
		DataC->alloc_from_coarse(NbrScale, Data.nx(), Data.ny(), Data.nz(), NbrDir2d, BlockOverlap);
		// build_64b => DataC->set_angle(0,8,-4);
		DataC->set_angle(0,8,-4);

	// Forward Transform
		DataC->transform(Data, TabBand, true);
	}
	else // We load the transformed cube
	{
/*		strcat(Name_Imag_In,".mr");
		DataC->CurTrans=CurTrans;
		DataC->BlockOverlap=BlockOverlap;
		DataC->set_stat(Extract_stat);
		DataC->set_3sigma(Use3sig);
		DataC->set_no_coarse(no_coarse);
		DataC->read(Name_Imag_In,TabBand,SB1D,&NormalizeInv);

		// Inverse normalisation if normalized
		if(NormalizeInv) DataC->normalize_self(TabBand,1, UseCubeSigma);
*/	}


// Statistic tools
	if(Noise_calib)
	{
		// variance estimation
//		DataC->noise_calibration(TabBand, Name_Imag_Out);

		// exact N_Sigma threshold equivalent
//		int N_Sigma=3;
//		DataC->calib_noise_nsig(TabBand, N_Sigma, Name_Imag_Out);
	}
//	if(Extract_stat) DataC->extract_stat(TabBand, Name_Imag_Out, UseCubeSigma);

// Thresholding
/*	if( SigmaNoise>=0 )
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
		{
			DataC->threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma, UseCubeSigma);
			for(int i=1;i<Niter;i++)
			{
				DataC->recons(TabBand, Data);
				DataC->transform(Data, TabBand, false);
				DataC->threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma, UseCubeSigma);
			}
		}
		else if(FilterType==FT_WIENER) DataC->wiener(TabBand, SigmaNoise, 3, UseCubeSigma);
		else if(FilterType==FT_FDR) DataC->fdr(TabBand, Alpha, SigmaNoise);
	}
*/

// Save the transform
	if(Output_mr)
	{
/*	// Normalize the transform if needed
		if(Normalize) DataC->normalize_self(TabBand, 0, UseCubeSigma);
	// Save the transform
		sprintf(filename,"%s.mr",Name_Imag_Out);
		DataC->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) DataC->normalize_self(TabBand, 1, UseCubeSigma);
*/	
	}
		
//sprintf(filename,"%s_trans.fits",Name_Imag_Out);
//writefltarr(filename, TabBand);

// Reconstruction & Save output
	if(Compute_recons) 
	{
		DataC->recons(TabBand, Recons);
		sprintf(filename,"%s_recons.fits",Name_Imag_Out);
		writefltarr(filename, Recons);
	}

	delete DataC;
	
	for(int i=0;i<DataC->nbr_block();i++)
	{
		for(int s=0;s<NbrScale;s++)
			delete [] TabBand[i][s];
		delete [] TabBand[i];
	}
	delete [] TabBand;
		
	for(int i=0;i<short_argc;i++)
		delete [] short_argv[i];
	delete [] short_argv;
	Data.~fltarray();
	Recons.~fltarray();
	stop=clock();
	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





