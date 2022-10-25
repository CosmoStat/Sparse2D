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
**    File:  mr3d_rcur.cc
**
*******************************************************************************
**
**    DESCRIPTION  RidCurvelet program
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "RCurvelet3D.h"
#include "GetLongOptions.h"

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

Bool Verbose=False;
bool Output_mr=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
bool Noise_calib=false;
bool Normalize=true;
bool NormalizeInv=false;
float BlockOverlap=0;
float SigmaNoise=-1;
float Alpha=0.01;
float NSigma=3;
bool Use3sig=false;
bool force4sigma=false;
bool UseCubeSigma=false;
bool no_coarse=false;
int Niter=1;
int NbrScale3D = 0; int NbrScale1D = 0;
int BlockSize = DEF_3DCUR_BLOCK_SIZE;
filter_type FilterType = FT_HARD;
type_wavelet3D TypeW3D = DEF_TYPE_W3D;//2
type_ridgelet3d_WTtrans CurTrans = DEF_RID3D_TRANS;//1

/***************************************/

static void usage(char *argv[])
{
	// int i;
	fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-n NbrScale3D]\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_3DCUR_NBR_SCALE);
	manline();

	fprintf(OUTMAN, "         [-N NbrScale1D]\n");
	fprintf(OUTMAN, "             default = floor(log2(BlockSize))\n");
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
	
	fprintf(OUTMAN, "         [--stat] Statistic extraction\n");
	fprintf(OUTMAN, "                in file 'NameOut_stat.dat'");
	manline();

	fprintf(OUTMAN, "         [--noise_calib][-S] Noise calibration\n");
	fprintf(OUTMAN, "                in file 'NameOut_noise_calib.dat'");
	manline();

	fprintf(OUTMAN, "         [--normalize][-c] Do NOT normalize the transform\n");
	fprintf(OUTMAN, "                in the .mr output file'");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] Set the coarse scale to zero\n");
	manline();

	fprintf(OUTMAN, "         [--overlap f][-O f]\n");
	fprintf(OUTMAN, "             BlockOverlap with overlapping size f*BlockSize\n");
	manline();

	fprintf(OUTMAN, "         [--BS Blocksize][-B BlockSize]\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_3DCUR_BLOCK_SIZE);
	manline();

	fprintf(OUTMAN, "         [-r 3D_Ridgelet_Type]\n");
	fprintf(OUTMAN, "             1:Orthonormal RT,(not implemented)\n");
	fprintf(OUTMAN, "             2:Pyramidal RT.(default)\n");
	manline();

    fprintf(OUTMAN, "         [-w 3D_Wavelet_Type]\n");
    fprintf(OUTMAN, "             1 : A Trou,\n");
    fprintf(OUTMAN, "            (2): Meyer.\n");
    fprintf(OUTMAN, "             3 : Poisson noise adapted Meyer.\n");
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

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;  
	int temp;
	
	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:N:vFRScOB:r:w:g:a:s:f:i:")) != -1) 
	{
		switch (c) 
		{
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale3D) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D scales parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'N': 
				if (sscanf(OptArg,"%d",&NbrScale1D) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 1D scales parameter: %s\n", OptArg);exit(-1);}
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
				Normalize = false;
				break;
			case 'O': 
				if (sscanf(OptArg,"%f",&BlockOverlap) != 1) 
				{fprintf(OUTMAN, "Error: bad Overlapping factor: %s\n", OptArg);exit(-1);}
				break;
			case 'B': 
				if (sscanf(OptArg,"%d",&BlockSize) != 1) 
				{fprintf(OUTMAN, "Error: bad Block size: %s\n", OptArg);exit(-1);}
				break;
			case 'r': 
				if (sscanf(OptArg,"%d",&temp) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D Ridgelet Type: %s\n", OptArg);exit(-1);}
				CurTrans=(type_ridgelet3d_WTtrans)temp;
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
			case 's': 
				if (sscanf(OptArg,"%f",&NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'a': 
				if (sscanf(OptArg,"%f",&Alpha) != 1) 
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
	map<string,string>::iterator it;

	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--BS");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>BlockSize; }
	it = opts.find("--overlap");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>BlockOverlap; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); Normalize=false; }
	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
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
	
// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;
		
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale1D = " <<  NbrScale1D    << endl;  
		cout << "NbrScale3D = " <<  NbrScale3D    << endl;  
		cout << "Wavelet_Type = " <<  (int)TypeW3D    << endl;  
		cout << "Ridgelet_Type = " <<  (int)CurTrans    << endl;  
		cout << "BlockSize = " <<  BlockSize    << endl;
		if(BlockOverlap>0) cout << "BlockOverlaping activated, factor=" << BlockOverlap  << endl;
		if(!Normalize) cout << "NOT normalizing the transform output" << endl;
		
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

#if USE_OMP_RC
	omp_set_nested(1);
#endif
	
	fltarray Data;
	fltarray *TabBand;
	fltarray Recons;
	fltarray Band;
	RCurvelet3D *DataC = new RCurvelet3D();
	char* filename;
	
// Forward transform
	if(Compute_forward)
	{
	// Input data
		Name_Imag_In = add_fits(Name_Imag_In);
		fits_read_fltarr(Name_Imag_In, Data);

	// RCurvelet initialisation
		DataC->CurTrans=CurTrans;
		DataC->BlockOverlap=BlockOverlap;
		DataC->set_no_coarse(no_coarse);
		DataC->init(Data.nx(), Data.ny(), Data.nz(), NbrScale3D, NbrScale1D, BlockSize, TypeW3D, Use3sig);

	// Forward Transform
		DataC->transform(Data, TabBand, true);
	}
	else // We load the transformed cube
	{
		Name_Imag_In = add_mr(Name_Imag_In);
		DataC->set_no_coarse(no_coarse);
		DataC->read(Name_Imag_In,TabBand,&NormalizeInv);
		
		// Inverse normalisation if normalized
		if(NormalizeInv) DataC->normalize_self(TabBand,1, UseCubeSigma);
	}
	
// Tests
	DataC->temp(TabBand);
	
// Get Band
//	DataC->get_band(0,0,TabBand, Band);
//	sprintf(filename,"%s_band.fits",Name_Imag_Out);
//	writefltarr(filename, Band);

// Statistic tools
	if(Noise_calib)
	{
		// variance estimation
		DataC->noise_calibration(TabBand, Name_Imag_Out);
		
		// exact N_Sigma threshold equivalent
		int N_Sigma=3;
		DataC->calib_noise_nsig(TabBand, N_Sigma, Name_Imag_Out);
	}
	if(Extract_stat) DataC->extract_stat(TabBand, Name_Imag_Out, UseCubeSigma, Normalize);
	
// Thresholding
	if( SigmaNoise>=0 )
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
		{
			if(Niter==1)
				DataC->threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma, UseCubeSigma);
			else
			{
				float lambda=5*SigmaNoise;
				

// support constraint (JL), equivatent to CFA with 1 transform
				float delta = lambda/(Niter-1);
				DataC->threshold(TabBand, SigmaNoise, NSigma, FT_HARD, force4sigma, UseCubeSigma);
				fltarray *TabBandO = new fltarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					TabBandO[s3]=TabBand[s3];
				DataC->threshold(TabBand, lambda, NSigma, FT_HARD, force4sigma, UseCubeSigma);
				
				for(int i=1;i<Niter;i++)
				{
					lambda-=delta;
					cerr<<"th="<<lambda<<endl;
					DataC->recons(TabBand, Data);
					DataC->transform(Data, TabBand, false);
					DataC->update_significant(TabBand, TabBandO, SigmaNoise, NSigma, force4sigma, UseCubeSigma);
					DataC->threshold(TabBand, lambda, NSigma, FT_SOFT, force4sigma, UseCubeSigma);
				}

	
				
/*
// stomp
				float delta = lambda/(Niter-1);

				fltarray *TabBandO = new fltarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					TabBandO[s3]=TabBand[s3];
				
				fltarray *TabErr = new fltarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					TabErr[s3]=TabBand[s3];
				
				intarray *support = new intarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					support[s3].resize(TabBand[s3].nx(),TabBand[s3].ny(),TabBand[s3].nz());
				
				fltarray xn(Data.nx(),Data.ny(),Data.nz());
				
				DataC->threshold(TabErr, lambda, NSigma, FT_HARD, force4sigma, UseCubeSigma);
				
				for(int i=1;i<Niter;i++)
				{
					// select significant in residual
					DataC->support_significant(support, TabBandO, SigmaNoise, NSigma, force4sigma, UseCubeSigma);
					DataC->select_significant(TabBand, support, TabBandO);
					
					DataC->recons(TabBand, xn);
					xn=Data-xn; // Error
					DataC->transform(xn, TabBand, false);
				}
*/
	
		
/*
// energy recovery : not much energy in detected coefficients for the residual
				float delta = lambda/(Niter-1);
				DataC->threshold(TabBand, SigmaNoise, NSigma, FT_HARD, force4sigma, UseCubeSigma);
				fltarray *TabErr = new fltarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					TabErr[s3]=TabBand[s3];
				fltarray *TabBandO = new fltarray[NbrScale3D];
				for(int s3=0;s3<NbrScale3D;s3++)
					TabBandO[s3]=TabBand[s3];
				
				DataC->threshold(TabBand, lambda, NSigma, FT_HARD, force4sigma, UseCubeSigma);
				
				fltarray xn(Data.nx(),Data.ny(),Data.nz());
				
				for(int i=1;i<Niter;i++)
				{
					lambda-=delta;
					cerr<<"th="<<lambda<<endl;
					DataC->recons(TabBand, xn);
					xn=Data-xn; // Error
					DataC->transform(xn, TabErr, false);
					DataC->select_significant(TabErr, TabBandO, SigmaNoise, NSigma, force4sigma, UseCubeSigma);
					for(int s3=0;s3<NbrScale3D;s3++)
						TabBand[s3]+=TabErr[s3];
					DataC->threshold(TabBand, lambda, NSigma, FT_SOFT, force4sigma, UseCubeSigma);
				}
*/



			}
		}
		else if(FilterType==FT_WIENER) DataC->wiener(TabBand, SigmaNoise, 3, UseCubeSigma);
		else if(FilterType==FT_FDR) DataC->fdr(TabBand, Alpha, SigmaNoise);
	}

// Get Band
//	DataC->get_band(0,0,TabBand, Band);
//	sprintf(filename,"%s_bandb.fits",Name_Imag_Out);
//	writefltarr(filename, Band);

// Save the transform
	if(Output_mr)
	{
	// Normalize the transform if not already
		if(Normalize) DataC->normalize_self(TabBand, 0, UseCubeSigma);
	// Save the transform
		filename = add_mr(Name_Imag_Out);
		DataC->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) DataC->normalize_self(TabBand, 1, UseCubeSigma);
	}

	//DataC->values_at(TabBand,(char*)"coeflist.dat",Name_Imag_Out);

// Reconstruction and save output
	if(Compute_recons) 
	{
		//	for(int s=0;s<DataC->NbrScale3D;s++){sprintf(filename,"%s_B%d.fits",Name_Imag_Out,s);writefltarr(filename, TabBand[s]);}
		DataC->recons(TabBand, Recons);
		filename = add_fits(Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	
// Free memory	
	delete [] TabBand;
	delete DataC;
	
	for(int i=0;i<short_argc;i++)
		delete [] short_argv[i];
	delete [] short_argv;
	
	Data.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}


