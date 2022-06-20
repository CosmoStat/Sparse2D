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
**    Date:  25/08/08
**    
**    File:  mr3d_rcur.cc
**
*******************************************************************************
**
**    DESCRIPTION  Meyer 3d wavelets program
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "MeyerWT3D.h"
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
bool Normalize=true;
bool NormalizeInv=false;
float SigmaNoise=-1;
float NSigma=3;
bool force4sigma=false;
int Niter=1;
int NbrScale = 4;
Bool Isotropic=False;
Bool NeedOdd=False;
Bool Extend=False;
filter_type FilterType = FT_HARD;

/***************************************/
 
/***************************************/
	
static void usage(char *argv[])
{
	// int i;
	fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-n NbrScale]\n");
	fprintf(OUTMAN, "         [-N NbrScale](alias)\n");
	manline();

	fprintf(OUTMAN, "         [-v] Verbose.\n");
	manline();

	fprintf(OUTMAN, "         [-I] Use Isotropic wavelets\n");
	manline();

	fprintf(OUTMAN, "         [-o] Force odd sizes\n");
	manline();

//	fprintf(OUTMAN, "         [-e] Extending the data by 4/3 (not working yet)\n");
//	manline();

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
	fprintf(OUTMAN, "                in the .mr output file or the statistics'");
	manline();

	fprintf(OUTMAN, "         [-g sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
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

	fprintf(OUTMAN, "         [--force4sigma] (not yet) Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at finest scale\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c;  

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"Ioen:N:vFRScg:s:f:i:")) != -1) 
	{
		switch (c) 
		{
			case 'I': 
				Isotropic = True;
				break;
			case 'o': 
				NeedOdd = True;
				break;
			case 'e': 
				Extend = True;
				break;
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{
					fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);
					exit(-1);
				}
				break;
			case 'N': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{
					fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);
					exit(-1);
				}
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
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
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
	bool tmp;
	GetLongOption(opts,"normalize",Normalize) && (Normalize = false);
	GetLongOption(opts,"stat",Extract_stat);
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"SigmaNoise",SigmaNoise);
	GetLongOption(opts,"NSigma",NSigma);
	GetLongOption(opts,"MRout",Output_mr);
	GetLongOption(opts,"noMR",No_mr);
	GetLongOption(opts,"hard",tmp)	&& (FilterType = FT_HARD);
	GetLongOption(opts,"soft",tmp)	&& (FilterType = FT_SOFT);
	GetLongOption(opts,"wiener",tmp)&& (FilterType = FT_WIENER);
	GetLongOption(opts,"force4sigma",force4sigma);
	GetLongOption(opts,"noise_calib",Noise_calib);

// Get short options
	filtinit(short_argc, short_argv);

// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
	
// Print parameters
	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(Compute_forward) cout << "Compute forward transform" << endl;
		if(Input_mr) cout << "Input MR file"<<endl;
		if(Compute_recons) cout << "Compute reconstruction" << endl;
		if(Output_mr) cout << "Output MR file"<<endl;
		if(Extend) cout << "Extending the data by 4/3"<<endl;
		if(Isotropic) cout << "Using Isotropic wavelets"<<endl;
		
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale << endl;  
		if(Normalize) 		cout << " Normalizing the transform" << endl;
		
		if(Noise_calib) 	cout << " Noise_calibration output" << endl;
		if(Extract_stat)	cout << " Statistics output" << endl;
		if(SigmaNoise>0)	cout << " Noise level = " <<  SigmaNoise    << endl;
		if(SigmaNoise>0)	cout << " NSigma Threshold = " <<  NSigma    << endl;
		if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
		if(force4sigma)		cout << " Force '4sigma lvl' at finest scale" << endl;
		cout << endl;  
	}

	fltarray Data;
	fltarray *TabBand;
	fltarray Recons;
	MEYER_WT3D *DataW = new MEYER_WT3D;
	char filename[64];
	
// Forward transform
	if(Compute_forward)
	{
	// Input data
		strcat(Name_Imag_In,".fits");
		fits_read_fltarr(Name_Imag_In, Data);
	}
	int Nx=Data.nx(),Ny=Data.ny(),Nz=Data.nz();

	if(Compute_forward)
	{
	// RCurvelet initialisation
		if(NeedOdd)
		{
			int NNx=2*(Nx/2)+1;
			int NNy=2*(Ny/2)+1;
			int NNz=2*(Nz/2)+1;
			if( Nx!=NNx || Ny!=NNy || Nz!=NNz )
			{
				fltarray data(NNx,NNy,NNz);
				int i,j,k;
				for(k=0;k<Nz;k++)
				for(j=0;j<Ny;j++)
				for(i=0;i<Nx;i++)
					data(i,j,k) = Data(i,j,k);
				for(k=Nz;k<NNz;k++)
				for(j=0 ;j<NNy;j++)
				for(i=0 ;i<NNx;i++)
					data(i,j,k) = Data(test_index_mirror(i,Nx), test_index_mirror(j,Ny), test_index_mirror(k,Nz));
				for(k=0 ;k<NNz;k++)
				for(j=Ny;j<NNy;j++)
				for(i=0 ;i<NNx;i++)
					data(i,j,k) = Data(test_index_mirror(i,Nx), test_index_mirror(j,Ny), test_index_mirror(k,Nz));
				for(k=0 ;k<NNz;k++)
				for(j=0 ;j<NNy;j++)
				for(i=Nx;i<NNx;i++)
					data(i,j,k) = Data(test_index_mirror(i,Nx), test_index_mirror(j,Ny), test_index_mirror(k,Nz));
				Data = data;
				writefltarr("out_data.fits", Data);
			}
		}
		DataW->init(NbrScale, Data.nx(), Data.ny(), Data.nz(), Extend, Isotropic, NeedOdd);
		
	// Forward Transform
		DataW->transform(Data, TabBand);
	}
	else // We load the transformed cube
	{
		strcat(Name_Imag_In,".mr");
		DataW->read(Name_Imag_In,TabBand,&NormalizeInv);
		
		// Inverse normalisation if normalized
		if(NormalizeInv) DataW->normalize_self(TabBand,1);
	}
	
// Tests
	DataW->test(TabBand);

// Statistic tools
	if(Noise_calib)
	{
		// variance estimation
		DataW->noise_calibration(TabBand, Name_Imag_Out);
	}
	if(Extract_stat) DataW->extract_stat(TabBand, Name_Imag_Out, Normalize);
	
// Thresholding
	if( SigmaNoise>=0 )
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
		{
			DataW->threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma);
			for(int i=1;i<Niter;i++)
			{
				DataW->recons(TabBand, Data);
				DataW->transform(Data, TabBand);
				DataW->threshold(TabBand, SigmaNoise, NSigma, FilterType, force4sigma);
			}
		}
		else if(FilterType==FT_WIENER) DataW->wiener(TabBand, SigmaNoise, 3);
	}
	
// Save the transform
	if(Output_mr)
	{//for(int i=0;i<NbrScale;i++) {sprintf(filename,"%s_B%i.fits",Name_Imag_Out,i);writefltarr(filename, TabBand[i]);}
	// Normalize the transform if not already
		if(Normalize) DataW->normalize_self(TabBand, 0);
	// Save the transform
		sprintf(filename,"%s.mr",Name_Imag_Out);
		DataW->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) DataW->normalize_self(TabBand, 1);
	}
	
	//DataW->values_at(TabBand,(char*)"coeflist.dat",Name_Imag_Out);
	
// Reconstruction and save output
	if(Compute_recons) 
	{
		//	for(int s=0;s<DataW->NbrScale3D;s++){sprintf(filename,"%s_B%d.fits",Name_Imag_Out,s);writefltarr(filename, TabBand[s]);}
		DataW->recons(TabBand,Recons);
		sprintf(filename,"%s.fits",Name_Imag_Out);
		if( Nx!=Data.nx() || Ny!=Data.ny() || Nz!=Data.nz() )
		{
				fltarray recons(Nx,Ny,Nz);
				for(int k=0;k<Nz;k++)
				for(int j=0;j<Ny;j++)
				for(int i=0;i<Nx;i++)
					recons(i,j,k) = Recons(i,j,k);
				Recons = recons;
		}
		writefltarr(filename, Recons);
	}
	
// Free memory	
	delete [] TabBand;
	delete DataW;
	
	for(int i=0;i<short_argc;i++)
		delete [] short_argv[i];
	delete [] short_argv;
	
	Data.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}





