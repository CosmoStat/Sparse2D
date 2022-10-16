
#include <time.h>

#include "PoissonAtrou3D.h"
#include "IM_IO.h"
#include "GetLongOptions.h"
#include "MGA_Inc.h"

char* Name_Imag_In; /* input file image */
char* Name_Imag_Out; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

#define DEF_NBR_ATROU_SCALE3D 4
/********************************/
//		Input variables			//
/********************************/

Bool Verbose=False;
bool Output_mr=false;
bool No_mr=false;
bool Compute_forward=false;
bool Compute_recons=false;
bool Extract_stat=false;
float SigmaNoise=0;
float NSigma=3;
int Niter=1;
int NbrScale3D = 3;
bool force4sigma=false;
bool Normalize=true;
bool NormalizeInv=false;
bool del_single=true;
bool no_fine=false;

filter_type FilterType = FT_HARD;
type_sb_filter wavelet_type = F_MALLAT_7_9;//F_HAAR

/***************************************/
/*
// list of filters
DEF_SB_FILTER F_MALLAT_7_9
enum type_sb_filter {	SB_UNKNOWN, 
						F_MALLAT_7_9, 
						F_DAUBE_4, 
						F_BI2HAAR, 
						F_BI4HAAR,
						F_ODEGARD_7_9,
						F_5_3,
						F_LEMARIE_1,
						F_LEMARIE_3,
						F_LEMARIE_5,
						F_USER,
						F_HAAR,
						F_3_5,
						F_4_4,
						F_5_3_DIV,
						F_MALLAT_9_7};
*/
static void usage(char *argv[])
{
	// int i;
	fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-n NbrScale3D]\n");
	fprintf(OUTMAN, "             default = %d\n",DEF_NBR_ATROU_SCALE3D);
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
/*	
	fprintf(OUTMAN, "         [--stat] Statistic extraction\n");
	fprintf(OUTMAN, "                in file 'NameOut_stat.dat'");
	manline();
*/
    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    	fprintf(OUTMAN, "              %d: %s \n",i, StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n\n", StringSBFilter ((type_sb_filter) wavelet_type));
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
	manline();

	fprintf(OUTMAN, "         [--force4sigma] Force thresholding at 4 sigma\n");
	fprintf(OUTMAN, "             at the finest scale\n");
	manline();

    exit(-1);
}

/*********************************************************************/

static void filtinit(int argc, char *argv[])
{
	int c;  
	int plop;
	
	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"n:N:vFRg:s:f:i:T:cP")) != -1) 
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
			case 'T': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad Filter Type: %s\n", OptArg);exit(-1);} else wavelet_type = (type_sb_filter) plop;
				break;
			case 'F': 
				Compute_forward = true;
				break;
			case 'R': 
				Compute_recons = true;
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 's': 
				if (sscanf(OptArg,"%f",&NSigma) != 1) 
				{fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);exit(-1);}
				break;
			case 'f': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad Filter Type: %s\n", OptArg);exit(-1);} else FilterType = (filter_type) plop;
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'c': 
				Normalize = false;
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

	it = opts.find("--stat");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Extract_stat; }
	it = opts.find("--verbose");
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop;	Verbose=(Bool)plop; }
	it = opts.find("--SigmaNoise");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>SigmaNoise; }
	it = opts.find("--NSigma");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>NSigma; }
	it = opts.find("--MRout");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>Output_mr; }
	it = opts.find("--noMR");
	if(it!=opts.end()){ istringstream ss(it->second); ss>>No_mr; }
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
	it = opts.find("--force4sigma");
	if(it!=opts.end()){ istringstream ss(it->second); force4sigma=true; }
	it = opts.find("--normalize");
	if(it!=opts.end()){ istringstream ss(it->second); Normalize=false; }
	it = opts.find("--keep-single");
	if(it!=opts.end()){ istringstream ss(it->second); del_single=false; }
	it = opts.find("--no-fine");
	if(it!=opts.end()){ istringstream ss(it->second); no_fine=true; }
	
// Get short options
	filtinit(short_argc, short_argv);

	if (Verbose == True)
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale3D = " <<  NbrScale3D    << endl;  
		cout << "Compute_Recons = " <<  Compute_recons    << endl;
		cout << "OutputMR = " << Output_mr << endl;
		cout << "Statistics = " << Extract_stat << endl;
		if(!Normalize) cout << "NOT normalizing the transform output" << endl;
		if(no_fine)			cout << " Set the fine scale to zero : " <<  Niter    << endl;
		if(SigmaNoise>0)
		{
								cout << "Filtering type : " << string_filter_type(FilterType) <<endl;
								cout << " Noise level = " <<  SigmaNoise    << endl;
								cout << " NSigma Threshold = " <<  NSigma    << endl;
			if(force4sigma)		cout << " Force '4sigma lvl' at finest scale" << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
		}
		cout << endl;  
	}

	fltarray Data;
	fltarray *TabBand;
	fltarray Recons;
	fltarray Band;
	char* filename;
	
	POISSONATROUS_3D_WT *atrou = new POISSONATROUS_3D_WT();
	
// Logical links between inputs
	if(!Compute_recons && !No_mr) Output_mr=true;
			
// forward transform
	if(Compute_forward)
	{
	// Input data
		Name_Imag_In = add_fits(Name_Imag_In);
		fits_read_fltarr(Name_Imag_In, Data);

	// UOWT initialisation
		atrou->set_no_fine(no_fine);
		atrou->alloc(TabBand, Data.nx(), Data.ny(), Data.nz(), NbrScale3D);

	// Forward Transform
		atrou->transform(Data,TabBand,NbrScale3D);
	}
	else // We load the transformed cube
	{
	// read data
		Name_Imag_In = add_mr(Name_Imag_In);
		atrou->read(Name_Imag_In,TabBand,&NormalizeInv);
	
	// Inverse normalisation if normalized
		if(NormalizeInv) atrou->normalize_self(TabBand,1);
	}
	
// Statistics
	if(Extract_stat) 
		atrou->extract_stat(TabBand,Name_Imag_Out);

// Thresholding
	if( SigmaNoise>0 )
	{
		if(FilterType==FT_HARD || FilterType==FT_SOFT)
			atrou->threshold(TabBand, SigmaNoise*NSigma, FilterType==FT_SOFT, false);
		else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"' not implemented yet";
		if(del_single) atrou->clean_single(TabBand, SigmaNoise);
	}

// Output transform
	if(Output_mr)
	{
	// Normalize the transform if needed
		if(Normalize) atrou->normalize_self(TabBand, 0);
	// Save the transform
		filename = add_mr(Name_Imag_Out);
		atrou->write(filename, TabBand, Normalize);
	// Inverse normalisation if normalized
		if(Normalize) atrou->normalize_self(TabBand, 1);
	}
	
// Compute and output reconstruction
	if(Compute_recons)
	{
		atrou->recons(TabBand,Recons);
		filename = add_fits(Name_Imag_Out);
		writefltarr(filename, Recons);
	}
	
// Free memory	
	delete [] TabBand;
	delete atrou;
//	delete SB1D;
	for(int i=0;i<short_argc;i++)
		delete [] short_argv[i];
	delete [] short_argv;
	Data.~fltarray();
	Recons.~fltarray();
	
	stop=clock();
	if(Verbose) cerr<<endl<<"Execution time = "<<(stop-start)/CLOCKS_PER_SEC<<endl;
	exit(0);
}



