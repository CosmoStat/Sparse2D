/******************************************************************************
 **                   Copyright (C) 2007 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Jean-Luc Starck - Jerome Bobin
 **
 **    Date:  11/30/07
 **
 **    File:  mr_gmca.cc
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Generalized Morphological Component Analysis
 **    -----------
 **
 **    Usage: mr_gmca options cube output
 **
 ******************************************************************************/
#include "DefMath.h"
#include "MatrixOper.h"  // defined in $TOOLS
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_Block.h"
#include "C_OMP.h"
#include "C_DL1D.h"
#include "C_ESTIMATE_REDSHIFT.h"
Bool Ov = True;
bool UseSparsity=true;
/****************************************************************************/
char Name_Dico[256];
char Name_Vect_in[256];
char Name_Shift_in[256];
char Name_Data_out[256];
char Name_RMS_in[256] = "";
double ErrorTarget = -1.;
extern int  OptInd;
extern char *OptArg;
dblarray Dico;
dblarray Image_out;
dblarray Image_in;
dblarray stacked_patches;
dblarray rms;
double lambda_step = log10(6003)-log10(6000)  ;
int shift_offset = 36;
int shift_range = 5;
int defaultSparsityTarget = 5;
int SparsityTarget = defaultSparsityTarget;
int OverlapNumber = 1;
extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose=False;
int Datanx,Datany,Datanz,Diconx,Dicony,Diconz,indData,indDico;
/****************************************************************************/
static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options input_vect input_dico shift_tab output_shift", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");
	fprintf(OUTMAN, "         [-S TargetSparsity]\n");
	fprintf(OUTMAN, "         [-E ErrorTarget, TargetSparsity ignored]\n");
	fprintf(OUTMAN, "         [-R RMS curved filename, flat by default]\n");
	fprintf(OUTMAN, "         [-N Testing 2N+1 shift values]\n");
	fprintf(OUTMAN, "         [-L Redshift quantization parameter]\n");
	vm_usage();
	manline();
	exit(-1);
}
/****************************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
	int c;
#ifdef LARGE_BUFF
	int VMSSize=-1;
	Bool OptZ = False;
	char VMSName[1024] = "";
#endif
	/* get options */
	while ((c = GetOpt(argc,argv,"OS:E:R:vzZ:N:L:")) != -1)
	{
		switch (c)
		{
		case 'S':
			/* -S <SparsityTarget> */
			if (sscanf(OptArg,"%d", &SparsityTarget) != 1)
			{
				fprintf(OUTMAN, "bad value: %s\n", OptArg);
				exit(-1);
			}
			break;
		case 'E':
			sscanf(OptArg,"%lf", &ErrorTarget);
			if (ErrorTarget > 0) UseSparsity = false;
			break;
		case 'R':
			sscanf(OptArg,"%s", &Name_RMS_in);
			break;
		case 'O': Ov = (Ov ==True) ? False: True; break;
		case 'v': Verbose = True; break;
		case 'N': sscanf(OptArg,"%d", &shift_range);
		case 'L': sscanf(OptArg,"%lf", &lambda_step);
#ifdef LARGE_BUFF
		case 'z':
			if (OptZ == True)
			{
				fprintf(OUTMAN, "Error: Z option already set...\n");
				exit(-1);
			}
			OptZ = True;
			break;
		case 'Z':
			if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
			{
				fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
				exit(-1);
			}
			if (OptZ == True)
			{
				fprintf(OUTMAN, "Error: z option already set...\n");
				exit(-1);
			}
			OptZ = True;
			break;
#endif
		case '?':
			usage(argv);
		}
	}
	/* get optional input file names from trailing
          parameters and open files */

	if (OptInd < argc) strcpy(Name_Vect_in, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Dico, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Shift_in, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Data_out, argv[OptInd++]);
	else usage(argv);
	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
	if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/****************************************************************************/
int main(int argc, char *argv[])
{
	/* Get command line arguments, open input file(s) if necessary */
	fitsstruct Header;
	char Cmd[512];
	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
	transinit(argc, argv);
	if (Verbose == True)
	{
		cout << "Dictionary Filename = " << Name_Dico << endl;
		cout << "Vector Filename in = " << Name_Vect_in  << endl;
		cout << "Redshift Filename in = " << Name_Shift_in  << endl;
		cout << "Redshift Filename out = " << Name_Data_out  << endl;
		if (UseSparsity  == True) cout << "Use Sparsity target = " << SparsityTarget << endl;
		else cout <<"Use Error target " << ErrorTarget << endl;
	}
	if (Verbose == True) cout << "Reading the data"<< endl;
	dblarray Dico, Data, Redshift;
	fits_read_dblarr (Name_Dico,Dico);
	int Nx = Dico.nx();
	int Ny = Dico.ny();
	int PatchSize = sqrt( double( Nx) );
	if (Verbose == True)
		if (Dico.naxis() == 2)
			cout << "Dictionary: " << Nx << " atoms of length " << Ny << endl;
	fits_read_dblarr (Name_Vect_in, Data);
	Nx = Data.nx();
	Ny = Data.ny();
	int Nz = Data.nz();
	if (Verbose == True)
		cout << "Input vector : Nx = " << Data.nx() << " Ny = " << Data.ny() << endl;
	fits_read_dblarr (Name_Shift_in, Redshift);
	if (Verbose == True)
		cout << "Redshift vector : " << Redshift.nx() << " redshift values " << endl;
	intarray Shift(Redshift.nx(),2*shift_range+1);
	for (int k=0;k<Redshift.nx();k++)
		for (int l=-shift_range;l<shift_range+1;l++)
		{
			Shift(k,l+shift_range) = shift_offset + round(log10(1+Redshift(k))/lambda_step) + l;
			// cout << dblShift(k,l) << " , " << Shift(k,l+shift_range) << endl;
		}
/*	fits_write_intarr("Shift.fits",Shift);
	cout << Shift.min() << " , " << Shift.max() << endl;
	dblarray dblShift(Shift.nx(),Shift.ny());
	for (int i=0;i<Shift.nx();i++)
		for (int j=0;j<Shift.ny();j++)
			dblShift(i,j) = Shift(i,j);
	fits_write_dblarr ("dblShift.fits",dblShift);*/
	if (Verbose == True)
		cout << "Shift vector : " << Shift.nx() << " series of " << Shift.ny() << " redshift values " << endl;
	// Fixing stop criterion for OMP
	if (UseSparsity == false) SparsityTarget = Dico.nx();
	else ErrorTarget = -1;
	// Reading RMS curve if given
	if (strlen(Name_RMS_in) == 0)
	{
		if (Verbose == True)
			cout << "No RMS curve specified, assuming flat noise" << endl;
		rms.alloc(Data.nx());
		for (int k=0;k<Data.nx();k++)
			rms(k) = 1;
	}
	else
		fits_read_dblarr (Name_RMS_in, rms);

	/****************************************************************************/
	// Calling main function
	C_ESTIMATE_REDSHIFT estimator(Dico);
	dblarray optimal_shift = estimator.estimate_redshift(Data,Shift,SparsityTarget,ErrorTarget,rms,Verbose);
	dblarray optimal_redshift(optimal_shift.nx());
	for (int i=0;i<optimal_shift.nx();i++)
		optimal_redshift(i) = pow(10,(optimal_shift(i) - shift_offset)*lambda_step) - 1;
	fits_write_dblarr (Name_Data_out, optimal_redshift);
	if (Verbose == True)
		cout << "Redshift estimation complete" << endl;
	exit(0);
}
/****************************************************************************/
