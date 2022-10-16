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
#include "C_EXTRACT_PATCHES.h"
Bool Ov = True;
bool UseSparsity=true;
/****************************************************************************/
char Name_Image_in[256];
char Name_Data_out[256];
extern int  OptInd;
extern char *OptArg;
int PatchSize,Npatch;
dblarray Data_out;
dblarray Image_in;
bool Verbose = False;
/****************************************************************************/
static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options Image Data \n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-v Verbose mode]\n");
	fprintf(OUTMAN, "         [-W Patch width]\n");
	fprintf(OUTMAN, "         [-N Number of patches]\n");
	verbose_usage();
	vm_usage();
	manline();
	exit(-1);
}
/****************************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
	int c;
	/* get options */
	while ((c = GetOpt(argc,argv,"W:N:v")) != -1)
	{
		switch (c)
		{
		case 'W':
			/* -W <Patch width> */
			sscanf(OptArg,"%d", &PatchSize);
			break;
		case 'N':
			/* -N <Number of patches> */
			sscanf(OptArg,"%d", &Npatch);
			break;

		case 'v': Verbose = True; break;
		case '?':
			usage(argv);
		}
	}
	/* get optional input file names from trailing
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Image_in, argv[OptInd++]);
	else usage(argv);
	if (OptInd < argc) strcpy(Name_Data_out, argv[OptInd++]);
	else usage(argv);
	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
/****************************************************************************/
int main(int argc, char *argv[]) 
{
	/* Get command line arguments, open input file(s) if necessary */
	char Cmd[512];
	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
	transinit(argc, argv);
	if (Verbose == True)
	{
		cout << "Image Filename in = " << Name_Image_in  << endl;
		cout << "Data Filename out = " << Name_Data_out  << endl;
	}
	if (Verbose == True) cout << "Reading the data"<< endl;
	dblarray Image_in;
	int Nx,Ny;
	fits_read_dblarr (Name_Image_in, Image_in);
	Nx = Image_in.nx();
	Ny = Image_in.ny();
	if (Verbose == True)
		cout << "Input image: Nx = " << Image_in.nx() << " Ny = " << Image_in.ny() << endl;



	/****************************************************************************/
	// Calling main function
	C_EXTRACT_PATCHES extractor(Image_in);
	dblarray Data_out = extractor.extract(PatchSize,Npatch);
	if (Verbose == true)
		cout << Data_out.ny() << "  patches of "<< Data_out.nx() << " pixels extracted." << endl;
	fits_write_dblarr (Name_Data_out, Data_out);
	exit(0);
}
