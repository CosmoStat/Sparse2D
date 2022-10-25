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
**    File:  mr3d_morpho.cc
**
*******************************************************************************
**
**    DESCRIPTION  Morphological operations on 3D volume : dilation/erosion
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

#include "Array.h"
#include "IM_IO.h"
#include "morpho3d.h"
#include "MGA_Inc.h"
#include "GetLongOptions.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

int Order = 1;

static void usage(char *argv[])
{
	// int i;
	fprintf(OUTMAN, "Usage: %s options in_cube result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [-n Order]\n");
	fprintf(OUTMAN, "             Order >0 : dilation");
	fprintf(OUTMAN, "             Order <0 : erosion");
	fprintf(OUTMAN, "             Default Order=1");
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
	while ((c = GetOpt(argc,argv,(char*)"n:")) != -1) 
	{
		switch (c) 
		{
			case 'n': 
				if (sscanf(OptArg,"%d",&Order) != 1) 
				{fprintf(OUTMAN, "Error: bad number of 3D scales parameter: %s\n", OptArg);exit(-1);}
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
	
	filtinit(argc, argv);
	
// Input data
	fltarray Data;
	strcat(Name_Imag_In,".fits");
	fits_read_fltarr(Name_Imag_In, Data);
	
	fltarray Out(Data.nx(),Data.ny(),Data.nz());
	int n = Order;
	if(n>0) morpho3d_dilation(Data,Out,2*n+1);
	else if(n<0) morpho3d_erosion(Data,Out,-2*n+1);
	
	char filename[64];
	sprintf(filename,"%s.fits",Name_Imag_Out);
	writefltarr(filename, Out);
}




