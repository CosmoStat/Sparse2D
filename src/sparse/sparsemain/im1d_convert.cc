/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  21/09/98
**    
**    File:  im1d_convert.cc
**
*******************************************************************************
**
**    DESCRIPTION   
**    ----------- 
**                 
**
**
******************************************************************************/


#include "Array.h"
#include "IM_Obj.h"
#include "IM1D_IO.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

int NSkip=0;
 
/****************************************************************************/

static void usage(char *argv[])
{
    
    fprintf(stdout, "Usage: %s options image output\n\n", argv[0]);
    fprintf(stdout, "   FITS format (.fits)  \n");
    fprintf(stdout, "   ACII format (.dat)  \n");
    fprintf(stdout, "   EXEL format (.csv)  \n");
    
    fprintf(stdout, "   where options =  \n");

    fprintf(stdout, "         [-s line_number_to_skip]\n");
    fprintf(stdout, "              Skip the s first lines of the input file.\n");

    fprintf(stdout, "\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"s:")) != -1) 
    {
	switch (c) 
        {
	  case 's':
 		if (sscanf(OptArg,"%d",&NSkip) != 1) 
                {
		    fprintf(stdout, "bad number of lines: %s\n", OptArg);
		    usage(argv);
		}
		break; 	   
	   case '?':
			usage(argv);
		}
	}
 
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(stdout, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}
 

/****************************************************************************/
 

int main(int argc, char *argv[])
{
    fltarray Dat;
    fltarray Mallat;
    /* Get command line arguments, open input file(s) if necessary */
    transinit(argc, argv);

    //  Dat.fits_read(Name_Imag_In);
    io_1d_read_data (Name_Imag_In,  Dat, NSkip);
    cout <<  endl << endl;

    cerr << "CONVERT " << String1DFormat(io_detect_1dformat(Name_Imag_In));
    cerr << " to " << String1DFormat(io_detect_1dformat(Name_Imag_Out)) << endl;

    io_1d_set_format(Name_Imag_Out);
    io_1d_write_data (Name_Imag_Out, Dat);

    exit(0);
} 

