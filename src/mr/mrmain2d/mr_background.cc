/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  mr_background.cc
**
*******************************************************************************
**
**    DESCRIPTION  substract the background to an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_background option image output
**        where options = 
**
**
**           [-n number_pixels]
**
**
**           [-w background]
**                if this option is set, the background image is created.
**                Default is no.
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_background.cc 3.1 96/05/02 CEA 1995 @(#)";
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_BGR[256];          /* background file name */
int Nbr_Plan=1;           /* number of scales */
type_transform Transform =  TM_PYR_MEDIAN; /* type of transform */
Bool WriteBGR = False;          /* write the background on the disk */

int NPix = 16;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
   
Bool Verbose = False;
int MedWind = 0;

/*********************************************************************/

static void usage(char *argv[])
{
 
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
   
    mrbgr_option_usage(NPix);
    manline();   
    fprintf(OUTMAN, "         [-W MedianWindowSize]\n");
    fprintf(OUTMAN, "             Median window size. \n");
    fprintf(OUTMAN, "             Default is 5.\n");
    manline();
    vm_usage();
    manline();
    verbose_usage();
      manline();    
     manline();    
     exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   
    /* get options */
    while ((c = GetOpt(argc,argv,"lW:n:w:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'l': Transform =  TM_TO_PYR; break;
 	   case 'v': Verbose = True; break;
	   case 'n':
		/* -n <Number_of_pixels> */
		if (sscanf(OptArg,"%d",&NPix) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of pixels: %s\n", OptArg);
		    exit(-1);
		}
		break;
  	   case 'W':
		/* -n <Number_of_pixels> */
		if (sscanf(OptArg,"%d",&MedWind) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad median window size: %s\n", OptArg);
		    exit(-1);
		}
		break;
	   case 'w':
		/* -w < background file name> */
		if (sscanf(OptArg,"%s",Name_BGR) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteBGR = True;
 		break;
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
 	    case '?': usage(argv); break;
	    default: usage(argv); break;
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
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    int  Ns;
    Ifloat Dat;
    int k;
    char Cmd[512];
    fitsstruct Header;
	
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
       cout << "Transform = " << StringTransform(Transform) << endl;
       cout << "Number of pixels = " << NPix << endl;
    }

    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Header.origin = Cmd;
    Ifloat Result (Dat.nl(), Dat.nc(), "Result Background");
 
    if ((NPix < 2) || (NPix >= Dat.nl()) || (NPix >= Dat.nc()))
    {
        printf ("Error: Bad number of pixels\n");
        exit(-1);
    }
    Ns = MIN(Dat.nl(), Dat.nc());
    
    if (Nbr_Plan == 1)
    {
        while (NPix < Ns) 
        {
           Nbr_Plan ++;
           Ns /= 2;
        }
    }
    if (Nbr_Plan < 3) Nbr_Plan = 3;
    if (Verbose == True) cout << "Number of scales = " << Nbr_Plan << endl;

    MultiResol MR_Data (Dat.nl(), Dat.nc(), Nbr_Plan, 
                        Transform, "MR_Transform");
    MR_Data.Border=I_MIRROR;
    if (MedWind > 0)
    {
       if (MedWind % 2 == 0) MedWind++;
       MR_Data.MedianWindowSize = MedWind;
    }
    MR_Data.transform (Dat);
    for (int b = 0; b < MR_Data.nbr_band()-1; b++) MR_Data.band(b).init();
    MR_Data.recons(Result);
    Dat -= Result;
    Header.bitpix = BP_FLOAT;
    io_write_ima_float(Name_Imag_Out, Dat, &Header);

    /* background image creation */
    if (WriteBGR == True)  
    {
        io_write_ima_float(Name_BGR, Result, &Header);
    }
    exit(0);
} 

