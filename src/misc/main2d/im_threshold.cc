/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/00
**    
**    File:  im_threshold.cc
**
*******************************************************************************
**
**    DESCRIPTION  image thresholding program
**    ----------- 
**                 
******************************************************************************/
   
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Rot.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool OptUp=False;
Bool OptAbs=False;
float ThresholdMin=0;
float ThresholdMax=0;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    manline();
    fprintf(OUTMAN, "         [-t ThresholdMin]\n");
    fprintf(OUTMAN, "             Threshold all values lower than ThresholdMin.\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    manline();
    fprintf(OUTMAN, "         [-T ThresholdMax]\n");
    fprintf(OUTMAN, "             Threshold all values larger than ThresholdMax.\n");
    fprintf(OUTMAN, "             Default is no thresholding.\n");

    manline();
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "             Compare the absolute value to the threshold.\n");
    fprintf(OUTMAN, "             Default is no.\n");
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
    while ((c = GetOpt(argc,argv,"at:T:vzZ")) != -1) 
    {
	switch (c) 
        { 
            case 'a': OptAbs=True;break;
            case 't': 
                if (sscanf(OptArg,"%f",&ThresholdMin) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad threshold: %s\n", OptArg);
                    exit(-1);
                    
                }
                break;
           case 'T': 
                if (sscanf(OptArg,"%f",&ThresholdMax) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad threshold: %s\n", OptArg);
                    exit(-1);
                    
                }
                OptUp=True;
                 break;
           case 'v': Verbose = True;break;
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
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/***************************************************************************/

int main(int argc, char *argv[])
{
    int k,i,j;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "ThresholdMin = " <<  ThresholdMin  << endl;
        if (OptUp == True) cout << "ThresholdMax = " <<  ThresholdMax  << endl; 
    }

    Ifloat Data;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    if (OptAbs == False)
    {
      for (i=0; i < Data.nl(); i++)
      for (j=0; j < Data.nc(); j++)
      {
        if (Data(i,j) < ThresholdMin) Data(i,j) = 0;
        else  if ((OptUp == True) && (Data(i,j) > ThresholdMax))
             Data(i,j) = ThresholdMax;
      }
    }
    else
    {
      for (i=0; i < Data.nl(); i++)
      for (j=0; j < Data.nc(); j++)
      {        
        float Val = Data(i,j);
        if (Val >= 0)  
        {
           if (Val < ThresholdMin) Data(i,j) =  0;
           else  if ((OptUp == True) && (Val > ThresholdMax))
                 Data(i,j) = ThresholdMax;
        }
        else 
        {
           Val = - Val;
           if (Val < ThresholdMin) Data(i,j) = 0;
           else  if ((OptUp == True) && (Val > ThresholdMax))
                                   Data(i,j) = -ThresholdMax;
        }
      }
    }
    io_write_ima_float (Name_Imag_Out, Data, &Header);
    exit(0);
}
