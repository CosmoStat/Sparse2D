/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  13/04/01
**    
**    File:  im1d_stf.cc
**
*******************************************************************************
**
**    DESCRIPTION  Short term Fourier transform 
**    -----------  
**                 
**
******************************************************************************/


#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "IM1D_IO.h"
#include "FFTN_1D.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
Bool Reverse = False;
float WinParam = 0.5;
int Step = 0;
int WindowSize = 1024;
type_std_win WindowType = W_HAMMING;
Bool DebugSTF = False;

type_spec TSpec = DEF_SPEC;
type_std_win WindowTypeFreq = W_HAMMING;

/***************************************************************************/
      
 
static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    spec_usage();
    manline();

    fprintf(OUTMAN, "         [-T type_of_window (time domain)]\n");
    for (int i = 0; i < NBR_STD_WIN ; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringSTDWin((type_std_win )i));
    fprintf(OUTMAN, "             default is %s.\n", StringSTDWin(WindowType));
    manline();

 //    fprintf(OUTMAN, "         [-F type_of_window (freq. domain)]\n");
//     fprintf(OUTMAN, "             Only for Choi-Williams distribution.\n");
//     fprintf(OUTMAN, "             default is %s.\n", StringSTDWin(WindowTypeFreq));
//     manline();
    
    fprintf(OUTMAN, "         [-w window_size (time domain)]\n");
    fprintf(OUTMAN, "             Window size. Default is 1024.\n");
    manline();

    fprintf(OUTMAN, "         [-W window_param (time domain)]\n");
    fprintf(OUTMAN, "             Window parameter. Default is 0.5.\n");
    manline();
    
    fprintf(OUTMAN, "         [-S Step]\n");
    fprintf(OUTMAN, "             Step between two points. Default is WindowSize/2.\n");
    manline();

    verbose_usage();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    /* get options */
    while ((c = GetOpt(argc,argv,"t:F:T:W:w:S:v")) != -1) 
    {
	switch (c) 
        {
 	  case 't':
	       if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "bad distribution type: %s\n", OptArg);
		    exit(-1);
		}
                if ((c <= 0) && (c > NBR_SPEC))   
                {
		    fprintf(OUTMAN, "bad distribution type: %s\n", OptArg);
	            usage(argv);
 		}
		TSpec = (type_spec ) (c-1);
                break;
  	  case 'T':
	       if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
		    exit(-1);
		}
                if ((c <= 0) && (c >  NBR_STD_WIN))   
                {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
	            usage(argv);
 		}
		WindowType = (type_std_win) (c-1);
                break;
 	  case 'F':
	       if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
		    exit(-1);
		}
                if ((c <= 0) && (c >  NBR_STD_WIN))   
                {
		    fprintf(OUTMAN, "bad window type: %s\n", OptArg);
	            usage(argv);
 		}
		WindowTypeFreq = (type_std_win) (c-1);
                break;
 	   case 'W':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%f",&WinParam) != 1) 
                {
		    fprintf(OUTMAN, "bad window parameter: %s\n", OptArg);
		    exit(-1);
		}
		break; 
	   case 'w':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%d",&WindowSize) != 1) 
                {
		    fprintf(OUTMAN, "bad window size: %s\n", OptArg);
		    exit(-1);
		}
		break;
	   case 'S':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%d",&Step) != 1) 
                {
		    fprintf(OUTMAN, "bad step: %s\n", OptArg);
		    exit(-1);
		}
		break;		
 	  case 'v': Verbose = True; break;
  	    case '?': usage(argv); break;
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
		usage(argv);
	}
}

 
/*********************************************************************/

int main(int argc, char *argv[])
{
    fltarray Data;
    fltarray TimeFreq;
    Icomplex_f STFDat;
    char Cmd[256];
    int  Nx,Nlt,Nct;
    ST_FFTN STF;
               
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    lm_check(LIC_M1D);
    filtinit(argc, argv);
   
    io_1d_read_data(Name_Imag_In, Data);
    // reform_to_1d(Data);
    // FitsHeader.origin = Cmd;
    Nx = Data.nx();
    if (WindowSize <= 0)  WindowSize = Nx / 4;
    else if (WindowSize >= Nx) WindowSize = Nx / 4;
    if (Step <= 0) Step = WindowSize  / 2;
         
    STF.alloc(Nx, WindowType, WinParam, WindowSize, Step);
    Nlt = STF.nl(Nx);
    Nct = STF.nc(Nx);
 
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name out = " << Name_Imag_Out << endl; 
       cout << "Distribution type = " << StringSPEC(TSpec) << endl;
       if (TSpec != SPEC_WIGNER_VILLE)
       {
          cout << "Window Type  = " <<  StringSTDWin(WindowType) << endl;
          cout << "Window size  = " << WindowSize << endl;
          cout << "Step  = " << Step << endl;
          cout << "Window parameter = " << WinParam << endl;
       }  
       cout << "Data size  = " << Nx << endl;
       cout << "Output distrib. size: Nl =  " << Nlt/2 << " Nc = " << Nct << endl;
    }
 
   switch(TSpec)
   {
      case SPEC_STF:
          STF.spectogram(Data, TimeFreq);
	  break;
      case SPEC_WIGNER_VILLE:
	  STF.wigner_wille(Data, TimeFreq);
	  break;
      case SPEC_CHOI_WILLIAMS:
	  STF.choi_williams(Data, TimeFreq, WindowTypeFreq);
	  break;
   }
   // cout << "Write " << Name_Imag_Out << TimeFreq.nl() << " " << TimeFreq.nc() << endl;
   io_1d_write_data(Name_Imag_Out, TimeFreq);
  
   exit(0);
}

