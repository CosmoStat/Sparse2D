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
 
/***************************************************************************/
      
 
static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options singal_in image_out\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_window]\n");
    for (int i = 0; i < NBR_STD_WIN ; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringSTDWin((type_std_win )i));
    fprintf(OUTMAN, "             default is %s.\n", StringSTDWin(WindowType));
    manline();
    
    fprintf(OUTMAN, "         [-w window_size]\n");
    fprintf(OUTMAN, "             Window size. Default is 1024.\n");
    manline();

    fprintf(OUTMAN, "         [-W window_param]\n");
    fprintf(OUTMAN, "             Window parameter. Default is 0.5.\n");
    manline();
    
    fprintf(OUTMAN, "         [-S Step]\n");
    fprintf(OUTMAN, "             Step between two points. Default is WindowSize/2.\n");
    manline();
        
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Reverse transform.\n");
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
    while ((c = GetOpt(argc,argv,"t:W:w:S:rvx")) != -1) 
    {
	switch (c) 
        {
	  case 'x':DebugSTF = True; break;
	  case 't':
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
 	  case 'r': Reverse = True; break;
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
    Ifloat STF_Spec;
    Icomplex_f STFDat;
    char Cmd[256];
    int i,j,Nx,Nlt,Nct;
    char Name_Imag_Out_Spec[256];
    ST_FFTN STF;
    complex_f *Buff_STF;
              
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    lm_check(LIC_M1D);
    filtinit(argc, argv);
    
    if (Reverse == False)
    {
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
       STFDat.alloc(Nlt, Nct, "STF");
    }
    else
    {
       //cout << "READ COMPLEX: " << Name_Imag_In << endl;
       io_read_ima_complex_f (Name_Imag_In, STFDat);
       //cout << "CF: " << STFDat.nl() << " " << STFDat.nc() << endl;
       Nlt = STFDat.nl();
       Nct = STFDat.nc();
       WindowSize = Nlt;
       if (Step <= 0) Step = WindowSize  / 2;
       Nx = Nct*Step; 
       //cout << "ALLOC" << Nx << " WindowSize = " << WindowSize << " Step = " << Step << endl;
       STF.alloc(Nx, WindowType, WinParam, WindowSize, Step);
       Data.alloc(Nx);
    }
    Buff_STF = STFDat.buffer();

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name out = " << Name_Imag_Out << endl;
       cout << "Window Type  = " <<  StringSTDWin(WindowType) << endl;
       cout << "Window size  = " << WindowSize << endl;
       cout << "Step  = " << Step << endl;
       cout << "Window parameter = " << WinParam << endl;
       cout << "Data size  = " << Nx << endl;
       cout << "STF size: Nl = " << Nlt << " Nc = " << Nct << endl;
    }
 
    if (DebugSTF == True)
    {
       float W=WinParam;
       fltarray Win(Nx);
       STD_Window STDW;
       STDW.hamming(Win, W);
       io_1d_write_data("xx_hamming.fits", Win);
       STDW.hanning(Win, W);
       io_1d_write_data("xx_hanning.fits", Win);
       STDW.gaussian(Win, W);
       io_1d_write_data("xx_gaussian.fits", Win);
       STDW.blackman(Win, W);
       io_1d_write_data("xx_blackman.fits", Win);
    }

    if (Reverse == False)
    {
       STF.transform(Data, Buff_STF);
       io_write_ima_complex_f (Name_Imag_Out, STFDat);
       STF_Spec.alloc(Nlt/2,Nct,"STF spec");
       for (i=0; i < Nlt/2; i++) 
       for (j=0; j < Nct; j++)
       {
          float Re = Buff_STF[i*Nct+j].real(); //  / (float) Nx;
          float Im = Buff_STF[i*Nct+j].imag(); //  / (float) Nx;
          STF_Spec(i,j) = Re*Re+Im*Im;
       }
       io_write_ima_complex_f (Name_Imag_Out, STFDat);
       sprintf(Name_Imag_Out_Spec,"%s_spec.fits", Name_Imag_Out);
       io_write_ima_float(Name_Imag_Out_Spec, STF_Spec);
    }
    else
    {
        STF.recons(Buff_STF, Data, True);
        io_1d_write_data(Name_Imag_Out, Data);
    }
    exit(0);
}

