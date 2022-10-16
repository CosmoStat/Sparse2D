/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  28/08/00
**    
**    File:  tend_est.cc
**
*******************************************************************************
**
**    DESCRIPTION  tendency estimation 
**    ----------- 
**                 
**
******************************************************************************/

#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "Usage.h"
#include "MatrixOper.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Imag_Out2[256];
char Name_Coeff[256];

 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
int NbrPixTendancy = 100;
int Nt=-1;

/***************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options input out_tend out_signal_no_tend \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-T WindowSize_for_tendency_estimation]\n");
    fprintf(OUTMAN, "             Default is %d.\n", NbrPixTendancy);
    manline();
    fprintf(OUTMAN, "         [-f FirstPixels]\n");
    fprintf(OUTMAN, "             Default is the input signal size.\n");
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
    while ((c = GetOpt(argc,argv,"T:f:v")) != -1) 
    {
	switch (c) 
        {
           case 'T':  
		if (sscanf(OptArg,"%d",&NbrPixTendancy) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
 		break;
           case 'f':  
		if (sscanf(OptArg,"%d",&Nt) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
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

	if (OptInd < argc) strcpy(Name_Imag_Out2, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    fltarray Data,TrueData;
    char Cmd[256];
    // extern softinfo Soft;

    lm_check(LIC_M1D);
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    filtinit(argc, argv);

    io_1d_read_data(Name_Imag_In, Data);   
    // reform_to_1d(Data);
    int Nx = Data.nx();   

    if (Nt < 1) Nt = Nx;
    else if (Nt > Nx) Nt = Nx;
    // FitsHeader.origin = Cmd;
 
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "    NbrPixTendancy = " << NbrPixTendancy <<  endl;
    }
  
   fltarray Tend(Nt);
   tendancy_est(Data, Tend, NbrPixTendancy, Nt);
   fltarray ST(Nt);
   for (int i = 0; i < Nt; i++) ST(i) = Data(i) - Tend(i);
   io_1d_write_data(Name_Imag_Out, Tend);
   io_1d_write_data(Name_Imag_Out2, ST);
   exit(0);
} 

