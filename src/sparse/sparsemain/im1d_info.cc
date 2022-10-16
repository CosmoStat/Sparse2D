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
**    File:  im1d_info.cc
**
*******************************************************************************
**
**    DESCRIPTION  return information about a time serie data set. 
**    ----------- 
**                 
******************************************************************************/

#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "Usage.h"
#include "MatrixOper.h"
#ifdef USELM
#include "Licence.h"
#endif
 
char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Imag_Out2[256];
char Name_Coeff[256];

 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
int NBShift = 10;
ar_order_detect_type AR_Order_Detect=AR_ORDER_BIC;
int BestNbrAR;
int MaxNbrAR=10;
Bool AutoCorr = False;

/***************************************************************************/

static void usage(char *argv[])
{
 
    fprintf(OUTMAN, "Usage: %s options input \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-a Nbr_of_Lags]\n");
    fprintf(OUTMAN, "             Calculate the autocorrelation function (AF)\n");
    fprintf(OUTMAN, "             with a given number of lags.\n");
    fprintf(OUTMAN, "             The output AF is saved in a file.\n");
    fprintf(OUTMAN, "                File Name is: autocor.\n");
    // fprintf(OUTMAN, "             Default number of lags is %d.\n", NBShift);
    manline();
    
    fprintf(OUTMAN, "         [-O Estimation_AR_Order_Method]\n");
    fprintf(OUTMAN, "             1: AIC \n");
    fprintf(OUTMAN, "             2: AICC\n");
    fprintf(OUTMAN, "             3: BIC\n");
    fprintf(OUTMAN, "             Default is %s method.\n", StringARDetectType(AR_Order_Detect));
    manline();
    
    fprintf(OUTMAN, "         [-M MaxAROrder]\n");
    fprintf(OUTMAN, "             Default is %d.\n", MaxNbrAR);
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
    while ((c = GetOpt(argc,argv,"a:O:M:v")) != -1) 
    {
	switch (c) 
        {
	   case 'O':  
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad order method: %s\n", OptArg);
	            exit(-1);
		}
                AR_Order_Detect = (ar_order_detect_type)(c-1);
		break;
	   case 'M':  
		if (sscanf(OptArg,"%d",&MaxNbrAR) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad maximum order value: %s\n", OptArg);
	            exit(-1);
		}
		break;
           case 'a':  
		if (sscanf(OptArg,"%d",&NBShift) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
		AutoCorr = True;
		strcpy(Name_Imag_Out, "autocor");
 		break;
	  case 'v': Verbose = True; break;
  	    case '?': usage(argv); break;
		}
	}
 
        /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
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

#ifdef USELM
    lm_check(LIC_M1D);
#endif
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    filtinit(argc, argv);

    io_1d_read_data(Name_Imag_In, Data);   
    // reform_to_1d(Data);
    int Nx = Data.nx(); 
    if (NBShift >= Nx-10)
    {
       cout << "Error: NBShift = " << NBShift << endl;
       cout << "       1 <  NBShift < " << Nx -10 << endl;    
    }
 
    // FitsHeader.origin = Cmd;
 
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       if (AutoCorr == True)
          cout << "File Name out = " << Name_Imag_Out << endl;
       cout << "AR Order Detection Method = " << StringARDetectType(AR_Order_Detect) << endl;
       if (AutoCorr == True)
          cout << "AF calculation: NBLags = " << NBShift <<  endl <<  endl;
    }
 
    cout << "  Np = " << Data.nx() << endl;
    float FluxIma =  Data.total();
    double Mean,Sigma,Skew,Curt;
    float Min,Max;
    moment4(Data.buffer(),Data.nx(), Mean,  Sigma, Skew,  Curt, Min,Max);
    cout << "  Min = " << Min << "  Max = " << Max << endl;
    cout << "  Mean = " << Mean << "  sigma = " << Sigma << endl;
    cout << "  Skew = " << Skew << "  Curtosis = " << Curt << endl;
    cout << "  Flux = " <<  FluxIma <<  "  Energy = " << (Data^2).total() << endl;

    fltarray BestARModel;
    AR_PREDICTION AR;
    AR.AR_Order_Detect = AR_Order_Detect;
    AR.MaxNbrAR = MaxNbrAR;
    AR.Verbose = Verbose;

    if (AutoCorr == True)
    {
        fltarray TabAuto(NBShift);
        autocor1d(Data, TabAuto, NBShift);
        io_1d_write_data(Name_Imag_Out, TabAuto);
    }
    AR.get_best_ar_model(Data, Data.nx(), BestNbrAR,  BestARModel);
    cout << "  Optimal AR order = " << BestNbrAR << endl;
    
    exit(0);
} 

