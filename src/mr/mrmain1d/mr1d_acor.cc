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

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "Usage.h"
#include "MR1D_Obj.h"
#include "NeuNet.h"
#ifdef USELM
#include "Licence.h"
#endif
#include "MatrixOper.h"
#include "MR1D_Predict.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Imag_Out2[256];
char Name_Coeff[256];

 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
static Bool Debug = False;

int NBShift = 10;
int NbrScale=5;           // number of scales
type_trans_1d Transform=TO1_PAVE_HAAR;    // type of transform

/***************************************************************************/

static void usage(char *argv[])
{
    // int i;

    fprintf(OUTMAN, "Usage: %s options input out_tend out_signal_no_tend \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    fprintf(OUTMAN, "             Default is 5.\n");
    manline();
    
    fprintf(OUTMAN, "         [-S Nbr_of_Shift]\n");
    fprintf(OUTMAN, "             Default is %d.\n", NBShift);
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
    while ((c = GetOpt(argc,argv,"n:S:xv")) != -1) 
    {
	switch (c) 
        {
	   case 'n':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
                if ((NbrScale <= 1) || (NbrScale > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    usage(argv);
		}
 		break;
	   case 'x': Debug = True; break; 
           case 'S':  
		if (sscanf(OptArg,"%d",&NBShift) != 1) 
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

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}

/****************************************************************************/

void mr1d_autocor(MR_1D & MR_Data, fltarray &TabAuto, int Np, int NBShift)
{
   int b,i,s;
   double XY,X2,Y2;
   if ((TabAuto.nx() != NBShift) || (TabAuto.ny() != MR_Data.nbr_scale()))
                                TabAuto.alloc(NBShift, MR_Data.nbr_scale());
   int Bord=1;
   for (b=0; b < MR_Data.nbr_scale(); b++)
   {
      TabAuto(0,b) = 1;
      for (s=Bord; s < NBShift; s++)
      {
         XY = X2 = Y2 = 0.;
         for (i=s; i < Np; i++)
         {
           XY += MR_Data(b,i) * MR_Data(b,i-s);
           X2 += MR_Data(b,i) * MR_Data(b,i);
           Y2 += MR_Data(b,i-s)*MR_Data(b,i-s);
         }
         TabAuto(s,b) = (float) ( XY / (sqrt(X2)*sqrt(Y2)));
      }
   }
}

/****************************************************************************/
 
 
int main(int argc, char *argv[])
{
    fltarray Data,TrueData;
    char Cmd[256];
    // int i,b,j;    
    extern type_1d_format IO_1D_Format;
    // extern softinfo Soft;

#ifdef USELM
    lm_check(LIC_M1D);
#endif

    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    filtinit(argc, argv);

    io_1d_read_data(Name_Imag_In, Data);   
    reform_to_1d(Data);
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
       cout << "File Name out = " << Name_Imag_Out << endl;
       cout << "    NBShift = " << NBShift <<  endl;
    }
    
    MR_1D MR_Data (Nx, Transform, "MR_Data", NbrScale);
    if (Transform !=TO1_PAVE_HAAR) MR_Data.Border = I_MIRROR;
    else MR_Data.Border = I_CONT;
    MR_Data.transform(Data);
    fltarray TabAuto(NBShift, NbrScale);
    mr1d_autocor(MR_Data, TabAuto, Nx, NBShift);
    if (IO_1D_Format == F1D_FITS)
            fits_write_fltarr(Name_Imag_Out, TabAuto);
    else io_write2d_ascii(Name_Imag_Out, TabAuto);
    exit(0);
} 

