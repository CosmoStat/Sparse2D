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
**    Date:  10/02/00
**    
**    File:  dct.cc
**
*******************************************************************************
**
**    DESCRIPTION  dct program
**    ----------- 
**                 
******************************************************************************/

#include "GlobalInc.h"
#include "IM_IO.h"
#include "Filter.h"
#include "SB_Filter1D.h"
#include "SB_Filter.h"


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int NbrScaleLine=-1;
int NbrScaleCol=-1;
Bool PNbrScaleLine = False;
Bool PNbrScaleCol = False;

sb_type_norm Norm = NORM_L2;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_border Bord = I_MIRROR;
Bool Reverse = False;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_data result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    sb_usage(SB_Filter);
    manline();
    
    fprintf(OUTMAN, "         [-x NbrScale_onXaxis]\n");
    fprintf(OUTMAN, "             Number of scales in the column direction.\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();
            
    fprintf(OUTMAN, "         [-y NbrScale_onYaxis]\n");
    fprintf(OUTMAN, "             Number of scales in the column direction.\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();
    
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Reverse transform.\n");
    manline();
        
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/
/*************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}
/*************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv,"rT:x:y:vzZ")) != -1) 
    {
	switch (c) 
        {
	   case 'r': Reverse = True; break; 
	   case 'T': 
 		SB_Filter = get_filter_bank(OptArg);
		break;	   
 	   case 'y': 
                if (sscanf(OptArg,"%d",&NbrScaleLine) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
		PNbrScaleLine = True;
		break;
           case 'x': 
                if (sscanf(OptArg,"%d",&NbrScaleCol) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
		PNbrScaleCol = True;
 	        break;
          case 'v': Verbose = True;break;
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
}
 
/*********************************************************************/
 
int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
#ifndef MRPOL    
    lm_check(LIC_MR4);
#else
    lm_check(LIC_POL);
#endif

    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;   
	if (PNbrScaleLine == True) cout << "Number of scales on Y-axis = " << NbrScaleLine    << endl; 
	if (PNbrScaleCol == True)  cout << "Number of scales on X-axis= " << NbrScaleCol    << endl;
	if (Reverse == True) cout << "Inverse Transform " << endl;
    }
 
    Ifloat Data;
    Ifloat Result;
    
    io_read_ima_float(Name_Imag_In, Data, &Header);

    float Min,Max;
    int Nx = Data.nl();
    int Ny = Data.nc();
    if (Verbose == True)
    {
       cout << "Nx = " << Nx <<  " Ny = " << Ny  <<  endl;
       cout << "Min = " << Data.min()  << "  Max = " << Data.max() <<  " Sigma = " << Data.sigma() << endl;
    }
   
    FilterAnaSynt FAS;
    // FAS.Verbose = Verbose;
    FAS.alloc(SB_Filter);
    SubBandFilter SBF(FAS,  Norm);
    LineCol LC;
    Bool UseMirror=False;
    LC.alloc(SBF, UseMirror);
    
    FilterAnaSynt *WT_SelectFilter = new FilterAnaSynt(SB_Filter);
    SubBandFilter *WT_SB1D = new SubBandFilter(*WT_SelectFilter, NORM_L2);
    // DirectionalLineCol LC1(*WT_SB1D);
    // LC1.NbrUndecimatedScaleCol = 0;
    // LC1.NbrUndecimatedScaleLine = 1;

    if (PNbrScaleLine == False)  NbrScaleLine = max_scale_number(Ny);
    if (PNbrScaleCol == False)  NbrScaleCol = max_scale_number(Nx);
    if (Verbose == True) cout << " NbrScaleY = " << NbrScaleLine << "  NbrScaleX = " << NbrScaleCol << endl; 
    
    if (Reverse == False) LC.transform(Data, NbrScaleCol, NbrScaleLine);
    else LC.recons(Data, NbrScaleCol, NbrScaleLine);
    io_write_ima_float(Name_Imag_Out, Data, &Header);
    
//     Ifloat Res;
//     if (Reverse == False) LC1.transform(Data, Res, NbrScaleCol, NbrScaleLine);
//     else LC1.recons(Data, Res, NbrScaleCol, NbrScaleLine);
//     
//     io_write_ima_float(Name_Imag_Out, Res, &Header);
    exit(0);
}
