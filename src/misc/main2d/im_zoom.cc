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
// #include "Filter.h"
#include "IM_Rot.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int NbrScaleLine=-1;
int NbrScaleCol=-1;
float ZoomVal=0.;
Bool PNbrScaleLine = False;
Bool PNbrScaleCol = False;

// sb_type_norm Norm = NORM_L2;
// type_sb_filter SB_Filter = F_MALLAT_7_9;
type_border Bord = I_MIRROR;
Bool PZoomVal = False;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-x Nbr_Pixels_in_x_dir]\n");
    fprintf(OUTMAN, "             Default is twice the number of the input image.\n");
    manline();
    fprintf(OUTMAN, "         [-y Nbr_Pixels_in_y_dir]\n");
    fprintf(OUTMAN, "             Default is twice the number of the input image.\n");
    manline();
    fprintf(OUTMAN, "         [-r RebinParameter]\n");
    fprintf(OUTMAN, "             Resize the input image by the value defined by RebinParameter. Default is 2.\n");
    manline();
    manline();
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;  

    /* get options */
    while ((c = GetOpt(argc,argv," x:y:r:vzZ")) != -1) 
    {
	switch (c) 
        {
	    case 'r': 
	       if (sscanf(OptArg,"%f",&ZoomVal) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
		PZoomVal = True;
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
       if ((PZoomVal == True) && ( (PNbrScaleLine == True) || (PNbrScaleCol == True)))
       {
          fprintf(OUTMAN, "Error: z option cannot be used with x or y option: %s\n", OptArg);
          exit(-1);
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
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;   
    }
 
    Ifloat Data;
    Ifloat Result;
    
    io_read_ima_float(Name_Imag_In, Data, &Header);
    int Nx = Data.nx();
    int Ny = Data.ny();
    if (PNbrScaleLine == False) NbrScaleLine = 2*Ny;
    if (PNbrScaleCol == False) NbrScaleCol = 2*Nx;
    if (PZoomVal == True)
    {
       NbrScaleLine = (int) (ZoomVal*Ny+0.5);
       NbrScaleCol =  (int) (ZoomVal*Nx+0.5);
    }
    
    Result.alloc(NbrScaleCol, NbrScaleLine);
    float ZoomX = (float) NbrScaleCol / Nx;
    float ZoomY = (float) NbrScaleLine / Ny;
    if (Verbose == True) 
    { 
        cout << "ZoomX = " << ZoomX << "  Zoomy " << ZoomY << endl;
	cout << "Output Image size: Nx  = " << NbrScaleCol << ", Ny = " << NbrScaleLine << endl;
    }
//    BSPLINE_DEC BSD;
//    BSD.SamplesToCoefficients(Data.buffer(), Data.nx(), Data.ny());
    
//     for (int i=0; i < NbrScaleLine; i++) 
//     for (int j=0; j < NbrScaleCol; j++) 
//     {
//        double x = j / ZoomX;
//        double y = i / ZoomY;
//        Result(j,i) = BSD.InterpolatedValue( Data, (double) x, (double) y);
//     }
    im_zoom(Data, Result);
    io_write_ima_float(Name_Imag_Out, Result);
    exit(0);
}
