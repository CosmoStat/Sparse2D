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
#include "IM_Rot.h"

char Name_Imag_In[256]; /* input file image */
char Name_List[256]; /* input file image */

char Name_List_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int NbrScaleLine=-1;
int NbrScaleCol=-1;
Bool PNbrScaleLine = False;
Bool PNbrScaleCol = False;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
     
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
    while ((c = GetOpt(argc,argv," x:y:vzZ")) != -1) 
    {
	switch (c) 
        {
//             case 'y': 
//                 if (sscanf(OptArg,"%d",&NbrScaleLine) != 1) 
//                 {
//                     fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
//                     exit(-1);
//                 }
// 		PNbrScaleLine = True;
// 		break;
//            case 'x': 
//                 if (sscanf(OptArg,"%d",&NbrScaleCol) != 1) 
//                 {
//                     fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
//                     exit(-1);
//                 }
// 		PNbrScaleCol = True;
//  	        break;
          case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
       

       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);
       if (OptInd < argc) strcpy(Name_List, argv[OptInd++]);
         else usage(argv);
	if (OptInd < argc) strcpy(Name_List_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
 
/*********************************************************************/

inline float bilin_inter(Ifloat &Data, float Xi, float Xj)
{
   int ii = (int) Xi;
   int jj = (int) Xj;
   int i1 = MIN(ii+1, Data.nl()-1);
   int j1 = MIN(jj+1, Data.nc()-1);
   
   float delta_x = (float) (Xj-jj);
   float delta_Val = Data (ii, j1) - Data (ii, jj);
   float ICol = Data (ii, jj)  + (delta_x * delta_Val);
   
   delta_x = (float) (Xi-ii);
   delta_Val = Data (i1, jj) - Data (ii, jj);
   float ILine = Data (ii, jj)  + (delta_x * delta_Val);
   return (ICol+ILine)/2.;

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
        cout << "File List Out = " << Name_List_Out << endl;   
    }
 
    fltarray Result;
    Bool Spline = True;
    if (Spline == True)
    {
       fltarray Data;
       fits_read_fltarr(Name_Imag_In, Data, &Header);
       BSPLINE_DEC BSD;
       BSD.SamplesToCoefficients(Data.buffer(), Data.nx(), Data.ny());
   
       fltarray List;
       fits_read_fltarr(Name_List, List);
       int Ni = List.nx();
       Result.alloc(Ni);
       for (int i=0; i < Ni; i++) 
         Result(i) = BSD.InterpolatedValue( Data, (double) List(i,0), (double) List(i,1));
    }
    else
    {
       Ifloat Data;
       io_read_ima_float(Name_Imag_In, Data, &Header);
       BSPLINE_DEC BSD;
       BSD.SamplesToCoefficients(Data.buffer(), Data.nx(), Data.ny());
   
       fltarray List;
       fits_read_fltarr(Name_List, List);
       int Ni = List.nx();
       Result.alloc(Ni);
       for (int i=0; i < Ni; i++) 
          Result(i) = bilin_inter( Data, (float) List(i,0), (float) List(i,1));
    }
    
    fits_write_fltarr(Name_List_Out, Result);
    exit(0);
}
