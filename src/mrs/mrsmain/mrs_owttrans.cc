/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  25/10/08
**    
**    File:  mrs_almtrans.cc
**
*******************************************************************************
**
**    DESCRIPTION:  2D Orthogonal Wavelet Transform on the sphere (each face separately)
**    -------------------
**
******************************************************************************
**
**    Modifications from mrs_owttrans.c : Original version only for testing that
**		forward and inverse transform were matching. 
**
**      27/09/11 : FCS : added InverseFlag to allow 
**		separate forward and backward transform , and choice of number of wavelet scales (not the full number).
**
******************************************************************************/

#include"MRS_Sparse.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
  
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

bool Verbose=false;
bool InverseFlag=false;
float SigmaNoise=0.;
 
int NbrScale=0;

/***************************************************************************/

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated.  \n");
  
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "n:vizZ")) != -1) 
    {
	switch (c)  
        { 
            case 'n':
 	     	      if (sscanf(OptArg,"%d",&NbrScale) != 1)  
                  {
		            fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
		            exit(-1); 
	            	}
 		          break;
 		    case 'i':
 	     	 	InverseFlag=true;
 		  		break;
     	    case 'v': Verbose = true; break;
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
/*********************************************************************/


int main(int argc, char *argv[])
{
    int k;
    double Min, Max;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == true)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (NbrScale > 0)  cout << "# NbrScale = " <<  NbrScale << endl;
     }

   Hdmap Map;
   Hdmap MapOut;
     
   Map.read(Name_Imag_In);
   if (Verbose == true)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == true) Map.info((char*) "Input MAP");
   } 

   int Nside = Map.Nside();
    
   C_OWT WT;
   WT.alloc(Nside, NbrScale, DEF_MRS_ORDERING);
   printf("NbrScale= %d\n",NbrScale);
   
   if(InverseFlag==false) {
   		WT.transform(Map);
   		(WT.WTTrans).write(Name_Imag_Out);
   } else {
     // Result.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
   		MapOut.alloc(Nside, DEF_MRS_ORDERING);
    	MapOut.fill(0.);
   		WT.WTTrans=Map;
   		WT.recons(MapOut);
   
   		MapOut.write(Name_Imag_Out);
   }  
   exit(0);
}

