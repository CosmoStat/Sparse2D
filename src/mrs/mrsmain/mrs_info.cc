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
**    File:  mrs_info.cc
**
*******************************************************************************
**
**    DESCRIPTION:  print information about a Healpix image
**    ------------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
 
/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map   \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

 
//    fprintf(OUTMAN, "         [-e Relaxation_parameter (in ]0,1] )]\n");
//    fprintf(OUTMAN, "             Default is 1.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptG = False;

    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "vzZ")) != -1) 
    {
	switch (c) 
        { 
    	    case 'v': Verbose = True; break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
      
        /* get optional input file names from trailing 
          parameters and open files */
    if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);
}
  
/*********************************************************************/

int main(int argc, char *argv[])
{
    int l,k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
    }
    
   Hfmap Map;
   Map.read(Name_Imag_In);
   //  float Min, Max;
   // Map.minmax(Min,Max);
   // Map.info((char*) "Input MAP");
   float *Dat = Map.buffer();
   int Np = Map.Npix();
   {
       double Mean,Sigma,Skew,Curt;
       float Min,Max;
       moment4(Dat, Np, Mean,  Sigma, Skew, Curt, Min,Max);
       double FluxIma = Mean*Np;
       float Energy = (float) (Sigma*Sigma+Mean*Mean)*Np;
       cout << "  Np = " << Np << endl;
       cout << "  Min = " << Min << "  Max = " << Max << endl;
       cout << "  Mean = " << Mean << "  Sigma = " << Sigma << endl;
       cout << "  Flux = " << FluxIma <<  "  Energy = " << Energy << endl;
  	   cout << "  Skew = " << Skew << "  Kurtosis = " << Curt << endl;
   }   
    exit(0);
}

/*********************************************************************/



