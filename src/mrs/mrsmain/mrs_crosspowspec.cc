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
**    DESCRIPTION:  ALM transformation
**    -------------------
**
******************************************************************************/

#include"HealpixClass.h"
#include"MatMask.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_In2[512]; /* input file image */
char SignalPowSpec_FileName [512];
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
 
Bool NormALM = False;
float  ZeroPadding=0.;
int Lmax = 0.;
Bool OptS_PS = False;
Bool CFIMA= False;

char MaskMatrix_FileName[512];
char Mask_FileName[512];
Bool UseMask = False;         
int MaskIter = 20;

Bool OptMat = False;
Bool OptMask = False;
Bool InputPowSpec = False;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map_1 in_map_2  out_CrossSpec_FileName \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 3*nside) \n", ALM_MAX_L);
    

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "Pi:m:M:l:p:vzZ")) != -1) 
    {
	switch (c) 
        { 
 	       case 'l':
 	     	      if (sscanf(OptArg,"%d",&Lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
	            	}
                break;
             case 'v': Verbose = True; break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 

        /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);
       if (OptInd < argc) strcpy(Name_Imag_In2, argv[OptInd++]);
       else usage(argv);
	  if (OptInd < argc) strcpy(SignalPowSpec_FileName, argv[OptInd++]);
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
    double Min, Max;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in 1 = " << Name_Imag_In << endl;
        cout << "# File Name in 2 = " << Name_Imag_In2 << endl;
		cout << "# Signal cross spectrum file name = " <<  SignalPowSpec_FileName << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
    }

   Hdmap Map, Map2;
   PowSpec SigPowSpec;
   int Lmax_Tot = 0;
   CAlmR  ALM, ALM2;
   dblarray MatMask;     // Matrix related to the Mask
   int Nside = 0;
   
        
   
   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP 1");
    }
    Map2.read(Name_Imag_In2);
    if (Verbose == True)
    {
        Map2.minmax(Min,Max);
        if (Verbose == True) Map2.info((char*) "Input MAP 2");
    }

    Nside = Map.Nside();
    if (Lmax_Tot <= 0) Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    if (Verbose == True) cout << "==> Use Lmax = " << Lmax_Tot << endl;
   
      // ALM.Verbose = Verbose;
    ALM.Norm = NormALM;
    ALM.UseBeamEff = (ZeroPadding  > 0) ? True: False;
    ALM2.Norm = ALM.Norm;
    ALM2.UseBeamEff = ALM.UseBeamEff;

      // fits_write_fltarr("beameff.fits", ALM.BeamEff);
   
    ALM.alloc(Nside, Lmax_Tot);
    if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);
    ALM2.alloc(Nside, Lmax_Tot);
    if (ALM2.UseBeamEff == True) ALM2.set_beam_eff(Lmax, Lmax_Tot);

    ALM.alm_trans(Map);
    ALM2.alm_trans(Map2);
    extract_crosspowspec(ALM, ALM2, SigPowSpec);

    if (Verbose == True)  cout << "ALM:  Lmax = " << ALM.Lmax() << " Mmax = " << ALM.Mmax() << endl;   
    mrs_write_powspec(SignalPowSpec_FileName,  SigPowSpec);    
   exit(0);
}

