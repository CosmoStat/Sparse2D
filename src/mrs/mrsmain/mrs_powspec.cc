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
    fprintf(OUTMAN, "Usage: %s options in_map out_PowSpec_FileName \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 3*nside) \n", ALM_MAX_L);
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-n]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");

    fprintf(OUTMAN, "         [-m Mask_FileName]\n");
    fprintf(OUTMAN, "             Healpix Mask File Name. By default, no mask.\n");

    fprintf(OUTMAN, "         [-M MatrixMask_FileName]\n");
    fprintf(OUTMAN, "             Matrix related to the mask. By default, it is computed from the mask.\n");
    
    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Number of iteration (in case of a mask). Default is %d.\n", MaskIter);

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
           case 'n':
 	     	       NormALM=True;
 		          break;
           case 'P': InputPowSpec = True; break;
	   case 'p':
 	     	      if (sscanf(OptArg,"%f",&ZeroPadding) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad zero padding  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
	       case 'l':
 	     	      if (sscanf(OptArg,"%d",&Lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
	            	}
                break;
	    case 'i':
 	     	      if (sscanf(OptArg,"%d",&MaskIter) != 1) 
                  {
		            fprintf(OUTMAN, "Error: Iteration number : %s\n", OptArg);
		            exit(-1);
	            	}
                break;
            case  'm':
                if (sscanf(OptArg,"%s", Mask_FileName) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                OptMask = True;
                UseMask = True;
                break;
            case  'M':
                if (sscanf(OptArg,"%s", MaskMatrix_FileName) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                OptMat = True;
                UseMask = True;
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
        cout << "# File Name in = " << Name_Imag_In << endl;
		cout << "# Signal power spectrum file name = " <<  SignalPowSpec_FileName << endl;
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
		if (NormALM == True)  cout << "# ALM Normalization " << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
    }

   Hdmap Map, Mask;
   PowSpec SigPowSpec;
   int Lmax_Tot = 0;
   CAlmR  ALM;
   dblarray MatMask;     // Matrix related to the Mask
   int Nside = 0;
   
   if (OptMat == True) 
   {
        fits_read_dblarr(MaskMatrix_FileName, MatMask);
        Lmax_Tot = MatMask.nx()-1;
        if (Verbose == True)
        {
            MatMask.info("MatMask");
        }
   }
      
      
   if (InputPowSpec == False)
   {
      Map.read(Name_Imag_In);
      if (Verbose == True)
      {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP");
      }

      Nside = Map.Nside();
      if (Lmax_Tot <= 0) Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
      if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
   
      // ALM.Verbose = Verbose;
      ALM.Norm = NormALM;
      ALM.UseBeamEff = (ZeroPadding  > 0) ? True: False;

      // fits_write_fltarr("beameff.fits", ALM.BeamEff);
   
     if (UseMask == False) 
     {
        ALM.alloc(Nside, Lmax_Tot);
        if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);
        mrs_alloc_powspec(SigPowSpec, Lmax_Tot);
        ALM.alm_trans(Map);
        ALM.alm2powspec(SigPowSpec);
        if (Verbose == True)  cout << "ALM:  Lmax = " << ALM.Lmax() << " Mmax = " << ALM.Mmax() << endl;   
     }
      
   }
   else
   {
       mrs_read_powspec(Name_Imag_In,  SigPowSpec);
       Lmax_Tot = Lmax = SigPowSpec.Lmax();
   }

   
   if (UseMask == True)
   {
      PowSpec ClSol;

      if (OptMask == True)
      {
         Mask.read(Mask_FileName);
      }
      else
      {
          cout << "Error: Mask file name is missing.  " <<   endl;
          exit(-1);
      }
      if (InputPowSpec == False)
      {
         for (int p=0; p<  Map.Npix(); p++) Map[p] *= Mask[p];
         ALM.alloc(Nside, Lmax_Tot);
         if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);
         ALM.alm_trans(Map);
         mrs_alloc_powspec(SigPowSpec, Lmax_Tot);
         ALM.alm2powspec(SigPowSpec);
         if (Verbose == True)  cout << "ALM:  Lmax = " << ALM.Lmax() << " Mmax = " << ALM.Mmax() << endl;   	 
      }
      

       iter_master_decconv(SigPowSpec,  ClSol,  Mask, MatMask,  MaskIter,  Lmax_Tot,  Verbose);
       mrs_write_powspec(SignalPowSpec_FileName,  ClSol);    
   }
   else  mrs_write_powspec(SignalPowSpec_FileName,  SigPowSpec);    
   exit(0);
}

