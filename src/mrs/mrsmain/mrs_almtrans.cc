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
// #include "cxxutils.h"
#include "datatypes.h"
// #include "openmp_support.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char SignalPowSpec_FileName [512];
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
 
Bool NormALM = False;
float  ZeroPadding=0.0;
int Lmax = 0.;
Bool OptS_PS = False;
Bool CFIMA= False;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-S Output_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Signal Power Spectrum file name . \n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    
//    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
//    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");
	
    fprintf(OUTMAN, "         [-T]\n");
    fprintf(OUTMAN, "             Write the output as a fits [*,*,2] array image.\n");
    fprintf(OUTMAN, "             with [0..Lmax, 0..Mmax, 0] for the real part and \n");
    fprintf(OUTMAN, "             with [0..Lmax, 0..Mmax, 1] for the imaginary part and \n");
    fprintf(OUTMAN, "             By default, it must be read in IDL with :  a = mrdfits(FileName, 1).\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "TNS:l:p:vzZ")) != -1) 
    {
	switch (c) 
        { 
		    case 'T':
 	     	       CFIMA=True;
 		          break;
           case 'N':
 	     	       NormALM=True;
 		          break;
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
 	     	case  'S':
  	              if (sscanf(OptArg,"%s", SignalPowSpec_FileName) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                OptS_PS = True;
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
/*
static void openmp_status()
  {
  if (openmp_enabled())
    {
    cout << "Application was compiled with OpenMP support," << endl;
    if (openmp_max_threads() == 1)
      cout << "but running with one process only." << endl;
    else
      cout << "running with up to " << openmp_max_threads()
           << " processes." << endl;
    }
  else
    cout << "Application was compiled without OpenMP support;" << endl
         << "running in scalar mode." << endl;
  }
 */
int main(int argc, char *argv[])
{
// openmp_status();

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
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
		if (NormALM == True)  cout << "# ALM Normalization " << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
		if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  SignalPowSpec_FileName << endl;
    }

   Hdmap Map;
   Hdmap RMSmap;
   Hdmap MaskMap;
    if (Verbose == True)  cout << " ALM read ... " << endl;
    
   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP");
   }

   int Nside = Map.Nside();
   int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);

   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
   
   CAlmR  ALM;
   // ALM.Verbose = Verbose;
   // if (Verbose == True)  cout << " ALM Alloc ... " << endl;
   ALM.Norm = NormALM;
   ALM.UseBeamEff = False; 
   ALM.alloc(Nside, Lmax_Tot);
   ALM.set_beam_eff(Lmax, Lmax_Tot);
   // fits_write_fltarr("beameff.fits", ALM.BeamEff);
  //  if (Verbose == True)  cout << " ALM transform ... " << endl;
   ALM.alm_trans(Map);

   if (Verbose == True)  cout << "ALM class:  Lmax = " << ALM.Lmax() << " Mmax = " << ALM.Mmax() << " MaxAbs = " <<  ALM.max_absalm() << endl;

   if (CFIMA == False)
   {
      ALM.write(Name_Imag_Out);
    // read the resultat in IDL:  a = mrdfits(FileName, 1)
   }
   else
   {
   
       dblarray A;
       A.alloc(ALM.Lmax()+1., ALM.Mmax()+1, 2);
       for (int l=0; l <= ALM.Lmax(); l++)
       for (int m=0; m <= l; m++) 
	   {
	       A(l,m,0) = ALM(l,m).real();
	       A(l,m,1) = ALM(l,m).imag();
	       // if ((l == 192) && (m == 10)) cout <<  A(l,m,0) << " " <<  A(l,m,1) << endl;
	   }
	 //  cout << "FITS W " << endl;
       fits_write_dblarr(Name_Imag_Out, A);
        //  cout << "FITS W " << endl;
   }
   
   if (OptS_PS == True)
   {
       PowSpec SigPowSpec(1, Lmax_Tot);
       ALM.alm2powspec(SigPowSpec);
   	   mrs_write_powspec( SignalPowSpec_FileName,  SigPowSpec);    
   }

   exit(0);
}

