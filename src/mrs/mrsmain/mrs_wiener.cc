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
**    File:  mrs_wiener.cc
**
*******************************************************************************
**
**    DESCRIPTION: Wiener filtering
**    ------------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Name_Mask[512]; 
char Name_RMSMap [512]; 
char SignalPowSpec_FileName [512];
char NoisePowSpec_FileName [512];
char Name_PowSpec_Out [512];

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
 #define  DEF_IWIENER_NBR_ITER  40
 int Max_Iter=DEF_IWIENER_NBR_ITER;
 
Bool UseRMSMap = False;   // Apply a denoising using an iterative Wiener filtering using a RMS map.
Bool WienerOK = False;        // Apply a denoising using an iterative Wiener filtering
Bool EstimPowSpec = False;
Bool All_WP_Band = False;

float  ZeroPadding=0.;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;

char WienerFilterName[512];
Bool WriteWienerFilter=False;
Bool UseNorm = False;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map out_map  [out_powspec]  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-g SigmaNoise]\n");
    fprintf(OUTMAN, "             Gaussian Noise standard deviation. Default is 0.\n");
    
    fprintf(OUTMAN, "         [-N Noise_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Input Noise Power Spectrum file name. Default is automatically estimated.\n");

    fprintf(OUTMAN, "         [-S Signal_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Input Signal Power Spectrum file name (for Wiener filtering). Default is automatically estimated. \n");

    fprintf(OUTMAN, "         [-R Noise_RmsMap]\n");
    fprintf(OUTMAN, "             Noise RMS map file name. Default is Gaussian noise\n");
         
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");
    fprintf(OUTMAN, "         [-w WienerFilter_FileName]\n");
    fprintf(OUTMAN, "             Write the Wiener filter to the disk.  Default is no.\n");
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");

         
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptG = False;

    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "tw:N:S:l:p:R:g:vzZ")) != -1) 
    {
	switch (c) 
        { 
        	  case 't': (UseNorm == True) ? False: True; break;
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
 	        case 'g':
 	     	      if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		            exit(-1);
	            	}
	            	OptG = True;
  		          break;
               case 'w':
     		      if (sscanf(OptArg,"%s", WienerFilterName) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
		           WriteWienerFilter = True;
   	        	break;
     	    case 'R':
     		       if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                UseRMSMap = True;
   	        	break;
   	       case  'N':
  	              if (sscanf(OptArg,"%s", NoisePowSpec_FileName) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                OptN_PS = True;
                WienerOK=True;
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
	
     if  ((OptG == False) && (UseRMSMap == False))
    {
     	  cerr << endl << endl;
          fprintf(OUTMAN, "Error: option -g and -R are not valid together. ...\n");
           exit(-1);
  	}

}
  
/*********************************************************************/

int main(int argc, char *argv[])
{
    int l,k;
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
        cout << "#  InpMap  File Name Out = " << Name_Imag_Out << endl;   
        if (UseRMSMap ==  True) cout << "# File Name Noise RMS = " <<   Name_RMSMap << endl;  
        else if (SigmaNoise > 0)  cout << "Sigma Noise = " <<   SigmaNoise << endl; 
         if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
        if (OptN_PS > 0)  cout << "# Noise power spectrum file name = " <<  NoisePowSpec_FileName << endl;
        if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  NoisePowSpec_FileName << endl;
    }
    
   Hdmap Map;
   Hdmap RMSmap;
   Hdmap MaskMap;

   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   Map.info((char*) "Input MAP");
   }
   
   if (UseRMSMap ==  True) 
   {
   	   RMSmap.read(Name_RMSMap);
   	   if (Verbose == True) RMSmap.info((char*) "RMSmap");
    }
       
   int Nside = Map.Nside();
   if (Lmax <= 0) Lmax = 3 * Nside;
   int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;

   CAlmR  ALM;  
   ALM.alloc(Nside, Lmax_Tot);
   ALM.Norm = UseNorm;
   ALM.set_beam_eff(Lmax, Lmax_Tot);
   if (ZeroPadding  > 0) ALM.UseBeamEff = True;
     
   // Get the noise power spectrum.
   PowSpec NoisePowSpec(1, Lmax_Tot);
   if (SigmaNoise > 0)  // case of White Gaussian noise
   {
   	   if (ALM.Norm == False)  SigmaNoise /= ALM.NormVal;
   	   double V2 = SigmaNoise*SigmaNoise;
       if (OptN_PS ==  False)  for (l=0; l <= NoisePowSpec.Lmax(); l++)   NoisePowSpec.tt(l) =  V2;
   }
   if (OptN_PS == True)
   {
   	   mrs_read_powspec(NoisePowSpec_FileName, NoisePowSpec);  
   	   if (ALM.Norm == True) 
   	        for (l=0; l <= NoisePowSpec.Lmax(); l++)  NoisePowSpec.tt(l) *= (ALM.NormVal*ALM.NormVal);
   }
   // cout << "==> NormVal = " << ALM.NormVal << endl;
   if (UseRMSMap ==  True) 
   {
   	   // void init_random (unsigned int Init);
   	   Hdmap MapNoise;
   	   MapNoise.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
           	   if (Verbose == True) RMSmap.info((char*) "RMSmap");
   	    for (int p=0; p < MapNoise.Npix() ; p++)  MapNoise[p] = get_random() * RMSmap[p];
       // if (Verbose == True) MapNoise.write((char *)"xx_mapnoise.fits");
       // if (Verbose == True) MapNoise.info((char*) "MapNoise");

        ALM.alm_trans(MapNoise);
        ALM.alm2powspec(NoisePowSpec);
   	}
    // if (Verbose == True)  mrs_write_powspec((char *)"xx_psnoise", NoisePowSpec);    
     
    //Signal noise power spectrum.
    PowSpec SigPowSpec(1, Lmax_Tot);
    ALM.alm_trans(Map);
    if (OptS_PS == True)
    {
   	   mrs_read_powspec(SignalPowSpec_FileName, SigPowSpec);
   	   if (ALM.Norm == True) 
   	        for (l=0; l <= SigPowSpec.Lmax(); l++)  SigPowSpec.tt(l) *= (ALM.NormVal*ALM.NormVal);
    }
    else  
    { 
       ALM.alm2powspec(SigPowSpec);
       for (int l=0; l < SigPowSpec.Lmax(); l++)  
       SigPowSpec.tt(l) = (SigPowSpec.tt(l) > NoisePowSpec.tt(l) ) ? (SigPowSpec.tt(l) - NoisePowSpec.tt(l) ) : 0.;
   	}
   	if (Verbose == True)  mrs_write_powspec((char *)"xx_pssig",  SigPowSpec);    

   
   double MaxAlm =ALM.max_absalm();
   if (Verbose == True)  cout << "==> Max Alm = " << MaxAlm << endl;
                 
   ALM.wiener(NoisePowSpec, SigPowSpec);
   ALM.alm_rec(Map);
   Map.write(Name_Imag_Out);
        if (Verbose == True)  Map.info( (char *) "  ==> Wiener filtering map ");

   if (WriteWienerFilter == True)
         fits_write_fltarr(WienerFilterName, ALM.WienerFilter);
   exit(0);
}
