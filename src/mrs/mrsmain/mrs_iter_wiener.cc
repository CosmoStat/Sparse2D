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
**    File:  mrs_iter_wiener.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Iterative wiener filtering
**    -------------------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Name_Mask[512]; 
char Name_RMSMap [512]; 
char SignalPowSpec_FileName [512];
char NoisePowSpec_FileName [512];

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
 #define  DEF_IWIENER_NBR_ITER  40
 int Max_Iter=DEF_IWIENER_NBR_ITER;
 
Bool UseRMSMap = False;
Bool UseMask = False;

float  ZeroPadding=0.1;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;

float Eps=1.;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map  out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
     
    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
    
    fprintf(OUTMAN, "         [-g SigmaNoise]\n");
    fprintf(OUTMAN, "             Gaussian Noise standard deviation. Default is 0.\n");
    
    fprintf(OUTMAN, "         [-N Noise_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Noise Power Spectrum file name . \n");

    fprintf(OUTMAN, "         [-S Signal_PowSpec_FileName]\n");
    fprintf(OUTMAN, "             Signal Power Spectrum file name . \n");

    fprintf(OUTMAN, "         [-R Noise_RmsMap]\n");
    fprintf(OUTMAN, "             Noise RMS map file name . Default is Gaussian noise\n");
    
    fprintf(OUTMAN, "         [-M MaskMap]\n");
    fprintf(OUTMAN, "             Mask map file name. Default is no mask. \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.1.\n");

    fprintf(OUTMAN, "         [-e Relaxation_parameter (in ]0,1] )]\n");
    fprintf(OUTMAN, "             Default is 1.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptG = False;

    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "e:N:S:l:p:R:M:g:i:vzZ")) != -1) 
    {
	switch (c) 
        { 
           case 'e':
 	     	      if (sscanf(OptArg,"%f",&Eps) != 1) 
                  {
		            fprintf(OUTMAN, "Error: relaxation parameter  value: %s\n", OptArg);
		            exit(-1);
	            	}
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
            case 'i':  
                   if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                   {
		            fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
		    		exit(-1);
		           } 
                   if (Max_Iter <= 0)   Max_Iter =  DEF_IWIENER_NBR_ITER;
 	        	break;
	        case 'g':
 	     	      if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		            exit(-1);
	            	}
	            	OptG = True;
 		          break;
    	    case 'R':
     		       if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                UseRMSMap = True;
  	        	break;
  	        case  'M':
  	              if (sscanf(OptArg,"%s", Name_Mask) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                UseMask = True;
  	        	break;
  	       case  'N':
  	              if (sscanf(OptArg,"%s", NoisePowSpec_FileName) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                OptN_PS = True;
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
	
	    if ((UseMask == False) && (OptG == False) && (UseRMSMap == False))
       {
     	  cerr << endl << endl;
          fprintf(OUTMAN, "Error: at least on of the option -g, -R or -M have to be set ...\n");
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
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (UseRMSMap ==  True) cout << "# File Name Noise RMS = " <<   Name_RMSMap << endl;  
        else if (SigmaNoise > 0)  cout << "Sigma Noise = " <<   SigmaNoise << endl; 
        if ( UseMask ==  True) cout << "# File Name Mask = " <<   Name_Mask << endl;  
        cout << "# Number of iterations = " <<  Max_Iter << endl;   
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
   	   if (Verbose == True) Map.info((char*) "Input MAP");
   }

   if (UseRMSMap ==  True) 
   {
   	   RMSmap.read(Name_RMSMap);
   	   if (Verbose == True) RMSmap.info((char*) "RMSmap");
    }
   
    if (UseMask  ==  True) 
   {
   	   MaskMap.read(Name_Imag_In);
   	   if (Verbose == True) MaskMap.info((char*) "MaskMap");
   }

   int Nside = Map.Nside();
   if (Lmax <= 0) Lmax = 3 * Nside;
   int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
   
   CAlmR  ALM;  
   // ALM.Verbose = Verbose;
   ALM.alloc(Nside, Lmax_Tot);
   ALM.Norm = True;
    if (Verbose == True)  cout << "ALM class:  Lmax = " << ALM.Lmax() << " Mmax = " << ALM.Mmax() << endl;
   
   for (int l=0; l < ALM.Lmax(); l++)
   for (int m=0; l < ALM.Mmax(); l++) ALM(l,m) = 0.;
   ALM.set_beam_eff(Lmax, Lmax_Tot);
   ALM.UseBeamEff = True;
   if (Verbose == True)      fits_write_fltarr((char *)"xx_beameff.fits", ALM.BeamEff);
   
   // Get the noise power spectrum.
   PowSpec NoisePowSpec(1, Lmax_Tot);
   if (SigmaNoise > 0)  // case of White Gaussian noise
   {
   	    double V2 = SigmaNoise*SigmaNoise;
       	if (OptN_PS ==  False)  for (l=0; l <= NoisePowSpec.Lmax(); l++)   NoisePowSpec.tt(l) =  V2;
   }
   if (OptN_PS == True)
   {
   	   mrs_read_powspec(NoisePowSpec_FileName , NoisePowSpec);  
   	   for (l=0; l <= NoisePowSpec.Lmax(); l++)  NoisePowSpec.tt(l) *= ALM.NormVal;
   }
    if (UseRMSMap ==  True) 
   {
   	   // void init_random (unsigned int Init);
   	   Hdmap MapNoise;
   	   MapNoise.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
   	    for (int p=0; p < MapNoise.Npix() ; p++)  MapNoise[p] = get_random() * RMSmap[p];
        if (Verbose == True)      MapNoise .write((char *)"xx_mapnoise.fits");
        ALM.alm_trans(MapNoise);
        ALM.alm2powspec(NoisePowSpec);
   	}
   if (Verbose == True)  mrs_write_powspec((char *)"xx_psnoise", NoisePowSpec);    
     
    //Signal noise power spectrum.
    PowSpec SigPowSpec(1, Lmax_Tot);
    ALM.alm_trans(Map);
    if (OptS_PS == True)
    {
   	   mrs_read_powspec(SignalPowSpec_FileName, SigPowSpec);
       for (l=0; l <= SigPowSpec.Lmax(); l++)  SigPowSpec.tt(l) *= ALM.NormVal;
    }
    else  
    {
       ALM.alm2powspec(SigPowSpec);
       for (int l=0; l < SigPowSpec.Lmax(); l++)  SigPowSpec.tt(l) = (SigPowSpec.tt(l) > NoisePowSpec.tt(l) ) ? (SigPowSpec.tt(l) - NoisePowSpec.tt(l) ) : 0.;
   	}
   	if (Verbose == True)  mrs_write_powspec((char *)"xx_pssig",  SigPowSpec);    

    Hdmap Result;
    Hdmap Resi;
   	Result.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
    Resi = Map;
   
   double MaxAlm =ALM.max_absalm();
   if (Verbose == True)  cout << "==> Max Alm = " << MaxAlm << endl;
   
  double MinRMS,MaxRMS;
  if (UseRMSMap ==  True)  
  {
       RMSmap.minmax(MinRMS,MaxRMS);
  }
  else MinRMS = MaxRMS = SigmaNoise;
   if (Verbose == True)  cout << "MinRMS  " <<  MinRMS <<   ",  MaxRMS = " << MaxRMS << endl;
  
  double FirstSoftThreshold = MaxAlm;
  double LastSoftThreshold = 0.;
  double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
  // double StepL = DeltaThreshold / (float) (Max_Iter -1);
  double Lambda = FirstSoftThreshold;
   
   PowSpec Iter_NoisePowSpec(1, Lmax_Tot);
   double Err=Map.rms();
   int NbrTotMask = 0;
   if (UseMask  ==  True)  for (int p=0; p<  Result.Npix(); p++)  if (MaskMap[p] > 0)  NbrTotMask ++;
   for (int i=0; i < Max_Iter; i++)
   { 
       char NN[256];
	   
      	  Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(4.*i/ (Max_Iter-1)));
   	      if (Verbose == True)  cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
   	      
   	      for (int p=0; p<  Result.Npix(); p++) 
   	      {
   	      	   // Compute the residual
   	           Resi[p] =  Map[p] - Result[p];
   	           if (UseRMSMap ==  True)  Resi[p] *=  MinRMS / RMSmap[p];
    	       if (UseMask  ==  True)  Resi[p] *= MaskMap[p];
    	       Err += Resi[p]*Resi[p];
   	           Resi[p] += Result[p];
   	       }
           if (UseMask  ==  True)   Err /= (float) (Result.Npix() - (Result.Npix() - NbrTotMask));
		   else  Err /= (float) Result.Npix();
           Err = sqrt(Err) / MinRMS;
            
           // Update Noise power spectrum
           for (l=0; l <= NoisePowSpec.Lmax(); l++)   Iter_NoisePowSpec.tt(l)  =  MinRMS*MinRMS  * (1 + Lambda);
             
           // Wiener Filtering 
		   Resi.info((char*) "BefW MAP");
           ALM.alm_trans(Resi);
   	       ALM.wiener(Iter_NoisePowSpec, SigPowSpec);
   	       ALM.alm_rec(Resi);
           Resi.info((char*) "AftW MAP");
		   ALM.WienerFilter.info("WF: ");
		   fits_write_fltarr(fitsname("xx_wien"), ALM.WienerFilter);
           sprintf(NN, "xx_it_%d.fits", i+1);
           // Update the solution
           for (int p=0; p<  Result.Npix(); p++)  Result[p] = (1-Eps)*Result[p] + Eps*Resi[p];
		   Result.write(NN);
    }
   
   Result.write(Name_Imag_Out);
    exit(0);
}
