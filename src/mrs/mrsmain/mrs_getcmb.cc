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
**    Date:  28/11/08
**    
**    File:  mrs_getcmb.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Get CMB realizations using a given cosmological model    
**    -----------
**
******************************************************************************/

#include"HealpixClass.h"
 
char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Cl_FileName [512];
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
 
Bool NormALM = False;
float  ZeroPadding=0.;
int Lmax = 0.;
 
Bool OptAllFile = True;
Bool OptCL = False;
fltarray Cl;
int Nside=256;
int NbrCMBMap=1;

#define DEF_RND_VAL 100
unsigned int InitRnd = DEF_RND_VAL;

float Fwhm_arcmin=5.;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-C Cl_FileName]\n");
    fprintf(OUTMAN, "             True Cl of the CMB. Default is \n");

    fprintf(OUTMAN, "         [-N NbrNaps]\n");
    fprintf(OUTMAN, "             Number of realization. Default is 1.\n");
    
    fprintf(OUTMAN, "         [-f Fwhm_arcmin]\n");
    fprintf(OUTMAN, "             Full Width at Half Maximum. Default is 5.\n");
    
    fprintf(OUTMAN, "         [-n Nside]\n");
    fprintf(OUTMAN, "             Nside value. Default is 256.\n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");
    
    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "             Default is 100. \n");
    
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "             Only one 2D output fits files (ring format). default is no.\n");
    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "aI:N:n:C:f:l:p:vzZ")) != -1) 
    {
	switch (c) 
        { 
        case 'a': OptAllFile = (OptAllFile == True) ? False: True;
                  break;
        case 'I':
 		        if (sscanf(OptArg,"%d",&c) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		            exit(-1);
		        }
                InitRnd = (unsigned int) c;
		        break;
         case 'N':
 	     	      if (sscanf(OptArg,"%d",& NbrCMBMap) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad number of maps: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
          case 'n':
 	     	      if (sscanf(OptArg,"%d",&Nside) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad nside  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;         
          case 'f':  
                  if (sscanf(OptArg,"%f",& Fwhm_arcmin) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad Fwhm value: %s\n", OptArg);
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
 	     	case  'C':
  	              if (sscanf(OptArg,"%s", Cl_FileName) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		           }
                OptCL  = True;
  	        	break;
  	    case 'v': Verbose = True; break;
        case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 

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
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
		if (OptCL == True)  cout << "# Cl file name = " <<  Cl_FileName << endl;
		if (NbrCMBMap > 1) cout << "# Number of realizations = " <<  NbrCMBMap << endl;
	    cout << "# Nside = " <<  Nside << endl;
	    if (InitRnd != DEF_RND_VAL)  cout << "# RND init = " <<  InitRnd << endl;
	    cout << "# Fwhm = " <<  Fwhm_arcmin << endl; 
    }
    
   
   Hdmap Map;
   Map.alloc(Nside); 
   int Nside = Map.Nside();
   if (Lmax <= 0) Lmax = 3 * Nside;
   int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
   
   PowSpec Cl(1, Lmax_Tot);
   if (OptCL == False)  
   {
   	  char *FN = getenv("ISAP");
   	sprintf(Cl_FileName, "%s/data/%s", FN, DEF_MRS_CL_LCDM);
   	 if (Verbose == True) cout << "==> Input Cl Filename = " << Cl_FileName << endl;
   }
   mrs_read_powspec(Cl_FileName, Cl);
   Fwhm_arcmin *= degr2rad/60;
   if (Fwhm_arcmin  > 0) Cl.smoothWithGauss(Fwhm_arcmin);

   CAlmR  ALM;
   // ALM.Verbose = Verbose;
   ALM.Norm = False; 
   ALM.UseBeamEff = (ZeroPadding  > 0) ? True: False;
   ALM.alloc(Nside, Lmax_Tot);
   if (ALM.UseBeamEff == True) ALM.set_beam_eff(Lmax, Lmax_Tot);

   
   planck_rng rng((int) InitRnd);
   
   if (NbrCMBMap == 1)
   {
      create_alm (Cl, ALM, rng);
      ALM.alm_rec(Map);
      if (Verbose == True)
      {
   	   if (Verbose == True) Map.info((char*) "Output MAP");
      }
      Map.write(Name_Imag_Out);
   }
   else
   {
   	  if (OptAllFile == True)
   	  {
   	  	  for (int i=0; i < NbrCMBMap; i++)
   	  	  {
   	         if (Verbose == True) cout << "==> CMB " << i+1 << endl;
   	  	  	 create_alm (Cl, ALM, rng);
             ALM.alm_rec(Map);
             char FN[512];
             sprintf(FN, "%s_%d", Name_Imag_Out,i+1);
             Map.write(FN);
   	  	  }
   	  }
   	  else
   	  {
   	  	 fltarray CMBTab;
   	  	 CMBTab.alloc(Map.Npix(), NbrCMBMap);
    	 for (int i=0; i < NbrCMBMap; i++)
   	  	 {
   	  	  	 create_alm (Cl, ALM, rng);
             ALM.alm_rec(Map);
             for (int p=0; p < Map.Npix(); p++) CMBTab(p,i) = Map[i];
   	  	 }
   	  	 fits_write_fltarr(Name_Imag_Out, CMBTab);
   	  }
   }
   exit(0);
}

