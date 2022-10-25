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
**    File:  mrs_inpainting.cc 
**
*******************************************************************************
**
**    DESCRIPTION:  Inpainting + Iterative wiener filtering
**    -------------------
**
******************************************************************************/

#include"HealpixClass.h"
#include"MRS_Sparse.h"
#include"MRS_Inp.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Name_Mask[512]; 
char Name_RMSMap [512]; 
char Name_Beam [512]; 

char SignalPowSpec_FileName [512];
char NoisePowSpec_FileName [512];
char Name_PowSpec_Out [512];
char MaskMatrix_FileName[512];
char RevMaskMatrix_FileName[512];
char Name_Alm_Out[512];
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0;
#define  DEF_IWIENER_NBR_ITER  40
#define  DEF_ITER_GAUSSIANITY 10
 
int Max_Iter=DEF_IWIENER_NBR_ITER;
Bool UseRMSMap = False;   // Apply a denoising using an iterative Wiener filtering using a RMS map.
Bool UseMask = True;          // Apply an inpainting
Bool ThresholdResi = True;   // Threshold the Alm of the residual instead of the Alm of the solution
Bool WienerOK = False;        // Apply a denoising using an iterative Wiener filtering
Bool HardThreshold = True;  // Apply an Alm iterative hard threshold for the inpainting
Bool SoftThreshold = False;
Bool Use_ZeroMeanConstraint = False; // Apply a zero mean constraint on the solution
Bool Use_WP_Constraint = False;   // Apply a variance constraint in the wavelet packet of the solution
Bool EstimPowSpec = False;
Bool All_WP_Band = False;
Bool IsotropyCst = False;
Bool ProjectConstraint = False;
Bool LogNormal = False;

float  ZeroPadding=0.0;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;
Bool UpdatePowSpec = False;

float Eps=1.;
Bool GaussianityConstraint = False;
float GaussCst_ProjectCst_Lambda = 5.;
Bool Acceleration=False;
Bool Pos=False;

Bool Analysis = True;
Bool OptMat = False;
Bool OldAnalysis = False;
Bool OptRevMat = False;
Bool AlmWrite=False;

type_inpainting InpaintMethod = DEF_INP_METHOD;
Bool Denoising = False;
int AccelationLevel=1;
Bool UseBeam=False;

Bool Use_Realization_Powspec = False;
char ReaPowSpec_FileName [512];
int NbrMasterIter = DEF_NBR_MASTER_ITER;
int NbrScale=0;
Bool MCA_InitCleabCMB=False;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map in_mask out_map  [out_powspec]  [out_alm] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-n Number_of_Scales]\n");
    fprintf(OUTMAN, "             Number of scales in the wavelet transform.  \n");

    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
        
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Maximal multipole for spherical harmonic transform.  \n");
    fprintf(OUTMAN, "             Default is automatically estimated according to the input map nside: min(3*Nside,%d). \n", ALM_MAX_L);
    
    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Apply a wavelet filtering to identify bad point sources and update the mask.\n");

    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptI = False;
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "fn:s:l:i:vzZ")) != -1) 
    {
	switch (c) 
        {   
	    case 'f':MCA_InitCleabCMB = True; break;
            case 'n':
                if (sscanf(OptArg,"%d",&NbrScale) != 1)  
                {
		            fprintf(OUTMAN, "Error: bad number of scales value: %s\n", OptArg);
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
                   OptI = True;
                   if (Max_Iter <= 0)   Max_Iter =  DEF_IWIENER_NBR_ITER;
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

    if (OptInd < argc) strcpy(Name_Mask, argv[OptInd++]);
         else usage(argv);
           
	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) 
	{
		   strcpy(Name_PowSpec_Out, argv[OptInd++]);
		   EstimPowSpec = True;
	}
    if (OptInd < argc) 
	{
        strcpy(Name_Alm_Out, argv[OptInd++]);
        AlmWrite = True;
	}
    

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
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    
    if (Verbose == True)
    { 
        cout << "# COMMAND: " << Cmd << endl;
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# InpMap  File Name Out = " << Name_Imag_Out << endl;   
        if (EstimPowSpec ==True)  cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if ( UseMask ==  True) cout << "# File Name Mask = " <<   Name_Mask << endl;  
        cout << "# Number of iterations = " <<  Max_Iter << endl;   
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
         
#ifdef MRSDEBUG                      
        cout << "DEBUG MODE " << endl;
#endif
    }
    
    
   CInpainting  INP;
   INP.Map.read(Name_Imag_In);
   INP.MaskMap.read(Name_Mask);
   // if (Verbose == True)  INP.MaskMap.info((char*) "Input Mask");

   for (int p=0; p < INP.Map.Npix(); p++)  
   {
      if (INP.MaskMap[p] != 1) 
      {
        INP.MaskMap[p] = 0;
      }
   }
 
   int Nside = INP.Map.Nside();
   int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    
   mrs_alloc_powspec(INP.PowSpecData, Lmax_Tot);
       
   if (Verbose == True) INP.MaskMap.info((char*) "MaskMap");
   INP.Verbose = Verbose;
   INP.HardThreshold = HardThreshold;
   INP.SoftThreshold = SoftThreshold;
   INP.ThresholdResi = ThresholdResi;
   INP.SigmaNoise = SigmaNoise;
   INP.Use_ZeroMeanConstraint = Use_ZeroMeanConstraint;
   INP.Use_WP_Constraint = Use_WP_Constraint;
   INP.All_WP_Band = All_WP_Band;
   INP.IsotropyCst = IsotropyCst;
   INP.ProjectConstraint = ProjectConstraint;
   INP.ZeroPadding = ZeroPadding;
   INP.Lmax = Lmax;
   INP.Eps=Eps;
     INP.Max_Iter = Max_Iter;
   INP.Denoising = Denoising;
   INP.UseRMSMap = UseRMSMap;
   INP.WienerFilter = WienerOK;
   INP.LogNormal = LogNormal;
   INP.Use_Realization_Powspec = Use_Realization_Powspec;
   INP.NbrMasterIter = NbrMasterIter;
   INP.InpaintMethod = InpaintMethod;
INP.Verbose= Verbose;

   Analysis = True;
   INP.Use_WP_Constraint = False;
   INP.MCA = True;
   INP.MCA_WT_NbrScale = NbrScale;
   INP.MCA_InitCleabCMB = MCA_InitCleabCMB;
   INP.analyse_alm_inpainting();
              
   INP.Result.write(Name_Imag_Out);

   char FN[512]; 
   sprintf(FN, "%s_alm.fits", Name_Imag_Out);
   ((INP.MCA_TabResult)[0]).write(FN);
   sprintf(FN, "%s_wt.fits", Name_Imag_Out);
   ((INP.MCA_TabResult)[1]).write(FN);
   if (MCA_InitCleabCMB == True)
   {
      sprintf(FN, "%s_mask.fits", Name_Imag_Out);
      INP.MaskMap.write(FN);
   }
   INP.ALM.Norm = False;
   if ((Analysis == True) && ((AlmWrite == True) || (EstimPowSpec == True))) 
   {
    // cout << "ALM ALM " << endl;
    INP.ALM.alm_trans((INP.MCA_TabResult)[0]);
   }
   
   if (AlmWrite == True)  INP.ALM.write(Name_Alm_Out);
   
   if (EstimPowSpec == True) 
   {
   // cout << "EstimPowSpec " << endl;
      PowSpec SigPowSpec(1, Lmax_Tot);
      INP.ALM.alm2powspec(SigPowSpec);
      mrs_write_powspec( Name_PowSpec_Out,  SigPowSpec);    
   }
   exit(0);
}

