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
Bool ThresholdResi = False;   // Threshold the Alm of the residual instead of the Alm of the solution
Bool WienerOK = False;        // Apply a denoising using an iterative Wiener filtering
Bool HardThreshold = True;  // Apply an Alm iterative hard threshold for the inpainting
Bool SoftThreshold = False;
Bool Use_ZeroMeanConstraint = True; // Apply a zero mean constraint on the solution
Bool Use_WP_Constraint = True;   // Apply a variance constraint in the wavelet packet of the solution
Bool EstimPowSpec = False;
Bool All_WP_Band = False;
Bool IsotropyCst = False;
Bool ProjectConstraint = False;
Bool LogNormal = False;

float  ZeroPadding=0.0;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;
Bool UpdatePowSpec = True;

float Eps=1.;
Bool GaussianityConstraint = False;
float GaussCst_ProjectCst_Lambda = 5.;
Bool Acceleration=False;
Bool Pos=False;

Bool Analysis = False;
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

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map in_mask out_map  [out_powspec]  [out_alm] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-m Inpainting method]\n");
    for (int m=0; m < NBR_INP_METHOD; m ++)
        fprintf(OUTMAN, "             %d: %s.\n", m+1, (char *) StringInpaintingMethod( (type_inpainting) m)); 
    fprintf(OUTMAN, "             default is %s\n", StringInpaintingMethod (InpaintMethod));

    
    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
        
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
    
    fprintf(OUTMAN, "         [-A AccelerationLevel]\n");
    fprintf(OUTMAN, "             Only valid with -m1 method. If set, the WP constraint is applied only once in N iterations, with N=AccelerationLevel. \n");
    fprintf(OUTMAN, "             When -A 0 is set, no WP constraint is used.  \n");
    
    fprintf(OUTMAN, "         [-g GaussianNoiseStd]\n");
    fprintf(OUTMAN, "             Apply also a denoising. GaussianNoiseStd is the noise standard deviation. \n");
    fprintf(OUTMAN, "             Default is no denoising. \n");
    fprintf(OUTMAN, "         [-r RMS_MapFileName]\n");
    fprintf(OUTMAN, "             Apply also a denoising. RMS_MapFileName is the filename of the map containing the RMS per pixel. \n");
    fprintf(OUTMAN, "             Default is no denoising. \n");
    fprintf(OUTMAN, "         [-S Prior_SignalPowspec]\n");
    fprintf(OUTMAN, "             Prior signal power spectrum for Wiener or the Cole denoising (i.e. Theoretical Cl).  Valid only if -g or -r is set (i.e. denoising mode).\n");
    fprintf(OUTMAN, "         [-s Prior_Realization_Powspec]\n");
    fprintf(OUTMAN, "             Prior on the output power spectrum. Constraint the solution to have a given power spectrum.\n");
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "            Replace the Cole denoising by the Wiener denoising. Valid only if -g or -r is set (i.e. denoising mode). \n");
    fprintf(OUTMAN, "         [-W]\n");
    fprintf(OUTMAN, "            Remove the variance constraint. \n");
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "            Remove the zero mean constraint. \n");

    
//    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
//    fprintf(OUTMAN, "             Default is %f.\n",ZeroPadding);

//    fprintf(OUTMAN, "         [-p]\n");
//    fprintf(OUTMAN, "             Add a positivity constraint.\n");
    
//    fprintf(OUTMAN, "         [-e Relaxation_parameter (in ]0,1] )]\n");
//    fprintf(OUTMAN, "             Default is 1.\n");
    
//    fprintf(OUTMAN, "         [-G]\n");
//    fprintf(OUTMAN, "             Add a Gaussianity constraint.\n");
    
//    fprintf(OUTMAN, "         [-C]\n");
//    fprintf(OUTMAN, "             Remove the zero mean constraint.\n");
    
//    fprintf(OUTMAN, "         [-W]\n");
//    fprintf(OUTMAN, "             Constraint on the variance of the wavelet packet decomposition.\n");
    
//    fprintf(OUTMAN, "         [-P GaussCst_ProjectCst_Lambda]\n");
//    fprintf(OUTMAN, "             Parameter for the Gaussianity constraint. Default is 5.\n");
    
//    fprintf(OUTMAN, "         [-I]\n");
//    fprintf(OUTMAN, "             Replace the L1 Alm minimization by an isotropic constraint.\n");
    
    fprintf(OUTMAN, "         [-M MatrixMask_FileName]\n");
    fprintf(OUTMAN, "             Matrix related to the mask. By default, it is computed from the mask.\n");

    fprintf(OUTMAN, "         [-R RevMatrixMask_FileName]\n");
    fprintf(OUTMAN, "             Reverse Matrix related to the mask. By default, it is computed from the mask.\n");
    
//    fprintf(OUTMAN, "         [-H]\n");
//    fprintf(OUTMAN, "             Replace the hard-thresholding by a soft-thresholding.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptI = False;
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "yNTwb:g:r:m:oM:rR:s:BS:ICWHGP:A:e:l:pg:i:vzZ")) != -1) 
    {
	switch (c) 
        {   
            case 'N': LogNormal = (LogNormal == True) ? False: True;  break;
            case 'w': WienerOK =  (WienerOK == True) ? False: True;  break;
	        case 'A':  if (sscanf(OptArg,"%d",&AccelationLevel) != 1) 
                       {
                          fprintf(OUTMAN, "Error: bad AccelationLevel: %s\n", OptArg);
                          exit(-1);
                       } 
                      Acceleration = True;  
                      break;
             case 'm':
				if (sscanf(OptArg,"%d",&c ) != 1) 
				{
					fprintf(OUTMAN, "Error: bad type of inpainting method: %s\n", OptArg);
					exit(-1);
 				}
                if ((c > 0) && (c <= NBR_INP_METHOD))  InpaintMethod = (type_inpainting) (c-1);
                else  
                {
					fprintf(OUTMAN, "Error: bad type of inpainting method: %s\n", OptArg);
					exit(-1);
				}
 				break;
            case 'o': OldAnalysis = True;
                break;
            case  'M':
                if (sscanf(OptArg,"%s", MaskMatrix_FileName) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                OptMat = True;
                break;
            case  'r':
                if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                Denoising = True;
                UseRMSMap = True;
                break;
            case  'b':
                if (sscanf(OptArg,"%s", Name_Beam) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                UseBeam = True;
                break;
            case  'R':
                if (sscanf(OptArg,"%s", RevMaskMatrix_FileName) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                OptRevMat = True;
                break;
            case 'T': ThresholdResi = (ThresholdResi == True) ? False: True;  break;
            case 'y': Analysis = (Analysis == True) ? False: True;  break;
            case 'B': ProjectConstraint = (ProjectConstraint == True) ? False: True;  break;
	        case  'S':
  	              if (sscanf(OptArg,"%s", SignalPowSpec_FileName) != 1) 
                      {
		              fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	  	           exit(-1);
		      }
                     //  IsotropyCst = True;
		      UpdatePowSpec = False;
		      OptS_PS = True;
              Denoising = True;
  	              break;
            case  's':
                if (sscanf(OptArg,"%s", ReaPowSpec_FileName) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
                    exit(-1);
                }
                Use_Realization_Powspec = True;
                Use_WP_Constraint = False;
                break;
	        case 'I': IsotropyCst = (IsotropyCst == True) ? False: True;  
                      if (IsotropyCst == True) ProjectConstraint = True; 
                      break;
	        case 'C': Use_ZeroMeanConstraint = (Use_ZeroMeanConstraint == True) ? False: True;  break;
         	case 'W': Use_WP_Constraint = (Use_WP_Constraint == True) ? False: True;  
        	           break;
        	case 'H': HardThreshold = (HardThreshold == True) ? False: True;
		          SoftThreshold = (HardThreshold == True) ? False: True;
        	          break;
         	case 'G':  GaussianityConstraint = (GaussianityConstraint == True) ? False: True;  
                           if (GaussianityConstraint == True) Use_WP_Constraint = True;
        	           break;
         	case 'P':  if (sscanf(OptArg,"%f",&GaussCst_ProjectCst_Lambda) != 1) 
                       {
		                  fprintf(OUTMAN, "Error: relaxation parameter  value: %s\n", OptArg);
		                exit(-1);
	                 	}
                       	   GaussianityConstraint = True; 
			            Use_WP_Constraint = True;
        	           break;
         	case 'a': All_WP_Band = (All_WP_Band == True) ? False: True;  
        	           break;
                case 'e':
 	     	      if (sscanf(OptArg,"%f",&Eps) != 1) 
                  {
		            fprintf(OUTMAN, "Error: relaxation parameter  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
                 case 'p': Pos=True; Use_ZeroMeanConstraint=False;break;
		 case 'O':
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
                   OptI = True;
                   if (Max_Iter <= 0)   Max_Iter =  DEF_IWIENER_NBR_ITER;
 	        	break;
	        case 'g':
 	     	      if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		            exit(-1);
	            	}
 	            	Denoising=True;
 		          break;
   	    case 'v': Verbose = True; break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 

    if (LogNormal == True) 
    {
       Use_ZeroMeanConstraint=False;
       Use_WP_Constraint=True;
    }
     if ((OptI == False) && (IsotropyCst == True)) Max_Iter =  DEF_ITER_GAUSSIANITY;

      
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
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# InpMap  File Name Out = " << Name_Imag_Out << endl;   
        if (EstimPowSpec ==True)  cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if ( UseMask ==  True) cout << "# File Name Mask = " <<   Name_Mask << endl;  
        cout << "# Number of iterations = " <<  Max_Iter << endl;   
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
        if (IsotropyCst == True)  cout << "# Use Isotropic constraint instead of l1 Alm minimization" << endl;
        if (Use_WP_Constraint == True) cout << "# Use WP constraint " << endl;
        if (GaussianityConstraint == True) cout << "# Use Gaussianity constraint " << endl;
        if (ProjectConstraint == True) cout << "# Use Projection constraint " << endl;
        if (HardThreshold == True) cout << "# Use HardThresholding " << endl;
        if (SoftThreshold == True) cout << "# Use SoftThresholding " << endl;
        //if (ThresholdResi == True) cout << "# Use residal thresholding " << endl;
        //else cout << "# Don't use residal thresholding " << endl;
        if (InpaintMethod == INP_L1_ALM_ANALYSIS_CST_VAR) cout << "# Use Analysis inpainting " << endl;
        else cout << "# Use Synthesis inpainting " << endl;
        cout << "Inpainting method: " << StringInpaintingMethod(InpaintMethod) << endl;
        if ((UseRMSMap == True) || (SigmaNoise > 0)) cout << "# Use denoising mode " << endl;
        if (UseRMSMap ==  True) cout << "# File Name Noise RMS = " <<   Name_RMSMap << endl;  
        else if (SigmaNoise > 0)  cout << "Sigma Noise = " <<   SigmaNoise << endl; 
        if (OptN_PS > 0)  cout << "# Noise power spectrum file name = " <<  NoisePowSpec_FileName << endl;
        if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  SignalPowSpec_FileName << endl;
        if (Acceleration == True)  cout << "# Acceleration: Level =  " << AccelationLevel << endl;
	if (Denoising == True)  cout << "# Denoising  " << endl; 
        
#ifdef MRSDEBUG                      
        cout << "DEBUG MODE " << endl;
#endif
    }
    
    
   CInpainting  INP;
   INP.Map.read(Name_Imag_In);
   // if (Verbose == True)  INP.Map.info((char*) "Input MAP");
   // double SigmaMap = INP.Map.sigma();
   // double CoefNorm = 500. / SigmaMap;
   // for (int p=0; p < INP.Map.Npix(); p++)  INP.Map[p] *= CoefNorm;
   
   INP.MaskMap.read(Name_Mask);
   // if (Verbose == True)  INP.MaskMap.info((char*) "Input Mask");

   for (int p=0; p < INP.Map.Npix(); p++)  
   {
      if (INP.MaskMap[p] != 1) 
      {
         INP.Map[p] = 0;
         INP.MaskMap[p] = 0;
      }
   }
 
   int Nside = INP.Map.Nside();
   int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
    
   mrs_alloc_powspec(INP.PowSpecData, Lmax_Tot);
   if (OptS_PS == True)
   {
        mrs_read_powspec(SignalPowSpec_FileName, INP.PowSpecSignal);
        if (Verbose == True) cout << " Powspec lmax = " << INP.PowSpecSignal.Lmax() << endl;
        INP.InputPowSpec = true;
        if (Lmax_Tot > INP.PowSpecSignal.Lmax())
        {
           cout << "Error: input signal power spec Lmax = " << INP.PowSpecSignal.Lmax() << endl;
           cout << "       It should be larger or equal to the chosen lmax = " << Lmax_Tot << endl;
           exit(-1);
        }
   }
   if (UseRMSMap == True)
   {
        INP.RMSMap.read(Name_RMSMap);
        if (Verbose == True)  INP.RMSMap.info((char*) "Input RMS Map");
   }
   if (Use_Realization_Powspec == True)
   {
        mrs_read_powspec(ReaPowSpec_FileName, INP.PowSpecRealization);
   }
   
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
   INP.GaussianityConstraint = GaussianityConstraint;
   INP.GaussCst_ProjectCst_Lambda = GaussCst_ProjectCst_Lambda;
   INP.Acceleration=Acceleration;
   INP.AccelationLevel = AccelationLevel;
   INP.Pos=Pos;
   INP.Max_Iter = Max_Iter;
   INP.Denoising = Denoising;
   INP.UseRMSMap = UseRMSMap;
   INP.WienerFilter = WienerOK;
   INP.LogNormal = LogNormal;
   INP.Use_Realization_Powspec = Use_Realization_Powspec;
   INP.NbrMasterIter = NbrMasterIter;
   INP.InpaintMethod = InpaintMethod;
    switch (InpaintMethod)
    {
       case INP_L1_ALM_ANALYSIS_CST_VAR:
            Analysis = True;
            if (AccelationLevel == 0) INP.Use_WP_Constraint = False;
            if (OldAnalysis == True) INP.old_analyse_alm_inpainting();
            else  INP.analyse_alm_inpainting();
            break;
       case INP_L1_ALM_SYNTHESIS:
              INP.synthesis_alm_inpainting();
            break;
      case INP_L1_ALM_SYNT_MASTER:
           if (OptMat == True) 
            {
                 fits_read_dblarr(MaskMatrix_FileName, INP.MatMask);
                 if (Verbose == True)
                 {
                    INP.MatMask.info("MatMask");
                 }
            }
/*
            if (OptRevMat == True) 
            {
                fits_read_dblarr(RevMaskMatrix_FileName, INP.RevMatMask);
                if (Verbose == True)
                {
                    INP.RevMatMask.info("RevMatMask");
                }
            }            
            INP.synthesis_alm_inpainting_with_master_decconv(); 
*/
  
            INP.synthesis_alm_inpainting();
            break;
         case INP_ALM_ISOTROP:
             INP.HardThreshold = INP.SoftThreshold = INP.Use_WP_Constraint = False;
             INP.IsotropyCst = True;
             INP.Denoising = False;
             if (OptMat == True) 
             {
                fits_read_dblarr(MaskMatrix_FileName, INP.MatMask);
                if (Verbose == True)
                {
                    INP.MatMask.info("MatMask");
                }
             }
                       
             // INP.WienerFilter = True;
            /* if (INP.InputPowSpec == False)
             {
                    cout << "Error: this method requires the input signal power spectrum.  " <<   endl;
                    exit(-1);
             } */
             INP.synthesis_alm_inpainting();
             break;
        case INP_L2_ALM:
            INP.L2Constraint=True;
            INP.HardThreshold =INP.SoftThreshold = INP.Use_WP_Constraint = False;
            INP.IsotropyCst = False;
            INP.Denoising = False;
            INP.WienerFilter = False;
            INP.Use_ZeroMeanConstraint=False;
         //   if (INP.InputPowSpec == False)
         //   {
         //       cout << "Error: this method requires the input signal power spectrum.  " <<   endl;
         //       exit(-1);
         //   }
            INP.synthesis_alm_inpainting();
            // INP.analyse_alm_inpainting();
            break;            
    }      
   
   INP.Result.write(Name_Imag_Out);

   INP.ALM.Norm = False;
   if ((Analysis == True) && ((AlmWrite == True) || (EstimPowSpec == True))) 
   {
    // cout << "ALM ALM " << endl;
    INP.ALM.alm_trans(INP.Result);
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

