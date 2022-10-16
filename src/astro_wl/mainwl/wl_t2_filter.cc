/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Sandrine PIRES
**
**    Date: 05/01/06
**    
**    File:  wl_t2_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION: Image filtering using the multi-scale entropy
**     In Variable alpha (FDR method) : reajusting the alpha/l 
**		(thanks to ROC curves) 
**       Options :
**          -n# number of scales
**          -F#  FirstScale used for detection (default is 0)
**          -g# sigma of image (default is automatically calculated)
**          -k  suppress isolated pixels
**          -I# number of loops in reconstruction iterative process 
**           (default is 10)
**          -s# Thresholding at k * SigmaNoise (default is k = 3.) 
**					or alpha in FDR (default is 0.05)
**          -f# Array (fits size n-1) of a different fix threshold per band in input
**          -C# Use a different fdr-threshold per band 
**          (default is cte alpha value = 0.05)
**          -c# Use a different alpha per band 
**           (default is c = variation of alpha = 2.)
**          -P Apply the positivty constraint
**          -e# kind of edges process (default is 0 : continous)
**          -r rms map is calculated
**          -R RMSMAP in input with -m5 (non uniform noise) or -m9 (no white noise)
**					-m# noise model
**
**  Last modification (12/07/2006) : add -f 
**		    
******************************************************************************/
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "CErf.h"
#include "CMem.h"
// #include "MR_Edge.h"
#include "MW_Filter.h"

/*char Name_W_Out[256];  output file name */
char Name_Imag_In[256]; /* input file image */
char Name_MR_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float alpha = 0.05; /* precision of fdr-thresholding - default value*/
float A_var = 2.0; /* variation of alpha for fdr-thresholding - default value*/
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
Bool UseNSigma = False;
Bool UseNSigmaFdr = False;
Bool UseNSigmaFdrl = False;
Bool UseNSigmaVar = False;
Bool UseRMSMap = False;
Bool SupIsol=False;             /* suppress isolated pixel in the support */
char Name_RMSMap[256]; 
char Name_ThreshArr[256]; 
float EpsilonPoisson = DEFAULT_EPSILON;
float RegulVal =1.;
float CvgParam= DEF_MEM_FILT_CONGER;
int NbrIter = DEF_MEM_FILT_MAX_ITER;
int niter = 10;
int bords = 0; /* continue edges*/
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);
#define PREC 0.000001
Bool  DataSNR = True;
Bool  DataEdge = False;
Bool  RegulPerScale = False;
Bool Verbose = False;
Bool Positiv = False;
Bool DetectOnlyPos = False;
int TypeOpt=3; // one alpha value per band;
int FirstScale = DEF_FIRST_DETECT_SCALE;
Bool KillLastScale = False;
type_threshold TypeThreshold = DEF_THRESHOLD;

/****************************************************************************/
/****************************************************************************/

void my_recons(MultiResol & MR_Data, MRNoiseModel & ModelData, int & FirstScale, Ifloat & Imarec)

/* Find a reversible reconstruction */

{		
   int Nbr_Band =  MR_Data.nbr_band()-1;
   int Nbr_Plan =  MR_Data.nbr_band();
   int i,j,b;
   int k;
   int Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   Ifloat I(Nl,Nc,"I");
   Ifloat Res(Nl,Nc,"Residu");
   MultiResol MR_Data_i(I.nl(),I.nc(),Nbr_Plan,Transform, "MRNoiseModel");
   if(bords == 0) MR_Data_i.Border=I_CONT;
   if(bords == 1) MR_Data_i.Border=I_MIRROR;
   if(bords == 2) MR_Data_i.Border=I_PERIOD;
   if(bords == 3) MR_Data_i.Border=I_ZERO;
   MultiResol MR_Data_res(I.nl(),I.nc(),Nbr_Plan,Transform, "MRNoiseModel");
   if(bords == 0) MR_Data_res.Border=I_CONT;
   if(bords == 1) MR_Data_res.Border=I_MIRROR;
   if(bords == 2) MR_Data_res.Border=I_PERIOD;
   if(bords == 3) MR_Data_res.Border=I_ZERO;
   //MR_Data.recons(I);
   for (b = 0; b < FirstScale; b++)
   {	
     for (i = 0; i < Nl; i++)
     {
       for (j = 0; j < Nc; j++)
       {
         MR_Data(b,i,j) = 0;
       }
     }
   }
   MR_Data.rec_adjoint(I);
   for (k = 0; k <  niter; k++)
   {
   MR_Data_i.transform(I);
   cout << "niter=" << k + 1  << endl;
   for (b = FirstScale; b < Nbr_Band; b++)
   {	
     for (i = 0; i < Nl; i++)
     {
       for (j = 0; j < Nc; j++)
       {
         if(ModelData(b,i,j) == True)
           MR_Data_res(b,i,j) = MR_Data(b,i,j) - MR_Data_i(b,i,j);
         else
           MR_Data_res(b,i,j) = 0;
         }
       }
     }
     //MR_Data_res.band(Nbr_Band).init();
     for (i = 0; i < Nl; i++)
     {
       for (j = 0; j < Nc; j++)
       {
         MR_Data_res(Nbr_Band,i,j) = MR_Data(Nbr_Band,i,j) - MR_Data_i(Nbr_Band,i,j);
       }
     }
     MR_Data_res.rec_adjoint(Res);
     for (i = 0; i < Nl; i++)
     {
       for (j = 0; j < Nc; j++)
       {
         if(Res(i,j)>0) I(i,j)=I(i,j)+Res(i,j);	
       }
     }
   }
   Imarec=I;
}   

/****************************************************************************/
/****************************************************************************/
static void usage(char *argv[])
{    
    int i;
    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "              default is automatically estimated\n");
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform\n");
    fprintf(OUTMAN, "              default is %d\n", Nbr_Plan);
    manline();
    
		fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-1; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
    manline();

    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "              Thresholding at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "              default is %f\n", N_Sigma);
    fprintf(OUTMAN, "              in FDR method is alpha value \n");
    fprintf(OUTMAN, "              default is %f\n", 0.05);
    manline();
     
		fprintf(OUTMAN, "         [-f ThreshArr]\n");
    fprintf(OUTMAN, "              Array of number_of_scale -1");
    fprintf(OUTMAN, "							 Thresholding at ThresArr[i] * SigmaNoise\n");
    manline();
 
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "               Apply the positivity constraint.\n");
    fprintf(OUTMAN, "              default is no.\n"); 
 		
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "              Remove the isolated pixels.\n");
    fprintf(OUTMAN, "              default is no.\n"); 
    manline();
    
    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "              Remove the smoothed plane.\n");
    fprintf(OUTMAN, "              default is no.\n"); 
    manline();
    
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "              Use a Thresholding at nsigma * SigmaNoise C=1 \n");
    fprintf(OUTMAN, "              Use a different fdr-threshold per band if C=2 \n");
    fprintf(OUTMAN, "              default is no.\n"); 
     
    fprintf(OUTMAN, "         [-c]\n");
    fprintf(OUTMAN, "              Use a different alpha per band .\n");
    fprintf(OUTMAN, "              default is 2.0 is the variation\n");
    fprintf(OUTMAN, "              default is no.\n"); 
     
    fprintf(OUTMAN, "         [-I niter]\n");
    fprintf(OUTMAN, "              Maximum number of iterations.\n");
    fprintf(OUTMAN, "              default is %d.\n", niter);  
    manline(); 

    fprintf(OUTMAN, "         [-e bords]\n");
    fprintf(OUTMAN, "              number of edges management.\n");
    fprintf(OUTMAN, "              default is %d.\n", bords);  
    manline(); 
      
    fprintf(OUTMAN, "         [-R RMSMap]\n");
    fprintf(OUTMAN, "              RMS Map.\n"); 
    fprintf(OUTMAN, "              If this Option is set, the noise model is automatically fixed to:.\n");
    fprintf(OUTMAN, "              Non-stationary additive noise. \n");
    manline();  
      
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "              rms map is calculated.\n");
    manline(); 
      
    first_detect_scale_usage();
    manline();
    kill_isol_pix_usage();
    manline();
    vm_usage();
    manline();   
         
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n"); 
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
  {
  int c;
   Bool TransfOpt=False;
   Bool NscaleOpt=False;
   #ifdef LARGE_BUFF
   int VMSSize=-1;
   Bool OptZ = False;
   char VMSName[1024] = "";
   #endif   
   /* get options */
   while ((c = GetOpt(argc,argv,"f:F:km:KI:e:g:n:s:vC:c:PR:rzZ:")) != -1) 
     {
     switch (c) 
       {
       case 'k': SupIsol = True; break;		 
       case 'K': KillLastScale = True; break;
       case 'v': Verbose = True; break;
       case 'P': Positiv = True; break;
       case 'g':
       /* -g <sigma_noise> */
       if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
         {
         fprintf(stderr, "Error: bad sigma noise: %s\n", OptArg);
         exit(-1);
         }
       Stat_Noise = NOISE_GAUSSIAN;
       break;   
       //Noise_Ima = 1.;
       case 'n':
       /* -n <Nbr_Plan> */
       if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
    	   {
		     fprintf(stderr, "Error: bad number of scales: %s\n", OptArg);
         exit(-1);
         }
       if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
         {
         fprintf(stderr, "Error: bad number of scales: %s\n", OptArg);
         fprintf(stderr, "1 < Nbr Scales <= %d\n", MAX_SCALE);
         exit(-1);
         }
       NscaleOpt = True;
       break; 
       case 'm':
       if (sscanf(OptArg,"%d",&c ) != 1) 
         {
         fprintf(stderr, "Error: bad type of noise: %s\n", OptArg);
         exit(-1);
         }
       if ((c > 0) && (c < NBR_NOISE)) 
         Stat_Noise = (type_noise) (c-1);
       else  
         {
         fprintf(stderr, "Error: bad type of noise: %s\n", OptArg);
         exit(-1);
         }
       break;
       case 's':
       /* -s <nsigma> */
       if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
         {
         fprintf(stderr, "Error: bad N_Sigma: %s\n", OptArg);
         exit(-1);
         }
       UseNSigma = True;
       if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
       break;
       case 'I':
		   /* -n <Nbr_Plan> */
		   if (sscanf(OptArg,"%d",&niter) != 1) 
    	   {
		     fprintf(stderr, "Error: bad number of iterations: %s\n", OptArg);
		     exit(-1);
			   }
		   break;
      case 'e':
		   /* -n <Nbr_Plan> */
		  if (sscanf(OptArg,"%d",&bords) != 1) 
    	  {
        fprintf(stderr, "Error: bad number of edges management: %s\n", OptArg);
        exit(-1);
        }
      break;
      case 'C':
      /* -f <type> type of filtering */
      if (sscanf(OptArg,"%d",&c ) != 1) 
        {
        fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
        exit(-1);
        }
      if ((c > 0) && (c <= NBR_THRESHOLD)) TypeThreshold = (type_threshold) (c-1);
      else  
        {
        fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
        exit(-1);
        }
      if (TypeThreshold == 0) 
        {
        UseNSigma = True;
        }    
      if (TypeThreshold == T_FDR)
        {
        //N_Sigma = 0.05;
        UseNSigma = True;
        UseNSigmaFdr = True;
        }         
      break;
      case 'c': 
      if (sscanf(OptArg,"%f",&A_var) != 1) 
        {
        fprintf(stderr, "Error: bad threshold: %s\n", OptArg);
        exit(-1);
        }
      UseNSigmaFdrl = True;
      break;	
      case 'R':
      /* -w < support file name> */
      if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
        {
        fprintf(stderr, "Error: bad file name: %s\n", OptArg);
        exit(-1);
        }
      //Stat_Noise = NOISE_NON_UNI_ADD;
      UseRMSMap = True;
      break;
      case 'f':
      /* -w < support file name> */
      if (sscanf(OptArg,"%s", Name_ThreshArr) != 1) 
        {
        fprintf(stderr, "Error: bad file name: %s\n", OptArg);
        exit(-1);
        }
      UseNSigmaVar = True;        
			UseNSigmaFdr = False;
			UseNSigmaFdrl = False;
			TypeThreshold = (type_threshold) 0;
      break;
      case 'r':
      /* -w < support file name> */
      Stat_Noise = NOISE_NON_UNI_ADD;
      UseRMSMap = False;
      break;
      case 'F':
      if (sscanf(OptArg,"%d", &FirstScale) != 1) 
        {
        fprintf(OUTMAN, "Error: bad first detection scale number: %s\n", OptArg);
        exit(-1);
        }
      FirstScale --;
      if (FirstScale < 0)
        {
        fprintf(OUTMAN, "Error: bad FirstScale parameter: %s \n", OptArg);
        fprintf(OUTMAN, "           FirstScale > 0\n");
        exit(-1);
        }
      break;
      #ifdef LARGE_BUFF
      case 'z':
      if (OptZ == True)
        {
        fprintf(OUTMAN, "Error: Z option already set...\n");
        exit(-1);
        }
      OptZ = True;
      break;
      case 'Z':
      if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
        {
        fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
        exit(-1);
        }
      if (OptZ == True)
        {
        fprintf(OUTMAN, "Error: z option already set...\n");
        exit(-1);
        }
      OptZ = True;
      break;
      #endif
      case '?': usage(argv); break;
      default: usage(argv); break;
 		  }
	  }
   /* get optional input file names from trailing 
          parameters and open files */
   if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
   else usage(argv);

   if (OptInd < argc) strcpy(Name_MR_Out, argv[OptInd++]);
   else usage(argv);
   /* make sure there are not too many parameters */
   if (OptInd < argc)
     {
     fprintf(stderr, "Too many parameters: %s ...\n", argv[OptInd]);
     exit(-1);
     }
   // Test inconsistencies in the option call 
   if (Stat_Noise == NOISE_EVENT_POISSON)
     {
     if ((TransfOpt == True) && (Transform != TO_PAVE_BSPLINE))
       {
       cerr << "WARNING: with this noise model, only the BSPLINE A TROUS can be used ... " << endl;
       cerr << "        Type transform is set to: BSPLINE A TROUS ALGORITHM " << endl;
       }
     Transform = TO_PAVE_BSPLINE;
     if (NscaleOpt != True) Nbr_Plan = DEF_N_SCALE;
    }
    if (UseRMSMap == True)
        {
           if ((Stat_Noise != NOISE_NON_UNI_ADD) && (Stat_Noise !=  NOISE_CORREL))
           {
              cerr << "Error: this noise model is not correct when RMS map option is set." << endl;
              cerr << "       Valid models are: " << endl;
              cerr << "        " << StringNoise(NOISE_NON_UNI_ADD) << endl;
              cerr << "        " << StringNoise(NOISE_CORREL) << endl;
              exit(-1);
           }
        } 
    if ((isotrop(Transform) == False) && ((Stat_Noise == NOISE_NON_UNI_ADD) ||
          (Stat_Noise  == NOISE_NON_UNI_MULT)))
      {
      cerr << endl << endl;
      cerr << "  Error: with this transform, non stationary noise models are not valid : " << StringFilter(FILTER_THRESHOLD) << endl;
      exit(-1);
      }
      #ifdef LARGE_BUFF
      if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
      #endif
}

/*********************************************************************/
int main(int argc, char *argv[])
{
    int s,i,j,k;
    Ifloat Imag;    
		Ifloat RmsMap;
    /* support image creation */
    fitsstruct Header;
    char Cmd[256];
    extern softinfo Soft;
    Soft.mr2();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR2);
     /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);       
    if (Verbose == True)
      {  
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Name in = " << Name_Imag_In << endl;
      cout << "File Name Out = " << Name_MR_Out << endl;
      cout << "Transform = " << StringTransform(Transform) << endl;
      cout << "Number of scales = " << Nbr_Plan << endl;
      cout << "Noise Type = " << StringNoise(Stat_Noise) << endl;
      if (Stat_Noise == NOISE_GAUSS_POISSON)
        {
        cout << "Type of Noise = POISSON" << endl;
        cout << "  Gain = " << PasCodeur << endl;
        cout << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
        cout << "  Read-out Mean = " << MeanGauss << endl;
        }
      if (Noise_Ima > 0) cout << "Sigma Noise = " << Noise_Ima << endl;
      cout << "N_Sigma = " << N_Sigma << endl;
    	}
    // Read the input image
    io_read_ima_float(Name_Imag_In, Imag, &Header);
    Header.origin = Cmd;
    // Noise modeling class intialization
    MRNoiseModel ModelData(Stat_Noise, Imag.nl(), Imag.nc(), Nbr_Plan, Transform);
		
    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima == 0) Noise_Ima = Imag.sigma();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
    	{
    	for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
    	}
    for (s=0; s < NbrBand; s++) ModelData.TabEps[s] = EpsilonPoisson;
    if (SupIsol == True) ModelData.SupIsol = True;
    //ModelData.FirstDectectScale = FirstScale;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss; 
		
		if (UseRMSMap == True)
      {
			ModelData.UseRmsMap = True;
      io_read_ima_float(Name_RMSMap, ModelData.RmsMap); 
      }
		
    MultiResol MR_Data(Imag.nl(), Imag.nc(),Nbr_Plan,Transform, "MRNoiseModel");
    if(bords == 0) MR_Data.Border=I_CONT;
    if(bords == 1) MR_Data.Border=I_MIRROR;
    if(bords == 2) MR_Data.Border=I_PERIOD;
    if(bords == 3) MR_Data.Border=I_ZERO;
	
    cout << "bords = " << MR_Data.Border << endl;  
    //ModelData.model(Imag, MR_Data);  
    ModelData.TypeThreshold=TypeThreshold;
    ModelData.FirstDectectScale = FirstScale;
    if (UseNSigmaFdrl  == True)
      {
      double puissance;
      double avar=double(A_var);
      float alpha_b;
      float alpha = N_Sigma;
      float alphaC;	 
			//double a;			
      for (i=0; i < NbrBand; i++) 
        {
        puissance = double((i+1)-(FirstScale+1));
        alpha_b = alpha * pow(avar,puissance);
        if (alpha_b > 0.5)  alpha_b = 0.5;
        alphaC = ABS(xerfc(0.5+(1-alpha_b)/2.));
        ModelData.NSigma[i]=alphaC;
        cout << "alpha_b = " << alpha_b << endl;
				//a = MR_Data.band_norm(i);
				// cout << "a = " << a << endl;
        }
       }
     if ((UseNSigma  == True)&&(UseNSigmaFdrl == False))
       {
       float alpha = N_Sigma;
       float alphaC = N_Sigma;
       if (UseNSigmaFdr  == True)
         {
         cout << "alpha = " << alpha << endl;  
         alphaC = ABS(xerfc(0.5+(1-alpha)/2.));
         }
       for (i=0; i < NbrBand; i++) 
			   {
			   ModelData.NSigma[i]= alphaC;
			   }
       }
     if (UseNSigmaVar == True)
       {  
       fltarray ThreshArr(Nbr_Plan);
       fits_read_fltarr(Name_ThreshArr, ThreshArr);
       for (i=0; i < Nbr_Plan-1; i++) 
         {			
         float Thresh;
         Thresh = ThreshArr(i);
         ModelData.NSigma[i]= Thresh;
         cout << "thresh = " << ModelData.NSigma[i] << endl;
         }
        }
    cout << "FirstScale = " << FirstScale << endl;
    ModelData.model(Imag, MR_Data);      
		//for (i=0; i < NbrBand; i++) cout << "noise = " << MR_Data.band_norm(i)<< endl;
    cout << "Noise_Ima = " << Noise_Ima << endl;  
    // regul. parameter intialization
    int Nb = NbrBand-1;
    fltarray TabAlpha(Nb);
    for (int b=0; b < Nb; b++) TabAlpha(b) = RegulVal;
    // Model for edges.
    MultiResol  *MR_Model = NULL; 	
    // apply the filtering
    mw_filter(MR_Data,  ModelData, TypeOpt, TabAlpha, DataSNR, 
        DataEdge, MR_Model, CvgParam, NbrIter, Positiv, Verbose);
    // write the MR support in a file
    //MR_Data.write("/home/pires/prog_WL/REAL/SPACE/HST_results/cosmos/Wev_temp");
    // Image reconstruction      
    // MR_Data.recons(Imag);
    if (KillLastScale == True)
      {
      int Nl = MR_Data.size_ima_nl();
      int Nc = MR_Data.size_ima_nc();
      for (i = 0; i < Nl; i++)
        {
        for (j = 0; j < Nc; j++)
          {
           MR_Data(Nbr_Plan-1,i,j) = 0;
          }
        }
      //MR_Data.recons(Imag);       
			my_recons(MR_Data, ModelData, FirstScale, Imag);
      }
    else
	    {
       my_recons(MR_Data, ModelData, FirstScale, Imag);
      }
    if (Positiv == True)  threshold(Imag);    
    //if (ModelData.TransImag == True) ModelData.im_invtransform(Imag);    
    // Save the result
    io_write_ima_float(Name_MR_Out,  Imag, &Header);
    exit(0);
  } 
