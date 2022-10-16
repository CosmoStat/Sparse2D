/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  26/01/02
**    
**    File:  cb_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  deconvolution using multiple basis
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MRBase.h"
#include "MR_Deconv.h"
#include "MR_Sigma.h"
#include "FFTN_2D.h"
#include "IM_Prob.h"
#include "MR_Sigma.h"
#include "IM_Deconv.h"
#include "FCur.h"

char Name_Imag_In[256];  // input file image  
char Name_Psf_In[256];   /* PSF */
char Name_Imag_Out[256]; //  output file name  
char Name_Imag_Start[256]; // First guess input solution 
char Name_Imag_ICF[256];   // ICF file name
char Name_Resi[256];            // residual file name 

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale1D = -1;
int NbrScale2D = 5;

float Noise_Ima=0.;
float N_Sigma=3.;
Bool PoissonNoise=False;
Bool AnscombeTrans=False;

Bool PositivSol=True;
int BlockSize=0;
int Nbr_Iter = 500;
 Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
int FirstDetectScale=0;
Bool KillLastScale = False;
Bool ThresholdRecIma = True;

Bool Filtering = True;

Bool TabBase[NBR_MRBASE];
type_mbase TabSelect[NBR_MRBASE];
int NbrBase = 0;

int WT_NbrUndecimatedScale = -1;
Bool WriteAllRec=False;
Bool UseHuberNorm=False;
float FirstSoftThreshold = 1;
float LastSoftThreshold = 1;
Bool RidKillLastScale = True;
Bool TV = True;


Bool WriteResi = False;         // write the residual 
Bool PsfMaxShift = True;        // shift the max of the PSF to the center of the image
Bool GaussConv=False;
Bool UseICF=False;
Bool UseGuess=False;
float RegulParam = 0.;
float Fwhm = 0.;
float Converg = 1.;                      // convergence parameter 
float Epsilon=DEFAULT_EPSILON_DECONV;/* convergence parameter */
type_deconv DecMethod=DEC_GRADIENT; // DEC_CITTERT;  
Bool PositivIma = True;
Bool Optim = False;
float SupRes = 0.;
Bool PosCoeff = False;

const int NBR_OK_METHOD = 3;
static type_deconv TabMethod[NBR_OK_METHOD] = {DEC_CITTERT,DEC_GRADIENT,DEC_LUCY};
const int NBR_OK_TRANS = 4;
static type_mbase TabTrans[NBR_OK_TRANS] = {MB_ATROU,MB_WT,MB_CUR,MB_FCUR};

Bool ExpDecreasingLambda = True;
int FCurNdir=32;
Bool OptBiWT=False;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image in_psf out_result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    
    fprintf(OUTMAN, "         [-d type_of_deconvolution]\n");
    for (int i = 0; i < NBR_OK_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringDeconv((type_deconv) TabMethod[i]));
    fprintf(OUTMAN, "              default is %s\n", StringDeconv((type_deconv)  DecMethod));
    
    fprintf(OUTMAN, "         [-t TransformSelection]\n");
    for (int i = 0; i < NBR_OK_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMRBase ((type_mbase) TabTrans[i]));

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the WT, the a trous, the PMT and the curvelet transform.\n");
    fprintf(OUTMAN, "             default is %d.\n", NbrScale2D);    
    // fprintf(OUTMAN, "             Number of ridgelet in the multi-ridgelet transform.\n");
    nsigma_usage(N_Sigma);
    gauss_usage();

   // fprintf(OUTMAN, "         [-b BlockSize]\n");
   // fprintf(OUTMAN, "             Block Size in the ridgelet transform.\n");
   //  fprintf(OUTMAN, "             default is image size. \n");    
 
    max_iter_usage(Nbr_Iter);   
    converg_param_usage(Epsilon);
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "            If set, strong negative coefficients.\n");
    fprintf(OUTMAN, "            are not detected. Default no.\n");    
 
       
    //poisson_noise_usage();
    //manline(); 

//     fprintf(OUTMAN, "         [-O]\n");
//     fprintf(OUTMAN, "             Supress the block overlapping. Default is no. \n");    
//     manline();

    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Poisson noise. Default is no. \n");    
 
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "             Minimize the Total Variation. \n"); 
    fprintf(OUTMAN, "             Default Regularization parameter for the TV method is %f.\n", RegulParam);
         
    //fprintf(OUTMAN, "         [-T]\n");
    //fprintf(OUTMAN, "             Minimize the Total Variation instead of the L1 norm. \n");    
    //manline();
    
//     fprintf(OUTMAN, "         [-C TolCoef]\n");
//     fprintf(OUTMAN, "              Default is 0.5. \n");    
//     manline();
    
    fprintf(OUTMAN, "         [-F First_Guess]\n");
    fprintf(OUTMAN, "              Input solution file.\n");
    write_residual_usage();        
    vm_usage();
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}
  
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"fl:L:Er:e:G:d:pPt:F:Oi:b:s:g:n:vzZ")) != -1) 
    {
	switch (c)
        {  
	   case 'f': OptBiWT=True; break;
	   case 'r':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Resi) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteResi = True;
 		break;
           case 'e': /* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Epsilon < 0) || (Epsilon > 1.)) 
                           Epsilon = DEFAULT_EPSILON_DECONV;
		break;
	    case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.1;
		 break;
	    case 'E': ExpDecreasingLambda = (ExpDecreasingLambda == True) ? False: True; break;
	    case 'd':
		/* -d <type> type of deconvolution */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
                    
		}        
		if ((c <= 0) || (c > NBR_OK_METHOD))
		{
		   fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	           exit(-1);
                    
		}                             
		DecMethod  = (type_deconv) (TabMethod[c-1]);
                break;
            case 'F':
 		if (sscanf(OptArg,"%s",Name_Imag_Start) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseGuess = True;
 		break; 
 	  // case 'T': TV = True; break;
          case 'L':               
                if (sscanf(OptArg,"%f",&FirstSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0)
                {
                   fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                   exit(-1);
                }
                break;
          case 'l':               
                if (sscanf(OptArg,"%f",&LastSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0)
                {
                   fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
                   exit(-1);
                }
                break;
          case 'P': PosCoeff = (PosCoeff == True) ? False: True; break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_MRBASE))
                {
                   if (TabBase[c-1] == True)
                   {
                       fprintf(OUTMAN, "Error: transform already selected: %s\n", OptArg);
                       exit(-1);
                   }
                   else
                   {
                      TabBase[c-1] = True;
                      TabSelect[NbrBase++] = TabTrans[c-1];
                   }
                }
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
           case 'i':
                if (sscanf(OptArg,"%d",&Nbr_Iter) != 1) 
                {
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'p': AnscombeTrans=True;
                     Noise_Ima = 1.;
                     break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
          case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
                    usage(argv);
                }
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
                break;
           case 'v': Verbose = True;break;
	   case 'n':
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D  > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
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
	 
	if (OptInd < argc) strcpy(Name_Psf_In, argv[OptInd++]);
         else usage(argv);
	 
	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/***************************************************************************/
enum  type_cb_weight {DEC_CB_NO_WEIGHT,  DEC_CB_WEIGHT_PROBA, DEC_CB_WEIGHT_SUPPORT};

class MBDeconv: public MRDeconv, public MRBase
{
   inline float get_resi(float Noise,float Level, float Coef, float Sol)
   {
       float P=0;
       float Resi = Coef - Sol;
       // float CoefTol=0.1;
       // float LevelTol = Noise*CoefTol;
       switch (TypeDataWeight)
       {
	    case DEC_CB_WEIGHT_PROBA:
	     {
	       float Vn = ABS(Coef) / (sqrt(2.)*Noise);
               if (Noise < FLOAT_EPSILON) P = 0;
	       else 
	       {
 	        if (Vn > 3.5) P = 0;
	        else P = (float) erfc (Vn);
  	       }
	       P = 1. - P;
	      }
 	      break;
	   case DEC_CB_WEIGHT_SUPPORT:
	      if (ABS(Coef) > Level) 
	      {
	         P = 1;
		 // if (ABS(Resi) < LevelTol) Resi = 0.;
              }
	      else P = 0;
	      if ((PosCoeff == True) && (Coef < 0)) P = 0;
 	      
	      if ((P == 0) && (ABS(Sol) > Level))
	      {
	         if ((PosCoeff == False) || (Sol > 0))
		 {
		    P = 1;
		    if (Sol > 0) Resi = Level - Sol;
		    else  Resi = -Level - Sol;
		 }  
 	      }
 	      break;
	    case DEC_CB_NO_WEIGHT:
	       P = 1;
	        break;
	}
	Resi *= P;
	return Resi;
   }
   public:
    type_cb_weight TypeDataWeight;
    MBDeconv():MRDeconv(),MRBase()  {TypeDataWeight = DEC_CB_WEIGHT_SUPPORT;
        SupRes=0; cout << "MRBase::N_Sigma = " << MRBase::N_Sigma << endl;}
    void cb_grad(Ifloat &Gradient);
    void cb_iter_deconv();
    void cb_deconv(Ifloat *FirstGuess, Ifloat *ICF);
    float SupRes;
    // void wt_regul(Ifloat &Obj, float Lambda, float Noise_Ima);
    void wt_proj_dec(Ifloat & ImaRec, float WT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
    void cur_proj_dec(Ifloat & ImaCur, float Cur_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);			

    void atrou_proj_dec(Ifloat & ImaAtrou, float AT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
    void fcur_proj_dec(Ifloat & ImaCur, float FCUR_SigmaNoise, float NSigma, float Lamba_Sigma);
};
  
/*****************************************************************************/
// void regulmr(Ifloat & Data, NbrBand, fltarray & TabRegul)

// void regulmr(Ifloat & Obj, Ifloat & Grad, float RegulVal)
// {
//    RegulIma CRG;
//    int i,j,s;
//    int Nl = Obj.nl();
//    int Nc = Obj.nc();
//    int NbrBand = 5;
//    Ifloat Buffer(Nl,Nc,"xx");
//    Ifloat Image(Nl,Nc,"xx");
// 
//    // TO_PAVE_LINEAR, TO_PAVE_BSPLINE, TO_PAVE_HAAR
//    type_transform Transform = TO_PYR_BSPLINE; // TO_PAVE_BSPLINE;
//    MultiResol MR_Data;
//    MR_Data.alloc (Nl, Nc, NbrBand, Transform, "mr");
//    MR_Data.transform(Obj);
//  
//    Nl = MR_Data.size_band_nl(NbrBand-1);  
//    Nc = MR_Data.size_band_nc(NbrBand-1);  
//    Buffer.alloc(Nl, Nc, "Buffer mr_iso_recons");
//    Image.resize (Nl, Nc);
//    Image = MR_Data.band(NbrBand-1);
//    for (s = NbrBand -2; s >= 0; s--)
//    {
//       Nl = MR_Data.size_band_nl(s);  
//       Nc = MR_Data.size_band_nc(s);   
//       Buffer.resize (Nl, Nc);
//       if ((s > 0) || (MR_Data.Set_Transform != TRANSF_SEMIPYR))
//       {
//          im_increase_size_2 (Image, Buffer, MR_Data.Border);
//          Image.resize (Nl, Nc);
//          for (i = 0; i < Nl; i++)
//          for (j = 0; j < Nc; j++)
//                         Image(i,j) = Buffer(i,j) + MR_Data(s,i,j);
//       }
//       else Image += MR_Data.band(s);
//       CRG.obj_regul(Image, Image, RegulVal);
//    }
//    for (i = 0; i < Obj.nl(); i++)
//    for (j = 0; j < Obj.nc(); j++)
//    {
//       Grad(i,j) -= Obj(i,j) - Image(i,j);
//    }
// }

/****************************************************************************/

// void MBDeconv::wt_regul(Ifloat &Obj, float Lambda, float Noise_Ima)
// {
//    int s,i,j; 
//    int LastScale = AT_NbrScale2D;
//    AWT.transform(Obj, AT_Resi, AT_NbrScale2D);
//    for (s=0; s < LastScale-1; s++)
//    {
//       float Norm = (s == AT_NbrScale2D-1) ? AWT.norm_band(s-1): AWT.norm_band(s);
//       float Noise =  AT_SigmaNoise*Norm;
//       float LWT =  Lambda*Noise;
//       for (i=0; i < MRDeconv::Nl; i++)
//       for (j=0; j < MRDeconv::Nc; j++)
//       {
//  	 (AT_Resi[s])(i,j) = soft_threshold((AT_Resi[s])(i,j),LWT) ;
//        }
//    }
//    AWT.recons(AT_Resi, Obj, AT_NbrScale2D);
//  
// //    MBSol->transform(Obj);
// //    for (s=0; s < WT_NbrBand-1; s++)
// //    {
// //       int Nlb = MBSol->size_band_nl(s);
// //       int Ncb = MBSol->size_band_nc(s);
// //       float Noise = TabMBSigma(s);
// //       float Level = TabMBLevel(s);
// //       float LWT =  Lambda*Noise;
// //       for (i=0; i < Nlb; i++)
// //       for (j=0; j < Ncb; j++)
// //       {
// //            float Coef = (*MBSol)(s,i,j);  
// // (*MBSol)(s,i,j) -=  Coef * LWT;
// //            (*MBSol)(s,i,j) = soft_threshold(Coef,LWT);
// //       }
// //    }
// //    MBSol->recons(Obj);
// }

/****************************************************************************/

void MBDeconv::atrou_proj_dec(Ifloat & ImaAtrou, float AT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma)
// Resi = IN: residual image
// AT_Resi = Ifloat[0..AT_NbrScale2D-1] = IN: buffer for WT computation
// AT_Trans = Ifloat[0..AT_NbrScale2D-1] = IN/OUT decomposition on
//                                           a trous algorithm
// ImaAtrou = OUT: reconstruction from AT_Trans
// AT_NbrScale2D = IN: number of scales
// AT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// AT_KillLastScale = IN: the last scale is taken into account
// Bord = IN: type of border management
{
   int s,i,j;
   int LastScale = AT_NbrScale2D;

   AWT.transform(ImaAtrou, AT_Resi, AT_NbrScale2D);
   
 
   for (s=0; s < LastScale-1; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      float Norm = (s == AT_NbrScale2D-1) ? AWT.norm_band(s-1): AWT.norm_band(s);
      float Noise =  AT_SigmaNoise*Norm;
      float Level = (s == AT_NbrScale2D-1) ? 0: NSig*Noise;
      // float LWT = (s == AT_NbrScale2D-1) ?  0: Lamba_Sigma*Noise;

      for (i=0; i < MRDeconv::Nl; i++)
      for (j=0; j < MRDeconv::Nc; j++)
      {
          float Coef = (AT_Trans[s])(i,j); // A trous coef
          // float Resi = Coef - (AT_Resi[s])(i,j);    // residual WT coef
  	 (AT_Resi[s])(i,j) = get_resi(Noise, Level, Coef, (AT_Resi[s])(i,j));
      }
   }
   s = LastScale-1; 
   for (i=0; i < MRDeconv::Nl; i++)
   for (j=0; j < MRDeconv::Nc; j++)
   {
      float Coef = (AT_Trans[s])(i,j); // A trous coef
      float Resi = Coef - (AT_Resi[s])(i,j);    // residual WT coef
      (AT_Resi[s])(i,j) = Resi;
   }
   AWT.recons(AT_Resi, MRDeconv::Resi, AT_NbrScale2D);
}

/***************************************************************************/

void MBDeconv::wt_proj_dec(Ifloat & ImaRec, float WT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma)
// Resi = IN: residual image
// WT_Resi = Ifloat[0..WT_NbrScale2D-1] = IN: buffer for WT computation
// WT_Trans = Ifloat[0..WT_NbrScale2D-1] = IN/OUT decomposition on WT
// ImaRec = OUT: reconstruction from WT_Trans
// WT_NbrScale2D = IN: number of scales
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   int s,i,j;
   Ifloat IPois;
   
   WT->transform(ImaRec, WT_Resi, WT_NbrScale2D, WT_TabDec);
   for (s=0; s < WT_NbrBand-1; s++)
   {
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      float Noise =  WT_SigmaNoise;
      float Level = (s == WT_NbrBand-1) ? 0: NSig*Noise;
      // float LWT = (s == WT_NbrBand-1) ? 0: Lamba_Sigma*Noise;

      for (i=0; i <  (WT_Trans[s]).nl(); i++)
      for (j=0; j <  (WT_Trans[s]).nc(); j++)
      {
         float Coef = (WT_Trans[s])(i,j); // A trous coef
         // float Resi = Coef - (WT_Resi[s])(i,j);    // residual WT coef
  	(WT_Resi[s])(i,j) = get_resi(Noise, Level, Coef, (WT_Resi[s])(i,j));
      }
   }
   s = WT_NbrBand-1; 
   for (i=0; i <  (WT_Trans[s]).nl(); i++)
   for (j=0; j <  (WT_Trans[s]).nc(); j++)
   {
      float Coef = (WT_Trans[s])(i,j); // A trous coef
      float Resi = Coef - (WT_Resi[s])(i,j);    // residual WT coef
      (WT_Resi[s])(i,j) = Resi;
   }
   WT->recons(WT_Resi, MRDeconv::Resi, WT_NbrScale2D, WT_TabDec);
}

/****************************************************************************/

void MBDeconv::fcur_proj_dec(Ifloat & ImaCur, float FCUR_SigmaNoise, float NSigma, float Lamba_Sigma)
{
   int s,b,i,j;
   Ifloat BandSol;
   Ifloat BandData;

   FCurSol->cur_trans(ImaCur);
   for (s=0; s <  FCurSol->nbr_scale()-1; s++)
   for (b=0; b <  FCurSol->nbr_band(s); b++)
   {
      FCurSol->get_band(s,b,BandSol);
      FCurData->get_band(s,b,BandData);
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      float Noise =  FCUR_SigmaNoise * FCurData->norm_band(s,b);
      float Level =  NSig*Noise;
      // float LWT = Lamba_Sigma*Noise;
 
      for (i=0; i <  BandSol.nl(); i++)
      for (j=0; j <  BandSol.nc(); j++)
      {
          // float Resi = BandData(i,j) - BandSol(i,j);    // residual WT coef
          float Coef = BandData(i,j); // A trous coef
          BandSol(i,j) = get_resi(Noise, Level, Coef, BandSol(i,j));     
      }
      FCurSol->put_band(s,b,BandSol);
   }
   s = FCurSol->nbr_scale()-1;
   FCurData->get_band(s,0,BandData);
   FCurSol->get_band(s,0,BandSol);
   for (i=0; i <  BandSol.nl(); i++)
   for (j=0; j <  BandSol.nc(); j++) BandSol(i,j) = BandData(i,j) - BandSol(i,j);
   FCurSol->put_band(s,0,BandSol);
   FCurSol->cur_recons(MRDeconv::Resi);
}

/****************************************************************************/

void MBDeconv::cur_proj_dec(Ifloat & ImaCur, float CurNoise, float NSigma, 
                      float Lamba_Sigma)
// Resi = IN: residual image
// CUR_Resi = Ifloat[0..AT_NbrScale2D-1] = IN: buffer for WT computation
// CUR_Trans = Ifloat[0..AT_NbrScale2D-1] = IN/OUT decomposition on
//                                           a trous algorithm
// ImaCur = OUT: reconstruction from the curvelet transform
// CUR_NbrScale2D = IN: number of scales
// CUR_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// Bord = IN: type of border management
{
   int s2d,s1d,i,j;
   int Nls,Ncs,Depi,Depj;
   float Coef,CoefResi;
 
   CUR.transform(ImaCur, CUR_Resi);
   for (s2d=0; s2d < CUR_NbrScale2D-1; s2d++)
   {
      // int NScale = CUR_Ridgelet[s2d].NbrScale;
      int NScale = CUR.nbr_rid_scale(s2d);
      for (s1d=0; s1d < NScale; s1d++)
      {
         float NSig = (s1d == 0) ? NSigma + 1: NSigma;
         float BS = ((CUR.TabBlockSize)(s2d) <= 0) ? MRBase::Nc: (CUR.TabBlockSize)(s2d);
         // float Norm = (s1d == NScale-1) ? CUR_Ridgelet[s2d].rid_norm(s1d-1) : CUR_Ridgelet[s2d].rid_norm(s1d);
         float Norm = (s1d == NScale-1) ? CUR.norm_band(s2d,s1d-1) : CUR.norm_band(s2d,s1d);
         float Noise = sqrt(BS)*CurNoise*Norm; 
         float Level = NSig*Noise;
         // float LRid = Lamba_Sigma*Noise;
         Nls =  CUR.size_nl(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nl(s1d);
         Ncs =  CUR.size_nc(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nc(s1d);
	 Depi = CUR.ipos(s2d, s1d); // CUR_Ridgelet[s2d].ipos(s1d);
         Depj = CUR.jpos(s2d, s1d); // CUR_Ridgelet[s2d].jpos(s1d);

         for (i=Depi; i < Depi+Nls; i++)
         for (j=Depj; j < Depj+Ncs; j++)
         {
            Coef = CUR_Trans(j,i,s2d);
            CoefResi = Coef - CUR_Resi(j,i,s2d);
    	    CUR_Resi(j,i,s2d) = get_resi(Noise, Level, Coef, CUR_Resi(j,i,s2d));
         }
      }
   }
   s2d = CUR_NbrScale2D-1; 
   for (i=0; i < MRBase::Nl; i++)
   for (j=0; j < MRBase::Nc; j++)
   { 
      float Coef = CUR_Trans(j,i,s2d);  
      float Resi = Coef - CUR_Resi(j,i,s2d);    
      CUR_Resi(j,i,s2d) = Resi;
   }
   
   CUR.recons(CUR_Resi, MRDeconv::Resi);
}

/****************************************************************************/

void MBDeconv::cb_iter_deconv()
{
    Ifloat Gradient(MRDeconv::Nl, MRDeconv::Nc, "Info");
    float Func, OldFunc=1e10, Sigma;
    int i,j,b,s,Iter=0;
    Bool Stop = False;
    float SigmaSignif,Cvgparam = IterCvg;
    float Sigma2N = Noise_Ima*Noise_Ima*MRDeconv::Nl*MRDeconv::Nc;
    if (NormFlux == True) FluxImag = flux(Imag);
    else FluxImag=0.;
    float DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    StepL = DeltaThreshold / (float) Nbr_Iter;
    Lambda = FirstSoftThreshold;
    Ifloat IPois;
    IPois.resize(MRBase::Nl, MRBase::Nc);

    // Initialization
    if (MRBase::Verbose == True) cout << "Initialization ... " << endl;
//     if (MBSol == NULL)
//     {
//         MB_SelectFilter = new FilterAnaSynt(F_MALLAT_7_9);
//         MB_SB1D = new SubBandFilter(*MB_SelectFilter, NORM_L2);
//         MBSol = new MIRROR_2D_WT(*MB_SB1D);
//         MBSol->alloc(MRBase::Nl, MRBase::Nc, MB_NbrScale2D);
//         MB_NbrBand = MBSol->nbr_band();
//     }
//     cout << "Noise manadgement .. " << endl;
//     CImaProb CP;
//     FFTN_2D FFT2D;
//     int Nl1 = 2*MRBase::Nl;
//     int Nc1 = 2*MRBase::Nc;    
//     Ifloat OutPsf(Nl1, Nc1, "OPsf");   
//     Icomplex_f Psf_cf(Nl1, Nc1, "OPsfcf");
//     dec_center_psf (Psf, OutPsf);
//     FFT2D.fftn2d(OutPsf, Psf_cf);
//     if (Noise_Ima < FLOAT_EPSILON)  
//                             Noise_Ima = detect_noise_from_med (Data);  
//     Ifloat NoiseData(MRBase::Nl,MRBase::Nc,"noise data");
//     im_noise_gaussian (NoiseData, Noise_Ima, 100);
//     dec_inverse(NoiseData, Psf_cf, NoiseData, 0.001);  
//     MBSol->transform(NoiseData);
//     TabMBLevel.alloc(MBSol->nbr_band());
//     TabMBSigma.alloc(MBSol->nbr_band());
//     Ifloat IBand;
//     float TabH=0.;
//     float TabV=0.;
//     float TabD=0.;
//     for (b=0; b < MBSol->nbr_band(); b++)
//     {
//        double LMin, LMax;
//        int Nlb = MBSol->size_band_nl(b);
//        int Ncb = MBSol->size_band_nc(b);
//        MBSol->get_band(IBand, b);
//        CP.set(IBand);
//        CP.find_gthreshold(MRBase::N_Sigma, LMin, LMax);
//        TabMBLevel(b) = MAX(ABS(LMin), LMax);
//        TabMBSigma(b) = sigma(IBand);
//        if (MBSol->get_indwt_band(b) == 0) 
//        {
//          if (TabMBSigma(b) > TabH) 
// 	 {
// 	     if ((TabH !=0) && (TabH < TabMBSigma(b)/5.)) TabH = 100.*TabMBSigma(b);
// 	     else TabH = TabMBSigma(b);
//  	 }
//          TabMBSigma(b) = TabH;
//        }
//        if (MBSol->get_indwt_band(b) == 1)  
//        {
//  	 if (TabMBSigma(b) > TabV) 
// 	 {
// 	     if ((TabV !=0) && (TabV < TabMBSigma(b)/5.)) TabV = 100.*TabMBSigma(b);
// 	     else TabV = TabMBSigma(b);
//  	 }
//          TabMBSigma(b) = TabV;
//        }
//        if (MBSol->get_indwt_band(b) == 2)
//        {
//  	 if (TabMBSigma(b) > TabD) 
// 	 {
// 	     if ((TabD !=0) && (TabD < TabMBSigma(b)/5.)) TabD = 100.*TabMBSigma(b);
// 	     else TabD = TabMBSigma(b);
//  	 }
//          else TabMBSigma(b) = TabD;
//        }
//        cout << " b = " << b+1 << " Sigma = " << TabMBSigma(b) << " Level = " << TabMBLevel(b) << endl;
//     }
//     cout << "End Noise manadgement .. " << endl;
          
    
    // Residual estimation
    psf_convol (Obj, Psf_cf, Imag_n);
    Result = Obj;
    // INFO_X(Imag_n, "Imag_n");
    // INFO_X(Obj, "obj");
    if (MRDeconv::Verbose == True) 
      cout << "Start iterating " << endl;
    Sigma = sigma(Data);
    do
    {
       Lambda -= StepL; 
       // Lambda = LastSoftThreshold + DeltaThreshold  *(1.-erf(2.5*Iter/ Nbr_Iter));
       if (Lambda < 0) Lambda = 0.;
       Result = Imag_n;
       // cout << "Iter " << Iter+1 <<  " LambdaTV = " << LambdaTV << endl;
        for (b=0; b < NbrBase; b++)
	// b = Iter % NbrBase;
        {
          if (AnscombeTrans == True) 
          {
             threshold(Result);
             noise_poisson_transform (Result, IPois);
          }
	  else IPois = Result;
    
          //if (MRBase::Verbose == True)
          //  cout << "  transform: " << StringMRBase(TabSelect[b]) << " Lambda = " << Lambda << endl;
	  switch(TabSelect[b])
          {
             case MB_PMT:
                   pmt_proj(IPois, DataSigmaNoise, MRBase::N_Sigma, Lambda); 
                  break;
            case MB_ATROU:
                   atrou_proj_dec(IPois, AT_SigmaNoise, MRBase::N_Sigma, Lambda);  
                  break;
             case MB_RID:
                   rid_proj(RID, RID_Trans, Result, DataSigmaNoise,MRBase::N_Sigma,  
                             RID.KillLastScale, RID_BlockSize, 
                             Lambda, RID_FirstDetectScale);   
                    break;
             case MB_PYRRID:
                   for (s = MRID_NbrRid-1; s >= 0; s--)
                  {
                      rid_proj(MRID[s], MRID_Trans[s], IPois, DataSigmaNoise,MRBase::N_Sigma,  
                                  MRID[s].KillLastScale,  MRID_TabBlockSize[s], 
                                  Lambda, MRID_FirstDetectScale);
                      Result = Imag_n;
                  }
                  break;
             case MB_WT:
                   wt_proj_dec(IPois,  DataSigmaNoise, MRBase::N_Sigma, Lambda);
                   break;
             case MB_COS:
                   cos_proj(IPois, DataSigmaNoise,  MRBase::N_Sigma, Lambda);
                  break;
             case MB_CUR:
                   cur_proj_dec(IPois, DataSigmaNoise, MRBase::N_Sigma, Lambda);
                  break;
	    case MB_FCUR:
                  fcur_proj_dec(IPois, DataSigmaNoise, MRBase::N_Sigma,  Lambda);
                 break;
           default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }
	  // if (TabSelect[b] !=  MB_PYRRID) Result = Imag_n;
       // INFO_X(MRDeconv::Resi, "Resi");
       
          if (AnscombeTrans == True)
          {
             noise_residu_estimation (Result, MRDeconv::Resi, MRDeconv::Resi);
          }
	  Result += MRDeconv::Resi;
       }
       for (i= 0; i < MRDeconv::Nl; i++)
       for (j= 0; j < MRDeconv::Nc; j++)  
          MRDeconv::Resi(i,j) = Result(i,j) - Imag_n(i,j);   

       if (DecMethod == DEC_LUCY)
       {
          for (i= 0; i < MRDeconv::Nl; i++)
          for (j= 0; j < MRDeconv::Nc; j++) 	  
          {
	      if (Imag_n(i,j) > FLOAT_EPSILON)  
	             Gradient(i,j) = MRDeconv::Resi(i,j) / Imag_n(i,j);
              else  Gradient(i,j) = 0.;
	  }
          psf_convol_conj (Gradient, Psf_cf);
          Gradient *= Obj;
       }
       else if (DecMethod == DEC_GRADIENT) 
                   psf_convol_conj (MRDeconv::Resi, Psf_cf, Gradient);
       else Gradient = MRDeconv::Resi;
       
       // if ((TotalVariation == True) && (Lambda > 0))
       //                          // regulmr(Obj, Gradient, Lambda); 
       //		      MRIm.obj_regul(Obj, Gradient, Lambda*Noise_Ima);
       // wt_regul(Obj, Lambda, Noise_Ima);
       
       for (i = 0; i < MRDeconv::Nl; i++) 
       for (j = 0; j < MRDeconv::Nc; j++) 
            Obj(i,j) += Cvgparam * Gradient(i,j);
	float TVParam = Lambda*LambdaTV;
	if ((TotalVariation == True) && (TVParam > 0))
		       MRIm.im_soft_threshold(Obj, Obj, TVParam*Noise_Ima);
	               //  MRIm.im_soft_iwt_threshold(Obj, Obj, TVParam*Noise_Ima);
                       // MRIm.im_soft_non_ortho_threshold(Obj, Obj, TVParam*Noise_Ima);
 
       if (PositivConstraint == True) threshold(Obj);
       if (SupRes > 0)
       {
          for (i= 0; i < MRDeconv::Nl; i++)
          for (j= 0; j < MRDeconv::Nc; j++) 
               if (Obj(i,j) > SupRes) Obj(i,j) = SupRes;
       }  
       if ((NormFlux == True) && (KillLastScale == False))
       {
          float FluxObj = flux(Obj);
	  FluxObj = FluxImag / FluxObj;
	  for (i= 0; i < MRDeconv::Nl; i++)
	  for (j= 0; j < MRDeconv::Nc; j++) Obj(i,j) *= FluxObj;
       }

       if (MRDeconv::Verbose == True) SigmaSignif = sigma(MRDeconv::Resi);	       
       psf_convol (Obj, Psf_cf, Imag_n);
       for (i= 0; i < MRDeconv::Nl; i++)
       for (j= 0; j < MRDeconv::Nc; j++) 
                MRDeconv::Resi(i,j) = Data(i,j) - Imag_n(i,j);
   
       Iter ++;
       
       if ((MRDeconv::Verbose == True) && ((Iter >= 0) && (Iter % 1 == 0)))
       {
          Func = fonctional() / (float) (MRDeconv::Nl*MRDeconv::Nc);
          printf("%d: SigmaResi = %f, SigmaSignif = %f, Func = %f, LambdaTV = %f\n", Iter, Sigma, SigmaSignif, Func, TVParam);
  	  // INFO_X(MRDeconv::Resi, "RESI");
	  // INFO_X(Gradient, "Gradient");
       }
       Sigma = sigma(MRDeconv::Resi);

//        if ((MRDeconv::Verbose == True) && (Iter % 100 == 0)) 
//        {
//           char ch[256];
// 	  sprintf(ch, "xx_sol%d.fits", Iter);
//           io_write_ima_float(ch, Obj);
//        }
       
       if (Iter == MaxIter) Stop = True;
       else if (EpsCvg > FLOAT_EPSILON)
       {       
          if (OldFunc < Sigma) Stop = True;
          if ((ABS(OldFunc - Sigma)) < EpsCvg) Stop = True;
       }
       OldFunc = Sigma;
    } while (Stop == False);
    
    if (MRDeconv::Verbose == True)
    {
       INFO_X(Obj, "Solution: ");
       INFO_X(MRDeconv::Resi, "Resi: ");
    }
}

/****************************************************************************/

void MBDeconv::cb_deconv(Ifloat *FirstGuess, Ifloat *ICF)
{
   Bool IterMethod = True;
   
   switch (DecMethod)
   {
       case DEC_INVERSE:
            IterMethod = False;
            StatNoise = NOISE_GAUSSIAN;
            MRDeconv::init_deconv(FirstGuess,  ICF);
	    dec_inverse (Imag, Psf_cf, Obj);
	    compute_resi();
            break;       
      case DEC_MAP:
             TypeInit = DEC_INIT_FLAT;
             OptimParam = False;
	     StatNoise = NOISE_POISSON;
 	     break;
       case DEC_LUCY:
             TypeInit = DEC_INIT_FLAT;
             StatNoise = NOISE_POISSON;
 	     break;	
      case DEC_CITTERT:
      case DEC_MEM:
      case DEC_MEM_MODEL:
             TypeInit = DEC_INIT_FLAT;
             NormFlux = True;
             OptimParam = False;
	     StatNoise = NOISE_GAUSSIAN;	
 	     break;
      case DEC_GRADIENT:
      case DEC_TIKHONOV:
             TypeInit = DEC_INIT_FLAT;
             StatNoise = NOISE_GAUSSIAN;	
 	     break;
      default:
                cerr << "mr_deconv: Not implemented in this procedure ... ";
                cerr << endl;
                exit (0);
                break;
   }
   if (AnscombeTrans == True)
                noise_poisson_transform (Data, Data);
   init_decomposition();
   if (AnscombeTrans == True)
                noise_inverse_poisson_transform (Data, Data);
		
   if (IterMethod == True)
   {
      MRDeconv::init_deconv(FirstGuess,ICF);
      cb_iter_deconv();
   }
   if (GaussConv == True)  convol_gauss(Obj, Fwhm);
   else if (UseICF == True) psf_convol(Obj, Ima_ICF, Obj);
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    int i,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR4);
    for (i = 0; i < NBR_MRBASE; i++) TabBase[i] = False;

    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl; 
	if (DecMethod == DEC_CITTERT)
	         cout << "Van Cittert iteration " << endl;
	else if (DecMethod == DEC_GRADIENT) 
	         cout << "Fixed Step Gradient " << endl; 
	else cout << "Lucy iteration " << endl; 
	
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        cout << "NbrScale2D = " << NbrScale2D << endl;
        if (FirstDetectScale > 0) cout << "FirstDetectScale = " << FirstDetectScale << endl;
        if (AnscombeTrans == True)
           cout << "Poisson noise " << endl; 
         cout << "NSigma = " << N_Sigma  << endl;  
	 if (RegulParam > 0) cout << "TV regularization, lambda = " << RegulParam << endl;
    }

    MBDeconv MB;  // Multiple base decomposition Class
    Ifloat Guess, Ima_ICF;

    io_read_ima_float(Name_Imag_In, MB.Data, &Header);
    io_read_ima_float(Name_Psf_In, MB.Psf);
    if (UseGuess == True) io_read_ima_float(Name_Imag_Start, Guess);
    if (UseICF == True) io_read_ima_float(Name_Imag_ICF, Ima_ICF);
        
    Header.origin = Cmd;
    if (NbrBase == 0)
    {
       NbrBase = 2;
       TabSelect[0] = MB_ATROU;
       TabSelect[1] = MB_FCUR;
    }
     
    if ((AnscombeTrans == False) && (Noise_Ima < FLOAT_EPSILON))
    {
       Noise_Ima = detect_noise_from_med (MB.Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    if (AnscombeTrans == True) Noise_Ima = 1.;

    // for (i = 0; i < NBR_MRBASE; i++) MB.TabBase[i] = TabBase[i];
    MB.SupRes = SupRes;
    MB.NbrBase = NbrBase;
    for (i = 0; i < NbrBase; i++) MB.TabSelect[i] = TabSelect[i];
    MB.Filtering = Filtering;
    MB.DataSigmaNoise = Noise_Ima;
    MB.Noise_Ima = Noise_Ima;
    MB.MRBase::N_Sigma = N_Sigma;
    MB.MRDeconv::N_Sigma = N_Sigma;
    MB.PoissonNoise = PoissonNoise;
    MB.MRBase::Verbose = Verbose;
    MB.MRDeconv::Verbose = Verbose;
    MB.Nbr_Iter = Nbr_Iter;
    MB.PositivRecIma = PositivSol;
    MB.FirstSoftThreshold = FirstSoftThreshold;
    MB.LastSoftThreshold = LastSoftThreshold;

    MB.RID.RidTrans = RidTrans;
    MB.RID.BlockOverlap = BlockOverlap;
    MB.RID_FirstDetectScale = FirstDetectScale;
    MB.RID.KillLastScale = RidKillLastScale;
    if (BlockSize > 0) MB.RID_BlockSize = BlockSize;

    MB.AT_NbrScale2D = NbrScale2D;
    MB.AT_PositivRecIma = True;
    MB.AT_KillLastScale  = KillLastScale;

    MB.PMT_NbrScale2D = NbrScale2D;
    MB.PMT_PositivRecIma = True;
    MB.PMT_KillLastScale  = KillLastScale;
 
    if (BlockSize > 0) MB.MRID_BlockSize = BlockSize;
    MB.MRID_NbrRid = NbrScale2D;
    MB.MRID_FirstDetectScale = FirstDetectScale;
    MB.MRID_BlockOverlap = BlockOverlap;
    MB.MRID_KillLastScale = RidKillLastScale;

    MB.WT_NbrScale2D = NbrScale2D;
    MB.WT_NbrUndecimatedScale = WT_NbrUndecimatedScale;
    MB.MB_NbrScale2D = NbrScale2D;
    MB.CUR_NbrScale2D = NbrScale2D;
    // MB.CUR_BlockSize = BlockSize;
    MB.CUR_BlockOverlap = BlockOverlap;
    MB.CUR_KillLastScale  = KillLastScale;
    MB.TotalVariation = TV;

    MB.FCUR_NbrDir = FCurNdir;
    MB.FCUR_NbrScale2D = NbrScale2D;
    
       
    MB.Imag = MB.Data;
    MB.PositivConstraint = PositivIma;
    MB.DecMethod = DecMethod;
    MB.PsfMaxShift = PsfMaxShift;
    // MB.Noise_Ima = Noise_Ima;
    // MB.RegulParam = RegulParam;
    MB.MaxIter = Nbr_Iter;
    MB.EpsCvg = Epsilon;
    MB.IterCvg = Converg;
    MB.GaussConv = GaussConv;
    MB.Fwhm = Fwhm;
    MB.OptimParam = Optim;
    Ifloat *Pt_G = NULL;
    if (UseGuess == True) Pt_G = &Guess;
    Ifloat *Pt_ICF = NULL;
    if (UseICF == True) Pt_ICF = &Ima_ICF;
    
    //DECONVOLUTION
    if (Verbose == TRUE) cout << " Start the deconvolution ... " << endl;
    MB.alloc();    
    MB.MRIm.ExpDecreasingLambda = ExpDecreasingLambda;
    // MB.MRIm.NbrScale=2;
    if (OptBiWT == True) MB.MRIm.TypeFilter = F_MALLAT_7_9;
    MB.LambdaTV = RegulParam;

    MB.cb_deconv(Pt_G, Pt_ICF);
    
     // MB.decomposition();
     io_write_ima_float(Name_Imag_Out, MB.Obj, &Header);
     MB.Data -= MB.Imag_n;
     if (WriteResi == True) io_write_ima_float(Name_Resi, MB.Data, &Header);
    exit(0);
}
