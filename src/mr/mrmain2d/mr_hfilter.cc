/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  01/02/00
**    
**    File:  mr_hfilter.cc
**
*******************************************************************************
**
**    DESCRIPTION  filter an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**  PB with decimated haar transform + laplacien regularization
**  -x , -d and -i options are not visible
**  -x simple thresholding
**  -d decimating (orthogonal wavelets)
**  -i number of iterations per scale
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Deconv.h"
#include "MR_HaarPoisson.h"
#include "MR_HaarPoissonFilter.h"
#include "NR.h"
#include "FFTN_2D.h"
#include "IM_Regul.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_BGR_In[256]; /* output file name */
char Name_PSF_In[256];

char Name_Write_Sup[256];          /* output support file name */
int Nbr_Plan=-1;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
// type_transform Transform = TO_MALLAT; /* type of transform */
type_transform Transform = TO_UNDECIMATED_MALLAT;
int Max_Iter = DEF_NBR_ITER_REC_PER_SCALE; /* Maximum number of iteration */
int MaxFilterIter =  DEF_NBR_ITER_FILTER;

type_filter Filter = FILTER_THRESHOLD;  /* default filtering method */
float Epsilon = DEFAULT_EPSILON_FILTERING; /* convergence parameter */
 Bool UseNSigma =False;
Bool PositivIma = DEF_POSITIV_CONSTRAINT;
Bool Decimated = False;

Ifloat BGR_Ima;

float Lambda=-1.;

float EpsilonPoisson;
Bool KeepPositivSup = False;
int FirstScale = DEF_FIRST_DETECT_SCALE;
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool EpsOpt=False;
Bool MaxIterOpt=False;
Bool PosOpt=False;
  
Bool Verbose=False;
Bool SoftThreshold = False;
Bool DeconvOpt = False;

Bool NscaleOpt=False;
Bool BGROpt = False;
Bool OptL = False;
Bool OptI = False;
Bool OptT = False;
float RegulParam = 0.1;

#define TEST 0
 

/****************************************************************************/

/****************************************************************************/

// void HFilter::deconv(Ifloat & Data, Ifloat & Psf, MultiResol &MR_Data, Ifloat &Result, float RegulParam)
// {   
//     FFTN_2D FFT2D;
//     int Iter,i,j;
//     Ifloat Resi, Grad;
//     float Mean, Delta;
//     Icomplex_f Psf_cf(Nlima,Ncima,"Psf cf");
//     Icomplex_f O_cf(Nlima,Ncima,"Psf cf");
//     FFT2D.fftn2d(Psf, Psf_cf);
//     float FluxPsf = 1. / Psf_cf(Nlima/2,Ncima/2).real();
//     float Cvgparam = 1.;
//     if (RegulParam > 0.5) Cvgparam = 1. / (2.*RegulParam );
//     
//     RegulIma CRG;
//     CRG.GradOperType = OPER_MARKOV_8; // {OPER_LAPLACIAN,OPER_DIR_LAPLACIAN,OPER_MARKOV_4,OPER_MARKOV_8};
//     Bool WaveResi = False;
//     
//     Resi.alloc(Nlima,Ncima,"resi");
//     Grad.alloc(Nlima,Ncima,"resi");
//     Mean = average(Data);
//     BGR_Ima.init(Mean);
//     if (WaveResi == True) get_mr_ima(BGR_Ima, NbrPlan);
//     // filter_noiter(MR_Data, Result);
//     // BGR_Ima = Result;
//     Result = BGR_Ima;
//     for (i=0; i < Nlima; i++)
//     for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
// 		   	
//     for (Iter=0; Iter < MaxFilterIter; Iter++)
//     {
//  	// Residual Filtering 
// 	if (WaveResi == True)
// 	{
// 	   MR_Data.transform(Resi);
//  	   threshold_all_scale(MR_Data);
// 	   MR_Data.recons(Resi);
// 	}
// 	
//  	// New solution
//         for (i = 0; i < Nlima; i++)
//         for (j = 0; j < Ncima; j++) 
//         {
//            if (BGR_Ima(i,j) > FLOAT_EPSILON) Resi(i,j) = Resi(i,j) / BGR_Ima(i,j);
//            else  Resi(i,j) = 0.;
//         }
// 	FFT2D.convolve_conj(Resi, Psf_cf);
//         Resi *= Result;
// 	
// 	if (RegulParam > 0.) 
// 	{
// 	   CRG.ima_regul(Result,Grad,RegulParam);
// 	   Resi -= Grad;
// 	}
// 	
// 	for (i = 0; i < Nlima; i++)
//         for (j = 0; j < Ncima; j++) 
//         {
//            Result(i,j) += Cvgparam*Resi(i,j);
// 	}
// 	threshold(Result);	
// 	// New image
// 	FFT2D.convolve(Result, Psf_cf, BGR_Ima);
//         threshold(BGR_Ima);
// 	for (i=0; i < Nlima; i++)
//         for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
// 	 
// 	if (WaveResi == True) get_mr_ima(BGR_Ima, NbrPlan);
//         if (Verbose == True) 
// 	{
//   	   Delta = sigma(Resi);
//            cout << "Filter. Iter " << Iter+1 << " Delta = " << Delta << endl << endl;
//         }
//     }   
// }

/****************************************************************************/

void HFilter::deconv (Ifloat & Data, Ifloat & Psf, MultiResol &MR_Data, Ifloat &Result, float RegulParam)
{   
    FFTN_2D FFT2D;
    int Iter,i,j,b,s;
    Ifloat Resi, Grad;
    float Mean, Delta;
    Icomplex_f Psf_cf(Nlima,Ncima,"Psf cf");
    Icomplex_f O_cf(Nlima,Ncima,"Psf cf");
    FFT2D.fftn2d(Psf, Psf_cf);
    float FluxPsf = 1. / Psf_cf(Nlima/2,Ncima/2).real();
    float Cvgparam = 1.;
    if (RegulParam > 0.5) Cvgparam = 1. / (2.*RegulParam );
    float CoefTol = 0.5 / (float) N_Sigma;
     double LambdaCoef,Level;
    float LastSoftThreshold = 0.;
    float FirstSoftThreshold = RegulParam;
    float DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
    float StepL = DeltaThreshold / (float) MaxFilterIter;
    Bool LucyMethod = False;
    
    RegulIma CRG;
    CRG.GradOperType = OPER_MARKOV_8; // {OPER_LAPLACIAN,OPER_DIR_LAPLACIAN,OPER_MARKOV_4,OPER_MARKOV_8};
    Bool WaveResi = True;
    
    Resi.alloc(Nlima,Ncima,"resi");
    Grad.alloc(Nlima,Ncima,"resi");
    Mean = average(Data);
    BGR_Ima.init(Mean);
    if (WaveResi == True) get_mr_ima(BGR_Ima, NbrPlan);
    // filter_noiter(MR_Data, Result);
    // BGR_Ima = Result;
    Result = BGR_Ima;
    for (i=0; i < Nlima; i++)
    for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
    MultiResol MR_Resi,MR_Ima;
     
    MR_Resi.alloc (Nlima, Ncima, MR_Data.nbr_scale(), 
                   MR_Data.Type_Transform, MR_Data.filter_bank(), NORM_L2);
    
    MR_Ima.alloc (Nlima, Ncima, MR_Data.nbr_scale(), 
                   MR_Data.Type_Transform, MR_Data.filter_bank(), NORM_L2);
		   	
    for (Iter=0; Iter < MaxFilterIter; Iter++)
    {
        // RegulParam = LastSoftThreshold + DeltaThreshold  *(1.-erf(2.5*Iter/ MaxFilterIter));
 	
	// Residual Filtering 
	if (WaveResi == True)
	{
           MR_Ima.transform(Resi);
           for (s = NbrPlan-2; s >= 0; s--)
           {
              NormCoef1 = 1. / (double) POW2(s+2);
              NormCoef2 = (double) POW2(2*(s+1));
             //cout << "LEVEL " << s << endl;
             for (b=3*s; b <= 3*s+2; b++)
             {
	        // cout << "BAND " << b << endl;
                for (i=0; i < MR_Data.size_band_nl(b); i++)
                for (j=0; j < MR_Data.size_band_nc(b); j++)
                {
                   float HaarCoeff = MR_Data(b,i,j);
                   LambdaCoef = (TabImaModel[s])(i,j);
                   Level = get_level(LambdaCoef, s);
		   // Non significant coefficient
		   // if (ABS(MR_Data(b,i,j)) < Level) MR_Ima(b,i,j) = 0;
		   if (MR_Data(b,i,j) < Level) MR_Ima(b,i,j) = 0;
		   else 
		   {
		      if (ABS(MR_Ima(b,i,j)) < Level*CoefTol) MR_Ima(b,i,j) = 0;
		   }
                }
	     }
	  }
 	  MR_Ima.recons(Resi);
 	}	
 	// New solution
	if (LucyMethod == True)
	{
          for (i = 0; i < Nlima; i++)
          for (j = 0; j < Ncima; j++) 
          {
           if (BGR_Ima(i,j) > FLOAT_EPSILON) Resi(i,j) = Resi(i,j) / BGR_Ima(i,j);
           else  Resi(i,j) = 0.;
          }
	  FFT2D.convolve_conj(Resi, Psf_cf);
          Resi *= Result;
	}
	else
	{
	   FFT2D.convolve_conj(Resi, Psf_cf);
	}
	 
	if (RegulParam > 0.) 
	{
	   CRG.ima_regul(Result,Grad,RegulParam);
	   Resi -= Grad;
	}
	if (RegulParam > 1.) Cvgparam = 1. / RegulParam;
	else Cvgparam = 1.;
	for (i = 0; i < Nlima; i++)
        for (j = 0; j < Ncima; j++) 
        {
           Result(i,j) += Cvgparam*Resi(i,j);
	}
	threshold(Result);	
	// New image
	FFT2D.convolve(Result, Psf_cf, BGR_Ima);
        threshold(BGR_Ima);
	for (i=0; i < Nlima; i++)
        for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
	 
	if (WaveResi == True) get_mr_ima(BGR_Ima, NbrPlan);
        if (Verbose == True) 
	{
  	   Delta = sigma(Resi);
           cout << "Filter. Iter " << Iter+1 << " Delta = " << Delta << endl << endl;
        }
    }   
}

/****************************************************************************/

void HFilter::jammal_deconv (Ifloat & Data, Ifloat & Psf, MultiResol &MR_Data, Ifloat &Result, float RegulParam)
{   
    FFTN_2D FFT2D;
    int Iter,i,j,b,s;
    Ifloat Resi, Grad;
    float Mean, Delta;
    Icomplex_f Psf_cf(Nlima,Ncima,"Psf cf");
    Icomplex_f O_cf(Nlima,Ncima,"Psf cf");
    FFT2D.fftn2d(Psf, Psf_cf);
    float FluxPsf = 1. / Psf_cf(Nlima/2,Ncima/2).real();
    float Cvgparam = 1.;
    float Level;
    double LambdaCoef;
    float CoefTol = 0.5;
    
    cout << "JAMMAL-BIJAOUI Deconvolution algorithm ... " << endl;
      
    RegulIma CRG;
    CRG.GradOperType = OPER_MARKOV_4; // {OPER_LAPLACIAN,OPER_MARKOV_4,OPER_MARKOV_8};
     
    Resi.alloc(Nlima,Ncima,"resi");
    Grad.alloc(Nlima,Ncima,"resi");
    Mean = average(Data);
    BGR_Ima.init(Mean);
    get_mr_ima(BGR_Ima, NbrPlan);
    // filter_noiter(MR_Data, Result);
    // BGR_Ima = Result;
    Result = BGR_Ima;
    for (i=0; i < Nlima; i++)
    for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
    MultiResol MR_Grad,MR_Ima;
    //cout << "ALLOC " << endl;
    MR_Grad.alloc (Nlima, Ncima, MR_Data.nbr_scale(), 
                   MR_Data.Type_Transform, MR_Data.filter_bank(), NORM_L2);
    
    MR_Ima.alloc (Nlima, Ncima, MR_Data.nbr_scale(), 
                   MR_Data.Type_Transform, MR_Data.filter_bank(), NORM_L2);
	
    //cout << "ITER " << endl;
    for (Iter=0; Iter < MaxFilterIter; Iter++)
    {
       CRG.ima_regul(Result,Grad,RegulParam);
       FFT2D.convolve(Grad, Psf_cf, Grad);
       MR_Grad.transform(Grad);	
       MR_Ima.transform(BGR_Ima);
       
       // cout << "DIFF " << endl;
       for (b=0; b < MR_Data.nbr_band(); b++)
       for (i = 0; i < MR_Data.size_band_nl(b); i++)
       for (j = 0; j < MR_Data.size_band_nc(b); j++) 
                                  MR_Ima(b,i,j) -= MR_Grad(b,i,j);
       for (s = NbrPlan-2; s >= 0; s--)
       {
           NormCoef1 = 1. / (double) POW2(s+2);
           NormCoef2 = (double) POW2(2*(s+1));
          //cout << "LEVEL " << s << endl;
           for (b=3*s; b <= 3*s+2; b++)
           {
	  // cout << "BAND " << b << endl;
              for (i=0; i < MR_Data.size_band_nl(b); i++)
              for (j=0; j < MR_Data.size_band_nc(b); j++)
              {
	      
                 float HaarCoeff = MR_Data(b,i,j);
                 LambdaCoef = (TabImaModel[s])(i,j);
                 Level = get_level(LambdaCoef, s);
		 // Non significant coefficient
		 if (ABS(MR_Data(b,i,j)) < Level)
		 {
		    if (MR_Ima(b,i,j) > Level) MR_Ima(b,i,j) = Level;
		    else if (MR_Ima(b,i,j) < -Level) MR_Ima(b,i,j) = -Level;
		 }
		 else
		 {
		    float LevelMin = MR_Data(b,i,j) - Level*CoefTol;
		    float LevelMax = MR_Data(b,i,j) + Level*CoefTol;
		    if (MR_Ima(b,i,j) < LevelMin) MR_Ima(b,i,j) = LevelMin;
		    else if (MR_Ima(b,i,j) > LevelMax) MR_Ima(b,i,j) = LevelMax;
		 }
              }
	   }
	}
	b = MR_Data.nbr_band()-1;
	for (i=0; i < MR_Data.size_band_nl(b); i++)
        for (j=0; j < MR_Data.size_band_nc(b); j++) 
	                                   MR_Ima(b,i,j) = MR_Data(b,i,j);
	      
	// cout << "RECONS " <<   endl;
	MR_Ima.recons(Grad);
	threshold(Grad);	
	
 	// New solution
        for (i = 0; i < Nlima; i++)
        for (j = 0; j < Ncima; j++) 
        {
           if (BGR_Ima(i,j) > FLOAT_EPSILON) Grad(i,j) /= BGR_Ima(i,j);
           else  Grad(i,j) = 0.;
        }
	FFT2D.convolve_conj(Grad, Psf_cf);
        Result *= Grad;
	threshold(Result);
	
 	// New image
	FFT2D.convolve(Result, Psf_cf, BGR_Ima);
        threshold(BGR_Ima);
	for (i=0; i < Nlima; i++)
        for (j=0; j < Ncima; j++) Resi(i,j) = Data(i,j) - BGR_Ima(i,j);
	get_mr_ima(BGR_Ima, NbrPlan); 
        if (Verbose == True) 
	{
  	   Delta = sigma(Resi);
           cout << "Filter. Iter " << Iter+1 << " Delta = " << Delta << endl << endl;
        }
    }   
}

/****************************************************************************/

// class HAAR_Poisson_Deconv: public MRDeconv {
//   protected:
//    void compute_resi();
//   public:
//      HAAR_Poisson_Deconv() {HF=NULL;init_param();};
//      HFilter *HF;
//      MultiResol *MR_Haar;
// };
// 
// void HAAR_Poisson_Deconv::compute_resi()
// {
//    int i,j;
//    //int Niter=1;
//    psf_convol (Obj, Psf_cf, Resi);
//    if (KeepImagn == True) Imag_n = Resi;
//    
//    // Filter the data using Imag_n as a new pdf
//    // Background image is Imag_n
//    Ifloat  ImagFilter(Nl, Nc, "Info");
//    // a verifier MR_Haar.transform(Imag);
//    BGR_Ima = Imag_n;
//    // TypeBGR = HPF_BGR_IMA;
//    // TypeThreshold = HPF_LAMBDA_GB;
//    // mr_haar_poisson_filter (*MR_Haar, ImagFilter, Niter);
// 
//    for (i=0; i< Nl; i++)
//    for (j=0; j< Nc; j++)
//       Resi(i,j) = ImagFilter(i,j) - Resi(i,j);
// }

/****************************************************************************/
  
type_haar_poisson_threshold  TypeThreshold = DEF_HAAR_POISSON_THRESHOLD;
type_background_model TypeBGR =  DEF_HAAR_POISSON_BGR_MODEL;

const int NBR_OK_METHOD = 3;
type_sb_filter  TabMethod[NBR_OK_METHOD] = {F_HAAR,F_BI2HAAR,F_BI4HAAR};
type_sb_filter HFilterBank = F_BI2HAAR;

/****************************************************************************/

static void usage(char *argv[])
{
    int i;
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
  
    nbr_scale_usage(Nbr_Plan);
    manline();
    
    nsigma_usage(N_Sigma);
    manline();
    
    //converg_param_usage(Epsilon);
    //convergp_param_usage(DEF_EPS_EVENT_FILTERING);
    //manline(); 
    //
    //detect_pos_usage();
    // manline();
    // prec_eps_poisson_usage(DEFAULT_EPSILON);
    // manline();
  
    first_detect_scale_usage();
    manline();
 
    fprintf(OUTMAN, "         [-M BackgroundModel]\n");
    for (i = 0; i <  NBR_HAAR_POISSON_BGR_MODEL; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,  StringHaarBgrModel ((type_background_model )i));
    fprintf(OUTMAN, "              default is %s\n", StringHaarBgrModel ((type_background_model ) TypeBGR));
    manline();    
 
    fprintf(OUTMAN, "         [-T ThresholdingMethod]\n");
    for (i = 0; i <  NBR_HAAR_POISSON_THRESHOLD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,    StringHaarThreshold((type_haar_poisson_threshold) i));
    fprintf(OUTMAN, "              default is %s\n",  StringHaarThreshold ((type_haar_poisson_threshold) TypeThreshold));
    manline();   
    
    fprintf(OUTMAN, "         [-h Haar_FilterBank]\n");
    for (i = 0; i < NBR_OK_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringSBFilter(TabMethod[i]));
    fprintf(OUTMAN, "              Default is %s\n", StringSBFilter(HFilterBank));
    manline();

    fprintf(OUTMAN, "         [-B FileName]\n");
    fprintf(OUTMAN, "             Background Image Model File Name.\n");
    manline();  

//     fprintf(OUTMAN, "         [-i NiterScale]\n");
//     fprintf(OUTMAN, "              Number of iteration per scale reconstruction.\n");
//     fprintf(OUTMAN, "              default is %d\n", Max_Iter);

    manline();  
    fprintf(OUTMAN, "         [-I NiterBgr]\n");
    fprintf(OUTMAN, "              Number of iteration for the iterative background estimation.\n");
    fprintf(OUTMAN, "              default is %d\n",  MaxFilterIter);
    manline(); 
    
    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "              Soft thresholding.\n");
    manline();  
    fprintf(OUTMAN, "         [-L LambdaValue]\n");
    fprintf(OUTMAN, "               Lambda Value for constant background.\n");
    manline();
    fprintf(OUTMAN, "         [-P PsfInFile]\n");
    fprintf(OUTMAN, "               Input Point Spread Function.\n");
    fprintf(OUTMAN, "               If set, a deconvolution is performed.\n");
    manline();
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "               Regularization parameter for the\n");
    fprintf(OUTMAN, "               deconvolution. Default is 0.1.\n");
    manline();
    vm_usage();
    manline();    
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
    // Bool FilterOpt=False;

#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     
    
    /* get options */
    while ((c = GetOpt(argc,argv,"G:Oh:pxdn:s:I:F:vzZ:L:M:T:SB:P:")) != -1) 
    {
	switch (c) 
        {
           case 'G':
 		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad regul paral value: %s\n", OptArg);
		    exit(-1);
		}
  		break; 
	    case 'O': Transform = TO_MALLAT; break;
           case 'h':
                if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of filter bank: %s\n", OptArg);
                    exit(-1);
                    
                }        
                if ((c <= 0) || (c > NBR_OK_METHOD))
                {
                   fprintf(OUTMAN, "Error: bad type of filter bank: %s\n", OptArg);
                   exit(-1);
                    
                }                             
                HFilterBank  = (type_sb_filter) (TabMethod[c-1]);
                break;
           case 'x': OptT = True;break;
           case 'd': Decimated = True; break;
	   case 'S': SoftThreshold = True; break;
           case 'T':
                if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of thresholding: %s\n", OptArg);
                    exit(-1);
                    
                }
                if ((c > 0) && (c <= NBR_HAAR_POISSON_THRESHOLD)) 
		          TypeThreshold = (type_haar_poisson_threshold) (c-1);
                else  
                {
                    fprintf(OUTMAN, "Error: bad type of thresholding: %s\n", OptArg);
                    exit(-1);
                }
                break;
          case 'M':
                if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of model: %s\n", OptArg);
                    exit(-1);
                    
                }
                if ((c > 0) && (c <= NBR_HAAR_POISSON_BGR_MODEL)) 
		          TypeBGR = (type_background_model) (c-1);
                else  
                {
                    fprintf(OUTMAN, "Error: bad type of thresholding: %s\n", OptArg);
                    exit(-1);
                }
                break;           
 	   case 'p': KeepPositivSup = True; break;
	   case 'B':  
	        if (sscanf(OptArg,"%s",Name_BGR_In) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad background file name: %s\n", OptArg);
		    exit(-1);
		}
 		BGROpt = True;
		break;
	   case 'P':
	         if (sscanf(OptArg,"%s",Name_PSF_In) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad background file name: %s\n", OptArg);
		    exit(-1);
		}
 		DeconvOpt = True;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		NscaleOpt = True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                UseNSigma =True;
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
		break;
	   case 'L':
 		if (sscanf(OptArg,"%f",&Lambda) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad lambda value: %s\n", OptArg);
		    exit(-1);
		}
               OptL =True;
 		break; 
	   case 'I':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&MaxFilterIter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad number of iterations: %s\n", OptArg);
		   exit(-1);
		}
                if (Max_Iter <= 0) Max_Iter = DEFAULT_MAX_ITER_FILTER;
                OptI=True;
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
		   exit(-1);
		}
                if (Max_Iter <= 0)   Max_Iter = DEFAULT_MAX_ITER_FILTER;
                MaxIterOpt=True;
		break;
	   case 'e': 
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Epsilon < 0) || (Epsilon > 1.)) 
                           Epsilon = DEFAULT_EPSILON_FILTERING;
		EpsOpt=True;
		break;
           case 'v': Verbose = True;break;
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
       if (DeconvOpt == True) TypeBGR = HPF_BGR_ITER;

       if ((TypeThreshold == HPF_PRESS) 
            && ((SoftThreshold == True) || (UseNSigma == True)))
       {
          cerr << "Error: -S and -s options not valid with PRESS filtering method ... " << endl;
	  exit(-1);
       }
       
       if ((BGROpt == False) && (TypeBGR == HPF_BGR_IMA))
       {
          cerr << "Error: a background image file name must be given ... " << endl;
	  exit(-1);
       }
       if ((BGROpt == True) && (TypeBGR != HPF_BGR_IMA)) TypeBGR = HPF_BGR_IMA;
       
       if ((OptL == True) && (TypeBGR != HPF_BGR_FLAT)) 
       {
          cerr << "Error: -L option only valid with flat background option (-M 1) ..." << endl;
	  exit(-1);
       }
       if ((OptI == True) && (TypeBGR != HPF_BGR_ITER)) 
       {
          cerr << "Error: -I option only valid with option (-M 4) ..." << endl;
	  exit(-1);
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
	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}
/****************************************************************************/

int main(int argc, char *argv[])
{
    int  k;
    Ifloat Data;
 
    /* support image creation */
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);
    
    EpsilonPoisson = (1. - erf((double) N_Sigma / sqrt((double) 2.)));
    if (TypeThreshold == HPF_LAMBDA_GB)
    {
        if ((EpsilonPoisson < TabHaar_MIN_EPS) || (EpsilonPoisson > TabHaar_MAX_EPS))
	{
	   cerr << "Error: EpsilonPoisson = " << EpsilonPoisson << endl;
	   cerr << "       with this thresholding method EpsilonPoisson must verify: " << endl;
	   cerr << "           1e-6 < EpsilonPoisson < 1e-2 " << endl << endl;
	}
    } 
	   
    io_read_ima_float(Name_Imag_In, Data, &Header);
    int N = MIN(Data.nl(), Data.nc());
    Header.origin = Cmd;   

    // Ifloat DT(Data.nl(), Data.nc(),"DT");
    // int number=-1;
//      for (i=0; i < Data.nl(); i++)
//      for (j=0; j < Data.nc(); j++)   
//             DT(i,j) = 0.381668; // poidev(Lambda,&number);
//      io_write_ima_float("tbgr.fits", DT);
 
    if (Nbr_Plan<=1) Nbr_Plan = (int) (log((double) N) / log(2.)) - 2;
    // cout << "TT = " << get_harr_poisson_threshold((double) 0., (double) 1e-3);
if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "N_Sigma = " << N_Sigma << " ==> Epsilon Poisson = " <<  EpsilonPoisson << endl;
    // cout << "Epsilon = " << Epsilon << endl;
    //if (KeepPositivSup == True) 
    //   cout << "Only positive wavelet coefficients are detected" << endl;
    if (FirstScale > 0)
       cout << "Start the detect at scale " << FirstScale+1 << endl;
    cout << "Thresholding type = " << StringHaarThreshold(TypeThreshold) << endl;
    cout << "Background modelization = " << StringHaarBgrModel(TypeBGR) << endl;
    if ((TypeBGR == HPF_BGR_FLAT) && (Lambda > 0))
                cout << "Lambda = " << Lambda << endl;
    if (Decimated  == True) 
    {
         cout << " Orthogonal Haar transform " << endl;
         cout << "Max_Iter per scale = " << Max_Iter << endl;
    }
    if (SoftThreshold == True) 
                cout << " Soft thresholding " << endl;
    else cout << " Hard thresholding " << endl<< endl<< endl;
    cout << "Filter bank = " << StringSBFilter(HFilterBank) << endl;
 }  
 
    Header.origin = Cmd;
    Ifloat Result (Data.nl(), Data.nc(), "Result Filtering");
    check_scale(Data.nl(),  Data.nc(), Nbr_Plan);
    if ((TypeBGR == HPF_BGR_FLAT) && (Lambda < 0))
    {
        Lambda = flux(Data) / (float)(Data.nl()*Data.nc());
        if (Verbose == True)
	  cout << "Estimated Lambda = " << Lambda << endl;
    }
    
    if (Decimated == True) Transform = TO_MALLAT;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    FAS.Verbose = Verbose;
    FAS.alloc(HFilterBank);
    PtrFAS = &FAS;
    MultiResol MR_Data;
    MR_Data.alloc (Data.nl(), Data.nc(), Nbr_Plan, Transform, PtrFAS, NORM_L2);
    MR_Data.transform(Data);

    HFilter HPoisFilter(TypeBGR, TypeThreshold, MR_Data, Name_BGR_In);
    HPoisFilter.Verbose = Verbose;
    HPoisFilter.FirstScale = FirstScale;
    HPoisFilter.SoftThreshold = SoftThreshold;
    HPoisFilter.KeepPositivSup = KeepPositivSup;
    if (OptI == True) HPoisFilter.MaxFilterIter = MaxFilterIter;
    if (MaxIterOpt == True) HPoisFilter.MaxIterPerScaleRec = Max_Iter;
    HPoisFilter.N_Sigma = N_Sigma;
    HPoisFilter.EpsilonPoisson = EpsilonPoisson;
    HPoisFilter.Lambda=Lambda;
     
    if (DeconvOpt == False) 
    {
       if (Verbose == True) cout << "Filtering ... " << endl;
       // HPoisFilter.test(MR_Data);
       if (OptT == False) HPoisFilter.filter(Data, MR_Data, Result);
       else
       {
          HPoisFilter.threshold_all_scale(MR_Data);
          MR_Data.recons(Result);
          threshold(Result);
       }
       Header.bitpix = BP_FLOAT;
       io_write_ima_float(Name_Imag_Out, Result, &Header);
    }
    else
    {
       Ifloat Psf;
       io_read_ima_float(Name_PSF_In, Psf);
       HPoisFilter.deconv(Data, Psf, MR_Data, Result, RegulParam);
       // HPoisFilter.jammal_deconv(Data, Psf, MR_Data, Result, RegulParam);
       Header.bitpix = BP_FLOAT; 
       io_write_ima_float(Name_Imag_Out, Result, &Header);
   }
    exit(0);
} 



//     float *LevelQ = new float[MR_Data.nbr_band()];
//     for (b=0; b < MR_Data.nbr_band()-1; b++)
//     { 
//        int s = b / 3; 
//        double NormCoef1 = 1. / (double) POW2(s+2);
//        double NormCoef2 = (double) POW2(2*(s+1));
//        double Lambdaj = NormCoef2*Lambda;
//        double z2 = N_Sigma*N_Sigma;
//        double z4 = z2*z2;
//        double Level = NormCoef1*(z2+sqrt(z4+4.*Lambdaj*z2));
//        LevelQ[b] = Level;
//        cout << "Band " << b+1 << " Level = " << Level << " coef = " << NormCoef1 << endl;
//        for (i=0; i < MR_Data.size_band_nl(b); i++)
//        for (j=0; j < MR_Data.size_band_nc(b); j++)
//        {
//        // MR_Data(b,i,j) /= NormCoef1;
//           if (KeepPositivSup == True)   
// 	  {
// 	     if (ABS(MR_Data(b,i,j)) < Level) MR_Data(b,i,j) = 0.; 
// 	  }
// 	  else if (MR_Data(b,i,j) < Level) MR_Data(b,i,j) = 0.;
//        }
//     }
//    MR_Data.write("xx_h");
