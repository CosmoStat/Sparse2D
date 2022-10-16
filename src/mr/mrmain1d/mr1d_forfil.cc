/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  28/08/00
**    
**    File:  mr_fcast.cc
**
*******************************************************************************
**
**    DESCRIPTION  Prediction using the undecimated Haar transform. 
**    ----------- 
**                 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "IM1D_IO.h"
#include "MatrixOper.h"
#include "MR1D_Predict.h"
#include "NeuNet.h"
#include "Usage.h"
#include "CErf.h"
#include "CMem.h"
#include "MR1D_NoiseModel.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];
char InfoFile[256];
char Name_Orig[256];

int NbrScale=5;           // number of scales
float NSigma=3.;   // number of sigma (for the noise) 
float SigmaNoise=0.;       // noise standard deviation
float NoiseAR = 0.;

type_noise Stat_Noise = DEFAULT_STAT_NOISE;   // type of noise
type_trans_1d Transform=TO1_PAVE_HAAR;    // type of transform
// type_trans_1d Transform=TO1_PAVE_B3SPLINE;
float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean  

Bool UseNSigma=False;
// Bool SetPlan = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
int Pos = -1;
 
Bool HardThreshold=False;
Bool EntrThreshold=False;
int NPredict=0;
int NbrEval=0;

Bool UseNoise=False;
extrapol_type  ExtrapolType = EX_MR_AR_PRED2; // EX_POLYNOM;
Bool Simu=False;

int NbrLearning =  0;

#define IS_AR4 1
#if IS_AR4
#define NBR_ARCoefSimu 4
float TabARSimu[NBR_ARCoefSimu] = {0.5,-0.5,-0.1,0.3};
#else
#define NBR_ARCoefSimu 1
float TabARSimu[NBR_ARCoefSimu] = {0.8};
#endif

int PredARNbr=-1;

int NpSimu=512;
unsigned int InitRnd = 100;
Bool NoLag=False;
Bool Stationary=True;
Bool OnlyNoise = False;
int MaxNbrAR = 10;
ar_order_detect_type AR_Order_Detect=AR_ORDER_BIC;

Bool MyTest = False;
Bool EvaluationMode=True;
int NbrPixTendancy = 100;
int NbrPixLastPred = 5;
int LastNPix=0;
int MinAR_OrderPerScale = 1;

pred_mirror_type TPredBord=P_POL2;
Bool UseB3=False;
Bool UseTend=False;
Bool DebugPred = False; 

#define TEST_B3SPLINE 0 
#define DEF_ITER 10

Bool WriteInfo = False;

float Noise_Ima=0.;       // noise standard deviation
Bool Orig = False;
float RegulWave=0.1;
int Max_Iter=1;      // maximum number of iterations
Bool UpDateEstimParam = False;
Bool RatioNoiseAR = False;

/*********************************************************************/

// static int max_scale_number (int Nc)
// {
//     int Nmin = Nc;
//     int ScaleMax;
// 
//     ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
//     return (ScaleMax);
// }
//  

/***************************************************************************/

static void usage(char *argv[])
{
   // int i;

    fprintf(OUTMAN, "Usage: %s options input output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
  
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    fprintf(OUTMAN, "             Default is 5.\n");
    manline();
    
    fprintf(OUTMAN, "         [-g SigmaNoise]\n");
    fprintf(OUTMAN, "             Noise standard deviation.\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();
    
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "             Regularization parameter.\n");
    fprintf(OUTMAN, "             Default is %f.\n", RegulWave);
    manline();
    fprintf(OUTMAN, "         [-A SigmaARInnovation]\n");
    fprintf(OUTMAN, "             AR innovation standard deviation.\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();

    fprintf(OUTMAN, "         [-a AR_Order]\n");
    fprintf(OUTMAN, "             AR oder used for the prediction. Default is automatically estimated.\n");
    manline();

    fprintf(OUTMAN, "         [-O Estimation_AR_Order_Method]\n");
    fprintf(OUTMAN, "             1: AIC \n");
    fprintf(OUTMAN, "             2: AICC\n");
    fprintf(OUTMAN, "             3: BIC\n");
    manline();

    fprintf(OUTMAN, "         [-h]\n");
    fprintf(OUTMAN, "             Non stationary signal. Default is stationary.\n");
    manline();

    fprintf(OUTMAN, "         [-p Number_of_Predict]\n");
    fprintf(OUTMAN, "             Number of prediction. Default is %d.\n", NPredict);
    manline();  
 
    fprintf(OUTMAN, "         [-m MinAROrder]\n");
    fprintf(OUTMAN, "             Default is %d.\n", MinAR_OrderPerScale);
    manline();
        
    fprintf(OUTMAN, "         [-M MaxAROrder]\n");
    fprintf(OUTMAN, "             Default is %d.\n", MaxNbrAR);
    manline();

    fprintf(OUTMAN, "         [-w InfoFileName]\n");
    fprintf(OUTMAN, "             Save some information in a file.\n");
    manline();
            
    fprintf(OUTMAN, "         [-L NPix]\n");
    fprintf(OUTMAN, "             Analyse the NPix last pixels. Default is all pixels.\n");
    manline();

    pred_border_usage(TPredBord);  // -B option
    manline();
    
    fprintf(OUTMAN, "        [-T Poly_Nbr_Pix]\n");
    fprintf(OUTMAN, "            Number of pixels used for the polynomial extrapolation.\n");
    fprintf(OUTMAN, "             Default is %d.\n",  NbrPixLastPred);
    manline();
    
    max_iter_usage(Max_Iter);
    manline();
        
    verbose_usage();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    /* get options */    // WN:SI: b X
    while ((c = GetOpt(argc,argv,"RA:Ui:G:g:L:B:T:xO:M:hla:p:m:n:v")) != -1) 
    {
	switch (c) 
        {
	   case 'R': RatioNoiseAR= True; break;
	   case 'U': UpDateEstimParam = True; break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
                   usage(argv);
		}
                if (Max_Iter <= 0)  Max_Iter = DEF_ITER;
		break;
           case 'w': 
	             if (sscanf(OptArg,"%s",InfoFile) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
	                exit(-1);
		     }
	             WriteInfo = True; 
	             break;
	   case 'x': DebugPred = True; break;
           case 'L':  
		if (sscanf(OptArg,"%d",&LastNPix) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
 		break;
	   case 'T':  
		if (sscanf(OptArg,"%d",&NbrPixLastPred) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
 		break;
           case 'X': MyTest = True; break;
	   case 'O':  
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad order method: %s\n", OptArg);
	            exit(-1);
		}
                AR_Order_Detect = (ar_order_detect_type)(c-1);
		break;
	   case 'm':  
		if (sscanf(OptArg,"%d",&MinAR_OrderPerScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad minimum order value: %s\n", OptArg);
	            exit(-1);
		}
		break;
	   case 'M':  
		if (sscanf(OptArg,"%d",&MaxNbrAR) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad maximum order value: %s\n", OptArg);
	            exit(-1);
		}
		break;
           case 'h': Stationary = False; break;
           case 'l': NoLag = True; break;
	   case 'I':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		    exit(-1);
		}
                InitRnd = (unsigned int) c;
		break;
 	   case 'S':Simu = True; break; 
	   case 'N':  
		if (sscanf(OptArg,"%d",&NpSimu) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of pointsn: %s\n", OptArg);
	            exit(-1);
		}
		break;
	   case 'f':
	   	if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad filtering method: %s\n", OptArg);
	            exit(-1);
		}
		if (c == 2) HardThreshold = True;
                else if (c == 3) EntrThreshold = True;
                else if (c != 1)
		{
		    fprintf(OUTMAN, "Error: bad  filtering method: %s\n", OptArg);
	            exit(-1);
		}
                if ((HardThreshold == True) || (EntrThreshold == True)) UseNoise = True;
		break;
	   case 'a':
	   	if (sscanf(OptArg,"%d",&PredARNbr) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
 		break;
	   case 'p':
	   	if (sscanf(OptArg,"%d",&NPredict) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number: %s\n", OptArg);
	            exit(-1);
		}
		if ( NPredict < 1)
		{
		    fprintf(OUTMAN, "Error: bad  number: %s\n", OptArg);
		    fprintf(OUTMAN, "        0 < x  \n");
	            exit(-1);
		}
		break;
           case 'B':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad border type: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_PRED_BORD)) 
                                       TPredBord  = (pred_mirror_type) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad border type: %s\n", OptArg);
	            exit(-1);
 		}
		break;		
	   case 'n':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
                if ((NbrScale <= 1) || (NbrScale > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    usage(argv);
		}
 		break;
        case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                Stat_Noise = NOISE_GAUSSIAN;
                break;
        case 'A':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&NoiseAR) != 1) {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
	case 'G':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&RegulWave) != 1) {
                    fprintf(OUTMAN, "bad regularization parameter: %s\n", OptArg);
                    usage(argv);
                }
                break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) 
                {
		    fprintf(OUTMAN, "bad NSigma: %s\n", OptArg);
		    usage(argv);
		}
                if (NSigma <= 0.) NSigma = 3.;
		UseNSigma = True;
		break;
	  case 'v': Verbose = True; break;
  	    case '?': usage(argv); break;
		}
	}

        /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
        else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
        else usage(argv);

	if (OptInd < argc) 
	{
	   strcpy(Name_Orig, argv[OptInd++]);
           Orig = True;
	}
	
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}
 
/*********************************************************************/

class MR_FIL_PRED: public MR_PRED  
{ 
    CMemWave MemWObj;
    void get_new_haar_coef(float NewVal, fltarray & Signal, int NbrScale, 
                           int Np, fltarray & TabHaar);
    void set_ar_model(fltarray & Signal, MR_1D & MR_Data, 
                      fltarray & ARModel, int User_AR_Order);
   public:
     float UserSigmaNoise;
     float UserSigmaARNoise;
     float SigmaNoise;
     float SigmaARNoise;
     float RegulWave;
     int MaxIter; // maximum number of iteration for the filtering
                  // Default is 1, i.e. no filtering
     float ConvgParam; // convergence parameter for the iterative filtering
     MR1DNoiseModel *NoiseModel;
     MR_FIL_PRED():MR_PRED() {NoiseModel=NULL;RegulWave=1.;MaxIter=1;SigmaARNoise=0.;
                              ConvgParam= 1e-4;UserSigmaNoise=0.;
			      UserSigmaARNoise=0.;SigmaNoise=0.;}
			      
     double pred_filter(fltarray &Signal, fltarray &TabSignalFil,
                              fltarray &TabSignalPred, MR_1D & MR_Data, 
                              int NbrPredict, int User_AR_Order);
     double ar2_pred_filter(fltarray & Signal, MR_1D & MR_Data, 
                            fltarray & TabFilter, MR_1D & MR_Fil,
                          fltarray &TabSignalPred, int User_AR_Order);

    double estime_error_pred(fltarray & Signal, 
           fltarray   &TabDelta,  MR_1D & MR_Data,  fltarray &ARModel, 
           fltarray &TabSignalPred, fltarray & NewFil, 
           MR_1D & MR_NewFil, fltarray & TabResi, fltarray & TabSigma, double & Err);

   double compute_predict_signal(fltarray & Signal, 
          MR_1D & MR_Data, fltarray & PredSignal, fltarray &ARModel,
	    fltarray &LastScale);
	       
    double estime_sigma_noise(fltarray & Signal, MR_1D & MR_Data,
            fltarray & TabFilter, MR_1D & MR_Fil, 
            fltarray &TabSignalPred, fltarray & TabAlpha, fltarray &ARModel,
	    fltarray &LastScale);
	    
	       float get_pred(MR_1D & MR_Data, int NusedScale, int Np, fltarray & ARModel);
    ~MR_FIL_PRED(){}
};


/*********************************************************************/

void MR_FIL_PRED::set_ar_model(fltarray & Signal, MR_1D & MR_Data, 
                               fltarray & ARModel, int User_AR_Order)
{
   int b,Step= (NoLag == False) ? 2: 1;
   int NbrScale = MR_Data.nbr_scale();
   int Nx= MR_Data.size_ima_np();
   // int NPredict = TabSignalPred.nx() - Nx;
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
    
   fltarray Scale(Nx);
     
   // Find the best AR order per scale
   if (AR_AutoDetect == True)
   {
      // cout << " NbrTraining = " << NbrTraining << endl;
      Step = (NoLag == False) ? 2: 1;
      for (b=0; b < NusedScale; b++)
      {
         int BestNbrAR;   
         MR_Data.scale(Scale, b);
         get_best_ar_model(Scale, NbrTraining, BestNbrAR, ARModel, Step,b);
         TabAR_OrderPerScale(b) = MAX(MinAR_OrderPerScale,BestNbrAR);
         if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
      }
   }
   else  for (b=0; b < NusedScale; b++) TabAR_OrderPerScale(b) = User_AR_Order;
 
   int NArParam = 0;
   for (b=0; b < NusedScale; b++) NArParam += TabAR_OrderPerScale(b);
   get_mr_ar2_model(Signal, MR_Data, NbrTraining, NArParam, ARModel);

   if (Verbose == True)
   {
       cout << "Total number of AR coefficients = " << NArParam << endl;
       if (Verbose == True) cout << "  Model AR = " << endl;
       int IndCoef = 0;
       for (b=0; b < NusedScale; b++)
       {
          cout << "    Scale " << b+1 << " NCoef = " <<  TabAR_OrderPerScale(b) << " Coef = ";
          for (int k=0; k<  TabAR_OrderPerScale(b); k++)  cout << " " << ARModel(IndCoef++);   
          cout << endl;
       }
   }
}

/*********************************************************************/

float MR_FIL_PRED::get_pred(MR_1D & MR_Data, int NusedScale, int Np, fltarray & ARModel)
{
   int Step = (NoLag == False) ? 2: 1;
   float Pred = 0.;  
   int IndCoef = 0;

   for (int b=0; b < NusedScale; b++)
   {
        for (int k=1; k<= TabAR_OrderPerScale(b); k++) 
              Pred += ARModel(IndCoef++) * MR_Data(b, Np-1-(Step*(k-1)));
        if ((NoLag == False) && (b != NbrScale-2)) Step *= 2;
   }
   return Pred;
}

/*********************************************************************/

void MR_FIL_PRED::get_new_haar_coef(float NewVal, fltarray & Signal, int NbrScale, int Np, fltarray & TabHaar)
// Signal = tab(0...Np-1)
{
   int i,b;
   int SupportWavelet=1;
   int Pos=0;
   TabHaar(0) = NewVal;
   
   // Calculate the scaling cofficient at each scale
   // cout << "II " << 0.5*(NewVal-Signal(Np-1)) << endl;
   for (b=1; b < NbrScale; b++)
   {
      float Sum=0.;
      for (i=0; i < SupportWavelet; i++)  Sum += Signal(MAX(0,Np-1-i-Pos));
      Sum /= (float) SupportWavelet;
      TabHaar(b) = (TabHaar(b-1) + Sum)*0.5;
      // cout << "Pos = " << Pos << " SupportWavelet " << SupportWavelet << " Np = " << Np << endl;
      Pos += SupportWavelet;
      SupportWavelet *=2;
   }
   // Calculate the wavelet cofficient at each scale
   for (b=0; b < NbrScale-1; b++) TabHaar(b) = TabHaar(b) - TabHaar(b+1);
}

/*********************************************************************/


/*********************************************************************/

double MR_FIL_PRED::estime_error_pred(fltarray & Signal, 
   fltarray   &TabAlpha,  MR_1D & MR_Data,  fltarray &ARModel, 
   fltarray &TabSignalPred, fltarray & NewFil, 
   MR_1D & MR_NewFil, fltarray & TabResi, fltarray & TabSig, double & ErrPredFil)
{
   int DistPred=0;
   int i,b; 
   int NbrScale = MR_Data.nbr_scale();
   int Nx= MR_Data.size_ima_np();
   int Np;
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
   fltarray TabHaar(NbrScale);
   fltarray TabHaarResi(NbrScale);
   fltarray TabWaveNorm(NbrScale);
   fltarray TabWiener(NbrScale);
   fltarray LastScale(Nx);
   fltarray ObjScale(Nx);
   fltarray Resi(Nx);
   fltarray ResiPred(Nx);
   double ErrPred = 0;
   double ResiVal, ErrResi = 0;
   ErrPredFil = 0;

   // cout << "Nx = " << Nx  << " signal.nx = " << Signal.nx() << endl;
   
   TabWaveNorm(0) = 1. / sqrt(2.);
   for (b=0; b < NbrScale; b++) 
   {
      TabWaveNorm(b) = (b == 0) ? TabWaveNorm(0): TabWaveNorm(b-1) * TabWaveNorm(0);
      MR_Data.scale(ObjScale,b);
      float SigmaObj = ObjScale.sigma();
      float Sig = SigmaNoise*TabWaveNorm(b);
      if (SigmaObj > Sig) SigmaObj = SigmaObj*SigmaObj - Sig*Sig;
      else SigmaObj = 0; 
      TabWiener(b) =  SigmaObj / (SigmaObj + Sig*Sig); 
      //cout << "Sig = " << Sig << " ObjScale.sigma() = " << ObjScale.sigma() << endl;
      //cout << "SigmaObj= " << SigmaObj << endl;
      //cout << "TW " << b+1 << " " << TabWiener(b) << endl;    
   }      
   TabResi.init();	
   TabSig.init();
   NewFil.init();
   int NpUsed = Nx - (int) TabAR_OrderPerScale.max();
   for (i= (int) TabAR_OrderPerScale.max(); i < Nx; i++)
   {
       Np = i; // Np = Nx - NbrEvaluate + i;
       // Do the MAR prediction	from the filtered data
       TabSignalPred(Np) = get_pred(MR_NewFil, NusedScale, Np, ARModel);
       if (Stationary == False) 
            TabSignalPred(Np) += extr_last_scale(LastScale, Np, DistPred);
       ResiVal = Signal(Np) - TabSignalPred(Np);
       ErrPred += ResiVal*ResiVal;
       // Get the Haar coefficients of the MAR prediction 
       get_new_haar_coef(TabSignalPred(Np), NewFil, NbrScale, Np, TabHaar);
       get_new_haar_coef(ResiPred(Np), ResiPred, NbrScale, Np, TabHaarResi);

          // Multiscale entropy filtering
      NewFil(Np) = 0;
      for (b=0; b < NbrScale-1; b++) 
      {
         float Wm = TabHaar(b);   // Predicted Model
	 float Wd = MR_Data(b,Np); // Data 
	 float Sig = SigmaNoise*TabWaveNorm(b); // noise stand. deviation
	                                       // Multiscale entropy filtering
         float SigAR = SigmaARNoise*TabWaveNorm(b);
	 // MR_NewFil(b,Np) *= TabWiener(b);
         MR_NewFil(b,Np) = MemWObj.filter(Wd, TabAlpha(b), Sig, SigAR, Wm);
	     // MR_NewFil(b,Np) = Wm;
	     // float ResiVal = Wd  - MR_NewFil(b,Np);
	 float ResiVal = TabHaarResi(b);
   	 // MR_NewFil(b,Np) = Wm + soft_threshold(Wd-Wm, TabWaveNorm(b)*RegulWave*SigmaARNoise);
	 NewFil(Np) += MR_NewFil(b,Np);
	 TabResi(b) += ResiVal*ResiVal;
	 TabSig(b) += Sig*Sig;
      }
      NewFil(Np) += MR_Data(NbrScale-1,Np);
      ResiPred(Np) = Signal(Np) - NewFil(Np);
      ErrResi += ResiPred(Np)*ResiPred(Np);
      ErrPredFil += (NewFil(Np) - TabSignalPred(Np)) * (NewFil(Np) - TabSignalPred(Np));
      MR_NewFil(NbrScale-1,Np) = MR_Data(NbrScale-1,Np);
          //TabFilter(Np) += MR_Obj(NbrScale-1,Np);
          //MR_Fil(NbrScale-1,Np) = MR_Obj(NbrScale-1,Np);
          //  TabFilter(Np) += TabHaar(NbrScale-1);
          //  MR_Fil(NbrScale-1,Np) = TabHaar(NbrScale-1);
   } // END for i= 1,Nx
   
   for (b=0; b < NbrScale-1; b++)
   {
      TabSig(b) = sqrt(TabSig(b) / (float) NpUsed);
      TabResi(b) = sqrt(TabResi(b) / (float) NpUsed);
   }
   ErrResi /= (float) NpUsed;    
   ErrPred /= (float) NpUsed;
   ErrPredFil /= (float) NpUsed;
   if (ErrPred > 0) ErrPred = sqrt(ErrPred);
   if (ErrResi > 0) ErrResi = sqrt(ErrResi);
   if (ErrPredFil > 0) ErrPredFil = sqrt(ErrPredFil);
   // double ErrorMin = ABS(SigmaNoise - ErrResi);
   return ErrPred;
}
   
/*********************************************************************/

double MR_FIL_PRED::compute_predict_signal(fltarray & Signal, 
          MR_1D & MR_Data, fltarray & PredSignal, fltarray &ARModel,
	    fltarray &LastScale)
{
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
   int i,Nx= MR_Data.size_ima_np();
   int NpUsed = Nx - (int) TabAR_OrderPerScale.max();
   float Resi, Err=0.;

   for (i= (int) TabAR_OrderPerScale.max(); i < Nx; i++)
   {
      int Np = i; // Np = Nx - NbrEvaluate + i;
      int DistPred=0;
      PredSignal(Np) = get_pred(MR_Data, NusedScale, Np, ARModel);
      if (Stationary == False) 
                PredSignal(Np) += extr_last_scale(LastScale, Np, DistPred);
      Resi  = Signal(Np) - PredSignal(Np);
      Err += Resi*Resi;
   }
   Err /= (float) NpUsed;  
   if (Err > 0) Err = sqrt(Err);
   else Err= 0.;
   return Err;  
}

/*********************************************************************/

double MR_FIL_PRED::estime_sigma_noise(fltarray & Signal, MR_1D & MR_Data,
            fltarray & TabFilter, MR_1D & MR_Fil, 
            fltarray &TabSignalPred, fltarray & TabAlpha, fltarray &ARModel,
	    fltarray &LastScale)
{
     // Estimation of the noise level
     int b,i; // Nx= MR_Data.size_ima_np();
     int NStepEval = 50;
     fltarray TabStepEval(NStepEval);
     int NbrScale = MR_Data.nbr_scale();
     // int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
     fltarray TabEval(NStepEval);
     fltarray TabResi(NbrScale);
     fltarray TabSig(NbrScale);
     float Err, SigmaMin = 0, ErrMin=0;
     float SigmaMax = Signal.sigma();
     float DeltaStep = SigmaMax / (float) (NStepEval-1);
     float VarData = Signal.sigma()*Signal.sigma();
     double SigmaMinPredFil = SigmaARNoise, ErrPredFil = 0.;
     fltarray TabAuto(1);
     fltarray Residual(TabSignalPred.nx());
     
     for (int It = 0; It < NStepEval; It++)
     {
       SigmaNoise = It * DeltaStep;
       for (b = 0; b < TabAlpha.nx(); b++) TabAlpha(b) = RegulWave;
       
       Err = estime_error_pred(Signal, TabAlpha,  MR_Data, ARModel, 
                    TabSignalPred,  TabFilter, MR_Fil, TabResi, TabSig, ErrPredFil);
       if (UserSigmaARNoise  <= 0) SigmaARNoise = ErrPredFil;
       float ResiFil = TabFilter.sigma();
       ResiFil = VarData - (ResiFil*ResiFil + SigmaNoise*SigmaNoise);
       float RR = ABS(VarData - (ErrPredFil*ErrPredFil+ SigmaNoise*SigmaNoise));
       // float RR = ResiFil;
       //cout << "  SigmaNoise = " << SigmaNoise << " SigmaAR = " << ErrPredFil;  
       //cout << "   J = " << RR  << endl;
         //cout << "Sigma Noise Est = " << sqrt(ErrPred*ErrPred - VarPredFil) << endl;   
       for (i= 0; i <Signal.nx(); i++)
                           Residual(i) = Signal(i) - TabSignalPred(i);
			   
       if ((It == 0) || (RR < ErrMin))
       {
           ErrMin = RR;
	   SigmaMin = SigmaNoise;
 	   SigmaMinPredFil = ErrPredFil;
       }
       // MR_Fil.transform(TabFilter);
       // set_ar_model(TabFilter, MR_Fil, ARModel, -1);
       // autocor1d(Residual, TabAuto, 4, MR_Data.size_ima_np());
       // cout << It+1 << "  SigmaNoise = " << SigmaNoise << "   Err = " << Err  <<  "  Corr = " << (ABS(TabAuto(1))+ ABS(TabAuto(2))+ABS(TabAuto(3))) <<  " RR = " << RR << endl;
     } 
     cout << "MIN Err ==> SigmaNoise = " << SigmaMin << " SigmaAR = " << SigmaMinPredFil << endl;
     SigmaNoise = SigmaMin;
     if (UserSigmaARNoise  <= 0) SigmaARNoise = SigmaMinPredFil;
     // MR_Fil.transform(Signal);
     // set_ar_model(Signal, MR_Fil, ARModel, -1);
     // float Err  = estime_error_pred(Signal, TabAlpha,  MR_Data, ARModel, 
     //               TabSignalPred,  TabFilter, MR_Fil, TabResi, TabSig, ErrPredFil);
     return SigmaMin;
}

/*********************************************************************/

double MR_FIL_PRED::ar2_pred_filter(fltarray & Signal, MR_1D & MR_Data,
                            fltarray & TabFilter, MR_1D & MR_Fil, 
                          fltarray &TabSignalPred, int User_AR_Order)
{
   int i,b;
   int NbrScale = MR_Data.nbr_scale();
   int Nx= MR_Data.size_ima_np();
   float Pred; 
   int NiterParam = 0;
   int NusedScale = (Stationary == True) ? NbrScale: NbrScale-1;
   int Nbr_Band = NbrScale;
   double ErrPredFil,ErrPred;
    
   fltarray  RegulMin(Nbr_Band); 
   fltarray  RegulMax(Nbr_Band);    
   fltarray  TabDelta(Nbr_Band);  
   fltarray  TabStep(Nbr_Band); 
   fltarray TabResi(NbrScale);
   fltarray TabResi1(NbrScale);
   fltarray TabDirBand(NbrScale);
   fltarray TabOldResi(NbrScale);
   fltarray TabSig(NbrScale);
   fltarray TabAlpha(NbrScale);
   // MR_1D MR_NewFil(Nx, MR_Data.Type_Transform, "MR_Data", MR_Data.nbr_scale());
   

   if (Debug == True) 
     cout << " AR2 PRED " << NbrEvaluate << endl;
   fltarray LastScale(Nx);
   fltarray Resi(Nx);
   fltarray ResiPred(Nx);
   fltarray ARModel;
   
   // Set the MAR model order from the NbrTraining first pixels
   set_ar_model(Signal, MR_Data, ARModel, User_AR_Order);

   //  Initialisazation 
  for (i=0; i< Nx; i++)
  {
      for (b=0; b < NbrScale; b++) MR_Fil(b,i) = MR_Data(b,i);
      TabSignalPred(i) =  TabFilter(i) = Signal(i);
  }
 
  // Do the evaluation 
  if (Stationary == False) MR_Data.scale(LastScale,NbrScale-1);

  for (b = 0; b < Nbr_Band; b++) TabAlpha(b) = RegulWave;
  if (UserSigmaARNoise  <= 0) 
      SigmaARNoise = compute_predict_signal(Signal,MR_Data,Resi,ARModel,LastScale);
  else SigmaARNoise = UserSigmaARNoise;
    
 
  if (UserSigmaNoise <= 0)
  {
     SigmaNoise = estime_sigma_noise(Signal,MR_Data,TabFilter,MR_Fil, 
                                  TabSignalPred,TabAlpha, ARModel, LastScale);
  }
  else SigmaNoise = UserSigmaNoise;   
                      
  // if (Verbose == True)
         cout << "SigmaNoise = " << SigmaNoise << " SigmaAR = " << SigmaARNoise << endl;
  

  do
  {
     for (b = 0; b < Nbr_Band; b++)
     {
       if (RatioNoiseAR == True)
            TabAlpha(b) = RegulWave * SigmaNoise / SigmaARNoise;
       else TabAlpha(b) = RegulWave;
     }
   
      ErrPred = estime_error_pred(Signal, TabAlpha,  MR_Data, ARModel, 
                    TabSignalPred,  TabFilter, MR_Fil, TabResi, TabSig, ErrPredFil);
      if (UserSigmaARNoise  <= 0) SigmaARNoise = ErrPredFil;
	// cout << "Prediction error = " << ErrPred << endl;
      // TabFilter.info();
          
     if (UpDateEstimParam == True)
     {
         if (Verbose == True) cout << " Update MAR coeff ... " << endl;
	 MR_Fil.transform(TabFilter);
         set_ar_model(TabFilter, MR_Fil, ARModel, User_AR_Order);
         if (UserSigmaNoise <= 0)
         {
             SigmaNoise = estime_sigma_noise(Signal,MR_Data,TabFilter,MR_Fil, 
                                  TabSignalPred,TabAlpha, ARModel, LastScale);
         }
     }
     NiterParam ++;
     if (Verbose == True)
      cout << NiterParam << ": Prediction error = " << ErrPred << " SigmaNoise = " << SigmaNoise << " SigmaARNoise " << SigmaARNoise << endl;
  } while  (NiterParam < MaxIter);
  // TabFilter.info();
  Resi = Signal;
  Resi -= TabFilter;         
  // if ((Verbose == True) && (DoFilterIter == True))
  //    for (b = 0; b <  Nbr_Band-1; b++) 
  //     cout << "band " << b+1 << " Optimal Regul. Param = " <<  TabAlpha(b) << endl;
  // float VarData = Signal.sigma();

  MR_Fil.transform(TabFilter);
  float ResiFil = compute_predict_signal(TabFilter,MR_Fil,Resi,ARModel,LastScale);

  cout << "Sigma Signal = " << Signal.sigma() << endl;
  cout << "Sigma Sol = " << TabFilter.sigma() << endl;
  cout << "Sigma Noise = " << SigmaNoise << endl;
  cout << "Sigma AR = " << SigmaARNoise << endl;
  cout << "Sigma (Data-Filter) = " << ResiFil << endl;
  cout << "Sigma (Data-Pred)   = " << ErrPred  << endl;
  Resi = TabFilter;
  Resi -= TabSignalPred;
  cout << "Sigma (Filter-Pred)   = " << Resi.sigma()  << endl;

    
  Pred = get_pred(MR_Fil, NusedScale, Nx, ARModel); 
  if (Stationary == False) Pred += extr_last_scale(LastScale, Nx, PredDistance);
 
  if (Debug == True) cout << "End AR 2 predict " << Pred << endl;
  return Pred;
}

/***************************************************************************/

double MR_FIL_PRED::pred_filter(fltarray &Signal, fltarray &TabSignalFil, 
                              fltarray &TabSignalPred, MR_1D & MR_Data, 
                              int NbrPredict, int User_AR_Order)
{
   double ErrPred=0.,Pred = 0.;
   int i,Nx=Signal.nx();
   TabSignalPred.alloc(NbrPredict+Nx);
   TabSignalFil.alloc(NbrPredict+Nx);
   MR_1D MR_FIL (NbrPredict+Nx, MR_Data.Type_Transform, "MR_Data", MR_Data.nbr_scale());

   for (i=0; i < Nx; i++) TabSignalPred(i) = Signal(i);
   if (EvaluationMode == True)
   {
      // By default, separate the serie into two parts: on part for the training
      // and one part the evaluation
      if (LastTraining  <= 0) LastTraining = Nx/2 - 1;
      if (FirstEvaluate <= 0) FirstEvaluate  = Nx/2;
      if (LastEvaluate  <= 0) LastEvaluate = Nx -1;
      NbrTraining = LastTraining-FirstTraining+1;
      NbrEvaluate = LastEvaluate-FirstEvaluate+1;
      for (i=FirstEvaluate; i <= LastEvaluate; i++) TabSignalPred(i) = 0.;
   }
   else
   {
      if (LastTraining <= 0) LastTraining = Nx - 1;
      NbrTraining = LastTraining-FirstTraining+1;
      LastEvaluate = -1;
      NbrEvaluate = FirstEvaluate = 0;
   }
   TabAR_OrderPerScale.alloc(MAX(MR_Data.nbr_scale(),1)); 
   if (AR_AutoDetect == False)
   {
      if (User_AR_Order < 0) 
      {
           cout << "Error: User_AR_Order is not correctly set ... " << endl;
           exit(-1);
      }
      if (User_AR_Order > MaxNbrAR) MaxNbrAR = User_AR_Order;
   }
   
   if (Verbose == True) 
   {
      if (AR_AutoDetect == False) cout << "ARNbr  = " << User_AR_Order  << endl;
      cout << "Np   = " << Nx  << endl;
      cout << "NbrTraining    = " << NbrTraining   << endl;
      cout << "NbrEvaluate     = " <<  NbrEvaluate   << endl;
      cout << "FirstTraining    = " << FirstTraining   << endl;
      cout << "FirstEvaluate     = " <<  FirstEvaluate   << endl;
      cout << "LastTraining     = " << LastTraining   << endl;
      cout << "LastEvaluate     = " <<  LastEvaluate   << endl;
      if (AR_AutoDetect == True) cout << "Automatic AR order estimation " << endl;
      if (Stationary == True) cout << "Stationary signal " << endl;
      if (AR_AutoDetect == True)      
        switch (AR_Order_Detect)
        {
         case AR_ORDER_AIC: 
                cout << "  Use AIC method criterion " << endl; break;
         case AR_ORDER_AICC: 
                cout << "  Use AICC method criterion " << endl; break;
         case AR_ORDER_BIC: 
                 cout << "  Use BIC method criterion " << endl; break;
         default: cout << "Error: unknown criterion ... " << endl;
                  exit(-1);
         }
   }
   int p;
   fltarray ARModel;
   int Nloop = MAX(1,NbrPredict);
   Bool MakeEval = EvaluationMode;
   
   MR_Data.transform(Signal);
   for (p=0; p < Nloop; p++)
   {
      PredDistance = p;      
      Pred = ar2_pred_filter(Signal, MR_Data, TabSignalFil, MR_FIL, 
                         TabSignalPred, User_AR_Order);
      if (((Nx+p-1) < 0) || ((Nx+p-1) >= TabSignalPred.nx()))
      {
	       cout << "Error: Nx = " << Nx << " p = " << p;
	       cout << " TabSignalPred.nx() = " << TabSignalPred.nx() << endl;
	       exit(-1);
      }
      if (NbrPredict > 0) TabSignalPred(Nx+p) = Pred;
    }
    if ((Verbose == True) && (Nloop == 1))
	          cout << "Predicted value = " << Pred << endl;
   EvaluationMode = MakeEval;
   
   // Error evaluation
//   double ErrPred=0.,Err,Pred0=0.;
//    if (EvaluationMode == True)
//    {
//      for (i=FirstEvaluate; i <= LastEvaluate; i++)
//      {
//         Err = Signal(i) - TabSignalPred(i);
//         ErrPred += Err*Err;
//         Pred0 += Signal(i)*Signal(i);
//      }
//      ErrPred /= (double) NbrEvaluate;
//      Pred0 /= (double) NbrEvaluate;
//      ErrPred = sqrt(ErrPred);
//      Pred0 = sqrt(Pred0);
//      if (Verbose == True) cout << "Pred Error for a 0 value prediction (i.e. max err.) = " << Pred0 << endl;
//    }

   return ErrPred;
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    fltarray Data,TrueData;
    char Cmd[256];
    int i,j;
    float Err=0.;
    
    // extern softinfo Soft;

    // lm_check(LIC_M1D);
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    // lm_check(LIC_MR2);
    filtinit(argc, argv);
    
    io_1d_read_data(Name_Imag_In, Data);
    double MeanData = Data.mean();
 
    if (DebugPred == True)
    { 
       cout << "Min(Data) = " << Data.min() << endl;
       cout << "Max(Data) = " << Data.max() << endl;
       cout << "Sigma(Data) = " << Data.sigma() << endl;
       cout << "Mean(Data) = " << MeanData  << endl;
    }
    if (NbrEval == 0)  NbrEval = Data.nx() / 2;
    int Nx = Data.nx();   
    reform_to_1d(Data);
    // FitsHeader.origin = Cmd;
    if ((LastNPix > 0) && (LastNPix < Nx))
    {
       fltarray Buff(LastNPix);
       for (i=0; i < Buff.nx(); i++) Buff(i) = Data(Nx-LastNPix+i);
       Data.reform(LastNPix);
       Data = Buff;
       Nx = LastNPix;
    }
    
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       if (NbrScale > 0) cout << "Number of scales = " << NbrScale << endl;
       cout << "Prediction method = " << StringExMethod(ExtrapolType) << endl;
       cout << "Number of values used for the prediction error evaluation = " << NbrEval << endl;
       if (Noise_Ima > 0) cout << "Noise = " << Noise_Ima << endl;
       if (NoiseAR > 0) cout << "NoiseAR = " << NoiseAR << endl;
       cout << "Regulariazation parameter = " << RegulWave << endl;
       
       if (UseTend == True)
           cout << "   Estimate the tendancy:  NbrPixTendancy = " << NbrPixTendancy <<  endl;
       if ((Stationary == True) && (UseTend == False))
           cout << "Stationary signal " << endl;
       if ((Stationary == True) && (UseTend == True))
           cout << "Stationary signal + tendancy" << endl; 
       if (Stationary == False) 
       {   
           cout << "Non stationary signal " << endl;
           cout << " Extrapolation type = " << StringPredBord(TPredBord) << endl;
	   if ((TPredBord == P_POL1) || (TPredBord == P_POL2))
	    cout << " Nbr of pixels used for the pol. extrapol. = " << NbrPixLastPred<< endl;
       }
       if (EvaluationMode == True)
            cout << "Evaluation mode" << endl;
       cout << endl;
    }
 
    MR_1D MR_Data (Nx, Transform, "MR_Data", NbrScale);
    if (Transform !=TO1_PAVE_HAAR) MR_Data.Border = I_MIRROR;
    else MR_Data.Border = I_CONT;
    if (NPredict > 0) EvaluationMode = False;
   fltarray TabSignalPred(Nx+NPredict);
   fltarray TabSignalFil(Nx+NPredict);
   
   MR_FIL_PRED MRP;
   MRP.AR_AutoDetect =  (PredARNbr < 0) ? True: False;
   MRP.NoLag = NoLag;
   MRP.MaxIter = Max_Iter;
   MRP.Stationary = Stationary;
   MRP.Verbose = Verbose;
   MRP.LastTraining = Nx - NbrEval - 1; 
   MRP.FirstEvaluate = Nx - NbrEval;
   MRP.LastEvaluate = Nx -1;
   MRP.MaxNbrAR = MaxNbrAR;
   MRP.EvaluationMode = EvaluationMode;
   MRP.NbrPixTendancy = NbrPixTendancy;
   MRP.NbrPixLastPred = NbrPixLastPred;
   MRP.MinAR_OrderPerScale = MinAR_OrderPerScale;
   MRP.TPredBord=TPredBord;
   MRP.UseB3SplineWT=UseB3;
   MRP.UseTend=UseTend;
   MRP.Debug = DebugPred;
   MRP.UserSigmaNoise=Noise_Ima;
   MRP.UserSigmaARNoise=NoiseAR;
   MRP.AR_Order_Detect = AR_Order_Detect;
   MRP.RegulWave = RegulWave;
   MR1DNoiseModel NoiseModel(Stat_Noise, Nx,  NbrScale, Transform);
   NoiseModel.SigmaDetectionMethod = SIGMA_MEDIAN_1D;
   
   if (Noise_Ima > FLOAT_EPSILON) NoiseModel.SigmaNoise = Noise_Ima;
   if (UseNSigma  == True)
        for (i=0; i < NbrScale; i++) NoiseModel.NSigma[i]=NSigma;
   // NoiseModel.CCD_Gain = CCD_Gain;
   // NoiseModel.CCD_ReadOutSigma = CCD_ReadOutSigma;
   // NoiseModel.CCD_ReadOutMean = CCD_ReadOutMean;
   NoiseModel.model(Data, MR_Data);
   // MR_Data.transform(Data);
   // cout << "TRANS OK" << endl;
   if ((Verbose == True) && (Noise_Ima <  FLOAT_EPSILON))
     cout << "Estimated Sigma noise = " <<   NoiseModel.SigmaNoise << endl;
   // NoiseModel.SigmaNoise = 0.;
   MRP.NoiseModel = &NoiseModel;
   // MRP.SigmaNoise = NoiseModel.SigmaNoise;
   //Err = MRP.make_prediction(ExtrapolType, Data, TabSignalPred, 
   //                          MR_Data, NPredict, PredARNbr);
   Err = MRP.pred_filter(Data, TabSignalFil,  TabSignalPred,  MR_Data, 
                          NPredict, PredARNbr);

//    if (EvaluationMode == True)
//    {
//       cout << endl;
//       cout << " Std deviation Prediction Err = " << Err << endl;
//       cout << endl;
//    }
   if (NPredict > 0)
   {
      cout << "Prediction = " << endl;
      for (i=0; i < NPredict; i++)
        cout << "+" << i+1 << " ==> " << TabSignalPred(Nx+i) << endl;
      cout << endl;
   }
   
   // cout << "Write results " << Name_Imag_Out << endl;
   // cout << "TabSignalFil " << TabSignalFil.nx() << " " << TabSignalFil.sigma() <<endl;
  
   io_1d_write_data(Name_Imag_Out, TabSignalFil);
  
   if (WriteInfo == True)
   {
      extern type_1d_format IO_1D_Format;
      fltarray TabRes(NbrScale+2);
      TabRes(0) = Err;
      TabRes(1) = NbrScale;
      if (ExtrapolType != EX_AR)
        for (j=0; j < NbrScale; j++) TabRes(j+2) = MRP.TabAR_OrderPerScale(j);
      else TabRes(2) = MRP.TabAR_OrderPerScale(0);
      if (IO_1D_Format == F1D_FITS) fits_write_fltarr(InfoFile, TabRes);
      else io_write2d_ascii(InfoFile, TabRes);
   }
    
   if (Orig == True) 
   {
      if (Verbose == True) cout << " Evaluation " << endl;
      TabSignalFil.info();
       int NbrEvaluate = MRP.LastEvaluate-MRP.FirstEvaluate+1;
      fltarray TabOrig;
      io_1d_read_data(Name_Orig, TabOrig);
      // Error evaluation
      double ErrPred=0.,Err;
      if (Verbose == True) 
        cout << "FirstEvaluate = " << MRP.FirstEvaluate << " Last   = " << MRP.LastEvaluate << endl;
      for (i=MRP.FirstEvaluate; i <= MRP.LastEvaluate; i++)
      {
        Err =  TabOrig(i) - TabSignalFil(i);
        ErrPred += Err*Err;
      }
      ErrPred /= (double) NbrEvaluate;
      ErrPred = sqrt(ErrPred);
      Data -= TabSignalFil;
      cout << "PRED RES: Nbr points for eval = " << NbrEvaluate << " ErrorPred   = " << ErrPred <<   "  Noise Estimation    = " << Data.sigma() <<  endl;

   } 
   exit(0);
} 

