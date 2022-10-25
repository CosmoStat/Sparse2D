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

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];
char InfoFile[256];

int NbrScale=5;           // number of scales
float NSigma=3.;   // number of sigma (for the noise) 
float SigmaNoise=0.;       // noise standard deviation
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

Bool WriteInfo = False;

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

    extrap_usage(ExtrapolType);
    manline();
  
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    fprintf(OUTMAN, "             Default is 5.\n");
    manline();

    fprintf(OUTMAN, "         [-a AR_Order]\n");
    fprintf(OUTMAN, "             AR oder used for the prediction. Default is automatically estimated.\n");
    manline();

//     fprintf(OUTMAN, "         [-l]\n");
//     fprintf(OUTMAN, "             Do not use lags between wavelet coefficients. Default is no.\n");
//     fprintf(OUTMAN, "             Only used with the global AR method. Default is no.\n");
//     manline();

    fprintf(OUTMAN, "         [-O Estimation_AR_Order_Method]\n");
    fprintf(OUTMAN, "             1: AIC \n");
    fprintf(OUTMAN, "             2: AICC\n");
    fprintf(OUTMAN, "             3: BIC\n");
    fprintf(OUTMAN, "            Default is 3.\n");
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

    fprintf(OUTMAN, "         [-w]\n");
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
    
#if TEST_B3SPLINE
    fprintf(OUTMAN, "         [-b]\n");
    fprintf(OUTMAN, "             Use a B3-spline WT instead of the Haar transform. Default is no.\n");
    manline();
#endif

#ifdef TEST_MODE
    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "             Simulation mode. An AR 4 is simulated for the input data.\n");
    fprintf(OUTMAN, "         [-N Number_of_pix]\n");
    fprintf(OUTMAN, "             Number of pixels in the simulation mode. Default is 512.\n");
    fprintf(OUTMAN, "         [-I RndVal]\n");
    fprintf(OUTMAN, "            Initial value for the random value generator. Only used in simulation mode.\n");
    fprintf(OUTMAN, "         [-W]\n");
    fprintf(OUTMAN, "             Replace the simulated AR 4 by a random noise. Only used in simulation mode.\n");
#endif

    manline();
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
    while ((c = GetOpt(argc,argv,"wL:tB:T:xO:M:hlP:a:p:m:n:v")) != -1) 
    {
	switch (c) 
        {
	   case 'w': 
// 	             if (sscanf(OptArg,"%s",InfoFile) != 1) 
//                      {
// 		        fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
// 	                exit(-1);
// 		     }
	             WriteInfo = True; 
	             break;
	   case 'x': DebugPred = True; break;
	   case 't': UseTend = True;break;
	   case 'b': Transform=TO1_PAVE_B3SPLINE;
	             UseB3 = True;
		     break; 
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
           case 'W': OnlyNoise = True; break;
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
	   case 'P':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of extrapolation: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_EX_METHOD)) 
                                        ExtrapolType  = (extrapol_type) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of extrapolation: %s\n", OptArg);
	            exit(-1);
 		}
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

        if ((UseB3 == True) && (ExtrapolType != EX_MR_AR_PRED2))
	{
	   cout << "Error: -b option is only valid when '-P 3' option is set ..."<<endl;
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
		usage(argv);
	}
}

/*********************************************************************/

void estime_param(fltarray & Data, float &SigmaNoiseData, float &SigmaNoise)
{
  int i,Np = Data.nx();
  fltarray Diff(Np-1);
  double E0, E1,E2;
  E0=E1=E2=0.;
  for (i=1; i < Np; i++) Diff(i-1) = Data(i) - Data(i-1);
  for (i=0; i < Np-1; i++)
  {
     E0 += Diff(i)*Diff(i);
     if (i > 0) E1 += Diff(i)*Diff(i-1);
     if (i > 1) E2 += Diff(i)*Diff(i-2);
  }
  E0 /= (double) (Np-1);
  E1 /= (double) (Np-2);
  E2 /= (double) (Np-3);

  cout << "E(0) = " << E0 << " E(1) = " << E1 << " E(2) = " << E2 << endl;
  SigmaNoiseData = -E1;
  SigmaNoise = E0 + 2*E1;

  cout << "Variance (noise_data) = " << SigmaNoiseData << endl;
  cout << "Variance (noise_signal) = " <<  SigmaNoise << endl;
  if (SigmaNoiseData > 0) SigmaNoiseData = sqrt(SigmaNoiseData);
  else SigmaNoiseData =0;
  if (SigmaNoise  > 0)  SigmaNoise = sqrt(SigmaNoise );
  else  SigmaNoise =0;

  cout << "Sigma (noise_data) = " << SigmaNoiseData << endl;
  cout << "Sigma (noise_signal) = " <<  SigmaNoise << endl;
  cout << "Est. AR(1) = " << E1 / E0 << endl;
}

/*********************************************************************/

float mksimu(fltarray & Data, fltarray & TrueData)
{
   int i,k;
   int Np = NpSimu;
   Data.alloc(Np);
   TrueData.alloc(Np);
   int Nar = NBR_ARCoefSimu;
   double x=0.;
   double Sigma = 1.;
   double SigmaNoise = 0.;

   init_random(InitRnd);
   cout << "Simulation: Np = " << Np << " AR(" << Nar << ")" << endl;
   cout << "AR = ";
   for (k=0; k < Nar; k++) cout << TabARSimu[k] << " ";
   cout << endl;


   Data.init();
   for (i=0; i < Nar;i++) TrueData(i) = Data(i) = (float) i/2.;
 
   for (i=Nar; i < Np;i++) 
   {
      for (k=0; k < Nar; k++) 
      {
         Data(i) += TabARSimu[k] *  Data(i-1-k);
      }
      x =   get_random();
      Data(i) +=   xerf(x) * Sigma;
      TrueData(i) = Data(i);
   }

   if (SigmaNoise > 0) for (i=0; i < Np;i++)  Data(i) += xerf(x) * SigmaNoise;

   if (OnlyNoise == True)
   {
      Data.init();
      init_random(InitRnd);
      for (i=0; i < Np;i++)  
      {
        x = get_random();
        TrueData(i) = Data(i) = xerf(x) * Sigma; //  * 0.01;
      } 
      for (i=Np-10; i < Np;i++) cout << " " << Data(i);
      cout << endl;
   }

   i = Np;
   float Sol=0.;
   for (k=0; k < Nar; k++) Sol +=  TabARSimu[k] * TrueData(i-1-k);
   if (Verbose == True) cout << "Ex. Sol = " << Sol << endl;

   fits_write_fltarr("xx.fits", Data);

   float ESigmaNoiseData,ESigmaNoise;
    if (Verbose == True) estime_param(Data, ESigmaNoiseData,ESigmaNoise);
   return Sol;
}

/*********************************************************************/
 
// static void write2d_ascii(char *Name_Dat_Out,  fltarray & Dat)
// {
//    FILE *FileOut=NULL;
//    int Np=0;
//    int Nl;	
//    int j,i=0;
// 
//    /* Read the size of the FileOut image */
//    if (FileOut != stdout) 
//    {
//       FileOut = fopen( Name_Dat_Out,"w");
//       if (FileOut == NULL) 
//       {
//         cout << "Error: cannot open file " <<  Name_Dat_Out  << endl;
//         exit(-1);
//       }
//    }
//    Np = Dat.nx();
//    Nl = Dat.ny();
//    // cout << "Np = " << Np << endl;
//     for (i=0; i < Np; i++)
//     {
//        for (j=0; j < Nl; j++) fprintf(FileOut,"%f  ", Dat(i,j));
//        fprintf(FileOut,"\n");
//     }
//     if (FileOut != stdout) fclose(FileOut);
// }
//  
/*********************************************************************/

void mk_all_tendancy(fltarray & Data, fltarray & Tend)
{
   int k;
   int Np = Data.nx();
   fltarray Tab(Np);
   cout << "mk_all_tendancy " << Np << " " << Tend.nx() << endl;
   
   tendancy_est(Data, Tend, NbrPixTendancy, Np);
   for (k=800; k<Np; k++) 
   {
      cout << " k = " << k << endl;
      tendancy_est(Data,  Tab, NbrPixTendancy, k);
      Tend(k-1) = Tab(k-1);
   }
   cout << "end " << endl;
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    fltarray Data,TrueData;
    char Cmd[256];
    int i,j;
    float TrueSol=0.;
    float Err=0.;  
    extern softinfo Soft;

#ifdef USELM
    lm_check(LIC_M1D);
#endif

    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    filtinit(argc, argv);
    
    // fits_read_fltarr(Name_Imag_In, Data);
    if (Simu == False) io_1d_read_data(Name_Imag_In, Data);
    else TrueSol = mksimu(Data,TrueData);
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
       cout << "Number of values used for the prediction error evaluation = " << NbrEval << endl << endl;
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
       if (UseB3 == True)
            cout << "Use the B3-spline WT instead of the Haar WT" << endl;   
       if (EvaluationMode == True)
            cout << "Evaluation mode" << endl;
    }
 
   if (MyTest == True)
   {      
     fltarray Pred(Nx);
     int N = (PredARNbr < 0) ? 128: PredARNbr;
     NN_PREDICTION NNP;
     NNP.TrainLWB = N;
     NNP.TrainUPB = Nx/2-1;
     NNP.TestLWB =  Nx/3;
     NNP.TestUPB = Nx/2;
     NNP.EvalLWB = Nx/2;
     NNP.EvalUPB = Nx - 1;
     NNP.alloc(N, Nx);
     // NNP.Verbose = True;
     double NNErr = NNP.eval_predict(Data.buffer(), Pred.buffer());  
     cout << "Sigma Err " << NNErr  << endl;;  
     exit(0);
   }
    

    MR_1D MR_Data (Nx, Transform, "MR_Data", NbrScale);
    if (Transform !=TO1_PAVE_HAAR) MR_Data.Border = I_MIRROR;
    else MR_Data.Border = I_CONT;
    
//     fltarray TT(Data.nx());
//     mk_all_tendancy(Data, TT);
//     fits_write_fltarr("xx_tt.fits", TT);
//     exit(-1);

//     if (ExtrapolType != EX_AR) 
//     {
//        MR_Data.transform(Data);
//        // wave_1d_trans(Data, MR_Data, Nx);
//        fits_write_fltarr("xx_trans.fits", MR_Data.image()); 
//     }

// 
//    // Filtering the wavelet coefficients
//    if (ExtrapolType != EX_AR)
//    {
//      if ((HardThreshold == True) || (EntrThreshold == True))
//      {
//         MR1DNoiseModel NoiseModel(Stat_Noise, Nx, NbrScale, Transform);
//         if (SigmaNoise > FLOAT_EPSILON) NoiseModel.SigmaNoise = SigmaNoise;
//         if (UseNSigma  == True)
//           for (i=0; i < NbrScale; i++) NoiseModel.NSigma[i]=NSigma;
//         NoiseModel.CCD_Gain = CCD_Gain;
//         NoiseModel.CCD_ReadOutSigma = CCD_ReadOutSigma;
//         NoiseModel.CCD_ReadOutMean = CCD_ReadOutMean;
//         // NoiseModel.SizeBlockSigmaNoise = SIZE_BLK_NOISE_ADD_NON_UN;
//   
//         NoiseModel.model(Data, MR_Data);
//         if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN))
//         cout << "Noise standard deviation = " << NoiseModel.SigmaNoise << endl;
// 
//         if (HardThreshold == True) NoiseModel.threshold(MR_Data);
//         else
//         {
//            fltarray TabAlpha(NbrScale);
//            float RegulParam=1.;
//            Bool DataSnr = True;
//            float Epsilon= DEF_EPS_CVG_MEMW1D_FILTER;   
//            int Max_Iter=DEF_MEMW1D_FILTER_MAX_ITER;
//            MR_1D  *MR_Model = NULL;
//            for (b=0; b < NbrScale; b++) TabAlpha(b) = 0.;
//            mw1d_filter (MR_Data, NoiseModel, TabAlpha,  RegulParam, 
//                    DataSnr, Epsilon,  Max_Iter, Verbose, MR_Model);  
//         }
//      }  
//      else  MR_Data.transform(Data);
//   }
   if (NPredict > 0) EvaluationMode = False;
   fltarray TabSignalPred(Nx+NPredict);
   MR_PRED MRP;
   MRP.AR_AutoDetect =  (PredARNbr < 0) ? True: False;
   MRP.NoLag = NoLag;
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
   MRP.AR_Order_Detect = AR_Order_Detect;
   float Inter;
   Err = MRP.make_prediction(ExtrapolType, Data, TabSignalPred, 
                             MR_Data, NPredict, PredARNbr, Inter);

   if (UseTend == True)
   {
     // fits_write_fltarr("xx_tend.fits", Tend);
      if (MRP.Tend.n_elem() > 0) 
             io_1d_write_data("xx_tend", MRP.Tend);
   }


   if (EvaluationMode == True)
   {
      cout << endl;
      cout << " Std deviation Prediction Err = " << Err << endl;
      cout << endl;
   }
   if (NPredict > 0)
   {
      cout << "Prediction = " << endl;
      for (i=0; i < NPredict; i++)
        cout << "+" << i+1 << " ==> " << TabSignalPred(Nx+i) << endl;
      cout << endl;
   }
   io_1d_write_data(Name_Imag_Out, TabSignalPred);
   sprintf(InfoFile, "%s_inf.dat", Name_Imag_Out);
   if (WriteInfo == True)
   {
      extern type_1d_format IO_1D_Format;
     fltarray TabRes(NbrScale+3);
     TabRes(0) = Err;
     TabRes(1) = NbrScale;
     TabRes(2) = Inter;
     if (ExtrapolType != EX_AR)
       for (j=0; j < NbrScale; j++) TabRes(j+3) = MRP.TabAR_OrderPerScale(j);
     else TabRes(3) = MRP.TabAR_OrderPerScale(0);
     
      if (IO_1D_Format == F1D_FITS) fits_write_fltarr(InfoFile, TabRes);
      else io_write2d_ascii(InfoFile, TabRes);
   }
   exit(0);
} 

