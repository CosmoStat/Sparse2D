/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  mr1d_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  spectrum filtering
**    ----------- 
**                 
**    Usage: mr1d_filter options image output
**
**   where options =
**
**         [-t type_of_multiresolution_transform]
**              1: linear wavelet transform: a trous algorithm
**              2: b1spline wavelet transform: a trous algorithm
**              3: b3spline wavelet transform: a trous algorithm
**              4: Morlet's wavelet transform
**              5: morphological median transform
**              6: mexican hat wavelet transform
**              7: french hat wavelet transform
**              8: pyramidal linear wavelet transform
**              9: pyramidal b3spline wavelet transform
**              10: pyramidal median transform
**         [-m type_of_noise]
**              1: NOISE_GAUSSIAN
**              2: NOISE_POISSON
**              .....
**         [-g sigma]
**             sigma = noise standard deviation
**         [-c gain sigma mean]
**             gain = gain of the CCD
**             sigma = read-out noise standard deviation
**             mean = read-out noise mean
**         [-n number_of_scales]
**             number of scales used in the multiresolution transform
**         [-s nsigma]
**              Thresholding at nsigma * SigmaNoise
**         [-i number_of_iterations]
**              Maximum number of iterations
**         [-e epsilon]
**              Convergence parameter
**         [-v]
**				Verbose 
**
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "MR1D_Filter.h"
#include "MR1D_NoiseModel.h"
 #include "MR_Threshold.h"


extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

char Name_Input[80];      /* input file image */
char Name_Output[80];     /* output file image */
char Name_MR[80];         /* filter Mr file name */
Bool WriteFilterMR = False;

int Nbr_Plan=0;           /* number of scales */

Bool UseNSigma=False;     /* use read sigma */

float NSigma=DEF_NSIG;   /* number of sigma (for the noise) */;
float SigmaNoise=0.;      
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_trans_1d Transform=TO1_PAVE_B3SPLINE;    /* type of transform */

float Epsilon= 1e-4;
int Max_Iter=DEF_ITER;
Bool SetPlan = False;

float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean 

Bool Verbose = False;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_1d_filter Filter = FIL_1D_THRESHOLD;  /* default filtering method */
float EpsilonPoisson = DEF_EPS_EVENT_1D_FILTERING;
Bool EpsOpt=False;
Bool MaxIterOpt=False;
Bool KillLastScale  = False;
char Name_RMSMap[256]; 
Bool UseRMSMap=False;
 
Bool SupIsol = False;  /* suppress isolated pixel in the support */
Bool KeepPositivSup = False;
int FirstScale=0;
int NiterClip = 1;
int SizeBlock = 25;
float RegulParam = 0.01;
Bool PositivIma = True;
type_threshold TypeThreshold = DEF_THRESHOLD;
Bool MultiVect = False;
int MedianWinSize = DEFAULT_MEDIAN_1D_WINSIZE;
extern int  DEF_MedianWinSize;  // see in MR1D_Obj.cc

/*********************************************************************/

int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.));
    return (ScaleMax);
}

/*********************************************************************/

static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options Signal1D_or_Image   Signal1D_or_Image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    filter_1d_usage(Filter);
    manline();
     
    get_threshold_usage(TypeThreshold);
    manline(); 
    
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_DIADIC_TRANS_1D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf1D((type_trans_1d)i));
    fprintf(OUTMAN, "             Default is %s \n", StringTransf1D(Transform));
    manline();

    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n", StringSBFilter ((type_sb_filter)  SB_Filter));
    manline();

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             Default is Gaussian noise.\n");    
    manline();

    gauss_usage();
    manline();  
     
    ccd_usage();
    manline();  
  
    prec_eps_poisson_usage(EpsilonPoisson);
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    manline();

    nsigma_usage(NSigma);
    fprintf(OUTMAN, "             Default is 2 for FDR detection method.\n"); 
    manline();
     
    max_iter_usage(Max_Iter);
    manline();
  
    converg_param_usage(Epsilon);
    manline();

    kill_last_scale_usage();
    manline();    
     
    detect_pos_usage();
    manline(); 
    
    kill_isol_pix_usage();
    manline();

    size_block_usage(SizeBlock);
    manline();
    
    sigma_clip_block_usage(NiterClip);
    manline();
 
    first_detect_scale_usage();
    manline();
    
    rms_noise_usage();
    manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
    manline();
    
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
    manline();
    fprintf(OUTMAN, "         [-W Median_Window_Size]\n");
    fprintf(OUTMAN, "             Window size in the pyramidal median transform. Default is %d.\n", DEFAULT_MEDIAN_1D_WINSIZE);
    manline();        
    verbose_usage();
    manline();
      
    manline();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    Bool TransfOpt=False;
    Bool FilterOpt=False;
    Bool NscaleOpt=False;
    Bool Optf=False;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"W:MR:C:PG:S:N:F:KkpE:T:f:t:m:g:c:n:s:i:e:vw")) != -1) {
	switch (c) {
	    case 'M': MultiVect = True; break;
            case 'P':
                PositivIma= (PositivIma == True) ? False: True;
                break;            
            case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.0;
		break;
        case 'W':
            /* -s <nsigma> */
            if (sscanf(OptArg,"%d",&MedianWinSize) != 1) 
            {
                fprintf(stderr, "bad Window Parameter: %s\n", OptArg);
                exit(-1);
            }
            if (MedianWinSize  < 0.) MedianWinSize = DEFAULT_MEDIAN_1D_WINSIZE;
            break;
            
           case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
		    exit(-1);
		}
		if (SizeBlock  < 2)
		{
		   fprintf(OUTMAN, "Error: bad  SizeBlock parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "              SizeBlock > 1\n");
		   exit(-1);
		}
		break; 
         case 'N':
		if (sscanf(OptArg,"%d",&NiterClip) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad it. number for the 3sigma clipping: %s\n", OptArg);
		    exit(-1);
		}
		if (NiterClip < 1)
		{
		   fprintf(OUTMAN, "Error: bad NiterClip parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "             NiterClip > 0\n");
		   exit(-1);
		}
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
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
	   case 'f':
		/* -f <type> type of filtering */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_1D_FILTERING)) Filter = (type_1d_filter) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of filtering: %s\n", OptArg);
	            exit(-1);
 		}
                FilterOpt = True;
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
 		break;
	    case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) {
		    fprintf (OUTMAN, 
                     "bad type of multiresolution transform: %s\n", OptArg);
			usage(argv);            
		}
		if ((c > 0) && (c <= NBR_DIADIC_TRANS_1D)) 
			Transform = (type_trans_1d) (c-1);
		else {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
			usage(argv);
 		}
                TransfOpt=True;
		break;
	case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
			exit(-1);
		}
		if ((c > 0) && (c <= NBR_NOISE)) Stat_Noise = (type_noise) (c-1);
		else {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
			exit(-1);
 		}
		break;
	case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&SigmaNoise) != 1) {
		    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
		    usage(argv);
		}
		Stat_Noise = NOISE_GAUSSIAN;
		break;
	case 'c':
		/* -c <gain sigma mean> */
		printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &CCD_Gain,
                   &CCD_ReadOutSigma, &CCD_ReadOutMean) <= 0) {
			fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
		Stat_Noise = NOISE_POISSON;
		SigmaNoise = 1.;
		break;
	case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
		if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE_1D)) {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    usage(argv);
		}
		NscaleOpt = SetPlan = True;
 		break;
	case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) {
		    fprintf(OUTMAN, "bad NSigma: %s\n", OptArg);
		    usage(argv);
		}
		if (NSigma <= 0.) NSigma = DEF_NSIG;
		UseNSigma = True;
		break;
	case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) {
		   fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
			usage(argv);
		}
		if (Max_Iter <= 0)  Max_Iter = DEF_ITER;
		MaxIterOpt=True;
                break;
	 case 'E':
		if (sscanf(OptArg,"%f",&EpsilonPoisson) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    exit(-1);
		}
 		if ((EpsilonPoisson <= 0) || (EpsilonPoisson >= 1)) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    fprintf(OUTMAN, "       0 < Precision < 1\n");
 		    exit(-1);
		}		
		break;
	case 'e':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) {
		   fprintf(OUTMAN, "bad convergence parameter: %s\n", OptArg);
		   usage(argv);
		}
		//  if ((Epsilon < 0) || (Epsilon > 1.))  Epsilon = DEF_EPS;
		EpsOpt=True;
		break;
        case 'p':
                KeepPositivSup=True;
                break;           
        case 'v':
		Verbose = True;
		break;
	case 'w':
	        WriteFilterMR = True;
		break;
         case 'R':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseRMSMap = True;
 		break;
	case 'k':
               SupIsol = True;break; 
        case 'K':
               KillLastScale = True;break;  
	case '?':
		usage(argv);
 	}
    }

	// Test inconsistencies in the option call 
 	if (Stat_Noise == NOISE_EVENT_POISSON)
 	{
 	   if ((TransfOpt == True) && (Transform != TO1_PAVE_B3SPLINE))
 	   {
 	      cerr << "WARNING: with this noise model, only the BSPLINE A TROUS can be used ... " << endl;
 
 	      cerr << "        Type transform is set to: BSPLINE A TROUS ALGORITHM " << endl;
 	   }
 	   Transform = TO1_PAVE_B3SPLINE;
  	   if (EpsOpt == False) Epsilon = DEF_EPS_EVENT_1D_FILTERING;
 	}

       // Eps option and MaxIter option are valid only with iterative methods
       if ((EpsOpt == True) || (MaxIterOpt == True))
       {
          if ((Filter != FIL_1D_ITER_THRESHOLD) && (Filter != FIL_1D_TV_CONSTRAINT) )
          {
             fprintf(OUTMAN, "Error: option -e and -i are not valid with non iterative filtering methods. ...\n");
	     fprintf(OUTMAN, "Filter method is %s\n", StringSBFilter ((type_sb_filter)  SB_Filter));
             exit(-1);
	  }
       }

       
       if ((UseNSigma == False) && (TypeThreshold == T_FDR)) 
       {
          NSigma = 2.;
	  UseNSigma = True;
       }
       
       if ((EpsOpt == False) && (UseNSigma == True))
          EpsilonPoisson = (1. - erf((double) NSigma / sqrt((double) 2.)));

 	if ((Transform != TU1_MALLAT) && (Transform != TO1_MALLAT) && (Optf == True))
	{
	   fprintf(OUTMAN, "Error: option -T is only valid with Mallat transform ... \n");
           exit(0);
	}
	     

	/* get optional input file names from trailing parameters and open files */
	if (OptInd < argc) strcpy(Name_Input, argv[OptInd++]);
	else usage(argv);

	if (OptInd < argc) strcpy(Name_Output, argv[OptInd++]);
	else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}

/*********************************************************************/

int main(int argc, char *argv[])
{
  int Nx,Ny=1,i;
  fltarray Data;
  fltarray MData,MResult;
  fitsstruct FitsHeader;

  // Get command line arguments, open input file(s) if necessary 
  lm_check(LIC_M1D);
  filtinit(argc, argv);
  char Cmd[256];
  Cmd[0] = '\0';
  for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
        
  // read the data
  if (MultiVect == True) 
  {
      fits_read_fltarr(Name_Input, MData, &FitsHeader);
      Nx = MData.nx();
      Ny = MData.ny();
      MResult.alloc(Nx,Ny);
      Data.alloc(Nx);
  }
  else  
  {
     io_1d_read_data(Name_Input, Data, &FitsHeader);
     reform_to_1d(Data);
     Nx = Data.n_elem();
  }
  FitsHeader.origin = Cmd;
   
  // out result
  fltarray Result (Nx);
  fltarray RMS;   
  Bool MultiRMS = False;

  // get nb Plan
  if (!SetPlan) Nbr_Plan = max_scale_number(Nx);

  // trace
  if (Verbose) 
  {
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Input << endl;
    cout << "File Name Out = " << Name_Output << endl;
    cout << "Transform = " << StringTransf1D(Transform) << endl;
    if ((Transform == TO1_MALLAT) || (Transform == TU1_MALLAT))
        cout << StringSBFilter(SB_Filter) << endl;
    cout << "Filtering Method = " << String1DFilter(Filter) << endl;
    cout << "Coef Detect Method = " <<  StringThreshold(TypeThreshold) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "Type of Noise = " << StringNoise(Stat_Noise) << endl;
    cout << "Sigma Noise = " << SigmaNoise << endl;
    cout << "NSigma = " << NSigma << endl;
    cout << "Epsilon = " << Epsilon << endl;
    cout << "Epsilon poisson = " << EpsilonPoisson << endl;  
    if (KeepPositivSup == True) 
       cout << "Only positive wavelet coefficients are detected" << endl;
    cout << "naxis = " << Data.naxis() << endl;
    cout << "Nx  = " << Nx << endl;
    if (MultiVect == True) cout << "Ny  = " << Ny << endl;
    cout << "Min = " << Data.min() << endl;
    cout << "Max = " << Data.max() << endl;
    cout << "sigma = " << Data.sigma() << endl;
    if (MedianWinSize != DEFAULT_MEDIAN_1D_WINSIZE) cout << "MedianWindowSize = " << MedianWinSize << endl;
   }

   DEF_MedianWinSize = MedianWinSize;
//    mr1d_filter (Data, Result, SigmaNoise, Transform,  Nbr_Plan,
//                  Stat_Noise,  NSigma, Epsilon,  Max_Iter);
//    if (Stat_Noise == NOISE_GAUSSIAN)
//       cout << "Noise Standard deviation = " << SigmaNoise << endl;


  // noise model class initialization
  FilterAnaSynt FAS;
  FilterAnaSynt *PtrFAS = NULL;
  if ((Transform == TO1_MALLAT) || (Transform == TU1_MALLAT))
  {
      FAS.Verbose = Verbose;
      FAS.alloc(SB_Filter);
      PtrFAS = &FAS;
  }
 
    if (UseRMSMap == True) 
    {
        io_1d_read_data(Name_RMSMap, RMS);
        if ((MultiVect == True) && (RMS.naxis() == 2))
        {
            if ( (MData.nx() != RMS.nx()) || (MData.ny() != RMS.ny()))
            {
                fprintf(OUTMAN, "Error: DATA and RMS inputs must have the same dimensions ...\n");
                fprintf(OUTMAN, "Data dim = %d, %d \n", MData.nx(),MData.ny());
                fprintf(OUTMAN, "RMS dim = %d, %d \n", RMS.nx(),RMS.ny());
                exit(-1);
            }
            MultiRMS = True;
        }
    }   


    
  for (int v=0; v < Ny; v++)
  {
     if (MultiVect == True)
     {
        for (i=0; i < Nx; i++)  Data(i) = MData(i,v);
	    Result.init();
     }
     MR1DNoiseModel NoiseModel;
     NoiseModel.alloc(Stat_Noise, Data.nx(), Nbr_Plan, Transform, PtrFAS);
     if (SigmaNoise > FLOAT_EPSILON) NoiseModel.SigmaNoise = SigmaNoise;
     if (UseNSigma  == True)
          for (i=0; i < Nbr_Plan; i++) NoiseModel.NSigma[i]=NSigma;
     for (i=0; i < Nbr_Plan; i++) NoiseModel.TabEps[i] = EpsilonPoisson;
     if (SupIsol == True) NoiseModel.SupIsol = SupIsol;
     NoiseModel.OnlyPositivDetect=KeepPositivSup;
     NoiseModel.FirstDectectScale = FirstScale;
     NoiseModel.SizeBlockSigmaNoise = SizeBlock;
     NoiseModel.NiterSigmaClip = NiterClip;
     NoiseModel.TypeThreshold = TypeThreshold;
     if (UseRMSMap == True) 
     {
       NoiseModel.UseRmsMap = True;
       NoiseModel.RmsMap.alloc(Nx);
       if ((MultiVect == True) && (MultiRMS ==  True))  for (i=0; i < Nx; i++)   NoiseModel.RmsMap(i) = RMS(i,v);
       else  for (i=0; i < Nx; i++)  NoiseModel.RmsMap(i) = RMS(i);
     }
     
     // Filtering class initialization
     MR1DFiltering CFilter(NoiseModel, Filter);
     if (MaxIterOpt == True)  CFilter.Max_Iter = Max_Iter;
     if (KillLastScale == True)
     { 
         CFilter.KillLastScale = KillLastScale;
         CFilter.PositivIma = False;
     }
     CFilter.WriteFilterMR = WriteFilterMR;
     if (WriteFilterMR == True) CFilter.NameSupport = strdup("FiterMr1d.mr");

     if (EpsOpt == True) CFilter.Epsilon = Epsilon;
     CFilter.PositivIma = PositivIma;
     // if (MaxIma == True) CFilter.MaxIma = True;
     // if (WindowOpt == True) CFilter.SmoothWindowSize = SmoothWindowSize;
     if (v==0) CFilter.Verbose = Verbose;
     if (Filter == FIL_1D_TV_CONSTRAINT) 
                   CFilter.RegulParam = RegulParam;
     CFilter.filter(Data, Result);
     
     if (MultiVect == True)  for (i=0; i < Nx; i++) MResult(i,v) = Result(i);
     CFilter.Verbose = False;
   }
   if (MultiVect == True) fits_write_fltarr(Name_Output, MResult, &FitsHeader);
   else io_1d_write_data(Name_Output, Result, &FitsHeader);
   exit(0);
} 

