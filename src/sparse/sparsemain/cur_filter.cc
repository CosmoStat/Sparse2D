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
**    Date:  03/02/00
**    
**    File:  cur_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet filtering program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Curvelet.h"
#include "Usage.h"
#include "IM_Regul.h"
#include "MR_SoftRegul.h"
#include "CErf.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
float N_Sigma2=0;
Bool Simu=False;
int BlockSize=DEF_CUR_BLOCK_SIZE;
int NbrScale1D = -1;
Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool PositivSol=True;
Bool ColTrans = False;
type_curvelet_block TypeBlock = DEF_TYPE_CUR_BLOCK;

char NoiseFile[256];  /*  list of input image names  */
Bool ReadNoiseIma = False;
Bool KillLastScale = False;
Bool MAD_Estim = False;
int FirstScale = 0;
Bool OptNSigma= False;
// #include"Simu.cc"
Bool NegativOK=True;
Bool Wiener=False;

#define NBR_NOISE_CUR 7

enum type_noise_cur {CUR_NOISE_GAUSSIAN, CUR_NOISE_POISSON, CUR_NOISE_MULT, 
                     CUR_NOISE_MULT_CORREL, CUR_NOISE_CORREL, 
                     CUR_NOISE_LOG_RAYLEIGH, CUR_NOISE_POISSON_FEW_EVENT};

inline char * StringTypeNoiseCur (type_noise_cur type)
{
    switch (type)
    {
     case CUR_NOISE_GAUSSIAN: 
      return ("Gaussian noise");
      break;
     case CUR_NOISE_POISSON: 
      return ("Poisson noise");
      break;
    case CUR_NOISE_MULT: 
      return ("Multiplicative noise");
      break;
    case CUR_NOISE_MULT_CORREL:
      return ("Multiplicative correlated noise");
      break;
     case CUR_NOISE_LOG_RAYLEIGH: 
      return ("Log transformed Rayleigh noise");
      break;
     case CUR_NOISE_CORREL: 
      return ("Stationary correlated noise");
      break;
     case CUR_NOISE_POISSON_FEW_EVENT: 
      return ("Poisson noise with few events");
      break;
    default:
      return ("Unknown noise");
    }
}

type_noise_cur Stat_Noise = CUR_NOISE_GAUSSIAN;
#define DEFAULT_MAX_ITER_CUR_FILTER 10
int Max_Iter = 0;
int NbrImage = 1;
float RegulParam = 0.;
unsigned int InitRnd = 100;
Bool FDR = False;
Bool Interval = False;
int WienerBlockSize=7;

/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

//     fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
//     for (int i = 0; i < NBR_RID_DECIMATED; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
//     fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
//    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    manline();

    fprintf(OUTMAN, "         [-N number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated. \n");    
    manline();
 
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is %d. \n", BlockSize);    
    manline();

//     fprintf(OUTMAN, "         [-B type_of_curvelet_block]\n");
//     for (int i = 0; i < NBR_CUR_TYPE_BLOCK; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i+1,  StringBlockType((type_curvelet_block) (i)));
//     fprintf(OUTMAN, "              Default is %s.\n",   StringBlockType(TypeBlock));
//     manline();

    fprintf(OUTMAN, "         [-m type of noise]\n");
    for (int i = 0; i < NBR_NOISE_CUR; i++)
      fprintf(OUTMAN, "              %d: %s \n",i+1,
	      StringTypeNoiseCur( (type_noise_cur)i ));
    manline();
         
    gauss_usage();
    manline();  
     
    nsigma_usage(N_Sigma);
    fprintf(OUTMAN, "             Default is 2 for FDR detection method.\n"); 
    manline();
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "             Coefficient detection using the FDR method \n"); 
    fprintf(OUTMAN, "             instead of the standard k-sigma approach. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Do not apply block overlapping. \n"); 
    fprintf(OUTMAN, "             By default, block overlapping is used. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Supress the positivity constraint. Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-I NoiseFileName]\n");
    manline();

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-i iteration]\n");
    fprintf(OUTMAN, "              number max of iterations\n");
    fprintf(OUTMAN, "              default is %d\n",DEFAULT_MAX_ITER_CUR_FILTER);
    manline();

    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
    
    fprintf(OUTMAN, "         [-F first_detection_scale]\n");
    fprintf(OUTMAN, "              First scale used for the detection\n");
    fprintf(OUTMAN, "              default is 1.\n");        
    
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
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"pCF:KM:G:i:m:I:B:Pt:Ob:W:s:g:n:N:vzZ")) != -1) 
    {
	switch (c) 
        { 
	 case 'p': NegativOK = False; break; 
	 case 'C': FDR = True; break;
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
	   case 'K': KillLastScale = True; break;
	   case 'M':
		/* -N <Number of image> */
		if (sscanf(OptArg,"%d",&NbrImage) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of image: %s\n", OptArg);
		    exit(-1);
		}
                if (NbrImage < 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of image: %s\n", OptArg);
		    fprintf(OUTMAN, "       0 < NbrImage \n");
 		    exit(-1);
		}
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
	   case 'm':
		/* -t <Stat_Noise> of noise */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "bad noise type: %s\n", OptArg);
	            usage(argv); 
		}
                if ((c > 0) && (c <= NBR_NOISE_CUR))  Stat_Noise = (type_noise_cur) (c-1);
                else  
                {
		    fprintf(stderr, "bad type of noise: %s\n", OptArg);
	            usage(argv);
 		}
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0)   Max_Iter = DEFAULT_MAX_ITER_CUR_FILTER;
		break;
           case 'I':
	        if (sscanf(OptArg,"%s",NoiseFile) != 1) 
                {
                    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
                    exit(-1);
                }
		ReadNoiseIma = True;
		Stat_Noise = CUR_NOISE_CORREL; 
 	        break;
	   // case 'C':  ColTrans =(ColTrans  == True) ? False: True; break; 
           case 'P': PositivSol = (PositivSol == True) ? False: True; break;
           case 't': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID_TRANS)) RidTrans = (type_ridgelet_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
	  case 'B': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of block: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_CUR_TYPE_BLOCK)) TypeBlock = (type_curvelet_block) (c-1);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of block: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
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
          case 'W': 
                /* -s <nsigma> */
                if (sscanf(OptArg,"%d",&WienerBlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad BlockSize: %s\n", OptArg);
                    usage(argv);
                }
		Wiener=True;
                if (WienerBlockSize <= 0.) WienerBlockSize = 7;
                break;
          case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
                    usage(argv);
                }
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		OptNSigma=True;
                break;
           case 'v': Verbose = True;break;
	   case 'n':
                /* -n <NbrScale2D> */
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
                    exit(-1);
                }
                break;
	   case 'N':
                 if (sscanf(OptArg,"%d",&NbrScale1D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
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
       if ((OptNSigma == False) && (FDR == True)) N_Sigma = 2;

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

  
/***************************************************************************/

// void cur_threshold (Curvelet &Cur, Ifloat * &Trans, float N_Sigma, float Noise_Ima, Bool Verbose, Ifloat * &TransNoise)
void cur_threshold (Curvelet &Cur, Ifloat * &Trans, float N_Sigma, float Noise_Ima, Bool Verbose, Bool Interval)
{
    int i,j,s1d,s2d,b;
    // char Name[256];  
    Bool BesselOK = False;
 
    for (int s2d=0; s2d < FirstScale; s2d++) Trans[s2d].init();
       
    for (b=0; b < Cur.nbr_band()-1; b++)
    {        
       float Sig = Cur.sigma_noise(b);
       Cur.get_scale_number(b, s2d, s1d);
       Ridgelet *Rid = Cur.get_ridgelet(s2d);
       // int BSize = Cur.TabBlockSize(s2d);
       Ifloat Band;
       Rid->get_scale(Trans[s2d], Band, s1d);
  
       int B2 = WienerBlockSize / 2;
       
//        if (Verbose == True)
//        {
//          cout << "Band " << s2d+1 << " " << s1d+1 << ": BSize BS = " << BSize << " Nl = " << Band.ny() << " Nc = " << Band.nx() << endl;
// 	 cout << "        SigmaNoise = " << Cur.sigma_noise(b) << " Nsigma = " << Cur.nsigma(b) << endl;
// 	 cout << "        MinLevel = " << Cur.min_level(b)  << " MaxLevel = " << Cur.max_level(b) << endl;
//        }
       if (BesselOK == True)
       {
          if (Verbose == True) cout << "BESSEL MODEL " << endl;
          void bessel_estimate_bayes(float *, int , float, Bool);
	  float *Ptr = Band.buffer();
	  int N = Band.nx()*Band.ny();
	  bessel_estimate_bayes(Ptr, N, Sig, True);
       }
       int Step=1;
       float  VarianceSignal=0.,VarianceData=0., VarianceNoise = Sig*Sig;
       for (i=0; i < Band.nl(); i++)
       for (j=0; j < Band.nc(); j++)
       {
           if (Wiener == True)
	   {
	      VarianceData = 0.;
              for (int k=i-B2*Step; k <= i+B2*Step; k+=Step) 
              for (int l=j-B2*Step; l <= j+B2*Step; l+=Step) VarianceData += Band(k,l,I_MIRROR)*Band(k,l,I_MIRROR);
              VarianceData /= (float)(WienerBlockSize*WienerBlockSize); 
	      VarianceSignal = MAX(0,VarianceData-VarianceNoise);
 	      Band(i,j) *= VarianceSignal  / VarianceData;
	   }
	   else
	   {
             float Coef = Band(i,j);
	     if ((Interval == False) || (ABS(Coef) > Sig))
	     {
                if ((Coef >= 0) && (Coef < Cur.max_level(b))) Band(i,j) = 0.; 
	        else  if ((Coef < 0) && ((NegativOK == False) ||  (Coef > Cur.min_level(b)))) Band(i,j) = 0.;
             }
	     else Band(i,j) = 0.; 
	   }
       }
       // Cur.put_band(Trans,s2d,s1d,Band);
       Rid->put_scale(Trans[s2d], Band, s1d);
    }
    if (KillLastScale == True)  Trans[NbrScale2D-1].init(); 
}

/***************************************************************************/

void cur_support_threshold (Curvelet &Cur, Ifloat * &Trans, Ifloat * &Support)
{
    int i,j,s1d,s2d,b;
   
    for (int s2d=0; s2d < FirstScale; s2d++) Trans[s2d].init();
    for (b=0; b < Cur.nbr_band()-1; b++)
    {        
       Cur.get_scale_number(b, s2d, s1d);
       Ridgelet *Rid = Cur.get_ridgelet(s2d);
       Ifloat Band,BandSupport;
       Rid->get_scale(Trans[s2d], Band, s1d);
       Rid->get_scale(Support[s2d], BandSupport, s1d);

       for (i=0; i < Band.nl(); i++)
       for (j=0; j < Band.nc(); j++) if (BandSupport(i,j) == 0) Band(i,j) = 0;
       Rid->put_scale(Trans[s2d], Band, s1d);
    }
    if (KillLastScale == True)  Trans[NbrScale2D-1].init(); 
}

/***************************************************************************/

int main(int argc, char *argv[])
{
    int i,j,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    // RegulIma RI;
    MR_Regul RI;
    RI.ExpDecreasingLambda = True;

    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
#ifndef MRPOL    
    lm_check(LIC_MR4);
#else
    lm_check(LIC_POL);
#endif

    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
        cout << "Noise type = " << StringTypeNoiseCur(Stat_Noise) << endl;
        cout << "NSigma = " << N_Sigma  << endl;  
	cout << "Type block = " << StringBlockType(TypeBlock) << endl;
	if (FDR == True) cout << "Detection method = FDR" << endl;
	else cout << "Detection method = K-Sigma" << endl;
    }
    RI.NbrScale = NbrScale2D;

    Ifloat Data;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    int Nl = Data.nl();
    int Nc = Data.nc();
    
    Header.origin = Cmd;
//     if ((is_power_of_2(Data.nl()) == False) || (is_power_of_2(Data.nc()) == False))
//     {
//            cout << "Error: image size must be power of two ... " << endl;
//            exit(-1);
//     }
    if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False))
    {
           cout << "Error: Block size must be a power of two ... " << endl;
           exit(-1);
    }
    if (BlockSize > Data.nc())
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }

    if ((Max_Iter == 0) && (RegulParam > 0)) Max_Iter = 10;
    
     // if (Simu == True) make_ima(Data);


    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    SB1D.setBorder(I_MIRROR);
    Curvelet Cur(SB1D);
    Cur.RidTrans = RidTrans;
    Cur.NbrScale2D = NbrScale2D;
    Cur.ColTrans = ColTrans;
    Cur.TypeBlock = TypeBlock;
    Cur.OddBlockSize=False;
    Cur.WPTrans = False;
    if (NbrScale1D >= 0)
    {
        Cur.NbrScale1D = NbrScale1D;
        Cur.GetAutoNbScale = False;
        if (NbrScale1D < 2) Cur.BlockOverlap = False;
        else Cur.BlockOverlap = BlockOverlap;
    }
    else Cur.BlockOverlap = BlockOverlap;
    if (Stat_Noise == CUR_NOISE_POISSON_FEW_EVENT) Cur.VarNorm =  True;
    else Cur.VarNorm = False;
    
    // Cur.Verbose = Verbose;
    if (Stat_Noise == CUR_NOISE_POISSON_FEW_EVENT) 
    {
        Cur.MSVST = True;
	Stat_Noise = CUR_NOISE_GAUSSIAN;
    }
    Cur.alloc(Data.nl(), Data.nc(), BlockSize);
    Ifloat NoiseData;
    switch (Stat_Noise)
    {
       case CUR_NOISE_CORREL:
               if (ReadNoiseIma == True)
	       {
                   io_read_ima_float(NoiseFile, NoiseData);
                   if ((NoiseData.nl() != Data.nl()) || (NoiseData.nc() != Data.nc()))
                   {
                      cout << "Error: the noise map must have the same size" << endl;
                      cout << "       as the input data ... " << endl;
	              exit(-1);
                   }             
                   Cur.set_noise_level(NoiseData, N_Sigma);
		   Noise_Ima = 1.;
              }
	      else 
	      { 
	         // The noise is estimated from the MAD of each band
		 // of the curvelet transform of the data
		 // if we do not iterate, this will be done in cur_threshold.
		 // Otherwise, we have to set the normalization coefficients
		 // in the curvelet transform
	         Noise_Ima = 1.;
	         MAD_Estim = True;
		 Cur.set_noise_level(Data, N_Sigma, True);
              }
	      break;
       case CUR_NOISE_POISSON:  
	        noise_poisson_transform (Data, Data);
                Noise_Ima = 1.;
		// Noise_Ima = detect_noise_from_med (Data);
		if (Verbose == True) 
		           cout << "Sigma Noise = " << Noise_Ima << endl;
		Cur.set_noise_model_gaussian(N_Sigma, 1.);
		break;
       case CUR_NOISE_GAUSSIAN: 
                if  ((Noise_Ima < FLOAT_EPSILON) && (Cur.MSVST == False))
                {
                   Noise_Ima = detect_noise_from_med (Data);
                   if (Verbose == True) 
		           cout << "Sigma Noise = " << Noise_Ima << endl;
                }     		
		if (Cur.MSVST == False) Cur.set_noise_model_gaussian(N_Sigma, Noise_Ima );
		else Cur.set_noise_model_gaussian(N_Sigma, 1.);
		
		// NoiseData.alloc(Data.nl(), Data.nc(), "noise");
                // im_noise_gaussian (NoiseData, Noise_Ima, InitRnd);
 		break;   
      case CUR_NOISE_MULT:
               noise_log_transform (Data, Data);
	       if (Noise_Ima < FLOAT_EPSILON)
               {
                   Noise_Ima = detect_noise_from_med (Data);
		   // Noise_Ima = Data.sigma();
                   if (Verbose == True) 
		           cout << "Sigma Noise = " << Noise_Ima << endl;
               }
	       Cur.set_noise_model_gaussian(N_Sigma, Noise_Ima);
 	       break; 
      case CUR_NOISE_MULT_CORREL:
               noise_log_transform (Data, Data);
	       Noise_Ima = 1.;
	       MAD_Estim = True;
	       Cur.set_noise_level(Data, N_Sigma, MAD_Estim);
	       break; 
      case CUR_NOISE_LOG_RAYLEIGH:
               NoiseData.alloc(Nl,Nc,"noise");
               im_noise_rayleigh(NoiseData,  NbrImage);
	       noise_log_transform (NoiseData, NoiseData);
	       io_write_ima_float("xx_logr.fits", NoiseData, &Header);
  	       Cur.set_noise_level(NoiseData, N_Sigma);
	       noise_log_transform (Data, Data);
	       Noise_Ima = 0.641;
	       break;
       default:
               cout << "Error: unknown noise model ... " << endl;
	       exit(-1);
          break;
    }

    Ifloat Filter(Data.nl(), Data.nc(), "filter");
    Ifloat *Trans = NULL;
    Ifloat *TransSup = NULL;
    // Ifloat *TransNoise = NULL;

    if (Max_Iter <= 1)
    {
       if (Verbose == True) cout << "Transformation ... " << endl;
//        fltarray T;
//        Cur.transform(Data,T);
//        Cur.recons(T,Data);
//        io_write_ima_float("pr.fits", Data);
//        exit(0);
       Cur.transform(Data,Trans);
       if (Verbose == True) cout << "Filtering ... " << endl;
       if (FDR == True)  Cur.get_fdr_detect_level(Trans,  N_Sigma, True);
       else  cur_threshold (Cur, Trans,  N_Sigma,  Noise_Ima, Verbose, Interval);
       if (Verbose == True) Cur.print_info_noise();
       if (Verbose == True) cout << "Reconstruction ... " << endl;
       Cur.recons(Trans,Filter);
    }
    else 
    {
      Bool MSVST = Cur.MSVST;
      if (MSVST == True)
      {
         Cur.transform(Data,TransSup);
	 if (FDR == True)  Cur.get_fdr_detect_level(TransSup,  N_Sigma, True);
         else  cur_threshold (Cur, TransSup,  N_Sigma,  Noise_Ima, Verbose, Interval);
         Cur.VarNorm =  True;
         Cur.MSVST = False;
      }
       
      Ifloat Resi(Nl, Nc, "Resi");
      Ifloat ImaGrad;
      float LambdaTV = RegulParam;
      Bool Verb= Verbose;
      // if (RegulParam > 0.5) alpha = 1. / (2.*RegulParam);
      Filter.init();
      Resi = Data;
      Interval = True;
      for (int Iter = 0; Iter < Max_Iter; Iter ++)
      {
          float LambdaTVIter = LambdaTV - LambdaTV / (float) (Max_Iter-1)  * Iter;
          if ((Verbose == True) && (RegulParam == 0))cout << "Iter " << Iter + 1 << ", Sigma Resi " << sigma(Resi) << endl;
	  else if ((Verbose == True) && (RegulParam > 0))
	     cout << "Iter " << Iter + 1 << ", Sigma Resi " << sigma(Resi) << " Lambda = " << LambdaTVIter << endl;
          Cur.transform(Resi,Trans);
	  if (MSVST == False)
          {
	     if ((FDR == True) && (Iter == 0)) Cur.get_fdr_detect_level(Trans, N_Sigma, False);
 	     cur_threshold (Cur, Trans,  N_Sigma,  Noise_Ima, Verb, Interval);
	  }
	  else cur_support_threshold (Cur,  Trans,TransSup);
	  Verb = False;
	  Cur.recons(Trans,Resi);
          Filter += Resi;
          if (LambdaTVIter > 0)  RI.im_soft_threshold(Filter, Filter, LambdaTVIter, Noise_Ima);
  	  for (i=0;i < Nl;i ++)
          for (j=0;j < Nc;j ++) Resi(i,j) = Data(i,j) - Filter(i,j);
      }
   }

 
    switch (Stat_Noise)
    {
      case CUR_NOISE_GAUSSIAN:
      case CUR_NOISE_CORREL:
	      break;
       case CUR_NOISE_POISSON:  
	        noise_inverse_poisson_transform (Filter, Filter);
 		break;
       case CUR_NOISE_LOG_RAYLEIGH:
       case CUR_NOISE_MULT:
       case CUR_NOISE_MULT_CORREL:
 		 noise_inv_log_transform(Filter, Filter);
 		break;
        default:
	   break;  
    }
    
    if (PositivSol == True) threshold(Filter);
    io_write_ima_float(Name_Imag_Out, Filter, &Header);
    exit(0);
}

