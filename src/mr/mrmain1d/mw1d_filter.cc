/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/20
**    
**    File:  mw1d_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  spectrum filtering
**    ----------- 
**                 
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"

#include "MR1D_Obj.h"  // RG
#include "MR1D_Filter.h" // RG for type_noise
#include "MR1D_NoiseModel.h"
 
#include "CErf.h"
#include "CMem.h"
#include "MW1D_Filter.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];

int Nbr_Plan=0;           // number of scales
float N_Sigma=DEF_NSIG;   // number of sigma (for the noise) 
float Noise_Ima=0.;       // noise standard deviation
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   // type of noise
type_trans_1d Transform=TO1_PAVE_B3SPLINE;    // type of transform
float Epsilon= DEF_EPS_CVG_MEMW1D_FILTER;     // convergence parameter
int Max_Iter=DEF_MEMW1D_FILTER_MAX_ITER;      // maximum number of iterations
float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean  

Bool UseNSigma=False;
Bool SetPlan = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool DataSnr = False;
float RegulParam=1.;
Bool Verbose = False;
Bool WriteCoeff = False;

/*********************************************************************/

static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}
 

/***************************************************************************/
      
 
static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options input output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_DIADIC_TRANS_1D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf1D((type_trans_1d)i));
    fprintf(OUTMAN, "             Default is %s \n", StringTransf1D(Transform));
    manline();

    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-1; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             Default is Gaussian noise.\n\n");    
   
    gauss_usage();
    manline();  
         
    ccd_usage();
    manline();    

    nsigma_usage(N_Sigma);
    manline();
     
    max_iter_usage(Max_Iter);
    manline();
  
    converg_param_usage(Epsilon);
    manline();
  
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    manline();
    
    fprintf( OUTMAN, "         [-G RegulParam]\n");
    fprintf( OUTMAN, "              Regularization parameter \n");
    fprintf( OUTMAN, "              default is %f\n", RegulParam);
    manline();
    
    fprintf( OUTMAN, "         [-D]\n");
    fprintf( OUTMAN, "              Alpha is modified using the data SNR.\n");
    fprintf( OUTMAN, "              default is no.\n");
     
    manline();
    fprintf( OUTMAN, "         [-w FilterCoefFileName]\n");
    fprintf( OUTMAN, "              Write to the disk the filtered wavelet coefficient.\n");
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
   
    /* get options */
    while ((c = GetOpt(argc,argv,"m:t:g:c:n:s:i:e:G:Dvw:v")) != -1) 
    {
	switch (c) 
        {
           case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_NOISE-1)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            usage(argv);
                    
		}
                if ((c > 0) && (c <= NBR_DIADIC_TRANS_1D)) 
                                        Transform = (type_trans_1d) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            usage(argv);
 		}
		break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
		    usage(argv);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		break;
            case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &CCD_Gain,
                                     &CCD_ReadOutSigma, &CCD_ReadOutMean) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima  = 1.;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    usage(argv);
		}
                SetPlan=True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
		    usage(argv);
		}
                if (N_Sigma <= 0.) N_Sigma = DEF_NSIG;
		UseNSigma = True;
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "bad Max_Iter: %s\n", OptArg);
				usage(argv);
		}
                if (Max_Iter <= 0)  Max_Iter = DEF_ITER;
		break;
	   case 'e':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "bad convergence parameter: %s\n", OptArg);
		   usage(argv);
		}
              //  if ((Epsilon < 0) || (Epsilon > 1.))  Epsilon = DEF_EPS;
		break;
	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    usage(argv);
		}
                if (RegulParam < 0.) RegulParam = 1.;
		break;
	  case 'v': Verbose = True; break;
	  case 'w' : 
	        if (sscanf(OptArg,"%s", Name_Coeff) != 1) 
                {
		   fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		   usage(argv);
		}
		WriteCoeff = True;
 		break;
 	    case 'D' : DataSnr = True;  break;
  	    case '?': usage(argv); break;
		}
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

int main(int argc, char *argv[])
{
    fltarray Dat;
    fitsstruct FitsHeader;
    char Cmd[256];
    int i,b,x;
    extern softinfo Soft;

    Soft.mr2();
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    lm_check(LIC_M1D);
    filtinit(argc, argv);
    
    io_1d_read_data(Name_Imag_In, Dat, &FitsHeader);
    reform_to_1d(Dat);
    FitsHeader.origin = Cmd;

    int Nx = Dat.nx();  
    fltarray Result (Nx);
    if (!SetPlan) Nbr_Plan = max_scale_number(Nx);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
       cout << "Transform = " << StringTransf1D(Transform) << endl;
       cout << "Number of scales = " << Nbr_Plan << endl;
       if (Stat_Noise == NOISE_GAUSSIAN)
            cout << "Type of Noise = GAUSSIAN" << endl;
       else cout << "Type of Noise = POISSON" << endl;
       cout << "Sigma Noise = " << Noise_Ima << endl;
       cout << "N_Sigma = " << N_Sigma << endl;
       cout << "Epsilon = " << Epsilon << endl;
       cout << "Max_Iter = " << Max_Iter << endl;
       cout << "naxis = " << Dat.naxis() << endl;
       cout << "Min = " << Dat.min() << endl;
       cout << "Max = " << Dat.max() << endl;
       cout << "sigma = " << Dat.sigma() << endl;
    }

  MR_1D MR_Data(Nx, Transform, "MR_Data", Nbr_Plan);
  MR_Data.Border = I_MIRROR;
 
  MR1DNoiseModel NoiseModel(Stat_Noise, Nx,  Nbr_Plan, Transform);
  if (Noise_Ima > FLOAT_EPSILON) NoiseModel.SigmaNoise = Noise_Ima;
  if (UseNSigma  == True)
        for (i=0; i < Nbr_Plan; i++) NoiseModel.NSigma[i]=N_Sigma;
  NoiseModel.CCD_Gain = CCD_Gain;
  NoiseModel.CCD_ReadOutSigma = CCD_ReadOutSigma;
  NoiseModel.CCD_ReadOutMean = CCD_ReadOutMean;
  NoiseModel.model(Dat, MR_Data);

  fltarray TabAlpha(Nbr_Plan);
  for (b=0; b < Nbr_Plan; b++) TabAlpha(b) = 0.;
  MR_1D  *MR_Model = NULL;

  mw1d_filter (MR_Data, NoiseModel, TabAlpha,  RegulParam, 
               DataSnr, Epsilon,  Max_Iter, Verbose, MR_Model);  

  
  if (WriteCoeff == True)
  { 
      extern type_1d_format IO_1D_Format;
      fltarray TabMR;
      TabMR = MR_Data.image();
      if (IO_1D_Format == F1D_FITS)
            fits_write_fltarr(Name_Coeff, TabMR);
      else io_write2d_ascii(Name_Coeff, TabMR);
  }
  MR_Data.recons (Result);
   
  // iterative reconstruction for pyramidal transform
  if (which_set_is_trans1d(Transform) == TRANS1_PYR)  
  {
      MR_1D MR_Rec(Nx, Transform, "MR_Data", Nbr_Plan);
      fltarray Resi (Nx);

      for (i=0; i < 4; i++)
      {
         MR_Rec.transform(Result);
	 for (b = 0; b < Nbr_Plan; b++) 
         for (x = 0; x < MR_Rec.size_scale_np(b); x++)
	    MR_Rec(b,x) = MR_Data(b,x) - MR_Rec(b,x);
	 MR_Rec.recons (Resi);
	 Result += Resi;
	 if (Verbose == True) cout << "Rec iter " << i+1 << " sigma resi = " << Resi.sigma() << endl;
      }
  }
  if (NoiseModel.TransImag == True) NoiseModel.signal_invtransform(Result);
   if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN))
          cout << endl << "Noise Standard deviation = " <<  NoiseModel.SigmaNoise << endl;

  io_1d_write_data(Name_Imag_Out, Result, &FitsHeader);
  exit(0);
} 

