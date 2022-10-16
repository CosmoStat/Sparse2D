/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  mr1d_detect.cc
**
******************************************************************************/
 

#include "Array.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "MR1D_Filter.h"
#include "NR.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
// #include "IM_Errors.h"
#include "MR1D_NoiseModel.h"
#include "MR1D_Segment.h"

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

char Name_Input[80]; /* input file image */
char Name_Output[80];
 
int NbrScale = 5;
int FirstScale= -1;
int LastScale= -1;
Bool OnlyPosDetect=False;
Bool OnlyNegDetect=False;
Bool DetectOnlyPosDetect=False;
Bool DetectOnlyNegDetect=False;

Bool UseNSigma=False;
Bool WriteAll = False;
int RecIterNumber = MAX_REC_ITER_OBJ1D;

float NSigma=DEF_NSIG;   /* number of sigma (for the noise) */;
float SigmaNoise=0.;
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */;
type_trans_1d Transform=TO1_PAVE_B3SPLINE;

float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean 
Bool Verbose = False;

#define SIZE_BLK_NOISE_ADD_NON_UN 20

/***********************************************************************/

static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN,"\n\nUsage: %s options in_signal out_signal\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    nbr_scale_usage(NbrScale);
    manline();
    
    nsigma_usage(NSigma);
    manline();   
  
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-2; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             Default is Gaussian noise.\n\n");    
    
    gauss_usage();
    manline();  
    
    ccd_usage();
    manline();    
    mr1ddetectr_option_usage(RecIterNumber);
    verbose_usage();
    manline();
    manline();
    manline();

    exit(-1);
}

 
/*********************************************************************/


/* GET COMMAND LINE ARGUMENTS */
static void init(int argc, char *argv[])
{
    int c;
    Bool SetIter = False;

    /* get options */
    while ((c = GetOpt(argc,argv,"g:n:s:f:l:m:c:aewi:MAEv")) != -1) 
    {
	switch (c) 
      	{
           case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_NOISE-2)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	 case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
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
                SigmaNoise = 1.;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);
		    exit(-1);
		}
                if (NSigma <= 0.) NSigma = DEF_NSIG;
                UseNSigma = True;
		break;
	 case 'n':
 		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((NbrScale <= 1) || (NbrScale > MAX_SCALE_1D)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
		    exit(-1);
		}
		break;
	 case 'f':
		if (sscanf(OptArg,"%d",&FirstScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad first scale: %s\n", OptArg);
		    exit(-1);
		}
		break;
	 case 'l':
		if (sscanf(OptArg,"%d",&LastScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad first scale: %s\n", OptArg);
		    exit(-1);
		}
		break;
	 case 'i':
		SetIter = True;
		if (sscanf(OptArg,"%d",&RecIterNumber) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad iteration number: %s\n", OptArg);
		    exit(-1);
		}
		break;
	 case 'a':
		OnlyNegDetect = True;
		break;
	 case 'e':
		OnlyPosDetect = True;
		break;
	 case 'M':
		Transform = TM1_PAVE_MEDIAN;
		break;
	 case 'A':
		DetectOnlyNegDetect = True;
		break;
	 case 'E':
		DetectOnlyPosDetect = True;
		break;
	 case 'w':
		WriteAll = True;
		break;
	 case 'v':
		Verbose = True;
		break;
	   case '?':
		usage(argv);
	}
   }

   if ( (Transform == TM1_PAVE_MEDIAN) && (SetIter == True))
   {
       fprintf(OUTMAN, "Error: iteration number option is not compatible with the median transform ...\n");
       exit(-1);
   }

   /* get input filename  */
   	if (OptInd < argc) strcpy(Name_Input, argv[OptInd++]);
   	else usage(argv);


   /* get output filename  */
   	if (OptInd < argc) strcpy(Name_Output, argv[OptInd++]);
         else usage(argv);

   /* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Error: oo many parameters: %d %d %s ...\n", 
                       OptInd, argc, argv[OptInd]);
		exit(-1);
	}
}

/************************************************************************/

int main(int argc, char *argv[])
{
   fltarray Data;
   fitsstruct FitsHeader;
   int i;
   char Cmd[256];
    extern type_1d_format IO_1D_Format;

   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
   /* get command line */
   lm_check(LIC_M1D);
   init(argc,argv);

    
   // read the data
  io_1d_read_data(Name_Input, Data, &FitsHeader);
  reform_to_1d(Data);
  FitsHeader.origin = Cmd;
 
  int Nx = Data.nx();  
  if ( (int) Transform < NBR_DIADIC_TRANS_1D)
  {
     int ScaleMax=iround(log((float)Nx/4.) / log(2.)+ 1.);
     if (NbrScale > ScaleMax)
     {
        cout << "Error: maximun Number of scales allowed: " << ScaleMax << endl;
        exit(-1);
     }  
  }
  fltarray Result(Nx);
 
  MR1DNoiseModel NoiseModel(Stat_Noise, Nx, NbrScale, Transform);
  if (SigmaNoise > FLOAT_EPSILON) NoiseModel.SigmaNoise = SigmaNoise;
  if (UseNSigma  == True)
          for (i=0; i < NbrScale; i++) NoiseModel.NSigma[i]=NSigma;
  if (DetectOnlyPosDetect == True) NoiseModel.OnlyPositivDetect = True;
  if (DetectOnlyNegDetect == True) NoiseModel.OnlyNegativDetect = True;
  NoiseModel.CCD_Gain = CCD_Gain;
  NoiseModel.CCD_ReadOutSigma = CCD_ReadOutSigma;
  NoiseModel.CCD_ReadOutMean = CCD_ReadOutMean;
  NoiseModel.SizeBlockSigmaNoise = SIZE_BLK_NOISE_ADD_NON_UN;
  
  MR_1D MR_Data(Nx, Transform, "MR_Data", NbrScale);
  MR_Data.Border = I_MIRROR;

  NoiseModel.model(Data, MR_Data);
  if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN))
      cout << "Noise standard deviation = " << NoiseModel.SigmaNoise << endl;
  NoiseModel.threshold(MR_Data);
  //MR_Data.recons( Data);
  //fits_write_fltarr("xx_filter",  Data);
 
  MR1D_Seg DataSeg;
  DataSeg.init(MR_Data);
  DataSeg.Verbose = Verbose;
  DataSeg.segment(MR_Data, NoiseModel);
 
  if (Verbose == True) DataSeg.print_obj_info();
  if (FirstScale >= 0) DataSeg.FirstScale = FirstScale;
  if (LastScale >= 0) DataSeg.LastScale =  LastScale;
  DataSeg.PosDetect = OnlyPosDetect;
  DataSeg.NegDetect = OnlyNegDetect;
  DataSeg.MaxRecIter = RecIterNumber;
  fltarray Tab_SignalRec, Tab_Add;
 
  DataSeg.recons(MR_Data, Tab_SignalRec, Tab_Add, Data, NoiseModel);
  if (Tab_SignalRec.n_elem() == 0) 
  {
     if (Verbose == True) cout << "No detection " << endl;
     Tab_SignalRec.alloc(1);
     Tab_SignalRec(0)=-1;
  }
  
  // MR_Data.rec_adjoint(Result);
  if (WriteAll == True)
  {
      // fits_write_fltarr("tabadd.fits", Tab_Add);
      io_1d_write_data("tabadd", Tab_Add, &FitsHeader);
      DataSeg.write_data_seg("tabseg");
  }
  if (IO_1D_Format == F1D_FITS)
       fits_write_fltarr(Name_Output, Tab_SignalRec);
  else io_write2d_ascii(Name_Output, Tab_SignalRec);  
  exit(0);
}


/************************************************************************/
