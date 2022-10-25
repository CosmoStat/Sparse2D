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
**    Date:  15/02/99
**    
**    File:  mr1d_pix.cc
**
*******************************************************************************
**
**    DESCRIPTION  analyse one pixel by multiresolution. Default is the 
**    -----------  last point.
**                 
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "IM1D_IO.h"
#include "MR1D_NoiseModel.h"
#include "Usage.h"

char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];

int NbrScale=0;           // number of scales
float NSigma=DEF_NSIG;   // number of sigma (for the noise) 
float SigmaNoise=0.;       // noise standard deviation
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   // type of noise
type_trans_1d Transform=TO1_PAVE_HAAR;    // type of transform
float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean  

Bool UseNSigma=False;
Bool SetPlan = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
int Pos = -1;
 
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

    fprintf(OUTMAN, "Usage: %s options input \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");


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

    nsigma_usage(NSigma);
    manline();
  
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    manline();

    fprintf(OUTMAN, "         [-x Position]\n");
    fprintf(OUTMAN, "             Position to analyse. Default is the last point.\n");
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
    while ((c = GetOpt(argc,argv,"x:m:g:c:n:s:v")) != -1) 
    {
	switch (c) 
        {
	   case 'x':
	   	if (sscanf(OptArg,"%d",&Pos ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad position: %s\n", OptArg);
	            exit(-1);
		}
		if (Pos < 2)
		{
		    fprintf(OUTMAN, "Error: bad  position: %s\n", OptArg);
		    fprintf(OUTMAN, "        1 < x <= N_pix \n");
	            exit(-1);
		}
		break;
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
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
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
                SigmaNoise  = 1.;
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
                SetPlan=True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) 
                {
		    fprintf(OUTMAN, "bad NSigma: %s\n", OptArg);
		    usage(argv);
		}
                if (NSigma <= 0.) NSigma = DEF_NSIG;
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
    fltarray Data;
    char Cmd[256];
    int i,b,x;
    // extern softinfo Soft;

    lm_check(LIC_M1D);
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    // lm_check(LIC_MR2);
    filtinit(argc, argv);
    
    // fits_read_fltarr(Name_Imag_In, Data);
    io_1d_read_data(Name_Imag_In, Data);

    reform_to_1d(Data);
    // FitsHeader.origin = Cmd;

    int Nx = Data.nx();  
    fltarray Result (Nx);
    if (!SetPlan) NbrScale = max_scale_number(Nx);
    int Last = Nx-1;
    if (Pos > 1) Last = Pos-1;
    if (Last > Nx-1)
    {
       cout << "Error: position must less than " << Nx << endl;
       exit(-1);
    }

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "Number of scales = " << NbrScale << endl;
       if (Stat_Noise == NOISE_GAUSSIAN)
            cout << "Type of Noise = GAUSSIAN" << endl;
       else cout << "Type of Noise = POISSON" << endl;
       cout << "Sigma Noise = " << SigmaNoise << endl;
       cout << "NSigma = " << NSigma << endl;
       cout << "naxis = " << Data.naxis() << endl;
       cout << "Min = " << Data.min() << endl;
       cout << "Max = " << Data.max() << endl;
       cout << "sigma = " << Data.sigma() << endl;
       cout << "Pixel position to analyse = " << Last + 1 << endl;
    }
 
 
MR_1D MR_Data (Nx, Transform, "MR_Data", NbrScale);
MR_Data.Border = I_MIRROR;
 
  MR1DNoiseModel NoiseModel(Stat_Noise, Nx, NbrScale, Transform);
  if (SigmaNoise > FLOAT_EPSILON) NoiseModel.SigmaNoise = SigmaNoise;
  if (UseNSigma  == True)
          for (i=0; i < NbrScale; i++) NoiseModel.NSigma[i]=NSigma;
  NoiseModel.CCD_Gain = CCD_Gain;
  NoiseModel.CCD_ReadOutSigma = CCD_ReadOutSigma;
  NoiseModel.CCD_ReadOutMean = CCD_ReadOutMean;
  // NoiseModel.SizeBlockSigmaNoise = SIZE_BLK_NOISE_ADD_NON_UN;
  
  NoiseModel.model(Data, MR_Data);
  if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN))
      cout << "Noise standard deviation = " << NoiseModel.SigmaNoise << endl;
  // NoiseModel.threshold(MR_Data);

 enum status_coeff {NO_DETECT,   // no detection
                    POS_DET,     // positive structure 
                    NEG_DET,     // negative structure
                    NEW_POS_DET, // new positive structure: upward transition 
                    NEW_NEG_DET, // new negative struct.: downward transition
                    END_POS_DET, // end of positive structure
                    END_NEG_DET  // end of negative structure
                   }; 

  for (b=0; b < NbrScale-1; b++)
  {
      status_coeff Status=NO_DETECT;
      int StartPos=0;
      int LastPos=0;
      float Coef = MR_Data(b,Last);
      Bool Signif = NoiseModel(b, Last);
      float Sig = NoiseModel.sigma(b, Last);
      float CoefBef = MR_Data(b,Last-1);
      Bool SignifBef = NoiseModel(b, Last-1);
 
      // test if the last coefficient is significant
      if (Signif == False) 
      {
          // if not, test if it corresponds to the end of a detected
	  // structure
          if (SignifBef == True)  
          {
             x = Last-1;
             while ((x > 0) && (NoiseModel(b, x) == True)  
                    && (CoefBef*MR_Data(b,x) > 0))  x--;
             Status = (CoefBef > 0) ? END_POS_DET: END_NEG_DET;
             StartPos = x+2;
             LastPos = Last;
          }
       }
       else
       {
          // if the previous coeff is not detected, it corresponds 
          // to a new detection
          if ((SignifBef == True) && (CoefBef*Coef > 0))
	  {
             x = Last-1;
             while ((x > 0) && (NoiseModel(b, x) == True)
                    && (Coef*MR_Data(b,x) > 0))  x--;
             Status = (Coef > 0) ? POS_DET: NEG_DET;
             StartPos = x+2;
             LastPos = Last+1;
          }
	  else  Status = (Coef > 0) ?  NEW_POS_DET: NEW_NEG_DET;
       }
 
      cout << "SCALE " << b+1;
      // cout << "  Coef = " << Coef << " Sigma = " << Sig << endl;    
      switch (Status)
      {
       case NO_DETECT: 
             cout << "  no detection "; 
             break;
       case POS_DET:
             cout << "  Positive significant structure [";
             cout  <<  StartPos << "," << Last+1 << "] ";
             break;
       case NEG_DET:
             cout << "  Negative significant structure [";
             cout <<  StartPos << "," << Last+1 << "] ";
             break;
       case NEW_POS_DET:
              cout << "  NEW Upward detection ";
             break;
       case NEW_NEG_DET :
              cout << "  NEW Downward transition: ";
             break;
       case END_POS_DET:
             cout << "  End of significant structure started at position " << StartPos;
             break;
       case END_NEG_DET:
             cout << "  End of significant structure started at position " << StartPos;
             break;
      }
      if (((Status == NEW_POS_DET) || (Status == NEW_NEG_DET)) &&
        ((Stat_Noise == NOISE_GAUSSIAN) || (Stat_Noise == NOISE_NON_UNI_ADD) 
          || (Stat_Noise == NOISE_UNDEFINED) || (Stat_Noise == NOISE_UNI_UNDEFINED))) 
                 cout << " SNR = " << ABS(Coef) / Sig << endl;
      else cout << endl;
   }
   exit(0);
} 

