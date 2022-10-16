/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  97/07/03
**    
**    File:  mr_proba.cc
**
*******************************************************************************
**
**    DESCRIPTION  Calculates the probability of each wavelet coefficient
**    -----------  not to be due to signal 
**                 
**    PARAMETRES    
**    ----------    
** 
**
**  
**
**
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"

char Name_Imag_In[256]; /* input file image */
char Name_Tab_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
Bool UseNSigma =False;

char Name_RMSMap[256]; 
Bool UseRMSMap=False;
float EpsilonPoisson = DEFAULT_EPSILON;

int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;

/****************************************************************************/

static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_TRANSFORM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                           StringTransform((type_transform)i));
    manline();  


    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "              default is automatically estimated\n");
    manline();  

    fprintf(OUTMAN, "         [-c gain,sigma,mean]\n");
    fprintf(OUTMAN, "             gain = gain of the CCD\n");
    fprintf(OUTMAN, "             sigma = read-out noise standard deviation\n");
    fprintf(OUTMAN, "             mean = read-out noise mean\n");
    fprintf(OUTMAN, "               noise = poisson + readout noise. default is no (Gaussian)\n");
    manline();  

    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
    manline();  

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform\n");
    fprintf(OUTMAN, "              default is %d\n", Nbr_Plan);
    fprintf(OUTMAN, "              default is %d in case of poisson noise with few events\n", DEF_N_SCALE);
    manline();  

    fprintf(OUTMAN, "         [-S SizeBlock]\n");
    fprintf(OUTMAN, "             Size of the  blocks used for local variance estimation.\n");
    fprintf(OUTMAN, "             default is %d \n", SizeBlock);
    manline();  

    fprintf(OUTMAN, "         [-N NiterSigmaClip]\n");
    fprintf(OUTMAN, "             iteration number used for local variance estimation.\n");
    fprintf(OUTMAN, "             default is %d \n", NiterClip);
    manline();  

    fprintf(OUTMAN, "         [-R RMS_Map_File_Name]\n");
    fprintf(OUTMAN, "              RMS Map.  If this Option is set, \n");
    fprintf(OUTMAN, "              the noise model is automatically fixed to:\n");
    fprintf(OUTMAN, "                 %s\n\n", StringNoise(NOISE_NON_UNI_ADD));
    manline();  

    vm_usage();
    manline();
    verbose_usage();    
    manline();
    
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
    while ((c = GetOpt(argc,argv,"m:t:g:c:n:S:N:R:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(stderr, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
                TransfOpt = True;
		break;
            case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_NOISE)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(stderr, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(stderr, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		break;
             case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(stderr, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(stderr, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(stderr, "bad number of scales: %s\n", OptArg);
		    fprintf(stderr, "1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		NscaleOpt = True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(stderr, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                UseNSigma =True;
                if ((N_Sigma <= 0.) || (N_Sigma > 100.)) 
                                        N_Sigma = DEFAULT_N_SIGMA;
		break;
 	   case 'R':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
		   fprintf(stderr, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseRMSMap = True;
 		break;
	 case 'E':
		if (sscanf(OptArg,"%f",&EpsilonPoisson) != 1) 
                {
		    fprintf(stderr, "Error: bad precision number: %s\n", OptArg);
		    exit(-1);
		}
		
		if ((EpsilonPoisson <MIN_EPSILON) || (EpsilonPoisson > MAX_EPSILON)) 
                {
		    fprintf(stderr, "Error: bad precision number: %s\n", OptArg);
		    fprintf(stderr, "%f <= Precision <= %f\n",MIN_EPSILON, MAX_EPSILON);
 		    exit(-1);
		}
		break; 
	 case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(stderr, "Error: bad block size: %s\n", OptArg);
		    exit(-1);
		}
		break; 
	 case 'N':
		if (sscanf(OptArg,"%d",&NiterClip) != 1) 
                {
		    fprintf(stderr, "Error: ad it. number for the 3sigma clipping: %s\n", OptArg);
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

	if (OptInd < argc) strcpy(Name_Tab_Out, argv[OptInd++]);
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
        if ((Stat_Noise == NOISE_CORREL) && (UseRMSMap != True))
        {
           cerr << endl << endl;
           cerr << "  Error: this noise model need a noise map (-R option) " << endl;
           exit(-1);
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
           // Stat_Noise = NOISE_NON_UNI_ADD;
        }
        if ((isotrop(Transform) == False)
              && ((Stat_Noise == NOISE_NON_UNI_ADD) ||
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


/****************************************************************************/

int main(int argc, char *argv[])
{
    int s,i,k;
    Ifloat Imag;

     /* support image creation */
    fitsstruct Header;
    char Cmd[256];
    extern softinfo Soft;

    Cmd[0] = '\0';
    Soft.mr2();
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR2);
    filtinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Tab_Out << endl;
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
   }

    io_read_ima_float(Name_Imag_In, Imag, &Header);
    Header.origin = Cmd;

    MRNoiseModel ModelData(Stat_Noise, Imag.nl(), Imag.nc(), 
                           Nbr_Plan, Transform);
    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
          for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
     for (s=0; s < NbrBand; s++) ModelData.TabEps[s] = EpsilonPoisson;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;

    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }

    MultiResol MR_Data(Imag.nl(), Imag.nc(),Nbr_Plan,Transform, "MRNoiseModel");
 
    MR_Data.transform(Imag);
    int NbrStatPerBand = 4; // moment of order 2,3,4 + multiscale entropy
    fltarray TabStat(NbrBand-1, NbrStatPerBand);   
    for (i=0; i < NbrBand-1; i++) 
    {
        int N = MR_Data.band(i).nl()*MR_Data.band(i).nc();
	double Mean, Sigma, Skew, Curt;
	float  Min, Max;
        moment4(MR_Data.band(i).buffer(),  N,  Mean,  Sigma, 
                Skew, Curt, Min, Max);
        TabStat(i, 0) = (float) Sigma;
	TabStat(i, 1) = (float) Skew;
	TabStat(i, 2) = (float) Curt;
    }
    
    ModelData.model(Imag, MR_Data);
    // calculate the probability p(w) of each wavelet coef. to be due to 
    // signal
    ModelData.prob_noise(MR_Data, True);
    if (Verbose == True) cout << "STAT = " << endl;
    for (int b=0; b< NbrBand-1; b++) 
    {
        int Nlb = MR_Data.band(b).nl();
	int Ncb = MR_Data.band(b).nc();
	double MeanProb=0.;
        for (i=0; i < Nlb; i++)
	for (int j=0; j < Ncb; j++) MeanProb += MR_Data(b,i,j);
	MeanProb /= (float)(Nlb*Ncb);
 	TabStat(b, 3) = MeanProb;
	
	if (Verbose == True)
            printf("  Band %d: Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, MEntrop = %5.3f\n",
	           b+1,TabStat(b,0),TabStat(b, 1),TabStat(b,2),MeanProb);
    }
    fits_write_fltarr(Name_Tab_Out, TabStat); 
    exit(0);
} 

