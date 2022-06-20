/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  97/09/03
**    
**    File:  mr_comb_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
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
#include "CErf.h"
#include "CMem.h"
#include "MR_Edge.h"
#include "MW_Filter.h"


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
Bool UseNSigma =False;

char Name_RMSMap[256]; 
Bool UseRMSMap=False;
float EpsilonPoisson = DEFAULT_EPSILON;

float RegulParam =1.;
float CvgParam=0.01;

int NbrIter = 20;

int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
   
Bool  DataSNR = False;
type_transform MemTransform = DEFAULT_MR_TRANS; /* type of transform */
int TypeOpt=DEF_MEM_ALPHA_OPT_SCALE;

#define MAX_NBR_COMBINED_FIL_METHOD 10
int NbrCombFil=0;
type_transform TabMethod[MAX_NBR_COMBINED_FIL_METHOD];
Bool FilterErode = False;
Bool FilterMem = False;

Bool Verbose = False;

/****************************************************************************/

static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_transform_for_thresholding]\n");
    for (i = 0; i < NBR_TRANSFORM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransform((type_transform)i));
     manline();       

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Filtering by an opening (erosion+dilation)\n");
    fprintf(OUTMAN, "             The structural element is a circle of size 3\n");
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
    for (i = 0; i < 4; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
     manline();       

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform\n");
    fprintf(OUTMAN, "              default is %d\n", Nbr_Plan);
     manline();       

    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "              Thresholding at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "              default is %f\n", N_Sigma);
    manline();       

    fprintf(OUTMAN, "         [-M type_transform]\n");
    fprintf(OUTMAN, "              Multiscale entropy method (MEM) using the transform\n");
    fprintf(OUTMAN, "              defined by type_transform.\n");
      manline();       
  
    fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "              MEM method:\n");
    fprintf(OUTMAN, "              Alpha is modified using the data SNR.\n");
    fprintf(OUTMAN, "              default is no.\n");
     manline();       

    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              MEM method: Regularization parameter.\n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
     manline();       
      
    fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "              Use default combined methods\n");
    manline();
    vm_usage();
    manline();       
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n"); 
    fprintf(OUTMAN, "\n");
    
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");
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
    while ((c = GetOpt(argc,argv,"m:t:g:c:n:s:G:dOM:DvzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
	   case 'D': DataSNR = True; 
	             FilterMem = True;
		     break;
	   case 'M':
	         /* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                      MemTransform = (type_transform) (c-1);
                else  
                {
		    fprintf(stderr, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}  
		FilterMem = True;
 		break;
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
 		TabMethod[NbrCombFil++] = Transform;
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
		    fprintf(stderr, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(stderr, "1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(stderr, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if ((N_Sigma <= 0.) || (N_Sigma > 100.)) 
                                        N_Sigma = DEFAULT_N_SIGMA;
 		break;
	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "Error: bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam < 0.) RegulParam = 1.;
		FilterMem = True;
 		break;
	 case 'd':
	         TabMethod[NbrCombFil++] = TM_PAVE_MEDIAN;
                 TabMethod[NbrCombFil++] = TO_MALLAT;
	         TabMethod[NbrCombFil++] = TO_HAAR;
	         TabMethod[NbrCombFil++] = TO_PAVE_BSPLINE;
                 FilterMem = True;
 		 break;
	case 'O': FilterErode = True; break;
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
	if (FilterMem == True) TabMethod[NbrCombFil++] = MemTransform;

        if ((NbrCombFil == 0) && (FilterErode == False))
	{
	   fprintf(stderr, "Warning: at least one filtering method should be selected ... \n");
	   fprintf(stderr, "         Default option (-d) is taken.\n");
           TabMethod[NbrCombFil++] = TM_PAVE_MEDIAN;
           TabMethod[NbrCombFil++] = TO_MALLAT;
	   TabMethod[NbrCombFil++] = TO_HAAR;
	   TabMethod[NbrCombFil++] = TO_PAVE_BSPLINE;
           FilterMem = True;
           TabMethod[NbrCombFil++] = MemTransform;
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
		fprintf(stderr, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/****************************************************************************/

 
void mr_cfm(Ifloat &Imag,  Ifloat &Result, Bool Verbose)
{
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   type_transform Trans;
   int i,s;
   Ifloat Sol(Nl,Nc,"Filter");
   // MRNoiseModel *ModelData;
   MultiResol MR_Data;
   MRNoiseModel ModelData;
   
   for (i=0; i < NbrCombFil; i++)
   {
      Trans = TabMethod[i];
      ModelData.alloc(Stat_Noise, Imag.nl(), Imag.nc(), Nbr_Plan, Trans);
      if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
      if (UseNSigma  == True)
          for (s=0; s < Nbr_Plan; s++) ModelData.NSigma[s]=N_Sigma;
      for (s=0; s < Nbr_Plan; s++) ModelData.TabEps[s] = EpsilonPoisson;
      ModelData.NiterSigmaClip = NiterClip;
      ModelData.SizeBlockSigmaNoise = SizeBlock;

      MR_Data.alloc(Imag.nl(), Imag.nc(),Nbr_Plan, Trans, "MRNoiseModel");
      ModelData.model(Imag, MR_Data);
      if ((i < NbrCombFil-1) || (FilterMem == False))
      {
         if (Verbose == True)
           cout << i+1 << ":Filtering using " << StringTransform(Trans) << endl;
         ModelData.threshold(MR_Data);
      }
      else 
      {
         if (Verbose == True)
           cout << i+1 << ":Filtering using MEM method " << endl;
	            // MEM Filtering    
	 int NbrBand = ModelData.nbr_band();
	 int Nb = NbrBand-1;
	 MultiResol  *MR_Model = NULL;
         fltarray TabAlpha(Nb);
         for (int b=0; b < Nb; b++) TabAlpha(b) = RegulParam;
         mw_filter(MR_Data, ModelData, TypeOpt, TabAlpha, DataSNR, 
	           False, MR_Model, CvgParam, NbrIter, True, Verbose);
      }
      MR_Data.recons(Sol);
      // threshold(Sol);
      Result += Sol;
      ModelData.free();
      MR_Data.free();
   }
   
   // Filtering by morpho math
   if (FilterErode == True) 
   {
      Ifloat Erode(Nl, Nc, "Erode");
      int Elstr_Size=3;
      morpho_cercle_erosion (Imag, Erode, Elstr_Size);
      morpho_cercle_dilation (Erode, Sol, Elstr_Size);
      Result += Sol;
      NbrCombFil++;
   }
   
   for (i=0; i < Nl*Nc; i++) Result(i) /= (float) NbrCombFil; 
   threshold(Result);
}
 
/****************************************************************************/

int main(int argc, char *argv[])
{
    int  k;
    Ifloat Imag;

     /* support image creation */
    fitsstruct Header;
    char Cmd[256];
    extern softinfo Soft;

    Soft.mr2();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR2);
    
    /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

    if (Verbose == True)
    { 
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
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
       cout << "N_Sigma = " << N_Sigma << endl;
    }

    io_read_ima_float(Name_Imag_In, Imag, &Header);
    Header.origin = Cmd;

    Ifloat Result(Imag.nl(), Imag.nc(), "Result");
    mr_cfm( Imag,   Result, Verbose);
    io_write_ima_float(Name_Imag_Out, Result, &Header);
    
    exit(0);
} 

