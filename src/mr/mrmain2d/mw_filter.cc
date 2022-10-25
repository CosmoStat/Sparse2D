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
**    Date:  98/01/13
**    
**    File:  mw_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION: Image filtering using the multi-scale entropy
**    ----------- 
**                 
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Filter.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "CErf.h"
#include "CMem.h"
#include "MR_Edge.h"
#include "MW_Filter.h"


char Name_Imag_In[256]; /* input file image */
char Name_MR_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_MR_TRANS; /* type of transform */
Bool UseNSigma =False;
Bool SupIsol=False;             /* suppress isolated pixel in the support */

char Name_RMSMap[256]; 
Bool UseRMSMap=False;
float EpsilonPoisson = DEFAULT_EPSILON;

float RegulVal =1.;
float CvgParam= DEF_MEM_FILT_CONGER;

int NbrIter = DEF_MEM_FILT_MAX_ITER;

int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
   
// #define WRITE_PARAM 1
#define PREC 0.000001

Bool  DataSNR = False;
Bool  DataEdge = False;
Bool  RegulPerScale = False;
Bool Verbose = False;
Bool Positiv = False;
Bool DetectOnlyPos = False;

int TypeOpt=DEF_MEM_ALPHA_CST;
int FirstScale = DEF_FIRST_DETECT_SCALE;


/****************************************************************************/

void mr_getmodel_from_edge(MultiResol & MR_Data, MultiResol &MR_Edge,
                           MRNoiseModel & ModelData)
/* Find a multiresolution model for edges:
    -- create an SNR edge image: ImaEdge
    -- threshold ImaEdge if ImaEdge < NSigma
    -- kill isolated pixel in ImaEdge
    -- in each band, 
        . threshold wavelet coefficient if there is no edge
        .average the wavelet coefficient
         with the two other values in the edge direction
   !! This routine works only if the a-trou algorithm is choosen
*/

{
   int Nbr_Band =  MR_Data.nbr_band()-1;
   int i,j,Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   Ifloat ImaEdge(Nl,Nc,"ImaEdge");
   Iint ImaAngle(Nl,Nc,"ImaAngle");

   mr_get_edge( MR_Data,  MR_Edge,  ImaEdge, ImaAngle,  ModelData);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)
   {
      ImaEdge(i,j) /= ModelData.NSigma[1];
      if (ImaEdge(i,j) > 1) ImaEdge(i,j) = 1.;
      else ImaEdge(i,j) = 0.;
   }
   for (i=0;i < Nl; i++)
   for (j=0;j < Nc; j++) 
       if (isolated_pixel(ImaEdge,i,j,I_MIRROR) == True) ImaEdge(i,j) = 0.;  

   for (int b = 0; b < Nbr_Band; b++) 
   {
      int Nlb = MR_Data.size_band_nl(b);
      int Ncb = MR_Data.size_band_nc(b);
      int Step = (MR_Data.band_to_scale(b) == b) ? b : 0;
  
      for (i = 0; i < Nlb; i++)
      for (j = 0; j < Ncb; j++)
      {
         if (ImaEdge(i,j) > FLOAT_EPSILON)
               MR_Edge(b,i,j) = val_contour_min(MR_Data.band(b), ImaEdge, ImaAngle,
                                              i,j, FLOAT_EPSILON, Step);
         else  MR_Edge(b,i,j) = 0.;
      }
   }
 /*
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++)

   for (int b = 0; b < Nbr_Band; b++) 
   {
      int Nlb = MR_Data.size_band_nl(b);
      int Ncb = MR_Data.size_band_nc(b);
      int Step = (MR_Data.band_to_scale(b) == b) ? b : 0;
Step=0;
      ImaEdge.resize(Nlb,Ncb);
      ImaAngle.resize(Nlb,Ncb);
      mr_band_edge(MR_Data, ImaEdge, ImaAngle, b, ModelData);
      // if (b == 0) io_write_ima_float("xx_edge.fits", ImaEdge);
      for (i = 0; i < Nlb; i++)
      for (j = 0; j < Ncb; j++)
      {
         int s = MR_Data.band_to_scale(b);
	 ImaEdge(i,j) /= ModelData.NSigma[s];
	 if (ImaEdge(i,j) > 1) ImaEdge(i,j) = 1.;
	 else ImaEdge(i,j) = 0.;
      }

      for (i=0;i < Nlb; i++)
      for (j=0;j < Ncb; j++) 
      {
         // supress isolated pixel
         if (isolated_pixel(ImaEdge,i,j,I_MIRROR) == True) ImaEdge(i,j) = 0.;       
      	 
	 // set the model edge
         if (ImaEdge(i,j) > FLOAT_EPSILON)
               MR_Edge(b,i,j) = val_contour_min(MR_Data.band(b), ImaEdge, ImaAngle,
                                              i,j, FLOAT_EPSILON, Step);
         else  MR_Edge(b,i,j) = 0.;
	 
	 
	 // MR_Edge(b,i,j) *= ImaEdge(i,j);                      
      }    
   }
*/
}



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
    for (i = 0; i < NBR_NOISE-1; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             default is Gaussian noise\n");
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform\n");
    fprintf(OUTMAN, "              default is %d\n", Nbr_Plan);
//     fprintf(OUTMAN, "              default is %d in case of poisson noise with few events\n", DEF_N_SCALE);
    manline();

    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "              Thresholding at nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "              default is %f\n", N_Sigma);
    manline();

//     fprintf(OUTMAN, "         [-E Epsilon]\n");
//     fprintf(OUTMAN, "             Epsilon = precision for computing thresholds\n");
//     fprintf(OUTMAN, "                       (only used in case of poisson noise with few events)\n");
//     fprintf(OUTMAN, "             default is %e \n", DEFAULT_EPSILON);

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
    fprintf(OUTMAN, "                 %s\n", StringNoise(NOISE_NON_UNI_ADD));
    manline();

    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              default is %f\n", RegulVal);
     manline();
   
    fprintf(OUTMAN, "         [-C CvgParam]\n");
    fprintf(OUTMAN, "              Convergence parameter \n");
    fprintf(OUTMAN, "              default is %f\n", CvgParam);
     manline();
   
    fprintf(OUTMAN, "         [-T Type_of_Regularization]\n");
    fprintf(OUTMAN, "              1: Use a fixed user Alpha value.\n");
    fprintf(OUTMAN, "              2: Estimate the optimal Alpha.\n");
    fprintf(OUTMAN, "              3: Estimate one  Alpha value per band.\n");
    fprintf(OUTMAN, "              default is 1.\n");
      manline();
  
//     fprintf(OUTMAN, "         [-A]\n");
//     fprintf(OUTMAN, "              Use the edge as a-priori model information.\n");
//     fprintf(OUTMAN, "              default is no.\n");
    
    fprintf(OUTMAN, "         [-D]\n");
    fprintf(OUTMAN, "              Alpha is modified using the data SNR.\n");
    fprintf(OUTMAN, "              default is no.\n");
    manline();
 
    fprintf(OUTMAN, "         [-i MaxIter]\n");
    fprintf(OUTMAN, "              Maximum number of iterations.\n");
    fprintf(OUTMAN, "              default is %d.\n", NbrIter);  
     manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "              Apply the positivity constraint.\n"); 
    manline();
    first_detect_scale_usage();
    manline();
    kill_isol_pix_usage();
    manline();
    vm_usage();
    manline();   
         
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "              Verbose\n"); 
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
    while ((c = GetOpt(argc,argv,"F:km:t:g:c:n:s:S:N:R:C:G:i:T:ADvpPzZ:")) != -1) 
    {
	switch (c) 
        {
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
	   case 'k':SupIsol = True;break;
	   case 'v': Verbose = True; break;
 	   case 'P': Positiv = True; break;
 	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "bad type of multiresolution transform: %s\n", OptArg);
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
                if ((c > 0) && (c < NBR_NOISE)) 
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
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
 		break;
	   case 'C':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&CvgParam) != 1) 
                {
		    fprintf(stderr, "Error: bad Cvg Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (CvgParam <= 0.) CvgParam = 0.1;
		break;
	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulVal) != 1) 
                {
		    fprintf(stderr, "Error: bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulVal < 0.) RegulVal = 1.;
		break;
           case 'T' :
                if (sscanf(OptArg, "%d", &TypeOpt) !=1)
                {
                    fprintf(stderr, "bad type of regularization %s\n", OptArg);
                    exit(-1);
                }
		if ((TypeOpt < 1) || (TypeOpt > 3))
		{
                   cerr << "Error: regularization type must be in [1,3] ... " << endl;
                   exit(-1);
                }
                break ;
            case 'A' : DataEdge = True; break;
	    case 'D' : DataSNR = True;  break;
 	    case 'i':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrIter) != 1) 
                {
		    fprintf(stderr, "Error: bad number of iterations: %s\n", OptArg);
		    exit(-1);
		}
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
// 	 case 'E':
// 		if (sscanf(OptArg,"%f",&EpsilonPoisson) != 1) 
//                 {
// 		    fprintf(stderr, "Error: bad precision number: %s\n", OptArg);
// 		    exit(-1);
// 		}
// 		
// 		if ((EpsilonPoisson <MIN_EPSILON) || (EpsilonPoisson > MAX_EPSILON)) 
//                 {
// 		    fprintf(stderr, "Error: bad precision number: %s\n", OptArg);
// 		    fprintf(stderr, "%f <= Precision <= %f\n",MIN_EPSILON, MAX_EPSILON);
//  		    exit(-1);
// 		}
// 		break; 
	 case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(stderr, "bad block size: %s\n", OptArg);
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

	if (OptInd < argc) strcpy(Name_MR_Out, argv[OptInd++]);
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

/*********************************************************************/
/****************************************************************************/

void iter_recons(MultiResol & MR_Data, MRNoiseModel & ModelData, int & FirstScale, 
                 Ifloat & Imarec);

void iter_recons(MultiResol & MR_Data, MRNoiseModel & ModelData, int & FirstScale, Ifloat & Imarec)

/* Find a reversible reconstruction */

{
   
   int Nbr_Band =  MR_Data.nbr_band()-1;
   int Nbr_Plan =  MR_Data.nbr_band();
   int i,j,b;
   int Niter=10;
   int k;
   int Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   Ifloat I(Nl,Nc,"I");
   Ifloat Res(Nl,Nc,"Residu");
   MultiResol MR_Data_i(I.nl(),I.nc(),MR_Data.nbr_scale(), Transform, "MRNoiseModel");
   MultiResol MR_Data_res(I.nl(),I.nc(),MR_Data.nbr_scale(), Transform, "MRNoiseModel");
 
   for (b = 0; b < FirstScale; b++) MR_Data.band(b).init();
   MR_Data.rec_adjoint(I);
   
   for (k = 0; k <  Niter; k++)
   {
   	MR_Data_i.transform(I);
     	for (b = 0; b < Nbr_Band; b++)
   	{	
   		for (i = 0; i < MR_Data.size_band_nl(b); i++)
                for (j = 0; j < MR_Data.size_band_nc(b); j++)
   		{
   		    if(ModelData(b,i,j) == True)
   				MR_Data_res(b,i,j) = MR_Data(b,i,j) - MR_Data_i(b,i,j);
   		    else MR_Data_res(b,i,j) = 0;
   		}
    	}
    	for (i = 0; i < MR_Data.size_band_nl(Nbr_Band); i++)
        for (j = 0; j < MR_Data.size_band_nc(Nbr_Band); j++)
    	   MR_Data_res(Nbr_Band,i,j) = MR_Data(Nbr_Band,i,j) - MR_Data_i(Nbr_Band,i,j);
 	MR_Data_res.rec_adjoint(Res);
	// MR_Data_res.recons(Res);
    	for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++)
   	{
   		// if (Res(i,j)>0) 
           I(i,j)=I(i,j)+Res(i,j);	
   	}
	// cout << " REC " << k+1 << " Max(rec) = " << I.max() << " Res.sigma = " << Res.sigma() << endl;
    }
    Imarec=I;
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

    Soft.mr2();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR2);
    
     /* Get command line arguments, open input file(s) if necessary */
    
    filtinit(argc, argv);

// CMemWave CM;
// fltarray TabTest(1000,8);
// CM.test_tab(TabTest, N_Sigma, DataSNR);
// fits_write_fltarr("tab.fits", TabTest);
// exit(-1);

    if (Verbose == True)
    {  
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_MR_Out << endl;
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
       if (FirstScale > 0)
          cout << "Start the detect at scale " << FirstScale+1 << endl;
    }

    // Read the input image
    io_read_ima_float(Name_Imag_In, Imag, &Header);
    Header.origin = Cmd;

    // Noise modeling class intialization
    MRNoiseModel ModelData(Stat_Noise, Imag.nl(), Imag.nc(), 
                           Nbr_Plan, Transform);
    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
          for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;
     for (s=0; s < NbrBand; s++) ModelData.TabEps[s] = EpsilonPoisson;
    if (SupIsol == True) ModelData.SupIsol = True;
    ModelData.FirstDectectScale = FirstScale;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;

    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }

    // Muliresolution object creation 
    MultiResol MR_Data(Imag.nl(), Imag.nc(),Nbr_Plan,Transform, "MRNoiseModel");
    if (Stat_Noise == NOISE_SPECKLE) ModelData.SigmaApprox = True;
    ModelData.model(Imag, MR_Data);

    // regul. parameter intialization
    int Nb = NbrBand-1;
    fltarray TabAlpha(Nb);
    for (int b=0; b < Nb; b++) TabAlpha(b) = RegulVal;
    
    // Model for edges.
    MultiResol  *MR_Model = NULL;
    if ((DataEdge == True) && (Transform == TO_PAVE_BSPLINE))
    {
        MR_Model = new MultiResol(Imag.nl(), Imag.nc(),Nbr_Plan,Transform,"Model");
        mr_getmodel_from_edge(MR_Data, *MR_Model, ModelData);
	(*MR_Model).write("xx_model.mr");
    }

    // apply the filtering
    mw_filter(MR_Data,  ModelData, TypeOpt, TabAlpha, DataSNR, 
	       DataEdge, MR_Model, CvgParam, NbrIter, Positiv, Verbose, True);
	
    // Image reconstruction      
    // MR_Data.recons(Imag);
    iter_recons(MR_Data, ModelData, FirstScale, Imag);
    if (Positiv == True)  threshold(Imag);
     
    if (ModelData.TransImag == True) ModelData.im_invtransform(Imag);

    // Save the result
    io_write_ima_float(Name_MR_Out,  Imag, &Header);
    
    exit(0);
} 

