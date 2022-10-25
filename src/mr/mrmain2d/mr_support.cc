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
**    File:  mr_support.cc
**
*******************************************************************************
**
**    DESCRIPTION  create the multiresolution support of an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_support option in_image out_file.mr
**        where options = 
**
**          [-t type_of_multiresolution_transform]
**                  1: linear wavelet transform: a trous algorithm 
**                  2: bspline wavelet transform: a trous algorithm 
**                  3: wavelet transform in Fourier space 
**                  4: morphological median transform 
**                  5: morphological minmax transform 
**                  6: pyramidal linear wavelet transform 
**                  7: pyramidal bspline wavelet transform 
**                  8: pyramidal wavelet transform in Fourier space: 
**                     wavelet =  between two resolutions 
**                  9: pyramidal wavelet transform in Fourier space: 
**                     wavelet = difference  between the square 
**                                                of two resolutions
**                 10: pyramidal median transform 
**                 11: morphological pyramidal minmax transform 
**                 12: pyramidal laplacian 
**                 13: decomposition on scaling function 
**                 14: Mallat's wavelet transform 
**                 15: G transform (morphological min-max algorithm 
**                 16: Feauveau's wavelet transform 
**                 17: Haar's wavelet transform 
**                 18: Feauveau's wavelet transform without undersampling 
**
**                 default is 2
**                 17 is not yet implemented
**
**           [-p]
**                Poisson noise
**                default is gaussian noise
**
**           [-g sigma]
**                Gaussian noise
**                  sigma = noise standard deviation 
**                by default, the noise is gaussian, and the standard 
**                devaition is automatically estimated. 
**
**           [-c gain,sigma,mean]
**                case of a CCD: noise = Poisson noise + read-out noise
**                  gain = CCD gain 
**                  sigma = standard deviation of the read-out noise
**                  mean = mean of the read-out noise
**                if this option is set, 
**                           Noise = Poisson + Gaussian read-out Noise
**                it is generally the case with the CCD.
**                Attention, these parameters must be separated by a comma 
**                without space. example: -c 0.133,1.733,0.
**                If mean, or sigma and mean are omitted, default values are 0.
**                gain can not be omitted. 
**
**           [-n number_of_scales]
**                number of scales used in the multiresolution transform
**                default is 4
**
**           [-s NSigma]
**                Thresolding at NSigma * SigmaNoise at each scale
**                default is 3
**
**           [-k]
**                kill isolated pixels in the multiresolution support
**                if the PSF is large, then isolated pixels are certainly 
**                residual noise, or cosmic rays, or artefacts. If this
**                option is set, we suppress these pixels in the support
**                default is no.
**
**           [-l]
**                dilate the multiresolution support
**                if this option is set, each scale is dilated by using 
**                the mathematical morphology operator. Usefull if artefacts
**                remain arround objects. 
**                default is no
**
**           [-w support_file_name]
**                if this option is set, an image is created from the 
**                multiresolution support. Default is no.
**
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_support.cc 3.1 96/05/02CEA 1994 @(#)";
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Support.h"
#include "MR_Sigma.h"

char Name_Imag_In[256];  /* input file image */
char Name_MR_Out[256]; /* output file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
float Fwhm=0.5;                   /* Full width at half maximum */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_TRANSFORM; /* type of transform */
Bool Dil=False;                 /* dilate the support */
Bool SupIsol=False;             /* suppress isolated pixel in the support */
char Name_Write_Sup[256];          /* output support file name */
Bool WriteSup = False;          /* write the support on the disk */
Bool UseNSigma =False;

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose = False;
sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
int NbrUndec = -1;                     /*number of undecimated scale */

/*********************************************************************/

static void usage(char *argv[])
{
 
    fprintf(OUTMAN, "Usage: %s options in_image out_mr_file\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
 
    transform_usage(Transform);
    manline();
    sb_usage(SB_Filter);
    manline();
    nbr_nbr_undec_usage(NbrUndec);
    manline();
    gauss_usage();
    manline();  
     
    poisson_noise_usage();
    manline(); 
      
    ccd_usage();
    manline();    

    nbr_scale_usage(Nbr_Plan);
    manline();
    
    nsigma_usage(N_Sigma);
    manline();

    kill_isol_pix_usage();
    manline();
   
    dilate_sup_usage();
    manline();

    support_file_usage();
    manline();
    vm_usage();
    manline();
    verbose_usage();    
    manline();
    manline();
 
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void supinit(int argc, char *argv[])
{
    Bool OptL = False, Optf = False;
    int c;
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif       
    /* get options */
    while ((c = GetOpt(argc,argv,"u:pklt:T:Lg:c:n:s:w:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'u':
 		if (sscanf(OptArg,"%d",&NbrUndec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
                    
		}
  		break;
	   case 'v': Verbose = True; break;
	   case 't':
		/* -t <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            usage(argv);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            usage(argv);
 		}
		break;
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
           case 'p':
                /* Poisson noise */
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
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
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "bad noise parameter: %s\n", OptArg);
		    usage(argv);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE);
		    usage(argv);
		}
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
		    usage(argv);
		}
		UseNSigma =True;
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		break;
            case 'k':
               /* kill i */
               SupIsol = True;
               break;
            case 'l':
               /* dilate the support */
                Dil = True;
               break;
 	   case 'w':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s",Name_Write_Sup) != 1) 
                {
		   fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		   usage(argv);
		}
                WriteSup = True;
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

       if ((isotrop(Transform) != True) &&
           ((SupIsol == True) || (Dil == True)) )
       {
          fprintf(OUTMAN, "Error: option -k and -l are not valid with non isotropic transform. ...\n");
		usage(argv);
       }
	if ((Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_MALLAT) && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with Mallat transform ... \n");
           exit(0);
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
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}

/*********************************************************************/

int main(int argc, char *argv[]) 
{
    int i,k;
    Ifloat Dat;
    int Nl, Nc;
    fitsstruct Header;
    char Cmd[256];

    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    supinit(argc, argv);
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "Transform = " << StringTransform(Transform) << endl;
       if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
       {
          cout << StringSBFilter(SB_Filter) << endl;
          if (Norm == NORM_L2) cout << "L2 normalization" << endl;
       }
       cout << "Number of scales = " << Nbr_Plan << endl;
       if (Transform == TO_UNDECIMATED_MALLAT) cout << "Number of undecimated scales = " <<  NbrUndec << endl;
    }

    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Nl = Dat.nl();
    Nc = Dat.nc();
    Header.origin = Cmd;
    check_scale(Dat.nl(), Dat.nc(), Nbr_Plan);

    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MRNoiseModel ModelData;
    ModelData.alloc(Stat_Noise, Dat.nl(),Dat.nc() ,Nbr_Plan, Transform, PtrFAS, Norm, NbrUndec);
    int NbrBand = ModelData.nbr_band();
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (SupIsol == True) ModelData.SupIsol = True;
    if (Dil == True) ModelData.DilateSupport = True;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;
    
    if (UseNSigma  == True)
          for (i=0; i < NbrBand; i++) ModelData.NSigma[i]=N_Sigma;

    ModelData.model(Dat);

    /* write the support */
    ModelData.write_support_mr(Name_MR_Out);

    /* support image creation */
    if (WriteSup) ModelData.write_support_ima(Name_Write_Sup);
    exit(0);
} 

/*********************************************************************/
