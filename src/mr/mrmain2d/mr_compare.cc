/******************************************************************************
**                   Copyright (C) 1995 CEA
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
**    File:  mr_compare.cc
**
*******************************************************************************
**
**    DESCRIPTION   Computes the correlation coefficients and the signal to
**    -----------   noise ratio between two images.  The comparison is done 
**                  in the wavelet space. 
**                  At each scale we compute the SNR and the
**                  correlation. We get two curves, 
**                        . (coeff_correlation - frequence)
**                        . (SNR - frequence)
**
**    PARAMETRES    
**    ----------
**         image 1
**
**         image 2 
**
**         Number of scales of the comparison 
**
**         Comparison parameter N_Sigma 
**            if N_Sigma > 0, the comparison is not done with
**            the full image, but only with the structure in the
**            wavelet space > N_Sigma * sigma(scale)
**
**         Input-Output correlation table name 
**
**         Input-Output signal to noise ratio table name
**
**
*******************************************************************************
**
**    PARAMETRES
**    ----------
**
**        Input_File: reference image
**
**        Image_1:  image to be compared to reference image
**                       
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
**
**           [-s nsigma ] : a thresholding is done at each scale at 
**                        nsigma * Sigma_Noise
**                        by default, nsigma = 3
**
**           [-n number of scales] : number of scales used by the 
**                                 multiresolution transform
**                                 by default number_of_scales = 4
**
**      
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_Support.h"
#include "MR_Sigma.h"

char File_Name_Imag[80];  /* input file image */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = DEFAULT_TRANSFORM; /* type of transform */

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

#define MAX_IMA 50
// #define MAX_BAND 30

char Tab_Name [MAX_IMA][80]; /* image file name array */
int ImageNbr = 0;

void print_im_cmp(FILE *FileDes, Ifloat &Imag_1,Ifloat &Imag_2, Ifloat &ImagErr,
                  int Ima=0, Bool UseSupport=False, int Scale=-1, 
		  Bool CmpWave = True);

float TabRMS[MAX_IMA][MAX_BAND+1];
float TabCorrel[MAX_IMA][MAX_BAND+1];
float TabSNR[MAX_IMA][MAX_BAND+1];
float TabSNRb[MAX_IMA][MAX_BAND+1];
float TabDiff[MAX_IMA][MAX_BAND+1];
float TabMin[MAX_IMA][MAX_BAND+1];
float TabMax[MAX_IMA][MAX_BAND+1];
float TabPSNR[MAX_IMA];

MRNoiseModel ModelData;
Bool Verbose = False;

/*********************************************************************/

static void usage(char *argv[])
{
 
    fprintf(OUTMAN, "Usage: %s options ref_image ima1 [ima2, [ima3, ...]]\n", argv[0]);
    fprintf(OUTMAN, "   where options = %s\n", argv[0]);

    isotrop_transform_usage(Transform);
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
    vm_usage();
    manline();
    verbose_usage();
    manline(); 
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void hcinit(int argc, char *argv[])
{
    int c;
    Bool OptMR=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"s:n:c::pg:t:vzZ:")) != -1) 
    {
	switch (c) 
        {
 	   case 'v': Verbose = True; break;
	   case 't':
		/* -t <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_ISOTROP_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
 		OptMR=True;
		break;
            case 'p':
                /* Poisson noise */
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
                OptMR=True;
               break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
                OptMR=True;
		break;
             case 'c':
		/* -c <gain sigma mean> */
 		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
                OptMR=True;
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan < 0) || (Nbr_Plan > MAX_BAND)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_BAND);
		    exit(-1);
		}
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma < 0.)  N_Sigma = DEFAULT_N_SIGMA;
                OptMR=True;
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
	if (OptInd < argc) 
             strcpy(File_Name_Imag, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) 
             strcpy(Tab_Name[0], argv[OptInd++]);
        else usage(argv);

        ImageNbr = 1;

	/*read the other image file names */
        while ((ImageNbr < MAX_IMA) && (OptInd < argc))
        {
             strcpy(Tab_Name[ImageNbr++], argv[OptInd++]);
        }
        
        if ((Nbr_Plan < 2) && (OptMR == True))
        {
           fprintf(OUTMAN, "Error: if Nbr_Plan < 2, t,p,g,c,s options are not valid ...\n");
           exit(-1);
        }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***************************************************/

void print_info_pict (FILE *File_Des, Ifloat &Pict, char *Name_Pict)
{
    float Min,Max;
    float Moy,Sigma;

    Sigma = sigma(Pict);
    Moy = average(Pict);
    Min = min(Pict);
    Max = max(Pict);

    fprintf(File_Des, "\n  Name : %s (%d,%d) \n",Name_Pict,Pict.nl(),Pict.nc());
    fprintf (File_Des, "     Minimum = %f, Maximum = %f\n",Min,Max);
    fprintf (File_Des, "     Mean    = %f, Sigma   = %f\n", Moy,Sigma);
}

/****************************************************************************/

void print_im_cmp(FILE *FileDes, Ifloat &Imag_1,Ifloat &Imag_2, 
                  Ifloat &ImagErr, int Ima, Bool UseSupport, 
                  int Scale, Bool WaveCmp)
/* computes information between two images

   Snr = variance (Imag_1) / variance de (Imag_1 - Imag_2)
   Snr_Db = 10 log10 (Snr)
   Correl = correlation des images Imag_1 et Imag_2
*/
{
    int i,j,s;
    int Nl = Imag_1.nl();
    int Nc = Imag_2.nc();
    float Snr,S1,S2,Pct;
    float Error;
    float Sum_X2= 0.,Sum_Y2= 0.,Sum_XY= 0.;
    float Moy = 0., Ecart = 0.;
    float Moy_Err = 0., Ecart_Err = 0.;
    float Moy_Abs_Err = 0., Ecart_Abs_Err = 0.;
    int Size = 0;
    float Snr_Db, Correl;

    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        if ( (UseSupport == False) || 
	     ((WaveCmp == True) && (ModelData(Scale,i,j) == True)))
	{
           /* Calcul Correlation */
           Sum_X2 += Imag_1 (i,j) * Imag_1 (i,j);
           Sum_Y2 += Imag_2 (i,j) * Imag_2 (i,j);
           Sum_XY += Imag_1 (i,j) * Imag_2 (i,j);

           /* Calcul of the standard deviation of Imag_1 */
           Moy += Imag_1 (i,j);
           Ecart += Imag_1 (i,j) * Imag_1 (i,j);

           /* Calcul of the standard deviation of the error */
           Error = ImagErr(i,j);
           Moy_Err += Error;
           Ecart_Err += Error * Error;

           Error = ABS(Error);
           Moy_Abs_Err += Error;
           Ecart_Abs_Err += Error*Error;
           Size ++;
	 }
    }
    Pct = (float) Size / (float)(Nl*Nc) * 100.;

    /* Correlation */
    Correl = Sum_XY / sqrt (Sum_X2 * Sum_Y2);

    /* variance of Imag_1 */
    Moy /= (float) Size;
    Ecart /= (float) Size;
    S1 = Ecart - Moy*Moy;

    /* variance of the error */
    Moy_Err /= (float) Size;
    Ecart_Err /= (float) Size;
    S2 = Ecart_Err - Moy_Err*Moy_Err;

    /* variance of the absolute error */
    Moy_Abs_Err /= (float) Size;
    Ecart_Abs_Err /= (float) Size;
    Ecart_Abs_Err = Ecart_Abs_Err - Moy_Abs_Err*Moy_Abs_Err;

    /* signal to noise ration */
    Snr = S1 / S2;

    /* signal to noise ration (dB) */
    Snr_Db = 10. * log10 (Snr);

    s = Scale+1;
    TabRMS[Ima][s] = sqrt(S2);
    TabCorrel[Ima][s] = Correl;
    TabSNR[Ima][s] = Snr;
    TabSNRb[Ima][s] = Snr_Db;
    TabDiff[Ima][s] = sqrt(Ecart_Abs_Err);
    TabMin[Ima][s] = min(Imag_1-Imag_2);
    TabMax[Ima][s] = max(Imag_1-Imag_2);
    if (Scale < 0)
    {
       if (TabRMS[Ima][s] > 0.) 
          TabPSNR[Ima] = 10. * log10(255.*255./ (TabRMS[Ima][s]*TabRMS[Ima][s]));
    }
    fprintf (FileDes, "\n");
   if (UseSupport == True)
     fprintf (FileDes, "     Signif = %5.2f %%\n", Pct);
   fprintf (FileDes, "     Correlation = %f\n", Correl);
   fprintf (FileDes, "     Sigma ( ABS(error) ) =  %f\n", TabDiff[Ima][s]);
   fprintf (FileDes, "     RMS =  %f\n", TabRMS[Ima][s]);
   fprintf (FileDes, "     Min(Error) =  %f\n", TabMin[Ima][s]);
   fprintf (FileDes, "     Max(Error) =  %f\n", TabMax[Ima][s]);
   fprintf (FileDes, "     SNR = variance(Ref) / variance(Error) = %f\n", Snr);
   fprintf (FileDes, "     SNRb = 10 log10(SNR) = %f dB\n", Snr_Db);
   if (Scale < 0) 
    if (TabRMS[Ima][s] > 0.)  
      fprintf (FileDes, "     PSNR = 10 log10(255^2/MSE) = %f dB\n\n\n",  TabPSNR[Ima]);
    else fprintf (FileDes, "     PSNR = Infinite (images are identical)\n\n\n");
}

/****************************************************************************/

void print_result(float Tab[MAX_IMA][MAX_BAND+1])
{
   char name[100];
   int i,j;
   
   fprintf (stdout,"            ");
     
   for (i = 0; i < ImageNbr; i++)
       fprintf (stdout,"%s   ", Tab_Name[i]);
   fprintf (stdout,"\n\n");  
   for (j = 0; j <= Nbr_Plan; j++)
   {
      if (j == 0) sprintf(name,"image  "); 
      else sprintf(name,"scale %d ", j);
      fprintf (stdout,"%s: ",name); 
      for (i = 0; i < ImageNbr; i++)
             fprintf (stdout,"%f ", Tab[i][j]);
       fprintf (stdout,"\n\n");
    }

}

/****************************************************************************/

int main (int argc, char *argv[])
{
    int i,s;
    int Nl,Nc;
    Ifloat Dat;
    Ifloat Imag_2;
    Bool UseWavelet=True;

    /* Get command line arguments, 
       open input file(s) if necessary */
    lm_check(LIC_MR1);
    hcinit(argc, argv);

    if (Nbr_Plan < 2)
    {
       UseWavelet = False;
       Nbr_Plan=0;
    }
    io_read_ima_float(File_Name_Imag, Dat);
    Nl = Dat.nl();
    Nc = Dat.nc();

    /* Variables declarations */
   MultiResol MR_Data (Nl, Nc, Nbr_Plan, Transform, "MR_Data");
   int NbrBand = MR_Data.nbr_band();
   Ifloat ImagGauss (Nl, Nc, "ImagGauss");
   Ifloat ImagErr (Nl, Nc, "ImagErr");
   
    
    if ((UseWavelet == True) && (N_Sigma > FLOAT_EPSILON))
    {
       ModelData.alloc(Stat_Noise,Nl,Nc,Nbr_Plan,Transform);
       if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
       if (N_Sigma > 0) for (i=0; i < NbrBand; i++) 
                             ModelData.NSigma[i]=N_Sigma;
       ModelData.CCD_Gain = PasCodeur;
       ModelData.CCD_ReadOutSigma = SigmaGauss;
       ModelData.CCD_ReadOutMean = MeanGauss;       
       ModelData.model(Dat, MR_Data);
       // if TransImag == true MR_Data contains the wavelet coefficient
       // of the transform image, and not from the image
       // then we compute the WT of the image
       if (ModelData.TransImag == True) MR_Data.transform (Dat, I_ZERO);
    }
    else if (UseWavelet == True) MR_Data.transform (Dat);
     

    /* MultiResol obj declaration for the image to compare */
    MultiResol MR_Comp (Nl, Nc, Nbr_Plan, Transform, "MR_Support");
    MultiResol MR_Err (Nl, Nc, Nbr_Plan, Transform, "MR_Err");

    /* on met une entete */
    fprintf (stdout,"   COMPARE IMAGES TO %s\n\n", File_Name_Imag);

    /* Information on the reference image */
    print_info_pict (stdout, Dat, File_Name_Imag);

    /* for each image to compare */
    for (i = 0; i < ImageNbr; i++)
    {
        io_read_ima_float(Tab_Name[i], Imag_2);

        /* Information on the image */
        print_info_pict (stdout, Imag_2, Tab_Name[i]);

        ImagErr = Dat - Imag_2;
        print_im_cmp (stdout, Dat, Imag_2, ImagErr, i, False);

        if (UseWavelet == True)
        {
            MR_Comp.transform (Imag_2);
            MR_Err.transform (ImagErr);

            /* compare scale by scale de tous les plans */
            for (s = 0; s < NbrBand-1; s++)
            {
                fprintf (stdout,"   band %d: \n", s+1);
                if (N_Sigma > FLOAT_EPSILON)
                       print_im_cmp (stdout, MR_Data.band(s), 
                                             MR_Comp.band(s), MR_Err.band(s), 
                                             i, True, s, True);     
                else print_im_cmp (stdout, MR_Data.band(s), 
                                           MR_Comp.band(s), MR_Err.band(s),
                                           i, False, s, True);  
            }
            s = NbrBand-1;
            fprintf (stdout,"   Band %d: \n", s+1);
            print_im_cmp (stdout, MR_Data.band(s), MR_Comp.band(s),
                                  MR_Err.band(s), i,False,s, False);
            fprintf (stdout,"\n================================\n\n");   
        }    
    }

    if ((UseWavelet == True) && (ImageNbr > 1))
    {
       fprintf (stdout,"\n\n================================\n\n");
       fprintf (stdout,"       RMS  comparison \n\n");
       print_result(TabRMS);

       fprintf (stdout,"\n\n================================\n\n");
       fprintf (stdout,"       SNR  comparison \n\n");   
       print_result(TabSNR );

       fprintf (stdout,"\n\n================================\n\n");   
       fprintf (stdout,"       SNR  (dB) comparison \n\n");
       print_result(TabSNRb);
        
       fprintf (stdout,"\n\n================================\n\n"); 
       fprintf (stdout,"       Correlation comparison \n\n");
       print_result(TabCorrel);
 
       fprintf (stdout,"\n\n================================\n\n"); 
       fprintf (stdout,"       Standard deviation of the difference \n\n"); 
 
       print_result(TabDiff);
    }
/*
print_result(TabMin);
print_result(TabMax);

*/
    exit(0);
}

/***************************************************/


