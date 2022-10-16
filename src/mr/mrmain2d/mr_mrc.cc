/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.4
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/02
**    
**    File:  mr_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  deconvolution of an image  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
**
**    USAGE: mr_deconv option image psf output
**        where options = 
**          [-d type_of_deconvolution]
**               type_of_deconvolution = 
**                  1: Deconvolution by Van Cittert's algorithm
**                  2: Deconvolution by Gradient's algorithm
**                  3: Deconvolution by division in Fourier space
**                  4: Deconvolution by Lucy's algorithm
**                  5: Deconvolution by CLEAN algorithm
**                  6: Deconvolution by multiresolution Van Cittert's algorithm
**                  7: Deconvolution by multiresolution Gradient's algorithm
**                  8: Deconvolution by multiresolution Lucy's algorithm 
**                  9: Deconvolution by multiresolution CLEAN algorithm 
** 
**                  default is 8
**                  5 and 9 are not yet implemented
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
**           [-e Epsilon]
**                Convergence parameter
**                default is 0.0001
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
**           [-a Boolean_image]
**                add a Boolean image to the multiresolution support
**                Boolean_image = image file name
**                if we have information about positions of stars, ...,
**                we can add this position to the support
**                default is no addition
**
**           [-f Fwhm]
**                Fwhm = full width at half maximum
**                only used if type_of_deconvolution in [3,5,9]
**
**           [-w support_file_name]
**                if this option is set, an image is created from the 
**                multiresolution support. Default is no.
**
**           [-r residual_file_name]
**                if this option is set, the residual is written to 
**                the file of name residual_file_name. By default, the
**                residual is not written.
**
**           t,p,g,c,n,s,e,k,l,a,w options are only used if
**           type_of_deconvolution in [6..9]
**  
**
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Deconv.h"
#include "MR_Deconv.h"
#include "MR_Sigma.h"
#include "MR_Sigma.h"

char Name_Imag_In[256];  /* input file image */
char Name_Psf_In[256];   /* PSF */
char Name_Imag_Out[256]; /* output file name */
char Name_Imag_Start[256]; // First guess input solution 
char Name_Imag_ICF[256];   // ICF file name

Bool WriteResi = False;         /* write the residual */
Bool PsfMaxShift = True;        /* shift the max of the PSF to the center of the image */
char Name_Resi[256];  

Bool UseNSigma =False;
Bool GaussConv=False;
Bool UseICF=False;
Bool UseGuess=False;

float Fwhm = 0.;
float Converg = 1.;                      // convergence parameter 
int Max_Iter=DEFAULT_MAX_ITER_DECONV;    // Maximum number of iteration 

int Nbr_Plan=DEFAULT_NBR_SCALE;   /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = TO_PAVE_FFT; /* type of transform */
type_deconv Deconv = DEC_MR_CLEAN; /* type of deconvolution */
float Epsilon=1e-4;               /* convergence parameter */
  
Bool KillLastScale =False;
Bool UseEnergy = True;
Bool CleanLastScale = False;

char Name_RMSMap[256]; 
Bool UseRMSMap=False;
Bool KeepPositivSup =  True;
Bool AddResi = True;
int CleanFirstScale=0;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose;
Bool PositivIma = True;
 
float RegulParam = 0.01;

const int NBR_OK_NOISE = 3;
static type_noise TabNoise[NBR_OK_NOISE] = {NOISE_GAUSSIAN, NOISE_UNI_UNDEFINED, NOISE_CORREL};

/*********************************************************************/

static void usage(char *argv[])
{
    int i;
    fprintf(OUTMAN, "Usage: %s options in_image in_psf out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
   
    gauss_usage();
    manline();
         
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_OK_NOISE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise)TabNoise[i]));
    fprintf(OUTMAN, "             default is %s\n",  StringNoise(Stat_Noise));
    manline();
      
    nbr_scale_usage(Nbr_Plan);
    manline();
 
    nsigma_usage(N_Sigma);
    manline();
     
    max_iter_usage(Max_Iter);
    manline();
  
    converg_param_usage(Epsilon);
    manline();
    
    //rms_noise_usage();
    //manline();
    fprintf(OUTMAN, "         [-R RMS_Map_File_Name]\n");
    fprintf(OUTMAN, "              RMS Map (only used with -m 2 and -m 3 options).\n");
    manline();
            
    fprintf(OUTMAN, "         [-f ICF_Fwhm]\n");
    fprintf(OUTMAN, "              Intrinsic correlation function.\n");
    fprintf(OUTMAN, "              Fwhm = Full Width at Half Maximum.\n");
    manline();

    fprintf(OUTMAN, "         [-G LoopGain]\n");
    fprintf(OUTMAN, "              Loop Gain parameter. \n");
    fprintf(OUTMAN, "              Default is %f\n", RegulParam);
    manline();

    fprintf(OUTMAN, "         [-E]\n");
    fprintf(OUTMAN, "              Use Dirty Map instead of Energy Dirty Map. \n");
    manline();

    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "              Apply CLEAN also on the last scale. \n");
    fprintf(OUTMAN, "              Default is Van-Cittert.\n");
    manline();
    
    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "              Do not add the residual to the CLEAN map. \n");
    manline();
    
    fprintf(OUTMAN, "         [-F FirstScaleToUse]\n");
    fprintf(OUTMAN, "             Default is 1. \n");
    manline();
    
    write_residual_usage();
    manline();    
    kill_last_scale_usage();
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
static void decinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    

    /* get options */
    while ((c = GetOpt(argc,argv,"R:Km:d:t:g:n:s:i:e:r:f:vzZ:G:ELAF:")) != -1) 
    {
	switch (c) 
        {
	  case 'A': AddResi=False; break;
	  case 'E': UseEnergy=False; break;
 	  case 'L': CleanLastScale=True; break;
 	  case 'R':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseRMSMap = True;
 		break;
            case 'K': KillLastScale = True; break;
            case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_OK_NOISE)) 
                                          Stat_Noise = (type_noise) (TabNoise[c-1]);
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
		    fprintf(OUTMAN, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	    case 'F':
 		if (sscanf(OptArg,"%d",&CleanFirstScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad first scale: %s\n", OptArg);
		    exit(-1);
		}
                CleanFirstScale --;
		if (CleanFirstScale < 0)
		{
		   fprintf(OUTMAN, "Error: bad first scale: %d\n", CleanFirstScale+1);
                   exit(-1);
		}
 		break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
 		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
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
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		UseNSigma =True;
                break;
	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0) Max_Iter = DEFAULT_MAX_ITER_DECONV;
		break;
	   case 'e':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Epsilon < 0) || (Epsilon > 1.)) 
                           Epsilon = DEFAULT_EPSILON_DECONV;
		break;
           case 'v': Verbose = True;break;
	   case 'r':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Resi) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteResi = True;
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
 	  case 'f':
		/* -f < Fwhm parameter> */
		if (sscanf(OptArg,"%f",&Fwhm) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Fwhm: %s\n", OptArg);
		   exit(-1);
		} 
		if (Fwhm > 0) GaussConv = True;
		if (UseICF == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
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
	if (CleanFirstScale >= Nbr_Plan)
	{
	    fprintf(OUTMAN, "Error: bad first scale ...\n");
	    fprintf(OUTMAN, "       FirstScale must be small than the number of scales.\n");
	    exit(-1);
	}
		
 	if (UseRMSMap == True)
 	{
 	   if ((Stat_Noise != NOISE_NON_UNI_ADD) && (Stat_Noise !=  NOISE_CORREL))
 	   {
 	      cerr << "Error: this noise model is not correct when RMS map option is set." << endl;
 	      cerr << "       Valid model is: " << endl;
 	      cerr << "        " << StringNoise(NOISE_CORREL) << endl;
	      exit(-1);
   	   }
  	}
 	   
	if ((Stat_Noise == NOISE_CORREL) && (UseRMSMap != True))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: this noise model need a noise map (-R option) " << endl;
           exit(-1);
  	}
	
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Psf_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    Ifloat DataB;
    int b, k;
    Ifloat Result;
    Ifloat Resi;
    fitsstruct Header;
    char Cmd[256];
    Ifloat Guess, Ima_ICF;
   
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    decinit(argc, argv);

if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    if (Stat_Noise == NOISE_GAUSSIAN)
    {
        cout << "Type of Noise = GAUSSIAN" << endl;
	if (Noise_Ima > 0) cout << "Sigma Noise = " << Noise_Ima << endl;  
    }
 
    cout << "Deconv = " << StringDeconv(Deconv) << endl;
    cout << "N_Sigma = " << N_Sigma << endl;
    cout << "Epsilon = " << Epsilon << endl;
    cout << "Max_Iter = " << Max_Iter << endl;
    cout << "Convergence paramter = " << Converg << endl;
    if (KillLastScale == True) cout << "Kill last scale " << endl;
    cout << "Fwhm = " << Fwhm << endl;
    if (WriteResi == True)
      cout << "Image file name : " << Name_Resi << endl;
}

    // read input image
    MRDeconv CDec;
    io_read_ima_float(Name_Imag_In, DataB, &Header);
    int Nl = DataB.nl();
    int Nc = DataB.nc();
    int Nl1,Nc1;
    Bool ModifSize=False;
    dec_line_column (Nl, Nl1);
    dec_line_column (Nc, Nc1);
    if (Nl1 < Nc1) Nl1 = Nc1;
    else Nc1 = Nl1;
    if ((Nl != Nl1) || (Nc != Nc1)) ModifSize = True;
    if (ModifSize)
    {
        CDec.Imag.alloc(Nl1, Nc1, "Ptr Imag");
        im_extend (DataB, CDec.Imag);
    }
    else CDec.Imag = DataB;
    
    io_read_ima_float(Name_Psf_In, CDec.Psf);
    Header.origin = Cmd; 
    CDec.UseMRCEnergy = UseEnergy;
    CDec.CleanLastScale = CleanLastScale; 
    CDec.KillLastScale = KillLastScale;
    CDec.PositivConstraint = PositivIma;
    CDec.DecMethod = Deconv;
    CDec.Noise_Ima = Noise_Ima;
    CDec.MaxIter = Max_Iter;
    CDec.EpsCvg = Epsilon;
    CDec.RegulParam = RegulParam;
    CDec.GaussConv = GaussConv;
    CDec.Fwhm = Fwhm;
    CDec.Verbose = Verbose;
    CDec.CleanFirstScale=CleanFirstScale;
    Ifloat *Pt_G = NULL;
    if (UseGuess == True) Pt_G = &Guess;
    Ifloat *Pt_ICF = NULL;
    if (UseICF == True) Pt_ICF = &Ima_ICF;
    
    //DECONVOLUTION
    if (Verbose == TRUE) cout << " Start the deconvolution ... " << endl;
    CDec.StatNoise = Stat_Noise;

    // noise model class initialization
    MRNoiseModel ModelData(Stat_Noise, CDec.Imag.nl(), CDec.Imag.nc(), 
                           Nbr_Plan, Transform);
    int NbrBand = ModelData.nbr_band();
    ModelData.OnlyPositivDetect = KeepPositivSup;
    if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
          for (b=0; b < NbrBand; b++) ModelData.NSigma[b]=N_Sigma;
    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }
    CDec.ModelData = &ModelData;    
    
    CDec.im_deconv(Pt_G, Pt_ICF);
    if (AddResi == True) CDec.Obj += CDec.Resi;
    
    // Write the results
    if (ModifSize == True) 
    {
       im_extract (CDec.Obj, DataB);
       io_write_ima_float(Name_Imag_Out, DataB, &Header);
    }
    else io_write_ima_float(Name_Imag_Out, CDec.Obj, &Header);
    if (WriteResi == True)
    {
       if (ModifSize == True) 
       {
           im_extract (CDec.Resi, DataB);
          io_write_ima_float(Name_Resi, DataB, &Header);
       }
       else  io_write_ima_float(Name_Resi, CDec.Resi, &Header);
    }   
    exit(0);
} 

