/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.12
**
**    Author: Jean-Luc Starck
**
**    Date:  98/05/20
**    
**    File:  mr_detect.cc
**
*******************************************************************************
**
**    DESCRIPTION  detect objects
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_detect option image output
**        where options = 
**
**           [-p]
**                Poisson noise
**                default is gaussian noise
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
**                normal unit Gain = e-/DN
**                normal unit RON = e-
**                 we normalize RON by RON_n = RON / Gain = unit (DN)
**                 
**           [-n number_of_scales]
**                number of scales used in the detect estimation
**
**
** file xx.mes:
**   for each object:
**         scale_number object_number 
**         meanx meany sigmax sigmay
**         Angle Max Flux Magnitude
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "MR_Filter.h"
#include "IM_VisTool.h"
#include "MR_VisTree.h"
#include "MR_ListObj.h"
#include "MR_Abaque.h"
#include "MR_MVM.h"

char Name_Imag_In[256];  /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Psf[256];      /* PSF file name */
char Name_BGR[256];          /* background file name */
int Nbr_Plan=5;           /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform =  TO_PAVE_BSPLINE; /* type of transform */
int Max_Iter = 1; /* Maximum number of iteration */
Bool WriteBGR = False;          /* write the detect on the disk */
float EpsilonPoisson = DEFAULT_EPSILON;
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
int FirstScale = DEF_FIRST_DETECT_SCALE;
int NPix = 16;
int LastScale=-1;
int MinEvent = DEF_MINEVENT_NUMBER;

Bool UseNSigma =False;
Bool WriteAllObj = False;
Bool WriteSimu = False;
Bool SubSegment = False;
Bool Verbose = False;
char Name_RMSMap[256]; 
Bool UseRMSMap=False;
Bool KillIsolObj=True;
Bool KeepPositivSup = True;

type_objrec_method ReconsMethod = DEF_RECOBJ_METHOD; 
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);
   
#define NB_ITER_REC_OBJ 5

int fin_rec = -1;
int deb_rec = 1;
int Nb_iter_rec = DEF_MAX_ITER_RECOBJ;
int NumObj = 0;
float Eps_ErrorRec = DEF_ERROR_RECOBJ; // parameter for convergence
float Fwhm = -1;
float ZoomPSF = 1;

Bool KeepBorderObj = False;
Bool Crea_RaDec_Table = False;
t_order Order = O_Obj;
float FluxMult = 1.;
Bool GetRmsCoeffGauss = True;

#define MVM_Max   1
#define MVM_Align 2

int TypeMVM = 1;
float DistMaxAl=1;

const int NBR_OK_TRANS = 5;

static type_transform TabTrans[NBR_OK_TRANS] = {TO_PAVE_BSPLINE,TO_SEMI_PYR,TO_PYR_BSPLINE,TM_TO_PYR,TM_TO_SEMI_PYR};

Bool InfoSubObj = False;
Bool UsePSF = False;
Bool Deconv = False;

t_mvm VisionModel=MVM_RUEBIJ;

Bool AperPhot = False;
float KSigmaAper = 3.;
int BGRSize = 16;
t_bgr BgrMethod = BGR_MR;   // BGR_IMA, BGR_MR, BGR_VALUE
float BackgroundValue = 0.;
Bool MVMDebug = False;
Bool WriteLarge=False;

Bool WriteObjFullSize = False;
Bool AnisotropAnalysis = False;
float AnisotropParam = 2.;
type_threshold TypeThreshold = DEF_THRESHOLD;

/*********************************************************************/

static void usage(char *argv[])
{
    int i;
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options are =  \n");
 
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_OK_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransform((type_transform)TabTrans[i]));
    fprintf(OUTMAN, "             default is %s\n",  StringTransform((type_transform)Transform));
    manline();

    fprintf(OUTMAN, "         [-V Multiscale_Vision_Model]\n");   
    for (i = 0; i < NBR_MVM; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMVM ((t_mvm)i));
    fprintf(OUTMAN, "             default is %s\n",StringMVM ((t_mvm) VisionModel));
    manline();
     
    nbr_scale_usage(Nbr_Plan);
    manline();  

    noise_usage();
    manline();  
  
    gauss_usage();
    manline();
    
    ccd_usage();
    manline();

    nsigma_usage(N_Sigma);
    manline(); 
  
    prec_eps_poisson_usage(DEFAULT_EPSILON);
    manline();

    min_event_usage(MinEvent);   
    manline();
        
    size_block_usage(SizeBlock);
    manline();
    
    sigma_clip_block_usage(NiterClip);
    manline();  

    first_detect_scale_usage();
    manline();

    rms_noise_usage();
    manline();    
   
    last_detect_scale_usage();
    manline();   
    
    mrdetect_option_usage(Nb_iter_rec,Eps_ErrorRec);
        
    fprintf(OUTMAN, "         [-o] \n");
    fprintf(OUTMAN, "             Sub-object analysis. Default is no. \n");
    manline();
    
    fprintf(OUTMAN, "         [-D] \n");
    fprintf(OUTMAN, "             Perform a deconvolution. Default is no. \n");
    manline();
    
    fprintf(OUTMAN, "         [-P PsfFileName] \n");
    fprintf(OUTMAN, "             PSF file name.\n");
    manline();
    
    fprintf(OUTMAN, "         [-f Fwhm] \n");
    fprintf(OUTMAN, "             Full Width at Half Maximum.\n");
    manline();
       
    fprintf(OUTMAN, "         [-O PSF_Sampling] \n");
    fprintf(OUTMAN, "             PSF over-sampling value.\n");
    manline();

    fprintf(OUTMAN, "         [-a BgrMethod] \n");
    fprintf(OUTMAN, "             Aperture photometry.\n");
    for (int i = 0; i <  NBR_BGR; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                          StringBgr((t_bgr) i));
    fprintf(OUTMAN, "              default is no aperture photometry.\n");
    manline();
        
    fprintf(OUTMAN, "         [-B BgrFileName] \n");
    fprintf(OUTMAN, "             Background image file name.\n");   
    manline();
    
    fprintf(OUTMAN, "         [-G BgrValue] \n");
    fprintf(OUTMAN, "             Constant background value.\n");  
    fprintf(OUTMAN, "             default is: %f.\n", BackgroundValue);
    manline();
    
    fprintf(OUTMAN, "         [-b BGR_Size] \n");
    fprintf(OUTMAN, "             Background image size for automatic background esitmation.\n");
    fprintf(OUTMAN, "             default is: %d.\n", BGRSize);
    manline();
    
    fprintf(OUTMAN, "         [-l KSigmaAperture] \n");
    fprintf(OUTMAN, "             Aperture photometry size parameter.\n");
    fprintf(OUTMAN, "             default is: %f.\n", KSigmaAper);
    manline();    
        
//    fprintf(OUTMAN, "\n         [-I Processus_Iteration_Number]\n");
//    fprintf(OUTMAN, "              Number of times the detection is repeated.\n");
//    fprintf(OUTMAN, "               by default, it is done only once.\n");
    
    obj_rec_methode_usage(ReconsMethod);
    manline();

    tab_order_usage(Order);
    manline();
    vm_usage();
    manline(); 
    verbose_usage();
    manline();
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
   Bool UseMethod = False;
   Bool EpsOpt = False;
   Bool OptB = False;
   Bool Optb = False;
   Bool OptL = False;
   Bool Opta = False;
     
    /* get options */
    while ((c =
    GetOpt(argc,argv,"T:I:XV:a:B:G:b:l:P:Dt:g:m:c:n:s:i:E:N:S:F:L:w:u:UM:e:R:vkpKC:A:rqd:ozZ:f:O:")) != -1) 
    {
	switch (c) 
        {
           case 'T':
		/* -f <type> type of filtering */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
	            exit(-1);
 		}
                if ((c > 0) && (c <= NBR_THRESHOLD)) TypeThreshold = (type_threshold) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of detection: %s\n", OptArg);
	            exit(-1);
 		}
 		break;
            case 'I':  
	            if (sscanf(OptArg,"%f", & AnisotropParam) != 1) 
                    {
		       fprintf(OUTMAN, "bad anisotrop param value: %s\n", OptArg);
	               exit(-1);
		    }
	            AnisotropAnalysis = (AnisotropAnalysis == True) ? False : True;
	            break;
	   case 'X':
                MVMDebug = (MVMDebug == True) ? False : True;
                break;
	   case 'G': 
	        if (sscanf(OptArg,"%f", &BackgroundValue) != 1) 
                {
		    fprintf(OUTMAN, "bad background value: %s\n", OptArg);
	            exit(-1);
		}
		AperPhot = True;
		OptL = True;
		break;
	   case 'B': 
	        if (sscanf(OptArg,"%s", Name_BGR) != 1) 
                {
		    fprintf(OUTMAN, "bad background file name: %s\n", OptArg);
	            exit(-1);
		}
		AperPhot = True;
		OptB = True;
		break;
	   case 'a': 
 	         if (sscanf(OptArg,"%d",&c) != 1) 
                 {
		    fprintf(OUTMAN, "Error: bad background parameter: %s\n", OptArg);
		    exit(-1); 
		 }           
		 if ((c > 0) && (c <= NBR_BGR)) {
		    BgrMethod = (t_bgr) (c-1);
		    AperPhot = True; 
		    Opta = True; 
                 } else if (c == 0) {
		 
		 } else  {
		    fprintf(OUTMAN, "Error: bad background parameter: %s\n", OptArg);
	            exit(-1);
 		 }
		break;
	   case 'b': 
	        if (sscanf(OptArg,"%d",&BGRSize) != 1) 
                {
		    fprintf(OUTMAN, "bad background image size: %s\n", OptArg);
	            exit(-1);
		}
		AperPhot = True;
		Optb = True;
		break;
	   case 'l':
	        if (sscanf(OptArg,"%f",&KSigmaAper) != 1) 
                {
		    fprintf(OUTMAN, "bad Ksigma value: %s\n", OptArg);
	            exit(-1);
		}
		AperPhot = True;
		break;
	   case 'V': 
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad vision model: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <=  NBR_MVM))  VisionModel = (t_mvm)(c-1);
                else  
                {
		    fprintf(OUTMAN, "bad vision model: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 'f':  
	        if (sscanf(OptArg,"%f", &Fwhm) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Fwhm: %s\n", OptArg);
		   exit(-1);
		}
		if (Fwhm <= 0)
		{
		   fprintf(OUTMAN, "Error: bad Fwhm ... (Fwhm > 0) \n");
		   exit(-1);
		}
		Deconv = True;
 		break; 
	   case 'D': Deconv = True; break;
	   case 'P':		
 		if (sscanf(OptArg,"%s", Name_Psf) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UsePSF = True;
 		break;
	   case 'o': InfoSubObj = True; break;
           case 'q': TypeMVM = MVM_Align; break;
           case 'd':
                if (sscanf(OptArg,"%f",&DistMaxAl) != 1) 
                {
		    fprintf(OUTMAN, "bad maximum distance: %s\n", OptArg);
	            exit(-1);
                    
		}
		if (DistMaxAl < 0.)
		{
		    fprintf(OUTMAN, "bad maximum distance: %s\n", OptArg);
	            exit(-1);
                    
		}
		TypeMVM = MVM_Align; 
  		break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_OK_TRANS)) 
                                        Transform = (type_transform) (TabTrans[c-1]);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	case 'O':
		if (sscanf(OptArg,"%f",&ZoomPSF) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad PSF Sampling: %s\n", OptArg);
		    exit(-1);
		}
		if (ZoomPSF < 1)
		{
		   fprintf(OUTMAN, "Error: bad PSF sampling ... \n");
		   exit(-1);
		}
		break;       
 	 case 'E':
		if (sscanf(OptArg,"%f",&EpsilonPoisson) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    exit(-1);
		}
		
		if ((EpsilonPoisson <MIN_EPSILON) || (EpsilonPoisson > MAX_EPSILON)) 
                {
		    fprintf(OUTMAN, "Error: bad precision number: %s\n", OptArg);
		    fprintf(OUTMAN, "       %f <= Precision <= %f\n",MIN_EPSILON, MAX_EPSILON);
 		    exit(-1);
		}
		EpsOpt = True;
		break; 
 	 case 'e':
		if (sscanf(OptArg,"%d",&MinEvent) != 1) 
                {
		    fprintf(OUTMAN, "Min number of events: %s\n", OptArg);
		    exit(-1);
		}
                if (MinEvent < 1)  
                {
		    fprintf(OUTMAN, "Error: bad number of events min: %s\n", OptArg);
		    exit(-1);
		}
		break;
	case 'u':
		if (sscanf(OptArg,"%f",&Eps_ErrorRec) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad error reconstruction value: %s\n", OptArg);
		    exit(-1);
		}
                break;
	 case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad size block: %s\n", OptArg);
		    exit(-1);
		}
		break; 
	 case 'N':
		if (sscanf(OptArg,"%d",&NiterClip) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad it. number for the 3sigma clipping: %s\n", OptArg);
		    exit(-1);
		}
		break; 
	 case 'F':
		if (sscanf(OptArg,"%d", &FirstScale) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad first detection scale number: %s\n", OptArg);
		   exit(-1);
		}
		deb_rec = FirstScale;
		FirstScale --;
		break;    
	 case 'L':
		if (sscanf(OptArg,"%d", &LastScale) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad last detection scale number: %s\n", OptArg);
		   exit(-1);
		}
		fin_rec = LastScale;
		LastScale --;
 		break;  		    
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                UseNSigma =True;
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
		break;
//  	   case 'I':
// 		/* -i < Number of iterations> */
// 		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
//                 {
// 		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
// 				exit(-1);
// 		}
//                 if (Max_Iter <= 0)  Max_Iter = 1;
// 		break; 
  	   case 'i':
		/* -i < Number of iterations per object reconstruction> */
		if (sscanf(OptArg,"%d",&Nb_iter_rec) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Nb_iter_rec: %s\n", OptArg);
				exit(-1);
		}
                if (Nb_iter_rec <= 0) Nb_iter_rec  = NB_ITER_REC_OBJ;
		break;        
           case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_NOISE)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
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
             case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &PasCodeur,
                                              &SigmaGauss, &MeanGauss) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
		break;
	   case 'n':
		/* -n <Number_of_pixels> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
		break;
	   case 'w': 
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad  writing  paramter: %s\n", OptArg);
		    exit(-1);
		}
		switch(c)
		{
		   case 1: WriteAllObj = True; break;
		   case 2: WriteSimu = True; break;
		   case 3: WriteAllObj = True; 
                           WriteSimu = True;
                           break;
                   case 4: WriteAllObj = WriteObjFullSize = True; break;
		   default: 
		    fprintf(OUTMAN, "Error: bad writing parameter: %s\n", OptArg);
		    exit(-1);
		    break;
		}
		break;
 	   case 'R':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_RMSMap) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseRMSMap = True;
 		break;
	   case 'U': SubSegment = True; break;
	   case 'M':
		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad reconstruction method: %s\n", OptArg);
		    exit(-1);
		}
                if ((c > 0) && (c <= NBR_RECOBJ_METHOD+1)) 
                                  ReconsMethod = (type_objrec_method) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad reconstruction method: %s\n", OptArg);
	            exit(-1);
 		}
		UseMethod = True;
		break;     
	   case 'v': Verbose = True; break;
	   case 'r': GetRmsCoeffGauss = False; break;
	   case 'k': KillIsolObj = (KillIsolObj == True) ? False: True; break;
	   case 'K': KeepBorderObj = (KeepBorderObj == True) ? False: True; break;
	   case 'C': 
	         if (sscanf(OptArg,"%d",&c) != 1) 
                 {
		    fprintf(OUTMAN, "Error: bad order parameter: %s\n", OptArg);
		    exit(-1); 
		 }           
		 if ((c > 0) && (c <= NBR_ORDER)) 
                                    Order = (t_order) (c-1);
                 else  
                 {
		    fprintf(OUTMAN, "Error: bad order parameter: %s\n", OptArg);
	            exit(-1);
 		 }
 		 Crea_RaDec_Table = True; 
 		 break;
 	   case 'A': 
 	        if (sscanf(OptArg,"%f",&FluxMult) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad flux parameter: %s\n", OptArg);
		    exit(-1);
		}
		break;
	   case 'p':
		 KeepPositivSup=False;
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

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
       
 	// Test inconsistencies in the option call 
 	if (Stat_Noise == NOISE_EVENT_POISSON)
 	{
 	   if (Transform != TO_PAVE_BSPLINE)
 	   {
 	      cerr << "WARNING: with this noise model, only the BSPLINE A TROUS can be used ... " << endl;
 
 	      cerr << "        Type transform is set to: BSPLINE A TROUS ALGORITHM " << endl;
 	   }
 	   Transform = TO_PAVE_BSPLINE;
        }
	
 	if ((BgrMethod != BGR_MR) && (Optb == True))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: bad background option -b with this bgr method. " << endl;
           exit(-1);
  	}
	if ((BgrMethod == BGR_IMA) && (OptB == False))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: a background image must be given (-B option) " << endl;
           exit(-1);
  	}
	if ((BgrMethod == BGR_VALUE) && (OptL == False))
	{
 	   cerr << endl << endl;
 	   cerr << "  Error: a background value must be given (-L option) " << endl;
           exit(-1);
  	}
	if ((OptL == True) && (OptB == True))
	{
	   cerr << endl << endl;
 	   cerr << "  Error: -B and -L options are exclusive. " << endl;
           exit(-1);
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
 	}
	
	if ((UsePSF == False) && (Deconv == True) && (ReconsMethod != GRAD_PSF_XMM))
	{
	   cerr << "ERROR: the PSF must be given for deconvolution " << endl;
	   cerr << "       using the -P option ... " << endl;
	   exit(-1);
	}   
	
	if ((UsePSF == True) && 
	     ((Transform == TM_TO_PYR) || (Transform == TM_TO_SEMI_PYR)))
	{
	   cerr << "ERROR: Deconvolution cannot be performed with this multiscale transform ... " << endl;
 	   exit(-1);
	}   
	if (UseMethod == True)
	{
	   if ((UsePSF == True) && (ReconsMethod != GRAD_PSF))
	   {
 	      cerr << "WARNING: for deconvolution, the reconstruction method must be: " << endl;
	      cerr << "         " << StringRecObj(GRAD_PSF) << endl;
	      ReconsMethod = GRAD_PSF;
	    } 
	}
	else if (UsePSF == True) ReconsMethod = GRAD_PSF;

	if ((AperPhot == True) && ((ReconsMethod == GRAD_PSF) || (ReconsMethod == GRAD_PSF_XMM)))
	{
	   cerr << "ERROR: aperture photometry cannot be used with deconvolution ... " << endl;
	   exit(-1);
	}
	
	if ((EpsOpt == False) && (UseNSigma == True))
            EpsilonPoisson = (1. - erff((double) N_Sigma / sqrt((double) 2.)));
		   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}

/*********************************************************************/


int main(int argc, char *argv[])
{
    int i,j,k,s;
    Ifloat Dat;
    fitsstruct Header, FitsIO_HD;
    char Cmd[256];
    ListObj2D Objects;
   char nom_gr[256];
   char nom_grobj[256];
   char nom_mes[256];
   char nom_mes_ani[256];
   char nom_psobj[256];
   char nom_psobj_ani[256];
   char nom_tex[256];
   char nom_tex_ani[256];
   char nom_radec[256];
   char nom_radec_ani[256];
   char nom_mes_sub[256];
   char nom_anisotrop[256];
   char nom_psobj_sub[256];
   char nom_tex_sub[256];
   char nom_radec_sub[256];
   int AddBorderY = (ReconsMethod != GRAD_PSF_XMM) ? 0:  (int) (XMM_PSF_NL / ZoomPSF / 2);
   int AddBorderX = (ReconsMethod != GRAD_PSF_XMM) ? 0:  (int) (XMM_PSF_NC / ZoomPSF / 2);
   extern Bool VisionCleanBord; // see MR_VisTree.cc
    
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);
   sprintf(nom_anisotrop, "%s%s", Name_Imag_Out, "_ani");
   sprintf(nom_gr, "%s%s", Name_Imag_Out, ".mrg");
   sprintf(nom_grobj, "%s%s", Name_Imag_Out, ".mro");
   
   sprintf(nom_mes, "%s%s", Name_Imag_Out, ".mes");
   sprintf(nom_mes_ani, "%s%s", Name_Imag_Out, "_ani.mes");
   
   sprintf(nom_psobj, "%s%s%s", Name_Imag_Out, "_obj", ".ps");
   sprintf(nom_psobj_ani, "%s%s%s", Name_Imag_Out, "_obj_ani", ".ps");

   sprintf(nom_tex, "%s%s", Name_Imag_Out, ".tex");
   sprintf(nom_tex_ani, "%s%s", Name_Imag_Out, "_ani.tex");
   
   sprintf(nom_radec, "%s%s%s", Name_Imag_Out, "_radec",".tex");
   sprintf(nom_radec_ani, "%s%s%s", Name_Imag_Out, "_radec_ani",".tex");

   sprintf(nom_mes_sub, "%s%s", Name_Imag_Out, "_sub.mes");
   sprintf(nom_psobj_sub, "%s%s%s", Name_Imag_Out, "_sub_obj", ".ps");
   sprintf(nom_tex_sub, "%s%s", Name_Imag_Out, "_sub.tex");
   sprintf(nom_radec_sub, "%s%s%s", Name_Imag_Out, "_sub_radec",".tex");
      
    if (Verbose == True)
    {
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
    cout << "Transform = " << StringTransform(Transform) << endl;
    cout << "Reconstruction method = " << StringRecObj(ReconsMethod) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "Noise Type = " << StringNoise(Stat_Noise) << endl;
    cout << "Sigma Noise = " << Noise_Ima << endl;
    cout << "N_Sigma = " << N_Sigma << endl;
    cout << "EpsilonPoisson = " << EpsilonPoisson << endl;

    if (Stat_Noise == NOISE_GAUSSIAN)
                             cout << "Type of Noise = GAUSSIAN" << endl;
    else
    {
      if (Stat_Noise == NOISE_POISSON)
      {
           cout << "Type of Noise = POISSON" << endl;
           cout << "  Gain = " << PasCodeur << endl;
           cout << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
           cout << "  Read-out Mean = " << MeanGauss << endl;
      }
    }
    if (Deconv == True) cout << "Restore deconvolved objects " << endl;
     if (AperPhot == True)
    {
        cout << "Aperture photometry" << endl;
	cout << "  Background estimation method = " << StringBgr(BgrMethod) << endl;
        switch (BgrMethod)
        {
            case BGR_VALUE: cout << "  Background Value = " << BackgroundValue << endl; break;
	    case BGR_IMA: cout << "  Background image file namme = " << Name_BGR << endl;break;
	    // case BGR_MR: cout << "  Automic background estimation " << endl;
	    default: break;
	}
     }
     cout << "Vision Model = " << StringMVM (VisionModel) << endl;;
   }

    io_read_ima_float(Name_Imag_In, Dat, &Header);
    Header.origin = Cmd;
    check_scale(Dat.nl(), Dat.nc(), Nbr_Plan);

    // Read Celestial coordinates with cfitsio routines
    if (Crea_RaDec_Table == True) 
    {
      if (cfitstio_read_celestial_coord(Name_Imag_In, &FitsIO_HD))
         Crea_RaDec_Table = False;
    }
    
    // Compute the thresholded wavelet coefficents  
    MultiResol MR_Data (Dat.nl(), Dat.nc(), Nbr_Plan, 
                        Transform, "MR_Transform");
    MR_Data.ModifiedATWT = True;
    MR_Data.Border = I_MIRROR;
    int NbrBand = MR_Data.nbr_band();
    if (Stat_Noise == NOISE_EVENT_POISSON) MR_Data.Border = I_MIRROR;
    else MR_Data.Border = I_CONT;
    // MR_Data.Border = I_MIRROR;
    MRNoiseModel ModelData(Stat_Noise, Dat.nl(), Dat.nc(), 
                           Nbr_Plan, Transform);
    // if (Transform == TM_TO_PYR) ModelData.DilateSupport = True;

   if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;
    if (UseNSigma  == True)
          for (i=0; i <  NbrBand; i++) 
    {
        ModelData.NSigma[i]=N_Sigma;
    }
    if (N_Sigma  == 111)
    {
        for (i=0; i < 4; i++)  ModelData.NSigma[i] = 5;
	ModelData.NSigma[4] = 4.;
	ModelData.NSigma[5] = 3.5;
 	for (i=6; i < NbrBand; i++) ModelData.NSigma[i] = 3.;
    }
    for (s=0; s <  NbrBand; s++) ModelData.TabEps[s] = EpsilonPoisson;
    // ModelData.FirstDectectScale = FirstScale;
    ModelData.MinEventNumber = MinEvent;
    ModelData.NiterSigmaClip = NiterClip;
    ModelData.SizeBlockSigmaNoise = SizeBlock;
    ModelData.OnlyPositivDetect = KeepPositivSup;
    ModelData.CCD_Gain = PasCodeur;
    ModelData.CCD_ReadOutSigma = SigmaGauss;
    ModelData.CCD_ReadOutMean = MeanGauss;
    ModelData.TypeThreshold=TypeThreshold;
    ModelData.SupIsol = True;
    VisionCleanBord = (KeepBorderObj == True) ? False : True;
    if (GetRmsCoeffGauss == False) ModelData.GetRmsCoeffGauss = False;
    
    if (UseRMSMap == True)
    {
       ModelData.UseRmsMap = True;
       io_read_ima_float(Name_RMSMap, ModelData.RmsMap);
    }
    if (UsePSF == True) 
    {
       io_read_ima_float(Name_Psf, Objects.Psf);
       AddBorderY =  (int) (Objects.Psf.nl() / ZoomPSF / 2);
       AddBorderX =  (int) (Objects.Psf.nc() / ZoomPSF / 2);
    }

    Ifloat Result(Dat.nl(), Dat.nc(), "res");
    ModelData.model(Dat, MR_Data);
    
    if ((Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN))
                 cout << "Sigma noise = " << ModelData.SigmaNoise << endl;
    if (Verbose == True) cout << "Graph creation ..." << endl;

    Foret F, *F_obj;
    MR2D_Seg MVM;
    type_filter Filter = FILTER_ITER_THRESHOLD;
    MRFiltering CFilter(ModelData, Filter);
    CFilter.Max_Iter = 10;
    CFilter.KillLastScale = True;
    CFilter.PositivIma = True;
    CFilter.Verbose = Verbose; 
    CFilter.IsDataModelled = True;
    CFilter.UseAdjointRec  = False;
    CFilter.Sup_Set = False;
    if (VisionModel == MVM_RUEBIJ)
    {
       if (Stat_Noise == NOISE_EVENT_POISSON)  MR_Data.transform(Dat);
//        else
//        {
//        cout << "FFB" << endl;
//        CFilter.Verbose = True;
//        CFilter.filter(Dat, Result);
// 	   cout << "FF1" << endl;
//         MR_Data.transform(Result);
// 	ModelData.threshold(MR_Data);
//        }
       
       ModelData.threshold(MR_Data);
       if (Verbose == True) cout << "Segmentation ... "  << endl;
       mr_make_graph(MR_Data, F, SubSegment, KillIsolObj);
       if (Verbose == True) cout << "Objects tree creation ... "  << endl;
       if (TypeMVM == MVM_Align)  mr_make_obj_tree(F, F_obj, False, True, DistMaxAl);
       else mr_make_obj_tree(F, F_obj, False, False);
    }
    else
    {
       int AddBorderY = (ReconsMethod != GRAD_PSF_XMM) ? 0:  (int) (XMM_PSF_NL / ZoomPSF / 2);
       int AddBorderX = (ReconsMethod != GRAD_PSF_XMM) ? 0:  (int) (XMM_PSF_NC / ZoomPSF / 2);
       int MinSizeObjX = 0;
       int MinSizeObjY = 0;

       for (s = 0; s < NbrBand-2; s++) ModelData.kill_isol(s);
       MVM.init(MR_Data);
       MVM.VisionModel = VisionModel;
       MVM.PercentFlux = 0.0;
       MVM.DistMin = 0;
       MVM.AddBorderX = AddBorderX;
       MVM.AddBorderY = AddBorderY;
       MVM.MinSizeObjX = MinSizeObjX;
       MVM.MinSizeObjY = MinSizeObjY;
       MVM.KillIsolatedObj = KillIsolObj;
       MVM.MVMDebug = MVMDebug;

       CFilter.filter(Dat, Result);
       if (Stat_Noise == NOISE_EVENT_POISSON) ModelData.Event_Image = Dat;  // ifloat(Dat, ModelData.Event_Image);

       MR_Data.transform(Result);
       for (s = 0; s < NbrBand-1; s++)
       {
           for (i=0; i < MR_Data.size_band_nl(s); i++)
           for (j=0; j < MR_Data.size_band_nc(s); j++) 
           {
              if (MR_Data(s,i,j) < ModelData.sigma(s,i,j)*ModelData.NSigma[s])    
                                                        MR_Data(s,i,j) = 0.;
           }
       }
       MVM.segment(MR_Data, ModelData);
       MR_Data.transform(Dat);
    }

   if (Verbose == True) 
   {
      cout << "Objects extraction ..." << endl;
      if (VisionModel == MVM_RUEBIJ)
      {
         cout << "Number of objects = " << F_obj->get_nb_obj()-F_obj->get_nb_s_obj() << endl;
         cout << "Number of sub-objects = "<< F_obj->get_nb_s_obj() <<endl;
      }
      else
        cout << "Number of objects = " <<  MVM.nbr_obj() << endl;
   }

   // object class initialization and objects extractions 
   Objects.NameImagIn = strdup(Name_Imag_In);
   Objects.FirstUseScale = FirstScale;  // deb_rec;
   Objects.LastUseScale = LastScale;    // fin_rec;
   Objects.Nb_iter_rec = Nb_iter_rec;
   Objects.KeepIma = False;
   Objects.ReconsMethod = ReconsMethod;
   Objects.ErrorRec = Eps_ErrorRec;
   Objects.Verbose = Verbose;
   Objects.FluxMult = FluxMult;
   Objects.AperPhot = AperPhot;
   Objects.BgrMethod  = BgrMethod;
   Objects.BackgroundValue = BackgroundValue;
   Objects.KSigmaAper = KSigmaAper;
   Objects.BGRSize = BGRSize;
   Objects.VisionModel = VisionModel;
   
   if (VisionModel == MVM_RUEBIJ) Objects.F_obj = F_obj;
   else  Objects.MVM = &MVM;

   if ((AperPhot == True) && (BgrMethod == BGR_IMA))
                         io_read_ima_float(Name_BGR,  Objects.BgrModel);  
   if (Fwhm > 0) Objects.Fwhm = Fwhm;
   if (ZoomPSF > 1) Objects.ZoomPSF  = ZoomPSF;   


   Objects.Deconv = Deconv;  
   Objects.IsotropSeparation = AnisotropAnalysis;
   Objects.AnitropicSeparationParam = AnisotropParam;
   Objects.create_list(Dat, ModelData, MR_Data, 
                       WriteAllObj, InfoSubObj, WriteObjFullSize);
   MR_Data.free();
   if (Crea_RaDec_Table == True)
   {
      double Ra, Dec;
      Object_2D *TabObj = Objects.TabObj;
      printf("crpixx = %f, crpixy = %f, crvalx = %f,crvaly = %f, cdeltx = %f, cdelty = %f\n", 
             FitsIO_HD.crpixx, FitsIO_HD.crpixy,FitsIO_HD.crvalx,FitsIO_HD.crvaly,FitsIO_HD.cdeltx,FitsIO_HD.cdelty);
      for (i = 0; i < Objects.Nb_TotalObj; i++)
      {
       xyad(&FitsIO_HD,  TabObj[i].PosX, TabObj[i].PosY, &Ra, &Dec);
       TabObj[i].Ra = Ra;
       TabObj[i].Dec = Dec;
      }
   }
        
  if (Verbose == True)  cout << "Write result " << endl;
  Header.bitpix = BP_FLOAT;

  io_write_ima_float(Name_Imag_Out, Objects.Ima_SumObj, &Header);
  if (Objects.IsotropSeparation == True) io_write_ima_float(nom_anisotrop, Objects.Ima_SumObjAnisotopic, &Header);

  // write all parameters in a ascii file
  mr_write_obj_ascii(Objects, nom_mes, True, False);
  if (Objects.IsotropSeparation == True) mr_write_obj_ascii(Objects, nom_mes_ani, True, False, False);
  
  // write the position of the object in a PS map
  mr_write_obj_ps(Objects,  nom_psobj, True, False);
  if (Objects.IsotropSeparation == True) mr_write_obj_ps(Objects,  nom_psobj_ani, True, False, False);
  
  // write object informations to a tex file 
  mr_write_obj_tex(Objects, nom_tex, True, False);
  if (Objects.IsotropSeparation == True) mr_write_obj_tex(Objects, nom_tex_ani, True, False, False);
  
  if (Crea_RaDec_Table == True)
  { 
     mr_write_radec_tex(Objects, nom_radec, Order, True, False);
     if (Objects.IsotropSeparation == True) mr_write_radec_tex(Objects, nom_radec_ani, Order, True, False, False);
  }
  
  if ((InfoSubObj == True) && (Objects.Nb_TotalSubObj > 0))
  {
      mr_write_obj_ascii(Objects, nom_mes_sub, True, True);
      mr_write_obj_ps(Objects,  nom_psobj_sub, True, True);
      mr_write_obj_tex(Objects, nom_tex_sub, True, True);
      if (Crea_RaDec_Table == True) 
          mr_write_radec_tex(Objects, nom_radec_sub, Order, True, True);
  }
   
  if (WriteSimu == True)
  {
      if (Verbose == True)  cout << "Draw objects" << endl;
      Ifloat DatAni;
      
      if (Objects.IsotropSeparation == True) DatAni = Dat;
      
      float Val = max(Dat) * 1.5;
      for (i = 0; i < Objects.Nb_TotalObj; i++)
        if ((Objects.IsotropSeparation == False) || (Objects.TabObj[i].Isotrop == True))
	                                                      Objects.TabObj[i].draw(Dat, False, Val);
      io_write_ima_float("xx_ellips", Dat, &Header);
       
      if (Objects.IsotropSeparation == True)
      {
         for (i = 0; i < Objects.Nb_TotalObj; i++)
            if (Objects.TabObj[i].Isotrop == False) Objects.TabObj[i].draw(DatAni, False, Val);
         io_write_ima_float("xx_ellips_ani", DatAni, &Header);
      }
      Dat.init();

      if (Verbose == True)        cout << "Simu objects" << endl;

      for (i = 0; i < Objects.Nb_TotalObj; i++)
          Objects.TabObj[i].simu(Dat);
      io_write_ima_float("xx_simu", Dat, &Header);
  }
  if (VisionModel == MVM_RUEBIJ) if (F_obj != NULL) delete F_obj;
  exit(0);
} 

