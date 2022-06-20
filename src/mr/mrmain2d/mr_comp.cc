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
**    File:  mr_comp.cc
**
*******************************************************************************
**
**    DECRIPTION  compression program
**    ---------- 
**
*******************************************************************************
**
**    PARAMETRES
**    ----------
**
**        Input_File: image file name to compress
**        [Output_File]: compressed image file
**                       by default, "Input_File.MRC"
**           [-v] : vebose
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
**           [-i number_of_iterations] : number of iterations to do the
**                                     compression by default, 0
**
**           [-s nsigma ] : a thresholding is done at each scale at 
**                        nsigma * Sigma_Noise
**                        by default, nsigma = 3
**
**           [-n number of scales] : number of scales used by the 
**                                 multiresolution transform
**                                 by default number_of_scales = 4
**           [-r]  : if -r option is given, the noise is compressed to
**                 with a step of sigma_noise/2. By default, the noise 
**                 is not conserved
**
**           [-q SignalQuantif] : 
**               The signal is quantified by Sq = S/(SignalQuantif*Sigma_Noise)
**               by default,SignalQuantif = 3
**
**           [-e NoiseQuantif]
**               The signal is quantified by Nq = S/(NoiseQuantif*Sigma_Noise)
**               by default, NoiseQuantif = 3
**
**           [-k]\n");
**               Keep isolated pixel in the support 
**               at the first scale. Default is no. 
**               If the PSF is on only one pixel, this 
**               option should be set\n");
**      
******************************************************************************/

// static char sccsid[] = "@(#)mr_comp.cc 3.1 96/05/02 CEA 1995 @(#)";

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"

char File_Name_Imag[256];  /* input file image */
char File_Name_Transform[256]; /* output file name */
int Nbr_Plan=DEFAULT_CMP_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
float Fwhm=0.5;                   /* Full width at half maximum */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = TM_PYR_MEDIAN; /* type of transform */
Bool SupIsol=True;             /* suppress isolated pixel in the support */
int MaxIter = DEFAULT_MAX_ITER_COMP;
float NoiseQuantif = DEFAULT_NOISE_QUANTIF;
float SignalQuantif = DEFAULT_SIGNAL_QUANTIF;
int KeepResi = 0;
int KeepFitsHeader = 0;
Bool NoBscale=False; // for integer coding. Fits image are not divided
                     // by the BScale fits keyword
Bool SupNeg=False;
Bool Verbose=False;
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

#define WRITE_PARAM 0
extern float BadPixalVal; 
extern Bool BadPixel;

type_comp Comp_Method = COMP_MRMEDIAN;
int WindowSize = 5;
int Elstr_Size = 5;
int Elstr_Shape = 1;
Bool UseMR_ForBgr = False;
int Npix_MR_Bgr = DEFAULT_NPIX_BGR;
Bool KillDiag0=False;

float CompressionRatio=-1.;
Bool UseBudget=False;
Bool NoiseInData = True;

int BlockSize=0;
Bool UseBlock=False;

/*********************************************************************/

static void usage(char *argv[])
{
   int i;
   
    fprintf(OUTMAN, "Usage: %s options in_image [out_file]\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-m Compression_Method]\n");
    for (i = 0; i <  NBR_COMP_METHOD; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringComp((type_comp)i));
    fprintf(OUTMAN, "              default is %s\n", StringComp((type_comp) Comp_Method));

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
   
    mrcomp_option_usage(Elstr_Size,  SignalQuantif,  NoiseQuantif, WindowSize);
    manline();    
    vm_usage();
    manline();           
    verbose_usage();
          
//    fprintf(OUTMAN, "         [-H]\n");
//    fprintf(OUTMAN, "              Kill the first diagonal scale. \n");
//    fprintf(OUTMAN, "              Only used with orthogonal transform.\n");
//    fprintf(OUTMAN, "              Default is No\n");
    manline();
    manline();
    manline(); 
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void hcinit(int argc, char *argv[])
{
    int c, L;
    char *Ptr;
    Bool OptInt = False;
    Bool Opte = False;
    Bool OptS = False;
    Bool OptD = False;
    Bool Optq = False;
    Bool Opts = False;
    Bool Optn = False;
    Bool OptNoise = False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif

    /* get options */
    while ((c = GetOpt(argc,argv,"m:vs:n:rc::pg:e:q:klfb:OWBPSD:MHR:NC:i:zZ:")) != -1) 
    {
	switch (c) 
        {
           case 'C': 
              if (sscanf(OptArg,"%d",&BlockSize ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
	            exit(-1);
                    
		}
		UseBlock = True;
                break;
           case 'R': 
                if (sscanf(OptArg,"%f", &CompressionRatio) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
		UseBudget = True;
		break;
           case 'N':NoiseInData = False;break;
           case 'H': KillDiag0=True;break;
	   case 'm':
		/* -f <type> type of compression */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of compression: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_COMP_METHOD)) Comp_Method = (type_comp) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of compression: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	    case 'v':
		/* Verbose flag -v */
		Verbose = True;
		break;
            case 'M':
                UseMR_ForBgr = True;
                break;
            case 'p':
                /* Poisson noise */
                Stat_Noise = NOISE_POISSON;
                Noise_Ima = 1.;
                if (OptNoise == True)
                {
                 fprintf(OUTMAN, "Error: g,p,c options cannot be used together ... \n");
                 exit(-1);
                }
                 OptNoise = True;
               break;
            case 'f':
                /* keep all keywords in fits header */
                KeepFitsHeader = 1;
               break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		if (OptNoise == True)
		{
                 fprintf(OUTMAN, "Error: g,p,c options cannot be used together ... \n");
                 exit(-1);
                }
                OptNoise = True;
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
                if (OptNoise == True)
                {
                 fprintf(OUTMAN, "Error: g,p,c options cannot be used together ... \n");
                 exit(-1);
		}
                OptNoise = True;
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
		Optn = True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)  N_Sigma = DEFAULT_N_SIGMA;
                Opts= True;
		break;
            case 'k':
               /* kill i */
               SupIsol = False;
               break;
            case 'l':
               /* save the residual */
               KeepResi = 2;
               break;
            case 'O':
               OptInt = True;
		break;
            case 'W':
               WindowSize = 3;
		break;
            case 'B':
                NoBscale = True;
 		break;
            case 'P':
               SupNeg = True;
		break;
	   case 'q':
		/* -c <SignalQuantif> */
		if (sscanf(OptArg,"%f",&SignalQuantif) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad SignalQuantif: %s\n", OptArg);
				exit(-1);
		}
		Optq = True;
		break;
	   case 'e':
		/* -c <NoiseQuantif> */
		if (sscanf(OptArg,"%f",&NoiseQuantif) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad NoiseQuantif: %s\n", OptArg);
				exit(-1);
		}
		Opte=True;
		break;
	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&MaxIter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad MaxIter: %s\n", OptArg);
				exit(-1);
		}
		break;           
	   case 'S':
                //NoBscale = True;
                // Comp_Method = COMP_O_MORPHO;
                Elstr_Shape = 0; //square for morphomate compression
                OptS = True;
                break;
            case 'D':
		/* -d <Elstr_Size> */
		if (sscanf(OptArg,"%d",&Elstr_Size) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad Elstr_Size: %s\n", OptArg);
		    exit(-1);
		}
                if ((Elstr_Size < 3.) || (Elstr_Size > 10.)) 
                                        Elstr_Size = DEFAULT_ELSTR_SIZE;
                OptD = True;
                break;
            case 'r':
               /* compress the residual */
               KeepResi = 1;
               break;
            case 'b':
		if (sscanf(OptArg,"%f",&BadPixalVal) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad Bad_Pixel value: %s\n", OptArg);
		    exit(-1);
		}
		BadPixel = True;
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

        strcpy (File_Name_Transform, File_Name_Imag);

	if (OptInd < argc) 
             strcpy(File_Name_Transform, argv[OptInd++]);

	L = strlen (File_Name_Transform);
	Ptr = File_Name_Transform;
	if (L < 5) strcat (File_Name_Transform, ".MRC");
	else if ((Ptr[L-1] != 'C') || (Ptr[L-2] != 'R') || (Ptr[L-3] != 'M')) 
               strcat (File_Name_Transform, ".MRC");

        // S,D valid only for morphology
        if (Comp_Method != COMP_MORPHO)  
        {
           if ( (OptS == True) || (OptD == True) || (UseMR_ForBgr == True))
           {
              fprintf(OUTMAN, "Error: S,D options are not valid with this compression method ...\n");
               exit(-1);
            }
         }
        if ((WindowSize == 3) && (Comp_Method != COMP_MRMEDIAN))
        {
           fprintf(OUTMAN, "Error: W options is not valid with this compression method ...\n");
           exit(-1);
        }
        
        // if -O option is set, compression method must be changed.
        if (OptInt == True)
        {
           if (Comp_Method == COMP_MRMEDIAN) Comp_Method = COMP_O_MRMEDIAN;
           else if (Comp_Method == COMP_MORPHO) Comp_Method = COMP_O_MORPHO;
           else if (Comp_Method == COMP_MIN_MAX) Comp_Method = COMP_O_MIN_MAX;
           else 
           {
		fprintf(OUTMAN, "Error: -O option is not valid with this compression method ...\n");
		exit(-1);
	   }
	}
	
	// if mathematical morphology, don't use scale the values
        if ((Comp_Method == COMP_O_MORPHO) || (Comp_Method == COMP_MIN_MAX))
         NoBscale = True;
         
        // if mathematical morphology, some options are not valid
        if (((Comp_Method == COMP_O_MIN_MAX) ||
             (Comp_Method == COMP_MORPHO) ||
             (Comp_Method == COMP_O_MORPHO)) && ( (KeepResi > 0) 
                                                || (Opte == True)
                                                || (Optq == True)
                                                || (NoiseInData == False)
                                                || (UseBudget==True)))
        {
           fprintf(OUTMAN, "Error: R,N,q,r,l,e options are not valid with this compression method ...\n");
           exit(-1);
        }
        if (((Comp_Method == COMP_MORPHO) ||
             (Comp_Method == COMP_O_MORPHO)) && ( (Optn==True) 
                                                || (Opts == True)  ))
        {
           fprintf(OUTMAN, "Error: s,n options are not valid with this compression method ...\n");
           exit(-1);
        }
        
	if ((UseBlock == True) && (BadPixel == True))
	{
	   cerr << "Error: bad pixel option is not valid using block compression ..." << endl;
	   exit(-1);
        }
         
/*	if ((UseBudget == False) && (NoiseInData == False))
	{
 		fprintf(OUTMAN, "Error: -N option is valid only if -R is set ...\n");
		exit(-1);
	}*/
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***************************************************/

int main(int argc, char *argv[])
{
    int k;
    char Cmd[256];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, 
       open input file(s) if necessary */
    lm_check(LIC_MR1);
    hcinit(argc, argv);
    MR_CompData CompIma;
 
    switch (Comp_Method)
    {       
       case COMP_HAAR: Transform = TO_HAAR; break;
       case COMP_MALLAT: Transform = TO_MALLAT;break;
       case COMP_FEAUVEAU: Transform = TO_FEAUVEAU;break;
       case COMP_MIN_MAX: Transform = TM_MIN_MAX; break;
       case COMP_MRMEDIAN: Transform = TM_PYR_MEDIAN;  break;
       case COMP_O_MRMEDIAN: Transform = TM_PYR_MEDIAN; break;
       case COMP_O_MORPHO:
       case COMP_MORPHO:
       default: break;
    }  
    if (UseBlock == True)
    {
         CompIma.UseBlock = True;
         CompIma.BlockSize = BlockSize;
    }
    CompIma.UseBudget = UseBudget;
    CompIma.CompressionRatio = CompressionRatio;
    CompIma.NoiseInData = NoiseInData;
    
    CompIma.KillDiag0 = KillDiag0;
    CompIma.SupNeg=SupNeg;
    CompIma.NoBscale=NoBscale;
    CompIma.Cmd=Cmd;
    CompIma.File_Name_Imag=File_Name_Imag;
    CompIma.File_Name_Transform=File_Name_Transform;
    CompIma.Nbr_Plan=Nbr_Plan;
    CompIma.N_Sigma= N_Sigma;
    CompIma.Noise_Ima = Noise_Ima;
    CompIma.Stat_Noise = Stat_Noise;
    CompIma.Transform = Transform;
    CompIma.SupIsol = SupIsol;
    CompIma.MaxIter = MaxIter;
    CompIma.NoiseQuantif = NoiseQuantif;
    CompIma.SignalQuantif = SignalQuantif;
    CompIma.KeepResi = KeepResi;
    CompIma.KeepFitsHeader = KeepFitsHeader;
    CompIma.Comp_Method = Comp_Method;
    CompIma.MedianWindowSize = WindowSize;
    CompIma.Elstr_Size = Elstr_Size;  
    CompIma.Elstr_Shape = Elstr_Shape;
    CompIma.UseMR_for_BGR = UseMR_ForBgr;
    CompIma.NpixBgr = Npix_MR_Bgr;
    if (Verbose == True) CompIma.Verbose = True;
    CompIma.compress();

    exit(0);
}


/***************************************************/

 

