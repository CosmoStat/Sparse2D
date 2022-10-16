/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Yanling Fang && Yves Bobichon && Jean-Luc Starck
**
**    Date:  97/12/04
**    
**    File:  mr_rfilter.cc
**
*******************************************************************************
**
**    DESCRIPTION  filter an image with Rayleigh noise
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_rfilter option image output
**    
**  
**
**
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_rfilter.cc 1.0 97/11/26 CEA 1997 @(#)";
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "MR_Filter.h"
   
char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
int Nbr_Plan=6;  /* default number of scales */
int Max_Iter = DEFAULT_MAX_ITER_FILTER; /* Maximum number of iteration */
int Nbr_Image=1;
Bool Verbose;

float Epsilon = DEFAULT_EPSILON; /* threshold  parameter */
float Converge = DEFAULT_EPSILON_FILTERING; /* convergence parameter */
float Ratio =10.;

Bool Enter_Epsilon = False;
Bool WindowOpt=False; 

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, char *opts);

type_transform Transform = TO_PAVE_BSPLINE; /* type of transform */

fitsstruct Header;
 
type_noise_speckle Stat_Noise = NOISE_LOG_RAYLEIGH;

static void usage(char *argv[])
{
int i;
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "\n");

    fprintf(OUTMAN, "         [-t type of noise]\n");
    for (i = 0; i < 3; i++)
      fprintf(OUTMAN, "              %d: %s \n",i,
	      StringTypeSpeckle( (type_noise_speckle)i ));
    fprintf(OUTMAN, "              default is %s\n",  StringTypeSpeckle(Stat_Noise));
    manline();


    nbr_scale_usage(Nbr_Plan);
    manline();
        
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-i iteration]\n");
    fprintf(OUTMAN, "              number max of iterations\n");
    fprintf(OUTMAN, "              default is %d\n",DEFAULT_MAX_ITER_FILTER);
    manline();

    manline();
    converg_param_usage(Epsilon);
    manline();

    fprintf(OUTMAN, "         [-E decision level]\n");
    fprintf(OUTMAN, "             statistical decision level used for all scales\n");
    fprintf(OUTMAN, "             if 0 : the value of Epsilon is set interactivally\n");
    fprintf(OUTMAN, "             for each scale (typically 0.01, 0.001, 0.0001)\n");
    fprintf(OUTMAN, "             default is %f for all scales\n",DEFAULT_EPSILON_FILTERING);
    manline();

    fprintf(OUTMAN, "         [-N number of images ]\n");
    fprintf(OUTMAN, "             if the input image is a multi-look image \n");
    fprintf(OUTMAN, "             i.e. input image is the average of N single-look images \n");
    fprintf(OUTMAN, "             default is 1\n");
    manline();

    fprintf(OUTMAN, "         [-r ratio\n");
    fprintf(OUTMAN, "             ratio=Epsilon2/Epsilon1 use for soft thresholding.\n");
    fprintf(OUTMAN, "             If set, soft thresholding is applied between the thresholds\n");
    fprintf(OUTMAN, "             defined by Epsilon1 and Epsilon2, and a hierarchical \n");
    fprintf(OUTMAN, "             thresholding is use for the first scale. \n");
    fprintf(OUTMAN, "             The value of Epsilon1 is set with -E option. \n");
    fprintf(OUTMAN, "             if ratio=1, standard hard thresholding is used for all scales\n");
    fprintf(OUTMAN, "             default ratio is 10\n");
    manline();
    fprintf(OUTMAN, "         [-u\n");
    fprintf(OUTMAN, "         Use the undecimated bi-orthogonal Wavelet transform\n");
    fprintf(OUTMAN, "         instead of the a trous isotropic Wavelet Transform.\n");
    vm_usage();
    manline();
    verbose_usage();
    manline();

 
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
    /* get options */
    while ((c = GetOpt(argc,argv,"uvt:n:r:i:e:E:N:dzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'u': Transform = TO_UNDECIMATED_MALLAT; break;
	   case 't':
		/* -t <Stat_Noise> of noise */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(stderr, "bad noise type: %s\n", OptArg);
	            usage(argv); 
		}
                if ((c >= 0) && (c < 3))  Stat_Noise = (type_noise_speckle) (c);
                else  
                {
		    fprintf(stderr, "bad type of noise: %s\n", OptArg);
	            usage(argv);
 		}
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
	   case 'N':
		/* -N <Number of image> */
		if (sscanf(OptArg,"%d",&Nbr_Image) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of image: %s\n", OptArg);
		    exit(-1);
		}
                if (Nbr_Plan < 1) 
                 {
		    fprintf(OUTMAN, "Error: bad number of image: %s\n", OptArg);
		    fprintf(OUTMAN, "       0 < Nbr_Image\n");
 		    exit(-1);
		}
		break;
 	   case 'i':
		/* -i < Number of iterations> */
		if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
				exit(-1);
		}
                if (Max_Iter <= 0)   Max_Iter = DEFAULT_MAX_ITER_FILTER;
		break;
	   case 'e': 
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Converge) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
                if ((Converge < 0) || (Converge> 1.)) 
                        Converge  = DEFAULT_EPSILON_FILTERING;
		break;
	 case 'E':
		if (sscanf(OptArg,"%f",&Epsilon) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad statistical decision level: %s\n", OptArg);
		    exit(-1);
		}
		if (Epsilon==0.) Enter_Epsilon=True;
		if ((Epsilon <0.) || (Epsilon> MAX_EPSILON)) 
                {
		    fprintf(OUTMAN, "Error: bad statistical decision level %s\n", OptArg);
		    fprintf(OUTMAN, "       %f <= Epsilon <= %f\n",MIN_EPSILON, MAX_EPSILON);
 		    exit(-1);
		}
		break;
	 case 'r':
		/* -r < ratio for soft thresholding> */

		if (sscanf(OptArg,"%f",&Ratio) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad Ratio parameter : %s\n", OptArg);
		    exit(-1);
		}
		
		if (Ratio <1) 
                {
		    fprintf(OUTMAN, "Error: bad Ratio parameter %s\n", OptArg);
		    fprintf(OUTMAN, "Ratio must be >= 1\n");
 		    exit(-1);
		}
		break;
            case 'v':
                /* Verbose flag -v */
                Verbose = True;
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
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif    
}


/****************************************************************************/

int main(int argc, char *argv[])
{
  int k;
  Ifloat Imag;
  fltarray Tab_Epsilon;
  char Cmd[256];
  
  Cmd[0]= '\0';
  for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
  
  /* Get command line arguments, open input file(s) if necessary */
  lm_check(LIC_MR1);
  filtinit(argc, argv);
  
   if (Verbose == True)
    {
      cout << endl << endl << "PARAMETERS: " << endl << endl;
      cout << "File Name in = " << Name_Imag_In << endl;
      cout << "File Name Out = " << Name_Imag_Out << endl;
      cout << "Transform = " << StringTransform(Transform) << endl;
      cout << "Number of scales = " << Nbr_Plan << endl;
      if (!Enter_Epsilon) cout << "Epsilon = " << Epsilon << endl;
      cout << "Ratio = " << Ratio << endl;
      cout << "Max iteration = " << Max_Iter << endl;
      cout << "Max Delta_Sigma = " <<  Converge << endl;
      cout << "Noise type = " << StringTypeSpeckle(Stat_Noise) << endl;
      cout << "Nbr of images = " << Nbr_Image   << endl;
    } 

  int NbrBand = number_band_per_resol(Transform)*(Nbr_Plan-1);
  Tab_Epsilon.alloc(NbrBand);
  if (Enter_Epsilon)
    {
      cout << "enter Epsilon value: " << endl;  
        for (int s = 0; s < NbrBand; s++)
	  {
	    if (Verbose == True) cout << "Epsilon[" << s << "] = ";
	    scanf("%f",&Tab_Epsilon(s));
		if ((Tab_Epsilon(s) <0.) || (Tab_Epsilon(s)> MAX_EPSILON)) 
                {
		    fprintf(OUTMAN, "Error: bad statistical decision level %f\n",Tab_Epsilon(s));
		    fprintf(OUTMAN, "       %f <= Epsilon <= %f\n",MIN_EPSILON, MAX_EPSILON);
 		    exit(-1);
		}
	    if (s>=4) 
	      {
		if (Verbose== True) cout << "Epsilon automatically set to " <<  Tab_Epsilon(s) << " for larger scales" << endl; 
		for (int s1 = s; s1 < (Nbr_Plan-1); s1++) 
		  Tab_Epsilon(s1)=Tab_Epsilon(s);
		s=Nbr_Plan-1;
	      }
	  }

    }
  else for (int s = 0; s < NbrBand; s++) Tab_Epsilon(s)=Epsilon;

  io_read_ima_float(Name_Imag_In,  Imag, &Header);
  Header.origin = Cmd;
  Ifloat Result (Imag.nl(),Imag.nc(), "Result");
  check_scale(Imag.nl(),  Imag.nc(), Nbr_Plan);

  if (Verbose == True) cout << "compute multiresolution noise model" << endl;
  // set noise model
  StatRayleigh NoiseModel(Stat_Noise, Nbr_Plan, Nbr_Image, Transform);
  // NoiseModel.write();
  NoiseModel.Verbose = Verbose;
  if (Verbose == True) cout << "Filtering ..." << endl;
  // image filtering
  mr_sfilter(Imag, Result, NoiseModel, Tab_Epsilon, Ratio, Max_Iter, Converge);

  if (Verbose == True) cout << "Output filtered image write in : " << Name_Imag_Out << endl ;
  io_write_ima_float(Name_Imag_Out, Result, &Header);
exit(0);
} 

/***********************************************************************************/


