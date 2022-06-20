/******************************************************************************
**                   Copyright (C) 1995 by CEA
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
**    File:  mr_sigma.cc
**
*******************************************************************************
**
**    DESCRIPTION  estimate the noise standard deviation of an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_sigma option image 
**        where options = 
**          [-m type_of_method]
**		   1: 3_sigma_clipping
**                 2: mediane + 3_sigma_clipping method
**                 3: wavelet + 3_sigma_clipping method
**                 4: multiresolution support method
**
**                 default is 2
**
**           [-t type_of_multiresolution_transform]
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
**
**          [-n number_of_scales]
**                 number of scales used in the multiresolution
**
**                 default is 4
**
**	    [-p gain]
**	      performs the standard deviation of the gaussian part of the noise.
**		   gain = gain used in the model of image
**		   it supposes the input image is of the following form :
**		   	Poisson_Noise + Zero_Mean_Gaussian_Noise
**
**		   default is no
**
**
**
**          t option is only used if type_of_method is [2-3]
**          n,w options are only used if type_of_method is [3]
**	    p option can't be used with all other options
**
**
**
*********************************************************************/

// static char sccsid[] = "@(#)mr_sigma.cc 3.1 96/05/02 CEA 1995 @(#)";

           
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_Sigma.h"
#include "MR_Support.h"

char Name_Imag_In[256]; /* input file image */
float Noise_Ima=0.0;              /* noise standard deviation */
type_transform Transform = TO_PAVE_BSPLINE; /* type of transform */
int Max_Iter = DEFAULT_MAX_ITER_METHOD; /* Maximum number of iteration */
int Nbr_Plan = DEFAULT_NBR_SCALE ; /* number of scales in Multiresolution */
type_method Method = DEFAULT_SIGMA_METHOD; /* default method for calculating noise sigma */
Bool Usesupport = True; /* Use the Multiresolution support */
Bool WriteSup = False; /* Write support on disk */
Bool Gaussian_Flag = True;
Bool Poisson_Flag = False; /* checks the use of p option  */
		 	   /* and realises exclusivity with all other options */

Bool SearchReadout  = False; 
Bool SearchGain = False; 

float Gain = 1.0;		
float ReadOutNoise=0.;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;

 
/**************************************************************************************/

static void usage(char *argv[])
{
 
  fprintf(OUTMAN, "Usage: %s options in_image \n\n", argv[0]);
  fprintf(OUTMAN, "   where options =  \n");
  
  sigma_method_usage(Method);
  manline();

   nbr_scale_usage(Nbr_Plan);
   manline();
    
   sigma_gain();
   manline();
  
//  fprintf(OUTMAN, "\n");
//  fprintf(OUTMAN, "         [-r read_out_noise]\n");
//  fprintf(OUTMAN, "              performs the gain \n");
//  fprintf(OUTMAN, "              Input image must of the form :\n");
//  fprintf(OUTMAN, "              Gain * Poisson_Noise +  Zero_Mean_Gaussian_Noise :\n");
//  fprintf(OUTMAN, "              default is no\n");

   vm_usage();
   manline();
   verbose_usage();    
   manline();
   manline();
   exit(0);
}
/*****************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sigminit(int argc, char *argv[])
{
    int c;
    Bool OptPlan=False;
 #ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
   
    /* get options */
    while ((c = GetOpt(argc,argv,"m:n:p:vzZ:")) != -1) 
    {
	switch (c) 
	  {
	  case 'v': Verbose = True; break;
	  case 'm':
	    /* -m <type> type of method */
	    	if  (sscanf(OptArg,"%d",&c ) != 1) 
                {
	    	    fprintf(OUTMAN, "Error: bad type of method: %s\n", OptArg);
	            exit(-1);
                    
		}
		if ((c > 0) && (c <= NBR_METHODS)) Method = (type_method) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of mehtod: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	  case 'n':
		/* -n <Nbr_Plan> */
		if ( (sscanf(OptArg,"%d",&Nbr_Plan) != 1) || (Poisson_Flag == True) )

                {
		    if (Poisson_Flag == True)
			fprintf(OUTMAN, "ERROR : p option must be used alone !\n");
		    else
			fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
	            exit(-1);
 
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
 		    exit(-1);
		}
		OptPlan=True;
		break;
	  case 'p':
		if (sscanf(OptArg,"%f",&Gain) != 1) usage(argv);
		else 
		{
		    SearchReadout=True;
		    Poisson_Flag = True;
		    Gaussian_Flag = False;
		}
		break;
	  case 'r':
		if (sscanf(OptArg,"%f",&ReadOutNoise) != 1) usage(argv);
		else 
		{
		    SearchGain =True;
		    Poisson_Flag = True;
		    Gaussian_Flag = False;
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

     /* make sure there are not too many parameters */
     if (OptInd < argc)
      {
	  fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
	  exit(-1);
      }
      
      if ((OptPlan == True) && (Method != METHOD_MR_SUPPORT))
      {
        fprintf(OUTMAN, "Error: -n option is not valid with this method ...\n");
	exit(-1);
      }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
      
 }

/***************************************************************************/

int main(int argc, char *argv[])
{
     Ifloat DataB;

     /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    sigminit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       if (Gaussian_Flag == True)
    	  cout << "Method = " << StringMethod(Method) << endl;
       else
       {
           if (SearchReadout == True)
	    cout << "estimation of the standard deviation of the gaussian noise" << endl;
          else cout << "estimation of the gain" << endl;
        }
    }

    io_read_ima_float(Name_Imag_In, DataB);

    if (Gaussian_Flag == True) 
    {
       Noise_Ima = get_noise(DataB, Method, DEFAULT_MAX_ITER_METHOD, Nbr_Plan);
       cout << "\nNoise standard deviation =  "<< Noise_Ima << "\n\n" << endl;
    } 
    else if (Poisson_Flag == True) 
    {
       if (SearchReadout == True) 
       {
          Noise_Ima = get_readout_noise(DataB, Gain, Method);
          cout << "\nReadout Noise standard deviation =  "<< Noise_Ima << "\n\n" << endl;
       }
       else if (SearchGain == True) 
       {
          Noise_Ima = get_gain(DataB, ReadOutNoise, Method);
          cout << "\nGain =  "<< Noise_Ima << "\n\n" << endl;
       }
   }
   exit(0);
}



