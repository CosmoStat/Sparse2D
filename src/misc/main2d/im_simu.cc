/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  im_simu.cc
**
*******************************************************************************
**
**    DESCRIPTION  transform an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: im_simu option image output
**        where options = 
**
**              [-p ]
**                Poisson noise
**                default is no
**
**              [-g sigma]
**                  Gaussian noise
**                  sigma = noise standard deviation 
**                  by default, the standard  deviation is set to 10. 
**
**		[-r PSF image]
**		    default is no
**
**		[-f full width at half-maximum]
**		     the instrumental response is gaussian and has the 
**		      same mean as image
**		      default is no
**
**
**		r and f options can't be used together
**              g and p options can't be used with r or f options
**
**
***************************************************************************/

// static char sccsid[] = "@(#)im_simu.cc 3.1 96/05/02 CEA 1996 @(#)";

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_IO.h"
#include "IM_Simu.h"
#include "IM_Deconv.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
float Noise_Ima = DEFAULT_NOISE_IMA; /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
Bool Instru_Resp = False; /* used for reading options; set to 1 if psf or f options are encountered */ 
Bool Gauss_Instru_Resp = False;
Bool Add_Noise = False;
char Name_IR_Image[80]; /* Instrumental Response */
float Gauss_IR_Width=0.0; /* width of the half-heigth gaussian instrumental response */
char Name_PSF_Out[256]; /* output file name */
Bool WritePSF=False;
float Gain=1.;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose = False;
unsigned int InitRnd = 100;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN, "\n");
    
    imsimu_option_usage();
    manline();
    vm_usage();
    manline();
    verbose_usage();
    manline(); 
    manline(); 
    fprintf(OUTMAN, "         r anf f options cannot be used together\n");
    fprintf(OUTMAN, "         g,p, and c options cannot be used together\n");
    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "\n");   
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void siminit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
   
    /* get options */
    while ((c = GetOpt(argc,argv,"p:g:r:f:c:w:G:vzZ:I:")) != -1) 
    {
	switch (c) 
        {
 	   case 'v': Verbose = True; break;
             case 'w': 
                if (sscanf(OptArg,"%s",Name_PSF_Out) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WritePSF = True;
                break;      
	     case 'p':
                /* Poisson noise */
                if (Add_Noise == False) Stat_Noise = NOISE_POISSON;
                else {
		    fprintf(OUTMAN, "\n\nError: only one type of noise can be selected: %s\n", OptArg);
		    exit(-1);
		}
                Noise_Ima = 1.0;
		Add_Noise = True;
		OptInd--;
               break;
	     case 'I':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		    exit(-1);
		}
                InitRnd = (unsigned int) c;
		break;
	     case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                if (Add_Noise == False) Stat_Noise = NOISE_GAUSSIAN;
                else {
		    fprintf(OUTMAN, "\n\nError: only one type of noise can be selected: %s\n", OptArg);
		    exit(-1);
		}
		Add_Noise = True;
		break;
		case 'G':
 		if (sscanf(OptArg,"%f",&Gain) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad gain: %s\n", OptArg);
		    exit(-1);
		}
		break;
	     case 'c':
		/* -c <sigma_noise> */
		if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
		    fprintf(OUTMAN, "\n\nError: bad sigma noise: %s\n", OptArg);
		    exit(-1);
		}
                if (Add_Noise == False) Stat_Noise = NOISE_GAUSS_POISSON;
                else {
		    fprintf(OUTMAN, "\n\nError: only one type of noise can be selected: %s\n", OptArg);
		    exit(-1);
		}
		Add_Noise = True;
		Stat_Noise = NOISE_GAUSS_POISSON;
		break;
	     case 'r':
		if (sscanf(OptArg,"%s",Name_IR_Image) != 1) 
                {
		    fprintf(OUTMAN, "\n\nError: bad PSF parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (Instru_Resp == False) Instru_Resp = True;
                else {
		    fprintf(OUTMAN, "\n\nError: -r and -f cannot be selected together: %s\n", OptArg);
		    exit(-1);
		}
		break;
	     case 'f':
		if (sscanf(OptArg,"%f",&Gauss_IR_Width) != 1) 
                {
		    fprintf(OUTMAN, "\n\nError: bad FWHM: %s\n", OptArg);
		    exit(-1);
		}
                if (Instru_Resp == False) Instru_Resp = True;
                else {
		    fprintf(OUTMAN, "\n\nError: -r and -f cannot be selected together: %s\n", OptArg);
		    exit(-1);
		}
                Gauss_Instru_Resp = True;
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

      if ((Instru_Resp == False) && (Add_Noise == False))
      {
          fprintf(OUTMAN, "\n\nError: noise addition option (pgc) or (and) convolution option has to selected \n");
          exit(-1);
      }
      
      if ((WritePSF == True) && (Gauss_Instru_Resp == False))
      {
          fprintf(OUTMAN, "\n\nError: -w is valid only if -f is set. \n");
          exit(-1);
      }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***********************************************************************/

void add_noise_poisson(Ifloat &DataB, Ifloat & Result, float Gain)
{
  int Nc,Nl;
  int i,j;
  int number=-1;


  Nc = Result.nc();
  Nl = Result.nl();

  for (i=0;i<Nl;i++)
    for (j=0;j<Nc;j++)
   	   Result(i,j) = poidev(DataB(i,j),&number) * Gain;
}

/***********************************************************************/

int main(int argc, char *argv[])
{
    int k,Nl,Nc;
    Ifloat DataB, Ir;
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    siminit(argc, argv);
    io_read_ima_float(Name_Imag_In, DataB, &Header);

    if (Verbose == True )
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
       if (Stat_Noise == NOISE_GAUSSIAN) 
               cout << "Type of Noise = GAUSSIAN" << endl;
       else if (Stat_Noise == NOISE_POISSON)  
               cout << "Type of Noise = POISSON" << endl;
       else
         cout << "Type of Noise = POISSON + GAUSSIEN" << endl;
       if (Noise_Ima > 0)
          cout << "Sigma Noise= " << Noise_Ima << endl;
       if ((Instru_Resp == True) && (Gauss_Instru_Resp == False))
	cout << "Name of Instrumental Response= " << Name_IR_Image << endl;
       else if ((Instru_Resp == True) && (Gauss_Instru_Resp == True))
	cout << "Gaussian Instrumental Response with width= " << Gauss_IR_Width << endl;
    }

    io_read_ima_float(Name_Imag_In, DataB, &Header);
    Nl = DataB.nl();
    Nc = DataB.nc();
    Header.origin = Cmd;
 
    Ifloat Result (DataB.nl(), DataB.nc(), "Result Transform");
    if (Instru_Resp == True)
    {
        /* we need square images */
//        dec_line_column (Nl, Nl1);
//        dec_line_column (Nc, Nc1);
//        if (Nl1 < Nc1) Nl1 = Nc1;
//        else Nc1 = Nl1;
//        if ((Nl != Nl1) || (Nc != Nc1)) 
//        {
//         fprintf(OUTMAN,"\n\nError: input image dimensions have to be power of 2\n");
//         exit(0);
//        }

       if (Gauss_Instru_Resp == False) 
       {
           io_read_ima_float(Name_IR_Image, Ir);
           // fft2d_conv(DataB, Ir, Result);
	   psf_convol (DataB, Ir, Result);
       }
       else 
       {
   	   Ir = im_gaussian(DataB.nl(),DataB.nc(), Gauss_IR_Width);
	   norm_flux(Ir);    
	   if (WritePSF == True) io_write_ima_float(Name_PSF_Out, Ir);
	   // fft2d_conv(DataB, Result, Result);
	   psf_convol (DataB, Ir, Result);
       }
       DataB = Result;
    }
 
    if (Add_Noise == True)
    {
    	switch (Stat_Noise)
    	{
    	   case NOISE_GAUSSIAN:
      		/* add a gaussian noise */
		im_noise_gaussian(Result, Noise_Ima, InitRnd);
     		Result += DataB;
      		break;
    	   case NOISE_POISSON:
      		/* add a poisson noise */
 		Result = DataB;
		im_noise_poisson (Result, 1., InitRnd);
      		break;
    	   case NOISE_GAUSS_POISSON:
     	 	/* add a poisson noise and then a gaussian noise */
 		Result = DataB;
		im_noise_poisson (Result, Gain, InitRnd);
		im_noise_gaussian(DataB, Noise_Ima, InitRnd);
     		Result += DataB;
       		break;
    	   default:
     		cerr << "Error: other king of noise are not generated by this routine ... " << endl;
     		exit(-1);
    	}
    }
	
    /* write result in output image */
    io_write_ima_float(Name_Imag_Out, Result, &Header);
    exit(0);
}





