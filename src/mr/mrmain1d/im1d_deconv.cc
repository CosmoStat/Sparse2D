/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  04/03/03
**    
**    File:  im1d_deconv.cc
**
*******************************************************************************
**
**    DESCRIPTION  deconvolution of an 1d signal by the standard deconvolution method 
**    ----------- 
**                 
**
******************************************************************************/
 
//#include "IM_Obj.h"
//#include "IM_IO.h"
//#include "IM_Sigma.h"
//#include "IM_Deconv.h"
//#include "IM_Errors.h" 
//#include "MREpsilon_Obj.h" 
//#include "MR_NoiseModel.h"
#include "MR1D_Deconv.h"

char Name_Signal_In[256];    // input file Signale 
char Name_Psf_In[256];       // PSF 
char Name_Signal_Out[256];   // output file name 
char Name_Signal_Start[256]; // First guess input solution  
char Name_Signal_ICF[256];   // ICF file name
char Name_Signal_Model[256]; // Model file name

Bool WriteResi = False;         // write the residual 
Bool PsfMaxShift = True;        // shift the max of the PSF to the center of thesignal
char Name_Resi[256];            // residual file name 
 
Bool GaussConv=False;
Bool UseICF=False;
Bool UseGuess=False;

float RegulParam = 0.1;
float Fwhm = 0.;
float Converg = 1.;                      // convergence parameter 
float Epsilon=DEFAULT_EPSILON_DECONV;    // convergence parameter 
int Max_Iter=DEFAULT_MAX_ITER_DECONV;    // Maximum number of iteration 

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
Bool Verbose;
Bool Optim = False;

type_deconv DecMethod=DEFAULT_STD_DECONV;

Bool PositivIma = True;
Bool UseModel =False;

/*********************************************************************/


static void usage(char *argv[])
{
    
    fprintf(OUTMAN, "Usage: %s options in_signal in_psf out_signal\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();
    fprintf(OUTMAN, "         [-d type_of_deconvolution]\n");
    for (int i = 0; i < NBR_STD_DECONV; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,StringDeconv((type_deconv)i));
    fprintf(OUTMAN, "              default is %s\n", StringDeconv((type_deconv)  DecMethod));
    manline();
      
    max_iter_usage(Max_Iter);
    manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Suppress the positivity constraint.\n");
    manline();
            
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "              Regularization parameter \n");
    fprintf(OUTMAN, "              default is %f\n", RegulParam);
    manline();

    fprintf(OUTMAN, "         [-C ConvergParam]\n");
    fprintf(OUTMAN, "              Convergence parameter. \n");
    fprintf(OUTMAN, "              default is %f\n", Converg);
    manline();
    
    fprintf(OUTMAN, "         [-f ICF_Fwhm]\n");
    fprintf(OUTMAN, "              Intrinsic correlation function.\n");
    fprintf(OUTMAN, "              Fwhm = Full Width at Half Maximum.\n");
    manline();
    
    fprintf(OUTMAN, "         [-I ICF_FileName]\n");
    fprintf(OUTMAN, "              Intrinsic correlation function file.\n");
    manline();
    
    fprintf(OUTMAN, "         [-F First_Guess]\n");
    fprintf(OUTMAN, "              Input solution file.\n");
    manline();
    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "              Optimization.\n");
    manline();
    fprintf(OUTMAN, "         [-M Model_Signale]\n");
    fprintf(OUTMAN, "              Input model for MEM method.\n");
    manline();
    write_residual_usage();
    manline();
    converg_param_usage(Epsilon);
    manline();    
    psf_not_center_usage();
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
    Bool OptG=False;
    Bool OptC=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif  
   
    /* get options */
    while ((c = GetOpt(argc,argv,"d:f:i:e:r:vC:G:OPSI:F:zZ:M:")) != -1) 
    {
	switch (c) 
        {
	   case 'd':
		/* -d <type> type of deconvolution */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_STD_DECONV)) DecMethod = (type_deconv) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of deconvolution: %s\n", OptArg);
	            exit(-1);
 		}
		break; 
          case 'f':
 		if (sscanf(OptArg,"%f",&Fwhm) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad ICF Fwhm: %s\n", OptArg);
		    exit(-1);
		}
		if (Fwhm > 0) GaussConv = True;
		if (UseICF == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
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
                if (Max_Iter <= 0) Max_Iter = DEFAULT_MAX_ITER_DECONV;
		break;
	   case 'C':
		/* -e < Convergence parameter> */
		if (sscanf(OptArg,"%f",&Converg) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad convergence parameter: %s\n", OptArg);
		   exit(-1);
		}
		OptC = True;
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
           case 'S': PsfMaxShift = False; break;
           case 'v': Verbose = True;break;
 	   case 'O' : Optim = True; break;
	   case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.1;
		OptG = True;
		break;	   
	  case 'r':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Resi) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteResi = True;
 		break;
	  case 'M':
		/* -r < residual file name> */
		if (sscanf(OptArg,"%s",Name_Signal_Model) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseModel = True;
 		break;
 	 case 'F':
 		if (sscanf(OptArg,"%s",Name_Signal_Start) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseGuess = True;
 		break; 
	  case 'I':
 		if (sscanf(OptArg,"%s",Name_Signal_ICF) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                UseICF = True;
	        if (GaussConv == True)
		{
		   fprintf(OUTMAN, "Error: -I and -f options are not compatible .. \n ");
		   exit(-1);
		}
 		break; 	   
            case 'P': PositivIma=False; break;  	   
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
       if ((UseModel == True) && (DecMethod != DEC_MEM_MODEL))
       {
          fprintf(OUTMAN, "Error: -M option is valid only with MEM method (Gull entropy) ...\n");
	  exit(-1);
       }
 
      
//        if ((Name_Psf_InOptG == True) && (DecMethod !=DEC_CLEAN)   
// 	     && (DecMethod != DEC_MEM_MODEL)  && (DecMethod != DEC_TIKHONOV))
//        {
//           fprintf(OUTMAN, "Error: -G option is not valid with this deconvolution method ...\n");
// 	  exit(-1);
//        }


	  
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Signal_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Psf_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Signal_Out, argv[OptInd++]);
         else usage(argv);

       if ((OptG == False) && (DecMethod !=DEC_CLEAN)   
	     && (DecMethod != DEC_MEM_MODEL)  && (DecMethod != DEC_TIKHONOV)
	     && (DecMethod != DEC_MARKOV)  && (DecMethod != DEC_MARKOV_LUCY)) {
           RegulParam = 0.;
       }
       if ((Optim == True) && (DecMethod != DEC_LUCY) && (DecMethod != DEC_TIKHONOV)
                      && (DecMethod != DEC_GRADIENT) && (DecMethod != DEC_MARKOV)
		      && (DecMethod != DEC_MARKOV_LUCY)) {
          fprintf(OUTMAN, "Error: -O option is not valid with this deconvolution method ...\n");
	  exit(-1);
       }
       if ((OptC == False) && (RegulParam > 1.))
                             Converg = 1. / (2. * RegulParam);
			     
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
	
	if ( (DecMethod == DEC_CLEAN) && (Max_Iter == DEFAULT_MAX_ITER_DECONV))
	   Max_Iter = DEFAULT_CLEAN_NITER;
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif	    
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    Ifloat Data;
    int k;
    fltarray Result;
    fltarray Resi;
    fltarray Psf, Psf1, Guess, Ima_ICF;
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR1);
    
     /* Get command line arguments, open input file(s) if necessary */
    decinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Signal_In << endl;
       cout << "File Name Out = " << Name_Signal_Out << endl;
       cout << "Method = " << StringDeconv((type_deconv) DecMethod) << endl;
       cout << "Epsilon = " << Epsilon << endl;
       cout << "Max_Iter = " << Max_Iter << endl;
       cout << "Cvg param iter = " << Converg << endl;
       cout << "RegulParam = " << RegulParam << endl;
       if (WriteResi == True)
          cout << "Signale residual file name : " << Name_Resi << endl;
     }

    /* read input Signale */
    MR1D_Deconv CDec;
    fits_read_fltarr (Name_Signal_In, CDec.Signal, &Header);
    fits_read_fltarr (Name_Psf_In, CDec.Psf);
cout <<  CDec.Signal.nx() << endl;  
cout <<  CDec.Psf.total() << endl;  
 
    if (UseGuess == True) fits_read_fltarr (Name_Signal_Start, Guess);
    if (UseICF == True) fits_read_fltarr (Name_Signal_ICF, Ima_ICF);
    if (UseModel == True) 
    {
       fits_read_fltarr (Name_Signal_Model, CDec.MemModel);
       if (CDec.MemModel.nx() != CDec.Signal.nx())
       {
          fprintf(OUTMAN, "Error: Signal and model must have the same size ...\n");
          exit(-1);
       }
    }
    //int Nx = Data.nx();
    Header.origin = Cmd;    
    CDec.UseModel = UseModel;
    CDec.PositivConstraint = PositivIma;
    CDec.DecMethod = DecMethod;
    CDec.PsfMaxShift = PsfMaxShift;
    // CDec.Noise_Ima = Noise_Ima;
    CDec.RegulParam = RegulParam;
    CDec.MaxIter = Max_Iter;
    CDec.EpsCvg = Epsilon;
    CDec.IterCvg = Converg;
    CDec.GaussConv = GaussConv;
    CDec.Fwhm = Fwhm;
    CDec.OptimParam = Optim;
    CDec.Verbose = Verbose;
    fltarray *Pt_G = NULL;
    if (UseGuess == True) Pt_G = &Guess;
    fltarray *Pt_ICF = NULL;
    if (UseICF == True) Pt_ICF = &Ima_ICF;
    
    //DECONVOLUTION
    if (Verbose == TRUE) cout << " Start the deconvolution ... " << endl;
    CDec.sig_deconv(Pt_G, Pt_ICF);

    fits_write_fltarr (Name_Signal_Out, CDec.Obj, &Header);
    if (WriteResi == True) fits_write_fltarr (Name_Resi, CDec.Resi, &Header);
    exit(0);
} 

