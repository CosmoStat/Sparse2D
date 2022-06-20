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
**    Date:  98/01/12
**    
**    File:  mr_edge.cc
**
*******************************************************************************
**
**    DESCRIPTION  detect the edges in an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_edge.cc 1.0 98/01/12 CEA 1998 @(#)";
  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "MR_Abaque.h"
#include "MR_Psupport.h"
#include "IM_Edge.h"
#include "MR_Edge.h"


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Write_Sup[256];          /* output support file name */
int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
float N_Sigma=DEFAULT_N_SIGMA;   /* number of sigma (for the noise) */
float Noise_Ima=0.;              /* noise standard deviation */
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   /* type of noise */
type_transform Transform = TO_PAVE_BSPLINE; 
Bool WriteSup = False;          /* write the support on the disk */
Bool UseNSigma =False;
 
char Name_RMSMap[256]; 
Bool UseRMSMap=False;

float EpsilonPoisson = DEFAULT_EPSILON;
int SizeBlock = DEFAULT_SIZE_BLOCK_SIG;
int NiterClip = 1;
 
extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool EpsOpt=False;
Bool MaxIterOpt=False;
Bool PosOpt=False;
Bool KillLastOpt = False;
Bool WindowOpt=False;    

Bool Verbose=False;
float Threshold = 0.;
Bool NoNoise=True;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
 //    gauss_usage();
//     manline();
//     
//     ccd_usage();
//     manline();
//     
//     noise_usage();
//     manline();

    nbr_scale_usage(Nbr_Plan);
//     nbr_scalep_usage(DEF_N_SCALE);
//     manline();
//     
//     nsigma_usage(N_Sigma);
//     manline();
//  
//     prec_eps_poisson_usage(DEFAULT_EPSILON);
//     manline();
// 
//     size_block_usage(SizeBlock);
//     manline();
//     
//     sigma_clip_block_usage(NiterClip);
//     manline();
//   
//     rms_noise_usage();
//     manline();
//  
//     fprintf(OUTMAN, "         [-L]\n");
//     fprintf(OUTMAN, "             Do not apply a noise modeling.\n");
//     manline();
//     
//     fprintf(OUTMAN, "         [-t MultiresolutionMethod]\n");
//     fprintf(OUTMAN, "             1: a trous algorithm + zero crossing.\n");
//     fprintf(OUTMAN, "             2: Wavelet transform modulus maxima.\n");
//     fprintf(OUTMAN, "             Default is 2.\n");
//     
//     manline();
//     fprintf(OUTMAN, "         [-w Reconstructed Edge Map]\n");
//     fprintf(OUTMAN, "             Reconstruct an edge map from the multiscale edges.\n");
//  
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
static void filtinit(int argc, char *argv[])
{
    int c;
    Bool NscaleOpt=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     
    
    /* get options */
//    while ((c = GetOpt(argc,argv,"m:g:c:n:s:e:w:E:S:N:R:t:LvzZ:")) != -1) 
    while ((c = GetOpt(argc,argv,"vn:")) != -1)
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
	    case 'L': NoNoise = False;break;
	    case 't':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad transform: %s\n", OptArg);
		    exit(-1);
		}
		if (c == 1) Transform = TO_PAVE_BSPLINE; 
		else if (c== 2) Transform = TO_DIADIC_MALLAT;
		else
		{
		    fprintf(OUTMAN, "Error: bad transform: %s\n", OptArg);
		    exit(-1);
		}
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
		NscaleOpt = True;
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
 	   case 'w':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s", Name_Write_Sup) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WriteSup = True;
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
		break; 
	 case 'S':
		if (sscanf(OptArg,"%d",&SizeBlock) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad block size: %s\n", OptArg);
		    exit(-1);
		}
		if (SizeBlock  < 2)
		{
		   fprintf(OUTMAN, "Error: bad  SizeBlock parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "              SizeBlock > 1\n");
		   exit(-1);
		}
		break; 
	 case 'N':
		if (sscanf(OptArg,"%d",&NiterClip) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad it. number for the 3sigma clipping: %s\n", OptArg);
		    exit(-1);
		}
		if (NiterClip < 1)
		{
		   fprintf(OUTMAN, "Error: bad NiterClip parameter: %s\n", OptArg);
		   fprintf(OUTMAN, "             NiterClip > 0\n");
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
 	   Transform = TO_PAVE_BSPLINE;
 	   if (NscaleOpt != True) Nbr_Plan = DEF_N_SCALE;
  	}

 	if (UseRMSMap == True)
 	{
 	   if ((Stat_Noise != DEFAULT_STAT_NOISE) &&
               (Stat_Noise != NOISE_NON_UNI_ADD))
 	   {
 	      cerr << "WARNING: noise model is automatically set to : " << endl;
 	      cerr << "        " << StringNoise(NOISE_NON_UNI_ADD) << endl;
 	      cerr << "         when RMS map option is set ... " << endl;
 	   }
 	   Stat_Noise = NOISE_NON_UNI_ADD;
 	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}


/****************************************************************************/

int main(int argc, char *argv[])
{
    int s,i,k;
    Ifloat Data;
    
    /* support image creation */
    fitsstruct Header;
    char Cmd[256];
 
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    lm_check(LIC_MR1);
     /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

 

if (Verbose == True)
{
    cout << endl << endl << "PARAMETERS: " << endl << endl;
    cout << "File Name in = " << Name_Imag_In << endl;
    cout << "File Name Out = " << Name_Imag_Out << endl;
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
    if ((Stat_Noise == NOISE_GAUSSIAN) && (Noise_Ima > FLOAT_EPSILON))
                cout << "Sigma Noise = " << Noise_Ima << endl;
    if (Stat_Noise ==  NOISE_EVENT_POISSON)
    {
       cout << "Epsilon Poisson = " <<  EpsilonPoisson << endl;
    }
    else cout << "N_Sigma = " << N_Sigma << endl;
       
    if (WriteSup == True)
       cout << "MR Edge file name : " << Name_Write_Sup << endl;
}

    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    int Nl = Data.nl();
    int Nc = Data.nc();
 
    // noise model class initialization
    MultiResol MR_Data (Nl,Nc, Nbr_Plan, Transform, "MR_Transform");
    MR_Data.Border = I_MIRROR;
       		 
    if (NoNoise == True)
    {
        MR_Data.transform(Data);
        mr_zero_cross_edge(MR_Data);
	MR_Data.write(Name_Imag_Out);
    }
    else 
    {
       MRNoiseModel ModelData(Stat_Noise, Nl,Nc, 
                           Nbr_Plan, Transform);    

       Ifloat ImaEdge(Nl,Nc, "Result Edges");
       Iint TCont(Nl,Nc, "Result type edges");
       // Compute the thresholded wavelet coefficents  
       MultiResol MR_Edge (Nl,Nc, Nbr_Plan, TO_PAVE_BSPLINE, "MR_Transform");        if (Noise_Ima > FLOAT_EPSILON) ModelData.SigmaNoise = Noise_Ima;

       if (UseNSigma  == True)
          for (i=0; i < Nbr_Plan; i++) ModelData.NSigma[i]=N_Sigma;
       for (s=0; s < Nbr_Plan; s++) ModelData.TabEps[s] = EpsilonPoisson;
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
      
      ModelData.model(Data, MR_Data);
      // if (UseNSigma == True) ModelData.threshold(MR_Data);

      mr_get_edge(MR_Data, MR_Edge, ImaEdge, TCont,  ModelData);
      if (UseNSigma == True) ModelData.threshold(MR_Edge);
      if (WriteSup == True) 
	         io_write_ima_float(Name_Imag_Out, ImaEdge);
      MR_Edge.write(Name_Imag_Out);
    }
    exit(0);
}
 
 
 
/*  PCA test
    cout << "Start edge detection ... " << MR_Data.nbr_band() << endl;
    if (Threshold > 0) ModelData.threshold(MR_Data);
 
   MR_Edge.write("edge.mr");
   int N_Image = Nbr_Plan-1;
   CorrelMatAna2D CMA(N_Image);
   CMA.compute(MR_Edge.tabband());
   CMA.print();
   CMA.transform(MR_Edge.tabband(), MR_Edge.tabband());
       
   // if (WriteEigen == True)
   {
      char NameEigen[256];
      for(k=0; k < N_Image; k++)
      {
         sprintf(NameEigen, "pca_%d",k+1);
         io_write_ima_float(NameEigen, MR_Edge.band(k));
      }
   }   
*/
