/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  16/03/2000
**    
**    File:  col_contrast.cc
**
*******************************************************************************
**
**    DESCRIPTION  color image enhancement program
**    ----------- 
**                 
******************************************************************************/
   
#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "SB_Filter.h"
#include "IM_Sigma.h"
#include "MR_Obj.h"
#include "ColRest.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

// #define MAX_SCALE 10
// #define  DEFAULT_N_SIGMA 3

Bool Verbose=False;
int Nbr_Plan = 4;
int Nbr_UndecScale = -1;
float N_Sigma = DEFAULT_N_SIGMA;
float Noise_Ima = 0;
Bool HardThreshold = True;
Bool MAD = False;

double QParam = 0.;
double CParam = 0.;  // must be the noise level
double PParam = 0.5; // P determined the degree of non lineraity
double MParam = 100.; // values larger than m are not modified

Bool UseFilter = False;    
Bool HistoEqual = False; 
Bool ImproveEdge = False;
Bool Retinex = False;
Bool MultiRetinex = False;
Bool Clip = True;
Bool ColSat = True;
float ClipVal=3.;
Bool AtrouRet = False;
Bool Norm = False;
Bool AtrouLog = False;
// http://www.webartz.com/fourcc/fccyvrgb.htm
// color information


/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
     
    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Performs a filtering before the enhancement.\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();     
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", Nbr_Plan);    
    manline();

    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "             Only used when is filtering is performed. \n");
    fprintf(OUTMAN, "             HardThres = nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             SoftThres = 0.5 * nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             default is 3.\n");
    manline();

    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             Only used when is filtering is performed.\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "             default is automatically estimated.\n");
    manline();

    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "             Color saturation. Default is true.\n");
    manline();
    
    fprintf(OUTMAN, "         [-c]\n");
    fprintf(OUTMAN, "             RGB clipping. Default is true.\n");
    manline();
        
    fprintf(OUTMAN, "         [-h]\n");
    fprintf(OUTMAN, "             Histogram equalization of the L component.\n");
    fprintf(OUTMAN, "             Default is no.\n");
    manline();
    
    fprintf(OUTMAN, "         [-e]\n");
    fprintf(OUTMAN, "             Multiscale Edge enhancement. Default is no.\n");
    manline();
    
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Retinex method. Default is no.\n");
    manline();

    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "             Multiscale Retinex method. Default is no.\n");
    manline();
    
    fprintf(OUTMAN, "         [-A]\n");
    fprintf(OUTMAN, "             A trous Multiscale Retinex method. Default is no.\n");
    manline();
            
    fprintf(OUTMAN, "         [-M M_parameter]\n");
    fprintf(OUTMAN, "             M Paramter. Default is %f.\n", MParam);
    manline();
    
    fprintf(OUTMAN, "         [-P P_parameter]\n");
    fprintf(OUTMAN, "             P Paramter. Default is %f.\n", PParam);
    manline();
    
    fprintf(OUTMAN, "         [-Q Q_parameter]\n");
    fprintf(OUTMAN, "             Q Paramter. Default is %f.\n", QParam);
    manline();
    
    fprintf(OUTMAN, "         [-C C_parameter]\n");
    fprintf(OUTMAN, "             C Paramter. Default is %f.\n", CParam);
    manline();
    
    fprintf(OUTMAN, "         [-K ClippingValue]\n");
    fprintf(OUTMAN, "             Clipping value. Default is 3.\n");
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
    float Val;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"AK:cSRrfheC:M:P:Q:fCSu:n:s:g:vzZ")) != -1) 
    {
	switch (c) 
        {  
 	  case 'A': AtrouRet = True; break;
	  case 'K': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&ClipVal) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad clipping param: %s\n", OptArg);
                    exit(-1);
                }
                if (MParam <= 0)  
                {
                    fprintf(OUTMAN, "Error: bad clipping param: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < MParam  \n");
                    exit(-1);
                }
                break;
	   case 'c': Clip = False; break;
	   case 'S': ColSat = False; break;
	   case 'r': Retinex = True; break;
	   case 'R': MultiRetinex = True; break;
	   case 'M': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad MParam: %s\n", OptArg);
                    exit(-1);
                }
		MParam = Val;
                if (MParam <= 0)  
                {
                    fprintf(OUTMAN, "Error: bad MParam: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < MParam  \n");
                    exit(-1);
                }
                break;
	   case 'C': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad CParam: %s\n", OptArg);
                    exit(-1);
                }
		CParam = Val;
                if (MParam <= 0)  
                {
                    fprintf(OUTMAN, "Error: bad CParam: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < CParam  \n");
                    exit(-1);
                }
                break;
	   case 'Q': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad QParam: %s\n", OptArg);
                    exit(-1);
                }
		QParam = Val;
                if ((QParam < -0.5) || (QParam > 0.5)) 
                {
                    fprintf(OUTMAN, "Error: bad QParam: %s\n", OptArg);
                    fprintf(OUTMAN, "-0.5 <= QParam <= 0.5 \n");
                    exit(-1);
                }
                break;
	   case 'P': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&Val) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad PParam: %s\n", OptArg);
                    exit(-1);
                }
		PParam = Val;
                if ((PParam < 0) || (PParam > 1)) 
                {
                    fprintf(OUTMAN, "Error: bad PParam: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < PParam < 1 \n");
                    exit(-1);
                }
                break;
           case 'e': ImproveEdge = True;; break;
	   case 'h': HistoEqual = True;; break;
           case 'f': UseFilter = True; break;
	   case 'n':
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE);
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
		CParam = Noise_Ima;
                break;
           case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
                    exit(-1);
                }
                if (N_Sigma < 0.)  N_Sigma = DEFAULT_N_SIGMA;
                break;
           case 'v': Verbose = True;break;
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

       
/***************************************************************************/

int main(int argc, char *argv[])
{
    int  k;
    char Cmd[512];
    extern softinfo Soft;

    Soft.mr3();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3); 
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;
	if ((UseFilter == True) || (ImproveEdge == True))
	{
          cout << "Number of scales = " << Nbr_Plan << endl;
	  if (Nbr_UndecScale >= 0)
	   cout << "Number of undecimated scales = " << Nbr_UndecScale << endl;
           if (HardThreshold == True) cout << "Filter = Hard thresholding "  << endl;
           else cout << "Filter = Soft thresholding "  << endl;
          if ((MAD == False) && (Noise_Ima > FLOAT_EPSILON))
                cout << "Sigma Noise = " << Noise_Ima << endl;
       }
       if (ImproveEdge == True) cout << "Wavelet Edge Enhancement" << endl;
       if (Retinex == True) cout << "Retinex method" << endl;
       if (MultiRetinex == True) cout << "Multiscale Retinex method" << endl; 
       if (HistoEqual == True) cout << "Histogram equalization of the Luminace" << endl;
       if (AtrouRet == True)  cout << "A trous Wavelet Retinex method" << endl;
       if (Clip == True)  cout << "Sigma Clipping" << endl;
       if (ColSat == True) cout << "Rescale between 0 and 255" << endl;
    }    
    if (Nbr_UndecScale < 0) Nbr_UndecScale = Nbr_Plan;
    
    fltarray Data;
    io_3d_read_data(Name_Imag_In, Data);
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    if (Nz != 3)
    {
        cout << "Error: input data are not a multichannel color image ... " << endl;
        exit(-1);
    }
    if (Verbose == True)
    { 
       cout << "naxis = " << Data.naxis() << endl;
       cout << "Nx = " << Nx << " Ny = "<< Ny << " Nz = " << Nz << endl;
       cout << "min = " << Data.min()  << " max = " << Data.max();
       cout << " sigma = " << Data.sigma() << endl;
    }
    if (Verbose == True) cout << "rgb_to_yuv  ... " << endl;

    fltarray NormTab;    
 
    if (Verbose == True) cout << "Transform" << endl;
    ColorRestore ColRest;
    ColRest.MAD = MAD;
    ColRest.Contrast_P_Param = PParam;
    ColRest.Contrast_Q_Param = QParam;
    ColRest.Contrast_C_Param = CParam;
    ColRest.Contrast_M_Param = MParam;
    ColRest.ClipVal = ClipVal;
    ColRest.Verbose = Verbose;
    if (UseFilter == True) 
         ColRest.col_filter(Data, Nbr_Plan, N_Sigma, Noise_Ima);

    if (ImproveEdge == True) 
        ColRest.enhance(Data, Nbr_Plan, N_Sigma, Noise_Ima);
	    
    if (HistoEqual == True) ColRest.equalize_histo_luminance(Data);
       // ColRest.equalize_histo(Data,2);   
       	
    if (Retinex == True)      ColRest.retinex(Data);
    if (MultiRetinex == True) ColRest.multiscale_retinex(Data);
    if (AtrouRet == True) ColRest.atrou_retinex(Data, Nbr_Plan, N_Sigma, Noise_Ima);
    // if (AtrouLog == True) ColRest.atrou_log(Data,  Nbr_Plan, False, Noise_Ima, 0.1);
    // ColRest.info(Data);
    if (Clip == True) ColRest.sature_clipping(Data);
    if (ColSat == True) ColRest.rescale_0_255(Data);
    io_3d_write_data(Name_Imag_Out, Data);
    exit(0);
}
