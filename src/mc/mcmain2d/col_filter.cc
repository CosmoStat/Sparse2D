/******************************************************************************
**                   Copyright (C) 2000 by CEA
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
**    File:  col_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  color image filtering program
**    ----------- 
**                 
******************************************************************************/
   
#include "Array.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_IO.h"
#include "SB_Filter.h"
#include "IM_Sigma.h"
#include "ColRest.h"

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

#define MAX_SCALE 10
#define  DEFAULT_N_SIGMA 3

Bool Verbose=False;
int Nbr_Plan = 4;
int Nbr_UndecScale = -1;
float N_Sigma = DEFAULT_N_SIGMA;
float Noise_Ima = 0;
Bool HardThreshold = True;
Bool MAD = False;
Bool ColSat  = True;


// http://www.webartz.com/fourcc/fccyvrgb.htm
// color information

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_catalogue result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
     
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", Nbr_Plan);    
    manline();

    fprintf(OUTMAN, "         [-u number_of_undecimated_scales]\n");
    fprintf(OUTMAN, "             Number of undecimated scales used in the  wavelet transform.\n");
    fprintf(OUTMAN, "             default is all. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-s nsigma]\n");
    fprintf(OUTMAN, "             HardThres = nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             SoftThres = 0.5 * nsigma * SigmaNoise\n");
    fprintf(OUTMAN, "             default is 3.\n");
    manline();

    fprintf(OUTMAN, "         [-S]\n");
    fprintf(OUTMAN, "             Use soft thresholding instead of hard thresholding.\n");
    fprintf(OUTMAN, "             default is False.\n");
    manline();

    fprintf(OUTMAN, "         [-g sigma]\n");
    fprintf(OUTMAN, "             sigma = noise standard deviation\n");
    fprintf(OUTMAN, "             default is automatically estimated.\n");
    manline();
    
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "             Correlated noise.\n");
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
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"TCSu:n:s:g:vzZ")) != -1) 
    {
	switch (c) 
        {  
	   case 'T': ColSat = False; break;
	   case 'C': MAD = True; break;
           case 'S': HardThreshold = False;
                 break;
          case 'u':
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%d",&Nbr_UndecScale) != 1) 
                {
                    fprintf(OUTMAN, "bad number of undecimated scales: %s\n", OptArg);
                    exit(-1);
                }
                if ((Nbr_UndecScale < 0) || (Nbr_Plan > MAX_SCALE)) 
                {
                    fprintf(OUTMAN, "bad number undecimated of scales: %s\n", OptArg);
                    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE);
                    exit(-1);
                }
                break;
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
 
/*********************************************************************/

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
        cout << "Number of scales = " << Nbr_Plan << endl;
	if (Nbr_UndecScale >= 0)
	   cout << "Number of undecimated scales = " << Nbr_UndecScale << endl;
        if (HardThreshold == True) cout << "Filter = Hard thresholding "  << endl;
        else cout << "Filter = Soft thresholding "  << endl;
        if ((MAD == False) && (Noise_Ima > FLOAT_EPSILON))
                cout << "Sigma Noise = " << Noise_Ima << endl;
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
    ColorRestore ColRest;
    ColRest.NbrUndec = Nbr_UndecScale; 
    // cout << "ColRest.NbrUndec = " <<  ColRest.NbrUndec << endl;
    ColRest.HardThreshold  = HardThreshold;
    ColRest.MAD = MAD;
    ColRest.Verbose = Verbose;
    ColRest.col_filter(Data, Nbr_Plan, N_Sigma, Noise_Ima);
    if (ColSat == True) ColRest.rescale_0_255(Data);
    io_3d_write_data(Name_Imag_Out , Data);
    exit(0);
}
