/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/01/02
**    
**    File:  mr_contrast.cc
**
*******************************************************************************
**
**    DESCRIPTION  contrast enhancement
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "FFTN_2D.h"
#include "IM_Math.h"
#include "MR_Contrast.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool Reverse=False; // Reverse transform
int NbrScale2D = DEFAULT_NBR_SCALE;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;

double QParam = 0.;
double CParam = 0.;  // must be the noise level
double PParam = 0.5; // P determined the degree of non lineraity
double MParam = 100.; // values larger than m are not modified

float ClipVal = 3.;
float LogParam = 0.1;
Bool PosIma = True;
float L_Coeff=1.;

type_contrast TContrast = DEF_CONTRAST;

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-m contrast_enhancement_method]\n");
    for (int i = 0; i < NBR_CONTRAST; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringContrast((type_contrast) (i)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringContrast(TContrast));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    manline();

    fprintf(OUTMAN, "         [-M M_parameter]\n");
    fprintf(OUTMAN, "             Coeff larger than M are not enhanced. \n");
    fprintf(OUTMAN, "             Default is %f.\n", MParam);
    manline();
    
    fprintf(OUTMAN, "         [-P P_parameter. P in ]0,1[)]\n");
    fprintf(OUTMAN, "             Determine the degree on non-linearity. \n");
    fprintf(OUTMAN, "             P must be in ]0,1[.  \n");
    fprintf(OUTMAN, "             Default is %f.\n", PParam);
    manline();
    
    fprintf(OUTMAN, "         [-Q Q_parameter]\n");
    fprintf(OUTMAN, "             Q > 0 ==> enhance less the darker part than \n");
    fprintf(OUTMAN, "                       the lighter part. \n"); 
    fprintf(OUTMAN, "             Q < 0 ==> enhance more the dark than \n");
    fprintf(OUTMAN, "                       the lighter part. \n");
    fprintf(OUTMAN, "             Q must be in [-0.5,0.5].  \n");
    fprintf(OUTMAN, "             Default is %f.\n", QParam);
    manline();
    
    fprintf(OUTMAN, "         [-C C_parameter]\n");
    fprintf(OUTMAN, "             Noise level. Default is %f.\n", CParam);
    manline(); 
    fprintf(OUTMAN, "         [-K ClippingValue]\n");
    fprintf(OUTMAN, "             Clipping value. Default is 3.\n");
    manline();

    fprintf(OUTMAN, "         [-G Param]\n");
    fprintf(OUTMAN, "             Parameter for for the Laplacian and \n");
    fprintf(OUTMAN, "             the log method.\n");
    fprintf(OUTMAN, "             Default is %f.\n", LogParam);
    manline();
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Remove positivity contraint.\n");
    manline();
    fprintf(OUTMAN, "         [-L Saturation]\n");
    fprintf(OUTMAN, "             Saturate the reconstructed image.\n");
    fprintf(OUTMAN, "             Coeff larger than L*Max are set to L*Max.\n");
    fprintf(OUTMAN, "             default is %5.2f. \n", L_Coeff);
    fprintf(OUTMAN, "             If L is set to 0, then no saturation is applied.\n");
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
    while ((c = GetOpt(argc,argv,"G:pm:L:K:M:Q:P:C:n:vzZ")) != -1) 
    {
	switch (c) 
        {  
	  case 'L': 
	        // if (sscanf(OptArg,"%f",&SParam) != 1) 
	        if (sscanf(OptArg,"%f",&L_Coeff) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad L_Coef param: %s\n", OptArg);
                    exit(-1);
                }

                if (L_Coeff < 0)
                {
                    fprintf(OUTMAN, "Error: bad L_Coef param: %s\n", OptArg);
                    fprintf(OUTMAN, "0 <= L_Coef <= 1  \n");
                    exit(-1);
                }
                break;
          case 'p': PosIma = False; break; 
          case 'm': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of contrast method: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_CONTRAST)) TContrast = (type_contrast) (c-1);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of contrast method: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'G': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&LogParam) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad Log param: %s\n", OptArg);
                    exit(-1);
                }
                if (LogParam <= 0)  
                {
                    fprintf(OUTMAN, "Error: bad Log param: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < MParam  \n");
                    exit(-1);
                }
                break;          
	 case 'K': 
               /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%f",&ClipVal) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad clipping param: %s\n", OptArg);
                    exit(-1);
                }
                if (ClipVal <= 0)  
                {
                    fprintf(OUTMAN, "Error: bad clipping param: %s\n", OptArg);
                    fprintf(OUTMAN, "0 < MParam  \n");
                    exit(-1);
                }
                break;
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
           case 'v': Verbose = True;break;
	   case 'n':
                /* -n <NbrScale2D> */
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
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
   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}
 
/*************************************************************************/
    
int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl; 
	cout << "Method = " << StringContrast(TContrast) << endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
    }
 
    Ifloat Data;
    fltarray Trans;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
 
    Noise_Ima = CParam;
    if (Noise_Ima < FLOAT_EPSILON)
    {
       Noise_Ima = detect_noise_from_bspline (Data);
                   //  detect_noise_from_med (Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    Contrast CT;
    CT.Noise_Ima = Noise_Ima;
    CT.Verbose = Verbose;    
    CT.Contrast_P_Param = PParam;
    CT.Contrast_Q_Param = QParam;
    CT.Contrast_C_Param = CParam;
    CT.Contrast_M_Param = MParam;
    CT.ClipVal = ClipVal;
    CT.LogParam =LogParam;
    CT.LaplacianParam=LogParam;
    int NbrBin = 255;
    float Max;
    if (L_Coeff != 0)  Max = max(Data);
    switch (TContrast)
    {
       case CT_LAPLACIAN:
            CT.laplacian(Data); break;
       case CT_ATROU_LOG_SGN:
            CT.atrou_log(Data, NbrScale2D, True); break;
       case CT_ATROU_LOG:
            CT.atrou_log(Data, NbrScale2D, False); break;
       case CT_HISTO:
            CT.histo_equal(Data, NbrBin);break;
       case CT_RET:
            CT.retinex(Data); break;
       case CT_MRET:
            CT.multiscale_retinex(Data); break;
       case CT_ATROU_RET:
             CT.atrou_retinex(Data, NbrScale2D); break;
       case CT_WT:
             CT.wt_enhance(Data,NbrScale2D); break;
       case CT_CLIPPING:
            CT.sature_clipping(Data); break;
       default:
            cout << "Error: not implemented method ... " << endl;
	    exit(-1);
    }
    if (PosIma == True) threshold(Data);
    if (L_Coeff != 0)
    {
       Max *= L_Coeff;
       for (int i = 0; i < Data.nl(); i++)
       for (int j = 0; j < Data.nc(); j++)
       {
           if (Data(i,j) > Max) Data(i,j) = Max;
       }
    }     
    io_write_ima_float(Name_Imag_Out, Data, &Header);
    exit(0);
}
