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
**    Date:  03/02/00
**    
**    File:  mbase.cc
**
*******************************************************************************
**
**    DESCRIPTION  Decomposition of an image on multiple bases
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Usage.h"
#include "MRBase.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale1D = -1;
int NbrScale2D = 4;

float Noise_Ima=0.;
float N_Sigma=4.;
Bool PoissonNoise=False;

Bool PositivSol=True;
int BlockSize=0;
int Nbr_Iter = DEF_MB_NBR_ITER;
 Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
int FirstDetectScale=0;
Bool KillLastScale = False;
Bool ThresholdRecIma = True;

Bool Filtering = True;

Bool TabBase[NBR_MRBASE];
type_mbase TabSelect[NBR_MRBASE];
int NbrBase = 0;

int WT_NbrUndecimatedScale = -1;
int WP_NbrUndecimatedScale = 2;

Bool WriteAllRec=False;
Bool UseHuberNorm=False;
float FirstSoftThreshold = DEF_MB_FIRST_SOFT_THRESHOLD;
float LastSoftThreshold = DEF_MB_LAST_SOFT_THRESHOLD;
Bool RidKillLastScale = False;
Bool TV = False;
Bool UsePsf = False;
char PSFFile[256];
float Eps= 1.e-3;
unsigned int InitRnd = 100;
float DetectCoefTol = 0.5;
float RegulParam = 0.1;
int CosBlockSize=0;
Bool CosWeightFirst = False;

int FCurNdir=32;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t TransformSelection]\n");
    for (int i = 0; i < NBR_USED_MRBASE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMRBase ((type_mbase) i));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the WT, the a trous, the PMT and the curvelet transform.\n");
    fprintf(OUTMAN, "             default is 4.\n");    
    // fprintf(OUTMAN, "             Number of ridgelet in the multi-ridgelet transform.\n");
    manline();
 
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is image size. \n");    
    // fprintf(OUTMAN, "             Starting Block Size in the multi-ridgelet transform.\n"); 
    // fprintf(OUTMAN, "             default is 8. \n");    
    manline();

    fprintf(OUTMAN, "         [-B DCT_BlockSize]\n");
    fprintf(OUTMAN, "             Local DCT block size.\n");
    fprintf(OUTMAN, "             By default, a global DCT is used. \n");    
    manline();


    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration. Default is%d.\n", Nbr_Iter);    
    manline();
 
    fprintf(OUTMAN, "         [-F FirstDetectionScale]\n");
    fprintf(OUTMAN, "             First detection scale in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is 1. \n");    
    manline();

//     fprintf(OUTMAN, "         [-k]\n");
//     fprintf(OUTMAN, "            Kill last scale in ridgelet, and multiscale ridgelet transform.\n");
//     fprintf(OUTMAN, "            default no.\n");    
//     manline();

    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "            Kill the last scale in the a trous algorithm and the curvelet.\n");
//    fprintf(OUTMAN, "            Kill the last scale in the a trous algorithm, the PMT, and the curvelet.\n");
    fprintf(OUTMAN, "            default no.\n");    
    manline();

    fprintf(OUTMAN, "         [-L FirstSoftThreshold]\n");
    fprintf(OUTMAN, "            First soft thresholding value.\n");
    fprintf(OUTMAN, "            default is %f.\n", FirstSoftThreshold);    
    manline();

    fprintf(OUTMAN, "         [-l LastSoftThreshold]\n");
    fprintf(OUTMAN, "            Last soft thresholding value..\n");
    fprintf(OUTMAN, "            default is %f.\n", LastSoftThreshold);    
    manline();

    fprintf(OUTMAN, "         [-u]\n");
    fprintf(OUTMAN, "             Number of undecimated scales in the WT.\n");
    fprintf(OUTMAN, "             default is all scales. \n");    
    manline();

    nsigma_usage(N_Sigma);
    manline();

    gauss_usage();
    manline();  
     
//    poisson_noise_usage();
//    manline(); 

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Supress the block overlapping. Default is no. \n");    
    manline();

    // fprintf(OUTMAN, "         [-P]\n");
    // fprintf(OUTMAN, "             Supress the positivity constraint. Default is no. \n");    
    // manline();
    
    fprintf(OUTMAN, "         [-G RegulParam]\n");
    fprintf(OUTMAN, "             Minimize the Total Variation instead of the L1 norm. \n"); 
    fprintf(OUTMAN, "             Default Regularization parameter for the TV method is %f.\n", RegulParam);
    manline();
    
    fprintf(OUTMAN, "         [-P PsfFile]\n");
    fprintf(OUTMAN, "             Apply a deconvolution using the PSF in the file PsfFile. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-e Eps]\n");
    fprintf(OUTMAN, "             Remove frequencies with |PP*| < Eps. \n");
    fprintf(OUTMAN, "             Default is 1e-3.\n");    
    manline();

    fprintf(OUTMAN, "         [-C TolCoef]\n");
    fprintf(OUTMAN, "              Default is 0.5. \n");    
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
    while ((c = GetOpt(argc,argv,"xB:G:C:I:e:P:Tl:L:u:t:KF:Oi:b:s:g:n:pvzZ")) != -1) 
    {
	switch (c)
        {  
	   case 'x': CosWeightFirst  = True;  break;
           case 'B':
                if (sscanf(OptArg,"%d",&CosBlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad cosinus transform block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'G':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
		    fprintf(stderr, "bad Regularization Parameter: %s\n", OptArg);
		    exit(-1);
		}
                if (RegulParam  < 0.) RegulParam = 0.1;
          case 'C':
	        if (sscanf(OptArg,"%f",&DetectCoefTol) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad tolerance coefficient: %s\n", OptArg);
                    exit(-1);
                }
                if (DetectCoefTol < 0)
                {
                   fprintf(OUTMAN, "Error: tolerance coefficient: %s\n", OptArg);
                   exit(-1);
                }
                break;
         case 'I':
 		if (sscanf(OptArg,"%d",&c) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number value: %s\n", OptArg);
		    exit(-1);
		}
                InitRnd = (unsigned int) c;
		break;	
         case 'e': /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&Eps) != 1) {
                fprintf(OUTMAN, "Error: bad Eps: %s\n", OptArg);
                exit(-1);
                }
                if (Eps <= 0.)  Eps = 0.001;
                break;
         case 'P':
	        if (sscanf(OptArg,"%s",PSFFile) != 1) 
                {
                    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
                    exit(-1);
                }
		UsePsf = True;
	        break;	
          case 'T': TV = (TV == True) ? False: True; break;
          case 'k': RidKillLastScale = True; break;
          case 'L':               
                if (sscanf(OptArg,"%f",&FirstSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0)
                {
                   fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                   exit(-1);
                }
                break;
          case 'l':               
                if (sscanf(OptArg,"%f",&LastSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold < 0)
                {
                   fprintf(OUTMAN, "Error: bad last soft threshold: %s\n", OptArg);
                   exit(-1);
                }
                break;
           // case 'P': PositivSol = False; break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_MRBASE))
                {
                   if (TabBase[c-1] == True)
                   {
                       fprintf(OUTMAN, "Error: transform already selected: %s\n", OptArg);
                       exit(-1);
                   }
                   else
                   {
                      TabBase[c-1] = True;
                      TabSelect[NbrBase++] = (type_mbase) (c-1);
                   }
                }
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'f': Filtering = (Filtering == True) ? False: True; break;
           case 'K': KillLastScale = True; break;
           case 'u': 
                if (sscanf(OptArg,"%d",&WT_NbrUndecimatedScale) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                    exit(-1);
                }
                if (WT_NbrUndecimatedScale < 0)  
                {
                   fprintf(OUTMAN, "Error: bad scale number: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'F': 
                if (sscanf(OptArg,"%d",&FirstDetectScale ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first detection scale: %s\n", OptArg);
                    exit(-1);
                }
                FirstDetectScale --;
                if (FirstDetectScale < 0)  
                {
                   fprintf(OUTMAN, "Error: bad first detection scale: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
           case 'i':
                if (sscanf(OptArg,"%d",&Nbr_Iter) != 1) 
                {
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'p': PoissonNoise=True;
                     Noise_Ima = 1.;
                     break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
		// FCurNdir = BlockSize;
                break;
           case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
          case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
                    usage(argv);
                }
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
                break;
           case 'v': Verbose = True;break;
	   case 'n':
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D  > MAX_SCALE)
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

/***************************************************************************/

int main(int argc, char *argv[])
{
    int i,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    for (i = 0; i < NBR_MRBASE; i++) TabBase[i] = False;

    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        cout << "NbrScale2D = " << NbrScale2D << endl;
        if (FirstDetectScale > 0) cout << "FirstDetectScale = " << FirstDetectScale << endl;
        if (PoissonNoise == True)
           cout << "Poisson noise " << endl; 
        cout << "NSigma = " << N_Sigma  << endl; 
	if ( TV == True) cout << "TV regularization: lambda = " <<  RegulParam << endl;
	else cout << "L_1 regularization: lambda = " <<  RegulParam << endl;
    }
    
    if (NbrBase == 0)
    {
       NbrBase = 2;
       TabBase[(int) MB_CUR] = True;    // Curvelet                      
       if (UsePsf == False) TabBase[(int) MB_WT] = True; // UWT
       else  TabBase[(int) MB_WTMIRROR] = True; 
       TabSelect[0] = (UsePsf == False) ? MB_WT: MB_WTMIRROR;
       TabSelect[1] = MB_CUR;
    }
    
    
    MRBase MB;  // Multiple base decomposition Class

    io_read_ima_float(Name_Imag_In, MB.Data, &Header);
    Header.origin = Cmd;

    // ridgelet_check_size_ima(MB.Data.nl(), MB.Data.nc(), BlockSize);
    
    if ((PoissonNoise == False) && (Noise_Ima < FLOAT_EPSILON))
    {
       Noise_Ima = detect_noise_from_med (MB.Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    if (PoissonNoise == True) Noise_Ima = 1.;

    // for (i = 0; i < NBR_MRBASE; i++) MB.TabBase[i] = TabBase[i];
    MB.NbrBase = NbrBase;
    for (i = 0; i < NbrBase; i++) MB.TabSelect[i] = TabSelect[i];
    MB.Filtering = Filtering;
    MB.DataSigmaNoise = Noise_Ima;
    MB.N_Sigma = N_Sigma;
    MB.PoissonNoise = PoissonNoise;
    MB.Verbose = Verbose;
    MB.Nbr_Iter = Nbr_Iter;
    MB.PositivRecIma = PositivSol;
    MB.FirstSoftThreshold = FirstSoftThreshold;
    MB.LastSoftThreshold = LastSoftThreshold;
    MB.DetectCoefTol = DetectCoefTol;

    MB.RID.RidTrans = RidTrans;
    MB.RID.BlockOverlap = BlockOverlap;
    MB.RID_FirstDetectScale = FirstDetectScale;
    // MB.RID.KillLastScale = RidKillLastScale;
    if (BlockSize > 0) MB.RID_BlockSize = BlockSize;

    MB.AT_NbrScale2D = NbrScale2D;
    MB.AT_PositivRecIma = True;
    MB.AT_KillLastScale  = KillLastScale;

    MB.PMT_NbrScale2D = NbrScale2D;
    MB.PMT_PositivRecIma = True;
    MB.PMT_KillLastScale  = KillLastScale;
 
    if (BlockSize > 0) MB.MRID_BlockSize = BlockSize;
    MB.MRID_NbrRid = NbrScale2D;
    MB.MRID_FirstDetectScale = FirstDetectScale;
    MB.MRID_BlockOverlap = BlockOverlap;
    // MB.MRID_KillLastScale = RidKillLastScale;

    MB.WT_NbrScale2D = NbrScale2D;
    MB.WT_NbrUndecimatedScale = WT_NbrUndecimatedScale;
    MB.MB_NbrScale2D = NbrScale2D;

    MB.CUR_NbrScale2D = NbrScale2D;
    // MB.CUR_BlockSize = BlockSize;
    MB.CUR_BlockOverlap = BlockOverlap;
    MB.CUR_KillLastScale  = KillLastScale;
    
    if (CosBlockSize > 0) MB.COS_BlockSize = CosBlockSize;
    MB.COS_Overlapping = BlockOverlap;
    MB.CosWeightFirst = CosWeightFirst;
     
    MB.WP_NbrUndecimatedScale = WP_NbrUndecimatedScale;
    MB.WP_NbrScale2D = NbrScale2D;
    MB.WP_PositivRecIma = False;
    MB.WP_Filter = F_MALLAT_7_9;
    MB.WP_KillLastScale = False;
    
    MB.FCUR_NbrDir = FCurNdir;
    MB.FCUR_NbrScale2D = NbrScale2D;
    
    MB.TotalVariation = TV;
    MB.LambdaTV = RegulParam;
    MB.alloc(); 
    if (UsePsf == True)
    {     
        Ifloat ImaPSF;
        io_read_ima_float(PSFFile, ImaPSF);
        MB.init_deconv(ImaPSF, Noise_Ima,  Eps,  InitRnd);
    }
    MB.decomposition();
    // MB.reconstruction(MB.Data);
    io_write_ima_float(Name_Imag_Out, MB.Result, &Header);
    MB.Data -= MB.Result;
    if (Verbose == True) io_write_ima_float("xx_resi.fits", MB.Data);

    exit(0);
}
