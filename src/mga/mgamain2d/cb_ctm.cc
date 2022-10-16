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
**    File:  cb_mca.cc
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
#include "MBase.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale1D = -1;
int NbrScale2D = 4;

float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
Bool PoissonNoise=False;

Bool PositivSol=True;
int BlockSize=0;
int Nbr_Iter = DEF_MB_NBR_ITER;
 Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
int FirstDetectScale=0;
Bool KillLastScale = False;
Bool ThresholdRecIma = True;

Bool Filtering = False;

Bool TabBase[NBR_MBASE];
type_mbase TabSelect[NBR_MBASE];
int NbrBase = 0;

int WT_NbrUndecimatedScale = 1;
Bool WriteAllRec=True;
Bool UseHuberNorm=False;
float FirstSoftThreshold = DEF_MB_FIRST_SOFT_THRESHOLD;
float LastSoftThreshold = DEF_MB_LAST_SOFT_THRESHOLD;
Bool RidKillLastScale = False;
Bool  UseNormL1 = False;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t TransformSelection]\n");
    for (int i = 0; i < NBR_MBASE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringMBase ((type_mbase) i));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the WT, the a trous, the PMT and the curvelet transform.\n");
    fprintf(OUTMAN, "             default is 4.\n");    
    fprintf(OUTMAN, "             Number of ridgelet in the multi-ridgelet transform.\n");
    manline();
 
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is image size. \n");    
    fprintf(OUTMAN, "             Starting Block Size in the multi-ridgelet transform.\n"); 
    fprintf(OUTMAN, "             default is 8. \n");    
    manline();

    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iteration. Default is%d.\n", Nbr_Iter);    
    manline();
 
    fprintf(OUTMAN, "         [-F FirstDetectionScale]\n");
    fprintf(OUTMAN, "             First detection scale in the (multi-) ridgelet transform.\n");
    fprintf(OUTMAN, "             default is 1. \n");    
    manline();

    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "            Kill last scale in ridgelet, and multiscale ridgelet transform.\n");
    fprintf(OUTMAN, "            default no.\n");    
    manline();

    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "            Kill last scale in the atrous algorithm, the PMT, and the curvelet.\n");
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

    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Apply the noise modeling (filtering).\n");
    fprintf(OUTMAN, "             default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-u]\n");
    fprintf(OUTMAN, "             Number of undecimated scales in the WT.\n");
    fprintf(OUTMAN, "             default is %d. \n", WT_NbrUndecimatedScale);    
    manline();

    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Write to the disk the reconstructed image.\n");
    fprintf(OUTMAN, "             from each decomposition. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-N]\n");
    fprintf(OUTMAN, "             Minimize the L1 norm.\n");
    fprintf(OUTMAN, "             Default is no. \n");    
    manline();
    
    nsigma_usage(N_Sigma);
    manline();

    gauss_usage();
    manline();  
     
    poisson_noise_usage();
    manline(); 

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Supress the block overlapping. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Supress the positivity constraint. Default is no. \n");    
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
    while ((c = GetOpt(argc,argv,"Nl:L:khwPu:t:fKF:Oi:b:s:g:n:pvzZ")) != -1) 
    {
	switch (c)
        {  
          case 'N': UseNormL1 = (UseNormL1  == True) ? False: True; break;
          case 'k': RidKillLastScale = True; break;
          case 'L':               
                if (sscanf(OptArg,"%f",&FirstSoftThreshold) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad first soft threshold: %s\n", OptArg);
                    exit(-1);
                }
                if (FirstSoftThreshold <= 0)
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
          case 'w': WriteAllRec = True; break;
          case 'P': PositivSol = False; break;
          case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_MBASE))
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
           case 'h': UseHuberNorm = (UseHuberNorm == True) ? False: True; break;
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
    for (i = 0; i < NBR_MBASE; i++) TabBase[i] = False;

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
        if (UseNormL1 == True) cout << "Use l1 optimazation " << endl;
    }

    MBase MB;  // Multiple base decomposition Class

    io_read_ima_float(Name_Imag_In, MB.Data, &Header);
    Header.origin = Cmd;

    // ridgelet_check_size_ima(MB.Data.nl(), MB.Data.nc(), BlockSize);
    
    if ((PoissonNoise == False) && (Noise_Ima < FLOAT_EPSILON))
    {
       Noise_Ima = detect_noise_from_med (MB.Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    if (PoissonNoise == True) Noise_Ima = 1.;

    // for (i = 0; i < NBR_MBASE; i++) MB.TabBase[i] = TabBase[i];
    MB.NbrBase = NbrBase;
    for (i = 0; i < NbrBase; i++) MB.TabSelect[i] = TabSelect[i];
    MB.Filtering = Filtering;
    MB.DataSigmaNoise = Noise_Ima;
    MB.N_Sigma = N_Sigma;
    MB.PoissonNoise = PoissonNoise;
    MB.Verbose = Verbose;
    MB.Nbr_Iter = Nbr_Iter;
    MB.PositivRecIma = PositivSol;
    MB.WT_PositivRecIma = PositivSol;
    MB.CUR_PositivRecIma = PositivSol;
    MB.MRID_PositivRecIma = PositivSol;
    MB.RID_PositivRecIma = PositivSol;
    MB.PMT_PositivRecIma =  PositivSol;
    MB.AT_PositivRecIma =  PositivSol;
    MB.COS_PositivRecIma = False;
    MB.FFT_PositivRecIma = False;
    MB.Bord = I_CONT;
    
    MB.UseHuberNorm = UseHuberNorm;
    MB.FirstSoftThreshold = FirstSoftThreshold;
    MB.LastSoftThreshold = LastSoftThreshold;
    MB.UseNormL1 = UseNormL1;

    MB.RID.RidTrans = RidTrans;
    MB.RID.BlockOverlap = BlockOverlap;
    MB.RID_FirstDetectScale = FirstDetectScale;
    MB.RID.KillLastScale = RidKillLastScale;
    if (BlockSize > 0) MB.RID_BlockSize = BlockSize;

    MB.AT_NbrScale2D = NbrScale2D;
    // MB.AT_PositivRecIma = True;
    MB.AT_KillLastScale  = KillLastScale;
    MB.AT_FirstDetectScale = FirstDetectScale;
    MB.AT_AdjointRec = False;
     
    
    MB.PMT_NbrScale2D = NbrScale2D;
    // MB.PMT_PositivRecIma = True;
    MB.PMT_KillLastScale  = KillLastScale;

 
    if (BlockSize > 0) MB.MRID_BlockSize = BlockSize;
    MB.MRID_NbrRid = NbrScale2D;
    MB.MRID_FirstDetectScale = FirstDetectScale;
    MB.MRID_BlockOverlap = BlockOverlap;
    MB.MRID_KillLastScale = RidKillLastScale;

    MB.WT_NbrScale2D = NbrScale2D;
    MB.WT_NbrUndecimatedScale = WT_NbrUndecimatedScale;

    MB.CUR_NbrScale2D = NbrScale2D;
    // MB.CUR_BlockSize = BlockSize;
    MB.CUR_BlockOverlap = BlockOverlap;
    MB.CUR_KillLastScale  = True; // KillLastScale;

     MB.alloc();
     MB.decomposition();
     MB.reconstruction(MB.Data);
     if (WriteAllRec == True) MB.write_allima();

    io_write_ima_float(Name_Imag_Out, MB.Data, &Header);
    exit(0);
}
