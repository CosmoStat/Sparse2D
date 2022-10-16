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
**    File:  ridgelet.cc
**
*******************************************************************************
**
**    DESCRIPTION  ridgelet filtering program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Ridgelet.h"
#include "Usage.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int Nbr_Plan = -1;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
Bool PoissonNoise=False;

Bool PositivSol=True;
Bool Simu=False;
int BlockSize=0;
int Nbr_Iter = 1;
float RegulParam = 0.2;
float CvgParam = 1.;
Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool OnlyHighFreq = False;
int FirstDetectScale=0;
Bool KillLastScale=False;
#ifdef DO_SIMU
#include"Simu.cc"
#endif

Bool WriteFilterCoef=False;
// Bool KillFluctu=False;
Bool ColTrans = False;

/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");


    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID_DECIMATED; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated.\n");    
    manline();

    fprintf(OUTMAN, "         [-h]\n");
    fprintf(OUTMAN, "             Apply the ridgelet transform only on the high frequencies.\n");
    fprintf(OUTMAN, "             default is no.\n");    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is image size. \n");    
    manline();

    fprintf(OUTMAN, "         [-F FirstDetectionScale]\n");
    fprintf(OUTMAN, "             First detection scale.\n");
    fprintf(OUTMAN, "             default is 1. \n");    
    manline();

//     fprintf(OUTMAN, "         [-i NbrIter]\n");
//     fprintf(OUTMAN, "             Number of iteration for the constraint reconstruction.\n");
//     fprintf(OUTMAN, "             default is 1 (no iteration). \n");    
//     manline();
//     
//     fprintf(OUTMAN, "         [-G RegulParam]\n");
//     fprintf(OUTMAN, "             Regularization parameter for the constraint reconstruction.\n");
//     fprintf(OUTMAN, "             default is 0.2 \n");    
//     manline();
//     
//     fprintf(OUTMAN, "         [-C ConvergParam]\n");
//     fprintf(OUTMAN, "            Convergence parameter.\n");
//     fprintf(OUTMAN, "             default is %f.\n", CvgParam);    
//     manline();

    nsigma_usage(N_Sigma);
    manline();

    gauss_usage();
    manline();  
     
    poisson_noise_usage();
    manline(); 

    //nsigma_usage(N_Sigma);
    //manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Do not apply block overlapping. \n"); 
    fprintf(OUTMAN, "             By default, block overlapping is used. \n");    
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
    while ((c = GetOpt(argc,argv,"CPkwF:ht:OG:i:b:Ss:g:n:pvzZ")) != -1) 
    {
	switch (c)
        {   
          // case 'k': KillFluctu = True; break;
          case 'C':  ColTrans =(ColTrans  == True) ? False: True; break; 
          case 'P': PositivSol = (PositivSol == True) ? False: True; break;
          case 'k': KillLastScale = True; break;
          case 'w': WriteFilterCoef = (WriteFilterCoef   == True) ? False: True; break; 
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
            case 'h': OnlyHighFreq  = (OnlyHighFreq  == True) ? False: True; break; 
            case 't': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID_TRANS)) RidTrans = (type_ridgelet_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
//            case 'C':
//                 if (sscanf(OptArg,"%f",&CvgParam) != 1) 
//                 {
//                     fprintf(OUTMAN, "bad Convergence parameter: %s\n", OptArg);
//                     exit(-1);
//                 }
//                 break;
           case 'G':
                if (sscanf(OptArg,"%f",&RegulParam) != 1) 
                {
                    fprintf(OUTMAN, "bad regularization parameter: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'i':
                if (sscanf(OptArg,"%d",&Nbr_Iter) != 1) 
                {
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'S': Simu=True;break;
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
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
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
      
       if ((OnlyHighFreq == true) && (BlockSize <= 0))
       {
           cout << "Error: option -h is valid only when block transform is activated ... " << endl;
           exit(-1);
       }
       if ((OnlyHighFreq == true) && (PoissonNoise == True))
       {
           cout << "Warning: option -h modify the poisson noise distribution ... " << endl;
           cout << "         Results may be not correct ... " << endl;
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
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        if (FirstDetectScale > 0) cout << "FirstDetectScale = " << FirstDetectScale+1 << endl;

        if (Nbr_Iter > 1)
        {
           cout << "Iterative reconstruction: nbr_iter = " << Nbr_Iter << endl;
           cout << "Regularization parameter = " << RegulParam << endl;
        }
        if (PoissonNoise == True)
           cout << "Poisson noise " << endl; 
         cout << "NSigma = " << N_Sigma  << endl;  
    }

    Ifloat Data;
    if (Simu == False) 
    {
    io_read_ima_float(Name_Imag_In, Data, &Header);
    // for (int q=0; q < Data.n_elem(); q++) Data(q) *= 100.;
    Header.origin = Cmd;
//     if ((is_power_of_2(Data.nl()) == False) || (is_power_of_2(Data.nc()) == False))
//     {
//            cout << "Error: image size must be power of two ... " << endl;
//            exit(-1);
//     }
    if (RidTrans == RID_FINITE) 
    {    
        CPRIME_NUMBER CPN;

        if ((BlockSize > 0) && ( CPN.is_prime_number((unsigned long) BlockSize) == False))
        {
           cout << "Error: Block size must be a prime number ... " << endl;
           cout << "       Previous prime number = " << CPN.previous_prime_number((unsigned long) BlockSize) << endl;
           cout << "       Next prime number = " << CPN.next_prime_number((unsigned long) BlockSize) << endl;
           exit(-1);
        }
    }
    if (BlockSize > Data.nc())
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }
    
    if ((PoissonNoise == False) && (Noise_Ima < FLOAT_EPSILON))
    {
       Noise_Ima = detect_noise_from_med (Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    }
#ifdef DO_SIMU
    if (Simu == True) 
    {
        Data.alloc(256,256,"simu");
        make_ima(Data);
    }
#endif
    
    // FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    FilterAnaSynt SelectFilter(F_HAAR);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    SB1D.Border = I_MIRROR;
    // SB1D.Border = I_PERIOD;

    Ridgelet RG(SB1D);
    RG.RidTrans = RidTrans;
    if (Nbr_Plan != DEF_RID_NBR_SCALE) 
    {
       RG.NbrScale = Nbr_Plan;
       RG.GetAutoNbScale = False;
       if (Nbr_Plan < 2)
       {
          RG.BlockOverlap = False;
       }
       else RG.BlockOverlap = BlockOverlap;
    }         
    else RG.BlockOverlap = BlockOverlap;
    RG.FirstDetectScale = FirstDetectScale;
    RG.KillLastScale = KillLastScale;
    RG.OnlyHighFreq = OnlyHighFreq;
    // if (OnlyHighFreq == True) RG.KillFluctu = KillFluctu;
    RG.NbrScale = Nbr_Plan;
    RG.Verbose = Verbose;
    RG.NbrFilterIter = Nbr_Iter;
    RG.CvgParam = CvgParam;
    RG.RegulParam = RegulParam;
    // In case of Poisson noise, we apply the variance stabilization    
    RG.VarianceStab=PoissonNoise;
    // RG.MadThreshold = False;
    RG.ColTrans = ColTrans;

    if (Verbose == True) cout << "Filtering ... " << endl;
    if (BlockSize < DEF_RID_MIN_BLOCK_SIZE) BlockSize = 0;

    if (WriteFilterCoef == False)
    {
       if (rid_class(RidTrans) != RID_CL_PAVE) 
          RG.filtering(Data,  Noise_Ima, N_Sigma, BlockSize);
        else RG.undec_wt_filtering(Data,  Noise_Ima, N_Sigma, BlockSize);
    }
    else
    {
       Ifloat Result;
       Ifloat Filter;
       int Nc = Data.nc();
       int Nl = Data.nl(); 
       if ((BlockSize > 0) && (OnlyHighFreq == True))
       {
          Ifloat ImaHigh(Nl,Nc, "ImaHigh");
          Filter.alloc(Nl,Nc, "Filter");
          RG.get_low_and_freq(Data, Filter, ImaHigh, BlockSize);
          Data = ImaHigh;
       }
       RG.transform(Data,Result,BlockSize);
       RG.thresholding(Result, Noise_Ima, N_Sigma);
       io_write_ima_float("xx_thres", Result);
       RG.recons(Result,Data,BlockSize);
       if ((BlockSize > 0) && (OnlyHighFreq == True)) Data += Filter;
    }

    if ((PositivSol == True) && (KillLastScale == False)) threshold(Data);
    io_write_ima_float(Name_Imag_Out, Data, &Header);
    exit(0);
}
