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
#include "Ridgelet3D.h"
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
Bool BlockOverlap=True;
type_ridgelet3d_WTtrans RidTrans = DEF_RID3D_TRANS;
int FirstDetectScale=0;
Bool KillLastScale=False;
 
Bool WriteFilterCoef=False;
  
/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");


    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID3D_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRid3DTransform((type_ridgelet3d_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRid3DTransform(RidTrans));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated.\n");    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is image size. \n");    
    manline();

    fprintf(OUTMAN, "         [-F FirstDetectionScale]\n");
    fprintf(OUTMAN, "             First detection scale.\n");
    fprintf(OUTMAN, "             default is 1. \n");    
    manline();

    nsigma_usage(N_Sigma);
    manline();

    gauss_usage();
    manline();  
     
    poisson_noise_usage();
    manline(); 

    nsigma_usage(N_Sigma);
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             No block overlapping. Default is no. \n");    
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
    while ((c = GetOpt(argc,argv,"PkF:t:Ob:s:g:n:pvzZ")) != -1) 
    {
	switch (c)
        {   
          // case 'k': KillFluctu = True; break;
          case 'P': PositivSol = (PositivSol == True) ? False: True; break;
          case 'k': KillLastScale = True; break;
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
             case 't': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID3D_TRANS)) RidTrans = (type_ridgelet3d_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
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
                if (N_Sigma < 0.) N_Sigma = DEFAULT_N_SIGMA;
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
    int i,j,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    // lm_check(LIC_MR3);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRid3DTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        if (FirstDetectScale > 0) cout << "FirstDetectScale = " << FirstDetectScale+1 << endl;

        if (PoissonNoise == True)
           cout << "Poisson noise " << endl; 
         cout << "NSigma = " << N_Sigma  << endl;  
    }

    fltarray Data;
    io_3d_read_data(Name_Imag_In, Data);
    // Header.origin = Cmd;
 
    if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False))
    {
        cout << "Error: Block size must be a power of two ... " << endl;
        exit(-1);
    }
  
 
    if (BlockSize > Data.nx())
    {
        cout << "Warning: Block size must lower than the cube size ... " << endl;
        cout << "         Block size is set to cube size " << endl;
        BlockSize = 0;
    }
    
//     if ((PoissonNoise == False) && (Noise_Ima < FLOAT_EPSILON))
//     {
//        Noise_Ima = detect_noise_from_med (Data);
//        if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
//     } 
    
    // FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    FilterAnaSynt SelectFilter(F_HAAR);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    SB1D.Border = I_MIRROR;
    // SB1D.Border = I_PERIOD;

    Ridgelet3D RG(SB1D);
    RG.RidTrans = RidTrans;
    if (Nbr_Plan != DEF_RID3D_NBR_SCALE) 
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
    RG.NbrScale = Nbr_Plan;
    RG.Verbose = Verbose;
     // In case of Poisson noise, we apply the variance stabilization    
    RG.VarianceStab=PoissonNoise;
    // RG.MadThreshold = False;
 
    if (Verbose == True) cout << "Filtering ... " << endl;
    if (BlockSize < DEF_RID3D_MIN_BLOCK_SIZE) BlockSize = 0;

    //RG.filtering(Data,  Noise_Ima, N_Sigma, BlockSize);
 
    if ((PositivSol == True) && (KillLastScale == False)) 
    {
       for (i=0; i < Data.nx(); i++)
       for (j=0; j < Data.ny(); j++)
       for (k=0; k < Data.nz(); k++) if (Data(i,j,k) < 0) Data(i,j,k) = 0.;
    }
     
    io_3d_write_data(Name_Imag_Out, Data);
    exit(0);
}
