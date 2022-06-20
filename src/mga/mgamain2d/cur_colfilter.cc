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
**    File:  cur_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet filtering program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Curvelet.h"
#include "Usage.h"
#include "IM3D_IO.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
Bool PositivSol=True;
Bool Simu=False;
int BlockSize=DEF_CUR_BLOCK_SIZE;
int NbrScale1D = -1;
Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;

/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

//     fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
//     for (int i = 0; i < NBR_RID_TRANS; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
//     fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
//     manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    manline();

    fprintf(OUTMAN, "         [-N number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated. \n");    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is 16. \n");    
    manline();

    gauss_usage();
    manline();  
     
    nsigma_usage(N_Sigma);
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Do not apply block overlapping. \n"); 
    fprintf(OUTMAN, "             By default, block overlapping is used. \n");    
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
    while ((c = GetOpt(argc,argv,"t:Ob:Ss:g:n:N:vzZ")) != -1) 
    {
	switch (c) 
        { 
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
           case 'S': Simu=True;break;
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
	   case 'N':
                 if (sscanf(OptArg,"%d",&NbrScale1D) != 1) 
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
    int  i,j,k,b;
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
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
        cout << "NSigma = " << N_Sigma  << endl;  
    }

    fltarray Data;
    io_3d_read_data(Name_Imag_In, Data, &Header);
    int Nx = Data.nx();
    int Ny = Data.ny();
    int Nz = Data.nz();
    if (Nz == 0) Nz = 1;

    if (Verbose == True) cout << "rgb_to_yuv  ... " << endl;
    rgb_to_yuv(Data);

    Header.origin = Cmd;
    if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False))
    {
           cout << "Error: Block size must be a power of two ... " << endl;
           exit(-1);
    }
    if (BlockSize > Nx)
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }

    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    SB1D.Border = I_MIRROR;
    Curvelet Cur(SB1D);
    Cur.RidTrans = RidTrans;
    Cur.NbrScale2D = NbrScale2D;
    if (NbrScale1D >= 0)
    {
        Cur.NbrScale1D = NbrScale1D;
        Cur.GetAutoNbScale = False;
        if (NbrScale1D < 2) Cur.BlockOverlap = False;
        else Cur.BlockOverlap = BlockOverlap;
    }
    else Cur.BlockOverlap = BlockOverlap;
    Cur.Verbose = Verbose;
    Cur.OddBlockSize=True;
    Cur.alloc(Ny, Nx, BlockSize);
     
    if (Verbose == True) cout << "Filtering ... " << endl;
    Ifloat Filter(Ny, Nx, "filter");

    for (b = 0; b < Nz; b++)
    {
       if (Verbose == True) cout << "Frame " << b+1 << endl;
       Ifloat Frame(Ny,Nx,"Frame");
       for (i=0; i < Ny; i++)
       for (j=0; j < Nx; j++) Frame(i,j) = Data(j,i,b);
       char NameBand[256];
       sprintf(NameBand, "band_%d.fits", b+1);

       if (Noise_Ima < FLOAT_EPSILON) 
       {
          Noise_Ima  = detect_noise_from_med (Frame);
          if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
       }
       Cur.filtering(Frame, Filter, Noise_Ima, N_Sigma);

       for (i=0; i < Ny; i++)
       for (j=0; j < Nx; j++) Data(j,i,b) = Filter(i,j);
    }
    if (Verbose == True) cout << "Write ... " << endl;
    yuv_to_rgb(Data);
    col_rescale(Data);
    //for (b = 0; b < Nz; b++)
    //for (i=0; i < Ny; i++)
    //for (j=0; j < Nx; j++) Data(j,i,b) =   saturation(Data(j,i,b), 255);
    io_3d_write_data(Name_Imag_Out, Data, &Header);
    exit(0);
}

