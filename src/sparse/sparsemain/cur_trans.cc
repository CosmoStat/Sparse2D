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
**    File:  cur_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet transform program
**    -----------
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Curvelet.h"

// Fabrice Poupon 2017/03/14
// For performance monitoring during code optimization
// This was only used on Unix-based systems
// Uncomment the line below to activate it
#define FP_CPU_TIMING
#ifdef FP_CPU_TIMING
#include <sys/time.h>
#endif

char Name_Imag_In[256];  // input file image
char Name_Imag_Out[256]; //  output file name

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool Reverse=False; // Reverse transform
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
int NbrScale1D = -1;
int BlockSize= DEF_CUR_BLOCK_SIZE;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool StatInfo = False;
Bool BlockOverlap = False;
Bool DebugTrans = False;
Bool WriteBand = False;
Bool MadAngleNorm= False;

/***************************************/

/***************************************/
static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
    manline();

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
    fprintf(OUTMAN, "             default is %d. \n", BlockSize);
    manline();

    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse Curvelet transform.\n");
    manline();

    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Print statistical information about.\n");
    fprintf(OUTMAN, "             each band. Default is no. \n");
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Use overlapping block. Default is no.  \n");
    manline();
    write_scales_x_band_usage();
    manline();

//     fprintf(OUTMAN, "         [-f]\n");
//     fprintf(OUTMAN, "             Filter each scan of the Radon transform.\n");
//     manline();
//
//     fprintf(OUTMAN, "         [-w]\n");
//     fprintf(OUTMAN, "             Filter width.\n");
//     fprintf(OUTMAN, "             Default is 100.\n");
//     manline();
//
//     fprintf(OUTMAN, "         [-s]\n");
//     fprintf(OUTMAN, "             Sigma paameter for the filtering.\n");
//     fprintf(OUTMAN, "             Default is 10.\n");
//     manline();

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
    while ((c = GetOpt(argc,argv,"MxXt:Oib:n:N:rvzZ")) != -1)
    {
	switch (c)
        {
          case 'M': MadAngleNorm= True; break;
	  case 'x': WriteBand= True; break;
	  case 'X': DebugTrans = True; break;
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
           case 'i': StatInfo= (StatInfo ==True) ? False: True;break;
           case 'O': BlockOverlap= (BlockOverlap==True) ? False: True;break;
           case 'r': Reverse=True;break;
           case 'v': Verbose = True;break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1)
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
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

       if ((StatInfo == True) && (Reverse == True))
       {
           fprintf(OUTMAN, "Error: i and r option are not compatible...\n");
           exit(-1);
       }

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/***************/

int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    Header.origin = Cmd;

#ifdef FP_CPU_TIMING
    struct timeval t0, t1, t2, t3, t4, t5, t6;
    gettimeofday( &t0, NULL );
#endif

#ifdef _OPENMP
    // The maximum number of threads returned by omp_get_max_threads()
    // (which is the default number of threads used by OMP in parallel
    // regions) can sometimes be far below the number of CPUs.
    // It is then required to set it in relation to the real number of CPUs
    // (the -1 is used to live one thread to the main process which ensures
    // better and more constant performances). - Fabrice Poupon 2013/03/09
    omp_set_num_threads( omp_get_num_procs() - 1 );
#endif

    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);
#ifdef FP_CPU_TIMING
    gettimeofday( &t1, NULL );
#endif

    if (Verbose == True)
    {
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;
        if (Reverse == True) cout << "Inverse Curvelet transform  " <<  endl;
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;
    }

    if (NbrScale2D < 1)
    {
        cout << "Error: number of scales must be larger than 1 ... " << endl;
        exit(-1);
    }
    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
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
    Cur.StatInfo = StatInfo;

#ifdef FP_CPU_TIMING
    gettimeofday( &t2, NULL );
#endif

    if (Reverse == False)
    {
       Ifloat Data;
       fltarray Result;
       io_read_ima_float(Name_Imag_In, Data, &Header);
#ifdef FP_CPU_TIMING
    gettimeofday( &t3, NULL );
#endif
       Header.origin = Cmd;
//        if ((is_power_of_2(Data.nl()) == False) || (is_power_of_2(Data.nc()) == False))
//        {
//            cout << "Error: image size must be power of two ... " << endl;
//            exit(-1);
//        }

       if (BlockSize > Data.nc())
       {
          cout << "Error: Block size must lower than the image size ... " << endl;
          exit(-1);
       }
       Cur.VarNorm = True;
       Cur.alloc(Data.nl(), Data.nc(), BlockSize);
       Cur.AngleNormalization=MadAngleNorm;

//        Ifloat *TabB;
//        Cur.transform(Data, TabB);
//        Data.init();
//        Cur.recons(TabB, Data);
//        // INFO_X(Data,"REC");
//        io_write_ima_float(Name_Imag_Out,Data);
//        exit(-1);

#ifdef FP_CPU_TIMING
    gettimeofday( &t4, NULL );
#endif
       Cur.transform(Data,Result);
#ifdef FP_CPU_TIMING
    gettimeofday( &t5, NULL );
#endif
       Cur.write(Name_Imag_Out,Result);
#ifdef FP_CPU_TIMING
    gettimeofday( &t6, NULL );
    double total=double(t6.tv_sec+t6.tv_usec/1000000.0-t0.tv_sec-t0.tv_usec/1000000.0);
    double tInit=double(t1.tv_sec+t1.tv_usec/1000000.0-t0.tv_sec-t0.tv_usec/1000000.0);
    double tBand=double(t2.tv_sec+t2.tv_usec/1000000.0-t1.tv_sec-t1.tv_usec/1000000.0);
    double tRead=double(t3.tv_sec+t3.tv_usec/1000000.0-t2.tv_sec-t2.tv_usec/1000000.0);
    double tAlloc=double(t4.tv_sec+t4.tv_usec/1000000.0-t3.tv_sec-t3.tv_usec/1000000.0);
    double tCompute=double(t5.tv_sec+t5.tv_usec/1000000.0-t4.tv_sec-t4.tv_usec/1000000.0);
    double tWrite=double(t6.tv_sec+t6.tv_usec/1000000.0-t5.tv_sec-t5.tv_usec/1000000.0);
    cout << "Overall execution time : " << total << "s" << endl;
    cout << "   initialization  : " << tInit << "s" << endl;
    cout << "   band allocation : " << tBand << "s" << endl;
    cout << "   data loading    : " << tRead << "s" << endl;
    cout << "   allocation      : " << tAlloc << "s" << endl;
    cout << "   computation     : " << tCompute << "s" << endl;
    cout << "   writing results : " << tWrite << "s" << endl;
#endif
       if (WriteBand == True)
       {
          char ch[256];
          char Prefix[256];
	  fltarray Band;
 	  io_strcpy_prefix(Prefix,  Name_Imag_Out);
          for (int b=0; b < Cur.nbr_band(); b++)
          {
              int N, s2d, s1d;
 	      double Mean, Sigma, Skew, Curt;
	      float  Min, Max;
	      Ifloat IBand;
              Cur.get_scale_number(b, s2d, s1d);
 	      Cur.get_band(Result, b, Band);
	      N = Band.nx()*Band.ny();
              moment4(Band.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
	      IBand.alloc(Band.buffer(), Band.ny(), Band.nx());
              sprintf (ch, "%s_band_%d", Prefix, b+1);
              cout << "Write " << ch << " Size = " << IBand.nl() << " " << IBand.nc() << endl;
              printf("  Band %d (%d,%d): Nl = %d, Nc = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Min = %5.3f, Max = %5.3f\n",
	           b+1, s2d+1, s1d+1, Band.ny(), Band.nx(), Sigma,
		   (float) Skew, (float) Curt, Min, Max);
             io_write_ima_float(ch, IBand);
	  }
       }
//        fits_write_fltarr(Name_Imag_Out, Result);
       if (DebugTrans == True)
       {
          Ifloat Rec(Data.nl(), Data.nc(), "REC");
          Cur.recons(Result,Rec, False);
          io_write_ima_float("xx_rec", Rec);
          Data -= Rec;
          INFO_X(Data, "Resi");
          io_write_ima_float("xx_resi", Data);
       }
    }
    else
    {
       fltarray Data;
       Ifloat Result,RidIma;
       Cur.read(Name_Imag_In,Data);
       Cur.VarNorm = True;

//        {
//           char ch[256];
//           char Prefix[256];
// 	  fltarray Band;
//  	  io_strcpy_prefix(Prefix,  Name_Imag_Out);
//           for (int b=0; b < Cur.nbr_band(); b++)
//           {
//               int N, Nx, Ny, s2d, s1d;
//  	      double Mean, Sigma, Skew, Curt;
// 	      float  Min, Max;
// 	      Ifloat IBand;
//               Cur.get_scale_number(b, s2d, s1d);
// 	      Ridgelet *Rid = Cur.get_ridgelet(s2d);
// 	      float *PtrRidIma= Data.buffer()+ s2d*Data.ny()*Data.nx();
// 	      RidIma.alloc(PtrRidIma, Data.ny(), Data.nx());
// 	      Rid->get_imablocksigma(RidIma, s1d, IBand);
//               sprintf (ch, "%s_band_%d", Prefix, b+1);
//               cout << "Write " << ch << " Size = " << IBand.nl() << " " << IBand.nc() << endl;
//               io_write_ima_float(ch, IBand);
// 	     }
//        }
       Cur.recons(Data,Result, False);
       Header.width = Result.nc();
       Header.height = Result.nl();
       Header.npix = Header.width*Header.height;
       Header.bitpix = -32;
       type_format FormatData = io_which_format(Name_Imag_Out);
       if ((FormatData == F_UNKNOWN) && (Cur.FormatInputImag != F_UNKNOWN))
                                           io_set_format(Cur.FormatInputImag);
       io_write_ima_float(Name_Imag_Out, Result, &Header);
    }
    exit(0);
}
