/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/02/00
**
**    File:  dct.cc
**
*******************************************************************************
**
**    DESCRIPTION  dct program
**    -----------
**
******************************************************************************/

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "Filter.h"
#include "SB_Filter1D.h"
#include "SB_Filter.h"
#include "MR_NoiseModel.h"
#include "IM_Sigma.h"


char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */

extern int  OptInd;
extern char *OptArg;

extern int GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int NbrScaleLine=-1;
int NbrScaleCol=-1;
Bool PNbrScaleLine = False;
Bool PNbrScaleCol = False;

sb_type_norm Norm = NORM_L2;
type_sb_filter SB_Filter =  F_MALLAT_7_9; // F_HAAR; // F_BI2HAAR;  // F_HAAR   F_MALLAT_7_9;
type_border Bord = I_MIRROR;
Bool Reverse = False;
Bool MAD_Estim = False;
Bool PositivSol=True;
Bool Wiener=False;
int WienerBlockSize=7;
Bool FDR = False;
float N_Sigma=DEFAULT_N_SIGMA;
Bool OptNSigma= False;
float Noise_Ima=0.;
Bool WriteBand= False;
Bool KillLastScale = False;
type_border T_BORDER = I_MIRROR; // border type
Bool MultNoise = False;
int FirstScale=0;
int FirstScale1=0;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_data result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    // sb_usage(SB_Filter);
    // manline();

    gauss_usage();
    manline();
    nsigma_usage(N_Sigma);
    manline();
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Multiplicative noise.\n");
    manline();
    fprintf(OUTMAN, "         [-m]\n");
    fprintf(OUTMAN, "             Correlated noise (use MAD).\n");
    manline();
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "             Poisson noise (use MAD).\n");
    manline();
    fprintf(OUTMAN, "         [-W WindowSize]\n");
    fprintf(OUTMAN, "             Multiresolution Wiener Filtering. Recommended value is 7. \n");
    manline();
    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "             Coefficient detection using the FDR method \n");
    fprintf(OUTMAN, "             instead of the standard k-sigma approach. \n");
    manline();

    fprintf(OUTMAN, "         [-x NbrScale_onXaxis]\n");
    fprintf(OUTMAN, "             Number of scales in the line direction..\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();

    fprintf(OUTMAN, "         [-y NbrScale_onYaxis]\n");
    fprintf(OUTMAN, "             Number of scales in the column direction.\n");
    fprintf(OUTMAN, "             Default is automatically calculated.\n");
    manline();
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "             Write the wavelet bands to the disk (name bands are xx_band_i_j.fits).\n");
    manline();
    fprintf(OUTMAN, "         [-B BorderType]\n");
    fprintf(OUTMAN, "             B = 0 : symmetric border (default). \n");
    fprintf(OUTMAN, "             B = 1 : periodic border.\n");
    manline();
    fprintf(OUTMAN, "         [-F first_detection_scale] \n");
    fprintf(OUTMAN, "             First scale used for the detection.\n");
    fprintf(OUTMAN, "             Default is 1.\n");

    fprintf(OUTMAN, "             B = 0 : symmetric border (default). \n");
    fprintf(OUTMAN, "             B = 1 : periodic border.\n");
    manline();
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose.\n");
    manline();
    exit(-1);
}

/*********************************************************************/
/*************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}
/*************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c,Ret;
     /* get options */
    while ((c = GetOpt(argc,argv,"L:F:MB:Kws:mW:CT:x:y:vzZ")) != -1)
    {
	switch (c)
        {
	  case 'M':MultNoise = True;break;
          case 'B':
	    int tbrd;
	    if (sscanf(OptArg, "%d", &tbrd) != 1)
	    {
	      cerr << "Bad or missing parameter " << OptArg << endl;
	      exit (-1);
	    }
	    if (tbrd == 0) T_BORDER = I_MIRROR;
	    else T_BORDER = I_PERIOD;
	    break;
	  case 'K':KillLastScale=True;
                break;
	  case 'w':WriteBand=True;
                break;
          case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1)
                {
                    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
                    usage(argv);
                }
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
		OptNSigma=True;
                break;
          case 'F':
                 if (sscanf(OptArg,"%d", &FirstScale) != 1)
                {
                    fprintf(OUTMAN, "bad First Scale: %s\n", OptArg);
                    usage(argv);
                }
		FirstScale --;
                if (FirstScale < 0.) FirstScale = 0;
                break;
	  case 'L':
                 if (sscanf(OptArg,"%d", &FirstScale1) != 1)
                {
                    fprintf(OUTMAN, "bad First Scale: %s\n", OptArg);
                    usage(argv);
                }
		FirstScale1 --;
                if (FirstScale1 < 0.) FirstScale1 = 0;
                break;
	  case 'm': MAD_Estim = True; break;
          case 'W':
                /* -s <nsigma> */
                Ret = sscanf(OptArg,"%d",&c);
		if (Ret == 1) WienerBlockSize = c;
		else if (Ret > 1)
                {
                    fprintf(OUTMAN, "bad BlockSize: %s\n", OptArg);
                    usage(argv);
                }
 		Wiener=True;
                if (WienerBlockSize <= 0.) WienerBlockSize = 7;
                break;
          case 'C': FDR = True; break;
 	  case 'T':
 		SB_Filter = get_filter_bank(OptArg);
		break;
 	   case 'y':
                if (sscanf(OptArg,"%d",&NbrScaleLine) != 1)
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
		PNbrScaleLine = True;
		break;
           case 'x':
                if (sscanf(OptArg,"%d",&NbrScaleCol) != 1)
                {
                    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
		PNbrScaleCol = True;
 	        break;
          case 'v': Verbose = True;break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	}
        if ((OptNSigma == False) && (FDR == True)) N_Sigma = 2;
        if ((OptNSigma == False) && (Wiener == True)) N_Sigma = 1.;

        if ((FDR == True) && (Wiener == True))
	{
           fprintf(OUTMAN, "Error: Wiener option cannot be used with FDR ...\n");
           exit(-1);
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
}

/*********************************************************************/

static double prob_noise(float Val, float Sig)
{
   double P=0;
   if (ABS(Val) < FLOAT_EPSILON) P = 1.;
   else
   {
      if (Sig < FLOAT_EPSILON) P = 0;
      else
      {
         double Vn = ABS(Val) / (sqrt(2.)*Sig);
         if (Vn > 3.5) P = 0;
         else P = (float) erfc (Vn);
       }
    }
    return P;
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

    /* Get command line arguments, open input file(s) if necessary */
#ifndef MRPOL
    lm_check(LIC_MR4);
#else
    lm_check(LIC_POL);
#endif
    filtinit(argc, argv);

    if (Verbose == True)
    {
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;
	if (PNbrScaleLine == True) cout << "Number of scales on Y-axis = " << NbrScaleLine    << endl;
	if (PNbrScaleCol == True)  cout << "Number of scales on X-axis= " << NbrScaleCol    << endl;
        cout << "NSigma = " << N_Sigma  << endl;
	if (FDR == True) cout << "Detection method = FDR" << endl;
	else cout << "Detection method = K-Sigma" << endl;
	if (Wiener == True) cout << "Multiscale Wiener Filtering: WienerBlockSize = " << WienerBlockSize << endl;
	if (FirstScale > 0) cout << "First Detection Scale = " << FirstScale << endl;
    }

    Ifloat Data;
    Ifloat Result;
    int FirstScaleLine = FirstScale;
    int FirstScaleCol = FirstScale1;

    io_read_ima_float(Name_Imag_In, Data, &Header);

    float Min,Max;
    int Nx = Data.nc();
    int Ny = Data.nl();
    if (Verbose == True)
    {
       cout << "Nx = " << Nx <<  " Ny = " << Ny  <<  endl;
       cout << "Min = " << Data.min()  << "  Max = " << Data.max() <<  " Sigma = " << Data.sigma() << endl;
    }

    FilterAnaSynt FAS;
    // FAS.Verbose = Verbose;
    FAS.alloc(SB_Filter);
    SubBandFilter SBF(FAS,  Norm);
    SBF.Border = T_BORDER;
    LineCol LC;
    Bool UseMirror=False;
    LC.alloc(SBF, UseMirror);

    FilterAnaSynt *WT_SelectFilter = new FilterAnaSynt(SB_Filter);
    SubBandFilter *WT_SB1D = new SubBandFilter(*WT_SelectFilter, NORM_L2);
    // DirectionalLineCol LC1(*WT_SB1D);
    // LC1.NbrUndecimatedScaleCol = 0;
    // LC1.NbrUndecimatedScaleLine = 1;

    if (PNbrScaleLine == False)  NbrScaleLine = max_scale_number(Ny);
    if (PNbrScaleCol == False)  NbrScaleCol = max_scale_number(Nx);
    if (Verbose == True) cout << " NbrScaleY = " << NbrScaleLine << "  NbrScaleX = " << NbrScaleCol << endl;

    if (MultNoise == True) noise_log_transform (Data, Data);

    if ((Noise_Ima < FLOAT_EPSILON) && (MAD_Estim == False))
    {
       Noise_Ima = detect_noise_from_med (Data);
       if (Verbose == True)
          cout << "Sigma Noise = " << Noise_Ima << endl;
    }


    Ifloat **TabTrans;

    LC.undec_transform (Data, TabTrans, NbrScaleCol, NbrScaleLine);
    for (int si = 0; si < NbrScaleCol; si++)
    for (int sj = 0; sj < NbrScaleLine; sj++)
      if ((si < FirstScaleCol) || (sj < FirstScaleLine))
      {
           if (Verbose == True) cout << "INIT  Scale " << si+1 << " " << sj+1  <<  endl;
          (TabTrans[si][sj]).init();
      }

    for (int si = 0; si < NbrScaleCol; si++)
    for (int sj = 0; sj < NbrScaleLine; sj++)
    {
       char Name[80];
       float NSigmaBand = ((si ==0) || (sj ==0)) ? N_Sigma + 1: N_Sigma;
       if (WriteBand == True)
       {
           sprintf(Name, "xx_band_%d_%d.fits", si+1,sj+1);
           io_write_ima_float(Name, TabTrans[si][sj], &Header);
       }

       if (MAD_Estim == True) Noise_Ima = detect_noise_from_mad(TabTrans[si][sj]);
       if (FDR == True)
       {
          dblarray PVal(Nx,Ny);
          for (int i=0; i < Ny; i++)
          for (int j=0; j < Nx; j++) PVal(j,i) = prob_noise((TabTrans[si][sj])(i,j), Noise_Ima);
   	  float Alpha = (1. - erf(N_Sigma / sqrt(2.)));
	  double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);
          NSigmaBand = ABS(xerfc(0.5+(1-PDet)/2.));
	  if((NSigmaBand < 5)||(NSigmaBand > 0)) NSigmaBand = ABS(xerfc(0.5+(1-PDet)/2.));
          else NSigmaBand = 5;
       }
       else NSigmaBand =N_Sigma;
       float Level = NSigmaBand * Noise_Ima;

       if (Verbose == True) cout << " Scale " << si+1 << " " << sj+1 << " N_Sigma = " << NSigmaBand << " Noise = " << Noise_Ima << " T = " << Level <<  endl;
       int Step=1;
       float  VarianceSignal=0.,VarianceData=0., VarianceNoise =  Noise_Ima*Noise_Ima*N_Sigma*N_Sigma;
       int B2 = WienerBlockSize / 2;

       if ((si != NbrScaleCol-1)  || (sj != NbrScaleLine-1))
       {
         for (int i=0; i < Ny; i++)
         for (int j=0; j < Nx; j++)
         {
	   if (Wiener == False)
	   {
	      if (ABS( (TabTrans[si][sj])(i,j)) < Level) (TabTrans[si][sj])(i,j) = 0.;
	   }
	   else
           {
               for (int k=i-B2*Step; k <= i+B2*Step; k+=Step)
               for (int l=j-B2*Step; l <= j+B2*Step; l+=Step)
		        VarianceData += (TabTrans[si][sj])(k,l,I_MIRROR)*(TabTrans[si][sj])(k,l,I_MIRROR);
               VarianceData /= (float)(WienerBlockSize*WienerBlockSize);
	       VarianceSignal = MAX(0,VarianceData-VarianceNoise);
 	       (TabTrans[si][sj])(i,j) *= VarianceSignal  / VarianceData;
	   }
         }
       }
       else if (KillLastScale == True)  (TabTrans[si][sj]).init();
         // INFO_X(TabTrans[si][sj], "ScaleI");
    }

    if (Verbose == True) cout << "Reconstruction ... " << endl;
    LC.undec_recons( TabTrans, Data, NbrScaleCol, NbrScaleLine);
    if (Verbose == True) cout << "Write ... " << endl;
    if (MultNoise == True) noise_inv_log_transform(Data, Data);
    io_write_ima_float(Name_Imag_Out, Data, &Header);

//     Ifloat Res;
//     if (Reverse == False) LC1.transform(Data, Res, NbrScaleCol, NbrScaleLine);
//     else LC1.recons(Data, Res, NbrScaleCol, NbrScaleLine);
//
//     io_write_ima_float(Name_Imag_Out, Res, &Header);
    exit(0);
}
