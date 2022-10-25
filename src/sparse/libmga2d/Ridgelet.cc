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
**    Date:  5/03/2000
**
**    File:  Ridgelet.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Ridgelet transform and reconstruction
**    -----------
**
******************************************************************************/

#include "Ridgelet.h"
#include "PrimeNumber.h"
// #include "LineCol.cc"
#include "FFTN_1D.h"
#include "MR1D_Obj.h"
#include "IM_Regul.h"
#include "MR_SoftRegul.h"
#include "Pyr1D.h"

#define DEGUG_RID 0

/*************************************************************/
// Ridgelet transform
/*************************************************************/

// Return true if the Fast Slant Stack Transform is used in a give ridgelet transform
const Bool UseFSS (type_ridgelet_WTtrans type)
{
   Bool Ret=False;
    switch (type)
    {
        case  RID_FSS:
        case  RID_FSS_PYR:
        case  RID_FSS_ORTHO:
        case  RID_UWT:
            Ret = True;
            break;
        default:
            Ret=False;
            break;
    }
    return Ret;
}

const char * StringRidTransform (type_ridgelet_WTtrans  type)
{
    switch (type)
    {
        case  RID_OWT:
	      return ("RectoPolar Ridgelet Transform using a standard bi-orthogonal WT");
        case  RID_PYR_FFT:
              return ("RectoPolar Ridgelet Transform using a FFT based Pyramidal WT");
	case RID_PYR_DIRECT:
              return ("RectoPolar Ridgelet Transform using a Pyramidal WT in direct space");
        case  RID_PAVE_FFT:
	      return ("Undecimated FFT based WT");
        case  RID_FINITE:
	      return ("Finite ridgelet transform");
        case  RID_FSS:
	      return ("Slant Stack Radon transform + FFT based pyramidal WT.");
	case  RID_FSS_PYR:
	      return ("Slant Stack Radon transformand + pyramidal WT in direct space");
	case  RID_FSS_ORTHO:
	      return ("Slant Stack Radon transformand + bi-orthogonal WT");
    case  RID_UWT:
            return ("Slant Stack Radon transformand + Undecimated Starlet WT");
	default:
	      return ("Undefined transform");
     }
}

type_ridgelet_WTclass rid_class(type_ridgelet_WTtrans Type)
{
    switch (Type)
    {
        case  RID_FINITE:
	      return  RID_CL_FINITE;
        case  RID_FSS_ORTHO:
	case  RID_OWT:
	      return  RID_CL_ORTHO;
        case  RID_PYR_FFT:
	case  RID_PYR_DIRECT:
	case  RID_FSS:
	case  RID_FSS_PYR:
              return  RID_CL_PYR;
        case  RID_UWT:
        case  RID_PAVE_FFT:
	      return RID_CL_PAVE;
        default:
	      return RID_CL_UNKNOW;
     }
}

/*************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

/****************************************************************************/

void Ridgelet::reset(SubBand1D *SB1D)
{
    extern type_format Format_Imag;
    FormatInputImag = Format_Imag;
    Ptr_SB1D = SB1D;
    VarianceStab=False;
    MadThreshold=False;
    WeightBefTrans=False;
    NbrFilterIter=1;
    RegulParam=0.2;
    CvgParam=1.;
    StatInfo=False;
    NbrScale=DEF_RID_NBR_SCALE;
    Verbose=False;
    BlockImaOverLap=0;
    BlockOverlap=True;
    GetAutoNbScale=True;
    RidTrans=DEF_RID_TRANS;
    RidClass = rid_class(RidTrans);
    WTinFFT=True;
    OnlyHighFreq = False;
    FirstDetectScale=0;
    KillLastScale=False;
    KillFluctu=False;
    ColTrans = False;
    WPTrans = False;
    SetTable = True;
    NoiseInvariantPerRotation = True;
    TabImaSigmaStat = NULL;
    VarNorm = False;
    AllocatedClass = False;
}

/*************************************************************************/

void Ridgelet::set_block_param_from_transform(int Nlt, int Nct, int BlockSize)
// intialize the variable for block reconstruction from
// the input ridgelet transform size
{
   int Nl,Nc;


   if (BlockSize == 0) BlockOverlap = False;

   // We need to calculate the size Nl,Nc of the original image
   // for setting the different parameters

   // The radon transform multiply be 2 the number of lines
   Nl = (rid_class(RidTrans) != RID_CL_FINITE) ? Nlt/2: Nlt-1;
   Nc = Nct;

   // For pyramidal transform, the redundancy of the WT is 2
   // (only for the column!)
   if ((GetAutoNbScale == True) || (NbrScale > 1))
   {
        if (rid_class(RidTrans) == RID_CL_PYR)
        {
           Nc /= 2;
           if ((RidTrans == RID_FSS) || (RidTrans == RID_FSS_PYR) || (RidTrans == RID_FSS_ORTHO)) Nc /= 2;
        }
        else if (RidTrans == RID_UWT)
        {
           if (GetAutoNbScale == True) NbrScale = max_scale_number(BlockSize);
           Nc /= NbrScale;  // UWT redundancy
           Nc /= 2; // FSS Radon redundancy
        }
   }

   // When overlapping, there is another redundancy of 2 in both direction
   if (BlockOverlap == True)
   {
      Nl /= 2;
      Nc /= 2;
   }

   // initialize the block parameters
   alloc(Nl, Nc, BlockSize);
}

/*************************************************************************/

void Ridgelet::alloc(int NlData, int NcData, int BlockSize)
{

 if ((AllocatedClass == False) || (NlIma != NlData) || (NcIma != NcData) || (BlockSize!= BS))
 {
    int RatioOverlap;
    CPRIME_NUMBER CPN;

    AllocatedClass = True;

    // Input image size
    NlIma = NlData;
    NcIma = NcData;
    RidClass = rid_class(RidTrans);
    if (RidClass != RID_CL_FINITE)
    {
       NlImaPow2 = next_power_of_2(NlIma);
       NcImaPow2 = next_power_of_2(NcIma);
    }
    else
    {
       NlImaPow2 = NlIma;
       NcImaPow2 = NcIma;
    }

   if (BlockSize <= 0)
   {
      Ncb = Nlb = 1;
      BlockOverlap = False;
      if (RidClass == RID_CL_FINITE)
      {
          NlImaPow2 = (int) CPN.next_prime_number((int) NlIma);
          NcImaPow2 = (int) CPN.next_prime_number((int) NcIma);
	      BS = BlockImaSize = MAX(NlImaPow2,NcImaPow2);
      }
      else
      {
         BlockImaSize = MAX(NlIma,NcIma);
	     BlockImaSize = next_power_of_2(BlockImaSize-1);
	      if ((RidTrans != RID_FSS)  && (RidTrans != RID_FSS_PYR) && (RidTrans != RID_FSS_ORTHO)  && (RidTrans != RID_UWT))
	      BlockImaSize ++;
         BS = BlockImaSize;
      }
      RatioOverlap = 1;
      BlockImaOverLap = 0;
   }
   else
   {
      // Set the the block sizes
      if (BlockOverlap == True)
      {
         BlockImaSize = (BlockSize+1) / 2;
         BS = BlockSize;
         BlockImaOverLap =  BlockSize/4;
         RatioOverlap = 2;
         Nlb =   NlIma / BlockImaSize;
         Ncb =   NcIma / BlockImaSize;
         if (Nlb*BlockImaSize < NlIma) Nlb++;
         if (Ncb*BlockImaSize < NcIma) Ncb++;
         if (RidClass == RID_CL_FINITE)
         {
             NlImaPow2 = Nlb*BlockImaSize;
             NcImaPow2 = Ncb*BlockImaSize;
         }
      }
      else
      {
         BS = BlockImaSize = BlockSize;
         BlockImaOverLap = 0;
         RatioOverlap = 1;
	     Nlb =  NlIma / BlockImaSize;
         Ncb =  NcIma / BlockImaSize;
	     if (Nlb*BlockImaSize < NlIma) Nlb++;
	     if (Ncb*BlockImaSize < NcIma) Ncb++;

         if (RidClass == RID_CL_FINITE)
         {
            NlImaPow2 = Nlb*BlockImaSize;
            NcImaPow2 = Ncb*BlockImaSize;
         }
      }
    }

   // Determine the number of scales to use
   if (GetAutoNbScale == True) NbrScale = max_scale_number(BS);

   if (Verbose == True) cout << "Number of scales = " << NbrScale << endl;
   if (NbrScale > 1) RidPixNbr = set_tab_pos(NbrScale, BS); // set the tables

    // Ridgelet image size
    // The FFT based Radon transform multiply by 2 the number of lines
    // If a pyramidal WT is used instead of an orthogonal one,
    // the number of columns are multiplied by 2.
    // Overlapping increase the ridgelet image size by the overlapping ratio
    if (RidClass != RID_CL_FINITE)
    {
       NlRid = 2*RatioOverlap*NlIma;
       NlBlockTransSize = 2*BS;
    }
    else
    {
       NlBlockTransSize = BS+1;
       NlRid = RatioOverlap*(NlImaPow2+Nlb);
    }
    NcBlockTransSize = BS;

    if (RidClass == RID_CL_PYR)
    {
        if ((GetAutoNbScale == True) || (NbrScale > 1))
        {
            NcBlockTransSize *=2;
        }
	    if ((RidTrans == RID_FSS) || (RidTrans == RID_FSS_PYR)) NcBlockTransSize *=2;
    }
    else if (RidTrans == RID_FSS_ORTHO) NcBlockTransSize *=2;
    else if (RidClass == RID_CL_PAVE)
    {
        NcBlockTransSize = 2 * BS * NbrScale;
    }

    NcRid = NcBlockTransSize*Ncb;
    NlRid = NlBlockTransSize*Nlb;
    set_tab_angle();
    TabImaSigmaStat = new Ifloat [NbrScale];

    if (RidTrans == RID_FINITE) RADON.alloc(BS, BS, RADON_FINITE);
    else
    {
       if ((RidTrans == RID_FSS) || (RidTrans == RID_FSS_PYR) || (RidTrans == RID_FSS_ORTHO) || (RidTrans == RID_UWT))
	    {
 	       if (Verbose == True) RADON.Verbose = True;
	       RADON.alloc(BS, BS, RADON_FSS);
           radon_transform_one_cst_block();
	       RADON.Verbose = False;
	    }
	    else RADON.alloc(BS, BS);
    }

    if (RidTrans == RID_UWT)
    {
       type_undec_filter U_Filter = U_B2SPLINE;
       Ptr_Undec_SB1D = new UndecSubBandFilter(U_Filter);
    }
    else Ptr_Undec_SB1D = NULL;

    if (Verbose == True)
    {
       cout << "RID INIT: " << endl;
       cout << "      Image size: Nl  = " <<  NlIma << " Nc = " <<  NcIma << endl;
       cout << "      Block size = " << BS << endl;
       cout << "      Ridgelet image size: Nl  = " <<  NlRid << " Nc = " <<  NcRid << endl;
       cout << "      Number of block: NbrBi = " << Nlb << " NbrBj = " << Ncb << endl;
       cout << "      Ridgelet Block size: Nlb = " << NlBlockTransSize << " Ncb = " << NcBlockTransSize << endl;
       cout << "      BlockImaSize = " << BlockImaSize << endl;
       cout << "      Overlapping size = " << BlockImaOverLap << endl;
       cout << "      NbrScale1D = " << NbrScale << endl;
    }
  }
}

/*************************************************************************/

void Ridgelet::block_transform(Ifloat &Ima, Ifloat &Transf)
// Ridgelet transform of image block per block
{
    int i,j;
    //int Nc = Ima.nc();
    //int Nl = Ima.nl();
    //alloc(Nl, Nc, BlockSize);
    Ifloat ImaBlock(BS, BS,"block");
    Ifloat ImaTransBlock(NlBlockTransSize,NcBlockTransSize,"blocktrans");

    if ((Transf.nl() != NlRid) || (Transf.nc() != NcRid))
                                          Transf.resize(NlRid, NcRid);
    // Loops on the block
    Bool UserVerbose = Verbose;
    if (UserVerbose == True)
    {
       cout << "RID Transform: Block Size = " <<  BS << " Nlb = " << Nlb << " Ncb = " << Ncb << endl;
       cout << "    Ridgelet image size  = " <<  NlRid << " " <<  NcRid  << endl;
       if (WPTrans == True)
         cout << "    Wavelet packets on the first scale " << endl;
    }
    // INFO_X(Ima, "Trans BAND");

    for (i = 0; i < Nlb; i++)
    for (j = 0; j < Ncb; j++)
    {
      // cout << "Block " << i << " " << j << endl;
        get_block_ima(i,j, Ima, ImaBlock);
      // cout << " Transform_one_block " << endl;
        transform_one_block(ImaBlock, ImaTransBlock);
	// io_write_ima_float("xx_b",  ImaTransBlock);

      // cout << " Put_block_trans  " << endl;
        put_block_trans(i,j, Transf, ImaTransBlock);
        if ((i == 0) && (j == 0)) Verbose = False;
    }
    Verbose = UserVerbose;
}

/***************************************/

void Ridgelet::block_recons(Ifloat &Transf, Ifloat &Ima)
// Image reconstruction from its ridgelet transform
// Block per block
{
    int i,j;
    Ifloat ImaBlock(BS, BS,"block");
    Ifloat ImaTransBlock(NlBlockTransSize,NcBlockTransSize,"blocktrans");

    if ((Ima.nl() != NlIma) || (Ima.nc() !=  NcIma)) Ima.resize(NlIma,NcIma);
    Ima.init();
    Bool UserVerbose = Verbose;
 // cout << "Block BS = " <<  BS << endl;
 // cout << "Block ImaTransBlock  = " <<  NlBlockTransSize << " " <<  NcBlockTransSize  << endl;
 // cout << "Ima  = " <<  NlIma << " " <<  NcIma  << endl;
    if (UserVerbose == True)
    {
       cout << "RID Transform: Block Size = " <<  BS << " Nlb = " << Nlb << " Ncb = " << Ncb << endl;
       cout << "    Ridgelet image size  = " <<  NlRid << " " <<  NcRid  << endl;
       cout << "    Reconstruted image  = " <<  NlIma << " " <<  NcIma  << endl;
    }

    for (i = 0; i < Nlb; i++)
    for (j = 0; j < Ncb; j++)
    {
       // if (j==0) cout << "Block " << i << " " << j << endl;
       get_block_trans(i, j, Transf,  ImaTransBlock);
       // cout << " get_block_trans " << endl;

       if (ImaTransBlock.energy() != 0)
       {
          recons_one_block(ImaTransBlock, ImaBlock);
          add_block_ima(i, j, Ima, ImaBlock);
       }
        //INFO_X(ImaTransBlock, "REC BAND BEF");
       //INFO_X(ImaBlock, "REC BAND REC");
       // cout << "  recons_one_block " << endl;

       if ((i==0) && (j==0)) Verbose = False;
    }
    // INFO_X(Ima, "REC BAND");

    Verbose = UserVerbose;
}

/***************************************/

void Ridgelet::variance_stab(Ifloat &Transf)
{
    int i,j;
    double C1 = 2.;
    double C2 = 3. / 8.;
    for (i = 0; i <  Transf.nl(); i++)
    for (j = 0; j <  Transf.nc(); j++)
    {
       if (Transf(i,j) > 0.)
               Transf(i,j) = (float)(C1*sqrt((double) Transf(i,j) + C2));
       else Transf(i,j) = 0.;
    }
}

/***************************************/

void Ridgelet::inv_variance_stab(Ifloat &Transf)
{
   int i,j;
   double C1 = 4.;
   double C2 = 3. / 8.;
   for (i = 0; i <  Transf.nl(); i++)
   for (j = 0; j <  Transf.nc(); j++)
            Transf(i,j) = (float)((double) Transf(i,j)*Transf(i,j)/C1 - C2);
}


/****************************************************************************/

int Ridgelet::set_tab_pos(int NbrScale, int N)
// Set the table position, for band position
{
    int s,Ns=N;
    int TotSize=0;
    if (TabSize.nx() != NbrScale) TabSize.reform(NbrScale);
    if (TabPos.nx() != NbrScale)  TabPos.reform(NbrScale);
    if ((RidTrans == RID_FSS)  || (RidTrans == RID_FSS_PYR) || (RidTrans == RID_FSS_ORTHO) || (RidTrans == RID_UWT) ) Ns *= 2;


    if (RidClass == RID_CL_PAVE)
    {
        for (s = 0; s < NbrScale; s++)
        {
            TabSize(s) = Ns;
            TotSize += TabSize(s);
        }
    }
    else
    {
    if ((RidClass != RID_CL_ORTHO) && (RidClass != RID_CL_FINITE))
    {
      for (s = 0; s < NbrScale; s++)
      {
         TotSize += Ns;
         TabSize(s) = Ns;
         Ns = (Ns+1)/2;
      }
    }
    else
    {
      for (s = 0; s < NbrScale; s++)
      {
         TabSize(s) = Ns/2;
         if (s != NbrScale-1)
         {
            TabSize(s) = Ns/2;
            Ns = (Ns+1)/2;
         }
         else TabSize(s) = Ns;
         TotSize += TabSize(s);
      }
    }
    }

    TabPos(NbrScale-1) = 0;
    for (s =  NbrScale-2; s >= 0; s--)
        TabPos(s) = TabPos(s+1) + TabSize(s+1);

   //  cout << "N = " << N << endl;
   if (Verbose == True)
    for (s = 0; s < NbrScale; s++)
      cout << "Scale " << s+1 << ", Pos = " << TabPos(s) << ", Size = " << TabSize(s) << ", End = " << TabPos(s)+TabSize(s)-1 << endl;

    return TotSize;
}

/****************************************************************************/

void Ridgelet::transform_one_block(Ifloat &Image, Ifloat &Transf)
{
    int i,j,Nlr,Ncr;
    Ifloat TransRadon;
    Icomplex_f TransRadonCF;

    if (Verbose == True)
      cout << "   Ridgelet transform of Block Image " << Image.nl() << " " << Image.nc() << endl;

    switch(RidTrans)
    {
        case RID_FINITE:
        case RID_UWT:
        case RID_OWT:
                 WTinFFT=False;
                 if (Ptr_SB1D == NULL)
                 {
                     cout << "Error: filter bank class not initialized ... " << endl;
                     exit(-1);
                 }
                 break;
       case RID_PYR_FFT:  WTinFFT=True; break;
       case RID_FSS:
       case RID_FSS_PYR:
       case RID_PYR_DIRECT:
       case RID_FSS_ORTHO:
            WTinFFT=False; break;
       case  RID_PAVE_FFT:
         cout << "Error: undecimated transforms are not implemented in this routine ... " << endl;
         exit(-1);
         break;
        default:
         cout << "Error: unknown ridgelet transform ... " << endl; exit(-1);
         break;
    }
    RidClass = rid_class(RidTrans);

    // Radon transform of the input image
    if (WTinFFT == False)
    {
       RADON.transform(Image, TransRadon);
       Nlr = TransRadon.nl();
       Ncr = TransRadon.nc();
    }
    else
    {
       RADON.fft_trans_cf(Image, TransRadonCF);
       Nlr = TransRadonCF.nl();
       Ncr = TransRadonCF.nc();
    }

#if DEGUG_RID
   if (WTinFFT == False) io_write_ima_float("xx_rad.fits", TransRadon);
#endif

    // Variance stablization of the Radon transform
    if (VarianceStab == True)
    {
        if (WTinFFT == False) variance_stab(TransRadon);
	    else
	    {
	       if ((TransRadon.nl() != Nlr) || (TransRadon.nl() != Nlr)) TransRadon.resize(Nlr,Ncr);
	       RADON.invfft_line(TransRadonCF, TransRadon);
	       variance_stab(TransRadon);
	       RADON.fft_line(TransRadon, TransRadonCF);
        }
    }
    else // Normalization of the Radon transform
    {
       // INFO_X(ImaCstRadonTransOneBlock, "ICRB");
       if ((RidTrans == RID_FSS) || (RidTrans ==  RID_FSS_PYR) || (RidTrans ==  RID_FSS_ORTHO) || (RidTrans ==  RID_UWT))
       {
          if ((TransRadon.nl() != ImaCstRadonTransOneBlock.nl()) || (TransRadon.nc() != ImaCstRadonTransOneBlock.nc()))
	      {
	         cout << "Error: pb initialization in Ridgelet::transform_one_block ... " << endl;
	         exit(-1);
	      }
          for (i=0; i < TransRadon.nl(); i++)
	      for (j=0; j < TransRadon.nc(); j++) TransRadon(i,j) /= ImaCstRadonTransOneBlock(i,j);
       }
    }
    // io_write_ima_float("xx_rad_trans.fits", TransRadon);


    // Apply the 1D wavelet transform on each line of the Radon transform
    if (NbrScale > 1)
    {
       int BNlRid = Nlr;
       // For non orthogonal WT transform, the redundancy is 2
       int BNcRid;
       if  ((RidClass == RID_CL_ORTHO) || (RidClass == RID_CL_FINITE)) BNcRid = Ncr;
       else
       {
           if  (RidClass == RID_CL_PAVE) BNcRid = NbrScale*Ncr;
           else BNcRid = 2*Ncr;
	  // if (RidTrans == RID_FSS) BNcRid *= 2;
       }
       if (Verbose == True)
       {
         cout << "   Radon block image size: Nl = " <<  Nlr<< " Nc = " <<  Ncr << endl;
         cout << "   Ridgelet block image size: Nl = " << BNlRid  << " Nc = " <<  BNcRid  << endl;
       }


       if ((Transf.nl() != BNlRid) || (Transf.nc() != BNlRid))
                                      Transf.resize(BNlRid,BNcRid);

        switch(RidTrans)
        {
            case RID_UWT:
                {
                    fltarray Line(Ncr);      // Radon line
                    fltarray WTLine(BNcRid); // Ridgelet line
                    // UndecSubBandFilter USF(*Ptr_SB1D);
                    PAVE_1D_WT UWT(*Ptr_Undec_SB1D);
                    for (i= 0; i < Nlr; i++)
                    {
                        for (j=0; j < Ncr; j++) Line(j) = TransRadon(i,j);
                                // if (i ==128) printf("L128BefTrans: NbrScale = %d,  Min = %f, Max = %f\n", NbrScale, Line.min(), Line.max());
                                // if (i ==128) fits_write_fltarr("xx_line128_trans1.fits", Line);
                        UWT.transform(Line, WTLine, NbrScale);
                                // if (i ==128) fits_write_fltarr("xx_line128_trans2.fits", WTLine);
                                // if (i ==128) printf("L128AfterTrans: Min = %f, Max = %f\n", WTLine.min(), WTLine.max());
                        for (int a=0; a < NbrScale; a++)
                        for (j=0; j < Ncr; j++) Transf(i,j+(NbrScale-1-a)*Ncr) = WTLine(j,a);
                    }
                            // io_write_ima_float("xx_trrid.fits", Transf);
                }
                break;
            case RID_OWT:       // Orthogonal wavelet transform
 	        case RID_FINITE:
            case RID_FSS_ORTHO:
              {
               fltarray Line(Ncr);      // Radon line
	           fltarray WTLine(BNcRid); // Ridgelet line
               if (ColTrans == False)
               {
                   Ortho_1D_WT WT1D(*Ptr_SB1D);
                   for (i= 0; i < Nlr; i++)
                   {
                      for (j=0; j < Ncr; j++) Line(j) = TransRadon(i,j);
                      // if (i ==128) printf("L128BefTrans: NbrScale = %d,  Min = %f, Max = %f\n", NbrScale, Line.min(), Line.max());
 		      // if (i ==128) fits_write_fltarr("xx_line128_trans1.fits", Line);
                      WT1D.transform(Line, WTLine, NbrScale);
 		      // if (i ==128) fits_write_fltarr("xx_line128_trans2.fits", WTLine);
                      // if (i ==128) printf("L128AfterTrans: Min = %f, Max = %f\n", WTLine.min(), WTLine.max());

                      for (j=0; j < BNcRid; j++) Transf(i,j) = WTLine(j);
                   }
		   // io_write_ima_float("xx_trrid.fits", Transf);
               }
               else
               {
                   LineCol WT2D(*Ptr_SB1D);
                   WT2D.transform(TransRadon, Transf, NbrScale);
               }
              }
	      break;
	    case RID_PYR_FFT:
	    case RID_FSS:
              {
                 WT1D_FFT WT1D(WT_PYR_FFT_DIFF_RESOL, NbrScale, Ncr);
                 int NPixTrans = WT1D.pyr_np();

	          // cout << "NPixTrans = " << NPixTrans << endl;
                 fltarray WTLine(NPixTrans);
	         cfarray LineCF(Ncr);   // Radon line of size Ncr

                 for (i= 0; i < Nlr; i++)
                 {
	    		 // cout << "i = " << i << endl;
             		// for (j=0; j < Nc; j++) Line(j) = TransRadon(i,j);
            		 // WT1D.transform(Line, WTLine, NbrScale);
	     	     if (RidTrans != RID_FSS) for (j=0; j < Ncr; j++) LineCF(j) = TransRadonCF(i,j);
	             else
	             {
 		         FFTN_1D F1D;
 		         for (j=0; j < Ncr; j++)  LineCF(j) = complex_f(TransRadon(i,j), 0.);
		         F1D.fftn1d(LineCF);
 	             }
	             WT1D.transform_cf(LineCF, WTLine, NbrScale);

                    if (WPTrans == True)
	            {
		       // Wavelet packets on the first Scale
	               int NpixB1 = WT1D.pyr_size(0);
	               int PosB1 = WT1D.pyr_pos(0);
		       int Nstep = NbrScale-2;
  	               fltarray Band1(NpixB1);
		       fltarray WP(NpixB1);
		       if (Nstep > 0)
		       {
	                  for (j=0; j < NpixB1; j++) Band1(j) = WTLine(PosB1+j);
		          wp1d_mallat_transform (Band1,WP, Nstep, *Ptr_SB1D);
		          for (j=0; j < NpixB1; j++) WTLine(PosB1+j) = WP(j);
		       }
	            }
                    for (j=0; j < NPixTrans; j++) Transf(i,j) = WTLine(j);
                    for (j=NPixTrans; j < BNcRid; j++) Transf(i,j) = 0.;
                 }
                 if (ColTrans == True)
                 {
                    LineCol WT2D(*Ptr_SB1D);
                    WT2D.transform(Transf, 0, NbrScale);
                 }
	      }
	      break;
        case RID_PYR_DIRECT:
 	    case RID_FSS_PYR:
	      {
               fltarray Line(Ncr);      // Radon line
	       fltarray WTLine(BNcRid); // Ridgelet line
               if (ColTrans == False)
               {
                   PYR1D_WT  WT1D(WT_PYR_BSPLINE,NbrScale,Ncr);
                   WT1D.ModifiedPWT = True;
                   for (i= 0; i < Nlr; i++)
                   {
                      for (j=0; j < Ncr; j++) Line(j) = TransRadon(i,j);
                      // if (i ==128) printf("L128BefTrans: NbrScale = %d,  Min = %f, Max = %f\n", NbrScale, Line.min(), Line.max());
 		      // if (i ==128) fits_write_fltarr("xx_line128_trans1.fits", Line);
                      WT1D.transform(Line, WTLine);
 		      // if (i ==128) fits_write_fltarr("xx_line128_trans2.fits", WTLine);
                      // if (i ==128) printf("L128AfterTrans: Min = %f, Max = %f\n", WTLine.min(), WTLine.max());

                      for (j=0; j < BNcRid; j++) Transf(i,j) = WTLine(j);
                   }
		   // io_write_ima_float("xx_trrid.fits", Transf);
               }
               else
               {
                   LineCol WT2D(*Ptr_SB1D);
                   WT2D.transform(TransRadon, Transf, NbrScale);
               }
              }
	      break;
 	      break;
	    default:
	      cout << "Error: this ridgelet transform is not defined for this routine decomposition routine... " << endl;
	      break;
        }
   }
   else // no wavelet decomposition
   {
      if ((RidTrans == RID_PYR_FFT) || (RidTrans ==  RID_PAVE_FFT))
      {
         TransRadon.resize(Nlr,Ncr);
         RADON.invfft_line(TransRadonCF, TransRadon);
         // INFO_X(TransRadon, "TR");
      }
      Transf = TransRadon;
   }
}


/****************************************************************************/

void Ridgelet::transform(Ifloat &Image, Ifloat &Transf, int BlockSize)
{
   alloc(Image.nl(), Image.nc(), BlockSize);
   block_transform(Image, Transf);
   if (StatInfo == True) info_stat(Transf);
}

/****************************************************************************/

void Ridgelet::transform(Ifloat &Image, Ifloat &Transf)
{
   block_transform(Image, Transf);
   if (StatInfo == True) info_stat(Transf);
}

/****************************************************************************/

void Ridgelet::recons_one_block(Ifloat &Transf, Ifloat &Image)
{
    int i,j;
    int Nl = Transf.nl();
    int Nc = Transf.nc();

    // cout << "NbrScale = " << NbrScale << endl;
    if (NbrScale < 2) RADON.recons(Transf, Image);
    else
    {
    switch(RidTrans)
    {
        case RID_FINITE:
        case RID_OWT:
              WTinFFT=False;
              if (Ptr_SB1D == NULL)
              {
                     cout << "Error: filter bank class not initialized ... " << endl;
                     exit(-1);
              }
              break;
        case RID_PYR_FFT:  WTinFFT=True;  break;
        case RID_FSS:
        case RID_FSS_ORTHO:
        case RID_FSS_PYR:
        case RID_UWT:
        case RID_PYR_DIRECT:
            WTinFFT=False; break;
        case RID_PAVE_FFT:
         cout << "Error: undecimated transforms are not implemented in this routine ... " << endl;
         exit(-1);
        default:
         cout << "Error: unknown ridgelet transform ... " << endl; exit(-1);
    }
    RidClass = rid_class(RidTrans);

    switch (RidTrans)
    {
      case RID_OWT:
      case RID_FINITE:
      case RID_FSS_ORTHO:
        {
           fltarray WTLine(Nc);
           fltarray RecLine(Nc);

           if (Verbose == True)
           {
             cout << "  Ridgelet block reconstruction ..." << endl;
             cout << "  Input block image size: Nl = " << Nl << " Nc = " << Nc << endl;
           }
           if (NbrScale < 2)
           {
             cout << "Error: number of scales must be larger than 2 ... " << endl;
             exit(-1);
           }
           if (Verbose == True) cout << "   Number of scales = " << NbrScale << endl;
#if DEGUG_RID
           io_write_ima_float("xx_rec_rid.fits", Transf);
#endif
          if (NbrScale > 1)
          {
            if (ColTrans == False)
            {
               Ortho_1D_WT WT1D(*Ptr_SB1D);
	       FFTN_1D F1D;

               for (i= 0; i < Nl; i++)
               {
                for (j=0; j < Nc; j++) WTLine(j) = Transf(i,j);
 		// if (i ==128) fits_write_fltarr((char*)"xx_line128_rec1.fits", WTLine);
        //       if (i ==128) printf((char*)"L128BefRecTrans:NbrScale = %d,  Min = %f, Max = %f\n", NbrScale, WTLine.min(), WTLine.max());
                WT1D.recons(WTLine, RecLine, NbrScale);
		// if (i ==128) fits_write_fltarr((char*)"xx_line128_rec2.fits", RecLine);
		// if (i ==128) printf((char*)"L128Rec: Min = %f, Max = %f\n", RecLine.min(), RecLine.max());
                for (j=0; j < Nc; j++)  Transf(i,j) = RecLine(j);
             }
           }
           else
           {
              LineCol WT2D(*Ptr_SB1D);
              WT2D.recons(Transf, NbrScale);
           }
          }
          if (VarianceStab == True) inv_variance_stab(Transf);
          else
          {
             if (RidTrans == RID_FSS_ORTHO)
             {
	        if ((Transf.nl() != ImaCstRadonTransOneBlock.nl()) || (Transf.nc() != ImaCstRadonTransOneBlock.nc()))
	        {
	           cout << "Error: pb initialization in Ridgelet::recons_one_block ... " << endl;
	           exit(-1);
	        }
                for (i=0; i < Transf.nl(); i++)
	        for (j=0; j < Transf.nc(); j++) Transf(i,j) *= ImaCstRadonTransOneBlock(i,j);
             }
          }
 #if DEGUG_RID
       // RADON.Verbose = Verbose;
       io_write_ima_float("xx_rec_rad.fits", Transf);
#endif
      // cout << "Radon method = " <<  StringRadon(RADON.radon_method()) << endl;

       RADON.recons(Transf, Image);
      }
      break;
    case RID_PYR_FFT:
    case RID_FSS:
     {
       Ifloat TransRadon;
       Icomplex_f TransRadonCF;
       cfarray RecLineCF(Nc/2) ;
                               // RecLine is one line of the Radon transform
                               // FFT WT1D has a factor 2 of redundancy
                               // The radon transform has half the number
			       // of columns than the Ridglet transform
       TransRadonCF.alloc(Nl,Nc/2,"trans radon");
       if (Verbose == True)
       {
          cout << "   Ridgelet block reconstruction ..." << endl;
          cout << "   Input block image size: Nl = " << Transf.nl() << " Nc = " << Transf.nc() << endl;
          cout << "   Number of scales = " << NbrScale << endl;
       }
       WT1D_FFT WT1D(WT_PYR_FFT_DIFF_RESOL, NbrScale, Nc/2);
       int NPixTrans = WT1D.pyr_np();
       fltarray WTLine(NPixTrans);
        if (ColTrans == True)
        {
            LineCol WT2D(*Ptr_SB1D);
            WT2D.recons(Transf, 0, NbrScale);
        }

       if (RidTrans == RID_FSS) TransRadon.alloc(Nl,Nc/2,"trans radon");
       for (i= 0; i < Nl; i++)
       {
          for (j=0; j < NPixTrans; j++) WTLine(j) = Transf(i,j);
          if (WPTrans == True)
	  {
	     // Wavelet packets on the first Scale
	     int NpixB1 = WT1D.pyr_size(0);
	     int PosB1 = WT1D.pyr_pos(0);
             int Nstep = NbrScale-2;
  	     fltarray Band1(NpixB1);
             fltarray WP(NpixB1);
	     if (Nstep > 0)
	     {
                 for (j=0; j < NpixB1; j++) WP(j) = WTLine(PosB1+j);
                 wp1d_mallat_rec (WP, Band1, Nstep, *Ptr_SB1D);
                 for (j=0; j < NpixB1; j++) WTLine(PosB1+j) = Band1(j);
	     }
	  }
 	  WT1D.recons_cf(WTLine, RecLineCF, NbrScale);
          if (RidTrans != RID_FSS) for (j=0; j < Nc/2; j++) TransRadonCF(i,j) = RecLineCF(j);
	  else
	  {
             FFTN_1D F1D;
	     F1D.fftn1d(RecLineCF, True);
 	     for (j=0; j < Nc/2; j++)  TransRadon(i,j) = RecLineCF(j).real();
 	  }
       }
#if DEGUG_RID
       // RADON.Verbose = Verbose;
       if (RidTrans == RID_FSS) io_write_ima_float("xx_rad_rec.fits", TransRadon);
#endif

       if (VarianceStab == True)
       {
           if (RidTrans == RID_FSS) inv_variance_stab(TransRadon);
	   else
	   {
	      TransRadon.alloc(Nl,Nc/2,"trans radon");
 	      RADON.invfft_line(TransRadonCF, TransRadon);
              inv_variance_stab(TransRadon);
	      RADON.fft_line(TransRadon, TransRadonCF);
	   }
       }

        //  cout << "Rec radon " << TransRadonCF.nl() << " " << TransRadonCF.nc() << endl;
        // RADON.recons(TransRadon, Image);
	// cout << "Radon method RID_FSS OR FFT = " <<  StringRadon(RADON.radon_method()) << endl;

       if (RidTrans != RID_FSS) RADON.fft_rec_cf(TransRadonCF, Image);
       else if (RidTrans == RID_FSS)
       {
          if ((TransRadon.nl() != ImaCstRadonTransOneBlock.nl()) || (TransRadon.nc() != ImaCstRadonTransOneBlock.nc()))
	  {
	     cout << "Error: pb initialization in Ridgelet::transform_one_block ... " << endl;
	     exit(-1);
	  }
          for (i=0; i < TransRadon.nl(); i++)
	  for (j=0; j < TransRadon.nc(); j++) TransRadon(i,j) *= ImaCstRadonTransOneBlock(i,j);
          // cout << "FSS Radon Rec " << TransRadon.nl() << " " << TransRadon.nc() << " " << Image.nl() << " " << Image.nc() << endl;
          // cout << "RADON: " << RADON.radon_nl() << " " << RADON.radon_nc() << endl;
	  // cout << "RADONI: " << RADON.imag_nl() << " " << RADON.imag_nc() << endl;
	   RADON.recons(TransRadon, Image);
        }
      }
      break;
    case RID_FSS_PYR:
    case RID_PYR_DIRECT:
        {
          Ifloat TransRadon(Nl,Nc/2);
	  fltarray RecLine(Nc/2);

           if (Verbose == True)
           {
             cout << "  Ridgelet block reconstruction ..." << endl;
             cout << "  Input block image size: Nl = " << Nl << " Nc = " << Nc << endl;
           }
           if (NbrScale < 2)
           {
             cout << "Error: number of scales must be larger than 2 ... " << endl;
             exit(-1);
           }
           if (Verbose == True) cout << "   Number of scales = " << NbrScale << endl;
#if DEGUG_RID
           io_write_ima_float("xx_rec_rid.fits", Transf);
#endif
          if (NbrScale > 1)
          {
            if (ColTrans == False)
            {
               PYR1D_WT WT1D(WT_PYR_BSPLINE, NbrScale, Nc/2);
	           WT1D.ModifiedPWT = True;
               int NPixTrans = WT1D.pyr_np();
               fltarray WTLine(NPixTrans);
               for (i= 0; i < Nl; i++)
               {
                for (j=0; j < NPixTrans; j++) WTLine(j) = Transf(i,j);
 		//if (i ==128) fits_write_fltarr("xx_line128_rec1.fits", WTLine);
        //        if (i ==128) printf("L128BefRecTrans:NbrScale = %d,  Min = %f, Max = %f\n", NbrScale, WTLine.min(), WTLine.max());
                WT1D.recons(WTLine, RecLine);
		//if (i ==128) fits_write_fltarr("xx_line128_rec2.fits", RecLine);
		//if (i ==128) printf("L128Rec: Min = %f, Max = %f\n", RecLine.min(), RecLine.max());
                for (j=0; j < Nc; j++)  TransRadon(i,j) = RecLine(j);
             }
           }
           else
           {
              LineCol WT2D(*Ptr_SB1D);
              WT2D.recons(Transf, NbrScale);
	      TransRadon = Transf;
           }
          }
          if (VarianceStab == True) inv_variance_stab(Transf);
          else if (RidTrans == RID_FSS_PYR)
          {
 	        if ((TransRadon.nl() != ImaCstRadonTransOneBlock.nl()) || (TransRadon.nc() != ImaCstRadonTransOneBlock.nc()))
	        {
	           cout << "Error: pb initialization in Ridgelet::recons_one_block ... " << endl;
	           exit(-1);
	        }
		for (i=0; i < TransRadon.nl(); i++)
	        for (j=0; j < TransRadon.nc(); j++) TransRadon(i,j) *= ImaCstRadonTransOneBlock(i,j);
          }
 #if DEGUG_RID
       // RADON.Verbose = Verbose;
       io_write_ima_float("xx_rec_rad.fits", TransRadon);
#endif
       // cout << "Radon method = " <<  StringRadon(RADON.radon_method()) << endl;
       RADON.recons(TransRadon, Image);
      }
 	      break;
    case RID_UWT:
        {
            int NcRadon = Nc/NbrScale;
            int NlRadon = Nl;
            Ifloat TransRadon(NlRadon, NcRadon);
            fltarray UWT_TransLine(NcRadon, NbrScale);
            fltarray RecLine(NcRadon);

            if (Verbose == True)
            {
                cout << "  Ridgelet block reconstruction ..." << endl;
                cout << "  Input block image size: Nl = " << Nl << " Nc = " << Nc << endl;
                cout << "  Radon block image size: Nl = " << NlRadon << " Nc = " << NcRadon << endl;
            }
            if (NbrScale < 2)
            {
                cout << "Error: number of scales must be larger than 2 ... " << endl;
                exit(-1);
            }
            if (Verbose == True) cout << "   Number of scales = " << NbrScale << endl;

            if (NbrScale > 1)
            {
                PAVE_1D_WT UWT(*Ptr_Undec_SB1D);
                // int NPixTrans = Transf.nc();
                for (i= 0; i < Nl; i++)
                {
                    // if (i ==0)  cout << "  Rec0 : Nl = " << UWT_TransLine.nx() << " Nc = " << UWT_TransLine.nc() << endl;

                    for (int a=0; a < NbrScale; a++)
                    for (j=0; j < NcRadon; j++) UWT_TransLine(j,a) = Transf(i, j+(NbrScale-1-a)*NcRadon);

                    //  if (i ==0)  cout << "  Rec1 : Nl = " << Transf.nx() << " Nc = " << Transf.nc() << endl;
                    UWT.recons(UWT_TransLine, RecLine, NbrScale);
                    // if (i ==0)  cout << "  Rec2 : Nl = " << RecLine.nx() << " Nc = " << Transf.nc() << endl;

                    for (j=0; j < NcRadon; j++)  TransRadon(i,j) = RecLine(j);
                }
            }
            // cout << "Radon method = " <<  StringRadon(RADON.radon_method()) << endl;


            if (VarianceStab == True) inv_variance_stab(TransRadon);
            else
            {
                if ((TransRadon.nl() != ImaCstRadonTransOneBlock.nl()) || (TransRadon.nc() != ImaCstRadonTransOneBlock.nc()))
                {
                    cout << "Error: pb initialization in Ridgelet::transform_one_block ... " << endl;
                    exit(-1);
                }
                for (i=0; i < TransRadon.nl(); i++)
                for (j=0; j < TransRadon.nc(); j++) TransRadon(i,j) *= ImaCstRadonTransOneBlock(i,j);
            }

            // io_write_ima_float("xx_rad_rec.fits", TransRadon);
            // cout << "Radon method = " <<  StringRadon(RADON.radon_method()) << endl;
            RADON.recons(TransRadon, Image);
        }
        break;

      default:
	      cout << "Error: this ridgelet transform is not defined for this routine reconstruction routine... " << endl;
	      break;
    } // endswitch
  } // endif
}

/****************************************************************************/

void Ridgelet::recons(Ifloat &Transf, Ifloat &Image, int BlockSize)
{
   set_block_param_from_transform(Transf.nl(), Transf.nc(), BlockSize);
   block_recons(Transf, Image);
}

/****************************************************************************/

void Ridgelet::recons(Ifloat &Transf, Ifloat &Image)
{
    if (AllocatedClass == False)
    {
       cout << "Error: The Ridgelet is not allocated ... " << endl;
       exit(-1);
    }
    if ((NlIma == 0) || (NcIma == 0))
    {
       cout << "Error: Ridgelet Class not initialized ... " << endl;
       exit(-1);
    }
    block_recons(Transf, Image);
}

/****************************************************************************/

void Ridgelet::recons(Ifloat &Transf, Ifloat &Image, int BlockSize, int InputImaNl, int InputImaNc)
{
   alloc(InputImaNl, InputImaNc, BlockSize);
   block_recons(Transf, Image);
}

/****************************************************************************/

void Ridgelet::info_stat(Ifloat &Transf)
{
    int s;
    double Mean, Sigma, Skew, Curt;
    float  Min, Max;
    Ifloat Buff(rid_nl(), rid_nc(), "Scale");

    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale; s++)
        {
           get_scale(Transf, Buff,s);
           // cout << "WT1D Band " << s+1 << " Sigma " << Buff.sigma() << endl;
           int N = Buff.nl()*Buff.nc();
           moment4(Buff.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
	   printf("  Band %d: Nl = %d, Nc = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f\n",
	           s+1, Buff.nl(), Buff.nc(), Sigma, Skew, Curt);
        }
    }
}

/****************************************************************************/

void Ridgelet::iter_recons(Ifloat & RidIma, Ifloat &Sol, float Lamda)
{
   int i,j,Iter;
   int Nl = RidIma.nl();
   int Nc = RidIma.nc();
   Ifloat RidImaAux(Nl,Nc,"Rid aux");
   // type_border Bord=I_CONT;
   // RegulIma RIM;
   int Nlima = Sol.nl();
   int Ncima = Sol.nc();
   MR_Regul RIM;
   RIM.ExpDecreasingLambda = True;
   RIM.NbrScale =  max_scale_number(MIN(Nlima,Ncima));
   if (Verbose == True)
      cout << "TV regul: Number of scales = " << RIM.NbrScale << " Lambda = " << RegulParam*Lamda << endl;

   Bool UserVerbose = Verbose;
   RidImaAux = RidIma;
   recons_one_block(RidImaAux, Sol);
   threshold(Sol);
   RIM.im_soft_threshold(Sol, Sol, RegulParam*Lamda);

   Verbose = False;
   float SigmaOld;
   float Sigma=Sol.sigma();
   float Delta=Sigma;
   //io_write_ima_float("xx_rid.fits", RidIma);
   //INFO_X(RidIma, "RIDIMA");
   for (Iter = 0; Iter < NbrFilterIter; Iter++)
   {
       transform_one_block(Sol,RidImaAux);
       for (i=0; i < Nl; i++)
       for (j=0; j < Nc; j++)
       {
           if (ABS(RidIma(i,j)) > FLOAT_EPSILON) RidImaAux(i,j) = RidIma(i,j);
       }
       recons_one_block(RidImaAux,Sol);
       RIM.im_soft_threshold(Sol, Sol, RegulParam*Lamda);
       threshold(Sol);
       SigmaOld=Sigma;
       Sigma=Sol.sigma();
       Delta = ABS(SigmaOld-Sigma) / Sigma;
       if (UserVerbose == True)
              printf("Iter %d: Delta = %f\n", Iter + 1, Delta);
       //io_write_ima_float("xx_sol.fits", Sol);
       //INFO_X(Sol, "SOL");
   }
   Verbose = UserVerbose;
}

/****************************************************************************/
//  RIDGELET DATA Structure management
/****************************************************************************/

float Ridgelet::get_coef(Ifloat &Trans, int i, int j, int s, int Bi, int Bj)
{
   int Indi = 0;
   int Indj = 0;

   if (NbrScale <= 1)
   {
       cout << "Error in get_coef: NbrScale = " << NbrScale << endl;
       exit(-1);
   }
   else
   {
       Indi = ipos(s, Bi) + i;
       Indj = jpos(s, Bj) + j;
   }
   return  Trans(Indi,Indj);
}

/****************************************************************************/

void Ridgelet::put_coef(Ifloat &Trans, int i, int j, float Val, int s, int Bi, int Bj)
{
   int Indi = 0;
   int Indj = 0;

   if (NbrScale <= 1)
   {
      cout << "Error in get_coef: NbrScale = " << NbrScale << endl;
      exit(-1);
   }
   else
   {
      Indi = ipos(s,Bi) + i;
      Indj = jpos(s,Bj) + j;
   }
   Trans(Indi,Indj) = Val;
}

/****************************************************************************/

void Ridgelet::get_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s)
{
   int i,j;
   int Nls = size_scale_nl(s);
   int Ncs = size_scale_nc(s);
   int Posi = ipos(s);
   int Posj = jpos(s);

   if ((ImaScale.nl() != Nls) || (ImaScale.nc() != Ncs))
                                             ImaScale.resize(Nls,Ncs);

   for (i=0; i < Nls; i++)
   for (j=0; j < Ncs; j++) ImaScale(i,j) = ImaTrans(Posi+i,Posj+j);
}


/****************************************************************************/

void Ridgelet::get_scale_without_bord(Ifloat &ImaTrans, int s, Ifloat &ImaScale)
{
   int i,j,Ind=0,bi,bj;
   // int Nls = size_scale_nl(s);
   // int Ncs = size_scale_nc(s);
   int NlOneBlock = rid_one_block_nl();
   int NcOneBlock = rid_size(s);
   int Posi = ipos(s);
   int Posj = jpos(s);
   int Nl1 = (rid_block_nl()-2)*NlOneBlock;
   int Nc1 = (rid_block_nc()-2)*NcOneBlock;
   // int Np = Nl1*Nc1;
   ImaScale.resize(Nl1,Nc1);

   for (i=0; i < NlOneBlock; i++)
   for (j=0; j < NcOneBlock; j++)
   {
     for (bi=1; bi < rid_block_nl()-1; bi++)
     for (bj=1; bj < rid_block_nc()-1; bj++)
     {
        int IndI = i+bi*NlOneBlock+Posi;
        int IndJ = j+bj*NcOneBlock+Posj;
        if ((IndI < 0) || (IndJ < 0)
          || (IndI >= ImaTrans.nl()) || (IndJ >= ImaTrans.nc()))
        {
	    printf("Error: IndI = %d, IndJ = %d, Nl = %d, Nc = %d\n",
	                 IndI,IndJ,ImaTrans.nl(),ImaTrans.nc());
	    exit(-1);
	}
        ImaScale(Ind++) = ImaTrans(IndI,IndJ);
     }
  }
  // if (Ind != Np)
  //{
  //   cout << "Min = " << ImaScale.min() << " Max = " << ImaScale.max() << endl;
  //   cout << "Size = " << ImaScale.n_elem() << " Sigma = " << ImaScale.sigma() << endl;
  //}
  // INFO_X(ImaScale, "ImaScale");
}

/****************************************************************************/

void Ridgelet::set_tab_angle()
{
   TabAngle.resize(2*BS);
   double D = (int)(BS/2);
   for (int t=0; t <= BS/2; t++)
   {
      double A = atan( (double) (t) / D);
      TabAngle(t) = A;
      TabAngle(BS-t) = PI/2. - A;
      TabAngle(BS+t) = PI/2. + A;
      if (t > 0) TabAngle(2*BS -t) = PI - A;
    }
}

/****************************************************************************/

float Ridgelet::get_angle_direction(int NumDirect, Bool InDegree)
{
    if ((NumDirect < 0) || (NumDirect >= 2*BS))
    {
        cout << "Error: bad direction number (" << NumDirect << "). Must be < " << 2*BS << endl;
        exit(-1);
    }
    if (InDegree == True) return ((float) (TabAngle(NumDirect) / PI * 180.));
    else return TabAngle(NumDirect);
}

/****************************************************************************/

int Ridgelet::get_num_direction(float A, Bool InDegree)
{
   float Angle = (InDegree == False) ? A: (float) (A / PI * 180.);
   float MinDist = ABS(TabAngle(0)-Angle);
   int IndMin = 0;
   for (int t=0; t < 2*BS; t++)
   {
      if (ABS(TabAngle(t)-Angle) < MinDist)
      {
         MinDist = ABS(TabAngle(t)-Angle);
	 IndMin = t;
      }
   }
   return IndMin;
}

/****************************************************************************/

void Ridgelet::get_scale_direction(Ifloat &ImaTrans, Ifloat &ImaScale, int s, float Angle)
{
   int Ind = get_num_direction(Angle);
   int Nls = rid_block_nl();
   int Ncs = size_scale_nc(s);
   if (Verbose == True)
   {
     cout << "Angle " << Angle << " Line Number = " << Ind+1 << " Exact Angle = " << TabAngle(Ind)*180./PI << endl;
     cout << "Size = " << Nls << " " << Ncs << endl;
   }
   get_scale_direction(ImaTrans,ImaScale, s, Ind);
}

/****************************************************************************/

void Ridgelet::get_scale_direction(Ifloat &ImaTrans, Ifloat &ImaScale, int s, int NumDirect)
{
   int i,j;
   int Nls = rid_block_nl();
   int Ncs = size_scale_nc(s);
   int Posi = ipos(s);
   int Posj = jpos(s);


   if ((ImaScale.nl() != Nls) || (ImaScale.nc() != Ncs))
                                             ImaScale.resize(Nls,Ncs);
   Posi+= NumDirect;
   for (i=0; i < Nls; i++)
   for (j=0; j < Ncs; j++)
     ImaScale(i,j) = ImaTrans(Posi+i*NlBlockTransSize,Posj+j);
}

/****************************************************************************/

void Ridgelet::put_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s)
{
   int i,j;
   int Nls = size_scale_nl(s);
   int Ncs = size_scale_nc(s);
   int Posi = ipos(s);
   int Posj = jpos(s);

   if ((ImaScale.nl() != Nls) || (ImaScale.nc() != Ncs))
   {
       cout << "Error: the scale image has not the correct size ... " << endl;
       exit(-1);
   }
   for (i=0; i < Nls; i++)
   for (j=0; j < Ncs; j++)  ImaTrans(i+Posi,Posj+j) = ImaScale(i,j);
}

/****************************************************************************/

void Ridgelet::get_blockvect(Ifloat &ImaTrans, int s,
                             int Indi, int Indj, fltarray & VectPix)
{
   int Ind,bi,bj;
   int NlOneBlock = rid_one_block_nl();
   int NcOneBlock = rid_size(s);

   VectPix.reform(rid_block_nbr());
   if (( Indi < 0) || (Indi >= NlOneBlock) || (Indj < 0) || (Indj >= NcOneBlock))
   {
      cout << "Error: bad pixel position in get_blockvect ... " << endl;
      cout << "        Posi = " << Indi <<  " Posj = " << Indj << endl;
      cout << "        NlOneBlock = " << NlOneBlock <<  " NcOneBlock = " << NcOneBlock << endl;
      exit(-1);
   }
   int Posi = ipos(s);
   int Posj = jpos(s);
   Ind = 0;
   for (bi=0; bi < rid_block_nl(); bi++)
   for (bj=0; bj < rid_block_nc(); bj++)
   {
      int IndI = Indi+bi*NlOneBlock+Posi;
      int IndJ = Indj+bj*NcOneBlock+Posj;
      if ((IndI < 0) || (IndJ < 0)
          || (IndI >= ImaTrans.nl()) || (IndJ >= ImaTrans.nc()))
      {
	    printf("Error: IndI = %d, IndJ = %d, Nl = %d, Nc = %d\n",
	                 IndI,IndJ,ImaTrans.nl(),ImaTrans.nc());
	    exit(-1);
      }
      VectPix(Ind++) = ImaTrans(IndI,IndJ);
   }
}

/****************************************************************************/

void Ridgelet::get_anglevect(Ifloat &ImaTrans, int s, int Indi,
                             fltarray & VectPix)
{
   int j,Ind,bi,bj;
   int NlOneBlock = rid_one_block_nl();
   int NcOneBlock = rid_size(s);
   if (( Indi < 0) || (Indi >= NlOneBlock))
   {
      cout << "Error: bad pixel position in get_blockvect ... " << endl;
      cout << "        Posi = " << Indi <<  endl;
      cout << "        NlOneBlock = " << NlOneBlock <<  endl;
      exit(-1);
   }

   int Posi = ipos(s);
   int Posj = jpos(s);
   //cout << "rid_block_nbr() = " << rid_block_nbr() << endl;
   //cout << "NcOneBlock = " << NcOneBlock << endl;

   VectPix.reform(rid_block_nbr()*NcOneBlock);
   Ind = 0;
   for (bi=0; bi < rid_block_nl(); bi++)
   for (bj=0; bj < rid_block_nc(); bj++)
     for (j=0; j < NcOneBlock; j++)
     {
        int IndI = Indi+bi*NlOneBlock+Posi;
        int IndJ = j+bj*NcOneBlock+Posj;
        if ((IndI < 0) || (IndJ < 0)
          || (IndI >= ImaTrans.nl()) || (IndJ >= ImaTrans.nc()))
        {
	    printf("Error: IndI = %d, IndJ = %d, Nl = %d, Nc = %d\n",
	                 IndI,IndJ,ImaTrans.nl(),ImaTrans.nc());
	    exit(-1);
	}
        VectPix(Ind++) = ImaTrans(IndI,IndJ);
     }
}

/****************************************************************************/

void Ridgelet::put_anglevect(Ifloat &ImaTrans, int s, int Indi,
                             fltarray & VectPix)
{
   int j,Ind,bi,bj;
   int NlOneBlock = rid_one_block_nl();
   int NcOneBlock = rid_size(s);
   if (( Indi < 0) || (Indi >= NlOneBlock))
   {
      cout << "Error: bad pixel position in get_blockvect ... " << endl;
      cout << "        Posi = " << Indi <<  endl;
      cout << "        NlOneBlock = " << NlOneBlock <<  endl;
      exit(-1);
   }

   int Posi = ipos(s);
   int Posj = jpos(s);
   //cout << "rid_block_nbr() = " << rid_block_nbr() << endl;
   //cout << "NcOneBlock = " << NcOneBlock << endl;

   Ind = 0;
   for (bi=0; bi < rid_block_nl(); bi++)
   for (bj=0; bj < rid_block_nc(); bj++)
     for (j=0; j < NcOneBlock; j++)
     {
        int IndI = Indi+bi*NlOneBlock+Posi;
        int IndJ = j+bj*NcOneBlock+Posj;
        if ((IndI < 0) || (IndJ < 0)
          || (IndI >= ImaTrans.nl()) || (IndJ >= ImaTrans.nc()))
        {
	    printf("Error: IndI = %d, IndJ = %d, Nl = %d, Nc = %d\n",
	                 IndI,IndJ,ImaTrans.nl(),ImaTrans.nc());
	    exit(-1);
	}
        ImaTrans(IndI,IndJ) = VectPix(Ind++);
     }
}

/****************************************************************************/

void Ridgelet::mad_normalize_per_angle(Ifloat &ImaTrans)
{
   float SigmaMad;

   for (int b=0; b < NbrScale; b++)
   {
      for (int i = 0; i < rid_one_block_nl(); i++)
      {
	 fltarray VectPix;
  	 get_anglevect(ImaTrans, b, i, VectPix);
	 // igmaMad = VectPix.sigma();
	 SigmaMad =  get_sigma_mad(VectPix.buffer(), VectPix.nx());
	 if (SigmaMad != 0)
	   for (int j=0; j < VectPix.nx(); j++) VectPix(j) /= SigmaMad;
	 put_anglevect(ImaTrans, b, i, VectPix);
      }
   }
}

/****************************************************************************/

void Ridgelet::mad_normalize_per_block(Ifloat &ImaTrans, Bool UseSigma)
{
   init_tab_sigma_block(ImaTrans,UseSigma);
   normalize_block(ImaTrans);
}

/****************************************************************************/

void Ridgelet::init_tab_sigma_block(Ifloat &ImaTrans, Bool UseSigma)
{
    for (int b=0; b < NbrScale; b++)
          get_imablocksigma(ImaTrans, b, TabImaSigmaStat[b], UseSigma);
}

/****************************************************************************/

void Ridgelet::get_imablocksigma(Ifloat &ImaTrans, int s, Ifloat &ImaSigmaStat, Bool UseSigma)
{
   int i,j,Ind,bi,bj;
   // int Nls = size_scale_nl(s);
   // int Ncs = size_scale_nc(s);
   int NlOneBlock = rid_one_block_nl();
   int NcOneBlock = rid_size(s);
   int Posi = ipos(s);
   int Posj = jpos(s);
   fltarray VectPix(rid_block_nbr());
   ImaSigmaStat.resize(NlOneBlock,NcOneBlock);
   ImaSigmaStat.init();
   for (i=0; i < NlOneBlock; i++)
   for (j=0; j < NcOneBlock; j++)
   {
     Ind = 0;
     for (bi=0; bi < rid_block_nl(); bi++)
     for (bj=0; bj < rid_block_nc(); bj++)
     {
        int IndI = i+bi*NlOneBlock+Posi;
        int IndJ = j+bj*NcOneBlock+Posj;
        if ((IndI < 0) || (IndJ < 0)
          || (IndI >= ImaTrans.nl()) || (IndJ >= ImaTrans.nc()))
        {
	    printf("Error: IndI = %d, IndJ = %d, Nl = %d, Nc = %d\n",
	                 IndI,IndJ,ImaTrans.nl(),ImaTrans.nc());
	    exit(-1);
	}
        VectPix(Ind++) = ImaTrans(IndI,IndJ);
     }
     if (UseSigma == True) ImaSigmaStat(i,j) = VectPix.sigma();
     else ImaSigmaStat(i,j) = get_sigma_mad(VectPix.buffer(), VectPix.nx());
   }
}

/****************************************************************************/

void Ridgelet::normalize_block(Ifloat &ImaTrans, Bool InverseNormalization)
{
   int i,j,bi,bj;

   for (int s=0; s < NbrScale; s++)
   {
      int NlOneBlock = rid_one_block_nl();
      int NcOneBlock = rid_size(s);
      int Posi = ipos(s);
      int Posj = jpos(s);
      for (i=0; i < NlOneBlock; i++)
      for (j=0; j < NcOneBlock; j++)
      {
         for (bi=0; bi < rid_block_nl(); bi++)
         for (bj=0; bj < rid_block_nc(); bj++)
         {
           int IndI = i+bi*NlOneBlock+Posi;
           int IndJ = j+bj*NcOneBlock+Posj;
	   float NormCoef = (TabImaSigmaStat[s])(i,j);
	   if (NormCoef > 0)
	   {
	      if (InverseNormalization == False) ImaTrans(IndI,IndJ) /= NormCoef;
              else ImaTrans(IndI,IndJ) *= NormCoef;
	   }
	 }
      }
   }
}

/****************************************************************************/



/****************************************************************************/

//               Ridgelet Block manadgement

/*************************************************************************/

float lap_weight(int BlockSize, int Overlap, int k, int B, int Nb)
// BlockSize = block image size
// Overlap = overlap size
// int k = pixel position in the block B
// Nb = Number of blocks
// return a weight value, following pixel position in the block
{
     float ValReturn=1.;
     if (Overlap != 0)
     {
        int Center = Overlap + BlockSize / 2;
        // int Size = BlockSize+2*Overlap;
        if ((B != 0) && (k < Center))
          ValReturn = (float) pow(sin((double) k / (double) Center*PI/2.), 2.);
          // ValReturn = (float) k / (float) Center;
        else if ((B != Nb-1) && (k >= Center))
          ValReturn = (float) pow(cos((double)(k-Center)/ (double) Center *PI/2.), 2.);
          // ValReturn =  (float) (Size - k) / (float) Center;
     }
     return ValReturn;
}

/*************************************************************************/

float frt_lap_weight(int BlockSize, int Overlap, int k, int B, int Nb)
// BlockSize = block image size
// Overlap = overlap size
// int k = pixel position in the block B
// Nb = Number of blocks
// return a weight value, following pixel position in the block
// Idem than previous one, but shifted by 1 when k < Center
{
     float ValReturn=1.;
     if ((Overlap != 0) && (k != BlockSize-1))
     {
        int Center = Overlap + BlockSize / 2;
        // int Size = BlockSize+2*Overlap;
        if ((B != 0) && (k < Center))
          ValReturn = (float) pow(sin((double) (k+1) / (double) Center*PI/2.), 2.);
          // ValReturn = (float) k / (float) Center;
        else if ((B != Nb-1) && (k >= Center))
          ValReturn = (float) pow(cos((double)(k-Center)/ (double) Center *PI/2.), 2.);
          // ValReturn =  (float) (Size - k) / (float) Center;
     }
     return ValReturn;
}

/*************************************************************************/

float Ridgelet::get_weight(int Bi, int Bj, int k, int l)
// return the weight value for pixel position k,l in the block Bi,Bj
{
    float ValReturn=1;

    if (BS  % 2 == 0)
      ValReturn = lap_weight(BlockImaSize, BlockImaOverLap, k, Bi, Nlb)*
                lap_weight(BlockImaSize, BlockImaOverLap, l, Bj, Ncb);
    else
    {
      ValReturn = frt_lap_weight(BlockImaSize, BlockImaOverLap, k, Bi, Nlb)*
                  frt_lap_weight(BlockImaSize, BlockImaOverLap, l, Bj, Ncb);
    }
    return ValReturn;
}

/*************************************************************************/

void Ridgelet::get_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Extract the block (Bi,Bj) from Image and put it in ImaBlock
{
   int k,l;
   int Depi = BlockImaSize*Bi - BlockImaOverLap;
   int Depj = BlockImaSize*Bj - BlockImaOverLap;

   if (WeightBefTrans == False)
   {
     for (k = 0; k < BS; k++)
     for (l = 0; l < BS; l++)
       ImaBlock(k,l) = Ima(Depi+k,Depj+l,I_MIRROR);
   }
   else
   {
     for (k = 0; k < BS; k++)
     for (l = 0; l < BS; l++)
       ImaBlock(k,l) = Ima(Depi+k,Depj+l,I_MIRROR) * get_weight(Bi,Bj,k,l);
   }
}

/*************************************************************************/

void Ridgelet::put_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Put the block (Bi,Bj) ImaBlock in Image
{
   int k,l;
   int Depi = BlockImaSize*Bi - BlockImaOverLap;
   int Depj = BlockImaSize*Bj - BlockImaOverLap;

   for (k = 0; k < BS; k++)
   for (l = 0; l < BS; l++)
         if ((Depi+k >= 0) && (Depi+k < NlIma)
             && (Depj+l >= 0) && (Depj+l < NcIma))
                               Ima(Depi+k,Depj+l) = ImaBlock(k,l);
}

/*************************************************************************/

void Ridgelet::add_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Add the block (Bi,Bj) ImaBlock in Image  with weighted values

{
   int k,l;
   int Depi = BlockImaSize*Bi - BlockImaOverLap;
   int Depj = BlockImaSize*Bj - BlockImaOverLap;

   if (WeightBefTrans == False)
   {
     for (k = 0; k < BS; k++)
     for (l = 0; l < BS; l++)
         if ((Depi+k >= 0) && (Depi+k < NlIma)
             && (Depj+l >= 0) && (Depj+l < NcIma))
                            Ima(Depi+k,Depj+l) += ImaBlock(k,l)*get_weight(Bi,Bj,k,l);
   }
   else
   {
     for (k = 0; k < BS; k++)
     for (l = 0; l < BS; l++)
         if ((Depi+k >= 0) && (Depi+k < NlIma)
             && (Depj+l >= 0) && (Depj+l < NcIma))
                            Ima(Depi+k,Depj+l) += ImaBlock(k,l);
   }
}

/*************************************************************************/

void Ridgelet::get_block_trans(int Bi, int Bj, Ifloat & Transf, Ifloat &  ImaTransBlock)
// Get a block from the ridgelet transformed image
{
   int s,i,j,k,l;
   int Depi = Bi*NlBlockTransSize;
   int Depj = Bj*NcBlockTransSize;

   if (NbrScale <= 1)
   {
      for (k = 0; k < NlBlockTransSize; k++)
      for (l = 0; l < NcBlockTransSize; l++)
             ImaTransBlock(k,l) = Transf(Depi+k, Depj+l);
   }
   else
   {
      for (s=0; s < NbrScale; s++)
      {
          // cout << "Scale = " << s+1 << endl;
          int NFirst= rid_pos(s);
          int Nsize = rid_size(s);
          int Posi = ipos(s,Bi);
          int NFirstTrans= jpos(s,Bj);

          // cout << " ImaTransBlock.nl() = " << ImaTransBlock.nl() << endl;
          // cout << "  NFirst = " << NFirst << endl;
          // cout << "  NFirstTrans  = " <<  NFirstTrans << endl;

          if ((RidTrans == RID_PYR_FFT)  && (VarNorm == True))
	  {
	     // int N = ImaTransBlock.nl();
             // float S2 = sqrt(sqrt(2.))-1.;

             for (i=0; i < ImaTransBlock.nl(); i++)
	     {
 	        float ParamNormAngle =  sqrt((float) BS); // *(1. + (1. - ABS( (i % (N/2)) - N /4.) / (N/4.))*S2);
                for (j=0; j < Nsize; j++)
		{
		  if (j+NFirst >= ImaTransBlock.nc())
	          {
	           cout << "ImaTransBlock : nc = " << ImaTransBlock.nc() << " j+NFirst= " << j+NFirst << " j = " << j << endl;
		   exit(-1);
                  }
		  if ((i+Posi >= Transf.nl()) || (NFirstTrans+j >= Transf.nc()))
	          {
	           cout << "Transf : Nl = " << Transf.nl() << endl;
	           cout << "         Nc = " << Transf.nc() << endl;
	           cout << "         Posi = " << Posi << endl;
	           cout << "         NFirstTrans = " << NFirstTrans << endl;
	           cout << "         i+Posi = " << i+Posi << endl;
	           cout << "         NFirstTrans+j = " << NFirstTrans+j << endl;
	           exit(-1);
                  }
		  ImaTransBlock(i,j+NFirst) = Transf(i+Posi,NFirstTrans+j)*ParamNormAngle;
		}
             }
	  }
	  else
	  {
             for (i=0; i < ImaTransBlock.nl(); i++)
             for (j=0; j < Nsize; j++)
	     {
	       if (j+NFirst >= ImaTransBlock.nc())
	       {
	        cout << "ImaTransBlock : nc = " << ImaTransBlock.nc() << " j+NFirst= " << j+NFirst << " j = " << j << endl;
		exit(-1);
               }
	       if ((i+Posi >= Transf.nl()) || (NFirstTrans+j >= Transf.nc()))
	       {
	        cout << "Transf : Nl = " << Transf.nl() << endl;
 		cout << "         Nc = " << Transf.nc() << endl;
		cout << "         Posi = " << Posi << endl;
		cout << "         NFirstTrans = " << NFirstTrans << endl;
		cout << "         i+Posi = " << i+Posi << endl;
		cout << "         NFirstTrans+j = " << NFirstTrans+j << endl;
	        exit(-1);
               }
               ImaTransBlock(i,j+NFirst) = Transf(i+Posi,NFirstTrans+j);
	     }
	  }
      }
   }
}

/*************************************************************************/

void Ridgelet::put_block_trans(int Bi, int Bj, Ifloat & Transf, Ifloat & ImaTransBlock)
// Put a block from the ridgelet transformed image
{
   int s,i,j,k,l;

   if ((Transf.nl() < ImaTransBlock.nl()) || (Transf.nc() < ImaTransBlock.nc()))
   {
       cout << "Block size < image size " <<  endl;
       cout << "Nbr line: " << ImaTransBlock.nl() << " " << Transf.nl() << endl;
       cout << "Nbr Col: " << ImaTransBlock.nc() << " " << Transf.nc() << endl;
       exit(-1);
   }
  // cout << "Nbr line block and ima: " << ImaTransBlock.nl() << " " << Transf.nl() << endl;
  // cout << "Nbr Col  block and ima: " << ImaTransBlock.nc() << " " << Transf.nc() << endl;

   if (NbrScale <= 1)
   {
     int Depi = Bi*NlBlockTransSize;
     int Depj = Bj*NcBlockTransSize;
     for (k = 0; k < NlBlockTransSize; k++)
     for (l = 0; l < NcBlockTransSize; l++)
            Transf(Depi+k, Depj+l) = ImaTransBlock(k,l);
   }
   else
   {
      for (s=0; s < NbrScale; s++)
      {
          int NFirst= rid_pos(s);
          int Nsize = rid_size(s);
          int Posi = ipos(s,Bi);
          int NFirstTrans= jpos(s,Bj);
          // float S2 = sqrt(sqrt(2.))-1.;
	  // int N=ImaTransBlock.nl();

          //cout << "Scale " << s+1 << " NFirst = " << NFirst << endl;
          //cout << "NFirst end  = " << NFirst  + Nsize -1  << endl;
          //cout << "NFirstTrans = " <<  NFirstTrans << endl;
          //cout << " END NFirstTrans " << NFirstTrans + Nsize -1 << endl;
          if ((RidTrans == RID_PYR_FFT) && (VarNorm == True))
	  {
             for (i=0; i < ImaTransBlock.nl(); i++)
	     {
	        float ParamNormAngle = sqrt((float) BS); // *(1. + (1. - ABS( (i % (N/2)) - N /4.) / (N/4.))*S2);
                for (j=0; j < Nsize; j++)
                            Transf(i+Posi,NFirstTrans+j) = ImaTransBlock(i,j+NFirst) / ParamNormAngle;
             }
	  }
	  else
	  {
	     for (i=0; i < ImaTransBlock.nl(); i++)
             for (j=0; j < Nsize; j++)
                                 Transf(i+Posi,NFirstTrans+j) = ImaTransBlock(i,j+NFirst);
 	  }
      }
   }
}

/*************************************************************************/
/*************************************************************************/
/****************************************************************************/

//               Ridgelet IO

/*************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 4) || (File_Name_In[L-1] != 'd')
                || (File_Name_In[L-2] != 'i')
                || (File_Name_In[L-3] != 'r')
                || (File_Name_In[L-4] != '.'))
    {
        strcat (File_Name_Out, ".rid");
    }
}

/****************************************************************************/
/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/

void Ridgelet::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nl", (long) NlIma,(char*)"NlIma",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Nc",(long) NcIma,(char*)"NcIma",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"BSize",(long) BS,(char*)"BlockImaSize",&status))
     PrintError( status );
  if ( ffpkyj(fptr, (char*)"Nlb", (long)Nlb,(char*)"Number of blocks Nlb",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Ncb",(long)Ncb,(char*)"Number of blocks Ncb",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Overlap",(long) ((BlockOverlap == True) ? 1: 0),(char*)"Block overlap",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"NbrScale",(long)NbrScale,(char*)"Number of scales",&status))
     PrintError( status );

  if ( ffpkyj(fptr, (char*)"Type_Tra", (long) RidTrans ,
                          (char*)StringRidTransform(RidTrans), &status))
     PrintError( status );
  if ( ffpkyj(fptr, (char*)"FormatIn",(long)FormatInputImag,(char*)"Format", &status))
      PrintError( status );
  if ( ffpkyj(fptr,(char*)"VarNorm",(long) ((VarNorm == True) ? 1: 0),(char*)"Var Norm",&status))
     PrintError( status );

}


/****************************************************************************/

void Ridgelet::write (char *Name, Ifloat &Ima)
/* new version with fits */
{
 char filename[256];
 fitsfile *fptr;
 int status;
 //int i,j,s;
 //float *Ptr;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[3];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 //long fpixels[3];
 //long lpixels[3];

/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */

/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 naxis = 2;
 naxes[0] = rid_nc();
 naxes[1] = rid_nl();

 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // write the header of the multiresolution file
  mr_io_fill_header(fptr);

  nelements = naxes[0] * naxes[1];
  if ( ffppre(fptr, group, firstpixel, nelements, Ima.buffer(), &status) )
              PrintError( status );

 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );
// cout << " end of write fits " << endl;

}
/****************************************************************************/

void Ridgelet::read (char *Name, Ifloat & Ima)
{
    // for fits
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT];
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
   // long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];

     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

   // inc[0]=1;  inc[1]=1; inc[2]=1;

#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif

    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) )
         PrintError( status );

    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

     nelements = naxes[0] * naxes[1];
     // cout << " begin to read " << endl;
     Ima.alloc(naxes[1], naxes[0], "read io");

    if (ffgkyj(fptr,"Nl", &mon_long, comment, &status)) PrintError( status );
    int Nli = (int) mon_long;
    if (ffgkyj(fptr,"Nc", &mon_long, comment, &status)) PrintError( status );
    int Nci = (int) mon_long;
    if (ffgkyj(fptr,"BSize", &mon_long, comment, &status)) PrintError( status );
    int Bi = (int) mon_long;
    if (ffgkyj(fptr,"NbrScale", &mon_long, comment, &status)) PrintError( status );
    NbrScale = (int) mon_long;
    GetAutoNbScale = False;
    if (ffgkyj(fptr,"Type_Tra", &mon_long, comment, &status)) PrintError( status );
    RidTrans = (type_ridgelet_WTtrans) mon_long;

    if ( ffgkyj(fptr,"Overlap",&mon_long, comment, &status))
     PrintError( status );
    BlockOverlap  = (mon_long == 0) ? False : True;

    if (ffgkyj(fptr,"FormatIn", &mon_long, comment, &status)) PrintError( status );
    FormatInputImag = (type_format)mon_long;
    if ( ffgkyj(fptr,"VarNorm",&mon_long, comment, &status))
     PrintError( status );
    VarNorm  = (mon_long == 0) ? False : True;

    // if (BlockOverlap == True) cout << "Block Over " << endl;
    // else cout << "NO Block Over " << endl;

    alloc(Nli,Nci,Bi);


#if DEBUG_IO
    cout << "Read in " << filename << endl;
    cout << "Nx = " << size_ima_nl () << endl;
    cout << "Ny = " << size_ima_nc () << endl;
    cout << "Nbr_Plan = " << nbr_scale () << endl;
    cout << "Type_Transform = " << RidTrans << " " <<
             StringRid3DTransform(RidTrans) << endl;
 #endif


    Ptr = Ima.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

  if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;

#if DEBUG_IO
    cout << "Read out " << filename << endl;
#endif
}

/****************************************************************************/
