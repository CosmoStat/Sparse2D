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
**    File:  RidFilter.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Ridgelet Filtering routine
**    -----------
**
******************************************************************************/

#include "Ridgelet.h"
#include "RidNoise.h"

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "MR1D_Filter.h"
#include "MR1D_NoiseModel.h"
 #include "MR_Threshold.h"

#define DEGUG_RID 0
#define NO_POISSON_FILTER 1

// In the case of the a FFT-based transform, the noise variance is
// keeped in the different bands. This table allows us to renormalize
// the noise standard deviation.

const float TabNormRid_PyrDiffResol[RID_MAX_SCALE] =
    {0.8546,0.215, 0.1265,0.097,0.0669,0.0483,0.0342,0.0242,0.0171,0.0121};


/****************************************************************************/

void Ridgelet::set_tab_norm_from_simu(Ifloat &Simu, int NbrScale)
{
    cout << "set_tab_norm_from_simu: Nl = " << Simu.nl() << " " << Simu.nc() << endl;
    TabNorm.reform(NbrScale);
    int i;
    Ifloat SimuRid;
    Bool OldBlockOverlap = BlockOverlap;


    transform(Simu,SimuRid);

    cout << "  alloc: NbrScale = " << NbrScale << endl;
    fltarray TabMin(NbrScale), TabMax(NbrScale);
    fltarray Tab_Nsig(NbrScale);
    Tab_Nsig.init(3.);
    StatRidNoise SRN;
    // SRN.Verbose = True;
    SRN.init(SimuRid, NbrScale, (*this));

    cout << "find_nsigma_threshold .. " << endl;
    SRN.find_nsigma_threshold(NbrScale, Tab_Nsig, TabMin, TabMax);
    for (i=0; i < NbrScale-1; i++)
    {
       float Coef = (RidTrans != RID_FSS) ? Tab_Nsig(i) * sqrt((float) Simu.nc()):  Tab_Nsig(i);
       TabMin(i) = ABS(TabMin(i))/Coef;
       TabMax(i) /= Coef;
       TabNorm(i) = MAX(TabMin(i),TabMax(i));
       cout <<  TabNorm(i) << ", ";
               // cout << i+1<< ": Snsig = " << Tab_Nsig(i) << " Min = " << TabMin(i) << " Max = " << TabMax(i) << endl;
    }
    cout << "OUT set_tab_norm_from_simu " << endl;
    BlockOverlap = OldBlockOverlap;
    SetTable = False;
}

/****************************************************************************/

void Ridgelet::set_tab_norm(int NbrScale)
// Set the noise standard deviation normalization
// In case of an orthognal transform, is is equal to 1 for all scale
{
    if (TabNorm.nx() != NbrScale) TabNorm.reform(NbrScale);
    int s;
    if ((RidTrans == RID_PYR_FFT) || (RidTrans == RID_PAVE_FFT) || (RidTrans == RID_FSS))
    {
       for (s = 0; s < NbrScale; s++)
          if (s < RID_MAX_SCALE) TabNorm(s) = TabNormRid_PyrDiffResol[s];
          else TabNorm(s) = TabNorm(s-1) / sqrt(2.);
    }
    else for (s = 0; s < NbrScale; s++) TabNorm(s) = 1.;
}

/****************************************************************************/

void Ridgelet::set_tab_norm_with_tab(int NbrScale, fltarray & Tab)
{
    if (TabNorm.nx() != NbrScale) TabNorm.reform(NbrScale);
    for (int s = 0; s < NbrScale; s++)
    {
       if (s < Tab.nx()) TabNorm(s) = Tab(s);
       else TabNorm(s) = 0.;
    }
}

/****************************************************************************/

void Ridgelet::ksig_threshold(Ifloat &Transf, float NoiseIma,
                              float N_Sigma, fltarray *TabN)
{
    int i,j,s;
    int Nl = Transf.nl();
    int NFirst,NLast;
    if (TabNorm.n_elem() < NbrScale-1) set_tab_norm(NbrScale);

    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale-1; s++)
        {
           float Nsig = (s == 0) ? N_Sigma + 1: N_Sigma;
           float Level = NoiseIma*Nsig*TabNorm(s);

           // float Level = NoiseIma*N_Sigma*TabNorm(s);
           NFirst= rid_pos(s);
           NLast = rid_pos(s) + rid_size(s);
	   // int Nc = NLast-NFirst;
           for (i=0; i < Nl; i++)
	   {
	       if (TabN != NULL) Level  = NoiseIma*Nsig*(*TabN)(i,s);
	      // float Param  =  1. - ABS( (i % (Nl/2)) - Nl /4.) / (Nl/4.);
	      // float CurLevel = Level*(1.+Param/10.);
              for (j=NFirst; j < NLast; j++)
               if ((s < FirstDetectScale) || (ABS(Transf(i,j)) < Level))
                                                            Transf(i,j) = 0.;
	   }
        }
	// cout << "Last scale " << endl;
	if (KillLastScale == True)
	{
	   s =  NbrScale-1;
           NFirst= rid_pos(s);
           NLast = rid_pos(s) + rid_size(s);
	   // cout << "First = " << NFirst << " Last = " << NLast << endl;
           for (i=0; i < Nl; i++)
           for (j=NFirst; j < NLast; j++) Transf(i,j) = 0.;
	}
    }
}

/****************************************************************************/


/****************************************************************************/

void Ridgelet::normalize_coef(Ifloat &Transf, Bool Inverse, float NoiseIma)
{
    int i,j,s;
    int Nl = Transf.nl();
    int NFirst,NLast;

    if (TabNorm.n_elem() < NbrScale-1) set_tab_norm(NbrScale);
    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale; s++)
        {
           int IS = (s == NbrScale-1) ? s-1: s;
           float NormVal = TabNorm(IS)*NoiseIma;
           if (Inverse == False) NormVal = 1. / NormVal;
           NFirst= rid_pos(s);
           NLast = rid_pos(s) + rid_size(s);
           for (i=0; i < Nl; i++)
           for (j=NFirst; j < NLast; j++) Transf(i,j) *= NormVal;
        }
    }
}

/****************************************************************************/

void Ridgelet::mad_threshold(Ifloat &Transf, float N_Sigma)
{
    int i,j,s;
    int Nl = Transf.nl();
    int NFirst,NLast;
    float Level,Sigma;
    float *Buff;

    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale-1; s++)
        {
           NFirst= rid_pos(s);
           NLast = rid_pos(s) + rid_size(s);
       //  cout << "MAD First = " << NFirst << " NLast = " << NLast <<  " size = " << rid_size(s) << endl;

           if (NoiseInvariantPerRotation == True)
	   {
              int Nelem = Nl*rid_size(s);
	      int Ind=0;
              Buff = new float [Nelem];
              for (i=0; i < Nl; i++)
              for (j=NFirst; j < NLast; j++) Buff[Ind++] = Transf(i,j);
              Sigma = get_sigma_mad(Buff, Nelem);
              Level = N_Sigma*Sigma;
              if (Verbose == True)
                  cout << "Band " << s+1 << " MAD = " << Sigma << " Level = " << Level << endl;
              for (i=0; i < Nl; i++)
              for (j=NFirst; j < NLast; j++)
                       if (ABS(Transf(i,j)) < Level) Transf(i,j) = 0.;
              delete [] Buff;
	   }
	   else
	   {
	      int Nelem = rid_size(s);
              Buff = new float [Nelem];
              cout << "MAD First = " << NFirst << " NLast = " << NLast <<  " size = " << rid_size(s) << endl;
	      for (i=0; i < Nl; i++)
	      {
	         int Ind=0;
	         for (j=NFirst; j < NLast; j++) Buff[Ind++] = Transf(i,j);
		 Sigma = get_sigma_mad(Buff, Nelem);
                 Level = N_Sigma*Sigma;
		 for (j=NFirst; j < NLast; j++)
                       if (ABS(Transf(i,j)) < Level) Transf(i,j) = 0.;
	      }
	      delete [] Buff;
	   }
        }
    }
}

/****************************************************************************/

void Ridgelet::filtering_one_block(Ifloat &Data, float NoiseIma, float N_Sigma, fltarray *TabN)
{
    Ifloat RidIma;
    float SigmaNoise = NoiseIma;

    // Ridgelet transform
    transform_one_block(Data,RidIma);
#if DEGUG_RID
    io_write_ima_float("xx_rid.fits",  RidIma);
#endif

    // we integrate on lines ==> the noise is multipled by sqrt(N)
    // in case of Gaussian noise
    if ((VarianceStab == False) && (RidTrans != RID_FSS))  SigmaNoise *= sqrt((float) Data.nc());

    if (SetTable == True)
    {
       // Thresholding
       if (TabN == NULL) set_tab_norm(NbrScale);
       else  if ((*TabN).naxis() == 1) set_tab_norm_with_tab(NbrScale, *TabN);
    }

    if (MadThreshold == True) mad_threshold(RidIma, N_Sigma);
    else ksig_threshold(RidIma, SigmaNoise, N_Sigma);

#if DEGUG_RID
    io_write_ima_float("xx_fil.fits", RidIma);
#endif
    // Reconstruction
    if (NbrFilterIter > 1) iter_recons(RidIma,Data,NoiseIma);
    else recons_one_block(RidIma,Data);
}

/****************************************************************************/

void Ridgelet::kill_fluctu(Ifloat &Ima, float NoiseIma, float N_Sigma,
                           int BlockSize)
{
    int i,j;
    int B = (BlockSize == 0) ? Ima.nc(): BlockSize;
    float ThresholdLevel = NoiseIma*N_Sigma / sqrt((float) B) / 2.;
    // cout << "Threshold = " << ThresholdLevel << endl;
    for (i = 0; i < Ima.nl(); i++)
    for (j = 0; j < Ima.nc(); j++)
              if (ABS(Ima(i,j)) <  ThresholdLevel) Ima(i,j) = 0.;
}

/****************************************************************************/

void Ridgelet::filtering(Ifloat &Ima, float NoiseIma, float N_Sigma,
                         int BlockSize, fltarray *TabN)
{
    int i,j;
    int Nc = Ima.nc();
    int Nl = Ima.nl();
    alloc(Nl, Nc, BlockSize);

    if ((BlockSize > 0) || (BlockSize != Nl) || (BlockSize != Nc))
    {
       Ifloat ImaBlock(BS, BS,"block");
       Ifloat Filter(Nl,Nc, "Filter");
       if (OnlyHighFreq == True)
       {
           Ifloat ImaHigh(Nl,Nc, "ImaHigh");
           get_low_and_freq(Ima, Filter, ImaHigh, BlockSize);
           Ima = ImaHigh;
           // io_write_ima_float("xx_low.fits", Filter);
           // io_write_ima_float("xx_high.fits",Ima);
       }
       else Filter.init();

       Bool UserVerbose = Verbose;
       // cout << "block Nlb = " << Nlb << " Ncb = " << Ncb << endl;
       for (i = 0; i < Nlb; i++)
       for (j = 0; j < Ncb; j++)
       {
          // if (Verbose == True)
          // cout << BlockSize << "Filter block " << i + 1 << " " << j+1 << endl;
           get_block_ima(i,j, Ima, ImaBlock);
           filtering_one_block(ImaBlock, NoiseIma, N_Sigma, TabN);
           if (KillFluctu == True)
                         kill_fluctu(ImaBlock, NoiseIma, N_Sigma, BlockSize);
           add_block_ima(i,j, Filter, ImaBlock);
           if ((i == 0) && (j == 0)) Verbose = False;
       }
       Verbose = UserVerbose;
       Ima=Filter;
   }
   else filtering_one_block(Ima, NoiseIma, N_Sigma, TabN);
}

/****************************************************************************/

void Ridgelet::radon_filtering_one_block(Ifloat &Data, float NoiseIma, float N_Sigma)
{
    Ifloat RadIma;
    float SigmaNoise = NoiseIma;

    // Ridgelet transform
    RADON.transform(Data, RadIma);
    // INFO_X(RadIma, "Radon image");
#if DEGUG_RID
    io_write_ima_float("xx_rad.fits", RadIma);
#endif

    // we integrate on lines ==> the noise is multipled by sqrt(N)
    // in case of Gaussian noise
    if (VarianceStab == False) SigmaNoise *= sqrt((float) Data.nc());

    // Thresholding
    int i,j;
    float Level = SigmaNoise*N_Sigma;
    // cout << "N_Sigma = " << N_Sigma << " Level = " << Level << endl;
    for (i=0; i <  RadIma.nl(); i++)
    for (j=0; j <  RadIma.nc(); j++)
            if (ABS(RadIma(i,j)) < Level)  RadIma(i,j) = 0.;

#if DEGUG_RID
    io_write_ima_float("xx_fil.fits", RadIma);
#endif

    // Reconstruction
   RADON.recons(RadIma, Data);
}

/****************************************************************************/

void Ridgelet::radon_filtering(Ifloat &Ima, float NoiseIma, float N_Sigma,
                         int BlockSize)
{
    int i,j;
    int Nc = Ima.nc();
    int Nl = Ima.nl();
    alloc(Nl, Nc, BlockSize);

    if (BlockSize > 0)
    {
       Ifloat Filter(Nl,Nc, "Filter");
       Ifloat ImaBlock(BS, BS,"block");

       if (OnlyHighFreq == True)
       {
           Ifloat ImaHigh(Nl,Nc, "ImaHigh");
           get_low_and_freq(Ima, Filter, ImaHigh, BlockSize);
           Ima = ImaHigh;
           // io_write_ima_float("xx_low.fits", Filter);
           // io_write_ima_float("xx_high.fits",Ima);
       }
       else Filter.init();

       Bool UserVerbose = Verbose;
       for (i = 0; i < Nlb; i++)
       for (j = 0; j < Ncb; j++)
       {
          get_block_ima(i,j, Ima, ImaBlock);
          radon_filtering_one_block(ImaBlock, NoiseIma, N_Sigma);
          add_block_ima(i,j, Filter, ImaBlock);
          if ((i == 0) && (j == 0)) Verbose = False;
       }
       Verbose = UserVerbose;
       Ima=Filter;
    }
    else radon_filtering_one_block(Ima, NoiseIma, N_Sigma);
}

/****************************************************************************/

static void smooth_bspline (Ifloat & Im_in, Ifloat &Im_out, type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    double Coeff_h0 = 3. / 8.;
    double Coeff_h1 = 1. / 4.;
    double  Coeff_h2 = 1. / 16.;
    Ifloat Buff(Nl,Nc,"Buff smooth_bspline");

    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
    {
       double Val = Coeff_h0 * (double) Im_in(i,j)
                 + Coeff_h1 * (  Im_in (i, j-Step, Type)
                               + Im_in (i, j+Step, Type))
                 + Coeff_h2 * (  Im_in (i, j-2*Step, Type)
                               + Im_in (i, j+2*Step, Type));
       Buff(i,j) = (float) Val;
    }

    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = Coeff_h0 * Buff(i,j)
                 + Coeff_h1 * (  Buff (i-Step, j, Type)
                               + Buff (i+Step, j, Type))
                + Coeff_h2 * (  Buff (i-2*Step, j, Type)
                               + Buff (i+2*Step, j, Type));
}

/****************************************************************************/

void Ridgelet::get_low_and_freq(Ifloat &Image,  Ifloat &Low, Ifloat &High, int BlockSize)
{
   // Calculate fits the number of scales used in the 1D WT
   int s,i,j,Ns = NbrScale;
   // if (GetAutoNbScale == True) Ns = max_scale_number(BlockSize);
   if (Ns < 4) Ns = 4;
   type_border Border=I_MIRROR;

   // cout << "Rad scale nbr = " << Ns << endl;
   for (i=0;i<Image.nl(); i++)
   for (j=0;j<Image.nc(); j++) High(i,j) =  Image(i,j);

   for (s = 0; s < Ns; s++)
   {
      smooth_bspline (High, Low, Border, s);
      if (s != Ns-1)
      {
         for (i=0;i<Image.nl(); i++)
         for (j=0;j<Image.nc(); j++) High(i,j) = Low(i,j);
      }
   }

   for (i=0;i<Image.nl(); i++)
   for (j=0;j<Image.nc(); j++) High(i,j) =  Image(i,j) -  Low(i,j);
}

/***************************************************************************/

void Ridgelet::undec_wt_filtering_one_block(Ifloat &Image, float NoiseIma, float N_Sigma)
{
    int i,j,s;
    Ifloat Transf;
    float SigmaNoise = NoiseIma;

    switch(RidTrans)
    {
        case RID_OWT:
        case RID_PYR_FFT:
              cout << "Error: decimated transform are not implemented in this routine  ... " << endl;
              exit(-1);
             break;
        case RID_UWT:
             if (Ptr_SB1D == NULL)
              {
                     cout << "Error: filter bank class not initialized ... " << endl;
                     exit(-1);
              }
              WTinFFT=False;
              break;
        case RID_PAVE_FFT:
               WTinFFT=True;
             break;
        default:
         cout << "Error: unknown ridgelet transform ... " << endl; exit(-1);
    }
    RidClass = rid_class(RidTrans);

    if (Verbose == True)  cout << "Ridgelet transform ..." << endl;
    RADON.transform(Image, Transf);
    if (Verbose == True)
    {
      cout << "Ouput image size: Nl = " << Transf.nl() << " Nc = " << Transf.nc() << endl;
      cout << "Number of scales = " << NbrScale << endl;
    }
    // we integrate on lines ==> the noise is multipled by sqrt(N)
    // in case of Gaussian noise
    if (VarianceStab == False) SigmaNoise *= sqrt((float) Image.nc());
    else variance_stab(Transf);
    float Level = SigmaNoise*N_Sigma;
    if (Verbose == True)
    {
       cout << "Detection level = " << Level << endl;
       cout << "Undecimated filtering: NbrScale = " << NbrScale  << endl;
    }
    // if (GetAutoNbScale == True) NbrScale = max_scale_number(Transf.nc());

    if (NbrScale > 1) // TabNormRid_PyrDiffResol
    {
       int Nl = Transf.nl();
       int Nc = Transf.nc();
       fltarray Line(Nc);
       fltarray WTLine(Nc, NbrScale);
       if (RidTrans == RID_UWT)
       {
#if NO_POISSON_FILTER
          PAVE_1D_WT WT1D(*Ptr_SB1D);
          for (i= 0; i < Nl; i++)
          {
             for (j=0; j < Nc; j++) Line(j) = Transf(i,j);
             WT1D.transform(Line, WTLine, NbrScale);
             for (s=0; s < NbrScale-1; s++)
             for (j=0; j < Nc; j++)
                 if ((s < FirstDetectScale) || (ABS(WTLine(j,s)) < Level))  WTLine(j,s) = 0.;
             if (KillLastScale == True)
                for (j=0; j < Nc; j++) WTLine(j,NbrScale-1) = 0.;
             WT1D.recons(WTLine, Line, NbrScale);
             for (j=0; j < Nc; j++) Transf(i,j) = Line(j);
          }
#else
	  FilterAnaSynt FAS;
          FilterAnaSynt *PtrFAS = NULL;
	  MR1DNoiseModel NoiseModel;
	  int Nbr_Plan=NbrScale;
	  float NSigma=N_Sigma;
	  float EpsilonPoisson = (1. - erf((double) NSigma / sqrt((double) 2.)));
	  type_trans_1d Transform=TO1_PAVE_B3SPLINE;
	  type_noise Stat_Noise = NOISE_EVENT_POISSON;
          NoiseModel.alloc(Stat_Noise, Nc, NbrScale, Transform, PtrFAS);
          for (i=0; i < Nbr_Plan; i++) NoiseModel.NSigma[i]=NSigma;
	  for (i=0; i < Nbr_Plan; i++) NoiseModel.TabEps[i] = EpsilonPoisson;
	  NoiseModel.SupIsol = False;
          NoiseModel.OnlyPositivDetect=True;
          NoiseModel.FirstDectectScale = FirstDetectScale;
          type_threshold TypeThreshold = DEF_THRESHOLD;
          NoiseModel.TypeThreshold = TypeThreshold;
	  type_1d_filter Filter = FIL_1D_ITER_THRESHOLD; //  FIL_1D_THRESHOLD;
	  MR1DFiltering CFilter(NoiseModel, Filter);
	  KillLastScale = True;
	  if (KillLastScale == True)
	  {
	     CFilter.KillLastScale = KillLastScale;
	     CFilter.PositivIma = False;
          }
	  CFilter.PositivIma = True;
	  fltarray Res(Nc);
	  for (i= 0; i < Nl; i++)
	  {
	      for (j=0; j < Nc; j++) Line(j) = Transf(i,j);
	      CFilter.filter(Line, Res);
	      for (j=0; j < Nc; j++) Transf(i,j) = Res(j);
	  }
#endif
       }
       else
       {
          WT1D_FFT WT1D(WT_PAVE_FFT, NbrScale, Nc);
          for (i= 0; i < Nl; i++)
          {
             for (j=0; j < Nc; j++) Line(j) = Transf(i,j);
             WT1D.transform(Line, WTLine, NbrScale);
             for (s=0; s < NbrScale-1; s++)
             {
                Level = TabNormRid_PyrDiffResol[s]*SigmaNoise*N_Sigma;
                for (j=0; j < Nc; j++)
                    if ((s < FirstDetectScale) || (ABS(WTLine(j,s)) < Level))  WTLine(j,s) = 0.;
             }
             if (KillLastScale == True)
                for (j=0; j < Nc; j++) WTLine(j,NbrScale-1) = 0.;
             WT1D.recons(WTLine, Line, NbrScale);
             for (j=0; j < Nc; j++) Transf(i,j) = Line(j);
          }
       }
    }
#if DEGUG_RID
    io_write_ima_float("xx_fil.fits", Transf);
#endif
    if (VarianceStab == True) inv_variance_stab(Transf);
    RADON.Verbose = Verbose;
    RADON.recons(Transf, Image);
}

/****************************************************************************/

void Ridgelet::undec_wt_filtering(Ifloat &Ima, float NoiseIma, float N_Sigma, int BlockSize)
{
    int i,j;
    int Nc = Ima.nc();
    int Nl = Ima.nl();
    alloc(Nl, Nc, BlockSize);

    if (BlockSize > 0)
    {
       Ifloat Filter(Nl,Nc, "Filter");
       Ifloat ImaBlock(BS, BS,"block");

       for (i = 0; i < Nlb; i++)
       for (j = 0; j < Ncb; j++)
       {
           get_block_ima(i,j, Ima, ImaBlock);
           undec_wt_filtering_one_block(ImaBlock, NoiseIma, N_Sigma);
           add_block_ima(i,j, Filter, ImaBlock);
       }
       Ima=Filter;
    }
    else undec_wt_filtering_one_block(Ima, NoiseIma, N_Sigma);
}

/****************************************************************************/

void Ridgelet::radon_transform_one_cst_block()
{
    int Nlr = RADON.radon_nl();
    int Ncr = RADON.radon_nc();
    //cout << "Radon Image Size = " << BS << " ==> " << Nlr << " " << Ncr << endl;
    Ifloat ImaBlock(BS, BS,"block");
    ImaBlock.init(1.);
    RADON.transform(ImaBlock, ImaCstRadonTransOneBlock);
    //cout << "Radon Trans Image Size = " << ImaCstRadonTransOneBlock.nl() << " " << ImaCstRadonTransOneBlock.nc() << endl;
    for (int i=0; i < ImaCstRadonTransOneBlock.nl(); i++)
    for (int j=0; j < ImaCstRadonTransOneBlock.nc(); j++)
    {
       // ImaCstRadonTransOneBlock(i,j) = (int) (ImaCstRadonTransOneBlock(i,j) + 0.5);
       if (ImaCstRadonTransOneBlock(i,j) < 1) ImaCstRadonTransOneBlock(i,j) = 1.;
       else ImaCstRadonTransOneBlock(i,j) = sqrt(ImaCstRadonTransOneBlock(i,j));
       // if (ImaCstRadonTransOneBlock(i,j) == 0) ImaCstRadonTransOneBlock(i,j) = 1.;
    }
    // io_write_ima_float("xx_cst.fits", ImaCstRadonTransOneBlock);
}

/****************************************************************************/

void Ridgelet::thresholding(Ifloat &Transf, float NoiseIma, float N_Sigma)
{
    int i,j,s;
    float SigmaNoise = NoiseIma;

    if (TabNorm.n_elem() < NbrScale-1) set_tab_norm(NbrScale);
    if ((VarianceStab == False) && (RidTrans != RID_FSS))  SigmaNoise *= sqrt((float) BS);

    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale-1; s++)
        {
           int Nls = size_scale_nl(s);
           int Ncs = size_scale_nc(s);
           int Posi = ipos(s);
           int Posj = jpos(s);

           float Nsig = (s == 0) ? N_Sigma + 1: N_Sigma;
           float Level =  SigmaNoise*Nsig*TabNorm(s);
           // float Level =  SigmaNoise*N_Sigma*TabNorm(s);
           for (i=0; i < Nls; i++)
           for (j=0; j < Ncs; j++)
              if ((s < FirstDetectScale) || (ABS(Transf(Posi+i,Posj+j)) < Level))
                                                      Transf(Posi+i,Posj+j) = 0.;
        }
	// cout << "Last scale " << endl;
	if (KillLastScale == True) init_scale(Transf, NbrScale-1);
    }
    else
    {
        int Nls = Transf.nl();
        int Ncs = Transf.nc();
        float Level =  SigmaNoise*N_Sigma;
        for (i=0; i < Nls; i++)
        for (j=0; j < Ncs; j++)
            if (ABS(Transf(i,j)) < Level) Transf(i,j) = 0.;
    }
}

/****************************************************************************/

void Ridgelet::init_scale(Ifloat &Transf, int s)
{
    int i,j;

    if (NbrScale <= s)
    {
         cout << "Error in init_scale: scale number is to high ... " << endl;
         cout << "        NbrScale =  " << NbrScale << " s = " << s << endl;
         exit(-1);
    }
    int Nls = size_scale_nl(s);
    int Ncs = size_scale_nc(s);
    int Depi = ipos(s);
    int Depj = jpos(s);
    for (i=Depi; i < Depi+Nls; i++)
    for (j=Depj; j < Depj+Ncs; j++) Transf(i,j) = 0.;
}

/****************************************************************************/
