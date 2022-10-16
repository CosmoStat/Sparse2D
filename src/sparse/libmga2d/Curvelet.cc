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
**    File:  Curvelet.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Curvelet transform and reconstruction
**    -----------
**
******************************************************************************/

#include "Curvelet.h"
#include "IM_IO.h"
#include "RidNoise.h"
#include "PrimeNumber.h"
#include "IM_Prob.h"
#include "DefFunc.h"

#define NBR_NOISE_2DSCALE 10
#define NBR_NOISE_1DSCALE 20
#define NBR_NOISE_FINITE_1DSCALE 20

#define CUR_DEBUG 0
#define CUR_THRESHOLD_ESTIMATION 0

// static float TabCurvNorm_sigma[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
//     {{0.7430,0.0497,0.0094,0.0039,0.0027},
//     {0.1919,0.0913,0.0252,0.0054,0.0015},
//     {0.0484,0.0763,0.0516,0.0166,0.0038},
//     {0.0108,0.0289,0.0480,0.0345,0.0117},
//     {0.0029,0.0073,0.0197,0.0327,0.0241}};

static float TabCurvNorm[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
    {{0.692084, 0.0488825, 0.00915464, 0.00390845, 0.00247424},
    {0.184628, 0.0853162, 0.0254304, 0.00558103, 0.00150869},
    {0.0489787, 0.0737917, 0.0490959, 0.016418, 0.00425767},
    {0.0118884, 0.0291999, 0.0459632, 0.0330577, 0.0122221},
    {0.00407141, 0.00783881, 0.0200084, 0.0305753, 0.0232205}};

static float TabCurvNormOrtho[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
    {{0.888171, 0.693059, 0.316226, 0.133076, 0.0792442},
    {0.0920643, 0.355331, 0.4478, 0.258705, 0.112331},
    {0.0133684, 0.0687863, 0.290654, 0.41791, 0.244424},
    {0.00712275, 0.017331, 0.0713842, 0.27431, 0.404305},
    {0.00264867, 0.00782153, 0.0244814, 0.0768586, 0.275698}};

static float TabCurvNormFinite[NBR_NOISE_2DSCALE][NBR_NOISE_FINITE_1DSCALE] =
    {{0.907205,   0.848266,  0.750752, 0.706939, 0.693558, 0.769639, 0.757993, 0.876104, 0.943324, 1.34013},
    {0.0575738,   0.219989,  0.308451, 0.285516, 0.278624, 0.297142, 0.32716,  0.402356, 0.487739, 0.691129},
    {0.00477841,  0.0294306, 0.126726, 0.205192, 0.189233, 0.233786, 0.239031, 0.299415, 0.282413, 0.358558},
    {0.00182506, 0.00365026, 0.0238334,0.104979, 0.160795, 0.158492, 0.213155, 0.297539, 0.241235, 0.373624},
    {0.000782209,0.00130493, 0.00312119,0.0230133,0.0992632,0.15913, 0.185009, 0.194789, 0.240767, 0.259313}};

static float TabCurvNormPoisson[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
    {{0.371955 , 0.0407213, 0.0359829, 0.0359829, 0.0359829},
    {0.095282, 0.0462022, 0.0164799, 0.0141048, 0.0141048},
    {0.0259683, 0.0364866, 0.025523, 0.0189549, 0.0189549},
    {0.00621935, 0.0140622, 0.0221193, 0.0245154, 0.0245154},
    {0.00624714, 0.00387116, 0.00893436,  0.027382, 0.027382}};

/****************************************************************************/

const char * StringBlockType (type_curvelet_block type)
{
    switch (type)
    {
        case BL_CONST:
	      return ("Constant block size");
        case BL_UP:
              return ("Block size is doubled at each resolution");
        case BL_UP2:
              return ("Block size is doubled at each two resolutions");
        case BL_DOWN:
	      return ("Block size is divided by 2 at each resolution");
        case BL_DOWN2:
	      return ("Block size is divided by 2 at each two resolutions");
        default:
	      return ("Undefined block type");
    }
}

/****************************************************************************/

float Curvelet::norm_band(int Inds2d, int Inds1d)
{
   int s2d = Inds2d;
   int s1d = Inds1d;
   switch (RidTrans)
   {
      case RID_OWT:
      case RID_UWT:
         if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
         if ((s1d < 0) || (s1d >= NBR_NOISE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
         break;
      case RID_PYR_FFT:
      case RID_PAVE_FFT:
      case RID_FSS:
         if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
         if ((s1d < 0) || (s1d >= NBR_NOISE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
         break;
      case RID_FINITE:
         if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
         if ((s1d < 0) || (s1d >= NBR_NOISE_FINITE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
         break;
     case RID_UNKNOWN:
     default:
         cout << "Error: unknown ridgelet transform (in Curvelet::norm_band) ... " << endl;
         exit(0);
         break;
   }
   return TabCurSigma(s1d,s2d);
}


/****************************************************************************/

void Curvelet::print_info_noise()
{
   for (int b=0; b < nbr_band()-1; b++)
   {
      int s2d,s1d;
      get_scale_number(b, s2d, s1d);
      int BSize = TabBlockSize(s2d);
      cout << "Band " << s2d+1 << " " << s1d+1 << ": BSize BS = " << BSize << " Nl = " << size_nl(s2d,s1d) << " Nc = " << size_nc(s2d,s1d)  << endl;
      cout << "        SigmaNoise = " << sigma_noise(b) << " Nsigma = " << nsigma(b) << endl;
      cout << "        MinLevel = " << min_level(b)  << " MaxLevel = " << max_level(b) << endl;
   }
}

/****************************************************************************/

void Curvelet::set_noise_model_gaussian_in_band(int b, float N_Sigma, float SigmaNoise)
{
    int s2d,s1d;
    get_scale_number(b, s2d,s1d);
    float BSize = TabBlockSize(s2d);
    float BSz = (MSVST == False) ? sqrt(BSize) : 1.;
    float Norm = (VarNorm == True) ? norm_band(s2d,s1d):  norm_band(s2d,s1d) * BSz;
    TabMinDetect(b) = -N_Sigma*Norm*SigmaNoise;
    TabMaxDetect(b) = N_Sigma*Norm*SigmaNoise;
    TabSigmaNoise(b) = Norm*SigmaNoise;
    TabNSigma(b) = N_Sigma;
}

/****************************************************************************/

void Curvelet::set_noise_model_gaussian(fltarray & TabN_Sigma, float SigmaNoise)
{
    int s2d,s1d;
    for (int b=0; b < nbr_band()-1; b++)
    {
       get_scale_number(b, s2d,s1d);
       set_noise_model_gaussian_in_band(b,  TabN_Sigma(b), SigmaNoise);
    }
}

/****************************************************************************/

void Curvelet::set_noise_model_gaussian(float N_Sigma,  float SigmaNoise)
{
    int s2d,s1d;
    for (int b=0; b < nbr_band()-1; b++)
    {
       get_scale_number(b, s2d,s1d);
       float Nsig = (b == 0) ? N_Sigma+1: N_Sigma;
       set_noise_model_gaussian_in_band(b, Nsig,  SigmaNoise);
    }
}

/****************************************************************************/

void Curvelet::set_noise_model_using_mad(Ifloat * &Trans, float N_Sigma)
{
    int s2d,s1d;
    for (int b=0; b < nbr_band()-1; b++)
    {
       Ifloat Band;
       get_scale_number(b, s2d,s1d);
       Ridgelet *Rid = get_ridgelet(s2d);
       Rid->get_scale(Trans[s2d], Band, s1d);
       float NoiseBand = detect_noise_from_mad(Band);
       TabMinDetect(b) = -N_Sigma*NoiseBand;
       TabMaxDetect(b) = N_Sigma*NoiseBand;
       TabSigmaNoise(b) = NoiseBand;
       TabNSigma(b) = N_Sigma;
    }
}


/****************************************************************************/

double Curvelet::prob_noise(float Val, int b)
{
   double P=0;
   if (TabCP != NULL)
   {
       P = TabCP[b].repartition(Val);
       if (Val > 0) P = 1. - P;
   }
   else
   {
      double Sig = TabSigmaNoise(b);
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
   }
   return P;
}

/****************************************************************************/

void Curvelet::get_fdr_detect_level(Ifloat * &Trans, float N_Sigma, Bool ThresholdBand)
{
    int i,j,b,s2d,s1d;
    for (b=0; b < nbr_band()-1; b++)
    {
        get_scale_number(b, s2d, s1d);
        Ridgelet *Rid = get_ridgelet(s2d);
        Ifloat Band;
        Rid->get_scale(Trans[s2d], Band, s1d);
	float Alpha = (1. - erf(N_Sigma / sqrt(2.)))*(s1d+s2d+1);
	if (Alpha > 0.5)  Alpha = 0.5;

       if (Verbose == True)
       {
         cout << "Band " << s2d+1 << " " << s1d+1 << " Nl = " << Band.ny() << " Nc = " << Band.nx() << endl;
	 cout << "        SigmaNoise = " <<  sigma_noise(b) << " Nsigma = " <<  nsigma(b) << endl;
	 cout << "        MinLevel = " <<  min_level(b)  << " MaxLevel = " <<  max_level(b) << endl;
       }

        // get_band(Trans, b, Band);
	dblarray PVal(Band.nc(),Band.nl());
        for (i=0; i < Band.nl(); i++)
        for (j=0; j < Band.nc(); j++) PVal(j,i) = prob_noise(Band(i,j) ,b);
         double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);
        float NSigmaBand = ABS(xerfc(0.5+(1-PDet)/2.));
	for (i=0; i < Band.nl(); i++)
        for (j=0; j < Band.nc(); j++)
	{
  	    if (PVal(j,i) > PDet)
	    {
 	       if ((Band(i,j) >= 0) && (Band(i,j) > TabMaxDetect(b))) TabMaxDetect(b) = Band(i,j);
               if ((Band(i,j) < 0) && (Band(i,j) < TabMinDetect(b))) TabMinDetect(b) = Band(i,j);
	       if (ThresholdBand == True) Band(i,j) = 0.;
 	    }
	}
	TabNSigma(b) = NSigmaBand;
	if (Verbose == True) cout << "      FDR: Level Alpha = " << Alpha << "   ==> Nsigma = " << NSigmaBand << " Min Level = " << TabMinDetect(b) << " Max level = " << TabMaxDetect(b) << endl;
        if (ThresholdBand == True) Rid->put_scale(Trans[s2d], Band, s1d);
  }
}


/****************************************************************************/

void Curvelet::set_noise_level(Ifloat &NoiseData, float N_Sigma, Bool UseMad)
{
    int b,s2d,s1d;
    double LMin, LMax;
    fltarray Band,CUR_Trans;
    Ifloat IBand;
    if (Verbose == True) cout << "Curvelet transform of the noise " << endl;
    transform(NoiseData, CUR_Trans);

    if (UseMad == False) TabCP = new CImaProb [nbr_band()];

    for (b=0; b < nbr_band()-1; b++)
    {
        get_scale_number(b, s2d, s1d);
        float Nsig = (s2d == 0) ? N_Sigma+1: N_Sigma;
        get_band(CUR_Trans, b, Band);
        if (UseMad == False)
        {
          IBand.alloc(Band.buffer(),Band.ny(), Band.nx());
          TabCP[b].set(IBand);
          TabCP[b].find_gthreshold(Nsig, LMin, LMax);
	  TabMinDetect(b) = LMin;
          TabMaxDetect(b) = LMax;
	  // cout << "Band " << b+1 << " LMin = " << LMin << " Lmax = " << LMax << endl;
          // TabCurLevel(b) = MAX(ABS(LMin), LMax);
          // TabCurSigma(s1d,s2d) =  sigma(IBand);
          TabSigmaNoise(b) = MAX(ABS(LMin), LMax) / Nsig;
	  if (VarNorm == False)  TabSigmaNoise(b) /= sqrt((float) TabBlockSize(s2d));
	  TabNSigma(b) = Nsig;
       }
       else
       {
          // cout << "MAD est: " << s1d << " " << s2d << " " << Band.n_elem() << endl;
          float SigmaMad = get_sigma_mad(Band.buffer(),Band.ny()*Band.nx());
          float Nsig = (s2d == 0) ? N_Sigma+1: N_Sigma;
	  TabMinDetect(b) = - Nsig*SigmaMad ;
          TabMaxDetect(b) = Nsig*SigmaMad;
          TabSigmaNoise(b) = SigmaMad;
	  TabNSigma(b) =  Nsig;
       }
    }
}

/****************************************************************************/

static void smooth_bspline (Ifloat & Im_in, Ifloat &Im_out, type_border Type, int Step_trou)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step = 1 << Step_trou;
    int TwoStep = Step << 1;
    double Coeff_h0 = 3. / 8.;
    double Coeff_h1 = 1. / 4.;
    double  Coeff_h2 = 1. / 16.;
    Ifloat Buff(Nl,Nc,"Buff smooth_bspline");

#ifdef _OPENMP
    #pragma omp parallel for private(i,j) shared(Buff,Im_in)
#endif
    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
    {
       Buff(i,j) = (float)( Coeff_h0 * (double) Im_in(i,j)
                 + Coeff_h1 * (  Im_in (i, j-Step, Type)
                               + Im_in (i, j+Step, Type))
                 + Coeff_h2 * (  Im_in (i, j-TwoStep, Type)
                               + Im_in (i, j+TwoStep, Type)) );
    }

#ifdef _OPENMP
    #pragma omp parallel for private(i,j) shared(Im_out,Buff)
#endif
    for (i = 0; i < Nl; i ++)
    for (j = 0; j < Nc; j ++)
       Im_out(i,j) = Coeff_h0 * Buff(i,j)
                 + Coeff_h1 * (  Buff (i-Step, j, Type)
                               + Buff (i+Step, j, Type))
                + Coeff_h2 * (  Buff (i-TwoStep, j, Type)
                               + Buff (i+TwoStep, j, Type));
}

/****************************************************************************/

// static void atrous2d(Ifloat &Image, Ifloat *TabBand, int NbrScale, type_border Bord)
// {
//    int s,i,j;
//    TabBand[0] = Image;
//    for (s = 0; s <  NbrScale-1; s++)
//    {
//         smooth_bspline (TabBand[s], TabBand[s+1],Bord,s);
//         for (i=0;i<Image.nl(); i++)
//         for (j=0;j<Image.nc(); j++) (TabBand[s])(i,j) -= (TabBand[s+1])(i,j);
//    }
// }

/****************************************************************************/

void Curvelet::alloc(int Nl, int Nc, int BS, Bool WeightBeforeTrans)
{
   int s;
   CurNl = 0;
   CurNc = 0;
   extern type_format Format_Imag;
   FormatInputImag = Format_Imag;

   if (UseFSS(RidTrans) == True) OddBlockSize=False;

   if (Verbose == True) cout << "Class allocation ... " << endl;
   tab_block_size_init(BS);
   if (Verbose == True) cout << "Ridgelet allocation ... " << endl;

   if (TabRidgelet != NULL) delete [] TabRidgelet;
   TabRidgelet = new Ridgelet[NbrScale2D-1];
   SizeTabRid = NbrScale2D-1;
   AllocClass = True;

   for (s = 0; s <  NbrScale2D-1; s++)
   {
      if (Verbose == True) cout <<  "INIT RIDGELET Scale " << s+1 << ", Block size = " << TabBlockSize(s) << endl;
      ridgelet_init( TabRidgelet[s]);
      TabRidgelet[s].alloc(Nl, Nc, TabBlockSize(s));
      TabRidgelet[s].WeightBefTrans = WeightBeforeTrans;
      TabRidgelet[s].VarNorm = VarNorm;
   }
   NlIma=Nl;NcIma=Nc;
   int Coef = (rid_class(RidTrans) == RID_CL_ORTHO) ? 1: 2;
   // int CoefOver = (BlockOverlap == False) ? 1: 2;
   if (TabBlockSize(0) > 0) // {CurNl=CoefOver*2*Nl;CurNc=CoefOver*Coef*Nc;}
   {
     for (s=0; s < NbrScale2D-1; s++)
     {
       if (CurNl < TabRidgelet[s].rid_nl()) CurNl = TabRidgelet[s].rid_nl();
       if (CurNc < TabRidgelet[s].rid_nc()) CurNc = TabRidgelet[s].rid_nc();
     }
   }
   else  {CurNl=2*Nl;CurNc=Coef*Nc;}

   if (Verbose == True)
   {
       cout << "CURVELET size: Nl = " << CurNl  << " Nc = " << CurNc << " NbrScale2D = " <<  NbrScale2D << endl;
   }

   TabCurSigma.alloc(NBR_NOISE_1DSCALE,NBR_NOISE_2DSCALE);
   TabMinDetect.alloc(nbr_band());
   TabMaxDetect.alloc(nbr_band());
   TabSigmaNoise.alloc(nbr_band());
   TabNSigma.alloc(nbr_band());
   for (int b=0; b < nbr_band()-1; b++)
   {
      int s2d, s1d;
       get_scale_number(b, s2d, s1d);
       switch (RidTrans)
       {
         case RID_OWT:
           if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
           if ((s1d < 0) || (s1d >= NBR_NOISE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
           TabCurSigma(s1d,s2d) = (MSVST == False) ? TabCurvNormOrtho[s2d][s1d] : TabCurvNormPoisson[s2d][s1d];
	   break;
         case RID_PYR_FFT:
         case RID_PAVE_FFT:
         case RID_UWT:
         case RID_FSS:
           if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
           if ((s1d < 0) || (s1d >= NBR_NOISE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
           if (MSVST == False) TabCurSigma(s1d,s2d) = TabCurvNorm[s2d][s1d];
	   else TabCurSigma(s1d,s2d) = TabCurvNormPoisson[s2d][s1d];
	   break;
         case RID_FINITE:
           if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
           if ((s1d < 0) || (s1d >= NBR_NOISE_FINITE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
           TabCurSigma(s1d,s2d) = TabCurvNormFinite[s2d][s1d];
	   break;
         case RID_UNKNOWN:
	 default:
         cout << "Error: unknown ridgelet transform (in Curvelet::norm_band) ... " << endl;
         exit(0);
         break;
      }
   }
}

/****************************************************************************/

void Curvelet::set_param(int Nl, int Nc)
// Set the curvelet transform size from the input image size
{
   int s;
   CurNl = 0;
   CurNc = 0;

   if (TabRidgelet == NULL)
   {
       TabRidgelet = new Ridgelet[NbrScale2D-1];
       for (s = 0; s <  NbrScale2D-1; s++)
       {
          ridgelet_init( TabRidgelet[s]);
          TabRidgelet[s].alloc(Nl, Nc, TabBlockSize(s));
         cout << endl << "PARAM Scale " << s+1 << "NlRid = " << TabRidgelet[s].rid_nl() <<  " NcRid = " << TabRidgelet[s].rid_nc() << endl;
       }
       TabRidgelet[0].WPTrans = WPTrans;
   }
   NlIma=Nl;NcIma=Nc;
   int Coef = (rid_class(RidTrans) == RID_CL_ORTHO) ? 1: 2;
   // int CoefOver = (BlockOverlap == False) ? 1: 2;
   if (TabBlockSize(0) > 0) // {CurNl=CoefOver*2*Nl;CurNc=CoefOver*Coef*Nc;}
   {
     for (s=0; s < NbrScale2D-1; s++)
     {
       if (CurNl < TabRidgelet[s].rid_nl()) CurNl = TabRidgelet[s].rid_nl();
       if (CurNc < TabRidgelet[s].rid_nc()) CurNc = TabRidgelet[s].rid_nc();
     }
   }
   else  {CurNl=2*Nl;CurNc=Coef*Nc;}

// cout << "CoefOver = " << CoefOver;
// cout << "  Coef= " <<  Coef;
// cout << "  Calc. CurNl= " <<  CurNl;
// cout << "  Calc. CurNc = " <<  CurNc << endl;
}

/****************************************************************************/

void Curvelet::set_param_from_trans(int Nlc, int Ncc)
// Set the image size from the the curvelet transform image size
{
   CurNl = Nlc;
   CurNc = Ncc;
   int Coef = (rid_class(RidTrans) == RID_CL_ORTHO) ? 1: 2;
   int CoefOver = (BlockOverlap == False) ? 1: 2;
   if ( TabBlockSize(0) > 0) { NlIma=Nlc/CoefOver/2; NcIma=Ncc/CoefOver/Coef;}
   else  { NlIma=Nlc/2; NcIma=Ncc/Coef;}
}

/****************************************************************************/

void Curvelet::reset(SubBand1D *SB1D)
{
   extern type_format Format_Imag;
   FormatInputImag = Format_Imag;
   WPTrans= False;
   AllocClass = False;
   SizeTabRid= 0;
   StatInfo = False;
   BlockOverlap = True;
   Ptr_SB1D = SB1D;
   Border=I_MIRROR;
   GetAutoNbScale= True;
   Verbose=False;
   NbrScale2D= DEF_CUR_NBR_SCALE_2D;
   NbrScale1D=-1;
   TabBlockSize.reform(MAX_SCALE);
   TabRidgelet = NULL;
   TabCP = NULL;
   ColTrans= False;
   RidTrans=RID_PYR_FFT;
   TypeBlock = DEF_TYPE_CUR_BLOCK;
   UserBlockSize=DEF_CUR_BLOCK_SIZE;
   OddBlockSize=True;
   AngleNormalization=False;
   tab_block_size_init(DEF_CUR_BLOCK_SIZE);
   VarNorm=False;
   MSVST=False;
}

/****************************************************************************/

void Curvelet::tab_block_size_init(int BS)
{
    CPRIME_NUMBER CPN;
    // Bool Odd = (BS % 2 == 0) ? False: True;
    int s, B = BS;
    // cout << " Type block = " << StringBlockType(TypeBlock) << endl;

    if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
    else B = (int) next_power_of_2((int) BS);

    switch (TypeBlock)
    {
       case BL_CONST:
           for (s=0; s < MAX_SCALE; s++) TabBlockSize(s) = B;
           break;
       case BL_UP:
           for (s=0; s < MAX_SCALE; s++)
           {
              TabBlockSize(s) = B;
              B *= 2;
 	      if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
	   }
	   break;
       case BL_UP2:
           for (s=0; s < MAX_SCALE; s++)
           {
              TabBlockSize(s) = B;
	      if (s % 2 == 0) B *= 2;
 	      // if ((B % 2 == 0) && (Odd == True)) B++;
	      if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
           }
	   break;
       case BL_DOWN:
           for (s=0; s < MAX_SCALE; s++)
           {
              TabBlockSize(s) = B;
              B /= 2;
 	      if (B < 8) B = 8;
	      // if ((B % 2 == 0) && (Odd == True)) B++;
	      if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
           }
	   break;
       case BL_DOWN2:
           for (s=0; s < MAX_SCALE; s++)
           {
              TabBlockSize(s) = B;
              if (s % 2 == 0)
	      {
	         B /= 2;
	         if (B < 8) B = 8;
	      }
	      // if ((B % 2 == 0) && (Odd == True)) B++;
	      if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
           }
	   break;
     }
    if ((RidTrans != RID_FINITE) && (OddBlockSize == True))
             for (s=0; s < MAX_SCALE; s++) TabBlockSize(s) += 1;
}

/****************************************************************************/

void Curvelet::ridgelet_init(Ridgelet &Rid)
{
    // Initialize the ridgelet transform
    if (Ptr_SB1D != NULL) Rid.reset(Ptr_SB1D);
    Rid.BlockOverlap = BlockOverlap;
    Rid.ColTrans = ColTrans ;
    Rid.StatInfo = StatInfo;
    Rid.GetAutoNbScale = GetAutoNbScale;
    Rid.RidTrans = RidTrans;
    Rid.Verbose = Verbose;
    Rid.WPTrans = WPTrans;
    if (GetAutoNbScale == False) Rid.NbrScale = NbrScale1D;
}

/****************************************************************************/

void Curvelet::transform(Ifloat &Image, Ifloat * & TabBand)
{
    int s;
    int Nl = Image.nl();
    int Nc = Image.nc();
    if (AllocClass == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    if ((Nl != NlIma) || (Nc != NcIma))
    {
       cout << "Error: the curvelet class is not initialized for this given image size ... " << endl;
       cout << "       Image size = " << Nl << " " << Nc << endl;
       cout << "       Expected size = " << NlIma << " " << NcIma << endl;
       exit(-1);
    }
    if (SizeTabRid != NbrScale2D-1)
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }
    Ifloat BandCoef(Nl,Nc," BandCoef");
    Ifloat BandSmooth(Nl,Nc," BandSmooth");
    Ifloat BandMSVST;

    BandCoef = Image;
    if ((MSVST == True) && (Verbose == True)) cout << " Multiscale Variance Stabilization ... " << endl;

    if (MSVST == True)
    {
        BandMSVST = Image;
        msvst_transform (BandMSVST, BandMSVST, 1);
    }
    if (Verbose == True)
     cout << endl << "Curvelet Image Transformation " <<  endl;
    if (TabBand == NULL) TabBand = new Ifloat [NbrScale2D];

    for (s = 0; s <  NbrScale2D-1; s++)
    {
        if (Verbose == True)
        {
           cout << endl;
           cout << " BAND " << s+1 << " Block size = " << TabBlockSize(s) << endl;
        }
        smooth_bspline(BandCoef, BandSmooth, Border, s);
	if (MSVST == True)
	{
	    BandCoef = BandMSVST;
	    msvst_transform (BandSmooth, BandMSVST, s+2);
	    BandCoef -= BandMSVST;
        }
	else BandCoef -= BandSmooth;
	// cout << "Scale " << s+1 << " sigma band = " << BandCoef.sigma() << " " <<  sqrt(msvst_var (2, s+1)) <<  endl;

        TabRidgelet[s].transform(BandCoef, TabBand[s]);
        if (AngleNormalization == True)
	        TabRidgelet[s].mad_normalize_per_block(TabBand[s]);
		// TabRidgelet[s].mad_normalize_per_angle(TabBand[s]);

	if (Verbose == True)
           cout << "Number of scale for the 1D WT = " << TabRidgelet[s].NbrScale << endl;

        BandCoef = BandSmooth;
    }
    s = NbrScale2D-1;
    TabBand[s].resize(Image.nl(), Image.nc());
    if (MSVST == True) TabBand[s] = BandMSVST;
    else TabBand[s] = BandSmooth;

//     Ifloat Rec(NlIma, NcIma);
//     recons(TabBand, Rec);
//     Image -= Rec;
//     INFO_X(Image, "Residu");
//     io_write_ima_float("xx_rec", Rec);
//     io_write_ima_float("xx_resi", Image);
}

/****************************************************************************/

void Curvelet::recons(Ifloat * TabBand, Ifloat &Image)
{
    int s,i,j;
    if (AllocClass  == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    int Nl = NlIma;
    int Nc = NcIma;
    Ifloat BandCoef(Nl,Nc,"aux");
    Ifloat BandTrans(CurNl,CurNc,"aux");

    if (Verbose == True)
     cout << "Reconstructed image size: Nl =  " << Nl  << " Nc = " << Nc << endl;

    if ((Image.nl() != Nl) || (Image.nc() != Nc)) Image.resize(Nl,Nc);
    s = NbrScale2D-1;
    Image = TabBand[s];

    for (s = NbrScale2D-2; s >= 0; s--)
    {
        int Nlb = TabRidgelet[s].rid_nl();
	int Ncb = TabRidgelet[s].rid_nc();
        if (Verbose == True) cout << endl << "REC Band " << s+1 << " Block size = " << TabBlockSize(s) << " Nlb = " << Nlb << " Ncb = " << Ncb << endl;
        TabRidgelet[s].recons(TabBand[s], BandCoef);
        if (Verbose == True) cout << "  Number of scale for the 1D WT = " << TabRidgelet[s].NbrScale << endl;
        for (i=0;i< Image.nl(); i++)
        for (j=0;j< Image.nc(); j++) Image(i,j) += BandCoef(i,j);
    }
    if ((MSVST == True) && (Verbose == True)) cout << " Inverse Multiscale Variance Stabilization ... " << endl;
    if (MSVST == True) msvst_inv_transform(Image, Image,1);

}

/****************************************************************************/

void Curvelet::transform(Ifloat &Image, fltarray & CurTrans)
{
    int s,i,j;
    int Nl = Image.nl();
    int Nc = Image.nc();

    if (AllocClass == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    if ((Nl != NlIma) || (Nc != NcIma))
    {
       cout << "Error: the curvelet class is not initialized for this given image size ... " << endl;
       cout << "       Image size = " << Nl << " " << Nc << endl;
       cout << "       Expected size = " << NlIma << " " << NcIma << endl;
       exit(-1);
    }
    if (SizeTabRid != NbrScale2D-1)
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }
    Ifloat BandCoef(Nl,Nc," BandCoef");
    Ifloat BandSmooth(Nl,Nc," BandSmooth");

    Ifloat BandTrans(CurNl,CurNc,"BandTrans ");
    if ((CurTrans.nx() != CurNc) || (CurTrans.ny() != CurNl)
        || (CurTrans.nz() != NbrScale2D))
    {
       if (CurTrans.n_elem() > 0) CurTrans.free();
       CurTrans.alloc(CurNc, CurNl, NbrScale2D);
    }
    BandCoef = Image;
    if (Verbose == True)
     cout << endl << "Curvelet Image Transformation " <<  endl;

    Ifloat BandMSVST;

    BandCoef = Image;
    if ((MSVST == True) && (Verbose == True)) cout << " Multiscale Variance Stabilization ... " << endl;

    if (MSVST == True)  msvst_transform (Image, BandMSVST, 1);

    for (s = 0; s <  NbrScale2D-1; s++)
    {
        if (Verbose == True)
        {
           cout << endl;
           cout << " BAND " << s+1 << " Block size = " << TabBlockSize(s) << endl;
        }
        smooth_bspline (BandCoef, BandSmooth, Border, s);
	if (MSVST == True)
	{
	    BandCoef = BandMSVST;
	    msvst_transform (BandSmooth, BandMSVST, s+2);
	    BandCoef -= BandMSVST;
        }
	else BandCoef-= BandSmooth;

        // transform the band
  	if (s != 0) BandTrans.init();
        TabRidgelet[s].transform(BandCoef, BandTrans);
	if ((Verbose == True) && (AngleNormalization == True)) cout << "  mad_normalize_per_block ... " << endl;
	if (AngleNormalization == True)
	   // TabRidgelet[s].mad_normalize_per_angle(BandTrans);
           TabRidgelet[s].mad_normalize_per_block(BandTrans);

	if (Verbose == True)
           cout << "Number of scale for the 1D WT = " << TabRidgelet[s].NbrScale << endl;
     // cout << " Ridgelet OK " << BandTrans.nl() << " " << BandTrans.nc() << endl;

     // INFO_X(BandTrans, "trans band");

        for (i=0;i<Nl; i++)
        for (j=0;j<Nc; j++)  BandCoef(i,j) = BandSmooth(i,j);

        int Nlb = BandTrans.nl();
        int Ncb = BandTrans.nc();
        for (i=0;i< Nlb; i++)
        for (j=0;j< Ncb; j++) CurTrans(j,i,s) = BandTrans(i,j);


#if CUR_THRESHOLD_ESTIMATION
      if (MSVST == True)
      {
           CImaProb CP;
           double LMin, LMax;
	   Ifloat IBand;
	   for (int i = 0; i < TabRidgelet[s].NbrScale; i++)
	   {
	      fltarray Band;
	      get_band(CurTrans,  s, i, Band);
	      IBand.alloc(Band.buffer(),Band.ny(), Band.nx());
              CP.set(IBand);
	      CP.find_gthreshold(3., LMin, LMax); cout << "    " << i+1 <<  " Min = " <<LMin / 3. << " Max = " << LMax / 3. << " Sigma band = " << Band.sigma() << endl;
	   }
	   cout << endl;
       }
       else
        // Threshold estimation: we use to determine the normalization table
        {
            cout << endl << "Band " << s+1 << endl;
            fltarray TabMin(TabRidgelet[s].NbrScale), TabMax(TabRidgelet[s].NbrScale);
            fltarray Tab_Nsig(TabRidgelet[s].NbrScale);
            Tab_Nsig.init(3.);
            StatRidNoise SRN;
            // SRN.Verbose = True;
            SRN.init(BandTrans, TabRidgelet[s].NbrScale, TabRidgelet[s]);
            SRN.find_nsigma_threshold(TabRidgelet[s].NbrScale, Tab_Nsig, TabMin, TabMax);
            for (i=0; i < TabRidgelet[s].NbrScale-1; i++)
            {
               float Coef =  sqrt((float) Nl) * Tab_Nsig(i);
               TabMin(i) = ABS(TabMin(i))/Coef;
               TabMax(i) /= Coef;
               float ErrMax = ABS((TabCurvNorm[s][i] - TabMax(i))) / TabMax(i) * 100.;
               float ErrMin = ABS((TabCurvNorm[s][i] - TabMin(i))) / TabMin(i) * 100.;
               cout << MAX(TabMin(i),TabMax(i)) << ", ";
               // cout << i+1<< ": Snsig = " << Tab_Nsig(i) << " Min = " << TabMin(i) << " Max = " << TabMax(i) << endl;
               // cout << "        " <<  " ErrMin = " <<  ErrMin << " ErrMax  = " <<  ErrMax << endl;
            }
        }
#endif
//         Bool NormalizeCoeff = True;
//         if (NormalizeCoeff == True)
//         {
// 	   fltarray Tab;
// 	   int Ns1D = TabRidgelet[s].NbrScale;
// 	   int bs = TabBlockSize(s);
//       cout << " SCALE    " << s+1 << " NbrScale1D = " <<  Ns1D << " BlockSize = " << TabBlockSize(s) << endl;
//            get_norm_coeff(TabBlockSize(s), s, Ns1D, Tab);
//  	   for (int s1d = 0; s1d < TabRidgelet[s].NbrScale; s1d++)
// 	   {
// 	      fltarray Band,Band1,Band2;
// 	      get_band(CurTrans,  s, s1d, Band);
// 	      fits_write_fltarr("xx_tab.fits", Tab);
// 	      Band1.alloc(Band.nx());
// 	      Band2.alloc(Band.nx());
// 	      for (i=0; i < Band.nx(); i++) Band1(i) = Band(i,0);
// 	      for (i=0; i < Band.nx(); i++) Band2(i) = Band(i,bs/2);
//               cout << "    " << s1d+1 <<  " Nl = " << Band.ny() << " Nc = " << Band.nx() << endl;
// 	      cout << "       Coef1 = " <<  Tab(0,s1d)*sqrt((float) TabBlockSize(s)) << " Coef2 = " <<  Tab(bs/2,s1d)*sqrt((float) TabBlockSize(s));
// 	      cout << "       Sigma band = " << Band.sigma();
// 	      cout << " Sigma band1 = " << Band1.sigma() << " Sigma band2 = " << Band2.sigma() << endl;
//            }
//            for (int s1d = 0; s1d < TabRidgelet[s].NbrScale; s1d++)
//            {
//                int NFirst = TabRidgelet[s].rid_pos(s1d);
//                int NLast = TabRidgelet[s].rid_pos(s1d) + TabRidgelet[s].rid_size(s1d);
//                for (i=0; i < BandTrans.nl(); i++)
// 	       for (j=NFirst; j < NLast; j++)
//                            CurTrans(j,i,s) /= Tab(i%TabBlockSize(s),s1d);
//  	   }
//        }

    }
    if (MSVST == True)
    {
#ifdef _OPENMP
    #pragma omp parallel for private(i,j) shared(CurTrans,BandMSVST)
#endif
       for (i=0;i<Nl; i++)
       for (j=0;j<Nc; j++)  CurTrans(j,i,NbrScale2D-1) = BandMSVST(i,j);
    }
    else
    {
#ifdef _OPENMP
    #pragma omp parallel for private(i,j) shared(CurTrans,BandSmooth)
#endif
       for (i=0;i<Nl; i++)
       for (j=0;j<Nc; j++)  CurTrans(j,i,NbrScale2D-1) = BandSmooth(i,j);
    }

}

/****************************************************************************/

void Curvelet::recons(fltarray & CurTrans, Ifloat &Image, Bool Init)
{
    int s,i,j;
    if (AllocClass  == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    if ((CurTrans.nx() != CurNc) || (CurTrans.ny() != CurNl))
    {
       cout << "Error: the curvelet class is not initialized for this given transform size ... " << endl;
       cout << "       Transform size = " << CurTrans.ny() << " " << CurTrans.nx() << endl;
       cout << "       Expected size = " << CurNl << " " << CurNc << endl;
       exit(-1);
    }
    if ((CurTrans.nz() != SizeTabRid+1) || (NbrScale2D != SizeTabRid+1))
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }

    int Nl = NlIma;
    int Nc = NcIma;
    Ifloat BandCoef(Nl,Nc,"aux");
    Ifloat BandTrans(CurNl,CurNc,"aux");

    if (Verbose == True)
     cout << "Reconstructed image size: Nl =  " << Nl  << " Nc = " << Nc << endl;

    if ((Image.nl() != Nl) || (Image.nc() != Nc)) Image.resize(Nl,Nc);

    for (i=0;i< Nl; i++)
    for (j=0;j< Nc; j++) Image(i,j) = CurTrans(j,i,NbrScale2D-1);

    for (s = NbrScale2D-2; s >= 0; s--)
    {
        int Nlb = TabRidgelet[s].rid_nl();
	int Ncb = TabRidgelet[s].rid_nc();
        if (Verbose == True) cout << endl << "REC Band " << s+1 << " Block size = " << TabBlockSize(s) << " Nlb = " << Nlb << " Ncb = " << Ncb << endl;
        BandTrans.resize(Nlb, Ncb);
        for (i=0;i< Nlb; i++)
        for (j=0;j< Ncb; j++) BandTrans(i,j) = CurTrans(j,i,s);

        TabRidgelet[s].recons(BandTrans, BandCoef);
        if (Verbose == True) cout << "  Number of scale for the 1D WT = " << TabRidgelet[s].NbrScale << endl;
   // cout << " Ridgelet OK " <<  BandCoef.nl() << " " <<  BandCoef.nc() << endl;

        for (i=0;i< Image.nl(); i++)
        for (j=0;j< Image.nc(); j++) Image(i,j) += BandCoef(i,j);
    }
    if ((MSVST == True) && (Verbose == True)) cout << " Inverse Multiscale Variance Stabilization ... " << endl;

    if (MSVST == True) msvst_inv_transform(Image, Image,1);
   //  cout << " RECONS OK " << endl;
}
/****************************************************************************/

void Curvelet::get_norm_coeff(int N, int ScaleNumber, int Ns1D, fltarray &Tab)
{
   cout << "get_norm_coeff " << N << endl;
   int s,i,j;
   float Sig=1.;
   unsigned int InitRnd = 10;

   Ifloat ImSimu(N,N,"ImSimu");
   im_noise_gaussian (ImSimu, Sig, InitRnd);
   Ifloat Result,ImaScale;
   type_transform Transform = TO_PAVE_BSPLINE;
   MultiResol MR_Data;
   MR_Data.alloc (N,N, ScaleNumber+2, Transform);
   MR_Data.transform(ImSimu);
   cout << "WT: get_norm_coeff " << endl;
   ImSimu = MR_Data.band(ScaleNumber);
   Ridgelet RidC(*Ptr_SB1D);
   RidC.Verbose = False;
   RidC.BlockOverlap = False;
   RidC.RidTrans = RidTrans;
   RidC.NbrScale = Ns1D;
   RidC.ColTrans = ColTrans;
   RidC.alloc(N,N,N);

   cout << "WT: get_norm_coeff " << RidC.NbrScale << endl;
   RidC.transform(ImSimu,Result,N);
   cout << "Trans size " << Result.nl() << " " << Result.nc() << " Bs = " << RidC.NbrScale << endl;
   cout << "WT: Ns1D " << RidC.NbrScale << endl;
   Tab.reform(Result.nl(),  RidC.NbrScale);

   for (s = 0; s < RidC.NbrScale; s++)
   {
       // cout << "    " << s+1 << " RidSize = " << RidC.rid_size(s) << endl;

       int NFirst = RidC.rid_pos(s);
       // int NLast = RidC.rid_pos(s) + RidC.rid_size(s);
       for (i=0; i < Result.nl(); i++)
       {
           float Sum=0.;
           for (j=0; j < RidC.rid_size(s); j++)
	   {
	      Sum += Result(i,NFirst+j)*Result(i,NFirst+j);
	   }
	   Sum /= (float) RidC.rid_size(s);
	   Tab(i,s) = sqrt(Sum) / sqrt((double) N);
       }
   }
   // cout << "Normalization coeff. for scale " <<ScaleNumber +1 << endl;
   Ifloat Band(Result.nl(), Result.nc(), "scale");
   s = 0;
   // for (s = 0; s < RidC.NbrScale; s++)
   {
      RidC.get_scale(Result,Band,s);
      cout << "SCALE " << ScaleNumber+1 << "    band " << s+1 << " Sigma = " << sqrt(energy (Band)) << endl;
      // for (i=0; i < Result.nl(); i++) printf(" %6.4f",  Tab(i,s));
   }
   cout << endl <<"End Normalization " << endl;
  // exit(-1);
}

/****************************************************************************/

void Curvelet::filtering(Ifloat &Image, Ifloat  & Filter, float NoiseIma, float N_Sigma)
{
    fltarray TN;
    int s,i,j;
    Ridgelet Rid;
    ridgelet_init(Rid);
    int Nl = Image.nl();
    int Nc = Image.nc();
    Ifloat BandCoef(Nl,Nc," BandCoef");
    Ifloat BandSmooth(Nl,Nc," BandSmooth");
    set_param(Nl,Nc);
    Ifloat BandTrans(CurNl,CurNc,"BandTrans ");
    //static double TabNorm[MAX_SCALE] =
    //{0.7538,  0.2677, 0.1659, 0.1141, 0.0794, 0.01, 0.01, 0., 0., 0.};

    BandCoef = Image;
    if (Verbose == True)
    {
        cout << "Transform size: Nl =  " << CurNl  << " Nc = " << CurNc  << " Nb = " << NbrScale2D << endl;
        for (s = 0; s < NbrScale2D-1; s++)
	   cout << "   Band " << s+1 << "  BlockSize = " << TabBlockSize(s) << endl;
    }

    Filter.init();
    TN.alloc(NBR_NOISE_1DSCALE);

    for (s = 0; s < NbrScale2D-1; s++)
    {
//         fltarray Tab;
// 	int Ns1D = Rid.NbrScale;
//         get_norm_coeff(TabBlockSize(s), s, Ns1D, Tab);

        if (Verbose == True) cout << "Band " << s+1 << " Block size = " << TabBlockSize(s) << endl;
        smooth_bspline(BandCoef, BandSmooth, Border, s);

        for (i=0;i<Image.nl(); i++)
        for (j=0;j<Image.nc(); j++)  BandCoef(i,j) -= BandSmooth(i,j);

        // transform the band
        if (s < NBR_NOISE_2DSCALE)
          for (i=0; i < NBR_NOISE_1DSCALE; i++) TN(i)= norm_band(s,i);
        else for (i=0; i < NBR_NOISE_1DSCALE; i++) TN(i)=0.;

        if (rid_class(RidTrans) != RID_CL_PAVE)
        {
           Rid.filtering(BandCoef, NoiseIma, N_Sigma, TabBlockSize(s), &TN);
        }
        else
        {
           Rid.undec_wt_filtering(BandCoef,  NoiseIma,  N_Sigma, TabBlockSize(s));
        }

#if CUR_DEBUG
        if (Verbose == True)
        {
           cout << "Number of scale for the 1D WT = " << Rid. NbrScale << endl;
           char NameBand[256];
           sprintf(NameBand, "cur_band_%d.fits", s+1);
           fits_write_fltarr(NameBand, BandCoef);
        }
#endif
        for (i=0;i<Image.nl(); i++)
        for (j=0;j<Image.nc(); j++)
        {
           Filter(i,j) += BandCoef(i,j);
           BandCoef(i,j) = BandSmooth(i,j);
        }
    }
    for (i=0;i<Image.nl(); i++)
    for (j=0;j<Image.nc(); j++)  Filter(i,j) += BandSmooth(i,j);
}

/****************************************************************************/

int Curvelet::nbr_band()
{
   int Nb = 1;
   for (int b=0; b < NbrScale2D-1; b++)  Nb+=TabRidgelet[b].NbrScale;
   return Nb;
}

/****************************************************************************/

void Curvelet::get_scale_number(int NumBand, int & s2d, int & s1d)
{
   int b=0;
   int NbrBand = nbr_band();
   s1d = s2d = 0;
   if ((NumBand < 0) || (NumBand >= NbrBand))
   {
       cerr << "Error: bad number of band in get_scale_number: " << NumBand << endl;
       cerr << "       0 <= NumBand < " << NbrBand << endl;
       exit(-1);
   }
   // cout << endl << "get_scale_number" << TabRidgelet[s2d].NbrScale << endl;
   while (b+TabRidgelet[s2d].NbrScale <= NumBand)
   {
      b += TabRidgelet[s2d].NbrScale;
      s2d ++;
      // cout << " b = " << b << endl;
   }
   if (s2d == NbrBand-1) s1d = 0;
   else s1d = NumBand - b;
   // if (s1d < 0)
   // {
   // cout << "ERROR: s1d = " << s1d << " s2d = " << s2d <<" NumBand = " << NumBand << " b = " << b << endl;
   // exit(-1);
   // }
   // cout << "NumBand = " << NumBand << " ==> s2d = " <<  s2d << ", s1d = " <<  s1d << endl;
}

/****************************************************************************/

void Curvelet::get_band(fltarray &Trans, int s2d, int s1d, fltarray &Band)
{
   int i,j,Posi = ipos(s2d, s1d);
   int Posj = jpos(s2d, s1d);
   int Nl =  size_nl(s2d, s1d);
   int Nc =  size_nc(s2d, s1d);
   Band.resize(Nc,Nl);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Band(j,i) = Trans(j+Posj, i+Posi, s2d);
}

/****************************************************************************/

void Curvelet::get_band_no_border(fltarray &Trans, int s2d, int s1d, fltarray &Band)
{
   int NlOneBlock = TabRidgelet[s2d].rid_one_block_nl();
   int NcOneBlock = TabRidgelet[s2d].rid_size(s1d);
   int i,j,Posi = ipos(s2d, s1d);
   int Posj = jpos(s2d, s1d);
   int Nl =  size_nl(s2d, s1d);
   int Nc =  size_nc(s2d, s1d);
   int Nls = Nl - 2*NlOneBlock;
   int Ncs = Nc - 2*NcOneBlock;
   Band.resize(Nls,Ncs);
   int Ind=0;
   for (i = NlOneBlock; i < Nl-NlOneBlock; i++)
   for (j = NcOneBlock; j < Nc-NcOneBlock; j++) Band(Ind++) = Trans(j+Posj, i+Posi, s2d);
   if (Ind != Nls*Ncs) cout << "BIG PB " << " Ind = " << Ind << " Nls*Ncs = " << Nls*Ncs << endl;
}

/****************************************************************************/

void Curvelet::put_band(fltarray &Trans, int s2d, int s1d, fltarray &Band)
{
    int i,j,Posi =  ipos(s2d, s1d);
    int Posj =  jpos(s2d, s1d);
    int Nl =  size_nl(s2d, s1d);
    int Nc =  size_nc(s2d, s1d);

    if ((Nl != Band.ny()) || (Nc != Band.nx()))
    {
       cerr << "Error: bad image size: " << Band.ny() << " " <<  Band.nx() << endl;
       cerr << "       Expected size: " << Nl << " " << Nc << endl;
       exit(-1);
    }
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)  Trans(j+Posj, i+Posi, s2d) = Band(j,i);
}

/****************************************************************************/

void Curvelet::get_band(fltarray &TabBand, int NumBand, fltarray &Band)
{
   int s2d, s1d;
   get_scale_number(NumBand, s2d, s1d);
   get_band(TabBand, s2d, s1d, Band);
}

/****************************************************************************/

void Curvelet::put_band(fltarray &TabBand, int NumBand, fltarray &Band)
{
   int s2d, s1d;
   get_scale_number(NumBand, s2d, s1d);
   put_band(TabBand, s2d, s1d, Band);
}

/****************************************************************************/
//
//                         Curvelet IO
//
/****************************************************************************/

/*************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 4) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'u')
		|| (File_Name_In[L-3] != 'c')
                || (File_Name_In[L-4] != '.'))
    {
        strcat (File_Name_Out, ".cur");
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

void Curvelet::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nl", (long) NlIma,(char*)"NlIma",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Nc",(long) NcIma,(char*)"NcIma",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"NlTrans",(long) CurNl ,(char*)"CurNl",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"NcTrans",(long) CurNc ,(char*)"CurNc",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"BSize",(long) UserBlockSize,(char*)"BlockImaSize",&status))
     PrintError( status );
  if ( ffpkyj(fptr, (char*)"TBlock", (long) TypeBlock,(char*)"Type of block",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"Overlap",(long) ((BlockOverlap == True) ? 1: 0),(char*)"Block overlap",&status))
     PrintError( status );
  if ( ffpkyj(fptr,(char*)"NScale2D",(long)NbrScale2D,(char*)"Number of scales (2D)",&status))
     PrintError( status );
  if (GetAutoNbScale == True)
  {
     if ( ffpkyj(fptr,(char*)"NScale1D",(long) 0,(char*)"Number of scales (1D)",&status))
      PrintError( status );
  }
  else
  {
     if ( ffpkyj(fptr,(char*)"NScale1D",(long)NbrScale1D,(char*)"Number of scales (1D)",&status))
      PrintError( status );
  }
  if ( ffpkyj(fptr, (char*)"Type_Tra", (long) RidTrans ,
                          (char*)StringRidTransform(RidTrans), &status))
    PrintError( status );
  if ( ffpkyj(fptr, (char*)"ColWT", (long) ColTrans, (char*)"ColWT", &status))
    PrintError( status );
  if ( ffpkyj(fptr, (char*)"FormatIn",(long)FormatInputImag,(char*)"Format", &status))
      PrintError( status );
}


/****************************************************************************/

void Curvelet::write (char *Name, fltarray &TabBand)
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
 naxis = 3;
 naxes[0] = TabBand.nx();
 naxes[1] = TabBand.ny();
 naxes[2] = TabBand.nz();

 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // write the header of the multiresolution file
  mr_io_fill_header(fptr);

  nelements = naxes[0] * naxes[1] * naxes[2];
  if ( ffppre(fptr, group, firstpixel, nelements, TabBand.buffer(), &status) )
              PrintError( status );

 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );
// cout << " end of write fits " << endl;

}
/****************************************************************************/

void Curvelet::read (char *Name, fltarray &TabBand)
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
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];

     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;

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

     nelements = naxes[0] * naxes[1] * naxes[2];
     // cout << " begin to read " << endl;
     TabBand.alloc(naxes[0], naxes[1], naxes[2]);

    if (ffgkyj(fptr,(char*)"Nl", &mon_long, comment, &status)) PrintError( status );
    int Nli = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Nc", &mon_long, comment, &status)) PrintError( status );
    int Nci = (int) mon_long;
    if (ffgkyj(fptr,(char*)"BSize", &mon_long, comment, &status)) PrintError( status );
    UserBlockSize  = (int) mon_long;
    if (ffgkyj(fptr,(char*)"NScale2D", &mon_long, comment, &status)) PrintError( status );
    NbrScale2D = (int) mon_long;
    if (ffgkyj(fptr,(char*)"NScale1D", &mon_long, comment, &status)) PrintError( status );
    if (mon_long > 0)
    {
       NbrScale1D = mon_long;
       GetAutoNbScale = False;
    }
    // NbrScale1D = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
    RidTrans = (type_ridgelet_WTtrans) mon_long;
    if (ffgkyj(fptr,(char*)"TBlock", &mon_long, comment, &status)) PrintError( status );
    TypeBlock = (type_curvelet_block) mon_long;
    if ( ffgkyj(fptr,(char*)"Overlap",&mon_long, comment, &status))
     PrintError( status );
    BlockOverlap  = (mon_long == 0) ? False : True;
    if (ffgkyj(fptr,(char*)"FormatIn", &mon_long, comment, &status)) PrintError( status );

    if (Verbose == True)
    {
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;
        if (BlockOverlap == False) cout << "No overlapping " << endl;
        if (UserBlockSize > 0)  cout << "BlockSize = " << UserBlockSize <<endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;
	cout << "Input image size = " << Nli << " " << Nci << endl << endl;
        if (BlockOverlap == True) cout << "Block Over " << endl;
        else cout << "NO Block Over " << endl;
     }
      VarNorm = True;

     alloc(Nli,Nci,UserBlockSize);
    FormatInputImag = (type_format)mon_long;

     Ptr = TabBand.buffer();
     if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

  if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;


#if DEBUG_IO
    cout << "Read out " << filename << endl;
#endif
}

/****************************************************************************/
