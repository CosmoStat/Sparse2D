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
**    DESCRIPTION Pyramidal Curvelet transform and reconstruction
**    -----------
**
******************************************************************************/

#include "Curvelet.h"
#include "IM_IO.h"
#include "RidNoise.h"
#include "PrimeNumber.h"
#include "IM_Prob.h"
#include "PCur.h"

#define NBR_NOISE_2DSCALE 5
#define NBR_NOISE_1DSCALE 5
#define NBR_NOISE_FINITE_1DSCALE 10

#define CUR_DEBUG 0
#define CUR_THRESHOLD_ESTIMATION 0

// static float TabCurvNorm_sigma[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
//     {{0.7430,0.0497,0.0094,0.0039,0.0027},
//     {0.1919,0.0913,0.0252,0.0054,0.0015},
//     {0.0484,0.0763,0.0516,0.0166,0.0038},
//     {0.0108,0.0289,0.0480,0.0345,0.0117},
//     {0.0029,0.0073,0.0197,0.0327,0.0241}};

/*
static float TabPCurvNorm[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
    {{0.692084, 0.0488825, 0.00915464, 0.00390845, 0.00247424},
    {0.184628, 0.0853162, 0.0254304, 0.00558103, 0.00150869},
    {0.0489787, 0.0737917, 0.0490959, 0.016418, 0.00425767},
    {0.0118884, 0.0291999, 0.0459632, 0.0330577, 0.0122221},
    {0.00407141, 0.00783881, 0.0200084, 0.0305753, 0.0232205}};

static float TabPCurvNormOrtho[NBR_NOISE_2DSCALE][NBR_NOISE_1DSCALE] =
    {{0.888171, 0.693059, 0.316226, 0.133076, 0.0792442},
    {0.0920643, 0.355331, 0.4478, 0.258705, 0.112331},
    {0.0133684, 0.0687863, 0.290654, 0.41791, 0.244424},
    {0.00712275, 0.017331, 0.0713842, 0.27431, 0.404305},
    {0.00264867, 0.00782153, 0.0244814, 0.0768586, 0.275698}};

static float TabPCurvNormFinite[NBR_NOISE_2DSCALE][NBR_NOISE_FINITE_1DSCALE] =
    {{0.907205,   0.848266,  0.750752, 0.706939, 0.693558, 0.769639, 0.757993, 0.876104, 0.943324, 1.34013},
    {0.0575738,   0.219989,  0.308451, 0.285516, 0.278624, 0.297142, 0.32716,  0.402356, 0.487739, 0.691129},
    {0.00477841,  0.0294306, 0.126726, 0.205192, 0.189233, 0.233786, 0.239031, 0.299415, 0.282413, 0.358558},
    {0.00182506, 0.00365026, 0.0238334,0.104979, 0.160795, 0.158492, 0.213155, 0.297539, 0.241235, 0.373624},
    {0.000782209,0.00130493, 0.00312119,0.0230133,0.0992632,0.15913, 0.185009, 0.194789, 0.240767, 0.259313}};
*/

/****************************************************************************/

float PyrCurvelet::norm_band(int Inds2d, int Inds1d)
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
      case RID_PYR_DIRECT:
      case RID_FSS_ORTHO:
      case RID_FSS_PYR:
         if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
         if ((s1d < 0) || (s1d >= NBR_NOISE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
         break;
      case RID_FINITE:
         if ((s2d < 0) || (s2d >= NBR_NOISE_2DSCALE)) s2d = NBR_NOISE_2DSCALE-1;
         if ((s1d < 0) || (s1d >= NBR_NOISE_FINITE_1DSCALE)) s1d = NBR_NOISE_1DSCALE-1;
         break;
     case RID_UNKNOWN:
         cout << "Error: unknown ridgelet transform (in PyrCurvelet::norm_band) ... " << endl;
         exit(0);
         break;
   }
   return 0; // TabCurSigma(s1d,s2d);
}


/****************************************************************************/

void PyrCurvelet::alloc(int Nl, int Nc, int BS)
{
   int b;
   extern type_format Format_Imag;
   FormatInputImag = Format_Imag;

   if (Verbose == True) cout << "Class allocation ... " << endl;
   tab_block_size_init(BS);

   if (Verbose == True) cout << "Wavelet allocation ... " << endl;
   type_sb_filter SB_Filter = F_MALLAT_7_9;
   FilterAnaSynt *PtrFAS = new FilterAnaSynt;
   if ((TransformWT2D == TO_MALLAT) || (TransformWT2D == TO_UNDECIMATED_MALLAT))
   {
      PtrFAS->Verbose = Verbose;
      PtrFAS->alloc(SB_Filter);
   }
   int NbrUndec = -1;
   sb_type_norm Norm = NORM_L2;
   MR_Data.alloc (Nl, Nc, NbrScale2D, TransformWT2D, PtrFAS, Norm, NbrUndec);
   NbrBand2D = MR_Data.nbr_band();
   MR_Data.ExactPyrRec = True;

   if (Verbose == True) cout << "Ridgelet allocation ... " << endl;
   if (TabRidgelet != NULL) delete [] TabRidgelet;
   else TabRidgelet = new Ridgelet[NbrBand2D-1];
   SizeTabRid = NbrBand2D-1;
   AllocClass = True;

   TabRidTrans = new Ifloat [NbrBand2D];
   for (b = 0; b <  NbrBand2D-1; b++)
   {
      int s = MR_Data.band_to_scale(b);
      if (Verbose == True) cout <<  "INIT RIDGELET Scale " << b+1 << "Block size = " << TabBlockSize(s) << endl;
      ridgelet_init( TabRidgelet[b]);
      TabRidgelet[b].alloc(MR_Data.size_band_nl(b), MR_Data.size_band_nl(b), TabBlockSize(s));
      TabRidTrans[b].alloc(TabRidgelet[b].rid_nl(),TabRidgelet[b].rid_nc(),"Scale Trans");
   }
   b = NbrBand2D-1;
   TabRidTrans[b].alloc(MR_Data.size_band_nl(b), MR_Data.size_band_nc(b),"Last Scale");

   NlIma=Nl;NcIma=Nc;

   NbrBandCur = 1;
   for (b=0; b < NbrBand2D-1; b++)  NbrBandCur+=TabRidgelet[b].NbrScale;

   TabCurSigma = new fltarray [NbrBand2D];
   for (b = 0; b < NbrBand2D-1; b++)
   {
      //int s2d, s1d;
      //get_scale_number(b, s2d, s1d);
      // cout << "Beg alloc " << TabRidgelet[s2d].NbrScale << " " << TabRidgelet[s2d].block_size() << endl;
      TabCurSigma[b].alloc(TabRidgelet[b].NbrScale, TabRidgelet[b].block_size());
   }
   // cout << "End alloc " << endl;
}

/****************************************************************************/

void PyrCurvelet::reset(SubBand1D *SB1D)
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
   ColTrans= False;
   RidTrans=RID_PYR_FFT;
   TypeBlock = DEF_TYPE_CUR_BLOCK;
   UserBlockSize=DEF_CUR_BLOCK_SIZE;
   OddBlockSize=True;
   // tab_block_size_init(DEF_CUR_BLOCK_SIZE);
   TransformWT2D = TO_PYR_BSPLINE;
   TabRidTrans= NULL;
   NbrBandCur= 0;
}

/****************************************************************************/

void PyrCurvelet::tab_block_size_init(int BS)
{
    CPRIME_NUMBER CPN;
    // Bool Odd = (BS % 2 == 0) ? False: True;
    int s, B = BS;
    // cout << " Type block = " << StringBlockType(TypeBlock) << endl;

    if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
    if (SetTransform(TransformWT2D) == TRANSF_MALLAT) B /= 2;
     if (SetTransform(TransformWT2D) == TRANSF_FEAUVEAU) B /= 2;

    for (s=0; s < MAX_SCALE; s++)
    {
        // cout << "TabBlockSize = " << s+1 << " " << B << endl;
	TabBlockSize(s) = B;
        if ((SetTransform(TransformWT2D) == TRANSF_PYR)  ||
	         (SetTransform(TransformWT2D) == TRANSF_MALLAT)||
	         (SetTransform(TransformWT2D) == TRANSF_FEAUVEAU))
        {
	    if  ((s!= 0) && (s % 1 == 0))   B /= 2;
        }
	else if ((SetTransform(TransformWT2D) == TRANSF_PAVE) ||
	         (SetTransform(TransformWT2D) == TRANSF_UNDECIMATED_MALLAT) ||
	         (SetTransform(TransformWT2D) ==  TRANSF_DIADIC_MALLAT))
        {
	   if (s % 2 == 0)   B *= 2;
        }
	else if (SetTransform(TransformWT2D) == TRANSF_SEMIPYR)
	              if (s != 0) B = iround((float) B / sqrt(2.));

       if (B < 8) B = 8;

       if (RidTrans == RID_FINITE)
                 B = (int) CPN.next_prime_number((int) BS);
    }
    if ((RidTrans != RID_FINITE) && (OddBlockSize == True))
    {
       for (s=0; s < MAX_SCALE; s++)
	     if (TabBlockSize(s) % 2 == 0) TabBlockSize(s) += 1;
    }
}

/****************************************************************************/

void PyrCurvelet::ridgelet_init(Ridgelet &Rid)
{
    // Initialize the ridgelet transform
    if (Ptr_SB1D != NULL) Rid.reset(Ptr_SB1D);
    Rid.BlockOverlap = BlockOverlap;
    Rid.ColTrans = ColTrans ;
    Rid.StatInfo = StatInfo;
    Rid.GetAutoNbScale = GetAutoNbScale;
    Rid.RidTrans = RidTrans;
    // Rid.Verbose = Verbose;
    Rid.WPTrans = WPTrans;
    if (GetAutoNbScale == False) Rid.NbrScale = NbrScale1D;
}

/****************************************************************************/

void PyrCurvelet::transform(Ifloat &Image)
{
    int b;
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
    if (SizeTabRid != NbrBand2D-1)
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }

    MR_Data.transform(Image);
    cout << " NbrBand2D = " << NbrBand2D  << endl;
    for (b = 0; b < NbrBand2D-1; b++)
    {
//        int s = MR_Data.band_to_scale(b);

        if (Verbose == True)
        {
           // cout << endl;
           cout << " BAND " << b+1 << " Block size = " << TabBlockSize(b) << " NbrScale WT1D = " << TabRidgelet[b].NbrScale << endl;
        }
        INFO_X(MR_Data.band(b), (char*)"BANDWT2D");
        TabRidgelet[b].transform(MR_Data.band(b), TabRidTrans[b]);
    }
    b = NbrBand2D-1;
    TabRidTrans[b] = MR_Data.band(b);
}

/****************************************************************************/

void PyrCurvelet::recons(Ifloat &Image)
{
    int b;
    char Name[256];
    if (AllocClass  == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    int Nl = NlIma;
    int Nc = NcIma;

    // if (Verbose == True)
     cout << "Reconstructed image size: Nl =  " << Nl  << " Nc = " << Nc << " NbrBand2D = " << NbrBand2D << endl;

    if ((Image.nl() != Nl) || (Image.nc() != Nc)) Image.resize(Nl,Nc);
    b = NbrBand2D-1;
    MR_Data.band(b) = TabRidTrans[b];
    for (b = NbrBand2D-2; b >= 0; b--)
    {
        int s = MR_Data.band_to_scale(b);
        // int Nlb = TabRidgelet[b].rid_nl();
	// int Ncb = TabRidgelet[b].rid_nc();
        if (Verbose == True)
	    cout << " RECBAND " << b+1 << " Block size = " << TabBlockSize(s) << " NbrScale WT1D = " << TabRidgelet[b].NbrScale << endl;
        sprintf(Name, "xx_rid_band_%d.rid", b+1);
	TabRidgelet[b].write(Name, TabRidTrans[b]);
	TabRidgelet[b].recons(TabRidTrans[b], MR_Data.band(b));
        INFO_X(MR_Data.band(b), (char*)"BANDWT2D");
    }
    MR_Data.recons(Image);
}

/****************************************************************************/

void PyrCurvelet::get_scale_number(int NumBand, int & s2d, int & s1d)
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
   if (NumBand == NbrBand-1) s2d = NbrBand2D-1;
   else
   {
      while (b+TabRidgelet[s2d].NbrScale <= NumBand)
      {
         b += TabRidgelet[s2d].NbrScale;
         s2d ++;
      // cout << " b = " << b << endl;
      }
      if (s2d == NbrBand2D-1) s1d = 0;
      else s1d = NumBand - b;
   }
   // cout << "NumBand = " << NumBand << " ==> s2d = " <<  s2d << ", s1d = " <<  s1d << endl;
}

/****************************************************************************/

void PyrCurvelet::get_band(int s2d, int s1d, Ifloat &Band)
{
    int i,j,Posi = ipos(s2d, s1d);
    int Posj = jpos(s2d, s1d);
    int Nl =  size_nl(s2d, s1d);
    int Nc =  size_nc(s2d, s1d);
    // cout << "  Nl " <<Nl << ", Nc " << Nc << " Pi" << Posi << " Pj" << Posj << endl;

    Band.resize(Nl,Nc);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++) Band(i,j) = TabRidTrans[s2d](i+Posi, j+Posj);
}

/****************************************************************************/

void PyrCurvelet::get_band(int NumBand, Ifloat &Band)
{
   int s2d, s1d;
   get_scale_number(NumBand, s2d, s1d);
   // cout << "GET band " << NumBand+1 << " ==> " << s2d+1 << " " << s1d+1 << endl;
   get_band(s2d, s1d, Band);
}

/****************************************************************************/

void PyrCurvelet::put_band(int s2d, int s1d, Ifloat &Band)
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
    for (j = 0; j < Nc; j++)  TabRidTrans[s2d](i+Posi, j+Posj) = Band(i,j);
}


/****************************************************************************/

void PyrCurvelet::put_band(int NumBand, Ifloat &Band)
{
   int s2d, s1d;
   get_scale_number(NumBand, s2d, s1d);
   put_band(s2d, s1d, Band);
}

/****************************************************************************/

int  PyrCurvelet::size_band_nl(int b)
{
   int s2d, s1d;
   get_scale_number(b, s2d, s1d);
   int Nl =  size_nl(s2d, s1d);
   return Nl;
}

/****************************************************************************/

int PyrCurvelet::size_band_nc(int b)
{
   int s2d, s1d;
   get_scale_number(b, s2d, s1d);
   int Nc =  size_nc(s2d, s1d);
   return Nc;
}

/****************************************************************************/

float & PyrCurvelet::operator() (int b, int i, int j)
{
   int s2d, s1d;
   get_scale_number(b, s2d, s1d);
   int Posi = ipos(s2d, s1d);
   int Posj = jpos(s2d, s1d);
   int Nl =  size_nl(s2d, s1d);
   int Nc =  size_nc(s2d, s1d);
   if ((i < 0) || (i >= Nl) || (j < 0) || (j >= Nc))
   {
      cout << "Error: bad coordinated ==> i = " << i << " j = " << j << endl;
      cout << "       Band size = " << Nl << " " << Nc  << endl;
      exit(-1);
   }
   return TabRidTrans[s2d](i+Posi, j+Posj);
}


/****************************************************************************/

void PyrCurvelet::get_norm_coeff(float N_Sigma)
{
   if (Verbose == True)
       cout << "get_norm_coeff " << N_Sigma << endl;
   float Sig=1.;
   unsigned int InitRnd = 10;

   Ifloat ImSimu(NlIma,NcIma,"ImSimu");
   im_noise_gaussian (ImSimu, Sig, InitRnd);
   transform(ImSimu);
   fltarray TabStat;
   // get_stat(TabStat);

   set_noise_level(N_Sigma);
}

/****************************************************************************/

void PyrCurvelet::get_norm_coeff(Ifloat &ImaNoise, float N_Sigma)
{
   if (Verbose == True)
       cout << "get_norm_coeff " << N_Sigma << endl;
   transform(ImaNoise);
   set_noise_level(N_Sigma);
}

/****************************************************************************/

int PyrCurvelet::block_size(int b)
{
   int s2d, s1d;
   get_scale_number(b, s2d, s1d);
   int BS = TabRidgelet[s2d].block_size();
   return BS;
}

/****************************************************************************/


void PyrCurvelet::set_noise_level(float N_Sigma)
{
    CImaProb CP;
//    int i,j;
	int b,s2d,s1d,Nlb,Ncb;
    double LMin, LMax;
    Ifloat Band;
    Ifloat LineBand;

    if (Verbose == True) cout << "Curvelet transform of the noise " << endl;
    for (b=0; b < nbr_band()-1; b++)
    {
       get_scale_number(b, s2d, s1d);
       int BS = block_size(b);
//       int NbrBlockNl = TabRidgelet[s2d].rid_block_nl();
       get_band(b, Band);
       Nlb = Band.nl();
       Ncb = Band.nc();
       cout <<  " Band " << b+1 << " Scale " << s2d+1 << " " << s1d  << " BS = " << BS << " RBS = " << TabRidgelet[s2d].block_size() << " Sigma = " << Band.sigma() << endl;
       CP.set(Band);
       CP.find_gthreshold(N_Sigma, LMin, LMax);
       float NoiseLevel = MAX(ABS(LMin), LMax) /  (N_Sigma * sqrt( (float) BS));
       // for (int bb=0; bb <BS; bb++) TabCurSigma[s2d](s1d,bb) = NoiseLevel;
        for (int bb=0; bb <BS; bb++) TabCurSigma[s2d](s1d,bb) = MAX(ABS(LMin), LMax) /  N_Sigma;
       if (Verbose == True)
             printf("Band(%2d,%2d): Noise Level = %f, Sig = %f, SigBand = %f\n", s2d+1, s1d+1,  TabCurSigma[s2d](s1d,0), NoiseLevel, Band.sigma());
     }
}

/****************************************************************************/

void PyrCurvelet::get_stat(fltarray &TabStat)
{
    int NbrStatPerBand = 5; // moment of order 2,3,4
    TabStat.alloc(nbr_band(), NbrStatPerBand);
    Ifloat Band;

    for (int b=0; b < nbr_band(); b++)
    {
        int N, s2d, s1d;
	double Mean, Sigma, Skew, Curt;
	float  Min, Max;
        get_scale_number(b, s2d, s1d);
	get_band(b, Band);
	N = Band.nl()*Band.nc();
        moment4(Band.buffer(), N, Mean, Sigma, Skew, Curt, Min, Max);
        TabStat(b, 0) = (float) Sigma;
	TabStat(b, 1) = (float) Skew;
	TabStat(b, 2) = (float) Curt;
        TabStat(b, 3) = (float) Min;
        TabStat(b, 4) = (float) Max;
        if (Verbose == True)
            printf("  Band %d (%d,%d): Nl = %d, Nc = %d, Sigma = %5.3f, Skew = %5.3f, Curt = %5.3f, Min = %5.3f, Max = %5.3f\n",
	           b+1, s2d+1, s1d+1, Band.nl(), Band.nc(), TabStat(b,0),
		   TabStat(b, 1),TabStat(b,2),TabStat(b, 3), TabStat(b, 4));
    }
}

/****************************************************************************/

void PyrCurvelet::threshold(float SigmaNoise, float N_Sigma)
{
    int i,j,b,s2d,s1d,Nlb,Ncb;
    // double LMin, LMax;
    Ifloat Band;
    if (Verbose == True) cout << "Noise thresholding ..." << endl;
    for (b=0; b < nbr_band()-1; b++)
    {
       get_scale_number(b, s2d, s1d);
       int s = MR_Data.band_to_scale(s2d);
       int BS = TabBlockSize(s);
       int Posi = ipos(s2d, s1d);
       int Posj = jpos(s2d, s1d);
       Nlb =  size_nl(s2d, s1d);
       Ncb =  size_nc(s2d, s1d);

       float Nsig = (s1d==0) ? N_Sigma+1: N_Sigma;
       //Nlb = size_band_nl(b);
       //Ncb = size_band_nc(b);
       // float Param  =  1. - ABS( (i % (Nl/2)) - Nl /4.) / (Nl/4.);
       // get_band(b, Band);
       for (i=0; i < Nlb; i++)
       for (j=0; j < Ncb; j++)
       {
           float Level = Nsig * SigmaNoise * TabCurSigma[s2d](s1d,i%BS); //  * sqrt((float) BS);
           if (ABS(TabRidTrans[s2d](i+Posi, j+Posj))  < Level)
	       TabRidTrans[s2d](i+Posi, j+Posj) = 0.;
       }
           // if (ABS(Band(i,j)) < Level) Band(i,j) = 0;
           //if (ABS((*this)(b,i,j)) < Level) (*this)(b,i,j) = 0;
       //put_band(b, Band);
       // if (Verbose == True)
       //    printf("Band(%2d,%2d) Nl = %d, Nc = %d: Threshold = %f\n", s2d+1, s1d+1, Nlb, Ncb, Level);
    }
}

/****************************************************************************/

void PyrCurvelet::wiener_filter(float SigmaNoise, int WienerBlockSize)
{
    int i,j,b,s2d,s1d,Nlb,Ncb;
    int B2 = WienerBlockSize / 2;
    // double LMin, LMax;
    Ifloat Band;
    if (Verbose == True) cout << "Wiener Noise Filtering ..." << WienerBlockSize <<  " Noise = " <<  SigmaNoise << endl;
    for (b=0; b < nbr_band()-1; b++)
    {
       get_scale_number(b, s2d, s1d);
       int s = MR_Data.band_to_scale(s2d);
       int BS = TabBlockSize(s);
       Nlb =  size_nl(s2d, s1d);
       Ncb =  size_nc(s2d, s1d);
       Band.resize(Nlb,Ncb);
       get_band(s2d, s1d, Band);
       printf("   Band(%2d,%2d) Nl = %d, Nc = %d: Max = %f\n", s2d+1, s1d+1, Nlb, Ncb, Band.max());
       float  VarianceSignal=0.,VarianceData=0., VarianceNoise;
       for (i=0; i < Nlb; i++)
       for (j=0; j < Ncb; j++)
       {
           VarianceData = 0.;
           for (int k=i-B2 ; k <= i+B2 ; k++)
           for (int l=j-B2 ; l <= j+B2 ; l++)
	       VarianceData += Band(k,l,I_MIRROR)*Band(k,l,I_MIRROR);
           VarianceData /= (float)(WienerBlockSize*WienerBlockSize);
           float NoiseLevel =  SigmaNoise * TabCurSigma[s2d](s1d,i%BS);
	   VarianceNoise = NoiseLevel*NoiseLevel;
           VarianceSignal = MAX(0,VarianceData-VarianceNoise);
	   if ((i==10) && (j==10)) printf(" VarianceNoise = %f, VarianceSignal = %f, VarianceData = %f, Coef = %f\n",
	        VarianceNoise,VarianceSignal,VarianceData,Band(i,j));
           Band(i,j) *= VarianceSignal  / VarianceData;
       }
       printf("Band(%2d,%2d) Nl = %d, Nc = %d: Max = %f\n", s2d+1, s1d+1, Nlb, Ncb, Band.max());
       put_band(s2d, s1d, Band);
   }
}

/****************************************************************************/
