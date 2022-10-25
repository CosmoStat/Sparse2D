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
**    File:  MBase.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Decomposition of an imaghe on multiple bases
**    -----------
**
******************************************************************************/

#include "MBase.h"

/*************************************************************/
// MBase transform
/*************************************************************/

const char *StringMBase (type_mbase  type)
{
    switch (type)
    {
         case MB_PIX:
	      return ("Pixel Basis");
        case MB_ATROU:
	      return ("A trous algorithm");
        case MB_WT:
              return ("bi-orthogonal WT with 7/9 filters");
        case MB_CUR:
              return ("Curvelet transform");
        case MB_RID:
	      return ("Ridgelet transform");
        case MB_COS:
	      return ("Cosinus transform");
        case MB_FFT:
	      return ("FFT transform");
        case MB_PYRRID:
	      return ("Multi-Ridgelet");
        case MB_PMT:
	      return ("Pyramidal Median transform");
        default:
	      return ("Undefined transform");
     }
}

/***********************************************************************/

void MBase::reset()  // Reset all internal variables
{
    // for (int i=0; i < NBR_MBASE; i++) TabBase[i] = False;
    Bord = I_MIRROR;
    Nl = Nc = NbrBase = 0;
    FirstSoftThreshold = DEF_MB_FIRST_SOFT_THRESHOLD;
    LastSoftThreshold =  DEF_MB_LAST_SOFT_THRESHOLD;
    Nbr_Iter = DEF_MB_NBR_ITER;
    Filtering = False;
    PositivRecIma = False;
    DataSigmaNoise = 0.;
    Verbose = False;
    N_Sigma = 0.;
    UseHuberNorm = False;
    UseNormL1= False;
    AT_NbrScale2D = 4;
    AT_PositivRecIma = PositivRecIma;
    AT_KillLastScale = False;
    AT_SigmaNoise = 1.;
    AT_Trans = NULL;
    AT_Resi = NULL;
    AT_FirstDetectScale = 0;
    AT_AdjointRec = False;

    PMT_NbrScale2D = 4;
    PMT_PositivRecIma = PositivRecIma;
    PMT_KillLastScale = False;
    PMT_Trans = NULL;
    PMT_Resi = NULL;

    RID_BlockSize = -1;
    RID_FirstDetectScale = 0;
    RID_PositivRecIma = PositivRecIma;

    MRID_NbrRid = 4;
    MRID_FirstDetectScale = 0;
    MRID_PositivRecIma = PositivRecIma;
    MRID_BlockSize = 8;
    MRID_BlockOverlap = False;
    MRID = NULL;
    MRID_TabNl = MRID_TabNc = NULL;
    MRID_Resi = NULL;
    MRID_Trans = NULL;
    MRID_TabImaRec= NULL;

    CUR_NbrScale2D = 4;
    CUR_BlockSize = 16;
    CUR_PositivRecIma = PositivRecIma;
    CUR_BlockOverlap = False;
    CUR_KillLastScale= False;
    // CUR_Resi = NULL;
    // CUR_Trans = NULL;
    // CUR_TabNl = CUR_TabNc = NULL;

    WT_NbrScale2D = 4;
    WT_SelectFilter = NULL;
    WT_SB1D = NULL;
    WT_TabDec = NULL;
    WT = NULL;
    WT_PositivRecIma = PositivRecIma;
    WT_NbrUndecimatedScale = 1;

    COS_PositivRecIma = False;
    FFT_PositivRecIma = False;
}

/***********************************************************************/

void MBase::alloc()
{
   int i,s;
   Nl = Data.nl();
   Nc = Data.nc();
   Resi.alloc(Nl,Nc,"resi");
   Resi = Data;
   StepL = (FirstSoftThreshold - LastSoftThreshold ) / (float) Nbr_Iter;
   Lambda = FirstSoftThreshold;

   for (i =0; i < NbrBase; i++)
   {

        TabImaRec[i].alloc(Nl,Nc,StringMBase(TabSelect[i]));
        if (Verbose == True)
        cout << "Selection: " << StringMBase(TabSelect[i]) << endl;
        switch (TabSelect[i])
        {
             case MB_ATROU:
                    AWT.alloc(AT_Trans, Nl, Nc, AT_NbrScale2D);
                    AWT.alloc(AT_Resi, Nl, Nc, AT_NbrScale2D);
                    if (AT_KillLastScale == True) AT_PositivRecIma = False;
		    AWT.AdjointRec = AT_AdjointRec;
		    AWT.Bord = Bord;
                   break;
             case MB_PMT:
                    PMT.alloc(PMT_Trans, Nl, Nc, PMT_NbrScale2D);
                    PMT.alloc(PMT_Resi, Nl, Nc, PMT_NbrScale2D);
                    if (PMT_KillLastScale == True) PMT_PositivRecIma = False;
                   break;
             case MB_RID:
                   RID.alloc(Nl, Nc, RID_BlockSize);
                   NlRid = RID.rid_nl();
                   NcRid = RID.rid_nc();
                   RID_Trans.alloc(NlRid, NcRid, "RidTrans");
                   RID.VarianceStab=PoissonNoise;
                   if (RID.KillLastScale == True) RID_PositivRecIma = False;
                   if (Verbose == True)
                   {
                     cout << " Ridgelet: " << " Nl = " << NlRid  << " Nc = " <<  NcRid << endl;
                     cout << " BlockSize = " << RID_BlockSize << endl;
                   }
                  break;
             case MB_PYRRID:
                 // cout << "ALLOC MRID" << Nl << " "  << Nc << endl;
                   if (MRID_NbrRid < 2)
                   {
                       cout << "Error: bad number of ridgelet transform ... " << endl;
                       exit(-1);
                   }
                   MRID_TabBlockSize = new int [MRID_NbrRid];
                   MRID_TabBlockSize[0] = MRID_BlockSize;
                   for (s=1; s < MRID_NbrRid; s++)
                              MRID_TabBlockSize[s] = 2*MRID_TabBlockSize[s-1];
                   MRID_Trans  = new Ifloat [MRID_NbrRid];
                   MRID_Resi = new Ifloat [MRID_NbrRid];
                   MRID = new Ridgelet [MRID_NbrRid];
                   MRID_TabNl = new int [MRID_NbrRid];
                   MRID_TabNc = new int [MRID_NbrRid];
                   MRID_TabImaRec = new Ifloat [MRID_NbrRid];
                   if (MRID_KillLastScale == True) MRID_PositivRecIma = False;
                   if (Verbose == True)
                        cout << "MRID: Nbr scales = " <<  MRID_NbrRid << endl;
                // cout << "ALLOC MRID" << Nl << " "  << Nc << endl;

                   for (s=0; s < MRID_NbrRid; s++)
                   {
                      char ch[80];
                      MRID[s].VarianceStab=PoissonNoise;
                      MRID[s].BlockOverlap = MRID_BlockOverlap;
                      MRID[s].KillLastScale = MRID_KillLastScale;
                      MRID[s].alloc(Nl, Nc, MRID_TabBlockSize[s]);
                      MRID_TabNl[s]  = MRID[s].rid_nl();
                      MRID_TabNc[s]  = MRID[s].rid_nc();
                      sprintf (ch, "band_%d", s+1);
                      MRID_Trans[s].alloc(MRID_TabNl[s], MRID_TabNc[s], ch);
                      MRID_Resi[s].alloc (MRID_TabNl[s], MRID_TabNc[s], ch);
                      MRID_TabImaRec[s].alloc(Nl, Nc, ch);
                      if (Verbose == True)
                      {
                          cout << " Rid " << s+1 << " Nl = " <<  MRID_TabNl[s];
                          cout << " Nc = " << MRID_TabNc[s]  << " BlockSize = ";
                          cout <<  MRID_TabBlockSize[s]  << endl;
                          if (MRID[s].BlockOverlap == False) cout << "no overlap " << endl;
                          else  cout << "overlap " << endl;
                       }
                   }

                  // cout << "ALLOC MRID_Buff" << Nl << " "  << Nc << endl;
                   MRID_Buff.alloc(Nl,Nc,"MRID_Buff");
                   break;
             case MB_WT:
                 if (WT_NbrScale2D  < 2)
                 {
                     cout << "Error: bad number of scales ... " << endl;
                     exit(-1);
                 }
                 WT_SelectFilter = new FilterAnaSynt(F_MALLAT_7_9);
                 WT_SB1D = new SubBandFilter(*WT_SelectFilter, NORM_L2);
                 WT = new HALF_DECIMATED_2D_WT(*WT_SB1D);
                 WT_TabDec = new Bool[WT_NbrScale2D];
                 if ((WT_NbrUndecimatedScale < 0) || (WT_NbrUndecimatedScale > WT_NbrScale2D))
                      WT_NbrUndecimatedScale = WT_NbrScale2D;
                 for (s=0; s < MIN(WT_NbrScale2D,WT_NbrUndecimatedScale); s++) WT_TabDec[s] = False;
                 for (s=WT_NbrUndecimatedScale; s < WT_NbrScale2D; s++)  WT_TabDec[s] = True;
                 WT_NbrBand = WT->alloc(WT_Resi, Nl, Nc, WT_NbrScale2D,WT_TabDec);
                 WT_NbrBand = WT->alloc(WT_Trans, Nl, Nc, WT_NbrScale2D,WT_TabDec);
                 if (Verbose == True)
                 {
                    cout << "WT: Nbr scales = " << WT_NbrScale2D;
                    cout << " Nbr band = " << WT_NbrBand << endl;
                    //for (s=0; s <  WT_NbrBand; s++)
                    //    cout << " Band " << s+1 << " Nl = " << WT_Resi[s].nl() << " Nc = " << WT_Resi[s].nc() << endl;
                 }
                 break;
             case MB_PIX:
                  Pix_Trans.alloc(Nl, Nc, "cos trans");
                 break;
             case MB_COS:
                  COS_Resi.alloc(Nl, Nc, "cos resi");
                  COS_Trans.alloc(Nl, Nc, "cos trans");
                 break;
             case MB_FFT:
                  FFT_Resi.alloc(Nl, Nc, "cos resi");
                  FFT_Trans.alloc(Nl, Nc, "cos trans");
                 break;
             case MB_CUR:
                   CUR.Border = Bord;
                   CUR.NbrScale2D = CUR_NbrScale2D;
		   CUR.BlockOverlap = CUR_BlockOverlap;

                   // CUR.tab_block_size_init(CUR_BlockSize);
                   if (CUR_KillLastScale == True) CUR_PositivRecIma = False;
                   if (CUR_NbrScale2D  < 2)
                   {
                       cout << "Error: bad number scales in the curvelet transform ... " << endl;
                       exit(-1);
                   }
		   CUR.alloc(Nl, Nc, CUR_BlockSize);
		   CUR_Trans.alloc(CUR.cur_nc(), CUR.cur_nl(), CUR_NbrScale2D);
		   CUR_Resi.alloc(CUR.cur_nc(), CUR.cur_nl(), CUR_NbrScale2D);

                   // CUR_Trans = new Ifloat [CUR_NbrScale2D];
                   //CUR_Resi  = new Ifloat [CUR_NbrScale2D];
                   // CUR_Ridgelet = new Ridgelet [CUR_NbrScale2D];
                   // CUR_TabNl = new int [CUR_NbrScale2D];
                   // CUR_TabNc = new int [CUR_NbrScale2D];
                   if (Verbose == True)
                       cout << "Curvelet: Nbr scales = "<< CUR_NbrScale2D << endl;
//                    for (s=0; s < CUR_NbrScale2D-1; s++)
//                    {
//                       char ch[80];
//                       CUR_Ridgelet[s].BlockOverlap = CUR_BlockOverlap;
//                       CUR_Ridgelet[s].alloc(Nl, Nc, (CUR.TabBlockSize)(s));
//                       int Ns = CUR_Ridgelet[s].NbrScale;
//                     // cout << "Scale2D " << s << ":  Scale1D = " << Ns << endl;
//                       fltarray TN(Ns);
//                       for (int b=0; b < Ns; b++) TN(b)= CUR.norm_band(s,b);
//
//                       CUR_Ridgelet[s].set_tab_norm_with_tab(Ns, TN);
//                       CUR_TabNl[s]  = CUR_Ridgelet[s].rid_nl();
//                       CUR_TabNc[s]  = CUR_Ridgelet[s].rid_nc();
//                     // cout << "ALLOC IMA " << CUR_TabNl[s] << " " << CUR_TabNc[s] << endl;
//
//                       sprintf (ch, "band_%d", s+1);
//                       CUR_Trans[s].alloc(CUR_TabNl[s], CUR_TabNc[s], ch);
//                       CUR_Resi[s].alloc (CUR_TabNl[s], CUR_TabNc[s], ch);
//                       if (Verbose == True)
//                       {
//                           cout << " Cur " << s+1 << " Nl = " <<  CUR_TabNl[s];
//                           cout << " Nc = " << CUR_TabNc[s]  << endl;
//                           cout << " BlockSize = " <<  CUR.TabBlockSize(s)  << endl;
//                           //if (CUR_Ridgelet[s].BlockOverlap == False) cout << "no overlap " << endl;
//                           //else  cout << "overlap " << endl;
//                        }
//                    }
//                    CUR_Trans[CUR_NbrScale2D-1].alloc(Nl, Nc, "last");
//                    CUR_Resi[CUR_NbrScale2D-1].alloc(Nl, Nc, "last");
                   break;

	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
        }
   }

   if (NbrBase  < 1)
   {
       cout << "Error: no selected transform ... " << endl;
       exit(-1);
   }

   if (Verbose == True)
            cout << "Number of selected bases = " << NbrBase << endl;

   if (Filtering == False)
   {
      N_Sigma = 0.;
      AT_SigmaNoise = DataSigmaNoise; //  = 1.;
   }
   else
   {
      AT_SigmaNoise = DataSigmaNoise;
      if (N_Sigma == 0) N_Sigma = 3.;
   }
   if (PoissonNoise == True)
   {
      DataSigmaNoise = 1.;
#ifndef CURALONE
      AT_SigmaNoise = flux(Data) / (float) (Nl*Nc);
#else
      AT_SigmaNoise = Data.total() / (float) (Nl*Nc);
#endif
      if (Verbose == True) cout << "Mean data = " << AT_SigmaNoise  << endl;
      AT_SigmaNoise = sqrt(AT_SigmaNoise);
      if (Verbose == True) cout << "sqrt(Mean) = " << AT_SigmaNoise  << endl;
   }
}

/****************************************************************************/

void MBase::free()
{
   //if (NbrBase > 0) delete [] TabImaRec;
   for (int b=0; b < NbrBase; b++)
   {
      switch(TabSelect[b])
      {
            case MB_ATROU:
                  AWT.free(AT_Resi, AT_NbrScale2D);
                  AWT.free(AT_Trans, AT_NbrScale2D);
                  break;
             case MB_PMT:
                  PMT.free(PMT_Resi, PMT_NbrScale2D);
                  PMT.free(PMT_Trans, PMT_NbrScale2D);
                  break;
             case MB_RID:
                   break;
             case MB_PYRRID:
                  if (MRID_NbrRid > 2)
                  {
                      delete []  MRID_Resi;
                      delete []   MRID_Trans;
                      MRID_NbrRid = 0;
                      delete [] MRID_TabNl;
                      delete [] MRID_TabNc;
                      delete [] MRID;
                      delete [] MRID_TabBlockSize;
                      delete [] MRID_TabImaRec;
                  }
                  break;
             case MB_WT:
                 if  (WT_NbrScale2D > 2)
                 {
                    WT->free(WT_Resi, WT_NbrScale2D);
                    WT->free(WT_Trans, WT_NbrScale2D);
                    delete WT_SelectFilter;
                    delete WT_SB1D;
                    delete WT;
                    delete [] WT_TabDec;
                    WT_NbrScale2D = 0;
                    WT_SelectFilter = NULL;
                    WT_SB1D = NULL;
                    WT_TabDec = NULL;
                    WT = NULL;
                 }
                 break;
             case MB_FFT:
             case MB_COS:
	     case MB_PIX:
                 break;
             case MB_CUR:
//                  if  (CUR_NbrScale2D > 2)
//                  {
//                      delete [] CUR_Resi;
//                      delete [] CUR_Trans;
//                      delete [] CUR_Ridgelet;
//                      CUR_NbrScale2D = 0;
//                      delete [] CUR_TabNl;
//                      delete [] CUR_TabNc;
//                  }
                 break;
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
      }
   }
   NbrBase = 0;
}

/****************************************************************************/

void MBase::fft_proj(Ifloat & ImaRec, float FFTNoise, float NSigma, float Lamba_Sigma)
// Resi = IN: residual image
// FFT_Resi = Ifloat = IN: buffer for COS computation
// FFT_Trans = Ifloat  = IN/OUT decomposition on cos
// ImaRec = OUT: reconstruction from cos. transf.
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   FFTN_2D FFT;
   int i,j;
   float Noise = FFTNoise * Nc / sqrt(2.);
   float HardT = NSigma*Noise;
   float SoftT = Lamba_Sigma*Noise;
//    double Flux = flux(Resi) / (double)(Nl*Nc);
//    for (i=0; i < Nl; i++)
//    for (j=0; j < Nc; j++) Resi(i,j) -= Flux;

   FFT.fftn2d(Resi, FFT_Resi, False);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
       float Resi_re = (FFT_Resi)(i,j).real();
       float Resi_im = (FFT_Resi)(i,j).imag();
       float Coef_re = (FFT_Trans)(i,j).real();
       float Coef_im = (FFT_Trans)(i,j).imag();
       float Up_Re, Up_Im;
       if (UseHuberNorm == False)
       {
          Up_Re = update_coef(Coef_re, Resi_re, HardT, SoftT, UseNormL1);
          Up_Im = update_coef(Coef_im, Resi_im, HardT, SoftT, UseNormL1);
       }
       else
       {
          Up_Re  = hubert_update_coef(Coef_re, Resi_re, HardT, SoftT, Noise);
          Up_Im =  hubert_update_coef(Coef_im, Resi_im, HardT, SoftT, Noise);
       }
       FFT_Resi(i,j) = FFT_Trans(i,j) = complex_f(Up_Re,Up_Im);
       // if ((i!=0) || (j!=0)) COS_Trans(i,j) = update_coef(Coef, Resi, HardT, SoftT);
       // else COS_Trans(i,j) += Resi;
   }
   ImaRec.init();
   FFT_Trans(Nl/2,Nc/2) = FFT_Resi(Nl/2,Nc/2) = complex_f(0,0);
   FFT.fftn2d(FFT_Resi, True);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++) ImaRec(i,j) = FFT_Resi(i,j).real();
   INFO_X(ImaRec, (char*)"rec fft");
}

/****************************************************************************/

void MBase::cos_proj(Ifloat & ImaRec, float CosNoise, float NSigma, float Lamba_Sigma)
// Resi = IN: residual image
// COS_Resi = Ifloat = IN: buffer for COS computation
// COS_Trans = Ifloat  = IN/OUT decomposition on cos
// ImaRec = OUT: reconstruction from cos. transf.
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   int i,j;
   float Noise = CosNoise * Nc / 2.;
   float HardT = NSigma*Noise;
   float SoftT = Lamba_Sigma*Noise;
//    double Flux = flux(Resi) / (double)(Nl*Nc);
//    for (i=0; i < Nl; i++)
//    for (j=0; j < Nc; j++) Resi(i,j) -= Flux;

   im_dct(Resi, COS_Resi, False);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
       float Resi = (COS_Resi)(i,j);
       float Coef = (COS_Trans)(i,j);
       if (UseHuberNorm == False)
            COS_Trans(i,j) = update_coef(Coef, Resi, HardT, SoftT, UseNormL1);
       else COS_Trans(i,j) = hubert_update_coef(Coef, Resi, HardT, SoftT, Noise);
       // if ((i!=0) || (j!=0)) COS_Trans(i,j) = update_coef(Coef, Resi, HardT, SoftT);
       // else COS_Trans(i,j) += Resi;
   }
   ImaRec.init();
   im_dct(COS_Trans, ImaRec, True);
   // INFO_X(ImaRec, "rec cos");
}

/****************************************************************************/

void MBase::pix_proj(Ifloat & ImaRec, float PixNoise, float NSigma, float Lamba_Sigma)
// Resi = IN: residual image
// Pix_Resi = Ifloat = IN: buffer for PIX computation
// COS_Trans = Ifloat  = IN/OUT decomposition on cos
// ImaRec = OUT: reconstruction from cos. transf.
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   int i,j;
   float Noise = PixNoise;
   float HardT = NSigma*Noise;
   float SoftT = Lamba_Sigma*Noise;

   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++)
   {
       float Res = Resi(i,j);
       float Coef = (Pix_Trans)(i,j);
       if (UseHuberNorm == False)
            Pix_Trans(i,j) = update_coef(Coef, Res, HardT, SoftT, UseNormL1);
       else Pix_Trans(i,j) = hubert_update_coef(Coef, Res, HardT, SoftT, Noise);
       // if ((i!=0) || (j!=0)) Pix_Trans(i,j) = update_coef(Coef, Res, HardT, SoftT);
       // else Pix_Trans(i,j) += Res;
   }
   ImaRec = Pix_Trans;
}

/****************************************************************************/

void MBase::wt_proj(Ifloat & ImaRec, float WT_SigmaNoise, float NSigma,
                        float Lamba_Sigma)
// Resi = IN: residual image
// WT_Resi = Ifloat[0..WT_NbrScale2D-1] = IN: buffer for WT computation
// WT_Trans = Ifloat[0..WT_NbrScale2D-1] = IN/OUT decomposition on WT
// ImaRec = OUT: reconstruction from WT_Trans
// WT_NbrScale2D = IN: number of scales
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   int s,i,j;

   // WT of the residual
  // INFO_X(Resi, "Resi");
  WT->transform(Resi, WT_Resi, WT_NbrScale2D, WT_TabDec);
  // cout << "WT " << WT_NbrBand << endl;

   for (s=0; s < WT_NbrBand; s++)
   {
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      if (NSigma == 0) NSig = 0;
      float Noise =  WT_SigmaNoise;
      float Level = (s == WT_NbrBand-1) ? 0: NSig*Noise;
      float LWT = (s == WT_NbrBand-1) ? 0: Lamba_Sigma*Noise;

      for (i=0; i <  (WT_Trans[s]).nl(); i++)
      for (j=0; j <  (WT_Trans[s]).nc(); j++)
      {
          float Resi = (WT_Resi[s])(i,j);    // residual WT coef
          float Coef = (WT_Trans[s])(i,j); // A trous coef
          if (UseHuberNorm == False)
             (WT_Trans[s])(i,j) = update_coef(Coef, Resi, Level, LWT, UseNormL1);
          else
             (WT_Trans[s])(i,j)= hubert_update_coef(Coef, Resi, Level, LWT,Noise);
      }
   }
  // cout << "REC" << endl;

   WT->recons(WT_Trans,  ImaRec, WT_NbrScale2D, WT_TabDec);
}

/****************************************************************************/

void MBase::atrou_proj(Ifloat & ImaAtrou, float AT_SigmaNoise, float NSigma,
                        float Lamba_Sigma)
// Resi = IN: residual image
// AT_Resi = Ifloat[0..AT_NbrScale2D-1] = IN: buffer for WT computation
// AT_Trans = Ifloat[0..AT_NbrScale2D-1] = IN/OUT decomposition on
//                                           a trous algorithm
// ImaAtrou = OUT: reconstruction from AT_Trans
// AT_NbrScale2D = IN: number of scales
// AT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// AT_KillLastScale = IN: the last scale is taken into account
// Bord = IN: type of border management
{
   int s,i,j;
   int LastScale = (AT_KillLastScale == True) ? AT_NbrScale2D-1: AT_NbrScale2D;

   // WT of the residual
   // INFO_X(Resi, "Resi");
   AWT.transform(Resi,AT_Resi, AT_NbrScale2D);
   // cout << "ATROU Lambda = " << Lamba_Sigma << endl;
   for (s=AT_FirstDetectScale; s < LastScale; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      if (NSigma == 0) NSig = 0;
      float Norm = (s == AT_NbrScale2D-1) ? AWT.norm_band(s-1): AWT.norm_band(s);
      float Noise =  AT_SigmaNoise*Norm;
      float Level = (s == AT_NbrScale2D-1) ? 0: NSig*Noise;
      float LWT = Lamba_Sigma*Noise;

      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)
      {
          float Resi = (AT_Resi[s])(i,j);    // residual WT coef
          float Coef = (AT_Trans[s])(i,j); // A trous coef
          // Coef += Resi - sgn(Coef)*LWT;
          if (UseHuberNorm == False)
             (AT_Trans[s])(i,j) = update_coef(Coef, Resi, Level, LWT, UseNormL1);
          else (AT_Trans[s])(i,j) = hubert_update_coef(Coef, Resi, Level, LWT, Noise);
      }
   }
   AWT.recons(AT_Trans, ImaAtrou, AT_NbrScale2D);
}

/****************************************************************************/

void MBase::pmt_proj(Ifloat & ImaPmt, float PMT_SigmaNoise, float NSigma,
                        float Lamba_Sigma)
// Resi = IN: residual image
// PMT_Resi = Ifloat[0..PMT_NbrScale2D-1] = IN: buffer for PMT computation
// PMT_Trans = Ifloat[0..PMT_NbrScale2D-1] = IN/OUT decomposition on
//                                           PMT algorithm
// ImaPmt = OUT: reconstruction from PMT_Trans
// PMT_NbrScale2D = IN: number of scales
// PMT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// PMT_KillLastScale = IN: the last scale is taken into account
// Bord = IN: type of border management
{
   int s,i,j;
   int LastScale = (PMT_KillLastScale == True) ? PMT_NbrScale2D-1: PMT_NbrScale2D;

   // WT of the residual
   // INFO_X(Resi, "Resi");
   PMT.transform(Resi,PMT_Resi, PMT_NbrScale2D, Bord);
   // cout << "SOFT" << endl;

   for (s=0; s < LastScale; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      if (NSigma == 0) NSig = 0;
      float Norm = (s == PMT_NbrScale2D-1) ? PMT.norm_band(s-1): PMT.norm_band(s);
      float Noise =  PMT_SigmaNoise*Norm;
      float Level = (s == PMT_NbrScale2D-1) ? 0: NSig*Noise;
      float LWT = Lamba_Sigma*Noise;
      int Nls = (PMT_Resi[s]).nl();
      int Ncs = (PMT_Resi[s]).nc();

      for (i=0; i < Nls; i++)
      for (j=0; j < Ncs; j++)
      {
          float Resi = (PMT_Resi[s])(i,j);    // residual WT coef
          float Coef = (PMT_Trans[s])(i,j); // A trous coef
          // Coef += Resi - sgn(Coef)*LWT;
          if (UseHuberNorm == False)
            (PMT_Trans[s])(i,j) = update_coef(Coef, Resi, Level, LWT, UseNormL1);
          else (PMT_Trans[s])(i,j) = hubert_update_coef(Coef, Resi, Level, LWT, Noise);
      }
   }
   PMT.recons(PMT_Trans, ImaPmt, PMT_NbrScale2D, Bord);
}

/****************************************************************************/

void MBase::cur_proj(Ifloat & ImaCur, float CurNoise, float NSigma,
                      float Lamba_Sigma)
// Resi = IN: residual image
// CUR_Resi = Ifloat[0..AT_NbrScale2D-1] = IN: buffer for WT computation
// CUR_Trans = Ifloat[0..AT_NbrScale2D-1] = IN/OUT decomposition on
//                                           a trous algorithm
// ImaCur = OUT: reconstruction from the curvelet transform
// CUR_NbrScale2D = IN: number of scales
// CUR_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// CUR_KillLastScale = IN: the last scale is taken into account
// Bord = IN: type of border management
{
   int s2d,s1d,i,j;
   int Nls,Ncs,Depi,Depj;
   float Coef,CoefResi;
   // WT of the residual

   // INFO_X(Resi, "Resi");
   // CUR.Verbose = True;
   // CUR.transform(Resi, CUR_Resi, CUR_Ridgelet);
   CUR.transform(Resi, CUR_Resi);
   for (s2d=0; s2d < CUR_NbrScale2D-1; s2d++)
   {
      int NScale = CUR.nbr_rid_scale(s2d);
      // int NScale = CUR_Ridgelet[s2d].NbrScale;
      for (s1d=0; s1d < NScale; s1d++)
      {
         float NSig = (s1d == 0) ? NSigma + 1: NSigma;
         if (NSigma == 0) NSig = 0;
         float BS = ((CUR.TabBlockSize)(s2d) <= 0) ? Nc: (CUR.TabBlockSize)(s2d);
         // float Norm = (s1d == NScale-1) ? CUR_Ridgelet[s2d].rid_norm(s1d-1) : CUR_Ridgelet[s2d].rid_norm(s1d);
         float Norm = (s1d == NScale-1) ? CUR.norm_band(s2d,s1d-1) : CUR.norm_band(s2d,s1d);
         float Noise = sqrt(BS)*CurNoise*Norm;
         float Level = NSig*Noise;
         float LRid = Lamba_Sigma*Noise;
	 // cout << "CUR " << s2d << "  " << s1d << " Hard = " << Level << " Soft = " << LRid << endl;
         Nls =  CUR.size_nl(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nl(s1d);
         Ncs =  CUR.size_nc(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nc(s1d);
	 Depi = CUR.ipos(s2d, s1d); // CUR_Ridgelet[s2d].ipos(s1d);
         Depj = CUR.jpos(s2d, s1d); // CUR_Ridgelet[s2d].jpos(s1d);

         for (i=Depi; i < Depi+Nls; i++)
         for (j=Depj; j < Depj+Ncs; j++)
         {
            Coef = CUR_Trans(j,i,s2d);  // (CUR_Trans[s2d])(i,j);
            CoefResi = CUR_Resi(j,i,s2d);  // (CUR_Resi[s2d])(i,j);
            if (UseHuberNorm == False)
                 CUR_Trans(j,i,s2d) = update_coef(Coef,CoefResi,Level,LRid, UseNormL1);
            else CUR_Trans(j,i,s2d) =  hubert_update_coef(Coef,CoefResi,Level,LRid,Noise);
         }
      }
 //      s1d = NScale-1;
//       Nls = CUR_Ridgelet[s2d].size_scale_nl(s1d);
//       Ncs = CUR_Ridgelet[s2d].size_scale_nc(s1d);
//       Depi = CUR_Ridgelet[s2d].ipos(s1d);
//       Depj = CUR_Ridgelet[s2d].jpos(s1d);
//
//       for (i=Depi; i < Depi+Nls; i++)
//       for (j=Depj; j < Depj+Ncs; j++)
//             (CUR_Trans[s2d]) (i,j) += (CUR_Resi[s2d])(i,j);
   }
   s2d = CUR_NbrScale2D-1;
   if (CUR_KillLastScale == False)
   {
       //Nls = CUR_Trans[s2d].nl();
       //Ncs = CUR_Trans[s2d].nl();
       //for (i=0; i < Nls; i++)
       //for (j=0; j < Ncs; j++)
       for (i=0; i < Nl; i++)
       for (j=0; j < Nc; j++)
       {
          float Norm = AWT.norm_band(s2d-1);
          float Noise = CurNoise*Norm;
          Coef = CUR_Trans(j,i,s2d);
          CoefResi = CUR_Resi(j,i,s2d);
          float Level = 0.;
          float LRid = Lamba_Sigma*CurNoise;
          if (UseHuberNorm == False)
               CUR_Trans(j,i,s2d) = update_coef(Coef,CoefResi,Level,LRid, UseNormL1);
          else CUR_Trans(j,i,s2d) = hubert_update_coef(Coef,CoefResi,Level,LRid,Noise);
       }
       // CUR_Trans[s2d] += CUR_Resi[s2d];
   }
   // CUR.recons(CUR_Trans, ImaCur, CUR_Ridgelet);
   CUR.recons(CUR_Trans, ImaCur);
}

/****************************************************************************/

void MBase::rid_proj(Ridgelet & RG, Ifloat &RidTrans, Ifloat & ImaRid,
                     float SigmaNoise, float NSigma, Bool KillLastScale,
                     int BlockSize, float Lamba_Sigma, int FirstScale)
// Ridgelet base project of the residual
// Resi = IN: residual image
//  RID_Resi = Ifloat = IN: buffer for WT computation
//  RidTrans = Ifloat = IN/OUT decomposition on a ridgelet base
//  ImaRid = OUT: reconstruction from  RidTrans
// SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// KillLastScale = IN: the last scale is taken into account
{
   int s,i,j;
   int Nc = Resi.nc();
   int NScale = RG.NbrScale;
   int LastScale = (KillLastScale == True) ? NScale-1: NScale;
   float BS = (BlockSize <= 0) ? Nc: BlockSize;
   float RidNoise = (PoissonNoise == True) ? 1. : sqrt(BS)*SigmaNoise;

   RG.transform(Resi, RID_Resi, BlockSize);
   RG.set_tab_norm(RG.NbrScale);
   for (s=FirstScale; s < LastScale; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      if (NSigma == 0) NSig = 0;
      float Noise = (s == NScale-1) ? RidNoise*RG.rid_norm(s-1) :
                                               RidNoise*RG.rid_norm(s);
      float Level = (s == NScale-1) ? 0: NSig*Noise;
      float LRid = Lamba_Sigma*Noise;

      int Nls = RG.size_scale_nl(s);
      int Ncs = RG.size_scale_nc(s);
      int Depi = RG.ipos(s);
      int Depj = RG.jpos(s);

      for (i=Depi; i < Depi+Nls; i++)
      for (j=Depj; j < Depj+Ncs; j++)
      {
          float Coef = RidTrans(i,j);
          float Resi = RID_Resi(i,j);
          if (UseHuberNorm == False)
            RidTrans(i,j) = update_coef(Coef,Resi,Level,LRid, UseNormL1);
          else
            RidTrans(i,j) = hubert_update_coef(Coef,Resi,Level,LRid, Noise);

      }
   }
//    if (KillLastScale == False)
//    {
//       s = NScale-1;
//       int Nls = RG.size_scale_nl(s);
//       int Ncs = RG.size_scale_nc(s);
//       int Depi = RG.ipos(s);
//       int Depj = RG.jpos(s);
//
//       for (i=Depi; i < Depi+Nls; i++)
//       for (j=Depj; j < Depj+Ncs; j++)
//             RidTrans(i,j) += RID_Resi(i,j);
//    }
   RG.recons(RidTrans, ImaRid, BlockSize, Nl, Nc);
}

/****************************************************************************/

void MBase::reconstruction(Ifloat &Result)
{
   int b;
   Result = TabImaRec[0];
   for (b=1; b < NbrBase; b++)  Result += TabImaRec[b];
   if (PositivRecIma == True) threshold(Result);
}

/****************************************************************************/

void MBase::get_residual()
{
   int i,j;

   // Get the new solution
   reconstruction(Resi);

   // Get the Residual
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Resi(i,j) = Data(i,j) - Resi(i,j);
}

/****************************************************************************/

void MBase::decomposition()
{
   int i,j,b,s;

    Resi = Data;
    for (int Iter=0; Iter < Nbr_Iter; Iter++)
    {
       Lambda -= StepL;
       if (Lambda < 0) Lambda = 0.;
       // Lambda = 1.;

       // if (Verbose == True)
       if (Verbose == True)
               cout << endl << "Iteration " << Iter+1 << " Lambda = " << Lambda << endl;

       for (b=0; b < NbrBase; b++)
       {
          switch(TabSelect[b])
          {
              case MB_PIX:
	        if (Verbose == True) cout << "PIX proj " << endl;
	        pix_proj(TabImaRec[b], DataSigmaNoise,  N_Sigma, Lambda);
                break;
             case MB_PMT:
                 if (Verbose == True) cout << "PMT proj " << endl;
                  pmt_proj(TabImaRec[b], DataSigmaNoise, N_Sigma, Lambda);
                  if (PMT_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
            case MB_ATROU:
                  if (Verbose == True) cout << "Atrou proj " << endl;
                  atrou_proj(TabImaRec[b], AT_SigmaNoise, N_Sigma, Lambda);
                  if (AT_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
             case MB_RID:
                  if (Verbose == True) cout << "Rid proj " << endl;
                  rid_proj(RID, RID_Trans, TabImaRec[b], DataSigmaNoise,N_Sigma,
                             RID.KillLastScale, RID_BlockSize,
                             Lambda, RID_FirstDetectScale);
                  if (RID_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
             case MB_PYRRID:
                  if (Verbose == True) cout << "Multiscale Rid proj " << MRID_NbrRid  << endl;
                  // s = MRID_NbrRid-1;
                  // MRID[s].get_low_and_freq(Resi, TabImaRec[b], MRID_Buff,  MRID_TabBlockSize[s]);
                  // Resi = MRID_Buff;
                  // INFO_X(Resi, "residual");
                  for (s = MRID_NbrRid-1; s >= 0; s--)
                  // for (s = 0; s < MRID_NbrRid; s++)
                  {
                      MRID_Buff = MRID_TabImaRec[s];
                      rid_proj(MRID[s], MRID_Trans[s], MRID_TabImaRec[s], DataSigmaNoise,N_Sigma,
                                  MRID[s].KillLastScale,  MRID_TabBlockSize[s],
                                  Lambda, MRID_FirstDetectScale);
                      for (i = 0; i < Nl; i++)
                      for (j = 0; j < Nc; j++)
                      {
                         float Delta = (MRID_TabImaRec[s])(i,j) - MRID_Buff(i,j);
                         (TabImaRec[b])(i,j) += Delta;
                         Resi(i,j) -= Delta;
                      }
                     // Resi -= MRID_Buff;
                     // cout << "AFTER " << MRID_Buff.nl() << endl;
                     // INFO_X(TabImaRec[b], "rec");
                     // INFO_X(Resi, "residual");
                    // get_residual();
                  }
                  if (MRID_PositivRecIma == True) threshold(TabImaRec[b]);
                  // cout << "END Multiscale Rid proj " << TabImaRec[b].nl() << endl;
                  break;
             case MB_WT:
                  if (Verbose == True) cout << "WT proj " << endl;
                  wt_proj(TabImaRec[b],  DataSigmaNoise, N_Sigma, Lambda);
                  if (WT_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
             case MB_FFT:
                  if (Verbose == True) cout << "FFT proj " << endl;
                  fft_proj(TabImaRec[b], DataSigmaNoise,  N_Sigma, Lambda);
                  if (FFT_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
             case MB_COS:
                  if (Verbose == True) cout << "COS proj " << endl;
                  cos_proj(TabImaRec[b], DataSigmaNoise,  N_Sigma, Lambda);
                  if (COS_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
             case MB_CUR:
                  if (Verbose == True) cout << "Curvelet proj " << endl;
                  cur_proj(TabImaRec[b], DataSigmaNoise, N_Sigma, Lambda);
                  if (CUR_PositivRecIma == True) threshold(TabImaRec[b]);
                  break;
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }

          // New residual calculation
          get_residual();
	  if (Verbose == True)
	  {
             // cout << " Base " << b+1 << ": Reconstruction" << endl;
             INFO_X(TabImaRec[b], (char*)"   REC");
             INFO_X(Resi, (char*)"   Residual");
	  }
        } // ENDFOR (b=0 ...)
    } // ENDFOR (Iter=...)
}

/****************************************************************************/

void MBase::write_allima()
{
   for (int b=0; b < NbrBase; b++)
   {
      switch(TabSelect[b])
      {
         case MB_PIX:
               io_write_ima_float((char*)"xx_pix.fits", TabImaRec[b]); break;
         case MB_ATROU:
               io_write_ima_float((char*)"xx_atrou.fits", TabImaRec[b]); break;
          case MB_RID:
               io_write_ima_float((char*)"xx_rid.fits", TabImaRec[b]);
               break;
          case MB_PYRRID:
               io_write_ima_float((char*)"xx_mrid.fits", TabImaRec[b]);
               for (int i=0; i < MRID_NbrRid; i++)
               {
                  char ch[256];
                  sprintf(ch, (char*)"xx_mrid_r%d.fits", i+1);
                  io_write_ima_float(ch, MRID_TabImaRec[i]);
               }
               break;
          case MB_CUR:
               io_write_ima_float((char*)"xx_cur.fits", TabImaRec[b]);
               break;
          case MB_COS:
              io_write_ima_float((char*)"xx_cos.fits", TabImaRec[b]); break;
          case MB_FFT:
              io_write_ima_float((char*)"xx_fft.fits", TabImaRec[b]); break;
          case MB_WT:
              io_write_ima_float((char*)"xx_wt.fits", TabImaRec[b]); break;
          case MB_PMT:
              io_write_ima_float((char*)"xx_pmt.fits", TabImaRec[b]); break;
	  default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
       }
   }
   io_write_ima_float((char*)"xx_resi.fits", Resi);
}

/****************************************************************************/
