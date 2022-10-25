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
**    File:  MRBase.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Decomposition of an imaghe on multiple bases
**    -----------
**
******************************************************************************/

#include "MRBase.h"
#include "IM_Prob.h"
#include "IM_Simu.h"
#include "IM_Sigma.h"
#include "FFTN_2D.h"
#include "IM_Deconv.h"
#include "MR_Sigma.h"
#include "IM_Regul.h"
#include "MR_SoftRegul.h"
#include "WPackets.h"

/*************************************************************/
// MRBase transform
/*************************************************************/

const char *StringMRBase (type_mbase  type)
{
    switch (type)
    {
        case MB_ATROU:
	      return ("A trous algorithm");
        case MB_WT:
              return ("bi-orthogonal WT with 7/9 filters");
        case MB_CUR:
              return ("Curvelet transform");
        case MB_FCUR:
              return ("Fast Curvelet transform");
	case MB_RID:
	      return ("Ridgelet transform");
        case MB_COS:
	      return ("Local Cosinus transform");
        case MB_PYRRID:
	      return ("Multi-Ridgelet");
        case MB_PMT:
	      return ("Pyramidal Median transform");
        case MB_WP:
	      return ("Wavelet Packet");
	case MB_WTMIRROR:
	      return ("Mirror Basis Wavelet transform");
	default:
	      return ("Undefined transform");
     }
}

/***********************************************************************/
 // static float mbr_update(float Coef, float CoefSol, float HardThres, float SoftThreshold)
//  {
//
//    float NewCoef=0;
//    float Resi = Coef - CoefSol;
//
//    if (ABS(Coef) >= HardThres)
//    {
//        if (ABS(CoefSol) < FLOAT_EPSILON) NewCoef = Coef;
//        else if (ABS(Resi) >= SoftThreshold)
//                 NewCoef = CoefSol + soft_threshold(Resi,SoftThreshold);
//    }
//    else if (CoefSol > HardThres) NewCoef = HardThres;
//    else if (CoefSol < -HardThres) NewCoef = -HardThres;
//
//    return (NewCoef);
// }

/***********************************************************************/

float MRBase::mbr_update(float CoefData, float CoefSol, float HardThres, float SoftThreshold, float Noise)
// L2 norm minimization with L1 norm penalty
{
//    float Resi = Coef - CoefSol;
//    if (ABS(Coef) < HardThres) Resi = 0;
//
//    return  soft_threshold(Resi,SoftThreshold);

   float NewCoef=CoefSol;
   float Resi = CoefData - CoefSol;
//   float Upper = CoefData+SoftThreshold;
//   float Lower = CoefData-SoftThreshold;
   if (ABS(CoefData) >= HardThres)
   // if (CoefData >= HardThres)
   {
       if (ABS(Resi) > Noise*DetectCoefTol)  NewCoef += Resi;
   }
   if (TotalVariation == True) return (NewCoef);
   else return (soft_threshold(NewCoef,SoftThreshold));
}
//    {
//        // NewCoef = CoefSol + soft_threshold(Resi,SoftThreshold);
//        // NewCoef = CoefSol + 0.5*Resi;
//        if (HardThres == 0)  NewCoef = CoefData;
//        else
//        {
//           if (NewCoef >  Upper) NewCoef = Upper;
// 	  else if (NewCoef < Lower ) NewCoef = Lower;
//        }
//        NewCoef = soft_threshold(NewCoef,SoftThreshold);
//    }
/***********************************************************************/

void MRBase::reset()  // Reset all internal variables
{
    // for (int i=0; i < NBR_MBASE; i++) TabBase[i] = False;
    Bord = I_MIRROR;
    Nl = Nc = NbrBase = 0;
    FirstSoftThreshold = DEF_MB_FIRST_SOFT_THRESHOLD;
    LastSoftThreshold =  DEF_MB_LAST_SOFT_THRESHOLD;
    Nbr_Iter = DEF_MB_NBR_ITER;
    Filtering = False;
    PositivRecIma = True;
    DataSigmaNoise = 0.;
    Verbose = False;
    N_Sigma = 0.;
    DetectCoefTol = 0.5;
    TotalVariation = False;
    LambdaTV = 0.1;

    AT_NbrScale2D = 4;
    AT_PositivRecIma = PositivRecIma;
    AT_KillLastScale = True;
    AT_SigmaNoise = 1.;
    AT_Trans = NULL;
    AT_Resi = NULL;

    PMT_NbrScale2D = 4;
    PMT_PositivRecIma = PositivRecIma;
    PMT_KillLastScale = True;
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
    MRID_Trans = NULL;

    CUR_NbrScale2D = 4;
    CUR_BlockSize = 16;
    CUR_PositivRecIma = PositivRecIma;
    CUR_BlockOverlap = False;
    CUR_KillLastScale= False;
    CUR_PositivRecIma = PositivRecIma;
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

    COS_PositivRecIma = PositivRecIma;
    COS_BlockSize = 32;
    COS_Overlapping = False;
    CosWeightFirst= False;

    WP = NULL;
    WPR = NULL;
    //WP_SelectFilter = NULL;
    //WP_SB1D = NULL;
    WP_PositivRecIma = False;
    WP_NbrScale2D = 4;
    WP_NbrUndecimatedScale = 2;
    WP_Filter = F_MALLAT_7_9;
    WP_NbrBand = 0;

    FCurData=NULL;
    FCurSol=NULL;
    FCUR_NbrDir = 16;
    FCUR_NbrScale2D = 0;
    FCUR_PositivRecIma = PositivRecIma;

    UseNoiseTab=False;
    MBData= NULL;
    MBSol= NULL;
    MB_NbrBand=0;
    MB_SelectFilter= NULL;
    MB_SB1D= NULL;
    MB_PositivRecIma= PositivRecIma;
    MB_NbrScale2D= 4;
}

/***********************************************************************/

void MRBase::alloc()
{
   int i,s;
   Nl = Data.nl();
   Nc = Data.nc();
   Resi.alloc(Nl,Nc,"resi");
   Result.alloc(Nl,Nc,"resi");
   // cout << "Image size: Nl = " << Nl << " " << " Nc = " << Nc << endl;
   Resi = Data;
   StepL = (FirstSoftThreshold - LastSoftThreshold ) / (float) Nbr_Iter;
   Lambda = FirstSoftThreshold;
   for (i =0; i < NbrBase; i++)
   {
       if (Verbose == True)
          cout << "Selection: " << StringMRBase(TabSelect[i]) << endl;
        switch (TabSelect[i])
        {
             case MB_ATROU:
                    AWT.alloc(AT_Trans, Nl, Nc, AT_NbrScale2D);
                    AWT.alloc(AT_Resi, Nl, Nc, AT_NbrScale2D);
                   break;
             case MB_PMT:
                    PMT.alloc(PMT_Trans, Nl, Nc, PMT_NbrScale2D);
                    PMT.alloc(PMT_Resi, Nl, Nc, PMT_NbrScale2D);
                   break;
             case MB_RID:
                   RID.alloc(Nl, Nc, RID_BlockSize);
                   NlRid = RID.rid_nl();
                   NcRid = RID.rid_nc();
                   RID_Trans.alloc(NlRid, NcRid, "RidTrans");
                   RID.VarianceStab=PoissonNoise;
                   // RID.KillLastScale = False;
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
                   MRID = new Ridgelet [MRID_NbrRid];
                   MRID_TabNl = new int [MRID_NbrRid];
                   MRID_TabNc = new int [MRID_NbrRid];
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
                      if (Verbose == True)
                      {
                          cout << " Rid " << s+1 << " Nl = " <<  MRID_TabNl[s];
                          cout << " Nc = " << MRID_TabNc[s]  << " BlockSize = ";
                          cout <<  MRID_TabBlockSize[s]  << endl;
                          if (MRID[s].BlockOverlap == False) cout << "no overlap " << endl;
                          else  cout << "overlap " << endl;
                       }
                   }
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
             case MB_WP:  {
	         // cout << "CB_WP" << endl;
	         if (WP_NbrScale2D  < 2)
                 {
                     cout << "Error: bad number of scales ... " << endl;
                     exit(-1);
                 }
                 FilterAnaSynt*  FAS = new FilterAnaSynt(WP_Filter);
                 SubBandFilter* WP1D = new SubBandFilter(*FAS, NORM_L2);
                 WP = new WPACKETS_2D(*WP1D);
                 WP->alloc(Nl, Nc, WP_NbrScale2D, WP_NbrUndecimatedScale);
 	         WPR = new WPACKETS_2D(*WP1D);
		 WPR->alloc(Nl, Nc, WP_NbrScale2D, WP_NbrUndecimatedScale);
		 if (Verbose == True) {
                    cout << "WP: Nbr scales = " << WP_NbrScale2D;
		    cout << " Nbr undecimated scales   = " << WP_NbrUndecimatedScale << endl;
                    cout << " Nbr band = " << WT_NbrBand << endl;
                  }
 		 }
                 break;
	     case MB_FCUR:
	          FCurData = new FCUR;
		  FCurSol = new FCUR;
		  FCurData->alloc_from_fine(FCUR_NbrScale2D,Nl,Nc,FCUR_NbrDir,False,False);
		  FCurSol->alloc_from_fine(FCUR_NbrScale2D,Nl,Nc,FCUR_NbrDir,False,False);
		  FCurData->get_norm_coeff(3.);
		 break;
	     case MB_COS:
                  COS_Resi.alloc(Nl, Nc, "cos resi");
                  COS_Trans.alloc(Nl, Nc, "cos trans");
                  LDCT.alloc(Nl, Nc, COS_BlockSize, COS_Overlapping, CosWeightFirst);
		  LDCT_Data.alloc(Nl, Nc, COS_BlockSize, COS_Overlapping, CosWeightFirst);
                  break;
             case MB_CUR:
                   CUR.Border = Bord;
                   CUR.NbrScale2D = CUR_NbrScale2D;
                   // CUR.tab_block_size_init(CUR_BlockSize);
		   // CUR.RidTrans = RID_FSS;
		   CUR.RidTrans = RID_PYR_FFT;
		   CUR.BlockOverlap = CUR_BlockOverlap;
                   if (CUR_NbrScale2D  < 2)
                   {
                       cout << "Error: bad number scales in the curvelet transform ... " << endl;
                       exit(-1);
                   }
		   CUR.alloc(Nl, Nc, CUR_BlockSize);
		   CUR_Trans.alloc(CUR.cur_nc(), CUR.cur_nl(), CUR_NbrScale2D);
		   CUR_Resi.alloc(CUR.cur_nc(), CUR.cur_nl(), CUR_NbrScale2D);


//                    CUR_Trans = new Ifloat [CUR_NbrScale2D];
//                    CUR_Resi  = new Ifloat [CUR_NbrScale2D];
//                    CUR_Ridgelet = new Ridgelet [CUR_NbrScale2D];
//                    CUR_TabNl = new int [CUR_NbrScale2D];
//                    CUR_TabNc = new int [CUR_NbrScale2D];
//                    if (Verbose == True)
//                        cout << "Curvelet: Nbr scales = "<< CUR_NbrScale2D << endl;
//                    for (s=0; s < CUR_NbrScale2D-1; s++)
//                    {
//                       char ch[80];
//                       CUR_Ridgelet[s].BlockOverlap = CUR_BlockOverlap;
//                       CUR_Ridgelet[s].set_block_param(Nl, Nc, (CUR.TabBlockSize)(s));
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
             case MB_WTMIRROR:
                 if (WT_NbrScale2D  < 2)
                 {
                     cout << "Error: bad number of scales ... " << endl;
                     exit(-1);
                 }
                 MB_SelectFilter = new FilterAnaSynt(F_MALLAT_7_9);
                 MB_SB1D = new SubBandFilter(*MB_SelectFilter, NORM_L2);
		 MBData = new MIRROR_2D_WT(*MB_SB1D);
		 MBSol = new MIRROR_2D_WT(*MB_SB1D);

                 MBData->alloc(Nl, Nc, MB_NbrScale2D);
                 MBSol->alloc(Nl, Nc, MB_NbrScale2D);
		 MB_NbrBand = MBSol->nbr_band();
                 if (Verbose == True)
                 {
                    cout << "MB: Nbr scales = " << MB_NbrScale2D;
                    cout << " Nbr band = " << MB_NbrBand << endl;
                 }
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

void MRBase::free()
{
   for (int b=0; b < NbrBase; b++)
   {
      switch(TabSelect[b])
      {
            case MB_FCUR:
	          delete FCurData;
		  delete FCurSol;
		  break;
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
                      delete []   MRID_Trans;
                      MRID_NbrRid = 0;
                      delete [] MRID_TabNl;
                      delete [] MRID_TabNc;
                      delete [] MRID;
                      delete [] MRID_TabBlockSize;
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
             case MB_WP:
	          delete WP;
		  delete WPR;
	         break;
             case MB_WTMIRROR:
	           MBData->free();
		   MBSol->free();
	           MB_NbrScale2D=0;
		   MB_SelectFilter = NULL;
                   MB_SB1D = NULL;
	           break;
             case MB_COS:
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
		 CUR_NbrScale2D = 0;
                 break;
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
      }
   }
   NbrBase = 0;
}
/****************************************************************************/

void MRBase::cos_proj(Ifloat & ImaRec, float CosNoise, float NSigma, float Lamba_Sigma)
// Resi = IN: residual image
// COS_Resi = Ifloat = IN: buffer for COS computation
// COS_Trans = Ifloat  = IN/OUT decomposition on cos
// ImaRec = OUT: reconstruction from cos. transf.
// WT_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
{
   int i,j;
   float Noise = CosNoise * LDCT.norm();
   int BS = LDCT.B2DTrans.block_size();
   Ifloat BlockCosImaRec(BS, BS, "blockima");
   Ifloat BlockCosImaData(BS, BS, "blockima");
   float HardT = NSigma*Noise;
   float SoftT = Lamba_Sigma*Noise;
//   float Coef, Level = NSigma*Noise;
   SoftT = 0.;
   float Mean = ImaRec.mean();
   LDCT.transform(ImaRec);

   for (int bi = 0; bi < LDCT.B2DTrans.nbr_block_nl(); bi++)
   for (int bj = 0; bj < LDCT.B2DTrans.nbr_block_nc(); bj++)
   {
       LDCT.B2DTrans.get_block_ima(bi,bj, LDCT.DCTIma, BlockCosImaRec);
       LDCT_Data.B2DTrans.get_block_ima(bi,bj, LDCT_Data.DCTIma, BlockCosImaData);
       for (i=0; i < BS; i++)
       for (j=0; j < BS; j++)
       {
           float CoefSol = BlockCosImaRec(i,j);
           float CoefData = BlockCosImaData(i,j);
           if ((i != 0) || (j!=0))
	    BlockCosImaRec(i,j) = mbr_update(CoefData, CoefSol, HardT, SoftT, Noise);
       }
       LDCT.B2DTrans.put_block_ima(bi,bj, LDCT.DCTIma, BlockCosImaRec);
   }

//    for (i=0; i < LDCT.nl(); i++)
//    for (j=0; j < LDCT.nc(); j++)
//    {
//       float CoefSol = LDCT(i,j);
//       float CoefData = LDCT_Data(i,j);
//       LDCT(i,j) = mbr_update(CoefData, CoefSol, HardT, SoftT, Noise);
//    }
   ImaRec.init();
   LDCT.recons(ImaRec);
   float Mean1 = ImaRec.mean();
   Mean /= Mean1;
   for (i=0; i < ImaRec.nl(); i++)
   for (j=0; j < ImaRec.nc(); j++) ImaRec(i,j) *= Mean;
}

/****************************************************************************/

void MRBase::wt_proj(Ifloat & ImaRec, float WT_SigmaNoise, float NSigma,
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

   WT->transform(Result, WT_Resi, WT_NbrScale2D, WT_TabDec);
   for (s=0; s < WT_NbrBand-1; s++)
   {
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      float Noise =  (UseNoiseTab == False) ? WT_SigmaNoise: TabWTSigma(s);
      float Level = (UseNoiseTab == False) ? NSig*Noise: TabWTLevel(s);
      float LWT = Lamba_Sigma*Noise;
      if (s == MB_NbrBand-1) Level = LWT = 0.;

      for (i=0; i <  (WT_Trans[s]).nl(); i++)
      for (j=0; j <  (WT_Trans[s]).nc(); j++)
      {
          float Resi = (WT_Resi[s])(i,j);    // residual WT coef
          float Coef = (WT_Trans[s])(i,j); // A trous coef
          (WT_Resi[s])(i,j) = mbr_update(Coef, Resi, Level, LWT, Noise);
      }
   }
   s = WT_NbrBand-1;
   for (i=0; i < (WT_Trans[s]).nl(); i++)
   for (j=0; j < (WT_Trans[s]).nc(); j++) (WT_Resi[s])(i,j) = (WT_Trans[s])(i,j);
   WT->recons(WT_Resi, ImaRec, WT_NbrScale2D, WT_TabDec);
}

/****************************************************************************/

void MRBase::fcur_proj(Ifloat & ImaRec, float FCUR_SigmaNoise, float NSigma, float Lamba_Sigma)
{
   int s,b,i,j;
   Ifloat BandSol;
   Ifloat BandData;

   FCurSol->cur_trans(Result);
   for (s=0; s <  FCurSol->nbr_scale()-1; s++)
   for (b=0; b <  FCurSol->nbr_band(s); b++)
   {
      FCurSol->get_band(s,b,BandSol);
      FCurData->get_band(s,b,BandData);
      float NSig = (s < 3) ? NSigma + 1: NSigma;
      float Noise =  FCUR_SigmaNoise * FCurData->norm_band(s,b);
      float Level =  NSig*Noise;
      float LWT = Lamba_Sigma*Noise;

      for (i=0; i <  BandSol.nl(); i++)
      for (j=0; j <  BandSol.nc(); j++)
      {
          float Resi = BandSol(i,j);    // residual WT coef
          float Coef = BandData(i,j); // A trous coef
          BandSol(i,j) = mbr_update(Coef, Resi, Level, LWT, Noise);
      }
      FCurSol->put_band(s,b,BandSol);
   }
   s = FCurSol->nbr_scale()-1;
   FCurData->get_band(s,0,BandSol);
   FCurSol->put_band(s,0,BandSol);

   FCurSol->cur_recons(ImaRec);
}

/****************************************************************************/

void MRBase::wp_proj(Ifloat & ImaRec, float WP_SigmaNoise, float NSigma,
                        float Lamba_Sigma)

{
   int s,i,j;

   WPR->transform(Result, WP_Resi);
   for (s=0; s < WP_TabTrans->nbr_band()-1; s++)
   {
      // float NSig = (s < 3) ? NSigma + 1: NSigma;
      float NSig =  NSigma;
      float Noise =  (UseNoiseTab == False) ? WP_SigmaNoise: TabWPSigma(s);
      float Level = (UseNoiseTab == False) ? NSig*Noise: TabWPLevel(s);
      float LPT = Lamba_Sigma*Noise;

      for (i=0; i <  WP_TabTrans->size_band_nl(s); i++)
      for (j=0; j <  WP_TabTrans->size_band_nc(s); j++)
      {
          float Resi = (*WP_Resi)(s,i,j);    // residual WT coef
          float Coef = (*WP_TabTrans)(s,i,j); // A trous coef
          (*WP_Resi)(s,i,j) = mbr_update(Coef, Resi, Level, LPT, Noise);
       }
   }
   s = WP_TabTrans->nbr_band()-1;
   for (i=0; i < WP_TabTrans->size_band_nl(s); i++)
   for (j=0; j < WP_TabTrans->size_band_nc(s); j++) (*WP_Resi)(s,i,j) = (*WP_TabTrans)(s,i,j);
   WPR->recons(WP_Resi, ImaRec);
}

/****************************************************************************/

void MRBase::mb_proj(Ifloat & ImaRec, float MB_SigmaNoise, float NSigma,
                     float Lamba_Sigma)
{
   int s,i,j;

   MBSol->transform(Result);
   for (s=0; s < MB_NbrBand-1; s++)
   {
      int Nlb = MBSol->size_band_nl(s);
      int Ncb = MBSol->size_band_nc(s);
      float NSig =  NSigma;
      float Noise = (UseNoiseTab == False) ? MB_SigmaNoise:TabMBSigma(s);
      float Level = (UseNoiseTab == False) ? NSig*Noise:TabMBLevel(s);;
      float LWT = Lamba_Sigma*Noise;
      if (s == MB_NbrBand-1) Level = LWT = 0.;

      for (i=0; i < Nlb; i++)
      for (j=0; j < Ncb; j++)
      {
         float Sol =  (*MBSol)(s,i,j);    // residual WT coef
         float Coef = (*MBData)(s,i,j); // A trous coef
         (*MBSol)(s,i,j) = mbr_update(Coef, Sol, Level, LWT, Noise);
      }
   }
   s = MB_NbrBand-1;
   for (i=0; i < MBSol->size_band_nl(s); i++)
   for (j=0; j < MBSol->size_band_nc(s); j++) (*MBSol)(s,i,j) = (*MBData)(s,i,j);
   MBSol->recons(ImaRec);
}

/****************************************************************************/

void MRBase::atrou_proj(Ifloat & ImaAtrou, float AT_SigmaNoise, float NSigma,
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
   int LastScale = AT_NbrScale2D;

   // WT of the residual
   // INFO_X(Resi, "Resi");
   AWT.transform(Result,AT_Resi, AT_NbrScale2D);
   // cout << "SOFT" << endl;

   for (s=0; s < LastScale-1; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
      float Norm = (s == AT_NbrScale2D-1) ? AWT.norm_band(s-1): AWT.norm_band(s);
      float Noise =  AT_SigmaNoise*Norm;
      float Level = (s == AT_NbrScale2D-1) ? 0: NSig*Noise;
      float LWT = (s == AT_NbrScale2D-1) ?  0: Lamba_Sigma*Noise;

      for (i=0; i < Nl; i++)
      for (j=0; j < Nc; j++)
      {
          float Resi = (AT_Resi[s])(i,j);    // residual WT coef
          float Coef = (AT_Trans[s])(i,j); // A trous coef
          (AT_Resi[s])(i,j) = mbr_update(Coef, Resi, Level, LWT, Noise);
       }
   }
   s = LastScale-1;
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++) (AT_Resi[s])(i,j) = (AT_Trans[s])(i,j);
   AWT.recons(AT_Resi, ImaAtrou, AT_NbrScale2D);
}

/****************************************************************************/

void MRBase::pmt_proj(Ifloat & ImaPmt, float PMT_SigmaNoise, float NSigma,
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
   int LastScale = PMT_NbrScale2D;

   // WT of the residual
   // INFO_X(Resi, "Resi");
   PMT.transform(Result, PMT_Resi, PMT_NbrScale2D, Bord);
   // cout << "SOFT" << endl;

   for (s=0; s < LastScale-1; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
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
          (PMT_Resi[s])(i,j) = mbr_update(Coef, Resi, Level, LWT, Noise);
      }
   }
   s = LastScale-1;
   for (i=0; i < (PMT_Resi[s]).nl(); i++)
   for (j=0; j < (PMT_Resi[s]).nc(); j++) (PMT_Resi[s])(i,j) = (PMT_Trans[s])(i,j);
   PMT.recons(PMT_Resi, ImaPmt, PMT_NbrScale2D, Bord);
}

/****************************************************************************/

void MRBase::cur_proj(Ifloat & ImaCur, float CurNoise, float NSigma,
                      float Lamba_Sigma)
// Resi = IN: residual image
// CUR_Resi = Ifloat[0..AT_NbrScale2D-1] = IN: buffer for WT computation
// CUR_Trans = Ifloat[0..AT_NbrScale2D-1] = IN/OUT decomposition on
//                                           a trous algorithm
// ImaCur = OUT: reconstruction from the curvelet transform
// CUR_NbrScale2D = IN: number of scales
// CUR_SigmaNoise = IN: noise standard deviation
// NSigma = IN: N sigma filtering
// Bord = IN: type of border management
{
   int b,s2d,s1d,i,j;
   int Nls,Ncs,Depi,Depj;
   float Coef,CoefResi;
   // WT of the residual

   // INFO_X(Resi, "Resi");
   // CUR.Verbose = True;
   CUR.transform(Result, CUR_Resi);

   for (b=0; b < CUR.nbr_band()-1; b++)
   {
      CUR.get_scale_number(b,s2d,s1d);
      int NScale = CUR.nbr_rid_scale(s2d);
      float NSig = (s1d == 0) ? NSigma + 1: NSigma;
      float BS = ((CUR.TabBlockSize)(s2d) <= 0) ? Nc: (CUR.TabBlockSize)(s2d);
      // float Norm = (s1d == NScale-1) ? CUR_Ridgelet[s2d].rid_norm(s1d-1) : CUR_Ridgelet[s2d].rid_norm(s1d);
      float Norm = (s1d == NScale-1) ? CUR.norm_band(s2d,s1d-1) : CUR.norm_band(s2d,s1d);
      float Noise = (UseNoiseTab == False) ? sqrt(BS)*CurNoise*Norm: TabCurSigma(b);
      float Level = (UseNoiseTab == False) ? NSig*Noise: TabCurLevel(b);
      float LRid = Lamba_Sigma*Noise;
      if (b == CUR.nbr_band()-1)   Level = LRid = 0.;
      Nls =  CUR.size_nl(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nl(s1d);
      Ncs =  CUR.size_nc(s2d, s1d); // CUR_Ridgelet[s2d].size_scale_nc(s1d);
      Depi = CUR.ipos(s2d, s1d); // CUR_Ridgelet[s2d].ipos(s1d);
      Depj = CUR.jpos(s2d, s1d); // CUR_Ridgelet[s2d].jpos(s1d);

      for (i=Depi; i < Depi+Nls; i++)
      for (j=Depj; j < Depj+Ncs; j++)
      {
         Coef = CUR_Trans(j,i,s2d);
         CoefResi = CUR_Resi(j,i,s2d);
         CUR_Resi(j,i,s2d) = mbr_update(Coef,CoefResi,Level,LRid, Noise);
      }
   }
   s2d = CUR_NbrScale2D-1;
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++) CUR_Resi(j,i,s2d) = CUR_Trans(j,i,s2d);
   CUR.recons(CUR_Resi, ImaCur);
}

/****************************************************************************/

void MRBase::rid_proj(Ridgelet & RG, Ifloat &RidTrans, Ifloat & ImaRid,
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

   RG.transform(Result, RID_Resi, BlockSize);
   RG.set_tab_norm(RG.NbrScale);

   for (s=FirstScale; s < LastScale; s++)
   {
      float NSig = (s == 0) ? NSigma + 1: NSigma;
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
          RID_Resi(i,j) = mbr_update(Coef,Resi,Level,LRid,Noise);
      }
   }
   RG.recons(RID_Resi, ImaRid, BlockSize, Nl, Nc);
}

/*****************************************************************************/

void MRBase::init_decomposition()
{
   int b,s;
    for (b=0; b < NbrBase; b++)
    {
      switch(TabSelect[b])
      {
       case MB_PMT: PMT.transform(Data, PMT_Trans, PMT_NbrScale2D, Bord);
                    break;
       case MB_ATROU:
                    AWT.transform(Data,AT_Trans, AT_NbrScale2D); break;
       case MB_RID: RID.transform(Data, RID_Trans, RID_BlockSize); break;
       case MB_PYRRID:
                  for (s = MRID_NbrRid-1; s >= 0; s--)
                      MRID[s].transform(Data, MRID_Trans[s], MRID_TabBlockSize[s]);
                  break;
       case MB_WT: WT->transform(Data, WT_Trans, WT_NbrScale2D, WT_TabDec);
                   break;
       case MB_WP: WP->transform(Data,WP_TabTrans);
                    break;
       case MB_COS: LDCT_Data.transform(Data);  break;
       case MB_CUR: CUR.transform(Data, CUR_Trans); break;
       case MB_WTMIRROR: MBData->transform(Data);break;
       case MB_FCUR:  FCurData->cur_trans(Data); break;
       default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
      }
    }
}

/*****************************************************************************/

void MRBase::decomposition()
{
   // RegulIma RIM;
   MR_Regul RIM;
//   int i,j;
   int b,s;
//   float Converg=1.;
   // RIM.NbrScale =
   // Ifloat ImaRec(Nl,Nc,"ima rec");
   Ifloat FilterData(Nl,Nc,"ima rec");
   Ifloat ResiData(Nl,Nc,"ima rec");
   Result.init();
   init_decomposition();
   RIM.ExpDecreasingLambda = True;
   RIM.NbrScale = 2; // WT_NbrScale2D;
   // RIM.TypeFilter= F_HAAR; // F_MALLAT_7_9;
    for (int Iter=0; Iter < Nbr_Iter; Iter++)
    {
       Lambda -= StepL;
       if (Lambda < 0) Lambda = 0.;
       // if ((TotalVariation == True) && (LambdaTV > 0.5)) Converg = 1. / (2. * Lambda);
       //else Converg = 1.;

       // if (Verbose == True)
       if (Verbose == True)
       {
          cout << endl << "Iteration " << Iter+1 << " Lambda = " << Lambda << endl;
       }
       FilterData = Result;
       for (b=0; b < NbrBase; b++)
       {
          switch(TabSelect[b])
          {
             case MB_PMT:
                 if (Verbose == True) cout << "PMT proj " << endl;
                  pmt_proj(Result, DataSigmaNoise, N_Sigma, Lambda);
                  break;
            case MB_ATROU:
                  if (Verbose == True) cout << "Atrou proj " << endl;
                  atrou_proj(Result, AT_SigmaNoise, N_Sigma, Lambda);
                  break;
             case MB_RID:
                  if (Verbose == True) cout << "Rid proj " << endl;
                  rid_proj(RID, RID_Trans, Result, DataSigmaNoise,N_Sigma,
                             RID.KillLastScale, RID_BlockSize,
                             Lambda, RID_FirstDetectScale);
                    break;
             case MB_PYRRID:
                  if (Verbose == True) cout << "Multiscale Rid proj " << MRID_NbrRid  << endl;
                  for (s = MRID_NbrRid-1; s >= 0; s--)
                  {
                      rid_proj(MRID[s], MRID_Trans[s], Result, DataSigmaNoise,N_Sigma,
                                  MRID[s].KillLastScale,  MRID_TabBlockSize[s],
                                  Lambda, MRID_FirstDetectScale);
                      // Result = ImaRec;
                  }
                  break;
             case MB_WT:
                  if (Verbose == True) cout << "WT proj " << endl;
                  wt_proj(Result,  DataSigmaNoise, N_Sigma, Lambda);
                   break;
             case MB_WP:
                  if (Verbose == True) cout << "WP proj " << endl;
                  wp_proj(Result,  DataSigmaNoise, N_Sigma, Lambda);
		  break;
             case MB_FCUR:
                  if (Verbose == True) cout << "Fast Curvelet proj " << endl;
		  fcur_proj(Result, DataSigmaNoise,  N_Sigma,   Lambda);
		  break;
	     case MB_COS:
                  if (Verbose == True) cout << "COS proj " << endl;
                  cos_proj(Result, DataSigmaNoise,  N_Sigma, Lambda);
                  break;
             case MB_CUR:
                 if (Verbose == True) cout << "Curvelet proj " << endl;
                  cur_proj(Result, DataSigmaNoise, N_Sigma, Lambda);
                  break;
	     case MB_WTMIRROR:
	          if (Verbose == True) cout << "WT Mirror proj " << endl;
		  mb_proj(Result, DataSigmaNoise, N_Sigma, Lambda);
		  break;
	     default:
               cout << "Error: not implemeted transform ... " << endl;
               exit(-1);
          }
           // if (TabSelect[b] !=  MB_PYRRID) Result = ImaRec;
	   // Result = ImaRec;
      } // ENDFOR (b=0 ...)

      if ((TotalVariation == True) && (LambdaTV > 0))
      {
//          // calculate the gradient
//          for (i=0; i < Result.nl(); i++)
//          for (j=0; j < Result.nc(); j++)
//  	    ResiData(i,j) = Result(i,j) - FilterData(i,j);
//
//  	 RIM.obj_regul(FilterData, ResiData, LambdaTV*DataSigmaNoise);
// 	               // laplacien_regul(Result, Lambda);
//
// 	 for (i=0; i < Result.nl(); i++)
// 	 for (j=0; j < Result.nc(); j++)
// 	        Result(i,j) = FilterData(i,j) + Converg*ResiData(i,j);
         RIM.im_soft_threshold(Result, Result, LambdaTV*DataSigmaNoise*Lambda, LambdaTV*DataSigmaNoise*Lambda);

	 if (PositivRecIma == True) threshold(Result);
         if (Verbose == True)
         {
 	     ResiData = Data - Result;
	     INFO_X(Result, (char*)"SOL");
	     INFO_X(ResiData, (char*)"RESI");
         }
       }  // end TV
       else if (PositivRecIma == True) threshold(Result);
    } // ENDFOR (Iter=...)
}

/****************************************************************************/

void MRBase::init_deconv(Ifloat & ImaPSF, float SigmaNoise, float Eps, int InitRnd)
{
   Ifloat NoiseData(Nl,Nc,"noise data");
   if (SigmaNoise < FLOAT_EPSILON)
                            SigmaNoise = detect_noise_from_med (Data);
   im_noise_gaussian (NoiseData, SigmaNoise, InitRnd);
   init_deconv(ImaPSF, NoiseData, Eps);
}

/****************************************************************************/

void MRBase::init_deconv(Ifloat & ImaPSF, Ifloat &NoiseData, float Eps)
{
   CImaProb CP;
   int b;
   double LMin, LMax;
   int Nlb, Ncb;
   fltarray Band;
   Ifloat IBand;
   FFTN_2D FFT2D;
   int Nl1 = Nl;
   int Nc1 = Nc;
   Ifloat OutPsf(Nl1, Nc1, "OPsf");
   Icomplex_f Psf_cf(Nl1, Nc1, "OPsfcf");
   dec_center_psf (ImaPSF, OutPsf);
   FFT2D.fftn2d(OutPsf, Psf_cf);
   if (Verbose == True) cout << "Inverse the PSF " << endl;
   dec_inverse(NoiseData, Psf_cf, NoiseData, Eps);
   io_write_ima_float("xx_inv_noise.fits", NoiseData);
   dec_inverse(Data, Psf_cf, Data, Eps);
   io_write_ima_float("xx_inv_data.fits", Data);
   UseNoiseTab = True;

   for (int ba=0; ba < NbrBase; ba++)
   {
       switch(TabSelect[ba])
       {
	  case MB_CUR:
             // Curvelet Transform of the inverse solution
             TabCurLevel.alloc(CUR.nbr_band());
             TabCurSigma.alloc(CUR.nbr_band());
	     if (Verbose == True) cout << "Curvelet transform of the noise " << endl;
             CUR.transform(NoiseData, CUR_Trans);
             for (b=0; b <  CUR.nbr_band(); b++)
             {
                 int s2d, s1d;
                 CUR.get_scale_number(b, s2d, s1d);
                 CUR.get_band(CUR_Trans, b, Band);
                 Nlb = Band.ny();
                 Ncb = Band.nx();
                 IBand.alloc(Band.buffer(),Nlb, Ncb);
                 // cout << "Band " << b+1 << " Size = " << IBand.nl() << " " << IBand.nc() << endl;
                 CP.set(IBand);
                 if (s1d == 0) CP.find_gthreshold(N_Sigma+1, LMin, LMax);
                 else CP.find_gthreshold(N_Sigma, LMin, LMax);
                 TabCurLevel(b) = MAX(ABS(LMin), LMax);
		 // TabCurSigma(b) =  sigma(IBand);
                 TabCurSigma(b) = TabCurLevel(b) / N_Sigma;
             }
	     break;
         case MB_WTMIRROR:
             // Mirror Wavelet Transform of the inverse solution
             if (Verbose == True) cout << "Mirror wavelet transform of the noise " << endl;
             MBSol->transform(NoiseData);
             TabMBLevel.alloc(MBSol->nbr_band());
             TabMBSigma.alloc(MBSol->nbr_band());
             for (b=0; b < MBSol->nbr_band(); b++)
             {
                  Nlb = MBSol->size_band_nl(b);
                  Ncb = MBSol->size_band_nc(b);
                  MBSol->get_band(IBand, b);
                  CP.set(IBand);
                  CP.find_gthreshold(N_Sigma, LMin, LMax);
                  TabMBLevel(b) = MAX(ABS(LMin), LMax);
                  // TabMBSigma(b) = sigma(IBand);
		  TabMBSigma(b) = TabMBLevel(b) / N_Sigma;
              }
	      break;
         case MB_WT:
	     if (Verbose == True) cout << "Bi-ortho wavelet transform of the noise " << endl;
             WT->transform(NoiseData, WT_Trans, WT_NbrScale2D, WT_TabDec);
	     TabWTLevel.alloc(WT_NbrBand);
             TabWTSigma.alloc(WT_NbrBand);
             for (b=0; b < WT_NbrBand; b++)
             {
                  Nlb = (WT_Trans[b]).nl();
                  Ncb = (WT_Trans[b]).nc();
                  CP.set(WT_Trans[b]);
                  CP.find_gthreshold(N_Sigma, LMin, LMax);
                  TabWTLevel(b) = MAX(ABS(LMin), LMax);
                  // TabWTSigma(b) = sigma(WT_Trans[b]);
		  TabWTSigma(b) = TabWTLevel(b) / N_Sigma;
              }
	      break;
         default:
	     cout << "Error: this transform is not implemented for the deconvolution ..." << endl;
	     exit(-1);
      }
   }
}

/****************************************************************************/
