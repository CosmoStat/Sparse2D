/******************************************************************************
**                   Copyright (C) 2000 by CEA + Santford University
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
**    File:  Ridgelet.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the Ridgelet transform
**    ----------- 
******************************************************************************/

#ifndef _MBASE_H_
#define _MBASE_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM_DCT.h"
#include "FFTN_2D.h"

#ifndef CURALONE
#include "MR_Obj.h"
#include "IM_Sigma.h"
#endif 

#include "Ridgelet.h"
#include "Curvelet.h"

#define NBR_MBASE 9

enum type_mbase {MB_PIX,MB_ATROU,MB_RID,MB_PYRRID,MB_WT,MB_COS,MB_CUR,MB_PMT,MB_FFT,MB_UNKNOWN};

const char * StringMBase (type_mbase type);

/***********************************************************************/

// inline float soft_threshold(float Val, float T)
// {
//    float Coef = Val;
//    if (ABS(Coef) < T) Coef = 0.;
//    else if (Coef > 0) Coef -= T;
//         else Coef += T;
//    return Coef;
// }
// 
// /***********************************************************************/
// 
// inline float hard_threshold(float Val, float T)
// {
//    float Coef = Val;
//    if (ABS(Coef) < T) Coef = 0.;
//    return Coef;
// }

/***********************************************************************/

inline float update_coef(float Coef, float Resi, float HardThres, float SoftThreshold, Bool UseNormL1)
// L2 norm minimization with L0 norm penalty
{
   float NewCoef=Coef;
   if (UseNormL1 == False)
   {
      if ((ABS(Resi) > HardThres) ||  (ABS(Coef) > FLOAT_EPSILON))
                       NewCoef += hard_threshold(Resi, SoftThreshold);
   }
   else
   {
      if ((ABS(Resi) >  HardThres) ||  (ABS(Coef) > FLOAT_EPSILON))
              NewCoef  += Resi;                   
      NewCoef = soft_threshold(NewCoef,  SoftThreshold);
   }
   return NewCoef;
}

/***********************************************************************/

inline float hubert_update_coef(float Coef, float Resi, float HardThres, float SoftThreshold, float Sigma)
// L2 norm minimization with L1 norm penalty
{
   float NewCoef=Coef;
   
   if ((ABS(Resi) >  HardThres) ||  (ABS(Coef) > FLOAT_EPSILON))
   {
       if (ABS(Resi) < 1.35*Sigma) NewCoef += Resi;
       else if (Resi > 0) NewCoef += 1.35 * Sigma;      
            else NewCoef -= 1.35 * Sigma;
       NewCoef = soft_threshold(NewCoef, SoftThreshold);
   }
   return NewCoef;
}

/***********************************************************************/

#define DEF_MB_FIRST_SOFT_THRESHOLD 0.5  // First soft thresholding level 
#define DEF_MB_LAST_SOFT_THRESHOLD  0.5  // Last soft thresholding level 

#define DEF_MB_NBR_ITER 10             // Default number of iterations
                                       // for the decomposition

/***********************************************************************/

class MBase {
   // General parameters
   int Nl;                          // input image size
   int Nc;

   float Lambda;                    // Soft threshold parameter
   float StepL;                     // Decreasing step per iteration 
                                    // applied to the Lambda value


   // Parameter for the atrous algorithm
   // ----------------------------------
   ATROUS_2D_WT AWT;      // Atrou class
   Ifloat *AT_Resi;       // Buffer for the residual transformation
   float AT_SigmaNoise;   // Noise standard deviation for the atrous algo.

   // Parameter for the PMT algorithm
   //-------------------------------- 
   PMT_2D PMT;         // PMT class
   Ifloat *PMT_Resi;   // Buffer for the residual transformation

   // Parameter for the Ridgelet transform
   // -------------------------------------
   Ifloat RID_Resi;    // Buffer for the ridgelet transform of the residual
   int NlRid;
   int NcRid;          // Ridgelet image size

   // Parameter for the multiblock Ridgelet transform
   //------------------------------------------------
   Ridgelet *MRID;         // Array of ridgelet class
   int *MRID_TabBlockSize; // block size for each ridgelet transform
   Ifloat *MRID_Resi;      // Buffer for the multi-ridgelet transform of the residual
   Ifloat MRID_Buff;       // Addition buffer
   int *MRID_TabNl;        // Ridgelet images size (number of lines)
   int *MRID_TabNc;        // Ridgelet images size (number of columns)
   Ifloat *MRID_TabImaRec; // Reconstructed images

   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   HALF_DECIMATED_2D_WT *WT; // Pointer to the wavelet transform class
   Ifloat *WT_Resi;          // Buffer for the wavelet transform of the residual
   int WT_NbrBand;           // Number of bands in the transform
   FilterAnaSynt *WT_SelectFilter; // Selected filter bank
   SubBandFilter *WT_SB1D;  // Pointer to the fiter bank decomposition class
   Bool *WT_TabDec;         // Table indicating scales to be decimated

   // Parameter for the cosinus transform
   // -----------------------------------
   Ifloat COS_Resi;  // Cosinus bases decomposition
   
   // Parameter for the Pixel Basis
   // -----------------------------------
   
   // Parameter for the FFT transform
   // ----------------------------------
   Icomplex_f FFT_Resi;  // FFT bases decomposition

   // Parameter for the curvelet transform
   // ------------------------------------
   fltarray CUR_Resi;   // Buffer for the curvelet transform of the residual     
   // Ridgelet *CUR_Ridgelet; // Table of ridgelet class
   // int *CUR_TabNl;     // Curvelet images size (number of lines)
   // int *CUR_TabNc;     // Curvelet images size (number of columns)


   void atrou_proj(Ifloat & ImaAtrou, float AT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the a trous transform

   void rid_proj(Ridgelet & RG, Ifloat &RidTrans, Ifloat & ImaRid,
                 float SigmaNoise, float NSigma, Bool KillLastScale,
                 int BlockSize, float Lamba_Sigma, int FirstScale);
   // Projection on the ridgelet transform

   void wt_proj(Ifloat & ImaRec, float WT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the orthogonal WT

   void cos_proj(Ifloat & ImaRec, float Noise, float NSigma, float Lamba_Sigma);
   // Projection on the cosinus transform

  void pix_proj(Ifloat & ImaRec, float PixNoise, float NSigma, float Lamba_Sigma);

   void fft_proj(Ifloat & ImaRec, float Noise, float NSigma, float Lamba_Sigma);
   // Projection on the cosinus transform

   void cur_proj(Ifloat & ImaCur, float Cur_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the curvelet transform

  void pmt_proj(Ifloat & ImaPmt, float PMT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the pyramidal median transform

   void get_residual();
   // Get the reisudal

   public:

   // ***************** PUBLIC MBASE PART  *********************
  
   // Bool TabBase[NBR_MBASE]; // TabBase[i] = True when bases number i is selected
   type_mbase TabSelect[NBR_MBASE]; // selected bases
   int NbrBase;                     // number of selected bases                    

   Ifloat Data;             // Input data
   Ifloat Resi;             // residual
   Ifloat TabImaRec[NBR_MBASE]; // reconstructed image for each bases
   Bool PositivRecIma;   // If true, impose positive constraint
                         // to the finale reconstructed image
   type_border Bord;     // border management
   int Nbr_Iter;         // Number of iterations for the decomposition
   Bool Filtering;       // If true, noise is represented in the bases
   float N_Sigma;        // N_Sigma fixes the thresholding level
   float DataSigmaNoise; // Noise standard deviation
   Bool PoissonNoise;    // True if the noise follows a Poisson distribution
   Bool Verbose;         // Verbose mode
   Bool UseHuberNorm;    // If true, the hubert norm is used instead
                         // of the Xhi2 for the minimization
   Bool UseNormL1;       // If true, the L1 norm is used instead of the L0

   float FirstSoftThreshold; // First thresholding value
   float LastSoftThreshold;  // Last thresholding value

   // Parameter for the atrous algorithm
   // ----------------------------------
   Ifloat *AT_Trans;      // atrous bases decomposition
   int AT_NbrScale2D;     // Number of scales for the atrous transform
   Bool AT_KillLastScale; // do not used the last scale for the atrous 
                          // decomposition
   Bool AT_PositivRecIma; // If true, impose positive constraint
                          // to the reconstructed image from
                          // the atrou algorithm
   int AT_FirstDetectScale; // First detection scale in the \`a trous transform
   Bool AT_AdjointRec;    // if true, the reconstruction is done by
                          // adjoint reconstruction
			  
   // Parameter for the PMT algorithm
   //-------------------------------- 
   Ifloat *PMT_Trans;      // PMT bases decomposition
   int PMT_NbrScale2D;     // Number of scales for the PMT transform
   Bool PMT_KillLastScale; // do not used the last scale for the PMT 
                           // decomposition
   Bool PMT_PositivRecIma; // If true, impose positive constraint
                           // to the reconstructed image from
                           // the PMT algorithm

   // Parameter for the Ridgelet transform
   // -------------------------------------
   Ridgelet RID;             // Ridgelet transform class
   Ifloat RID_Trans;         // Ridgelet decomposition
   int RID_BlockSize;        // Ridgelet block size
   int RID_FirstDetectScale; // First detection scale in the ridgelet transform
   Bool RID_PositivRecIma;   // If true, impose positive constraint  

   // Parameter for the multiblock Ridgelet transform
   //------------------------------------------------
   Ifloat *MRID_Trans;        // atrous bases decomposition
   int MRID_BlockSize;        // block size for the first ridgelet transform
   int MRID_NbrRid;           // Number of ridgelet transform
   int MRID_FirstDetectScale; // First detection scale
   Bool MRID_BlockOverlap;    // If true, block overlapping
   Bool MRID_PositivRecIma;   // If true, impose positive constraint
   Bool MRID_KillLastScale;   // If true, kill last scale

   // Parameter for the curvelet transform
   // ------------------------------------
   Curvelet CUR;            // Curvelet class
   int  CUR_BlockSize;      // block size at the first scale
   Bool CUR_BlockOverlap;   // block overlapping
   Bool CUR_PositivRecIma;  // If true, impose positive constraint
   int  CUR_NbrScale2D;     // Number of scales in the a trou algorithm
   Bool CUR_KillLastScale;  // If true, last scale is killed
   fltarray CUR_Trans;       // Curvelet bases decomposition

   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   int WT_NbrUndecimatedScale; // Number of undecimated scale
   int WT_NbrScale2D;          // Number of scales
   Bool WT_PositivRecIma;      // If true, impose positive constraint
   Ifloat *WT_Trans;           // Wavelet bases decomposition

   // Parameter for the cosinus transform
   // -----------------------------------
   Ifloat COS_Trans;          // Cosinus bases decomposition
   Bool COS_PositivRecIma;    // If true, impose positive constraint

   // Parameter for the Pixel basis
   // -----------------------------------
   Ifloat Pix_Trans;          // Pixel bases decomposition
   
   // Parameter for the FFT transform
   // -----------------------------------
   Icomplex_f FFT_Trans;          // Cosinus bases decomposition
   Bool FFT_PositivRecIma;    // If true, impose positive constraint


   // Routines
   // --------
   MBase ()  {reset();}   // Constructor
   void reset();          // Reset all internal variables
   void alloc();          // Allocate ayy memory

   void decomposition (); // Start the decomposition
                          // alloc MUST have been called before
                           
   void reconstruction(Ifloat &Result);
                          // Return the reconstruction image from the 
                          // decompositions

   void write_allima();   // Write the reconstruction of each decompositon
                          // to the disk

   void free();           // deallocation of the class
   ~MBase(){free(); reset();}
};

/***********************************************************************/

#endif
