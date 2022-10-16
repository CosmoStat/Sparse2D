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

#ifndef _MRBASE_H_
#define _MRBASE_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM_DCT.h"
#include "WT_Mirror.h"
#include "IM_DCT.h"
#include "IM_DCT.h"
#include "WPackets.h"
#include "FCur.h"

#ifndef CURALONE
#include "MR_Obj.h"
#include "IM_Sigma.h"
#endif 

#include "Ridgelet.h"
#include "Curvelet.h"

#define NBR_MRBASE 9
#define NBR_USED_MRBASE 6

enum type_mbase {MB_ATROU, MB_WT, MB_RID, MB_CUR, MB_COS, MB_FCUR, MB_WP, MB_WTMIRROR, MB_PYRRID, MB_PMT, MB_UNKNOWN};
const char * StringMRBase (type_mbase type);

/***********************************************************************/

#define DEF_MB_FIRST_SOFT_THRESHOLD 1.  // First soft thresholding level 
#define DEF_MB_LAST_SOFT_THRESHOLD  0.  // Last soft thresholding level 

#define DEF_MB_NBR_ITER 10             // Default number of iterations

/***********************************************************************/

class MBDeconv;
class MRBase {
   friend class MBDeconv;
   // General parameters
   int Nl;                          // input image size
   int Nc;

   float Lambda;                    // Soft threshold parameter
   float StepL;                     // Decreasing step per iteration 
                                    // applied to the Lambda value
   Bool UseNoiseTab;                // In case of deconvolution
                                    // UseNoiseTab is to set to True
				    // It means that the noise std deviation
				    // in each band must be read from the
				    // specific relative arrays.

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
   fltarray TabWTLevel;
   fltarray TabWTSigma;
   
   // Parameter for the mirror wavelet transform
   //-----------------------------------------------
   MIRROR_2D_WT *MBData; // Pointer to the wavelet transform class
   MIRROR_2D_WT *MBSol; // Pointer to the wavelet transform class
   int MB_NbrBand;           // Number of bands in the transform
   FilterAnaSynt *MB_SelectFilter; // Selected filter bank
   SubBandFilter *MB_SB1D;  // Pointer to the fiter bank decomposition class
   fltarray TabMBLevel;
   fltarray TabMBSigma;
       
   // Parameter for the cosinus transform
   // -----------------------------------
   LOCAL_DCT2D LDCT; // Local cosinus transform class
   LOCAL_DCT2D LDCT_Data; // Local cosinus transform class
   Ifloat COS_Resi;  // Cosinus bases decomposition


   // Parameter for the curvelet transform
   // ------------------------------------
   fltarray CUR_Resi;   // Buffer for the curvelet transform of the residual     
   fltarray TabCurLevel;
   fltarray TabCurSigma;
   // Ridgelet *CUR_Ridgelet; // Table of ridgelet class
   // int *CUR_TabNl;     // Curvelet images size (number of lines)
   // int *CUR_TabNc;     // Curvelet images size (number of columns)

    // Parameter for the Packets WT
   // -----------------------------------
   WPACKETS_2D *WP;         // Wavelet Packet Class
   WPACKETS_2D *WPR;         // Wavelet Packet Class
   FilterAnaSynt *WP_SelectFilter; // Selected filter bank
   SubBandFilter *WP_SB1D;  // Pointer to the fiter bank decomposition class
   int WP_NbrBand;           // Number of bands in the transform
   WPTransf_2D* WP_Resi;        
   fltarray TabWPLevel;
   fltarray TabWPSigma;
   
   // Parameter for the Fast Curvelet Transform
   // -----------------------------------
   FCUR *FCurData; // Pointer the Fast curvelet transform CLASS
   FCUR *FCurSol; // Pointer the Fast curvelet transform CLASS
   
   
   
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
   
   void wp_proj(Ifloat & ImaRec, float WP_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the orthogonal WP
   
   void cos_proj(Ifloat & ImaRec, float Noise, float NSigma, float Lamba_Sigma);
   // Projection on the cosinus transform

   void cur_proj(Ifloat & ImaCur, float Cur_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the curvelet transform

  void pmt_proj(Ifloat & ImaPmt, float PMT_SigmaNoise, float NSigma, 
                        float Lamba_Sigma);
   // Projection on the pyramidal median transform

   void mb_proj(Ifloat & ImaRec, float MB_SigmaNoise, float NSigma, 
                     float Lamba_Sigma);
   // Projection on the mirror basis WT		     

   void  fcur_proj(Ifloat & ImaRec, float FCUR_SigmaNoise, float NSigma, float Lamba_Sigma);

		     
  float mbr_update(float CoefData, float CoefSol, float HardThres, 
                   float SoftThreshold, float Noise);
   // Calculate the update coefficient of the transform, taking into acount:
   // CoefData: IN = coefficient of the input data
   // CoefSol: IN = coefficient of the solution
   // HardThres: IN =  detection level for CoefData
   // SoftThreshold: IN =  Regularization parameter
   // Noise: IN =  Noise level
 
  
   public:
   void init_decomposition();

   // ***************** PUBLIC MRBASE PART  *********************
  
   // Bool TabBase[NBR_MRBASE]; // TabBase[i] = True when bases number i is selected
   type_mbase TabSelect[NBR_MRBASE]; // selected bases
   int NbrBase;                     // number of selected bases                    

   Ifloat Data;             // Input data
   Ifloat Resi;             // residual
   Ifloat Result;             // Input data
   Bool PositivRecIma;   // If true, impose positive constraint
                         // to the finale reconstructed image
   type_border Bord;     // border management
   int Nbr_Iter;         // Number of iterations for the decomposition
   Bool Filtering;       // If true, noise is represented in the bases
   float N_Sigma;        // N_Sigma fixes the thresholding level
   float DataSigmaNoise; // Noise standard deviation
   Bool PoissonNoise;    // True if the noise follows a Poisson distribution
   Bool Verbose;         // Verbose mode
   float FirstSoftThreshold; // First thresholding value
   float LastSoftThreshold;  // Last thresholding value
   float DetectCoefTol; // coefficient of the transformed solution must
                        // satisfy the constraint:
		        //  | c_s - c_d | < DetectCoefTol*SigmaNoise
			//  when | c_d | > DetectionLevel
			// c_s is the coefficient of the transformed solution
			// c_d is the coefficient of the transformed data
   Bool TotalVariation; // If true, the total variation is minimized instead
                        // of the L_1 norm of the multiscale coefficient
   float LambdaTV;      // Regularization parameter for the TV method.
   
   // Parameter for the atrous algorithm
   // ----------------------------------
   Ifloat *AT_Trans;      // atrous bases decomposition
   int AT_NbrScale2D;     // Number of scales for the atrous transform
   Bool AT_KillLastScale; // do not used the last scale for the atrous 
                          // decomposition
   Bool AT_PositivRecIma; // If true, impose positive constraint
                          // to the reconstructed image from
                          // the atrou algorithm

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

   // Parameter for the mirror wavelet transform
   //-----------------------------------------------
   int MB_NbrScale2D;          // Number of scales
   Bool MB_PositivRecIma;      // If true, impose positive constraint
   
    // Parameter for the mirror wavelet transform
   //-----------------------------------------------
   int FCUR_NbrScale2D;          // Number of scales
   int FCUR_NbrDir;
   Bool FCUR_PositivRecIma;      // If true, impose positive constraint
   
   // Parameter for the cosinus transform
   // -----------------------------------
   Ifloat COS_Trans;          // Cosinus bases decomposition
   Bool COS_PositivRecIma;    // If true, impose positive constraint
   int  COS_BlockSize;
   Bool COS_Overlapping;
   Bool CosWeightFirst;

   
    // Parameter for the Packets WT
   // -----------------------------------
   int WP_NbrScale2D;          // Number of scales
   int WP_NbrUndecimatedScale; // Number of undecimated scale
   Bool WP_PositivRecIma;    // If true, impose positive constraint
   Bool WP_KillLastScale;    
   type_sb_filter WP_Filter;
   WPTransf_2D* WP_TabTrans;        
   fltarray WP_TabNoise;    // Noise level per band
      
   // Routines
   // --------
   MRBase ()  {reset();}   // Constructor
   void reset();          // Reset all internal variables
   void alloc();          // Allocate the memory

   void decomposition (); // Start the decomposition
                          // alloc MUST have been called before
                           
   void free();           // deallocation of the class

   void init_deconv(Ifloat & ImaPsf, Ifloat &NoiseData, float Eps=0.001);
   void init_deconv(Ifloat & ImaPSF, float SigmaNoise, float Eps=0.001, int InitRnd=100);
   // Divide the data by the PSF in the Fourier space, and estimate
   // the noise detection level in the different transforms.

   ~MRBase(){free(); reset();}
};

/***********************************************************************/

#endif
