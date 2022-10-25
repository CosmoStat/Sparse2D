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

#ifndef _MCA_H_
#define _MCA_H_

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Usage.h"
#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM_DCT.h"
#include "Ridgelet.h"
#include "Curvelet.h" 
#include "IM_Regul.h"
#include "IM_Rot.h"
#include "IM1D_IO.h"
#include "MDCT.h"
#include "WT_Mirror.h"
#include "WPackets.h"
#include "PCur.h"
#include "FCur.h"

// #include "FFT2D_CircleWT1d.h"

#define NBR_TOT_MCA 9
#define NBR_MCA 8
#define DEF_MCA_FIRST_SOFT_THRESHOLD 50.  // First soft thresholding level 
#define DEF_MCA_LAST_SOFT_THRESHOLD  3.  // Last soft thresholding level 
#define NBR_THRES_DEC 3 //Nb of type of decreasing threshold

#define DEF_MCA_NBR_ITER 1000             // Default number of iterations

enum type_mcabase {MCA_ATROU,MCA_WT,MCA_RID,MCA_CUR, MCA_COS,MCA_WP, MCA_PCUR, MCA_FCUR,
                   MCA_MCOS,MCA_PIXEL, MCA_DDWT, MCA_FFT2D_WT1D};

enum type_threshold_decrease {MCA_THRES_LINEAR, MCA_THRES_EXP, MCA_THRES_MOM};

/***********************************************************************/

inline const char *StringMCABase (type_mcabase  type)
{
    switch (type)
    {
        case MCA_ATROU:
	      return ("A trous algorithm"); 
        case MCA_WT: 
              return ("bi-orthogonal WT with 7/9 filters"); 
        case MCA_CUR: 
              return ("Curvelet transform 02"); 
        case MCA_PCUR: 
              return ("Curvelet transform 05 (Pyramidal Construction)"); 
	case MCA_FCUR: 
              return ("Curvelet transform 04 (Fast Transform)");
	case MCA_RID:
	      return ("Ridgelet transform");
 	case MCA_FFT2D_WT1D:
	      return ("FFT2D-WT1D-PowSpectrum");
	case MCA_DDWT:
	      return ("Multi-Direction WT");
	case MCA_COS:
	      return ("Local Discrete Cosinus Transform");
	case MCA_MCOS:
	      return ("Multiscale Local Discrete Cosinus Transform");
	case MCA_PIXEL:
	      return ("Pixel basis");
	case MCA_WP:
	      return ("Wavelet Packet basis");
	default: 
	      return ("Undefined transform");
     }
}

/***********************************************************************/

class MCA {
    // General parameters
   int Nl;                          // input image size
   int Nc;
   int NlRid;
   int NcRid;          // Ridgelet image size

   float Lambda;                    // Soft threshold parameter
   float StepL;                     // Decreasing step per iteration 
                                    // applied to the Lambda value

   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   HALF_DECIMATED_2D_WT *WT; // Pointer to the wavelet transform class
   int WT_NbrBand;           // Number of bands in the transform
   FilterAnaSynt *WT_SelectFilter; // Selected filter bank
   SubBandFilter *WT_SB1D;  // Pointer to the fiter bank decomposition class
   Bool *WT_TabDec;         // Table indicating scales to be decimated
   Ifloat *WT_TransCorrel;  // Wavelet bases decomposition used for correlation constraint

   // Parameter for the cosinus transform
   // -----------------------------------
   LOCAL_DCT2D LDCT; // Local cosinus transform class
   Ifloat COS_Resi;  // Cosinus bases decomposition
   LOCAL_DCT2D LDCT_Correl;  // Cosinus bases decomposition

   // Parameter for the multiscale cosinus transform
   // -----------------------------------
   MDCT M_DCT;

   ATROUS_2D_WT AWT;      // Atrou class
   Ifloat *AT_TransCorrel;      // atrous bases decomposition

    // Parameter for the curvelet 
   // -----------------------------------
   Ifloat  *CUR_Correl;   // Curvelet bases decomposition
   
    // Parameter for the pyramidal curvelet 
   // -----------------------------------
   PyrCurvelet PCUR_Correl;   // Curvelet bases decomposition
   
   // Parameter for the fast curvelet transform 
   // -----------------------------------
   FCUR *FCUR_Correl; // Pointer the Fast curvelet transform CLASS
    
    // Parameter for the Ridgelet Transform 
   // -----------------------------------
   Ifloat RID_TransCorrel;         // Ridgelet decomposition

    // Parameter for the Packets WT
   // -----------------------------------
   WPACKETS_2D *WP;         // Wavelet Packet Class
   FilterAnaSynt *WP_SelectFilter; // Selected filter bank
   SubBandFilter *WP_SB1D;  // Pointer to the fiter bank decomposition class
   int WP_NbrBand;           // Number of bands in the transform
   WPACKETS_2D *WPC;               // WP class used for correlation constraint
   WPTransf_2D* WP_TabTransCorrel; // WP transform used for correlation constraint
   
   LineCol LC;  // Directional WT
   FilterAnaSynt *DWT_SelectFilter; // Selected filter bank
   SubBandFilter *DWT_SB1D;  // Pointer to the fiter bank decomposition class

   // CircleWT Fft2dWt1d;   // FFT2D-WT1D around zero on power spectrum class

   float atrou_max(Ifloat & ImaAtrou);
    float atrou_init(Ifloat & ImaAtrou, Bool SetMAD);
    void atrou_proj(Ifloat & ImaAtrou, float Lamba_Sigma);
   // Projection on the a trous transform

   void rid_proj(Ifloat & ImaRid, float Lamba_Sigma);
   // Projection on the ridgelet transform

   float wt_max(Ifloat & ImaWT);
   float wt_init(Ifloat & ImaWT, Bool SetMAD);
   void wt_proj(Ifloat & ImaRec,  float Lamba_Sigma);
   // Projection on the orthogonal WT

  float wp_max(Ifloat & ImaWP);
  float wp_init(Ifloat & ImaWP, Bool SetMAD);
   void wp_proj(Ifloat & ImaRec,  float Lamba_Sigma);
   // Projection on the orthogonal Wavelet Packets

   float cos_max(Ifloat & ImaCos);
   float cos_init(Ifloat & ImaCos);
   void cos_proj(Ifloat & ImaRec, float Lamba_Sigma);
   // Projection on the cosinus transform

   void mcos_proj(Ifloat & ImaRec, float Lamba_Sigma);
   // Projection on the multiscale cosinus transform

   float cur_max(Ifloat & ImaCur);
   float cur_init(Ifloat & ImaCur, Bool SetMAD);
   void cur_proj(Ifloat & ImaCur,  float Lamba_Sigma);
   // Projection on the curvelet transform

   float pcur_max(Ifloat & ImaCur);
   float pcur_init(Ifloat & ImaCur, Bool SetMAD);
   void pcur_proj(Ifloat & ImaCur, float Lamba_Sigma);
   // Projection on the pyramidal curvelet transform

   float fcur_max(Ifloat & ImaFCur);
   float fcur_init(Ifloat & ImaCur, Bool SetMAD);
   void fcur_proj(Ifloat & ImaCur, float Lamba_Sigma);
   // Projection on the fast curvelet transform

   void dwt_proj(Ifloat & Sol, float Lamba_Sigma, float Angle);
   // Projection on the directional WT

   // void fft2dwt1d_proj(Ifloat & ImaPS,  float Lamba_Sigma);
   // Projection on the curvelet transform

   float pixel_max(Ifloat & ImaPix);
   float pixel_init(Ifloat & ImaPix, Bool SetMAD);
   void pixel_proj(Ifloat & ImaPix, float Lamba_Sigma);
   float PIX_TabNoise;
   // Projection on the pixel basis
		     
   inline float mca_update(float CoefSol, float Threshold, float Noise);
   // Calculate the update coefficient of the transform, taking into acount:
   // CoefSol: IN = coefficient of the solution
   // Threshold: IN =  Regularization parameter 

   Ifloat ImaCorrel;  // Used if the correlation constraint is True
   public:
   void init_decomposition();

   // ***************** PUBLIC MCA PART  *********************
	 Bool Bounded;   //Min-Max Bounded
   Bool SigmaBounded;   //Sigma Bounded
   Bool WriteAll;   // Write intermediate results on the disk (every 10 iterations)
   Bool UseZoom;    // Zoom mode
   int ZoomFactor; // Zoom Factor
   Ifloat UnZoomData; // Original  Data (without interpolation)
   int TV_NbrScale2D; 
   Bool Linear;  
   type_threshold_decrease  T_Decreasing;
   
   // Bool TabBase[NBR_MCA]; // TabBase[i] = True when bases number i is selected
   type_mcabase TabSelect[NBR_MCA]; // selected bases
   int NbrBase;                     // number of selected bases                    

   Ifloat Data;             // Input data
   Ifloat Resi;    
   Ifloat TabImaRec[NBR_MCA]; // reconstructed image for each bases
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
   float LambdaTV;      // Total Variation regularization parameter
   Bool UseNormL1;      // If true, the L1 norm is used instead of the L0
   float NSigmaSoft;
   Bool  UseCorrelConstraint; // if true use the correlation constraint for the separation
   float CorrelConstraint;    // correlation constraint parameter
   
   // Parameter for the atrous algorithm
   // ----------------------------------
   Ifloat *AT_Trans;      // atrous bases decomposition
   Ifloat *AT_Trans1;      // atrous bases decomposition
   Ifloat *AT_Trans2;      // atrous bases decomposition
   int AT_NbrScale2D;     // Number of scales for the atrous transform
   Bool AT_KillLastScale; // do not used the last scale for the atrous 
                          // decomposition
   Bool AT_PositivRecIma; // If true, impose positive constraint
                          // to the reconstructed image from
                          // the atrou algorithm
   fltarray AT_TabNoise;    // Noise level per band

    // Parameter for the Ridgelet transform
   // -------------------------------------
   Ridgelet RID;             // Ridgelet transform class
   Ifloat RID_Trans;         // Ridgelet decomposition
   int RID_BlockSize;        // Ridgelet block size
   int RID_FirstDetectScale; // First detection scale in the ridgelet transform
   Bool RID_PositivRecIma;   // If true, impose positive constraint  

   // Parameter for the curvelet transform
   // ------------------------------------
   Curvelet CUR;            // Curvelet class
   int  CUR_BlockSize;      // block size at the first scale
   Bool CUR_BlockOverlap;   // block overlapping
   Bool CUR_PositivRecIma;  // If true, impose positive constraint
   int  CUR_NbrScale2D;     // Number of scales in the a trou algorithm
   Bool CUR_KillLastScale;  // If true, last scale is killed
   Ifloat  *CUR_Trans;      // Curvelet bases decomposition
   fltarray CUR_TabNoise;   // Noise level per band

   // Parameter for the Fast curvelet transform
   // ------------------------------------
   FCUR *FCur;            // Curvelet class
   int  FCUR_NDIR;      // block size at the first scale
   Bool FCUR_PositivRecIma;  // If true, impose positive constraint
   int  FCUR_NbrScale2D;     // Number of scales in the a trou algorithm
   Bool FCUR_KillLastScale;  // If true, last scale is killed
   fltarray FCUR_TabNoise;   // Noise level per band
   
   // Parameter for the pyramidal curvelet transform
   // ------------------------------------
   PyrCurvelet PCUR;            // Curvelet class
   int  PCUR_BlockSize;      // block size at the first scale
   Bool PCUR_BlockOverlap;   // block overlapping
   Bool PCUR_PositivRecIma;  // If true, impose positive constraint
   int  PCUR_NbrScale2D;     // Number of scales in the a trou algorithm
   Bool PCUR_KillLastScale;  // If true, last scale is killed
   fltarray PCUR_TabNoise;   // Noise level per band
   
   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   int WT_NbrUndecimatedScale; // Number of undecimated scale
   int WT_NbrScale2D;          // Number of scales
   Bool WT_PositivRecIma;      // If true, impose positive constraint
   Ifloat *WT_Trans;           // Wavelet bases decomposition
   fltarray WT_TabNoise;    // Noise level per band

   // Parameter for the cosinus transform
   // -----------------------------------
   Ifloat COS_Trans;          // Cosinus bases decomposition
   Bool  COS_PositivRecIma;    // If true, impose positive constraint
   int   COS_BlockSize;
   Bool  COS_Overlapping;
   float COS_Sensibility;
   float COS_Min;
   Bool  COS_WeightFirst;
   Bool  COS_Isotrop;  

      // Parameter for the multiscale cosinus transform
   // -----------------------------------
   Bool MCOS_PositivRecIma;    // If true, impose positive constraint
   int  MCOS_FirstBlockSize;
   Bool MCOS_Overlapping;
   int MCOS_NbrScale2D; 

   // Parameter for the directional WT
   // -----------------------------------
   Ifloat *DWT_TabImaRec; // Reconstructed images
   int DWT_NbrAngle;
   fltarray DWT_TabAngle;
   int DWT_NbrScale2Di;
   int DWT_NbrScale2Dj;
   int DWT_PositivRecIma;
   Ifloat DWT_Trans;
   Bool UseMad;
   Bool RemoveLastScale;
   Ifloat LastScale;
   Bool UseMask; // if true, some data are missing
   Ifloat MaskedData; // MaskedData(x,y) = 0 ==> missing data
                      // MaskedData(x,y) = 1 ==> good data

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
   MCA ()  {reset();}   // Constructor
   void reset();          // Reset all internal variables
   void alloc(int Nli, int Nci, Ifloat * ); // Allocate the memory

   void reconstruction(Ifloat &Result);
                          // Return the reconstruction image from the 
                          // decompositions

   void make_residual(Ifloat &Data, Ifloat &Resi);
   void decomposition (Ifloat &Result); 
                          // Start the decomposition
                          // alloc MUST have been called before
   void write_allima(char *FileName, fitsstruct &Header);
   void write_allima(char *FileName, int It);
   void free();           // deallocation of the class

   ~MCA(){free(); reset();}
};

/***********************************************************************/

#endif
