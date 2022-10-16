/******************************************************************************
**                   Copyright (C) 2000 by CEA + Santford University
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Hubert Druesne, Philippe Querre
**
**    Date:  02/06/2003 
**    
**    File:  cb_mca1d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : 
**    ----------- 
******************************************************************************/

#ifndef _MCA1D_H_
#define _MCA1D_H_


//const float DEF_MCA1D_FIRST_SOFT_THRESHOLD = 50.;// First soft thresholding level 
const float DEF_MCA1D_LAST_SOFT_THRESHOLD  = 3.; // Last soft thresholding level 

const int DEF_MCA1D_NBR_ITER = 30;             // Default number of iterations
const int MAX_MGA_SCALE_1D = 100;
const int MAX_SMOOTH_REMOVE = 3;

enum type_mca1dbase { MCA1D_ATROU=0,
                      MCA1D_WT=1,
                      MCA1D_COS=2,
                      MCA1D_MCOS=3,
		      MCA1D_ALDCT=4};		   
const int  NBR_MCA1D = 5;

/***********************************************************************/

inline char* StringMCA1DBase (type_mca1dbase  type) {

    switch (type)
    {
        case MCA1D_ATROU:
	      return ("A trous algorithm"); break;
        case MCA1D_WT: 
              return ("Half decimated WT");  break;
	case MCA1D_COS:
	      return ("Local Discrete Cosinus Transform"); break;
	case MCA1D_MCOS:
	      return ("Multiscale Local Discrete Cosinus Transform"); break;
	case MCA1D_ALDCT:
	      return ("Adaptive Local Discrete Cosinus Transform"); break;
	default: 
	      return ("Undefined transform"); break;
     }
}
   
   extern float Tab1DSignifLevel[MAX_MGA_SCALE_1D];
   //float mr1d_detect_noise (fltarray& Signal, int NbrScale1D=5, 
   //                         float FirstSoftThreshold=3, int Niter=10, 
   //                         int NiterClip=3);
                            
/***********************************************************************/

class MCA1D {

   // General parameters
   float Lambda;                    // Soft threshold parameter
   float StepL;                     // Decreasing step per iteration 
                                    // applied to the Lambda value
   
   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   HALF_DECIMATED_1D_WT *WT; // Pointer to the wavelet transform class
   int WT_NbrBand;           // Number of bands in the transform
   FilterAnaSynt *WT_SelectFilter; // Selected filter bank
   SubBandFilter *WT_SB1D;  // Pointer to the fiter bank decomposition class
   Bool *WT_TabDec;         // Table indicating scales to be decimated

   // Parameter for the cosinus transform
   // -----------------------------------
   LOCAL_DCT1D LDCT; // Local cosinus transform class
   fltarray COS_Resi;  // Cosinus bases decomposition
   
   // Parameter for the multiscale cosinus transform
   // -----------------------------------
   MDCT1D M_DCT;

   // Parameter for the adaptated cosinus transform
   // -----------------------------------
   ALDCT AL_DCT;

   // Parameter for the Atrou transform
   // -----------------------------------   
   ATROUS_1D_WT AWT;      

   void atrou_proj(fltarray& SigRec, float Lamba_Sigma, int Iter);
   // Projection on the a trous transform   

   void wt_proj(fltarray& SigRec,  float Lamba_Sigma, int Iter);
   // Projection on the orthogonal WT

   void cos_proj (fltarray& SigRec, float Lamba_Sigma, int Iter);
   // Projection on the cosinus transform

   void mcos_proj (fltarray& SigRec, float Lamba_Sigma, int Iter);
   // Projection on the multiscale cosinus transform
		     
   void MCA1D::aldct_proj(fltarray& SigRec, float Lamba_Sigma, int Iter);
   // Projection on the adaptated cosinus transform

   inline float mca1d_update (float CoefSol, float Threshold, float Noise);
   // Calculate the update coefficient of the transform, taking into acount:
   // CoefSol: IN = coefficient of the solution
   // Threshold: IN =  Regularization parameter 

   inline bool mca1d_is_on_support (float CoefSol, float Threshold, float Noise);
   
   public:
   //void init_decomposition();
   //void threshold (fltarray& Signal, float T=0);
   float compute_first_lambda (fltarray& Signal);
   void  remove_smooth_plane(fltarray& Signal, int NbSmoothRemove);

   // ***************** PUBLIC MCA1D PART  *********************
  
   // Bool TabBase[NBR_MCA1D]; // TabBase[i] = True when bases number i is selected
   type_mca1dbase TabSelect[NBR_MCA1D]; // selected bases
   int NbrBase;                     // number of selected bases                    
   int Nx;                          // input signal size
   fltarray Smooth;			    

   fltarray Data;             // Input data
   fltarray Resi;    
   fltarray TabSigRec[NBR_MCA1D]; // reconstructed signal for each bases
   Bool PositivRecSig;   // If true, impose positive constraint
                         // to the finale reconstructed signal
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
   
   // Parameter for the atrous algorithm
   // ----------------------------------
   int AT_FirstDetectScale;
   fltarray *AT_Trans;      // atrous bases decomposition
   int AT_NbrScale1D;     // Number of scales for the atrous transform
   Bool AT_KillLastScale; // do not used the last scale for the atrous 
                          // decomposition
   Bool AT_PositivRecSig; // If true, impose positive constraint
                          // to the reconstructed signal from
                          // the atrou algorithm

   // Parameter for the orthogonal wavelet transform
   //-----------------------------------------------
   int WT_NbrUndecimatedScale; // Number of undecimated scale
   int WT_NbrScale1D;          // Number of scales
   Bool WT_PositivRecSig;      // If true, impose positive constraint
   fltarray *WT_Trans;           // Wavelet bases decomposition

   // Parameter for the cosinus transform
   // -----------------------------------
   fltarray COS_Trans;          // Cosinus bases decomposition
   Bool COS_PositivRecSig;      // If true, impose positive constraint
   int  COS_BlockSize;
   Bool COS_Overlapping;
   float COS_Sensibility;
   float COSMin;
   Bool CosWeightFirst;
   
   // Parameter for the multiscale cosinus transform
   // -----------------------------------
   Bool MCOS_PositivRecSig;    // If true, impose positive constraint
   int  MCOS_FirstBlockSize;
   Bool MCOS_Overlapping;
   int MCOS_NbrScale1D;
   int MCOS_FirstDetectScale;
   
   // Parameter for the adaptated cosinus transform
   // -----------------------------------
   int ALDCT_NbrScale1D;
   fltarray *TabALDCT;      // ALDCT bases decomposition
   float ALDCT_Sensibility;
   int ALDCT_InfoCost;

   // Routines
   // --------
   MCA1D ()  {reset();}   // Constructor
   void reset();          // Reset all internal variables
   void alloc(int Nx, fltarray* Tab); // Allocate the memory

   void reconstruction(fltarray& Result);
                          // Return the reconstruction signal from the 
                          // decompositions

   void make_residual(fltarray& Data, fltarray& Resi);
   void decomposition (fltarray& Result); 
                          // Start the decomposition
                          // alloc MUST have been called before
   void write_allima(char *FileName, fitsstruct &Header);
   void write_allima(char *FileName, int It);
   void free();           // deallocation of the class
   inline void infoVerbose ();

   ~MCA1D(){free(); reset();}
};


#endif
