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
**    File:  Curvelet.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _CURVELET_H_
#define _CURVELET_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter.h"
#include "Ridgelet.h"
#include "MR_Obj.h"
#include "IM_Prob.h"

#define NBR_CUR_TYPE_BLOCK 5
enum type_curvelet_block {BL_CONST,BL_UP,BL_UP2,BL_DOWN, BL_DOWN2};
// the block size can be modified at each scale of the wavelet transform
// BL_CONST ==> the block size is always the same
// BL_UP ==> the block size is doubled at each resolution level
//           (from high frequency to low freq)
// BL_UP2 ==> the block size is doubled at each two resolutions
// BL_DOWN ==> the block size is divided by 2 at each resolution level
// BL_DOWN2==> the block size is divided by 2at each two resolutions 

#define DEF_CUR_NBR_SCALE_2D 4 // Default number of scale in the a trous algo
#define DEF_CUR_BLOCK_SIZE 16  // Default block size at the first scale
#define DEF_TYPE_CUR_BLOCK BL_UP2
const char * StringBlockType (type_curvelet_block type);

/***********************************************************************/


class Curvelet {
   SubBand1D *Ptr_SB1D; // pointer to the filter bank to use
                        // in case of othogonal 1D WT

   int NlIma,NcIma;     // input image size
   int CurNl,CurNc;     // Curvelet band size
   int UserBlockSize;   // First scale chosen Block size
   void set_param(int Nl, int Nc); // Set the curvelet transform size from the 
                                   // input image size
 
   void set_param_from_trans(int Nlc, int Ncc); 
                                  // Set the image size from the 
                                 // the curvelet transform image size

   void ridgelet_init(Ridgelet &Rid);
   // Initialize the Ridgelet transform with curvelet parameters

   Ridgelet *TabRidgelet; // One ridgelet transform per scale
   
   Bool AllocClass; // True if the CLASS has been allocated
   int SizeTabRid;  // Size of the TabRidgelet table
                    
   void mr_io_fill_header(fitsfile *fptr);

   void tab_block_size_init(int BS);
                         // Initialise TabBlockSize starting with 
                         // a block size = BS at the finest scale
   
   fltarray TabCurSigma; // Noise Standard deviation in a given band
                         // for a realization of the noise 
                        
   CImaProb *TabCP; // For noise modeling, TabCP[b] is pdf modelisation of the noise
   fltarray TabMaxDetect;  // Detection level per band for positive coefficient
   fltarray TabMinDetect;  // Detection level per band for negative coefficient
   fltarray TabSigmaNoise;  // noise level per scale
   fltarray TabNSigma;  // Nsigma detection level

   public:
   inline float & max_level(int b) { return TabMaxDetect(b);}
   inline float & min_level(int b) { return TabMinDetect(b);}
   inline float & nsigma(int b) { return TabNSigma(b);}
   inline float & sigma_noise(int b) { return TabSigmaNoise(b);}
   void print_info_noise();

   Bool  VarNorm;      // If true, the curvelet coefficients are normalized by sqrt(BLOCK_SIZE)
   Bool MSVST;         // if True, the Multiscale variance stabilization is applied during
                       // the transform and the reconstruction.
		       // Only used for Poisson noise denoising.

   type_format FormatInputImag;   // data format of the input image
   type_curvelet_block TypeBlock;  // Type of curvelet blocks
   type_ridgelet_WTtrans RidTrans; // Ridgelet transform type

   Bool BlockOverlap;  // If True, Overlapped blocks are used.
   Bool ColTrans;      // If true, apply an orthogonal 1DWT on column
   Bool GetAutoNbScale; // If true the number of scale is automatically
                        // calculated by:
                        // NbrScale = fix( log( (3N/4) / log(2)));
  
   Bool Verbose;        // Verbose mode  (default is no).

   Bool StatInfo;       // If true, statistical information of each band 
                        // is printed (default is no).

   type_border Border;  // Border used in the 2D a trous WT

   Bool WPTrans;        // Wavelet packets are applied on the first scale
                        // of the ridgelet transform

   int NbrScale2D;      // Number of scales in the 2D a trous WT

   int NbrScale1D;     // Number of scales used by the 1D wavelet transform
                       // if NbrScale = 0 then Ridgelet transform = Radon transform
                       // if NbrScale < 0 the number of scales is automatically
                       // calculated by:
                       // NbrScale = fix( log( (3N/4) / log(2)));

   intarray TabBlockSize; // Block size used in the Curvelet transform
                          // TabBlockSize(s) = Block size of scale s
                          // s = 0..MAX_SCALE-1

   Bool OddBlockSize;     // If true block size is an odd number
   Bool AngleNormalization; // if true, then renormalize the curvelet coeff
                           // so that the std of the coeff. for a given angle
			   // is equal to 1.

   void reset(SubBand1D *SB1D=NULL); // Reset all parameters to default

   Curvelet(SubBand1D &SB1D) { reset(&SB1D);} // Constructor
   Curvelet() { reset();}
   
   void alloc (int Nl, int Nc, int BS, Bool WeightBeforeTrans=False);
   // initialize the CLASS for a given image size and 
   // a given block size.
   
   float norm_band(int s2d, int s1d); 
   // Normalisation parameter to be applied to the noise for a 
   // a given scale s2d of the a trous algorithm, and a given scale
   // s1d of the ridgelet transform

   void transform(Ifloat &Image,  fltarray &TabBand);
   // Apply the curvelet transform and store the result in TabBand

   void recons(fltarray &TabBand, Ifloat &Image, Bool Init=True);
   // Reconstruct an image from its curvelet transform
   // If the Curvelet class has already been initialized with a given
   // input image size, then Init should be set to False

   void transform(Ifloat &Image,  Ifloat * & TabBand);
   // Apply the curvelet transform and store the result in TabBand
   // A each scale, the Ridgelet Class defined by TabRid is used
   // if TabRid == NULL, use a default one

   void recons(Ifloat *TabBand, Ifloat &Image);
   // Reconstruct an image from its curvelet transform
   // A each scale, the Ridgelet Class defined by TabRid is used
   // if TabRid == NULL, use a default one

   void filtering(Ifloat &Image, Ifloat & Filter, float NoiseIma, float N_Sigma);
   // Image filtering by the curvelet transform
   
   int nbr_rid_scale(int s2d) {return TabRidgelet[s2d].NbrScale;}
   int cur_nl() {return CurNl;}
   int cur_nc() {return CurNc;}
   int ipos(int s2d, int s1d) {return TabRidgelet[s2d].ipos(s1d);}
   int jpos(int s2d, int s1d){return TabRidgelet[s2d].jpos(s1d);}
   int size_nl(int s2d, int s1d) {return TabRidgelet[s2d].size_scale_nl(s1d);}
   int size_nc(int s2d, int s1d) {return TabRidgelet[s2d].size_scale_nc(s1d);}
   Ridgelet * get_ridgelet(int s2d) {return TabRidgelet+s2d;}
                    // return the ridgelet class used at a given scale.
   int nbr_band(); 
   // return the number of band of the curvelet transform
   //         NbrBand = total( TabRidgelet[0:NbrScale2D-2].NbrScale ) + 1

   void get_band(fltarray &TabBand, int s2d, int s1d, fltarray &Band);
   void put_band(fltarray &TabBand, int s2d, int s1d, fltarray &Band);
   void get_band(fltarray &TabBand, int NumBand, fltarray &Band);
   void put_band(fltarray &TabBand, int NumBand, fltarray &Band);
   void get_scale_number(int NumBand, int & s2d, int & s1d);
   void get_norm_coeff(int N, int ScaleNumber, int Ns1D, fltarray &Tab);
   
   void set_noise_level(Ifloat &NoiseData, float N_Sigma=3., Bool UseMad=False);
   // Calculate the detection level in the curvelet bands using a realization 
   // of the noise
   // NoiseData = noise realization
   // Output ==> TabMinDetect, TabMaxDetect and TabCP are initiatized
   
   void  set_noise_model_gaussian(fltarray & TabN_Sigma, float SigmaNoise);
   // Calculate the detection level in the curvelet bands using a Gaussian noise model
   // Output ==> TabMinDetect, TabMaxDetect and TabSigmaNoise are initialized
   
   void  set_noise_model_gaussian(float N_Sigma,  float SigmaNoise);
   // Calculate the detection level in the curvelet bands using a Gaussian noise model
   // Output ==> TabMinDetect, TabMaxDetect and TabSigmaNoise are initialized
   // The same value of NSigma is used in all bands except the fist scale where
   // NSigma = N_Sigma + 1 is used
   
   void set_noise_model_using_mad(Ifloat * &Trans, float N_Sigma);
   // Calculate the detection level in the curvelet bands using a Correlated Gaussian noise model
   // Output ==> TabMinDetect, TabMaxDetect and TabSigmaNoise are initialized
   // The MAD is used to estimate the standard deviation in the bands.
   
   double prob_noise(float Val, int b);
   // Return the probability that the value Val is due to the noise.

   void get_fdr_detect_level(Ifloat * &Trans, float N_Sigma, Bool ThresholdBand=False);
   // Estimate the detection level using the FDR (False Discovery Rate Method).
   //  Output ==> TabMinDetect, TabMaxDetect are calculated
   //  the routine set_noise_model_gaussian or set_noise_model_using_mad or set_noise_level
   //  must be called before.
   
   void set_noise_model_gaussian_in_band(int b, float N_Sigma, float SigmaNoise);
   // Calculate the detection level in a given curvelet band using a Gaussian noise model
   // Output ==> TabMinDetect(b), TabMaxDetect(b) and TanSigmaNoise(b) are initialized
   
   void get_band_no_border(fltarray &Trans, int s2d, int s1d, fltarray &Band);

   void read(char *Name, fltarray &TabBand);
   // Read a ridgelet transform and initialize the class
    
   void write(char *Name, fltarray &TabBand);
   // write the ridgelet transform in the FITS format

   ~Curvelet(){Ptr_SB1D = NULL; if (TabCP != NULL) delete [] TabCP;}
};

/***********************************************************************/
#endif
