/******************************************************************************
**                   Copyright (C) 2002 by CEA
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

#ifndef _PCURVELET_H_
#define _PCURVELET_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter.h"
#include "Ridgelet.h"
#include "MR_Obj.h"
#include "Curvelet.h"

/***********************************************************************/


class PyrCurvelet {
   int NlIma,NcIma;     // input image size
   int UserBlockSize;   // First scale chosen Block size
   
   MultiResol MR_Data;     // 2D WT used in the Curvelet Transform
   int NbrBand2D;          // Number of bands in the 2D WT  
    
   Ifloat *TabRidTrans;    // Array of ridgelet transform
   Ridgelet *TabRidgelet;  // Array of ridgelet Class
                           // One ridgelet transform per scale
   int SizeTabRid;         // Size of the TabRidgelet Class table
   
   SubBand1D *Ptr_SB1D; // pointer to the filter bank to use
                        // in case of othogonal 1D WT

   int NbrBandCur;      // Number of bands in the curvelet transform
   
   void ridgelet_init(Ridgelet &Rid);
   // Initialize the Ridgelet transform with curvelet parameters
   
   Bool AllocClass; // True if the CLASS has been allocated
                    
   void mr_io_fill_header(fitsfile *fptr);

   void tab_block_size_init(int BS);
                         // Initialise TabBlockSize starting with 
                         // a block size = BS at the finest scale
   
   fltarray *TabCurSigma; // Noise Standard deviation in a given band
                          // for a realization of the noise 
                        
   public:
   type_transform TransformWT2D;  // type of 2D Wavelet transform. 

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
   

   void reset(SubBand1D *SB1D=NULL); // Reset all parameters to default

   PyrCurvelet(SubBand1D &SB1D) { reset(&SB1D);} // Constructor
   PyrCurvelet() { reset();}

   void alloc (int Nl, int Nc, int BS);
   // initialize the CLASS for a given image size and 
   // a given block size.
   
   float norm_band(int s2d, int s1d); 
   // Normalisation parameter to be applied to the noise for a 
   // a given scale s2d of the a trous algorithm, and a given scale
   // s1d of the ridgelet transform

   void transform(Ifloat &Image);
   // Apply the curvelet transform  

   void recons(Ifloat &Image);
   // Reconstruct an image from its curvelet transform
 
   //void filtering(Ifloat &Image, Ifloat & Filter, float NoiseIma, float N_Sigma);
   // Image filtering by the curvelet transform

   int nbr_band() const{return NbrBandCur;}
   // return the number of band of the curvelet transform
   //         NbrBand = total( TabRidgelet[0:NbrBand2D-2].NbrScale ) + 1
      
   int nbr_rid_scale(int s2d) 
     {return (s2d == NbrScale2D-1) ? 1: TabRidgelet[s2d].NbrScale;}
   // return the number of band at a give scale 
   
   int ipos(int s2d, int s1d) 
      const{return (s2d == NbrScale2D-1) ? 0 : TabRidgelet[s2d].ipos(s1d);}
   int jpos(int s2d, int s1d)
      const{return (s2d == NbrScale2D-1) ? 0 : TabRidgelet[s2d].jpos(s1d);}
   int size_nl(int s2d, int s1d) 
      const{return (s2d == NbrScale2D-1) ? MR_Data.size_band_nl(NbrBand2D-1) : TabRidgelet[s2d].size_scale_nl(s1d);}
   int size_nc(int s2d, int s1d) 
      const{return (s2d == NbrScale2D-1) ? MR_Data.size_band_nc(NbrBand2D-1) : TabRidgelet[s2d].size_scale_nc(s1d);}
   int size_band_nl(int b);
   int size_band_nc(int b);
   
   void get_band(int s2d, int s1d, Ifloat &Band);
   void put_band(int s2d, int s1d, Ifloat &Band);
   void get_band(int NumBand, Ifloat &Band);
   void put_band(int NumBand, Ifloat &Band);
   void get_scale_number(int NumBand, int & s2d, int & s1d);
   int block_size(int NumBand);

   void get_norm_coeff(float N_Sigma=3.);
   // Simulate a Gaussian noise with sigma=1, apply the curvelet trans.
   // and call set_noise_level to store the detection thresholds in
   // TabCurSigma.
   
   void get_norm_coeff(Ifloat &ImaNoise, float N_Sigma);
   // Apply the curvelet transform to ImaNoise  
   // and call set_noise_level to store the detection thresholds in
   // TabCurSigma.

   void set_noise_level(float N_Sigma=3.);
   // Estimate the detection levels from the transformed coefficients.
   // The transform coefficients should be the coeff relative to a
   // noise distribution.
   
   float norm(int b, int i, int j);
   
   // void read(char *Name, fltarray &TabBand);
   // Read a curvelet transform and initialize the class
   // void write(char *Name, fltarray &TabBand);
   // write the curvelet  transform in the FITS format

   void get_stat(fltarray &TabStat);
   // Calculates some statistics about the coefficients
   //           TabStat(b, 0) = (float) Sigma;
   //	        TabStat(b, 1) = (float) Skew;
   //	        TabStat(b, 2) = (float) Curt;
   //           TabStat(b, 3) = (float) Min;
   //           TabStat(b, 4) = (float) Max; 
   
   float & operator() (int NumBand, int i, int j);
 
   void threshold(float SigmaNoise, float N_Sigma=3.);
   
   void mix_filtering(Ifloat &Image);
   void wiener_filter(float SigmaNoise, int WienerBlockSize=7);
   void mix_threshold_one_block(int s2d, Ifloat &ImaBlock, Ifloat &RidIma, 
                               float Nsig, float SigmaNoise);
	  
   ~PyrCurvelet(){Ptr_SB1D = NULL;}
};

/***********************************************************************/
#endif
