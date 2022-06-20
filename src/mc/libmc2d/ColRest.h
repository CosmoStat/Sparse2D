
/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  may 2001
**    
**    File:  ColRest.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for color image restoration
**    ----------- 
**                 
******************************************************************************/

#ifndef _COLREST_H_
#define _COLREST_H_

class ColorRestore
{
  // Filter one band of the wavelet transform
  void col_filter_band (Ifloat & Band, float  Noise, float NSigma);

  public:
  
   // Contrast function, following K. Vande Velde (IEEE, 1999).
   // y(x) = (m/c)^p if |x| < c
   // y(x) = (m/|x|)^p if c <= |x| < m
   // y(x) = 1 if x >= m
   inline float contrast_function(float x)
   {
       float ValRet=0.;
       if (x > FLOAT_EPSILON)
       {
          if (ABS(x) < Contrast_C_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/Contrast_C_Param), Contrast_P_Param);
          else if (ABS(x) < Contrast_M_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/ABS(x)), Contrast_P_Param);
          else ValRet = 1.;
       }
       return ValRet;
   }
   
  ColorRestore() {HardThreshold=True;Contrast_C_Param=0;Contrast_M_Param=100;
                  Contrast_P_Param=0.3; Contrast_Q_Param = 0; MAD = False;
		  NbrUndec=-1;Verbose=False;ClipVal=3.;}
  

  // Parameter for the contrast function		       
  double Contrast_P_Param; // determine the degree of non-linearity
                           // p must be in ]0,1[
  double Contrast_Q_Param; // q must be in [-0.5,0.5] 
                           // q > 0, enhance less the darker part than
			   //        the lighter part
			   // q < 0 ==> enhance more the dark than lighter part
  double Contrast_M_Param; // Transform coefficients larger than m are not
                           // modified
  double Contrast_C_Param; // correspond to the noise level

  // Parameter for the sigma clipping
  float ClipVal;
  
  void retinex(fltarray &Data);
  // RETINEX contrast enhancement
  
  void multiscale_retinex(fltarray &Data);
  // Multiscale retinex contrast enhancement
  
  void equalize_histo(fltarray &Data, int BandNumber=2);
  // Equalize histogram in the HSV system
  // by default, the band number 2 (intensity) is used.  

  void equalize_histo_luminance(fltarray &Data);
  // Histogram Equalization of the luminance.
  
  // Parameter for the filtering of color image
  // float SigmaNoise;
  Bool Verbose;
  Bool HardThreshold;  // If true, a hard threshold is applied when filtering
                       // otherwise, use a soft-thresholding
  Bool MAD;            // If true, assume correlated noise and use
                       // the Mediane Absolute Deviation for the noise
		       // standard deviation estimation.
  int NbrUndec;        // Number of decimated scale in the WT
  
  void atrou_retinex(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima=0);
  // Contrast enhancement using a trous wavelet transform
  
  void col_filter(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima=0);
  // Color image filtering
  
  void info(fltarray &Data);
  
  void rescale_0_255(fltarray &Data);
  // Rescale the Data between 0 and 255
  
  void enhance(fltarray &Data, int Nbr_Plan, float N_Sigma, float Noise_Ima);
  // edge enhancement using the diadyc wavelet transform
  
  void sature_clipping(fltarray &Data);
  void sature_clipping(Ifloat &Data);
  // Apply a sigma clipping to the data
  
  // void atrou_clipping(fltarray &Data, int Nbr_Plan);
  
  void atrou_log(fltarray &Data, int Nbr_Plan, Bool Sign=False, 
                float Noise_Ima=0, float LogParam=0.1);
  // Apply the log-wavelet tranformation.
};

#endif
