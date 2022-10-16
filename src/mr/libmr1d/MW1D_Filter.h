/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  28/07/98 
**    
**    File:  MW1D_Filter.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#ifndef __MW1DFILTER__
#define __MW1DFILTER__

const float DEF_MAX_MEM1D_ALPHA = 200.0;
const int DEF_MEMW1D_FILTER_MAX_ITER = 10;
const float DEF_EPS_CVG_MEMW1D_FILTER = 1e-4;

/****************************************************************************/ 

#if _NOISE1D
void mw1d_fm1 (MR_1D &MR_Data, 
	     INoiseModel1d & NoiseModel, 
	     fltarray & TabAlpha, CMemWave & MemWObj, Bool UseDataSNR,
             MR_1D *MR_Model);
#else
void mw1d_fm1 (MR_1D &MR_Data, 
	     MR1DNoiseModel & NoiseModel, 
	     fltarray & TabAlpha, CMemWave & MemWObj, Bool UseDataSNR, 
             MR_1D *MR_Model);	    
#endif           
/* Filters the wavelet coefficient using the multiscale entropy method:
   MR_Data : in out = wavelet coefficients
   NoiseModel: in = data noise modeling
   TabAlpha: in = regularization parameter (one per scale)
   CMemWave: in = multiscale entropy class
   MR_Model: in = Model of the solution (= NULL if no model)
   UseDataSNR: in = True if the regularization is varying with the SNR
*/


/****************************************************************************/ 

#if _NOISE1D
void  mw1d_filter (MR_1D &MR_Data,  INoiseModel1d &  NoiseModel, 
	      fltarray &TabAlpha, float Alpha, Bool UseDataSNR,
	    float ConvgParam, int MaxIter,  Bool Verbose, MR_1D *MR_Model);
#else
void mw1d_filter (MR_1D & MR_Data, MR1DNoiseModel & NoiseModel,  
              fltarray &TabAlpha, float Alpha, Bool UseDataSNR,
	    float ConvgParam, int MaxIter,  Bool Verbose, MR_1D *MR_Model);
#endif
	    
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively at scale separately
   (by dichotomy) in order to have a residu compatible with the noise
   at all the scales
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
*/	
 
/****************************************************************************/ 

#endif




