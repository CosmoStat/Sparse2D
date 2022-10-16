/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Philippe Querre
**
**    Date:  96/05/02 
**    
**    File:  MR_mcFilter.h
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


#ifndef __MWMCFILTER__
#define  __MWMCFILTER__

#include "CErf.h"
#include "CMem.h"

const int DEF_MEM_ALPHA_CST = 1;
const int DEF_MEM_ALPHA_OPT = 2;
const int DEF_MEM_ALPHA_OPT_SCALE = 3;
const float DEF_MAX_MEM_ALPHA = 50.0;

const int DEF_MEM_FILT_MAX_ITER = 20;
const float DEF_MEM_FILT_CONGER = 0.01;

/****************************************************************************/ 

void mw_mcfilter(int CurrentBand, 
        MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        int TypeFilter, fltarray &TabAlpha, Bool DataSNR=False, 
	Bool DataModel=False, MultiResol *MR_Model=NULL,
	float ConvgParam = DEF_MEM_FILT_CONGER, int MaxIter=DEF_MEM_FILT_MAX_ITER, 
	Bool PositivConst=False, Bool Verbose=False);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TypeFilter = in: type of filtering
             = DEF_MEM_ALPHA_CST  ==> alpha is fixed 
	     = DEF_MEM_ALPHA_OPT  ==> alpha is globally optimized
	     = DEF_MEM_ALPHA_OPT_SCALE ==> alpha is optimized at each scale
	     
   TabAlpha = in-out: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   DataModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   DataSNR  = in: true if alpha parameter are modified using the SNR
                    of the data
   PositivConst = in: apply the positivy constraint

  Warning: works well, when each wavelet can be modelized by a Gaussian		    
*/	

/****************************************************************************/ 
 
void mw_mcfm1(int CurrentBand, 
            MultiResol & MR_Data, MRNoiseModel & NoiseModel, 
            fltarray &TabAlpha,  CMemWave & MemWObj, 
	    Bool UseModel,  MultiResol *MR_Model, Bool UseDataSNR);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
*/


/****************************************************************************/
		    	
void mw_mcfm3(int CurrentBand,
        MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool Verbose);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively at scale separately
   (by dichotomy) in order to have a residu compatible with the noise
   at all the scales
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
*/	

#endif




