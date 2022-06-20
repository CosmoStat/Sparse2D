/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/01/20
**    
**    File:  MW1D_Filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  spectrum filtering
**    ----------- 
**                 
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
 
#include "MR1D_NoiseModel.h"

#include "CErf.h"
#include "CMem.h"
#include "MW1D_Filter.h"

/*********************************************************************/
 
void mw1d_fm1 (MR_1D &MR_Data, 
	     MR1DNoiseModel & NoiseModel, 
	     fltarray & TabAlpha, CMemWave & MemWObj, Bool UseDataSNR,
             MR_1D *MR_Model)	    
/* Filters the wavelet coefficient using the multiscale entropy method:
   MR_Data : in out = wavelet coefficients
   NoiseModel: in = data noise modeling
   TabAlpha: in = regularization parameter (one per scale)
   CMemWave: in = multiscale entropy class
   MR_Model: in = Model of the solution (= NULL if no model)
   UseDataSNR: in = True if the regularization is varying with the SNR
*/
{
    int s,i;
    int NbrScale = MR_Data.nbr_scale();
    float Model;

    for (s=0; s < NbrScale-1; s++)
    for (i=0; i <  MR_Data.size_scale_np(s); i++)
    {
       float Sigma = NoiseModel.sigma(s,i);
       float Coef = MR_Data(s,i);
       float RegulParam = TabAlpha(s);
       
       // use the input model for the wavelet coefficients
       if (MR_Model == NULL) Model = 0.;
       else  Model = (*MR_Model)(s,i);
	  
       // Modify Alpha using the SNR of the data
       if (UseDataSNR == True)
       {
	  // SNR = MAX [ ABS(data) / (N*Sigma), 1 ]
    
          float AlphaP = ABS(Coef / (Sigma* NoiseModel.NSigma[s]));
          if (AlphaP > 1) AlphaP = 1;
	     
	  // if AlphaP = 0, J = h_n ==> Solution = model
	  // when RegulParam < 0, the filtering routine set to
	  // the model the solution
          if (AlphaP < FLOAT_EPSILON) RegulParam = -1;
	  else RegulParam *= (1-AlphaP)/AlphaP;
       }
       MR_Data(s,i) = MemWObj.filter(Coef,  RegulParam, Sigma);
    }   
}

/*********************************************************************/

void mw1d_filter (MR_1D & MR_Data, MR1DNoiseModel & NoiseModel,  
              fltarray &TabAlpha, float Alpha, Bool UseDataSNR,
	    float ConvgParam, int MaxIter,  Bool Verbose, MR_1D *MR_Model)
	    
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
{
   CMemWave MemWObj;
   int NbrScale = NoiseModel.nbr_scale();
   int Nbr_Band = NbrScale-1;
   int Np = MR_Data.size_ima_np();
   MR_1D MR_Sol (Np, MR_Data.Type_Transform, "MR_Data", NbrScale);

   fltarray Data(Np);
   fltarray  Ima(Np);
   double  SigmaNoise;
   int b,i,Iter=0;
   fltarray  RegulMin( Nbr_Band ); 
   fltarray  RegulMax( Nbr_Band );    
   fltarray  TabDelta( Nbr_Band );    
   
   // we need the raw data in order to be able to calculate the residual
   // image: resi = data - solution
   MR_Data.recons(Data);
   
   // initialization
   for (i = 0; i <  MR_Sol.size_scale_np(Nbr_Band); i++) 
                             MR_Sol(Nbr_Band,i) = MR_Data (NbrScale-1,i);
 
  // initialization
   for (b = 0; b < Nbr_Band; b++)
   {
      RegulMin(b) = 0;                 // minimum alpha value
      RegulMax(b) = DEF_MAX_MEM1D_ALPHA; // maximum alpha value
      if (b == 0) TabAlpha(b)= 5;
      else TabAlpha(b)= 1.;
      TabDelta(b) = 0.;
   }   
   
   do 
   {
      Iter ++;
      if (Verbose == True) cout << "Iter " << Iter << endl;
      
      // copy the raw wavelet coef. into the solution array
      for (b = 0; b < Nbr_Band; b++) 
      for (i = 0; i < MR_Sol.size_scale_np(b); i++) MR_Sol(b,i) = MR_Data(b,i);
      
      // RegulParam = alpha parameter new calculation
      for (b = 0; b <  Nbr_Band; b++) 
                 TabAlpha(b) = (RegulMin(b) + RegulMax(b))/2;
      
      // correct the wavelet coefficients
      mw1d_fm1(MR_Sol, NoiseModel,TabAlpha, MemWObj, False, MR_Model);
      
      // new estimation of the Alpha paramters
      for (b = 0; b <  Nbr_Band; b++) 
      {
         SigmaNoise = 0.;
	 float T=0.;
         for (i = 0; i <  MR_Sol.size_scale_np(b); i++)
         {
	    float Coef = MR_Data(b,i);
	    float SigmaCoef = NoiseModel.sigma(b,i);
	    float Resi = Coef - MR_Sol(b,i);
	    
           float AlphaP = ABS(Coef / (SigmaCoef*3));
           if (AlphaP > 1) AlphaP = 0;
	   else AlphaP = 1 - AlphaP;
	   AlphaP=1.;
   	   T +=  AlphaP;
	  
	    float ExpectVariance =  SigmaCoef*SigmaCoef;
	    if (ExpectVariance > FLOAT_EPSILON)
 	          SigmaNoise += AlphaP*(Resi*Resi) / ExpectVariance;
            else SigmaNoise += AlphaP*(1 + Resi*Resi);
	 }
	 SigmaNoise = sqrt(SigmaNoise/T);
	 //if (SigmaNoise <  SigmaNoise*(1- Tolerance))  RegulMin(b) = TabAlpha(b);
	 //else RegulMax (b) = TabAlpha(b);
	 
	 //if (SigmaNoise >= 1 + (1-Tolerance)*SigmaNoise) 
	 //                           RegulMax (b) =  TabAlpha(b);
	 //else if (SigmaNoise < 1 - Tolerance*SigmaNoise	)  RegulMin(b) = TabAlpha(b); 		    

	 if (SigmaNoise >= 1)  RegulMax (b) =  TabAlpha(b);
	 else  RegulMin(b) = TabAlpha(b);      
  
	 TabDelta(b) = ABS(SigmaNoise - 1.);
 	 if (Verbose == True) 
	     cout << "   band " << b+1 << " Normalized sigma "   << SigmaNoise
	           << " Alpha = " <<  TabAlpha(b) << "  Delta = " << TabDelta(b)
		   << endl;      
      }
     // convergence test
   } while( (TabDelta.max() > ConvgParam) && (Iter < MaxIter));

   // the calculate alpha is 
   for (b = 0; b <  Nbr_Band; b++) 
                 TabAlpha(b) = (RegulMin(b) + RegulMax(b))/2*Alpha;    
   if (Verbose == True)
     for (b = 0; b <  Nbr_Band; b++) 
      cout << "band " << b+1 << " Optimal Regul. Param = " <<  TabAlpha(b) << endl;
    
   // final correction of the wavelet coefficients
   mw1d_fm1(MR_Data, NoiseModel, TabAlpha, MemWObj,  UseDataSNR, MR_Model);

}

/***************************************************************************/
      
 
 
