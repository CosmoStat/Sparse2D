/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: %V%
**
**    Author: Philippe Querre
**
**    Date:  02/10/30 
**    
**    File:  MW_mcFilter.cc
**
*******************************************************************************
**
**    DESCRIPTION  Image filtering using the multiscale entropy
**    ----------- 
**                 
**
******************************************************************************/

#include "IM_Obj.h"
// #include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_NoiseModel.h"
#include "CErf.h"
#include "CMem.h"
#include "MW_mcFilter.h"



/****************************************************************************/ 
void mw_mcfilter(int CurrentBand, 
        MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        int TypeFilter, fltarray &TabAlpha, Bool DataSNR, 
	Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool PositivConst, Bool Verbose)
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

  Warning: works well, when each wavelet can be modelized by a Gaussian		    
*/	
{
   CMemWave MemWObj;

   switch (TypeFilter)
   {
      case DEF_MEM_ALPHA_CST:
        mw_mcfm1(CurrentBand, MR_Data, NoiseModel,TabAlpha, MemWObj, DataModel, MR_Model, DataSNR);
 	break;
      case DEF_MEM_ALPHA_OPT:
        cerr << "type of filtering nor used in this case ... " << endl;
	exit(-1);
	break;      
        //mw_fm2(CurrentBand, MR_Data, NoiseModel,  TabAlpha, MemWObj, DataSNR, 
	//          DataModel, MR_Model, ConvgParam, MaxIter, PositivConst, Verbose);
        //break;  
       case DEF_MEM_ALPHA_OPT_SCALE:
        mw_mcfm3(CurrentBand, MR_Data, NoiseModel,  TabAlpha, MemWObj, DataSNR, 
	          DataModel, MR_Model, ConvgParam, MaxIter, Verbose);
        break; 
      default: 
        cerr << "Error: unknow type of filtering ... " << endl;
	exit(-1);
	break;
   }
}
 

/****************************************************************************/
void mw_mcfm1(int CurrentBand,
            MultiResol & MR_Data, MRNoiseModel & NoiseModel, 
            fltarray &TabAlpha,  CMemWave & MemWObj, 
	    Bool UseModel,  MultiResol *MR_Model, Bool UseDataSNR)
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
    
{
   float AlphaP,Val;
   int i,j;
   float Model=0.;

      int s = MR_Data.band_to_scale(CurrentBand);
      
      // for all wavelet coefficients
      for (i = 0; i < MR_Data.size_band_nl(CurrentBand); i++)
      for (j = 0; j < MR_Data.size_band_nc(CurrentBand); j++)
      {
	  float Sigma = NoiseModel.sigma(CurrentBand,i,j);
          float RegulParam = TabAlpha(CurrentBand);
	  Val = MR_Data(CurrentBand,i,j);
	  
	  // use the input model for the wavelet coefficients
	  if ((UseModel != True) || (MR_Model == NULL)) Model = 0.;
	  else  Model = (*MR_Model)(CurrentBand,i,j);
 	  
	  // Modify Alpha using the SNR of the data
          if (UseDataSNR == True)
          {
	     // SNR = MAX [ ABS(data) / (N*Sigma), 1 ]
             AlphaP = ABS(Val / (Sigma* NoiseModel.NSigma[s]));
             if (AlphaP > 1) AlphaP = 1;
	     
	     // if AlphaP = 0, J = h_n ==> Solution = model
	     // when RegulParam < 0, the filtering routine set to
	     // the model the solution
             if (AlphaP < FLOAT_EPSILON) RegulParam = -1;
	     else RegulParam *= (1-AlphaP)/AlphaP;
          }
	  
	  // correct the wavelet coefficient
 	  MR_Data(CurrentBand, i,j) = MemWObj.filter(Val, RegulParam, Sigma, Model);
       }
}


/****************************************************************************/
void mw_mcfm3(int CurrentBand,
        MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool Verbose)
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
   PositivConst = in: apply the positivy constraint
*/	
{
   int Nbr_Band = NoiseModel.nbr_band()-1;
   int Nbr_Plan = NoiseModel.nbr_scale();
   int Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   type_transform Transform = NoiseModel.type_trans();
   MultiResol MR_Sol(Nl,Nc,Nbr_Plan,Transform,"MR_Sol");
   Ifloat Data(Nl,Nc, "data");
   Ifloat Ima(Nl,Nc, "Ima");
   double  SigmaNoise;
   int i,j, Iter=0;
   float  RegulMin; 
   float  RegulMax;    
   float  TabDelta;    
   float  AlphaUser = TabAlpha(0);
   //float  Alpha;

   // initialization
   MR_Sol.band(Nbr_Band) = MR_Data.band(Nbr_Band);

   RegulMin = 0;                 // minimum alpha value
   RegulMax = DEF_MAX_MEM_ALPHA; // maximum alpha value
   if (CurrentBand == 0) TabAlpha(CurrentBand) = 5;
   else TabAlpha(CurrentBand)= 1.;
   TabDelta = 0.;
   
   do 
   {
      Iter ++;
      if (Verbose == True) cout << "Iter " << Iter << endl;
      
      // copy the raw wavelet coef. into the solution array
     MR_Sol.band(CurrentBand) = MR_Data.band(CurrentBand);
      
      // RegulParam = alpha parameter new calculation
      TabAlpha(CurrentBand) = (RegulMin + RegulMax)/2;
      
      // correct the wavelet coefficients
      mw_mcfm1(CurrentBand, MR_Sol, NoiseModel,TabAlpha, MemWObj, 
               DataModel, MR_Model,DataSNR);
      
      // new estimation of the Alpha paramters
         int Nlb = MR_Data.size_band_nl(CurrentBand);
	 int Ncb = MR_Data.size_band_nc(CurrentBand);
         SigmaNoise = 0.;
         for (i = 0; i < Nlb; i++)
         for (j = 0; j < Ncb; j++)
         {
  	    float Resi = MR_Data(CurrentBand,i,j) - MR_Sol(CurrentBand,i,j);
 	    float SigmaCoef = NoiseModel.sigma(CurrentBand,i,j);
	    float ExpectVariance =  SigmaCoef*SigmaCoef;
	    if (ExpectVariance > FLOAT_EPSILON)
 	          SigmaNoise += (Resi*Resi) / ExpectVariance;
            else SigmaNoise += 1 + Resi*Resi;
 	 }
	 SigmaNoise = sqrt(SigmaNoise/(float) (Nlb*Ncb));
	 if (SigmaNoise >= 1)  RegulMax =  TabAlpha(CurrentBand);
	 else  RegulMin = TabAlpha(CurrentBand);      
  
	 TabDelta = ABS(SigmaNoise - 1.);
	 if (Verbose == True)
 	   cout << "   band " << CurrentBand+1 << " Normalized sigma " << 
	           SigmaNoise << " Alpha = " <<  TabAlpha(CurrentBand) << 
		   "  Delta = " << TabDelta << endl;      

   // convergence test
   } while( (TabDelta > ConvgParam) && (Iter < MaxIter));
   
   //***********
   //avant , tabDelta.max() > COnver... => PB ?????
   //***********

   // the calculate alpha is 
   TabAlpha(CurrentBand) = (RegulMin + RegulMax)/2*AlphaUser;    

   if (Verbose == True)
      cout << "band " << CurrentBand+1 << " Optimal Regul. Param = " <<  
              TabAlpha(CurrentBand) << endl;
    
   // final correction of the wavelet coefficients
   mw_mcfm1(CurrentBand, MR_Data, NoiseModel, TabAlpha, MemWObj, 
          DataModel, MR_Model,DataSNR);
}
