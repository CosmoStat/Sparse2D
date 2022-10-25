
/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  31.8.99
**    
**    File: MR_HaarPoissonFilter.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for poisson-data filtering using Haar transform
**    -----------  
**                 
******************************************************************************/

#ifndef _HFILTER_H_
#define _HFILTER_H_


/*****************************************************************************/
    
#define NBR_HAAR_POISSON_THRESHOLD 4
enum type_haar_poisson_threshold {HPF_LAMBDA_K_1,HPF_LAMBDA_K_2,HPF_LAMBDA_GB,HPF_PRESS};
#define DEF_HAAR_POISSON_THRESHOLD HPF_LAMBDA_K_1

#define NBR_HAAR_POISSON_BGR_MODEL 4
enum type_background_model {HPF_BGR_FLAT, HPF_BGR_MR, HPF_BGR_IMA, HPF_BGR_ITER};
#define DEF_HAAR_POISSON_BGR_MODEL  HPF_BGR_MR

#define TabHaar_MIN_EPS  1e-6 
#define TabHaar_MAX_EPS  1e-2 

#define DEF_NBR_ITER_REC_PER_SCALE 1
#define DEF_NBR_ITER_FILTER 4


inline const char *StringHaarThreshold  (type_haar_poisson_threshold  type)
{
    switch (type)
    {
        case  HPF_LAMBDA_K_1:
              return ("Kolaczyk-Dixon Threshold ");break;
        case  HPF_LAMBDA_K_2: 
              return ("Kolaczyk-Dixon Upper bound Threshold");break;
        case  HPF_LAMBDA_GB: 
              return ("Jammal-Bijaoui Threshold ");break;
	case  HPF_PRESS: 
              return ("PRESS-optimal filter");break;
        default:
              return ("Undefined sub-band filters");
              break;
    }
}
inline const char *StringHaarBgrModel(type_background_model type)
{
    switch (type)
    {
        case  HPF_BGR_FLAT:
              return ("Flat background");break;
        case  HPF_BGR_MR:
              return ("Multiresolution background estimation");break;
        case  HPF_BGR_IMA: 
              return ("Background image model");break;
	case HPF_BGR_ITER:  
	      return ("Iterative background estimation");break;
        default: 
              return ("Undefined sub-band filters");
              break;
    }
}


double get_harr_poisson_threshold(double Lambda, double Eps);
// return the threshold level using HPF_LAMBDA_GB method
// Lambda = number of count per pixel
// Eps = probability of false detection


class HFilter {
 void laplacien_recons(MultiResol &MR_Step, int s, Ifloat &Result);
 // void mr_ortho_regul_ima_rec(MultiResol &MR_Step, int s, Ifloat &Result, int Step=1);
 // iterative reconstruction to the previous resolution
 // MR_Step = mutliresolution transform with two scales
 // if 	TabLevel == NULL then LevelScale is the detection level
 //                           for all the scale
 // else TabLevel[w] is the detection level for the wth coefficient

  void threshold_scale(MultiResol &MR_Data,  Ifloat &LowResol, int s);
  // Threshold a band using one level per coefficient 
  // or the same level if TabLevel == NULL
  
  // variable use with HPF_LAMBDA_K_1 and HPF_LAMBDA_K_2
  double NormCoef1;
  double NormCoef2;  
  double Np;	    
 
  double get_level(float Lambda, int Scale);
  // return the threshold level for a given scale and a given lambda

  // void get_scale_level(MultiResol &MR_Step, int s, float *TabLevel);
  // get the threshold level using the low resolution scale as a background 
  // model: MR_Step.band(3) = background model
  // Norm = normalization coefficient when L2 normalization is used to
  //        get an unormalized coefficient
  // s = scale 
  // TabLevel : out = level table 
  
  // void get_model_scale_level(int Scale, Ifloat & ImaModel, MultiResol &MR_Step, float *TabLevel);
  // get the threshold level using an image model given by ImaModel
  // Scale: in = current scale 
  // ImaModel: in =  model image
  // MR_Step: in = 2 scales (4 bands) Haar transform at a given resolution 
  
  void get_mr_ima(Ifloat &Ima, int NbrScale);
  // set TabImaModel, the background model at each scale, from the
  // the background image Ima
  // Ima: in = background image model

  Ifloat *TabImaModel; // TabImaModel[s] contains the background model
                       // at scale s
		       // TabImaModel[s].nl() = Nl / 2^(s+1)
		       // TabImaModel[s].nc() = Nc / 2^(s+1)
		        
  type_background_model TypeBGR; // type background model
  type_haar_poisson_threshold TypeThreshold; // Type of thresholding
  
  int NbDetect;                  // Number of detected coeff.
                                 // Verbose usage only
  int Nlima;
  int Ncima;                     // image size to filter
  int NbrPlan;                   // number of scales
  int NbrBand;                   // number of bands
//   float *TabHLevel;              // level table  
//                                  // TabHLevel[w] = threshold level of the wth
// 				 // wavelet coefficient in a given scale
  Bool DecimatedTransform;       // True if it is a decimated haar transform

  void filter_noiter (MultiResol &MR_Data, Ifloat &Result);
         // filtering only for methods != HPF_BGR_ITER and not PRESS-optimal         

  void press_filter (Ifloat & Data, MultiResol &MR_Data, Ifloat &Result);
         // Filtering using the PRESS-optimal filter

  void alloc(type_background_model TBGR, 
              type_haar_poisson_threshold TThreshold, 
	      MultiResol &MR_Data, char *NameFileBgr);

 public:
  Bool Verbose;          
  int MaxIterPerScaleRec; // number of iterations for the regularisation
                          // scale per scale
  int FirstScale;      // threshold coef of scale < FirstScale
  Bool KeepPositivSup; // if true, threshold negative coeff.
  Bool SoftThreshold; // if apply a softthresholding instead of a hard pme
  int MaxFilterIter;  // number of iterations for HPF_BGR_ITER method
  double EpsilonPoisson; // Confidence interval for the detection
                         // for method HPF_LAMBDA_GB
  float N_Sigma;         // Confidence interval for the detection
                         // for method HPF_LAMBDA_K_1 and K_2
  Ifloat BGR_Ima;        // Background image model
  double Lambda;         // Number of count per pixel
                         // used only with HPF_BGR_FLAT	method	 
  int haar_norm_coef(int b, int Nb)
  // return the normalization coefficient when L2 normalization is used
  // An unormalized coefficient u_j = w_j * 2^{j+1}
  {
    int s = (b == Nb-1) ? (b-1) / 3: b/3;
    return POW2(s+1);
  }

  HFilter(type_background_model TBGR, type_haar_poisson_threshold TThreshold, 
	   MultiResol &MR_Data, char *NameFileBgr)
	  { alloc(TBGR, TThreshold, MR_Data, NameFileBgr); }
	    
  void filter (Ifloat & Data, MultiResol &MR_Data, Ifloat &Result);
        // filter an image
	// Data = input data
	// MR_Data = its Haar wavelet transform
	// Result: out = result
  void threshold_all_scale(MultiResol &MR_Data);

  void deconv(Ifloat & Data, Ifloat & Psf, MultiResol &MR_Data, Ifloat &Result,
              float RegulParam=0.1);
  void jammal_deconv(Ifloat & Data, Ifloat & Psf, MultiResol &MR_Data, Ifloat &Result,
              float RegulParam=0.1);
  void test(MultiResol &MR_Data);
   
  ~HFilter()
  { 
     NbrPlan=0;
     // if (TypeBGR != HPF_BGR_FLAT) delete [] TabHLevel;
     if ((TypeBGR == HPF_BGR_IMA) || (TypeBGR == HPF_BGR_ITER)) 
                                  delete [] TabImaModel;
  }
};


#endif
