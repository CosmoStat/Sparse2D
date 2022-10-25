
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
**    Date:  1.3.99
**    
**    File:  MR_CorrNoise.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for  correlated noise manadgement
**    ----------- 
**                 
******************************************************************************/

#ifndef _CCRADNOISE_H_
#define _CCRADNOISE_H_

#define SIZE_CRAD_HISTO 1024 /* number of bins for histogram */

class RadNoiseMap {
  float StepHisto,MinHisto,MaxHisto;
  fltarray Histo;
  fltarray Repart;
  fltarray Tab_Bin;

public:
  Bool Verbose;  
  float Sigma;

  RadNoiseMap();

  void set(Ifloat &NoiseMap);

  void find_threshold(float Proba, float & ThresholdMin,  float &ThresholdMax);
       // return the thresholds (min, max) corresponding 
       // to a Epsilon value 
 
   void find_gthreshold(float Nsigma, float & ThresholdMin,  float &ThresholdMax);
       // return the thresholds (min, max) corresponding 
       // to a Nsigma detection
 
  float prob(float Val); 
       // return the probability to have a  coefficient of value Val
   
  float repartition(float Val);
       // return the  integrated probability to have a  coefficient
       // lower than Val  (Repartition Function)
    
  void write(); 
       // write histogram, repartition function and bin value 
  
  ~RadNoiseMap(){;}
};

#endif
