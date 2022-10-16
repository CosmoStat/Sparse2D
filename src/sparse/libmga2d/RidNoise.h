
/******************************************************************************
**                   Copyright (C) 200 by CEA + Stanford university
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  06/04/00
**    
**    File:  RidNoise.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for  correlated noise manadgement in the Ridgelet 
**    -----------  transform
**                 
******************************************************************************/

#ifndef _CRIDNOISE_H_
#define _CRIDNOISE_H_

#define RID_SIZE_HISTO 1024 /* number of bins for histogram */

#include "Ridgelet.h"

class StatRidNoise {
  fltarray StepHisto,MinHisto,MaxHisto;
  fltarray TabSigma;
  int NbrBand;
  fltarray Histo;
  fltarray Repart;
  fltarray Tab_Bin;

public:
  void init(Ifloat & RidTrans, int Nscale,  Ridgelet & Rid);

  Bool Verbose;
  int nbr_band() const {return  NbrBand;}
  float sigma_band(int Band)  const {return TabSigma(Band);}
  
  StatRidNoise() {NbrBand=0;Verbose = False;};
  void find_threshold(int NbrMaxBand, fltarray &Tab_Proba, fltarray &ThresholdMin, fltarray &ThresholdMax);
       // return the thresholds (min, max) corresponding 
       // to a Epsilon value for each scale

  void find_nsigma_threshold(int NbrMaxBand, fltarray &Tab_Proba, fltarray &ThresholdMin, fltarray &ThresholdMax);
       // return the thresholds (min, max) corresponding 
       // to a NSigma (equivalent Gaussian) detection value for each scale
 
  float prob(int Band,  float Val); 
       // return the probability to have a wavelet
       // coefficient of value Val
   
  float repartition(int Band, float Val);
       // return the  integrated probability to have a wavelet coefficient
       // lower than Val  (Repartition Function)
    
  void write(char *StepFileName); 
       // write histogram, repartition function and bin value 
  
  ~StatRidNoise(){;}
};

#endif
