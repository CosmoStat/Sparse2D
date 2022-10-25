/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre, Hubert Druesne
**
**    Date:  6/08/03
**    
**    File:  IM1D_ALDCT.h
**
**    Modification history :
**
******************************************************************************/

#ifndef _1DALDCT_H_
#define _1DALDCT_H_

#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM1D_Block.h"
#include "IM1D_Dct.h"

//#include "IM_Obj.h"
//#include "IM1D_Dct.h"
//#include "SB_Filter1D.h"

/************************************************************************/

class BASIS{
   public:
      BASIS () {}
      intarray levels;
      fltarray* contents;
      int indice;
      void alloc(int NbrBand);
      ~BASIS () {}
};

void sig_inv (fltarray& Trans, fltarray& Sig, int BlockSize, int Nx, 
                      Bool Overlap, Bool WeightFirst);


class  ALDCT {
   
      LOCAL_DCT1D AL_ldct;
      float update (float CoefSol, float Threshold, float SoftLevel);
      float compute_base (fltarray* parents, int ind_left, int ind_right,
                          int level, int Dir, int infocost);   
      inline float norm() { return 1.;}
public:
   
      float SigmaNoise;
      float Sensibility;
      int NbScale;
      int NbPts;
      bool Verbose;
      bool Write;    
      bool UseNormL1;
      float NSigmaSoft;
        
      BASIS base;
      ALDCT () {}
      void alloc (fltarray* & TabALDCT, int Nx, int NbrBand);
      void transform (fltarray& Signal, fltarray* TabALDCT); 
      void threshold (fltarray* TabALDCT, float NSigma, int IterNumber);  
      void recons(fltarray& SigRec);
      float cost_L1(fltarray& Sig);
      float cost_L2(fltarray& Sig);
      float ml2logl2(fltarray& Sig);
      float getAbsMaxTransf (fltarray* TabALDCT);
      void getBestBasis (fltarray* TabALDCT, fltarray& BestBasis, int Infocost);
      ~ALDCT(){}

};

#endif
