/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.1
**
**    Author: 08/23/99
**
**    Date:  01/03/02
**    
**    File:  wth_MC2D.h
**
*******************************************************************************/

#ifndef _WTH_MC2D_H
#define _WTH_MC2D_H

#include "IM_Obj.h"
#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_NoiseModel.h" 
#include "mc_filter.h"




class WTH_MC2D : public mc_filter {

private:
   int       i_Length;
   fltarray  o_wthMc2d;
   FilterAnaSynt o_FAS;
  
public:

   WTH_MC2D () {init();};
   WTH_MC2D (int Length) {
      init();
      filt_Alloc(Length);
   };
   void init() {
      o_FAS.Verbose = False;
      o_FAS.alloc(F_MALLAT_7_9);   
   }
   void filt_Alloc (int Length);
   ~WTH_MC2D ();
   mc_filter* filt_GetpFilter() {return this;}
   
   //void filt_Compute (MultiResol *TabMR, MRNoiseModel *TabMRNoise);
   
   void filt_Transform (MultiResol  *TabMRin, MultiResol *TabMRout,
                        MRNoiseModel *TabNoiseIn=NULL, 
			MRNoiseModel *TabNoiseOut=NULL);  
   
   void filt_InvTransform (MultiResol *TabMRin , MultiResol *TabMRout, Bool KillLastScale=True);
   
   void filt_TransfSignal (MultiResol *TabMRin , MultiResol *TabMRout,
                           MRNoiseModel *TabNoiseIn=NULL);
   
   void filt_Threshold  (MultiResol  *TabMRin_out, MRNoiseModel *TabNoiseOut);
   
   void filt_ComputeNoise  (MultiResol *TabMRin, MRNoiseModel *TabNoiseIn, 
                           MRNoiseModel *TabNoiseOut); 
   
   //void filt_Print ();
			   
};

#endif
