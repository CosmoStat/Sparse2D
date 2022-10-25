/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.2
**
**    Author: 09/09/99
**
**    Date:  01/02/09
**    
**    File:  mc_wthfilter.h
**
*******************************************************************************/

#ifndef _MC_WTHFILTER_H
#define _MC_WTHFILTER_H

#include "mc2d_com.h"
#include "IM_Noise.h"
#include "MR_Psupport.h"
#include "mc_usage.h"
#include "wth_MC2D.h"

class wth_Param : public to_Param {
public:
   type_transform e_Transform;   /* type of transform */   
   sb_type_norm e_Norm;          /* type norm */
   type_sb_filter e_SBFilter;    /* SB Filter type */  
   FilterAnaSynt* po_FAS;        /* filter analisys and synthesis */
   int i_NbPlan;                 /* number of scales */  
   int i_NbBand;                 /* number of band */
   int i_NbIter;                 /* maximumnumber of iteration */
   Bool e_Verbose;               /* verbose exec */
   Bool e_WriteTransf;           /* write the Transf vector */
   Bool e_TraceParamIn;          /* Trace parametres In */
   Bool e_UseReconsAdjoint;      /* use rec_adjoint in Mr2d classes */  
   // noise model
   Bool e_WithNoiseModel;
   type_noise e_TypeNoise;
   float f_NoiseIma;
   float f_Gain;
   float f_SigmaGauss;
   float f_MeanGauss;  
   float f_NSigmaMr2d; 
   float f_NSigmaWth; 
   Bool e_WriteSupportMr2d;
   Bool e_WriteSupportWth; 
   type_lift e_LiftingTrans;
   Bool e_PositivImag;  
   Bool KillLastScale;  
public:
   wth_Param ();
   virtual void operator () (int argc, char *argv[]); /* read parameters */
private:
   virtual void read_image (char *argv[]); /* read images in */
   virtual void usage (char *argv[]);      /* trace usage */
   virtual void trace ();                  /* trace parameters */    
};


class wth_Result : public to_Result {
public:   
   MRNoiseModel* po_TabNoiseModel;
   WTH_MC2D      o_Mc2dWth;
   wth_Param*    po_Param;
public:
   wth_Result (wth_Param* ppo_Param);
   ~wth_Result ();   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute ();     /* curent loop work */
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   to_Param* res_GetpParam ();      /* get to_Param pointeur */
private:
   void res_InitNoiseModelWth ();
   void res_InitNoiseModelMr2d ();
};

class wth_Iter : public to_Iter {
public:
   wth_Result* po_Result;
   wth_Iter (wth_Result* ppo_Result) : to_Iter () {po_Result=ppo_Result;};
   virtual Bool iter_End ();
private:
   virtual void iter_Residu ();
   virtual to_Result* iter_getpResult();
};


#endif
