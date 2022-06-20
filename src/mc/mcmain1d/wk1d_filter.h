/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.3
**
**    Author: 01/23/01
**
**    Date:  01/02/09
**    
**    File:  mc_pca.h
**
*******************************************************************************/

#ifndef _MC1D_PCAFILTER_H
#define _MC1D_PCAFILTER_H

#include "IM_Noise.h"
#include "mc1d_com.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "MR1D_Pca.h"
#include "mc_usage.h"
//#include "Filter.h"


class pca_Param : public to_Param1d {
public:
   type_trans_1d e_Transform;     /* type of transform */
   sb_type_norm e_Norm;           /* type norm */
   type_sb_filter e_SBFilter;     /* SB Filter type */
   FilterAnaSynt* po_FAS;         /* filter analisys and synthesis */
   Bool e_WithNoiseModel;         /* noise model is used */
   type_noise e_TypeNoise;        /* type of noise */
   float f_NoiseSpectre;          /* noise standard deviation */
   float f_NSigmaMr1d;            /* Thresolding at NSigma * SigmaNoise */
   float f_NSigmaPCA;             /* Thresolding at NSigma * SigmaNoise */
   Bool e_PositivSpectre;         /* positiv constraint for image */
   Bool e_MaxSpectre;             /* max value for an image is 256 */
   type_border e_Bord;            /* Border type */
   int i_NbrUndec;                /* nb undecimated scale */
   int i_NbPlan;                  /* number of scales */
   int i_NbBand;                  /* number of band */
   int i_FirstScale;              /* first detetction scale */
   Bool e_Verbose;                /* verbose exec */
   Bool e_NormCorrelMatrix;       /* normalize correl mat */
   Bool e_NoDisplay;              /* no display on screen */
   Bool e_WriteTransf;            /* write the Transf vector */
   Bool e_WriteEigenOnDisk;       /* write Eigen on disk */
   Bool e_TraceParamIn;           /* Trace parametres In */
   Bool e_WriteCorrelMatrix;      /* Write correl matrix */   
   Bool e_WriteSupportMr1d;       /* write support Multiresolution */
   Bool e_WriteSupportPCA;        /* write support PCA */
   Bool e_SupIsolPixel;           /* supress iso pixel in support */
   Bool e_OnlyPositivDetect;      /* detect only on positiv pixel */
   Bool e_RemoveLastScale;        /* remove last scale */
   Bool e_OptimSoftThreshold;     /* call optim proc with soft threshold */
   Bool e_UseReconsAdjoint;       /* use rec_adjoint in Mr2d classes */
   Bool e_TresholdInMrSpace;      /* threshold in mr spcae with Sigma tab */
   CorrelComputeType1d e_CorrelComp;/* correlation compute type */   
   char tc_ImportCorMatName[356]; /* matrix correlation name */
   Ifloat* po_ImportCorMat;       /* correlation matrix */   
   CorrelNoiseType1d e_CorrelNoise; /* correlation noise type */
   int ti_NbrUseEigen[MAX_NB_BAND];
   int tti_TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND];
   fitsstruct ts_Header;          /* Header fits */
   char tc_Cmd[356];              /* command */
public:
   pca_Param ();
   virtual void operator () (int argc, char *argv[]); /* read parameters */
private:
   void read_spectre (char *argv[]); /* read images in */
   void usage (char *argv[]);      /* trace usage */
   void trace ();                  /* trace parameters */    
};


class pca_Result : public to_Result1d {
public:
   MR1DNoiseModel* po_TabNoiseModel;
   PCA_MR1D        o_Mr1dPca;
   pca_Param*      po_Param;
   fltarray        MrSigmaLevel;
   int             LoopNumber;
                   // use in soft threshold case
   fltarray        FirstPcaCoef;   
   fltarray        ResidualPcaCoef;
public:
   pca_Result (pca_Param* ppo_Param);
   ~pca_Result ();   
   void res_FirstCompute ();       /* first loop work */ 
   void res_CurrentCompute ();     /* curent loop work */
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   to_Param1d* res_GetpParam ();   /* get to_Param pointeur */
private:
   void res_InitNoiseModelPca ();
   void res_InitNoiseModelMr1d ();
   void res_SoftThreshold ();
   void res_SoftInitAlloc ();
   void res_ComputeCommonSupport ();
};

class pca_Iter : public to_Iter1d {
public:
   pca_Result* po_Result;
   pca_Iter (pca_Result* ppo_Result) : to_Iter1d () {po_Result=ppo_Result;};
   virtual Bool iter_End ();
private:
   virtual void iter_Residu ();
   virtual to_Result1d* iter_getpResult();
};

class pca_SoftIter : public to_SoftIter1d {
public:
   pca_Result* po_Result;
   pca_SoftIter (pca_Result* ppo_Result) : to_SoftIter1d () {
      po_Result=ppo_Result;};
   virtual Bool iter_End ();
private:
   virtual void iter_ActionOnSol();
   virtual to_Result1d* iter_getpResult();
};


#endif
