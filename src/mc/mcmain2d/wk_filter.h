/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.4
**
**    Author: 02/19/01
**
**    Date:  01/02/23
**    
**    File:  wk_filter.h
**
*******************************************************************************/

#ifndef _MC_PCAFILTER_H
#define _MC_PCAFILTER_H

#include "mc2d_com.h"
#include "IM_Noise.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "MR_PCA.h"
#include "MR_Psupport.h"
#include "mc_usage.h"
#include "IM3D_IO.h"

class pca_Param : public to_Param {
public:
   type_transform e_Transform;   /* type of transform */   
   sb_type_norm e_Norm;          /* type norm */
   type_sb_filter e_SBFilter;    /* SB Filter type */  
   FilterAnaSynt* po_FAS;        /* filter analisys and synthesis */
   Bool e_WithNoiseModel;        /* noise model is used */
   Bool e_ModNoiseModelPCA;      /* modify PCA noise model */
   Bool e_ModNoiseModelMr2d;     /* modify Mr2d noise model */
   type_noise e_TypeNoise;       /* type of noise */
   Bool e_PositivImag;           /* positiv constraint for image */
   Bool e_MaxImag;               /* max value for an image is 256 */
   CorrelComputeType e_CorrelComp;/* correlation compute type */   
   char tc_ImportCorMatName[356]; /* matrix correlation name */
   Ifloat* po_ImportCorMat;      /* correlation matrix */
   CorrelNoiseType e_CorrelNoise; /* correlation noise type */
   Bool e_NormCorrelMatrix;      /* normalize correl mat */
   float f_NoiseIma;             /* noise standard deviation */
   float f_Gain;                 /* CCD gain */
   float f_SigmaGauss;           /* CCD stand deviation of the read-out noise */
   float f_MeanGauss;            /* CCD mean of the read-out noise */
   float f_NSigmaMr2d;           /* Thresolding at NSigma * SigmaNoise */
   float f_NSigmaPCA;            /* Thresolding at NSigma * SigmaNoise */
   float f_ModNsigmaMr2d;        /* mod Mr2d Thresolding at Nsigma * SigmaNoise */
   float f_ModNsigmaPCA;         /* mod Pca Thresolding at Nsigma * SigmaNoise */
   int i_NbPlan;                 /* number of scales */   
   int i_NbBand;                 /* number of band */
   int i_NbIter;                 /* maximumnumber of iteration */
   Bool e_Verbose;               /* verbose exec */
   Bool e_WriteTransf;           /* write the Transf vector */
   Bool e_WriteEigenOnDisk;      /* write eigen vectors on disk */
   Bool e_WriteSupportMr2d;      /* write support Multiresolution */
   Bool e_WriteSupportPCA;       /* write support PCA */
   Bool e_SupIsolPixel;          /* supress iso pixel in support */
   Bool e_DilateSupportPCA;      /* dilate Support PCA */
   Bool e_DestroyRings;          /* destroy rings in pca noise support */
   Bool e_NormEigen;             /* Nomalise Eigen */
   Bool e_TraceParamIn;          /* Trace parametres In */
   Bool e_PCALinearThreshold;    /* Linear Threshold */
   Bool e_Mr2dLinearThreshold;   /* Linear Threshold */
   Bool e_UseReconsAdjoint;      /* use rec_adjoint in Mr2d classes */
   Bool e_WriteCorrelMatrix;     /* Write correl matrix */
   int ti_NbrUseEigen[MAX_NB_BAND];
   int tti_TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND];
public:
   pca_Param ();
   virtual void operator () (int argc, char *argv[]); /* read parameters */
private:
   virtual void read_image (char *argv[]); /* read images in */
   virtual void usage (char *argv[]);      /* trace usage */
   virtual void trace ();                  /* trace parameters */    
};


class pca_Result : public to_Result {
public:
   MRNoiseModel* po_TabNoiseModel;
   PCA_MR2D      o_Mr2dPca;
   pca_Param*    po_Param;
public:
   pca_Result (pca_Param* ppo_Param);
   ~pca_Result ();   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute ();     /* curent loop work */
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   to_Param* res_GetpParam ();      /* get to_Param pointeur */
private:
   void res_InitNoiseModelPca ();
   void res_InitNoiseModelMr2d ();
   void res_ModSupport (float f_ModNSigma);
   void res_DestroyRingsInPcaSup ();
};

class pca_Iter : public to_Iter {
public:
   pca_Result* po_Result;
   pca_Iter (pca_Result* ppo_Result) : to_Iter () {po_Result=ppo_Result;};
   virtual Bool iter_End ();
private:
   virtual void iter_Residu ();
   virtual to_Result* iter_getpResult();
};


#endif
