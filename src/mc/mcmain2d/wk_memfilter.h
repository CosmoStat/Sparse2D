/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.3
**
**    Author: 02/29/00
**
**    Date:  01/02/09
**    
**    File:  mc_mwfilter.h
**
*******************************************************************************/

#ifndef _MC_MWFILTER_H
#define _MC_MWFILTER_H

#include "mc2d_com.h"
#include "IM_Noise.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "MR_PCA.h"
#include "MR_Psupport.h"
#include "mc_usage.h"
//#include "XXX_MC2D.h"

class mw_Param : public to_Param {
public:
   type_transform e_Transform;   /* type of transform */   
   sb_type_norm e_Norm;          /* type norm */
   type_sb_filter e_SBFilter;    /* SB Filter type */  
   FilterAnaSynt* po_FAS;        /* filter analisys and synthesis */
   Bool e_WithNoiseModel;        /* noise model is used */
   type_noise e_TypeNoise;       /* type of noise */
   Bool e_PositivImag;           /* positiv constraint for image */
   Bool e_MaxImag;               /* max value for an image is 256 */   
   Bool e_NormCorrelMatrix;      /* normalize correl matrix */
   CorrelComputeType e_CorrelComp;/* correlation compute type */
   char tc_ImportCorMatName[356]; /* matrix correlation name */
   CorrelNoiseType e_CorrelNoise; /* correlation noise type */   
   Ifloat* po_ImportCorMat;      /* correlation matrix */
   float f_NoiseIma;             /* noise standard deviation */
   float f_Gain;                 /* CCD gain */
   float f_SigmaGauss;           /* CCD stand deviation of the read-out noise */
   float f_MeanGauss;            /* CCD mean of the read-out noise */
   float f_NSigmaMr2d;           /* Thresolding at NSigma * SigmaNoise */
   float f_NSigmaPCA;            /* Thresolding at NSigma * SigmaNoise */
   int i_NbPlan;                 /* number of scales */  
   int i_NbBand;                 /* number of band */
   int i_NbIter;                 /* maximumnumber of iteration */
   Bool e_Verbose;               /* verbose exec */
   Bool e_WriteTransf;           /* write the Transf vector */
   Bool e_WriteEigenOnDisk;      /* write eigen vectors on disk */
   Bool e_WriteSupportMr2d;      /* write support Multiresolution */
   Bool e_WriteSupportPCA;       /* write support PCA */
   Bool e_NormEigen;             /* Nomalise Eigen */
   Bool e_TraceParamIn;          /* Trace parametres In */   
   Bool e_UseReconsAdjoint;      /* use rec_adjoint in Mr2d classes */
   int ti_NbrUseEigen[MAX_NB_BAND];
   int tti_TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND];
   
   float f_CvgParam;             /* convergence parametre */
   float f_RegulVal;             /* regularisation parametre */
   int i_TypeOpt;                /* optimisation type */
   Bool e_DataEdge;              /* model for edge */
   Bool e_DataSNR;               /* use data SNR */
   Bool e_UseRMSMap;             /* use RMS Map */
   char tc_NameRMSMap[256];      /* name of RMS Map File */
   int i_SizeBlock;              /* size block for sigma clipping */
   int i_NiterClip;              /* number iter for sigma clipping */
   int i_NbIterEntrop;           /* number oe iteration for entrop programm */
   Bool e_NscaleOpt;             /*  */
   Bool e_TransfOpt;             /*  */
#ifdef LARGE_BUFF
   Bool e_OptZ;                  /* large buffer */
   int i_VMSSize;                /* large buffer */
   char tc_VMSName[1024];        /* large buffer */
#endif   
     
public:
   mw_Param ();
   virtual void operator () (int argc, char *argv[]); /* read parameters */
private:
   virtual void read_image (char *argv[]); /* read images in */
   virtual void usage (char *argv[]);      /* trace usage */
   virtual void trace ();                  /* trace parameters */    
};


class mw_Result : public to_Result {
public:   
   MRNoiseModel* po_TabNoiseModel;
   PCA_MR2D      o_Mr2dPca;
   mw_Param*     po_Param;
public:
   mw_Result (mw_Param* ppo_Param);
   ~mw_Result ();   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute ();     /* curent loop work */   
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   to_Param* res_GetpParam ();      /* get to_Param pointeur */
private:
   void res_InitNoiseModelPca ();
   void res_InitNoiseModelMr2d ();
   //void res_ModSupport (float f_ModNSigma);
};



class mw_Iter : public to_Iter {
public:
   mw_Result* po_Result;
   mw_Iter (mw_Result* ppo_Result) : to_Iter () {po_Result=ppo_Result;};
   virtual Bool iter_End ();
private:
   virtual void iter_Residu ();
   virtual to_Result* iter_getpResult();
};

#endif
