/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: 10/30/02
**
**    Date:  02/10/30
**    
**    File:  wk_trec.h
**
*******************************************************************************/

#ifndef _MC_PCA_REC_H
#define _MC_PCA_REC_H

#include "mc2d_com.h"
#include "IM_Noise.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "MR_PCA.h"
#include "MR_Psupport.h"
#include "mc_usage.h"

class pca_Param : public to_Param {
public:
   type_transform e_Transform;   /* type of transform */
   sb_type_norm e_Norm;           /* type norm */
   type_sb_filter e_SBFilter;     /* SB Filter type */
   FilterAnaSynt* po_FAS;         /* filter analisys and synthesis */
   type_border e_Bord;            /* Border type */
   int i_NbrUndec;                /* nb undecimated scale */
   int i_NbPlan;                 /* number of scales */   
   int i_NbBand;                 /* number of band */
   Bool e_Verbose;               /* verbose exec */
   Bool e_NormCorrelMatrix;       /* normalize correl mat */
   Bool e_NoDisplay;              /* no display on screen */
   Bool e_WriteTransf;           /* write the Transf vector */
   Bool e_WriteEigenOnDisk;      /* write eigen vectors on disk */
   Bool e_TraceParamIn;          /* Trace parametres In */
   Bool e_WriteCorrelMatrix;      /* Write correl matrix */   
   CorrelComputeType e_CorrelComp;/* correlation compute type */   
   char tc_ImportCorMatName[356]; /* matrix correlation name */
   Ifloat* po_ImportCorMat;       /* correlation matrix */
   int ti_NbrUseEigen[MAX_NB_BAND];
   int tti_TabKill[MAX_NBR_PCA_IMA][MAX_NB_BAND];
public:
   pca_Param ();
   void operator () (int argc, char *argv[]); /* read parameters */
private:
   void read_image (char *argv[]); /* read images in */
   void usage (char *argv[]);      /* trace usage */
   void trace ();                  /* trace parameters */    
};


class pca_Result : public to_Result {
public:
   PCA_MR2D      o_Mr2dPca;
   pca_Param*    po_Param;
public:
   pca_Result (pca_Param* ppo_Param);
   ~pca_Result ();   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute () {};  /* curent loop work */
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   to_Param* res_GetpParam ();     /* get to_Param pointeur */
};

#endif
