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

#ifndef _MC1D_PCA_H
#define _MC1D_PCA_H

#include "mc1d_com.h"
#include "IM_Noise.h"
#include "CPca.h"
#include "IM_Pca.h"
#include "MR1D_Pca.h"
#include "mc_usage.h"




class pca_Param : public to_Param1d {
public:
   type_trans_1d e_Transform;     /* type of transform */
   sb_type_norm e_Norm;           /* type norm */
   type_sb_filter e_SBFilter;     /* SB Filter type */
   FilterAnaSynt* po_FAS;         /* filter analisys and synthesis */
   type_border e_Bord;            /* Border type */
   int i_NbrUndec;                /* nb undecimated scale */
   int i_NbPlan;                  /* number of scales */
   int i_NbBand;                  /* number of band */
   Bool e_Verbose;                /* verbose exec */
   Bool e_NormCorrelMatrix;       /* normalize correl mat */
   Bool e_NoDisplay;              /* no display on screen */
   Bool e_WriteTransf;            /* write the Transf vector */
   Bool e_WriteEigenOnDisk;       /* write Eigen on disk */
   Bool e_TraceParamIn;           /* Trace parametres In */
   Bool e_WriteCorrelMatrix;      /* Write correl matrix */   
   CorrelComputeType1d e_CorrelComp;/* correlation compute type */   
   char tc_ImportCorMatName[356]; /* matrix correlation name */
   Ifloat* po_ImportCorMat;       /* correlation matrix */
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
   PCA_MR1D      o_Mr1dPca;
   pca_Param*    po_Param;
public:
   pca_Result (pca_Param* ppo_Param);
   ~pca_Result ();   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute () {};  /* curent loop work */
   void res_Recons () {};          /* inv reconstruct process */
   void res_Write () {};           /* write 3D out result */
   to_Param1d* res_GetpParam ();     /* get to_Param pointeur */
};


#endif
