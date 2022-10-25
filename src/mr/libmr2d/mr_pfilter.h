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
**    Date:  00/07/13
**    
**    File:  mc_pcafilter.h
**
*******************************************************************************/

#ifndef _MR_PFILTER_H
#define _MR_PFILTER_H

#include "mr_com.h"
#include "IM_Noise.h"
#include "MR_Psupport.h"
#include "mr_FiltUsage.h"
#include "MR_Filter.h"


class pfilt_Param {
public:

   Iint          o_EventImage;
   fitsstruct    s_Header; 
   Ifloat        o_Abaque;  
   Iint          o_EventCorImag;
   Ifloat        o_CorImag; 
   
   char     tc_NameIn[256];      /* input name */
   char     tc_NameOut[256];     /* output name */
   Ifloat   o_Imag;              /* data in */   
   int      i_NbLin,             /* data in dimension */
            i_NbCol, 
	    i_NbIter;  
	    
   int i_NbScale;                /* number of scales */ 
   int i_MinNumberEvent;         /* number min of event */
   int i_FirstScale;             /* first scale detection */
   float f_ConvParam;            /* convergence parametre */
   char tc_NameAbaque[80];       /* file name abaque */  
   char tc_CorectedImag[80];     /* corrected (transmission) image */ 
   type_filter e_FilterType;     /* filter type */           
   Bool e_UseAdjoint;            /* recons with adjoint */
   Bool e_DetectOnlyPos;         /* detect only positiv structure */
   Bool e_PositivImag;            /* image must be positive  */
   Bool e_DilSupport;            /* dilate support */
   Bool e_DelIsolPix;            /* delete isolates pixels */
   Bool e_KillLastScale;         /* kill last scale */
   Bool e_WriteAll;              /* write MR transf and threshold coef */
   Bool e_OptimSoftThreshold;    /* uses soft iter filtering method */
   Bool e_Ascii;                 /* ascii file name */
   Bool e_Imag;                  /* image file name */
   Bool e_CorImag;               /* corected image given */
   Bool e_TwoStep;               /* corected image given */
   Bool e_OptZ;                  /* used memory */
   type_border e_Border;         /* border type */
   int i_VMSSize;                /* memory size */
   char tc_VMSName[1024];        /* memory directory */
   Bool e_Trace;                 /* trace param in */
   char tc_Cmd[256];             /* command prog */
   Bool e_Verbose;               /* verbose exec */

public:
   pfilt_Param ();
   void operator () (int argc, char *argv[]); /* read parameters */
private:
   void usage (char *argv[]);      /* trace usage */
   void trace ();                  /* trace parameters */ 
   void read_image_and_abaque ();  /* read image and abaque */
};


class pfilt_Result {
public:
   Ifloat        o_DataOut;  
   MultiResol*   po_Mr2d;        
   MRNoiseModel  o_NoiseModel;
   pfilt_Param*  po_Param;
   float         SigmaNoise;
   float         CurSigma;
   float         OldSigma;
   int             LoopNumber;
   fltarray        FirstCoef;   
   fltarray        ResidualCoef;
public:
   pfilt_Result (pfilt_Param* ppo_Param);
   ~pfilt_Result () {};   
   void res_FirstCompute ();       /* first loop work */  
   void res_CurrentCompute ();     /* curent loop work */
   void res_Recons ();             /* inv reconstruct process */
   void res_Write ();              /* write 3D out result */
   pfilt_Param* res_GetpParam ();  /* get to_Param pointeur */
   
   void res_SoftInitAlloc ();
   void res_SoftThreshold ();
   
};

class pfilt_Iter  { 
public:  
   Ifloat     o_FilteredImag;
   Ifloat     o_Init;
   int        NbIter;
   pfilt_Result* po_Result;
public:
   pfilt_Iter (pfilt_Result* ppo_Result) {NbIter=1;po_Result=ppo_Result;};
   void iter_Begin ();
   Bool iter_End ();   
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~pfilt_Iter() {};
private:
   virtual void iter_Residu ();
   virtual pfilt_Result* iter_getpResult();
};


class pfilt_SoftIter  {  
public:   
   int        NbIter;
   pfilt_Result* po_Result;
   pfilt_SoftIter (pfilt_Result* ppo_Result) {
      NbIter=1;po_Result=ppo_Result;};
   void iter_Begin ();  
   virtual Bool iter_End ();   
   virtual void operator++ ();
   virtual void operator++ (int);
   virtual ~pfilt_SoftIter() {};
private:
   virtual void iter_ActionOnSol();
   virtual pfilt_Result* iter_getpResult();
};

#endif
