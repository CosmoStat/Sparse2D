/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Hubert Druesne, Philippe Querre
**
**    Date:  24/007/03 
**    
**    File:  MDCT1D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale DCT
**    -----------  
**                 
******************************************************************************/
 

#ifndef _1DMDCT_H_
#define _1DMDCT_H_

#include "GlobalInc.h"
#include "IM1D_IO.h"
#include "IM1D_Block.h"
#include "IM1D_Dct.h"
#include "IM1D_ALDCT.h" 
#include "SB_Filter1D.h" 

/************************************************************************/

class MDCT1D {
	   fltarray *TabTransWT;
	   fltarray **TabTransMDCT;
           ATROUS_1D_WT AWT;
	   int NbrScale;
	   inline void test_scale(int s)
	   {
	       if ((s < 0) || (s >= NbrScale))
	       {
	          cout << "Error: bad scale index. s = " << s << endl;
		  exit(-1);
	       }
	   }
	   inline void test_indi(int s, int i)
	   {
	       if ((i < 0) || (i >= TabDCT[s].nx()))
	       {
	          cout << "Error: bad i index. i = " << i  << endl;
		  exit(-1);
	       }
	   }
	   void reset()
	   {
	      TabTransWT = NULL; TabTransMDCT = NULL; 
	      TabDCT = NULL;
	      NbrScale = 0;
	      BlockOverlap=False;
	      Bord = I_MIRROR;Verbose=True;
	      Write=False;
	   }
           
           float update (float CoefSol, float Threshold, float SoftLevel);
        public:
        
          int FirstDetectScale;
          float SigmaNoise;
          float COSMin;
          float COS_Sensibility;
          bool Write;
          float NSigmaSoft;
          bool UseNormL1;
          bool OnlyPositivDetect;
        
	  LOCAL_DCT1D *TabDCT;
 	  int nbr_scale() {return NbrScale;}
	  int nx(int s) {test_scale(s);return TabDCT[s].nx();}
	  fltarray & band(int s) {test_scale(s);return TabDCT[s]._DCTSig;}
          // fltarray * trans_buffer()   { return TabTransMDCT;}
 	  inline float & operator() (int s, int i) 
	                        {test_scale(s);test_indi(s,i);
				return (TabDCT[s])(i);}
	  type_border Bord;
	  Bool BlockOverlap;
	  intarray TabBlockSize;
	  Bool Verbose;
	  MDCT1D () {reset();}
 	  void alloc(int Nx, int Nbr_Plan, int FirstBlockSize=8);
          void transform(fltarray &Signal);
          float norm_band(int s);
	  void recons(fltarray &Sig);
	  void threshold(float Noise, float NSigma);
          void KillScaleNotUsed (int FirstDetScale);
          void KillLastScale ();
          void threshold (float NSigma, int IterNumber);
          float getAbsMaxTransf ();
 	  void free() 
	         {if (TabTransMDCT != NULL) delete [] TabTransMDCT;
 		  if (TabDCT != NULL) delete [] TabDCT;
		  if (NbrScale != 0) AWT.free(TabTransWT, NbrScale);
	         }
 
	  ~MDCT1D(){ free();}
};

/************************************************************************/

#endif
