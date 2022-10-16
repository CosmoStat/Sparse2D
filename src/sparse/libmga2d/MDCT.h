/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  21/04/02 
**    
**    File:  MDCT.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale DCT
**    -----------  
**                 
******************************************************************************/
 

#ifndef _MDCT_H_
#define _MDCT_H_

#include "IM_Obj.h"
#include "IM_DCT.h"
#include "SB_Filter.h"

/************************************************************************/

class MDCT {
	   Ifloat *TabTransWT;
	   Ifloat **TabTransMDCT;
	   LOCAL_DCT2D *TabDCT;
           ATROUS_2D_WT AWT;
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
	       if ((i < 0) || (i >= TabDCT[s].nl()))
	       {
	          cout << "Error: bad i index. i = " << i  << endl;
		  exit(-1);
	       }
	   }
	   inline void test_indj(int s, int j)
	   {
	       if ((j < 0) || (j >= TabDCT[s].nc()))
	       {
	          cout << "Error: bad j index. j = " << j  << endl;
		  exit(-1);
	       }
	   }
	   void reset()
	   {
	      TabTransWT = NULL; TabTransMDCT = NULL; 
	      TabDCT = NULL;
	      NbrScale = 0;
	      BlockOverlap=False;
	      Bord = I_MIRROR;Verbose=False;
	   }
        public:
 	  int nbr_scale() {return NbrScale;}
	  int nl(int s) {test_scale(s);return TabDCT[s].nl();}
	  int nc(int s) {test_scale(s);return TabDCT[s].nc();}
	  Ifloat & band(int s) {test_scale(s);return TabDCT[s].DCTIma;}
          // Ifloat * trans_buffer()   { return TabTransMDCT;}
 	  inline float & operator() (int s, int i, int j) 
	                        {test_scale(s);test_indj(s,j);test_indi(s,i);
				return (TabDCT[s])(i,j);}
	  type_border Bord;
	  Bool BlockOverlap;
	  intarray TabBlockSize;
	  Bool Verbose;
          MDCT () {reset();}
 	  void alloc(int Nl, int Nc, int Nbr_Plan, int FirstBlockSize=8);
          void transform(Ifloat &Image);
          float norm_band(int s);
	  void recons(Ifloat &Imag);
	  void threshold(float Noise, float NSigma);
	  void free() 
	         {if (TabTransMDCT != NULL) delete [] TabTransMDCT;
 		  if (TabDCT != NULL) delete [] TabDCT;
		  if (NbrScale != 0) AWT.free(TabTransWT, NbrScale);
	         }
          ~MDCT(){ free();}
};

/************************************************************************/

#endif
