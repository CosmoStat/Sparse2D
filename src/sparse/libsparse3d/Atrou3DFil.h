/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  18/09/98 
**    
**    File:  SB_Filter.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _Atrou3DFil_H_
#define _Atrou3DFil_H_

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Atrou3D.h"
#include "Mr3d_FewEvent.h"
#include "CUBE_Segment.h"

/************************************************************************/

class ATROUS_3D_FIL: public ATROUS_3D_WT {
         public:
 	  // Remove in the support isolated pixels
	  // i.e. if Support[b](i,j,k) > 0 and all neighboors are equal
	  // to zero, then  Support[b](i,j,k) is set to zero.
 	  void No_Single_Point(int Nx,int Ny, int Nz,intarray * & Support,
	  	               FewEventPoisson FEP,int Nbr_Scale, Bool Verbose=False);
		
          void threshold(intarray *& Support,  int NbPlan);
 	  void Filtering_Structure(intarray *& Support,int LastScale,Bool InfoBand,
	  	Bool RemoveBorder, int FistScaleDetect);
          void wttransform(fltarray & Cube, int Nbr_Plan,
	  	Bool WriteRes,Bool InfoBand, char Name_Cube_Out[]);
	  void wtrecons(fltarray & Cube, int Nbr_Plan, 
	  	Bool InfoBand,Bool AddLastScale, Bool Adjoint);
          fltarray *Wavelet_Coef;
};
 
/************************************************************************/
 
void free(fltarray * TabBand, int Nbr_Plan);
void free(intarray * TabBand, int Nbr_Plan);
void alloc_array (intarray * & Nb_Event, int Nx, int Ny, int Nz, int NbrBand);
void alloc_array (fltarray * & TabBand, int Nx, int Ny, int Nz, int NbrBand);
void init (fltarray * & TabBand, int NbrBand, float val=0.);

#endif
