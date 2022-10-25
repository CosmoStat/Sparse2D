/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  Curvelet.h
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

#ifndef _CURVELETCON_H_
#define _CURVELETCON_H_

#include "IM_Obj.h"
#include "IM_Radon.h"
#include "SB_Filter.h"
#include "Ridgelet.h"
#include "MR_Obj.h"
#include "MR_Contrast.h"
#include "Curvelet.h"

/***********************************************************************/

class CurContrast: public Contrast
{
   inline float cur_contrast_function_velde(float x)
   {
       float ValRet=0.;
       if (x > FLOAT_EPSILON)
       {
          if (ABS(x) < Contrast_C_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/Contrast_C_Param), Contrast_P_Param);
          else if (ABS(x) < Contrast_M_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/ABS(x)), Contrast_P_Param);
          else ValRet = 1.;
       }
       return ValRet;
   }
   inline float cur_contrast_function(float x)
   {
       float ValRet=1.;
       float C2 = Contrast_C_Param*2.;
       float DiffC = Contrast_C_Param;
       
       if (x > FLOAT_EPSILON)
       {
          if (x < Contrast_C_Param) 
	  {
	    if (Filtering== True) ValRet = 0;
	    else ValRet = 1.;
	  }
          else if (x < C2) 
	    ValRet = (x - Contrast_C_Param)/DiffC  * 
	          pow( (double) (Contrast_M_Param/Contrast_C_Param), Contrast_P_Param)
		  +  (C2-x)/DiffC;
 	  else if (x < Contrast_M_Param) 
	    ValRet = pow( (double) (Contrast_M_Param/x), Contrast_P_Param);
          else if (Contrast_S_Param==0.) ValRet = 1.;
	  else ValRet = pow( (double) (Contrast_M_Param/x), Contrast_S_Param);
       }
       return ValRet;
  }
   public:
  CurContrast(): Contrast() {Contrast_M_ParamCoef=0.5;Contrast_S_Param=0.;
              NSigmaLow=5.;NSigmaUp=20.;UseSigmaUp=False;
	      L_Coeff=1.;Filtering=False;UseVeldeFunction=False;}
  // Parameter for the contrast function	
  Bool UseVeldeFunction;   // Use the Velde enhancement function	       
  double NSigmaLow;        // curvelet coefficient between NSigmaLow and 
  double NSigmaUp;         // NSigmaUp are amplified
  double Contrast_S_Param; // Saturation parameter
  Bool Filtering;
  Bool UseSigmaUp;
  float L_Coeff;          //  Coeff larger than L*Max are set to L*Max.     
  double Contrast_M_ParamCoef;// Transform coefficients larger than  
                              //     MaxBandCoef/Contrast_M_ParamCoef    
			      //  are not modified. M must be in ]0,1]
  void cur_enhance(Ifloat &Data, Curvelet & Cur, Bool PosIma=True);
  void curcol_enhance_luminance(fltarray &Data, Curvelet & Cur, Bool UseYUV=False);
  void curcol_enhance(fltarray &Data, Curvelet & Cur, Bool UseYUV=False);
};

/***********************************************************************/
#endif
