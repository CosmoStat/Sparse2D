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
**    Date:  18/09/98 
**    
**    File:  SB_Filter.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  Line Column Multiscale  decomposition
**    ----------- 
**
******************************************************************************/

#ifndef _BORDLINECOL_H_
#define _BORDLINECOL_H_

#include "SB_Filter1D.h"
#include "BordLineCol.h"
#include "IM_Rot.h"

/****************************************************************************/
// 2D decimated  BordLineCol transformation
/*************************************************************/

class  BordLineCol {

       SubBand1D *Ptr_SB1D_LINE; // 1D subband decomposition along the lines
	   SubBand1D *Ptr_SB1D_COL;  // 1D subband decomposition along the columns
       SubBand1D *Ptr_SB1D;  // equal to Ptr_SB1D_LINE
	   
       Rotation Rot;
       inline void test_class()
       {
          if (Ptr_SB1D == NULL)
	  {
	     cout << "Error in BordLineCol: Class not allocated ... " << endl;
	     exit(-1);
	  }
       }
	   
        public:
        Bool Mirror;
	BordLineCol () {Ptr_SB1D = Ptr_SB1D_LINE = NULL; Ptr_SB1D_COL = NULL; Mirror=False;}
	BordLineCol (SubBand1D &SB1D) { alloc(SB1D);}
	void alloc (SubBand1D &SB1D, Bool UseMirror=False);
	void alloc(SubBand1D &SB1D_Line, SubBand1D &SB1D_Col, Bool UseMirror=False);

         // Constructor: SB1D = subband method to use

        SubBand1D *get_subband_method() {return Ptr_SB1D;}
	    SubBand1D *get_line_subband_method() {return Ptr_SB1D_LINE;}
        SubBand1D *get_col_subband_method() {return Ptr_SB1D_COL;}

        // return the subband decomposition method
        void transform_bord_V0_V1(Ifloat &Data);
        void transform_bord_V1_V0(Ifloat &Data);
        void recons_bord_V0_V1(Ifloat & Data);
        void recons_bord_V1_V0(Ifloat & Data);

        void transform(Ifloat &Data, int NbrScaleLine, int NbrScaleCol=-1);
        void recons(Ifloat &Data, int NbrScaleLine, int NbrScaleCol=-1);
        void transform(Ifloat &Data, Ifloat &Result, 
                               int NbrScaleLine, int NbrScaleCol=-1);
        void recons(Ifloat &Trans, Ifloat &Data, 
                               int NbrScaleLine, int NbrScaleCol=-1);

        void transform_one_step_line (Ifloat & Data, int Nl, int Nc, int Dep=0);
        void recons_one_step_line (Ifloat & Data, int Nl, int Nc, int Dep=0);
        void transform_one_step_col (Ifloat & Data, int Nl, int Nc, int Dep=0);
        void recons_one_step_col(Ifloat & Data, int Nl, int Nc, int Dep=0);
 
        void directional_transform(Ifloat &Data, Ifloat &DataTrans, 
	             float AngleRot,  int NbrScaleLine, int NbrScaleCol, 
	             Bool RadianUnitAngle=True);
        void directional_recons(Ifloat &DataTrans, Ifloat &Data, float AngleRot, 
                                         int NbrScaleLine, int NbrScaleCol, 
					 Bool RadianUnitAngle=True);
					 	
        void undec_transform_one_step_line (Ifloat & Data, Ifloat & Low, Ifloat & High,  int Step);
        void undec_recons_one_step_line (Ifloat & Low, Ifloat & High, Ifloat & Data,  int Step);
        void undec_transform_one_step_col (Ifloat & Data, Ifloat & Low, Ifloat & High,  int Step);
        void undec_recons_one_step_col(Ifloat & Low, Ifloat & High, Ifloat & Data,  int Step);
        void undec_transform (Ifloat &Data, Ifloat ** & TabTrans, int NbrScaleLine, int NbrScaleCol=-1);
	void undec_transform(Ifloat &Data, Ifloat * & TabTrans, int NbrScaleLine, int NbrScaleCol);

        void undec_recons(Ifloat ** & TabTrans, Ifloat &Data, int NbrScaleLine, int NbrScaleCol=-1);
        void undec_recons(Ifloat * & TabTrans, Ifloat &Data, int NbrScaleLine, int NbrScaleCol=-1);
	~BordLineCol (){Ptr_SB1D = NULL;}
};

/************************************************************************/

class  DirectionalBordLineCol {
       SubBand1D *Ptr_SB1D;
       Rotation Rot;
       inline void test_class()
       {
          if (Ptr_SB1D == NULL)
	  {
	     cout << "Error in BordLineCol: Class not allocated ... " << endl;
	     exit(-1);
	  }
       }
       void reset()
       {
          Ptr_SB1D = NULL; 
	  Mirror=False;
          NbrUndecimatedScaleLine=0;
	  NbrUndecimatedScaleCol=0;
       }
        Bool is_decimated_line_scale(int s) 
	    {Bool Ret = (s >= NbrUndecimatedScaleLine) ? True: False;
	     return Ret;}
	Bool is_decimated_col_scale(int s) 
	    {Bool Ret = (s >= NbrUndecimatedScaleCol) ? True: False;
	     return Ret;}    
       void transform_col (Ifloat & Data, int Nl, int Nc, int NScale, int Dep);
       public:
        Bool Mirror;
	int NbrUndecimatedScaleLine;
	int NbrUndecimatedScaleCol;
	DirectionalBordLineCol () {reset();}
        DirectionalBordLineCol (SubBand1D &SB1D) { alloc(SB1D);}
	void alloc (SubBand1D &SB1D) {reset();Ptr_SB1D = &SB1D;}
        // Constructor: SB1D = subband method to use

        SubBand1D *get_subband_method() {return Ptr_SB1D;}
        // return the subband decomposition method
  
        void transform(Ifloat & Data, Ifloat & Trans, int NbrScaleLine, int NbrScaleCol);
        void recons(Ifloat & Trans, Ifloat & Data, int NbrScaleLine, int NbrScaleCol);

	~DirectionalBordLineCol (){Ptr_SB1D = NULL;}
};

/************************************************************************/

#endif
