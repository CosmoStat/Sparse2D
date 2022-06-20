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
 
#ifndef _WTBORD_H_
#define _WTBORD_H_

#include "SB_Filter1D.h"
#include "BordLineCol.h"
#include "IM_Rot.h"

void tfd1D_bl_V0(float *Input,float *Output, int Np); /* transformee avec bords libres avec splines lineaires */
void tfd1Dinv_bl_V0(float *Input,float *Output, int Np); /* transformee avec bords libres avec splines lineaires */
void tfd1D_bl_V1(float *Input,float *Output, int Np); /* transformee avec bords libres avec splines quadratique */
void tfd1Dinv_bl_V1(float *Input,float *Output, int Np);/* transformee avec bords libres avec splines quadratique */

/***********************   Ondelettes a div nulle   **********************************************/

void tfodab2d(int N, float *c1,float *c2,float *ddn1,float *ddn2); /* transformee directe */
void tfodab2dinv(int N, float *c1,float *c2,float *ddn1,float *ddn2); /* transformee inverse */

/***********************   Ondelettes gradient   **********************************************/

void tfogab2d(int N, float *c1,float *c2,float *dg1,float *dg2); /* transformee directe */
void tfogab2dinv(int N, float *c1,float *c2,float *dg1,float *dg2); /* transformee inverse */
void qi10(float *u,float *v);
void qi01(float *u,float *v);
void rec10(float *u,float *v);
void rec01(float *u,float *v);
#endif
