/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  IM_Morpho.h
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

#ifndef _MORPHO3D_H_
#define _MORPHO3D_H_

#include "GlobalInc.h"
#include "IM_IO.h"
#include <cmath>


void morpho3d_erosion (fltarray &Imag1, fltarray& Imag2, int Window_Size = 3);
void morpho3d_dilation (fltarray &Imag1, fltarray &Imag2, int Window_Size = 3);

#endif

