/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: A. Woiselle
**
**    Date:  June, 10th 2008
**    
**    File:  IM_Prob.cc
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise management
**    ----------- 
**                 
******************************************************************************/

#include "Arr_Prob.h"

void CArrProb::set(fltarray &NoiseArr, bool sym)
{
   CImaProb::set(NoiseArr.buffer(),  NoiseArr.n_elem(), sym);
}
