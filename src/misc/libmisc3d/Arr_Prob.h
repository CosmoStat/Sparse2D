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
**    File:  Arr_Prob.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for correlated noise manadgement
**    ----------- 
**                 
******************************************************************************/

#ifndef _CARRPROB_H_
#define _CARRPROB_H_

#include "GlobalInc.h"
#include "IM_Prob.h"

class CArrProb : public CImaProb
{
	public:
		CArrProb(){};
		~CArrProb(){};
		void set(fltarray &NoiseMap, bool sym=false);
};

#endif
