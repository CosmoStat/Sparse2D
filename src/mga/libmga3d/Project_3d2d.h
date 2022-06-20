/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  27/05/2009
**    
**    File:  Project_3d2d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  Makes the projection on one direction : 
**    -----------  Partial Radon transform on only one direction, set by
**                 set_angle(int _a, int _b, int _c).
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _PROJECT3D2D_H_
#define _PROJECT3D2D_H_

#include "IM3D_PartialRadon.h"

class Project_3d2d : public PartialRadon3D {

private:
	
	int a_ref, b_ref, c_ref;
		
	// Override Comp_Plane to treat the single case (a,b,c)==(a_ref,b_ref,c_ref)
//	void comp_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage, Bool Inverse);
	void comp_odd_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage, Bool Inverse);
	void comp_FFT2D (fltarray& PartRad3DImage, cfarray& FPartRad3DImage);
public:
	Project_3d2d() {reset(); a_ref = 0; b_ref = 0; c_ref = 8;}
	~Project_3d2d() {reset();}
	void set_angle(int _a, int _b, int _c) {a_ref=_a; b_ref=_b; c_ref=_c;}
	void alloc();

};

#endif
