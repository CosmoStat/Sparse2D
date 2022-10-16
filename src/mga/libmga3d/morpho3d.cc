/******************************************************************************
**                   Copyright (C) 1994 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck and Simona Mei
**
**    Date:  96/06/13
**    
**    File:  IM_Morpho.cc
**
*******************************************************************************
**
**    DESCRIPTION  morphological tools for image processing
**    -----------  
**
******************************************************************************* 
**
** morpho3d_dilation (fltarray &Imag1, fltarray& Imag2, int Window_Size)
**
** Imag2 = erosion(Imag1) with a square filter 3x3
** erosion in gray level
**
******************************************************************************* 
**
** morpho_dilation (Ifloat &Imag1, Ifloat &Imag2)
**
** Imag2 = dilation(Imag1) with a square filter 3x3
** dilatation in gray level
** 
*******************************************************************************/

 
// static char sccsid[] = "@(#)IM_Morpho.cc 3.2 96/06/13 CEA @(#)";


#include "morpho3d.h"


// ***************************************************************************

void morpho3d_dilation (fltarray &Imag1, fltarray& Imag2, int Window_Size)
{
	int i,j,k,l,m,n;
	int Nx = Imag1.nx();
	int Ny = Imag1.ny();
	int Nz = Imag1.nz();
	int Window2 = (Window_Size - 1) / 2;
	Imag2.init(0.F);

	for (i = 0; i < Nx; i++) 
	for (j = 0; j < Ny; j++) 
	for (k = 0; k < Nz; k++) 
	{
		if(Imag1(i,j,k) > 0.01)
			for (l = i - Window2; l <= i + Window2; l++)
			for (m = j - Window2; m <= j + Window2; m++)
			for (n = k - Window2; n <= k + Window2; n++)
				if (l>=0 & m>=0 & n>=0 & l<Nx & m<Ny & n<Nz)
					Imag2(l,m,n) = 1.F;
	}
}

// ***************************************************************************
 
void morpho3d_erosion (fltarray &Imag1, fltarray &Imag2, int Window_Size)
{
	int i,j,k,l,m,n;
	int Nx = Imag1.nx();
	int Ny = Imag1.ny();
	int Nz = Imag1.nz();
	int Window2 = (Window_Size - 1) / 2;
	Imag2.init(1.F);

	for (i = 0; i < Nx; i++) 
	for (j = 0; j < Ny; j++) 
	for (k = 0; k < Nz; k++) 
	{
		if(Imag1(i,j,k) < 0.01)
			for (l = i - Window2; l <= i + Window2; l++)
			for (m = j - Window2; m <= j + Window2; m++)
			for (n = k - Window2; n <= k + Window2; n++)
				if (l>=0 & m>=0 & n>=0 & l<Nx & m<Ny & n<Nz)
					Imag2(l,m,n) = 0.;
	}
}
