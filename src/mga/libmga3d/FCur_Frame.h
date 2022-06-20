
#ifndef _FCURTFRAME_H
#define _FCURTFRAME_H

#include "Array.h"
#include "IM_IO.h"
#include "IM_Obj.h"
#include "Usage.h"
#include "PCur.h"
#include "FCur.h"
#include "writefits3d.h"
#include "MGA_Inc.h"

/***********************************************************************/


class FCur_Frame
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNt;
	
	int NbrScale;
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	Bool RealBand;					// If true, the redundancy inside the curvelet transform is used to remove
									// the imaginary part.Then the imaginary part equals to zero.
	Bool ExtendWT;
	
	FCUR* TabCur;
	
//Methods
	
public:
	FCur_Frame();
	~FCur_Frame();
	
// Size and other getters
	inline int real() { return RealBand;}

	void alloc_from_coarse(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=False, Bool IsotropWT=False, Bool RealCur=True);
	void alloc_from_fine(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT, Bool IsotropWT, Bool Real);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Data, bool allocTB);
	
// Reconstruct a cube from its curvelet transform
	void recons(fltarray &Data);
	
// Statistic and information tools
	void get_norm_coeff(float NSigma);

// Filtering methods
	void threshold(float SigmaNoise, float NSigma, filter_type FilterType);

// IO tools
//	void write(char *Name, fltarray ** TabBand, bool Normalize);
//	void read(char *Name, fltarray ** &TabBand, bool *NormalizeInv);

//	void temp(fltarray** TabBand);
};

 
/***********************************************************************/

#endif
