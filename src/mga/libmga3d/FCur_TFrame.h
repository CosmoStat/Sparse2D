
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
#include "MeyerWT1D.h"

/***********************************************************************/


class FCur_TFrame
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNt;
	
	int NbrScale;
	int NbrScale1D;
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	Bool RealBand;					// If true, the redundancy inside the curvelet transform is used to remove
									// the imaginary part.Then the imaginary part equals to zero.
	Bool ExtendWT;
	
	FCUR* TabCur;
	
	int kill_time_fine;
	
	fltarray TabSigma;
//Methods
	
public:
	FCur_TFrame();
	~FCur_TFrame();
	
// Size and other getters
	inline int real() { return RealBand;}

	void alloc_from_coarse(int Nbr_Scale, int _NbrScale1D, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=False, Bool IsotropWT=False, Bool RealCur=True);
	void alloc_from_fine(int Nbr_Scale, int _NbrScale1D, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=False, Bool IsotropWT=False, Bool RealCur=True);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Data, fltarray ***&TabBand, bool allocTB);
	
// Reconstruct a cube from its curvelet transform
	void recons(fltarray ***&TabBand, fltarray &Data);
	
// Statistic and information tools
	void get_norm_coeff(float NSigma);
	void noise_calib(fltarray ***&TabBand, char* Outname);
	void extract_stat(fltarray ***&TabBand, char* Outname);
	inline int nbr_bands(int s) {return TabCur[0].nbr_band(s);}
	
// Filtering methods
	void threshold(fltarray ***&TabBand, float SigmaNoise, float NSigma, filter_type FilterType);
	void set_kill_time_fine(int nf) { kill_time_fine=nf; }
// IO tools
//	void write(char *Name, fltarray ** TabBand, bool Normalize);
//	void read(char *Name, fltarray ** &TabBand, bool *NormalizeInv);

//	void temp(fltarray** TabBand);
};

 
/***********************************************************************/

#endif
