
#ifndef _FCURSFrame_H
#define _FCURSFrame_H

#include "Array.h"
#include "IM_IO.h"
#include "IM_Obj.h"
#include "Usage.h"
#include "PCur.h"
#include "FCur.h"
#include "writefits3d.h"
#include "MGA_Inc.h"

/***********************************************************************/


class FCur_SFrame // Fast Curvelets Stack Frame
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNt;
	
	int NbrScale;
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)

	FCUR* TabCur;
	
	fltarray TabSigma;
	Ifloat TabStat;
	bool threshold_coarse;

//Methods
	
public:
	FCur_SFrame();
	~FCur_SFrame();
	
// Size and other getters
	void alloc_from_coarse(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=False, Bool IsotropWT=False, Bool RealCur=True);
	void alloc_from_fine(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=False, Bool IsotropWT=False, Bool RealCur=True);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Data, Ifloat **&TabBand, bool allocTB);
	
// Reconstruct a cube from its curvelet transform
	void recons(Ifloat **&TabBand, fltarray &Data);
	
// Statistic and information tools
	void get_norm_coeff(float NSigma);
	void noise_calib(Ifloat **&TabBand, char* Outname);
	void extract_stat(Ifloat **&TabBand, char* Outname);
	inline int nbr_bands(int s) {return TabCur->nbr_band(s);}
	
// Filtering methods
	void threshold(Ifloat **&TabBand, float SigmaNoise, float NSigma, filter_type FilterType);

// Filtering properties
	inline void set_threshold_coarse(bool th) { threshold_coarse=th; }
	
// IO tools
//	void write(char *Name, fltarray ** TabBand, bool Normalize);
//	void read(char *Name, fltarray ** &TabBand, bool *NormalizeInv);

//	void temp(fltarray** TabBand);
};

 
/***********************************************************************/

#endif
