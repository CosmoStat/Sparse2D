
#ifndef _DOWT_H
#define _DOWT_H

#include <fstream>
#include <vector>
#include "GlobalInc.h"
#include "MR3D_Obj.h"
#include "IM_IO.h"
#include "MGA_Inc.h"

#include "Array.h"

class DOWT : public Ortho_3D_WT
{
// Base elements
	int DataNx;				// Input image size
	int DataNy;	
	int DataNz;
	int NbrScale;
	SubBandFilter *SB1D;
	
public:

	DOWT();
	DOWT(SubBandFilter *SB);
	~DOWT();
	
// initialize and allocate the transform
	void init(int nx, int ny, int nz, int _NbrScale);
	
// Size and other getters
	inline int nbr_band() { return 7*NbrScale-6;}
	int size_band_nx(int s, int b);
	int size_band_ny(int s, int b);
	int size_band_nz(int s, int b);
	int start_band_nx(int s, int b);
	int start_band_ny(int s, int b);
	int start_band_nz(int s, int b);

// Filtering methods
	void threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void wiener(fltarray &TabBand, float Noiselvl, int LocalBS);
	void fdr(fltarray &TabBand, float Alpha, float SigmaNoise);
	void stein_block_threshold(fltarray &TabBand, float SigmaNoise);

// IO tools
	void write(char *Name, fltarray &TabBand, bool Normalize);
	void read(char *Name, fltarray &TabBand, bool *NormalizeInv);

};
void coef_import(DOWT *dwt, vector< vector< fltarray* > > &vTabBand, fltarray &TabBand);
void coef_export(DOWT *dwt, fltarray &TabBand, vector< vector< fltarray* > > &vTabBand, bool alloc=false);

// Functions using ordered vectors for the bands
void dwt3d_clear(vector< vector< fltarray* > > &C);
void dwt3d_transform(fltarray &Data, vector< vector< fltarray* > > &vTabBand, int NbrScale3D, int wavelet_type);
void dwt3d_recons(vector< vector< fltarray* > > &vTabBand, fltarray &Data, int wavelet_type);
void dwt3d_threshold(vector< vector< fltarray* > > &vTabBand, float threshold, filter_type FilterType);

// Functions using a simple fltarray as coefficients
void dwt3d_transform(fltarray &Data, fltarray &TabBand, int NbrScale3D, int wavelet_type);
void dwt3d_recons(fltarray &TabBand, fltarray &Data, int NbrScale3D, int wavelet_type);
void dwt3d_threshold(fltarray &TabBand, int NbrScale3D, float threshold, filter_type FilterType);
void dwt3d_filter(fltarray &Data, fltarray &Recons, int NbrScale3D, float threshold, filter_type FilterType, int wavelet_type);

/***********************************************************************/

#endif
