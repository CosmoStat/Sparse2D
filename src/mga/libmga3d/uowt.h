
#ifndef _UOWT_H
#define _UOWT_H

#include <fstream>
#include <vector>
#include "GlobalInc.h"
#include "MR3D_Obj.h"
#include "IM_IO.h"
#include "MGA_Inc.h"
#include "FloatTrans.h"

#include "Array.h"

class UOWT : public PAVE_3D_WT, public FloatTrans
{
protected:
// Base elements
	int DataNx;				// Input image size
	int DataNy;	
	int DataNz;
	int NbrScale;
	SubBandFilter *SB1D;
	
// Filtering properties
	bool no_fine;
	bool no_coarse;
	bool positivity;
	
public:

	UOWT(SubBandFilter *SB);
	~UOWT();
	
	fltarray *m_TabBand; // for the floating version : buffer pointing to the float*

// initialize and allocate the transform
	void init(fltarray * &TabBand, int nx, int ny, int nz, int _NbrScale);
	void alloc(fltarray& in, int _NbrScale);// for floating version
	
// Size and other getters
	inline int get_nbands() { return 7*NbrScale-6;}
	inline int get_sizebands() { return DataNx*DataNy*DataNz;}
	inline int size_transform() { return get_nbands()*get_sizebands();}

// Statistic tools
	void extract_stat(fltarray * TabBand);
	float get_maxabs(fltarray * TabBand);

// Transform / Recons alias
	void transform(fltarray &in, fltarray* TabBand) {PAVE_3D_WT::transform(in, TabBand, NbrScale);}
	void transform(fltarray &in, fltarray* TabBand, int ns){PAVE_3D_WT::transform(in, TabBand, ns);}
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(fltarray* TabBand, fltarray &out) {PAVE_3D_WT::recons(TabBand, out, NbrScale);}
	void recons(fltarray* TabBand, fltarray &out, int ns){PAVE_3D_WT::recons(TabBand, out, ns);}
	void recons(float* in, fltarray& out);

// Filtering methods
	void threshold(fltarray * TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void threshold_on_mask(fltarray * TabBand, fltarray * TabBand_mask, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
	void hard_threshold(float* in, float lvl, bool threshold_coarse=false);
	void wiener(fltarray * TabBand, float Noiselvl, int LocalBS);
	void fdr(fltarray * TabBand, float Alpha, float SigmaNoise);
	void stein_block_threshold(fltarray * TabBand, float SigmaNoise);
	void apply_positivity(fltarray & Data);

// Filtering properties
	inline void set_no_fine(bool nf) { no_fine=nf; }
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	inline void set_positivity(bool p) { positivity=p; }
	
// IO tools
	void write(char *Name, fltarray * TabBand);
	void read(char *Name, fltarray * &TabBand);

};
void uwt3d_clear(vector< vector< fltarray* > > &C);
void uwt3d_transform(fltarray &Data, vector< vector< fltarray* > > &vTabBand, int NbrScale3D, int wavelet_type);
void uwt3d_recons(vector< vector< fltarray* > > &vTabBand, fltarray &Data, int wavelet_type);
void uwt3d_threshold(vector< vector< fltarray* > > &vTabBand, float threshold, filter_type FilterType);
void uwt3d_filter(fltarray &Data, fltarray &Recons, int NbrScale3D, float threshold, filter_type FilterType, int wavelet_type);
 
/***********************************************************************/

#endif
