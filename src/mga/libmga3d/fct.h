
#ifndef _FASTCUR3D_H
#define _FASTCUR3D_H

#ifdef WINDOWS
	#define USE_OMP_FC3D 0
#else
	#define USE_OMP_FC3D 1
#endif

#if USE_OMP_FC3D
	#include <omp.h>
#endif

#include <fstream>
#include "GlobalInc.h"
#include "IM_IO.h"
#include "FFTN.h"
#include "FFTN_3D.h"
#include "MeyerWT3D.h"
#include "IM_Prob.h"
#include "IM_Math.h"
#include "Arr_Prob.h"
#include "MGA_Inc.h"

class FCurvelet3D : public MEYER_WT3D
{
protected:
// Base elements
	int DataNx;						// Input image size
	int DataNy;	
	int DataNz;
	Bool ModifSize;					// temp variable, to make the data odd-sized
	int NewNx;						// Input cube sizes must be odd. If it is not, we etend it by one in each requiered direction
	int NewNy;
	int NewNz;
	int NbrScale;					// Number of scales
	int NbrBand;					// total Number of Band 
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	type_extraction extract_type;	// Type of wedge extraction : backward, forward or both
	intarray TabNbrAnglePerScale;	// Number of directions per scale
	intarray TabNbrBandPerScale;	// Number of bands per scale
	double fw(double x);			// wedge angular window
	Bool RealBand;					// If true, the redundancy inside the curvelet transform is used to remove
									// the imaginary part.Then the imaginary part equals to zero.

	bool tabsigma;					// Specifies if the exact normalizing coefficients are known
	fltarray TabSigmaNoise;			// Contains the noise level at each scale. used only if TSN(0)>=0
	bool eval_stat;					// to extract statistics
	bool th_coarse;					// to threshold the coarse scale
	
	intarray *TabSizeNx;  // Band size 
	intarray *TabSizeNy;
	intarray *TabSizeNz;
	
// Filtering properties
	bool no_norm;					// If set, then no normalization is done
	bool no_fine;					// kills the finest scale
	bool no_coarse;					// kills the coarsest scale
	bool use_min, use_max;			// saturates the output between min_value and max_value
	float min_value, max_value;
	bool lowmem;					// low memory algorithm (filtering)
	fltarray ** TabBand_norm;		// external band normalization (for mask correlations)
	
//Methods
	void ifft_tabcube(cfarray ** & TabCF_Cube, fltarray ** & Tab_Cube, Bool Alloc);
	void get_size();
	float get_noise_lvl(float sigma, int s);
	virtual void get_wedges(cfarray * &TabWT, fltarray ** TabBand); // writes in TabBand only if RealCur
	void put_wedges(cfarray * &TabWT, fltarray ** TabBand); // reads TabBand only if RealCur
//	void scale_into_wedge(cfarray& Scale, cfarray*& TabCFBand, int nd, double L1, double L2, double L3, bool backward=false);
	void get_single_wedge(cfarray * &TabWT, fltarray ** Band, int s, int b); // writes in TabBand only if RealCur; normalizes the wedge if !no_norm
	void put_single_wedge(cfarray * &TabWT, fltarray ** Band, int s, int b); // reads TabBand only if RealCur; normalizes the wedge if !no_norm
	void scale_into_single_wedge(cfarray& Scale, cfarray & FBand, int b, int nd, double L1, double L2, double L3, bool backward=false);
	void meyer_transform(fltarray &Data); // Transforms data (may extend it to NewN.) into Tabcf_WT_Band
	void extern_band_norm(fltarray & Band, fltarray & Band_norm, bool backward=false);
	
	
	bool pointer;
	float* address;

public:
	inline void set_pointer(float* a){address=a; pointer=true;}
	float* get_pointer_of_band(int scale, int band);
	int get_size_of_transform();
	
	FCurvelet3D();
	~FCurvelet3D();
	
	fltarray TabSigma;		// Noise normalizing coeficients
	dblarray TabStat;

// Size and other getters
	inline int real() { return RealBand;}
	inline int nbr_angle(int s) {return TabNbrAnglePerScale(s);}
	inline int use_NAngle(int s) {return max( (real()) ? nbr_angle(s) / 2: nbr_angle(s) , 1 );}
	inline int nbr_band(int s) { return  TabNbrBandPerScale(s);}
	inline int nbr_tot_band() { return  (int) TabNbrBandPerScale.total();}
	inline int size_band_nx(int s, int b) { return TabSizeNx[s](b);}
	inline int size_band_ny(int s, int b) { return TabSizeNy[s](b);}
	inline int size_band_nz(int s, int b) { return TabSizeNz[s](b);}
	void get_band(int s, int b, fltarray **TabBand, fltarray &band);
	// float operator() (int s, int b, int i, int j, int k);

// initialize the class for given cube sizes, scales, block size, wavelet extention (for aliasing purposes), 
//		a flag to force 4sigma at finest scale, one for isotropic wavelets, and one to state that the input data is real (not complex).
// 		NbrDir2d = total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	void alloc_from_coarse(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT=True, Bool IsotropWT=False, Bool RealCur=True, type_extraction extract_type=DEF_TYPE_EXTRACTION);
	void alloc_from_fine(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool ExtendWT, Bool IsotropWT, Bool Real);
	void alloc_with_tab(int Nbr_Scale, int Nx, int Ny, int Nz, intarray & TabDir, Bool ExtendWT=True, Bool IsotropWT=False, Bool RealCur=True);

// Apply the curvelet transform and store the result in TabBand
	void cur_trans(fltarray &Data, fltarray ** &TabBand, bool allocTB);
	
// Reconstruct a cube from its curvelet transform
	void cur_recons(fltarray ** TabBand, fltarray &Data);
	
// Statistic and information tools
	inline void set_lowmem(bool lm) { lowmem=lm; }
	inline bool get_lowmem() { return lowmem; }
	bool isset_tabsigma() {return tabsigma;}
	void estim_normalization(char* Name_Imag_Out);
	void get_norm_coeff(float N_Sigma, char* Name_Imag_Out);
	void calib_noise_nsig(fltarray ** TabBand, float N_Sigma, char* Name_Imag_Out);
	//void noise_calibration(fltarray **TabBand, char* Outname);
	void extract_stat(fltarray **TabBand, char* Outname, bool not_centered=false);
	void normalize_band(fltarray &Band, int s, int b, bool inverse=false);
	void redundancy();
	inline void set_stat(bool stat) { eval_stat=stat; }
	inline void set_threshold_coarse(bool th) { th_coarse=th; }

// Filtering methods
	void filter(fltarray &Data, fltarray &Recons, float SigmaNoise, float NSigma, filter_type FilterType, 
					float Alpha, int LocalBS, float Rho, float Kmin, float Kmean, float Lmax, bool force4sigma=false);
	void threshold(fltarray **TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void threshold_single(fltarray * Band, int s, int b, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);
	void soft_threshold(fltarray **TabBand, float lvl, bool threshold_coarse=false);
	void hard_threshold(fltarray **TabBand, float lvl, bool threshold_coarse=false);
	void enhance(fltarray **TabBand, float Rho, float Kmin, float Kmean, float Lmax);
	void enhance_single(fltarray *Band, int s, int b, float Rho, float Tmin, float Tmean, float Lmax);
	void wiener(fltarray ** &TabBand, float Noiselvl, int LocalBS);
	void wiener_single(fltarray * Band, int s, int b, float Noiselvl, int LocalBS);
	void fdr(fltarray ** &TabBand, float Alpha, float SigmaNoise);
	void fdr_single(fltarray * Band, int s, int b, float Alpha, float SigmaNoise);
	void stein_block_threshold(fltarray ** &TabBand, float SigmaNoise);
	void stein_block_threshold_single(fltarray * Band, int s, int b, float SigmaNoise);

// Filtering properties
	inline void set_no_norm() {no_norm=true; }	// to prevent from normalizing (must be called after alloc)
	inline void set_no_fine(bool nf) { no_fine=nf; }
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	inline void set_positivity(bool p) { min_value=0; if(p) use_min=true; }
	inline void set_min(float mini) { min_value=mini; use_min=true; }
	inline void set_max(float maxi) { max_value=maxi; use_max=true; }
	void set_tabsigma(fltarray &TS);// Heritate TabSigma 
	void set_tabsigmanoise(int s, float sig) { TabSigmaNoise(s)=sig; }
	void set_extern_norm(fltarray ** TB) { TabBand_norm = TB;}
	
// IO tools
	void write(char *Name, fltarray ** TabBand);
	void read(char *Name, fltarray ** &TabBand);
	void write_nohead(char *Name, fltarray ** TabBand);
	void read_nohead(char *Name, fltarray ** &TabBand);

	void temp(fltarray** TabBand);
};

 
/***********************************************************************/

#endif
