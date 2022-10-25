/******************************************************************************
**                   Copyright (C) 2008 by CEA 
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  2008
**    
**    File:  BCurvelet3D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the 3D BeamCurvelet transform
**    ----------- 
******************************************************************************/

#ifndef _BCURVELET3D_H_
#define _BCURVELET3D_H_

#ifdef WINDOWS
	#define USE_OMP_BC 0
#else
	#define USE_OMP_BC 0
#endif

#if USE_OMP_BC
	#include <omp.h>
#endif

#include <vector>
#include "IM_IO.h"
#include "IM3D_Radon.h"
#include "SB_Filter.h"
#include "Linelet3D.h"
#include "MR_Obj.h"
#include "IM_Prob.h"
#include "MeyerWT3D.h"
#include "PoissonMeyerWT3D.h"
#include "Arr_Prob.h"
#include "MGA_Inc.h"

#define DEF_BCUR_NBR_SCALE3D 3;

/***********************************************************************/

class BCurvelet3D
{
// pointer to the filter bank to use in case of othogonal 1D WT
	SubBand1D *Ptr_SB1D;

// Base elements
	int NxCube, NyCube, NzCube;			// input image size
	Linelet3D *TabBeamlet;				// One Beamlet transform per scale
	Bool AllocClass;					// True if the CLASS has been allocated
	bool InitClass;						// True if the CLASS has been initialised
	SubBandFilter* SB1D;				// Filter Bank Class pointer
	type_linelet3d_WTtrans RidClass;	// WT can be non redundant or pyramidal
	
// Statistic elements
	fltarray TabSigma;					// Copy of the noise normalising coeficients
	fltarray* CubeSigma;				// Copy of the noise normalising coeficients
	intarray TabSizeW;					// Sizes of the wavelet scales
	bool eval_stat;
	bool no_recons;
	dblarray * TabStat;					// NScale3D arrays (block number, scale1d, stat number)
										// and last array (scale3d,scale1d,stat number)
	bool no_coarse;
	bool keep_energy;
	
//Methods
	void reset();
	void beamlet_init(Linelet3D &Rid);	// Initialize the Beamlet transform with curvelet parameters
	void alloc();
	void dealloc();
	void deinit();
	void init_recons(fltarray* &TabWaveBand, fltarray* &TabBand);
	void tab_block_size_init();			// Initialise TabBlockSize with BS at finest scale (s=0); BlockSize is odd for linelets
	void wavelet_transform(fltarray & Cube, fltarray * & TabWaveBand, bool alloc);
	void wavelet_recons(fltarray * & TabWaveBand, fltarray & Cube);
	
public:
	BCurvelet3D();
	~BCurvelet3D();

// Parameters
	type_linelet3d_WTtrans CurTrans;	// Beamlet3D transform type
	type_wavelet3D TypeW3D;				// Wavelet3D transform type
	type_pmeyer3D Type_PMW3D;			// Meyer Wavelet Poisson type
	float BlockOverlap;					// If >0, Overlapping blocks are used.
	Bool GetAutoNbScale;				// If true the number of scale is automatically NbrScale = fix( log( (3N/4) / log(2)));
	type_border Border;					// Border used in the 3D a trous WT
	type_3dcurvelet_block TypeBlock;
	bool Use_3sigma;					// Use the 3-sigma equivalent normalization coefficients
	int NbrScale3D;						// Number of scales in the 3D a trous WT
	int BlockSize;						// Block size used at finest scale
	int* TabBlockSize;					// Block size used in the Curvelet transform
										// TabBlockSize(s) = Block size of scale s
										// s = 0..MAX_SCALE-1
	
// initialize the CLASS for given cube sizes, scales and block size.
	void init(int Nx, int Ny, int Nz, int _NbrScale3D, int BS, type_wavelet3D _TypeW3D, SubBandFilter* _SB1D);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Cube, fltarray* & TabBand, bool TBalloc=true);
	void transform(fltarray * & TabWaveBand, fltarray* & TabBand, bool TBalloc=true);

// Reconstruct a cube from its curvelet transform
	void recons(fltarray* &TabBand, fltarray &Cube, bool Cubealloc=true);
	void recons(fltarray* &TabBand, fltarray* &TabWaveBand, bool TWBalloc=true);
	
// Information on bands
	int get_nbscale2d(int s3) {return TabBeamlet[s3].get_NbScale();};
	int get_nbband(int s3) {return TabBeamlet[s3].nbr_band();};
	int get_nbplane(int s3) {return TabBeamlet[s3].get_LinNbPlane();};
	int get_nbradplane(int s3) {return TabBeamlet[s3].get_RadNbPlane();};
	void get_band(int s3, int s1, fltarray *TabBand, fltarray &band);
	// Positions and length of 2d subbands
	int s2d_x0(int s3, int s);
	int s2d_y0(int s3, int s);
	int s2d_xn(int s3, int s);
	int s2d_yn(int s3, int s);
	int nbr_block(int s3) {return TabBeamlet[s3].get_NbBlock();};
	inline int nbr_block_nx (int s3) {return TabBeamlet[s3].get_NbBlock_nx();};
	inline int nbr_block_ny (int s3) {return TabBeamlet[s3].get_NbBlock_ny();};
	inline int nbr_block_nz (int s3) {return TabBeamlet[s3].get_NbBlock_nz();};
	
// Statistic and information tools
	int scale2d(int s3) {return TabBeamlet[s3].get_NbScale();};
	int num_block(int s3, int Bx, int By, int  Bz) {return TabBeamlet[s3].num_block(Bx,By,Bz);}	
	void extract_stat(fltarray *TabBand, char* Outname, bool UseCubeSigma=false, bool normalize=true);
	void values_at(fltarray *TabBand, char * filename, char* Outname); // use after normalize_self, to extract values at filename's coef list
	void noise_calibration(fltarray *TabBand, char* Outname);
	void normalize(fltarray *TabBand, fltarray *TabBandNorm, bool UseCubeSigma=false);
	void normalize_self(fltarray *TabBand, bool inverse=0, bool UseCubeSigma=false);
	void calib_noise_nsig(fltarray * TabBand, float N_Sigma, char* Name_Imag_Out);
	
// Filtering methods
	void threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false, bool UseCubeSigma=false);
	int get_BC_pos(int kx,int ky,int bx,int by,int LocalBS,int N);
	void wiener(fltarray *TabBand, float noise_lvl, int LocalBS, bool UseCubeSigma=false);
	void fdr(fltarray * &TabBand, float Alpha, float SigmaNoise);
	void filter(fltarray & Cube, fltarray & Recons, float SigmaNoise, float NSigma, filter_type FilterType, 
				float Alpha, int LocalBS, bool force4=false, bool UseCubeSigma=false, char* Outname=NULL);

// Filtering properties
	inline void set_3sigma(bool sig) { Use_3sigma=sig; }
	inline void set_stat(bool stat) { eval_stat=stat; }
	inline void set_no_recons(bool _no_recons) { no_recons=_no_recons; }
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	inline void set_keep_energy(bool nc) { keep_energy=nc; }

// IO tools
	void write(char *Name, fltarray * & Tab_WCube, bool Normalize);
	void read(char *Name, fltarray * & Tab_WCube, SubBandFilter * SB1D, bool *NormalizeInv);
	
	void temp(fltarray *TabBand);
	
	friend void bcur_transform(fltarray &Data, vector< fltarray* > &vTabBand, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D, type_pmeyer3D Type_PMW3D);
	friend void bcur_recons(vector< fltarray* > &vTabBand, fltarray &Recons, int Nx, int Ny, int Nz, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D, type_pmeyer3D Type_PMW3D);
};
void bcur_clear(vector< fltarray* > &C);
void bcur_transform(fltarray &Data, vector< fltarray* > &vTabBand, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D=W3D_MEYER, type_pmeyer3D Type_PMW3D=DEF_TYPE_PMW3D);
void bcur_recons(vector< fltarray* > &vTabBand, fltarray &Recons, int Nx, int Ny, int Nz, int NbrScale3D, int BlockSize, float BlockOverlap, type_wavelet3D W3D=W3D_MEYER, type_pmeyer3D Type_PMW3D=DEF_TYPE_PMW3D);
void bcur_threshold(vector< fltarray* > &vTabBand, float threshold, filter_type FilterType);
void bcur_filter(fltarray &Data, fltarray &recons, int NbrScale3D, int BlockSize, float BlockOverlap, float lvl, filter_type FilterType=FT_HARD, type_wavelet3D W3D=W3D_MEYER, type_pmeyer3D Type_PMW3D=DEF_TYPE_PMW3D);

/***********************************************************************/
#endif

