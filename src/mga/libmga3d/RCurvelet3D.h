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
**    Date:  10/03/2008
**    
**    File:  RCurvelet3D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition of the 3D RidCurvelet transform
**    ----------- 
******************************************************************************/

#ifndef _RCURVELET3D_H_
#define _RCURVELET3D_H_

#ifdef WINDOWS
	#define USE_OMP_RC 0
#else
	#define USE_OMP_RC 0
#endif

#if USE_OMP_RC
	#include <omp.h>
#endif

#include "IM_IO.h"
#include "IM3D_Radon.h"
#include "SB_Filter.h"
#include "Ridgelet3D.h"
#include "MR_Obj.h"
#include "IM_Prob.h"
#include "MeyerWT3D.h"
#include "PoissonMeyerWT3D.h"
#include "Arr_Prob.h"
#include "MGA_Inc.h"
#include "Atrou3D.h"

#define DEF_NBR_RCUR_SCALE3D 3

/***********************************************************************/

class RCurvelet3D
{
// pointer to the filter bank to use in case of othogonal 1D WT
	SubBand1D *Ptr_SB1D;

// Base elements
	int NxCube, NyCube, NzCube;			// input image size
	int *CurNl, *CurNc;					// Curvelet band sizes
	int *ScaleNx, *ScaleNy, *ScaleNz;	// Wavelet band sizes
	Ridgelet3D *TabRidgelet;			// One ridgelet transform per scale
	Bool AllocClass;					// True if the CLASS has been allocated
	bool InitClass;						// True if the CLASS has been initialised
	type_ridgelet3d_WTclass RidClass;	// WT can be non redundant or pyramidal
	
// Statistic elements
	fltarray TabSigma;					// Copy of the noise normalising coeficients
	fltarray** CubeSigma;				// Copy of the noise normalising coeficients
	dblarray TabStat;
	
	bool no_coarse;
//Methods
	void reset();
	void alloc(fltarray * & TabWaveBand);
	void dealloc();
	void deinit();
	void ridgelet_init(Ridgelet3D &Rid, int BBS);	// Initialize the Ridgelet transform with curvelet parameters
	void tab_block_size_init();				// Initialise TabBlockSize with BS at finest scale (s=0)
	void wavelet_transform(fltarray & Cube, fltarray * & TabWaveBand, bool alloc=true);
	void wavelet_recons(fltarray * & TabWaveBand, fltarray & Cube);
	
public:
	RCurvelet3D();
	RCurvelet3D(RCurvelet3D *original);
	~RCurvelet3D();
	
// Parameters
	type_ridgelet3d_WTtrans CurTrans;	// Ridgelet3D transform type
	type_wavelet3D TypeW3D;				// Wavelet3D transform type
	float BlockOverlap;		// If True, Overlapped blocks are used.
	Bool GetAutoNbScale;	// If true the number of scale is automatically NbrScale = fix( log( (3N/4) / log(2)));
	type_border Border;		// Border used in the 3D a trous WT
	type_3dcurvelet_block TypeBlock;
	Bool WPTrans;			// Wavelet packets are applied on the first scale of the ridgelet transform
	bool Use_3sigma;
	int NbrScale3D;		// Number of scales in the 3D a trous WT
	int NbrScale1D;		// Number of scales used by the 1D wavelet transform
	int BlockSize;		// Block size used at finest scale
	int* TabBlockSize;	// Block size used in the Curvelet transform
						// TabBlockSize(s) = Block size of scale s
						// s = 0..MAX_SCALE-1
	
// initialize the class for given cube sizes, scales, block size, wavelet type, and a flag to force 3sigma at all scales
	void init(int Nx, int Ny, int Nz, int _NbrScale3D, int _NbrScale1D, int BS, type_wavelet3D _TypeW3D, bool use3sig=false);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Cube, fltarray* & TabBand, bool TBalloc=true);
	void transform(fltarray * & TabWaveBand, fltarray* & TabBand, bool TBalloc=true);
	
// Reconstruct a cube from its curvelet transform
	void recons(fltarray* &TabBand, fltarray &Cube, bool Cubealloc=true);
	void recons(fltarray* &TabBand, fltarray* &TabWaveBand, bool TWBalloc=true);// default : TWB is not allocated
	
// Test of wavelet transform and reconstruction
	void test_wavelet(fltarray &Cube, fltarray &Recons);
	
// Test of ridgelet transform and reconstruction
	void test_ridgelet(fltarray &Cube, fltarray &Recons);
	
// Size and position of the bands, GetBand
	int get_nbscale1d(int s3) {return TabRidgelet[s3].NbrScale;};
	int cur_nl(int s) {return CurNl[s];};
	int cur_nc(int s) {return CurNc[s];};
	int ipos(int s3d, int s1d) {return TabRidgelet[s3d].ipos(s1d);};
	int jpos(int s3d, int s1d) {return TabRidgelet[s3d].jpos(s1d);};
	int jpos(int s3d, int s1d, int B) {return TabRidgelet[s3d].jpos(s1d,B);};
	int size_nl(int s3d, int s1d) {return TabRidgelet[s3d].size_scale_nl(s1d);};
	int size_nc(int s3d, int s1d) {return TabRidgelet[s3d].size_scale_nc(s1d);};
	int size_nc(int s3d, int s1d, int B) {return TabRidgelet[s3d].rid_size(s1d);};
	int nbr_block(int s3d)    { return TabRidgelet[s3d].nbr_block();};
	inline int nbr_block_nx(int s3d) { return TabRidgelet[s3d].nbr_block_nx();}
	inline int nbr_block_ny(int s3d) { return TabRidgelet[s3d].nbr_block_ny();}
	inline int nbr_block_nz(int s3d) { return TabRidgelet[s3d].nbr_block_nz();}
	int num_block(int s3, int Bx, int By, int  Bz){ return (Bz*(nbr_block_nx(s3)*nbr_block_ny(s3))+By*nbr_block_nx(s3)+Bx);}
	void get_band(int s3, int s1, fltarray *TabBand, fltarray &band);
	
// Statistic and information tools
	void noise_calibration(fltarray *TabBand, char* Outname);
	void extract_stat(fltarray *TabBand, char* Outname, bool UseCubeSigma=false, bool normalize=true);
	void values_at(fltarray *TabBand, char * filename, char* Outname); // use after normalize_self, to extract values at filename's coef list
	void normalize(fltarray *TabBand, fltarray *TabBandNorm, bool UseCubeSigma=false);
	void normalize_self(fltarray *TabBand, bool inverse=false, bool UseCubeSigma=false);
	void calib_noise_nsig(fltarray * TabBand, float N_Sigma, char* Name_Imag_Out);
	void update_significant(fltarray *TabBand, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4sigma=false, bool UseCubeSigma=false);
	void select_significant(fltarray *TabBand, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4sigma=false, bool UseCubeSigma=false);
	void select_significant(fltarray *TabBand, intarray *support, fltarray *TabBandO);
	void support_significant(intarray *TabBand, fltarray *TabBandO, float SigmaNoise, float NSigma, bool force4sigma=false, bool UseCubeSigma=false);

// Filtering methods
	void threshold(fltarray *TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false, bool UseCubeSigma=false);
	int get_RC_pos(int kv,int kh,int by,int bx,int LocalBS,int N);
	void wiener(fltarray *TabBand, float noise_lvl, int LocalBS, bool UseCubeSigma=false);
	void fdr(fltarray * &TabBand, float Alpha, float SigmaNoise);
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	
// IO tools
	void write(char *Name, fltarray * & Tab_WCube, bool Normalize);
	void read(char *Name, fltarray * & Tab_WCube, bool *NormalizeInv);
	
	void temp(fltarray *TabBand);
};

/***********************************************************************/
#endif

