
#ifndef _BFASTCUR3D_H
#define _BFASTCUR3D_H

#ifdef WINDOWS
	#define USE_OMP_BFC3D 0
#else
	#define USE_OMP_BFC3D 1
#endif

#if USE_OMP_BFC3D
	#include <omp.h>
#endif

#include <fstream>
#include <vector>
#include "IM3D_Block.h"
#include "fct.h"
#include "FloatTrans.h"

/***********************************************************************/

//v1.2 : update version when modifying. To know when to update the matlab version
// Default parameters set in the reset() method
class BFCurvelet3D_params
{
public:
	BFCurvelet3D_params(){reset_params();};
	BFCurvelet3D_params(fltarray a);
	~BFCurvelet3D_params(){};
	void import_params(BFCurvelet3D_params Q);
	void stat_import(BFCurvelet3D_params Q);
	void reset_params();

// Transform properties
	int Nx;							// Input image x size
	int Ny;							// Input image y size
	int Nz;							// Input image z size
	int NbrScale;					// Number of 3D scales
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	int BlockSize;					// Size of the blocks (default 1 full block)
	bool lapped;					// Doesn't work yet
	Bool BlockOverlap;				// Overlapping blocks
	type_extraction extract_type;	// Type of wedge extraction : backward, forward or both
	Bool RealData;					// If true, the redundancy inside the curvelet transform is used to remove

// Memory properties
	bool lowmem;					// Low memory (blockwize or directionwize) transform
	bool reuse_input_data;			// The output is stored in the input data (memory savings)

// Filtering properties
	bool no_norm;					// If set, then no normalization is done
	bool no_fine;					// The finest scale is set to 0
	bool no_coarse;					// The coarse scale is set to 0
	bool threshold_coarse;			// Treshold the coarse scale
	bool Positivity;				// The reconstruction is positive
	bool use_min;					// Use lower value for reconstruction
	bool use_max;					// Use higher value for reconstruction
	float min_value;				// Lower value for reconstruction
	float max_value;				// Higher value for reconstruction

	filter_type FilterType;			// Type of filtering
	float SigmaNoise;				// Noise level
	float NSigma;					// NSigma * SigmaNoise thresholding
	float Alpha;					// FDR parameter
	int WienerBS;					// Wiener Block Size
	float Rho;						// Contrast enhancement parameter : non linearity
	float Kmin;						// " Enhance from Kmin*SigmaNoise
	float Kmean;					// " Enhance cutoff Kmean*SigmaNoise
	float Lmax;						// " Enhance to Lmax*max(band); if(Lmax<0) Lmax*SigmaNoise
	bool force4sigma;
	
	void print();
};


class BFCurvelet3D : public FloatTrans, protected BFCurvelet3D_params 
{
// Parameters
//	BFCurvelet3D_params Param;
	
// Block management
	Block3D B3D;
	int BSx,BSy,BSz;
	void get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
		{B3D.get_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
	void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.put_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}
	void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.add_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}

// Curvelets
	FCurvelet3D *TabCur;
	
// Statistics
	bool tabsigma;					// Automatically set when the exact normalizing coef. are known
	dblarray *TabStat;				// Statistics of the bands

// Filtering properties
	bool eval_stat;					// Estimate the statistics of the bands
	fltarray *** extern_norm;		// Set an extern normalization (like for correlating the coef with those of a mask)
	
// vector form (float*) of the coefficients
	bool pointer;					// use the 'pointer mode' : the whole transform is stored on contiguous memory space (one float*)
	float* address;					// address of the transform
	fltarray *** local_TabBand;		// list of addresses mapped on the transform

public:
	void set_pointer(float* a);
	float *get_pointer(){return address;}
	float* get_pointer_of_block(int B);
	float* get_pointer_of_band(int B, int scale, int band);
	void update_TabBand(float *in);

	BFCurvelet3D();
	~BFCurvelet3D();
	
// Size and other getters
	inline int real() { return RealData;}
	inline int nbr_angle(int s) { return TabCur->nbr_angle(s);}
	inline int nbr_band(int s) { return TabCur->nbr_band(s);}
	inline int nbr_tot_band() { return 1;}
	inline int nbr_block_nx() { return B3D.nbr_block_nx();}
	inline int nbr_block_ny() { return B3D.nbr_block_ny();}
	inline int nbr_block_nz() { return B3D.nbr_block_nz();}
	inline int nbr_block()    { return B3D.nbr_block();}
	inline int block_size()   { return B3D.block_size();}
	inline int get_num_block(int Bi, int Bj, int Bk) { return Bk*nbr_block_nx()*nbr_block_ny()+Bj*nbr_block_nx()+Bi;};
	inline int size_band_nx(int s, int b) { return TabCur->size_band_nx(s,b);}
	inline int size_band_ny(int s, int b) { return TabCur->size_band_ny(s,b);}
	inline int size_band_nz(int s, int b) { return TabCur->size_band_nz(s,b);}
	inline bool isset_tabsigma() { return tabsigma;};
	void get_band(int s, int b, fltarray ***TabBand, fltarray &band);
	int size_transform();

// initialize the class for given parameters (see fct.h)
	void alloc_from_coarse(fltarray data, BFCurvelet3D_params &P);
	void alloc_from_coarse(BFCurvelet3D_params &P);
	void dealloc(fltarray*** &TabBand);

// Apply the curvelet transform and store the result in TabBand. 
	// First call MUST have alloc(TB) set to true to create TB or local_TB
	void transform(fltarray &Data, fltarray *** &TabBand, bool allocTB=false);
	// allocTB : creates the TabBand arrays, and if(pointer) allocates the float* address
	// !allocTB : if pointer : sets the address to TabBand (which MUST have been allocated in 'pointer' mode)
	void transform(fltarray &Data, float *&out, bool alloc=false);
	// alloc : transform(TB) allocates out and creates local_TabBand which follows address=out
	// !alloc : sets the address to out and thus local_TabBand will follow address (out and TB_loc must have been allocated ealier)
	
// Reconstruct a cube from its curvelet transform
	void recons(fltarray *** TabBand, fltarray &Data);
	void recons(float * in, fltarray &out);
	// Updates address and local_TB to in; apply reconstruction
	
// Statistic and information tools
	void estim_normalization(char* Name_Imag_Out);
	void extract_stat(fltarray ***TabBand, char* Outname);

// Filtering methods
	void filter(fltarray &Data, fltarray &Recons, BFCurvelet3D_params &P, char* Outname=NULL);
	void threshold(fltarray ***TabBand, BFCurvelet3D_params &P);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
	void hard_threshold(float* in, float lvl, bool threshold_coarse=false);
	void enhance(fltarray ***TabBand, BFCurvelet3D_params &P);
	void wiener(fltarray *** &TabBand, BFCurvelet3D_params &P);
	void fdr(fltarray *** &TabBand, BFCurvelet3D_params &P);
	void stein_block_threshold(fltarray *** &TabBand, BFCurvelet3D_params &P);

// Filtering properties
	double get_max_coef(int s);
	double get_max_coef();
	inline void set_stat(bool stat) { eval_stat=stat; }
	void set_tabsigmanoise(int s, float sig);
	void set_extern_norm(fltarray *** TB) {extern_norm=TB;}
	void put_extern_norm();
	
// IO tools
	void write(char *Name, fltarray *** TabBand);
	void read(char *Name, fltarray *** &TabBand, BFCurvelet3D_params &P);

	void temp(fltarray*** TabBand);
	
// Friend functions
	friend int fct3d_normalize(BFCurvelet3D_params &P, char filename[256]);
};
void fct3d_clear(vector< vector< vector< fltarray* > > > &C);
int fct3d_transform(fltarray &Data, vector< vector< vector< fltarray* > > > &vTabBand, BFCurvelet3D_params &params);
int fct3d_recons(vector< vector< vector< fltarray* > > > &vTabBand, fltarray &Data, BFCurvelet3D_params &params);
int fct3d_filter(fltarray &Data, fltarray &Recons, BFCurvelet3D_params &P);
int fct3d_threshold(vector< vector< vector< fltarray* > > > &vTabBand, BFCurvelet3D_params &params);
int fct3d_normalize(BFCurvelet3D_params &P, char filename[256]);
/***********************************************************************/

#endif

