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
**    Date: 24th Oct. 2008
**    
**    File: IM3D_DCT.h
**
**    Modification history:
**
******************************************************************************
**
**    DESCRIPTION 3D Block DCT class
**    -----------  
**                 
******************************************************************************/

#ifndef _IM3D_DCT_H_
#define _IM3D_DCT_H_

#include "IM_IO.h"
#include "IM_Obj.h"
#include "IM3D_Block.h"
#include "MGA_Inc.h"
//#include "ooura_dct.h"

/***********************************************************************/

class IM3D_DCT
{
// Base elements
	int NxCube, NyCube, NzCube;			// input image size
	int BlockSize;
	bool InitClass;						// True if the CLASS has been initialised
	Block3D B3D;
	bool AllocClass;
	double ***Block3;
	double *Block1;
	double TabSigma;					// Noise level for a given Block size
	
// dct parameters
	int n1d;
	int* dct_ip;
	double* dct_w;
	bool lapped;
	int skip_order;

//Methods
	void reset();
	void alloc();
	void dealloc();
	void transform_one_block(fltarray & Block, fltarray & BlockTrans);
	void transform_one_block(fltarray & Block);
	void recons_one_block(fltarray & Block, fltarray & BlockTrans);
	void recons_one_block(fltarray & Block);
	
// Block manipulation
	void get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
		{B3D.get_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
	void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.put_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}
	void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.add_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}

public:
	IM3D_DCT();
	~IM3D_DCT();
	
// Parameters
	Bool BlockOverlap;		// If True, Overlapped blocks are used.
	type_border Border;		// Border used in the 3D a trous WT
	
// initialize the class for given cube sizes, block size
	void init(int Nx, int Ny, int Nz, int _BlockSize, Bool _BlockOverlap);
	inline void set_lapped(bool l) {lapped=l;}
// Apply the block dct transform and store the result in TabBand
	void transform(fltarray & Cube, fltarray & TabBand);
	
// Reconstruct a cube from its block dct transform
	void recons(fltarray &TabBand, fltarray &Cube);
	
// Properties of the blocks
	inline int nbr_block_nx() { return B3D.nbr_block_nx();}
	inline int nbr_block_ny() { return B3D.nbr_block_ny();}
	inline int nbr_block_nz() { return B3D.nbr_block_nz();}
	inline int nbr_block()    { return B3D.nbr_block();}
	
// Statistic and information tools
	void normalize_self(fltarray TabBand, bool inverse=false);
	void noise_calibration(fltarray &TabBand, char* Outname);
	void extract_stat(fltarray &TabBand, char* Outname);

// Filtering methods
	void threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD);
	void wiener(fltarray &TabBand, float noise_lvl, int LocalBS);
	void fdr(fltarray &TabBand, float Alpha, float SigmaNoise);
	void set_skip_order(int i) { skip_order = i;}
	
// IO tools
	void write(char *Name, fltarray & TabBand, bool Normalize);
	void read(char *Name, fltarray & TabBand, bool *NormalizeInv);
	
	void temp(fltarray &TabBand);
};

/***********************************************************************/
#endif

