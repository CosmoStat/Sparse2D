
#ifndef _DIRDCT3D_H
#define _DIRDCT3D_H

#include <fstream>
#include "GlobalInc.h"
#include "IM_IO.h"
#include "FFTN.h"
#include "FFTN_3D.h"
#include "IM_Prob.h"
#include "IM_Math.h"
#include "Arr_Prob.h"
#include "MGA_Inc.h"
#include "IM3D_Block.h"

#include "IM_DCT.h"
#include "IM3D_PartialRadon.h"
#include "IM3D_IO.h"
#include "SB_Filter1D.h"
#include "WT1D_FFT.h"
#include "IM_Obj.h"
#include "fftw3.h"

class DirDCT3D
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNz;
	Bool ModifSize;
	int NewNx;	// Input cube sizes must be odd. If it is not, we etend it by one in each requiered direction
	int NewNy;
	int NewNz;
	int BlockSize;
	
	intarray TabSize;
	
	PartialRadon3D PRad3D;
	Block3D B3D;
	
// Filtering properties
	bool no_fine;
	bool no_coarse;
	
//Methods
	void transform_one_block(fltarray & Block, fltarray & BlockTrans);
	void transform_one_block(fltarray & Block);
	void recons_one_block(fltarray & Block, fltarray & BlockTrans);
	void recons_one_block(fltarray & Block);
	void get_size();
	
// Block manipulation
	void get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
		{B3D.get_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
	void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.put_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}
	void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.add_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}

public:
	DirDCT3D();
	~DirDCT3D();
	
	fltarray TabSigma;		// Copy of the noise normalising coeficients

// Size and other getters
	void init(int BS);
	void alloc();

// Properties of the blocks
	inline int nbr_block_nx() { return B3D.nbr_block_nx();}
	inline int nbr_block_ny() { return B3D.nbr_block_ny();}
	inline int nbr_block_nz() { return B3D.nbr_block_nz();}
	inline int nbr_block()    { return B3D.nbr_block();}
	
// Apply the curvelet transform and store the result in Trans
	void transform(fltarray &Data, fltarray &Trans);
	
// Reconstruct a cube from its curvelet transform
	void recons(fltarray &TabBand, fltarray &Data);

// Filtering methods
	void threshold(fltarray &TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD);

// Filtering properties
	inline void set_no_fine(bool nf) { no_fine=nf; }
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	
// IO tools
	void write(char *Name, fltarray &Trans, bool Normalize);
	void read(char *Name, fltarray &Trans, bool *NormalizeInv);

	void temp(fltarray &TabBand);
};


/***********************************************************************/

#endif
