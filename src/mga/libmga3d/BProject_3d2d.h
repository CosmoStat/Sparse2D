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
**    Date:  27/05/2009
**    
**    File:  BProject_3d2d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  Project_3d2d blockwise, followed by 2D fast curvelets
**    -----------  
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _BProject_H_
#define _BProject_H_

#include <fstream>
#include "IM_Obj.h"
#include "Project_3d2d.h"
#include "IM3D_IO.h"
#include "WT1D_FFT.h"
#include "IM3D_Block.h"
#include "MGA_Inc.h"
#include "FCur.h"
#include "SB_Filter.h"

enum type_BProject_3d2d_WTclass{LIN3D_CL_UNKNOW, 
							LIN3D_CL_PYR,
							LIN3D_CL_ORTHO};

#define NBR_LIN3D_TRANS 2
enum type_BProject_3d2d_WTtrans{LIN3D_UNKNOWN,
							LIN3D_OWT,       // standard bi-orthogonal WT
							LIN3D_PYR_FFT};  // Pyramidal FFT transform
#define DEF_LIN3D_TRANS LIN3D_OWT

/***********************************************************************/

class BProject_3d2d
{
// Base elements
	int DataNx;	// Input image size
	int DataNy;	
	int DataNz;
	int BlockSize;
	int BSx,BSy,BSz;
	bool BlockOverlap;
	int NbrScale;
	int NbrDir2d;					// Total number of directions at the coarsest resolution in 2d equivalent (minimum 8 : 2 per quad in 2d)
	type_extraction extract_type;	// Type of wedge extraction : backward, forward or both
	Bool RealBand;					// If true, the redundancy inside the curvelet transform is used to remove
									// the imaginary part.Then the imaginary part equals to zero.
	
	Block3D B3D; 
	FCUR* TabCur;
Ifloat* TabFrame;
	Project_3d2d *Project;
	type_BProject_3d2d_WTtrans LinTransf;
	
	bool tabsigma;
	dblarray *TabStat;

// Block manipulation
	void get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock, Bool Weight=False)
		{B3D.get_block_cube(Bi,Bj,Bk,Cube,CubeBlock,Weight);}
	void put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.put_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}
	void add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
		{B3D.add_block_cube(Bi,Bj,Bk,Cube,CubeBlock);}

// Filtering properties
	bool eval_stat;
	bool no_fine;
	bool no_coarse;
	bool threshold_coarse;
	bool positivity;
	bool use_min, use_max;
	float min_value, max_value;
	
public:
	BProject_3d2d();
	~BProject_3d2d();
	
// Size and other getters
//	inline int nbr_angle(int s) { return TabCur[0].nbr_angle(s);}
	inline int nbr_band(int s) { return TabCur[0].nbr_band(s);}
	inline int nbr_tot_band() { return 1;}
	inline int nbr_block_nx() { return B3D.nbr_block_nx();}
	inline int nbr_block_ny() { return B3D.nbr_block_ny();}
	inline int nbr_block_nz() { return B3D.nbr_block_nz();}
	inline int nbr_block()    { return B3D.nbr_block();}
	inline int block_size()   { return B3D.block_size();}
	inline int get_num_block(int Bi, int Bj, int Bk) { return Bk*nbr_block_nx()*nbr_block_ny()+Bj*nbr_block_nx()+Bi;};
	inline bool isset_tabsigma() { return tabsigma;};
//	void get_band(int s, int b, Ifloat ***TabBand, Ifloat &band);

// initialize the class for given parameters (see fct.h)
	void alloc_from_coarse(int Nbr_Scale, int Nx, int Ny, int Nz, int NbrDir2d, Bool BlockOverlap);
	void set_angle(int _a, int _b, int _c) {for(int i=0;i<nbr_block();i++) (Project[i]).set_angle(_a,_b,_c);}
	void dealloc(Ifloat*** &TabBand);

// Apply the curvelet transform and store the result in TabBand
	void transform(fltarray &Data, Ifloat *** &TabBand, bool allocTB=true);
	
// Reconstruct a cube from its curvelet transform
	void recons(Ifloat *** TabBand, fltarray &Data);
	
// Statistic and information tools
	inline int nbr_bands(int s) {return TabCur->nbr_band(s);}
	void get_norm_coeff(float N_Sigma, fltarray*** &TabBand, char* Name_Imag_Out, bool &allocTB);
	void noise_calibration(fltarray ***TabBand, char* Outname);
	void extract_stat(fltarray ***TabBand, char* Outname);
	void normalize_self(fltarray ***TabBand, bool inverse=false);
	void set_BlockOverlap (bool B) {BlockOverlap=B;}
	void set_LinTransf (type_BProject_3d2d_WTtrans L) {LinTransf=L;}

// Filtering methods
	void threshold(Ifloat ***TabBand, float SigmaNoise, float NSigma, filter_type FilterType=FT_HARD, bool force4sigma=false);

// Filtering properties
	inline void set_no_fine(bool nf) { no_fine=nf; }
	inline void set_no_coarse(bool nc) { no_coarse=nc; }
	inline void set_positivity(bool p) { min_value=0; if(p) use_min=true; }
	inline void set_min(float mini) { min_value=mini; use_min=true; }
	inline void set_max(float maxi) { max_value=maxi; use_max=true; }
	void set_threshold_coarse(bool th);
	
// IO tools
//	void write(char *Name, fltarray *** TabBand, bool Normalize);
//	void read(char *Name, fltarray *** &TabBand, bool *NormalizeInv);

};

 
/***********************************************************************/

#endif

